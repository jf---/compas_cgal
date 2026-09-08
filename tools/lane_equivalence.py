"""Per-test equivalence between a suite run and a frozen behavioural baseline.

A representation change must move **no** test outcome. Summary counts cannot
state that: two tests swapping outcomes leaves ``tests``, ``failures``,
``errors`` and ``skipped`` byte-identical, and a suite that collected nothing at
all still exits 0. This module therefore compares **by test identity**, and
refuses to report equivalence over an empty or malformed collection.

The frozen record is a distilled JSON manifest rather than the raw
``--junitxml`` document, because the raw document lives under the gitignored
``build/`` tree and cannot survive a clean checkout. The manifest carries its
provenance (the commit, the carrier state, the exact command) alongside the
identity/outcome map, and its ``summary`` block is cross-checked against that
map on load, so hand-editing an outcome to silence a diff fails loudly instead
of passing quietly.

Stage 2 of the number-type coherence work is the first user: the manifest at
``tests/build/stage2_station_lane_baseline.json`` was captured at ``45c5cc8a``
with ``StationEventSource2`` still carrying decimal text, and the inversion must
reproduce it exactly. Nothing here is stage-specific — a later stage freezes its
own manifest and passes ``--baseline``.

```bash
# run the baseline's own suites and diff the result
pixi run -e default python -m tools.lane_equivalence --run

# diff a junit report that already exists
pixi run -e default python -m tools.lane_equivalence --current build/stage2/after.xml
```
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Dict
from typing import List
from typing import Mapping
from typing import Tuple

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BASELINE = REPOSITORY_ROOT / "tests" / "build" / "stage2_station_lane_baseline.json"
DEFAULT_CURRENT = REPOSITORY_ROOT / "build" / "stage2" / "after.xml"

REQUIRED_BASELINE_KEYS = ("provenance", "suites", "summary", "outcomes")
REQUIRED_PROVENANCE_KEYS = (
    "stage",
    "commit",
    "captured",
    "carrier_state",
    "command",
    "pytest_summary",
    "source",
)
REQUIRED_SUMMARY_KEYS = ("tests", "failures", "errors", "skipped")

# pytest's own vocabulary, not the JUnit tag names: a `<failure>` child is a
# `failed` test and a `<error>` child is an `error`, so a report read back here
# reads the way the summary line does.
PASSED = "passed"
FAILED = "failed"
ERROR = "error"
SKIPPED = "skipped"
OUTCOMES = (PASSED, FAILED, ERROR, SKIPPED)

# JUnit child tag -> outcome. Order matters when a testcase carries more than
# one: an execution error outranks an assertion failure, which outranks a skip.
TAG_OUTCOME_PRECEDENCE = ((ERROR, "error"), (FAILED, "failure"), (SKIPPED, "skipped"))


class MalformedLaneJUnitError(Exception):
    """The JUnit document is malformed or lacks testcase identity."""


class MalformedLaneBaselineError(Exception):
    """The frozen baseline manifest does not satisfy its complete schema."""


class EmptyLaneCollectionError(Exception):
    """The JUnit report contains no testcase.

    A pytest run that collects nothing exits 0, so an empty report must never be
    read as agreement with the baseline.
    """


class LaneSuiteRunError(Exception):
    """The baseline's own suites could not be run."""


@dataclass(frozen=True)
class LaneSummary:
    """The four JUnit tallies, kept together because they are read together."""

    tests: int
    failures: int
    errors: int
    skipped: int

    @classmethod
    def of(cls, outcomes: Mapping[str, str]) -> LaneSummary:
        """Derive the tallies an outcome map implies.

        Args:
            outcomes: Test identity to outcome.

        Returns:
            The summary a JUnit report over exactly these outcomes would carry.
        """
        values = list(outcomes.values())
        return cls(
            tests=len(values),
            failures=values.count(FAILED),
            errors=values.count(ERROR),
            skipped=values.count(SKIPPED),
        )


@dataclass(frozen=True)
class LaneBaseline:
    """A frozen per-test behavioural record with the provenance to read it by.

    Attributes:
        path: Where the manifest was loaded from.
        provenance: Commit, capture time, carrier state and exact command.
        suites: The test files the record covers, repository-relative.
        summary: The tallies recorded at capture.
        outcomes: Test identity to outcome, one entry per collected test.
    """

    path: Path
    provenance: Mapping[str, str]
    suites: Tuple[str, ...]
    summary: LaneSummary
    outcomes: Mapping[str, str]


@dataclass(frozen=True)
class LaneEquivalence:
    """The three ways a run can differ from its baseline, kept apart.

    A count of changes is not a diagnosis. A test that disappeared, a test that
    appeared, and a test that changed verdict have different causes and
    different repairs, so they are never summed into one number.

    Attributes:
        only_in_baseline: Identities the baseline has and the run does not.
        only_in_current: Identities the run has and the baseline does not.
        outcome_changes: (identity, baseline outcome, current outcome) triples.
    """

    only_in_baseline: Tuple[str, ...]
    only_in_current: Tuple[str, ...]
    outcome_changes: Tuple[Tuple[str, str, str], ...]

    @property
    def is_identical(self) -> bool:
        """Whether every test in the run has the baseline's identity and verdict."""
        return not (self.only_in_baseline or self.only_in_current or self.outcome_changes)

    def report(self) -> str:
        """Render the diff, one finding per line, most diagnostic first."""
        lines: List[str] = []
        for identity in self.only_in_baseline:
            lines.append(f"disappeared: {identity}")
        for identity in self.only_in_current:
            lines.append(f"appeared: {identity}")
        for identity, was, now in self.outcome_changes:
            lines.append(f"outcome-changed: {identity}: {was} -> {now}")
        if not lines:
            lines.append("per-test outcomes identical to the baseline")
        lines.append(f"ids only in baseline: {len(self.only_in_baseline)}, ids only in current: {len(self.only_in_current)}, outcome changes: {len(self.outcome_changes)}")
        return "\n".join(lines)


def read_junit_outcomes(junit: Path) -> Dict[str, str]:
    """Read a pytest ``--junitxml`` report into an identity/outcome map.

    Args:
        junit: Path to the report.

    Returns:
        Test identity (``classname::name``) to outcome.

    Raises:
        MalformedLaneJUnitError: The report is unreadable, is not JUnit, lacks
            testcase identity, or repeats an identity — a repeat would silently
            drop a test from the comparison.
        EmptyLaneCollectionError: The report contains no testcase.
    """
    try:
        root = ET.parse(junit).getroot()
    except OSError as error:
        raise MalformedLaneJUnitError(f"{junit}: cannot read JUnit XML: {error}") from error
    except ET.ParseError as error:
        raise MalformedLaneJUnitError(f"{junit}: invalid XML: {error}") from error
    if root.tag not in ("testsuite", "testsuites"):
        raise MalformedLaneJUnitError(f"{junit}: root must be <testsuite> or <testsuites>, got <{root.tag}>")
    outcomes: Dict[str, str] = {}
    for case in root.iter("testcase"):
        classname = case.get("classname")
        name = case.get("name")
        if not classname or not name:
            context = ET.tostring(case, encoding="unicode")
            raise MalformedLaneJUnitError(f"{junit}: testcase requires non-empty classname and name: {context}")
        identity = f"{classname}::{name}"
        if identity in outcomes:
            raise MalformedLaneJUnitError(f"{junit}: duplicate testcase identity: {identity}")
        tags = {child.tag for child in case}
        outcome = PASSED
        for candidate, tag in TAG_OUTCOME_PRECEDENCE:
            if tag in tags:
                outcome = candidate
                break
        outcomes[identity] = outcome
    if not outcomes:
        raise EmptyLaneCollectionError(f"{junit}: contains no testcase; a run that collects nothing still exits 0")
    return outcomes


def read_baseline(manifest: Path) -> LaneBaseline:
    """Load and fully validate a frozen baseline manifest.

    The ``summary`` block is cross-checked against ``outcomes``, so an edited
    verdict cannot pass as a record.

    Args:
        manifest: Path to the JSON manifest.

    Returns:
        The validated baseline.

    Raises:
        MalformedLaneBaselineError: Any schema or internal-consistency breach.
    """
    try:
        payload: object = json.loads(manifest.read_text(encoding="utf-8"))
    except OSError as error:
        raise MalformedLaneBaselineError(f"{manifest}: cannot read baseline: {error}") from error
    except json.JSONDecodeError as error:
        raise MalformedLaneBaselineError(f"{manifest}: invalid JSON: {error}") from error
    if not isinstance(payload, dict):
        raise MalformedLaneBaselineError(f"{manifest}: root must be an object, got {type(payload).__name__}")
    if set(payload) != set(REQUIRED_BASELINE_KEYS):
        missing = sorted(set(REQUIRED_BASELINE_KEYS) - set(payload))
        unexpected = sorted(set(payload) - set(REQUIRED_BASELINE_KEYS))
        raise MalformedLaneBaselineError(f"{manifest}: root keys differ: missing={missing}, unexpected={unexpected}")

    provenance = payload["provenance"]
    if not isinstance(provenance, dict) or set(provenance) != set(REQUIRED_PROVENANCE_KEYS):
        raise MalformedLaneBaselineError(f"{manifest}: provenance must be an object with exactly {sorted(REQUIRED_PROVENANCE_KEYS)}")
    for key, value in provenance.items():
        if not isinstance(value, str) or not value.strip():
            raise MalformedLaneBaselineError(f"{manifest}: provenance.{key} must be a non-empty string")

    raw_suites = payload["suites"]
    if not isinstance(raw_suites, list) or not raw_suites:
        raise MalformedLaneBaselineError(f"{manifest}: suites must be a non-empty list")
    for index, suite in enumerate(raw_suites):
        if not isinstance(suite, str) or not suite.strip():
            raise MalformedLaneBaselineError(f"{manifest}: suites[{index}] must be a non-empty string")

    raw_summary = payload["summary"]
    if not isinstance(raw_summary, dict) or set(raw_summary) != set(REQUIRED_SUMMARY_KEYS):
        raise MalformedLaneBaselineError(f"{manifest}: summary must be an object with exactly {sorted(REQUIRED_SUMMARY_KEYS)}")
    for key, value in raw_summary.items():
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise MalformedLaneBaselineError(f"{manifest}: summary.{key} must be a non-negative non-boolean integer")

    raw_outcomes = payload["outcomes"]
    if not isinstance(raw_outcomes, dict) or not raw_outcomes:
        raise MalformedLaneBaselineError(f"{manifest}: outcomes must be a non-empty object")
    for identity, outcome in raw_outcomes.items():
        if "::" not in identity:
            raise MalformedLaneBaselineError(f"{manifest}: outcome key {identity!r} is not a classname::name identity")
        if outcome not in OUTCOMES:
            raise MalformedLaneBaselineError(f"{manifest}: outcomes[{identity!r}] is {outcome!r}, not one of {list(OUTCOMES)}")

    summary = LaneSummary(**raw_summary)
    derived = LaneSummary.of(raw_outcomes)
    if summary != derived:
        raise MalformedLaneBaselineError(f"{manifest}: summary {summary} disagrees with the outcome map {derived}; the record has been edited")

    return LaneBaseline(
        path=manifest,
        provenance=provenance,
        suites=tuple(raw_suites),
        summary=summary,
        outcomes=dict(raw_outcomes),
    )


def compare(baseline: Mapping[str, str], current: Mapping[str, str]) -> LaneEquivalence:
    """Diff two identity/outcome maps, per test.

    Args:
        baseline: The frozen record.
        current: The run under test.

    Returns:
        The three difference sets, each sorted by identity.
    """
    only_in_baseline = tuple(sorted(set(baseline) - set(current)))
    only_in_current = tuple(sorted(set(current) - set(baseline)))
    outcome_changes = tuple((identity, baseline[identity], current[identity]) for identity in sorted(set(baseline) & set(current)) if baseline[identity] != current[identity])
    return LaneEquivalence(only_in_baseline, only_in_current, outcome_changes)


def run_baseline_suites(baseline: LaneBaseline, junit: Path) -> None:
    """Re-run the baseline's own suites, writing ``junit``.

    Delegates to ``pixi run -e default pytest``, the project's one declared way
    to run the suite, so the editable-rebuild guard the task carries applies
    here too. The subprocess's exit status is deliberately ignored — reds are
    read from the report, and a run that collected nothing would exit 0.

    Args:
        baseline: The record naming the suites to run.
        junit: Where to write the report.

    Raises:
        LaneSuiteRunError: pytest wrote no report at all.
    """
    junit.parent.mkdir(parents=True, exist_ok=True)
    if junit.exists():
        junit.unlink()
    command = [
        "pixi",
        "run",
        "-e",
        "default",
        "pytest",
        "--",
        *baseline.suites,
        "-n",
        "auto",
        "-q",
        f"--junitxml={junit}",
    ]
    print(f"running: {' '.join(command)}", flush=True)
    subprocess.run(command, cwd=REPOSITORY_ROOT, check=False)
    if not junit.exists():
        raise LaneSuiteRunError(f"{junit}: pytest produced no JUnit report; the run did not start")


def check(baseline: LaneBaseline, current: Path) -> LaneEquivalence:
    """Compare a JUnit report against a frozen baseline, per test.

    Args:
        baseline: The validated frozen record.
        current: Path to the run's JUnit report.

    Returns:
        The per-test difference.

    Raises:
        MalformedLaneJUnitError: The report is malformed.
        EmptyLaneCollectionError: The report contains no testcase.
    """
    return compare(baseline.outcomes, read_junit_outcomes(current))


def main(argv: List[str] | None = None) -> int:
    """Run the command-line per-test equivalence check."""
    parser = argparse.ArgumentParser(description="Diff a suite run against a frozen per-test baseline.")
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--current", type=Path, default=DEFAULT_CURRENT)
    parser.add_argument(
        "--run",
        action="store_true",
        help="run the baseline's own suites into --current first",
    )
    args = parser.parse_args(argv)
    try:
        baseline = read_baseline(args.baseline)
    except MalformedLaneBaselineError as error:
        print(f"malformed-baseline: {error}", file=sys.stderr)
        return 2
    try:
        if args.run:
            run_baseline_suites(baseline, args.current)
        equivalence = check(baseline, args.current)
    except LaneSuiteRunError as error:
        print(f"suite-run-failed: {error}", file=sys.stderr)
        return 2
    except MalformedLaneJUnitError as error:
        print(f"malformed-junit: {error}", file=sys.stderr)
        return 2
    except EmptyLaneCollectionError as error:
        print(f"empty-collection: {error}", file=sys.stderr)
        return 2
    print(f"baseline: {baseline.provenance['commit']} ({baseline.provenance['carrier_state']}), {baseline.summary.tests} tests")
    print(equivalence.report())
    return 0 if equivalence.is_identical else 1


if __name__ == "__main__":
    sys.exit(main())
