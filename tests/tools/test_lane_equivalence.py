"""The per-test equivalence check must fail when an outcome moves.

The perturbations here are the reason the check is worth running. Each one is a
way a representation change could alter behaviour while leaving the pytest
summary line untouched — a renamed test, a flipped verdict, two tests swapping
verdicts — and each must be caught. `test_a_swap_is_invisible_to_the_summary`
states the case for comparing by identity at all: it asserts, in the same test,
that the four JUnit tallies are equal while the per-test comparison fails.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping

import pytest

from tools.lane_equivalence import DEFAULT_BASELINE
from tools.lane_equivalence import EmptyLaneCollectionError
from tools.lane_equivalence import LaneSummary
from tools.lane_equivalence import MalformedLaneBaselineError
from tools.lane_equivalence import MalformedLaneJUnitError
from tools.lane_equivalence import check
from tools.lane_equivalence import compare
from tools.lane_equivalence import main
from tools.lane_equivalence import read_baseline
from tools.lane_equivalence import read_junit_outcomes

BASELINE_COMMIT = "45c5cc8aad3885c166282618b22dc0c2fd2533ba"
BASELINE_FAILURES = (
    "tests.adaptive.test_generator::test_task13f_full_continuation",
    "tests.adaptive.test_generator::test_real_active_family_stops_at_unresolved_exact_event",
)
OUTCOME_TAG = {"failed": "failure", "error": "error", "skipped": "skipped"}


def _junit(outcomes: Mapping[str, str]) -> str:
    """Render an identity/outcome map as a pytest-shaped JUnit document."""
    summary = LaneSummary.of(outcomes)
    cases = []
    for identity, outcome in outcomes.items():
        classname, name = identity.split("::", 1)
        body = "" if outcome == "passed" else f'<{OUTCOME_TAG[outcome]} message="m">detail</{OUTCOME_TAG[outcome]}>'
        cases.append(f'<testcase classname="{classname}" name="{name}" time="0.1">{body}</testcase>')
    tallies = f'errors="{summary.errors}" failures="{summary.failures}" skipped="{summary.skipped}" tests="{summary.tests}"'
    return f'<?xml version="1.0" encoding="utf-8"?><testsuites name="pytest tests"><testsuite name="pytest" {tallies}>{"".join(cases)}</testsuite></testsuites>'


def _write_junit(directory: Path, name: str, outcomes: Mapping[str, str]) -> Path:
    path = directory / name
    path.write_text(_junit(outcomes), encoding="utf-8")
    return path


@pytest.fixture()
def baseline_outcomes() -> Mapping[str, str]:
    return {
        "tests.adaptive.test_circle_oracle::test_alpha": "passed",
        "tests.adaptive.test_circle_oracle::test_beta": "passed",
        "tests.adaptive.test_generator::test_gamma": "failed",
    }


def test_the_committed_baseline_is_the_pre_inversion_record() -> None:
    baseline = read_baseline(DEFAULT_BASELINE)
    assert baseline.provenance["commit"] == BASELINE_COMMIT
    assert baseline.summary == LaneSummary(tests=145, failures=2, errors=0, skipped=0)
    assert len(baseline.outcomes) == 145
    assert sorted(identity for identity, outcome in baseline.outcomes.items() if outcome != "passed") == sorted(BASELINE_FAILURES)
    assert len(baseline.suites) == 7


def test_the_baseline_summary_is_checked_against_its_own_outcome_map(tmp_path: Path) -> None:
    """Editing a recorded verdict to silence a diff must fail loudly."""
    payload = json.loads(DEFAULT_BASELINE.read_text(encoding="utf-8"))
    payload["outcomes"][BASELINE_FAILURES[0]] = "passed"
    tampered = tmp_path / "tampered.json"
    tampered.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(MalformedLaneBaselineError, match="the record has been edited"):
        read_baseline(tampered)


def test_a_junit_report_reads_back_as_the_baseline_it_was_built_from(tmp_path: Path) -> None:
    baseline = read_baseline(DEFAULT_BASELINE)
    current = _write_junit(tmp_path, "after.xml", baseline.outcomes)
    assert read_junit_outcomes(current) == dict(baseline.outcomes)
    assert check(baseline, current).is_identical


def test_every_outcome_vocabulary_member_survives_the_round_trip(tmp_path: Path) -> None:
    outcomes = {
        "m::a": "passed",
        "m::b": "failed",
        "m::c": "error",
        "m::d": "skipped",
    }
    assert read_junit_outcomes(_write_junit(tmp_path, "v.xml", outcomes)) == outcomes


def test_a_renamed_test_is_reported_as_a_disappearance_and_an_appearance(baseline_outcomes: Mapping[str, str]) -> None:
    current = dict(baseline_outcomes)
    current["tests.adaptive.test_circle_oracle::test_alpha_renamed"] = current.pop("tests.adaptive.test_circle_oracle::test_alpha")
    equivalence = compare(baseline_outcomes, current)
    assert not equivalence.is_identical
    assert equivalence.only_in_baseline == ("tests.adaptive.test_circle_oracle::test_alpha",)
    assert equivalence.only_in_current == ("tests.adaptive.test_circle_oracle::test_alpha_renamed",)
    assert equivalence.outcome_changes == ()


def test_a_flipped_outcome_is_reported(baseline_outcomes: Mapping[str, str]) -> None:
    current = dict(baseline_outcomes)
    current["tests.adaptive.test_generator::test_gamma"] = "passed"
    equivalence = compare(baseline_outcomes, current)
    assert not equivalence.is_identical
    assert equivalence.outcome_changes == (("tests.adaptive.test_generator::test_gamma", "failed", "passed"),)


def test_a_swap_is_invisible_to_the_summary(baseline_outcomes: Mapping[str, str]) -> None:
    """The whole reason the comparison is per test and not per count."""
    current = dict(baseline_outcomes)
    current["tests.adaptive.test_circle_oracle::test_alpha"] = "failed"
    current["tests.adaptive.test_generator::test_gamma"] = "passed"
    assert LaneSummary.of(current) == LaneSummary.of(baseline_outcomes)
    equivalence = compare(baseline_outcomes, current)
    assert not equivalence.is_identical
    assert len(equivalence.outcome_changes) == 2


def test_an_empty_collection_is_never_agreement(tmp_path: Path) -> None:
    empty = tmp_path / "empty.xml"
    empty.write_text('<?xml version="1.0"?><testsuites><testsuite name="pytest" errors="0" failures="0" skipped="0" tests="0" /></testsuites>', encoding="utf-8")
    with pytest.raises(EmptyLaneCollectionError):
        read_junit_outcomes(empty)


def test_a_repeated_identity_is_rejected(tmp_path: Path) -> None:
    """A duplicate would silently drop one of the two from the comparison."""
    repeated = tmp_path / "repeated.xml"
    repeated.write_text(
        '<?xml version="1.0"?><testsuites><testsuite name="pytest"><testcase classname="m" name="a" /><testcase classname="m" name="a" /></testsuite></testsuites>',
        encoding="utf-8",
    )
    with pytest.raises(MalformedLaneJUnitError, match="duplicate testcase identity"):
        read_junit_outcomes(repeated)


def test_a_document_that_is_not_junit_is_rejected(tmp_path: Path) -> None:
    wrong = tmp_path / "wrong.xml"
    wrong.write_text("<coverage><line/></coverage>", encoding="utf-8")
    with pytest.raises(MalformedLaneJUnitError, match="root must be"):
        read_junit_outcomes(wrong)


def test_a_testcase_without_identity_is_rejected(tmp_path: Path) -> None:
    anonymous = tmp_path / "anonymous.xml"
    anonymous.write_text('<?xml version="1.0"?><testsuite name="pytest"><testcase name="a" /></testsuite>', encoding="utf-8")
    with pytest.raises(MalformedLaneJUnitError, match="non-empty classname and name"):
        read_junit_outcomes(anonymous)


def test_the_cli_separates_agreement_divergence_and_malformation(tmp_path: Path) -> None:
    baseline = read_baseline(DEFAULT_BASELINE)
    identical = _write_junit(tmp_path, "identical.xml", baseline.outcomes)
    assert main(["--baseline", str(DEFAULT_BASELINE), "--current", str(identical)]) == 0

    moved = dict(baseline.outcomes)
    moved[BASELINE_FAILURES[0]] = "passed"
    diverged = _write_junit(tmp_path, "diverged.xml", moved)
    assert main(["--baseline", str(DEFAULT_BASELINE), "--current", str(diverged)]) == 1

    broken = tmp_path / "broken.xml"
    broken.write_text("not xml at all", encoding="utf-8")
    assert main(["--baseline", str(DEFAULT_BASELINE), "--current", str(broken)]) == 2
