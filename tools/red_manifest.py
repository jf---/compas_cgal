"""Invariant I3: compare deliberate suite reds with their manifest."""

from __future__ import annotations

import argparse
import json
import re
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import List
from typing import Literal

REQUIRED_KEYS = ("match", "count", "reason", "closes_with")
ViolationKind = Literal[
    "unexpected-red",
    "expected-red-went-green",
    "count-mismatch",
    "overlapping-red",
]
VIOLATION_KIND_ORDER: dict[ViolationKind, int] = {
    "expected-red-went-green": 0,
    "count-mismatch": 1,
    "overlapping-red": 2,
    "unexpected-red": 3,
}


class MalformedManifestError(Exception):
    """The manifest does not satisfy the complete red-manifest schema."""


class MalformedJUnitError(Exception):
    """The JUnit document is malformed or lacks testcase identity."""


@dataclass(frozen=True)
class ManifestEntry:
    """One validated expected-red ownership rule."""

    match: str
    count: int
    reason: str
    closes_with: str
    pattern: re.Pattern[str]

    @classmethod
    def build(cls, raw: object, index: int) -> ManifestEntry:
        """Validate and compile one manifest entry.

        Args:
            raw: Decoded JSON value for the entry.
            index: Zero-based entry position for diagnostics.

        Returns:
            A fully validated ownership rule.

        Raises:
            MalformedManifestError: The entry violates the schema.
        """
        if not isinstance(raw, dict):
            raise MalformedManifestError(f"expected_red[{index}] must be an object, got {type(raw).__name__}")
        actual_keys = set(raw)
        required_keys = set(REQUIRED_KEYS)
        if actual_keys != required_keys:
            missing = sorted(required_keys - actual_keys)
            unexpected = sorted(actual_keys - required_keys)
            raise MalformedManifestError(f"expected_red[{index}] keys differ: missing={missing}, unexpected={unexpected}")

        match = _nonempty_manifest_text(raw["match"], index, "match")
        count = raw["count"]
        if isinstance(count, bool) or not isinstance(count, int) or count <= 0:
            raise MalformedManifestError(f"expected_red[{index}].count must be a positive non-boolean integer")
        reason = _nonempty_manifest_text(raw["reason"], index, "reason")
        closes_with = _nonempty_manifest_text(raw["closes_with"], index, "closes_with")
        try:
            pattern = re.compile(match)
        except re.error as error:
            raise MalformedManifestError(f"expected_red[{index}].match is not a valid regex: {error}") from error
        return cls(match, count, reason, closes_with, pattern)


def _nonempty_manifest_text(value: object, index: int, key: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise MalformedManifestError(f"expected_red[{index}].{key} must be a non-empty string")
    return value


@dataclass(frozen=True)
class ManifestViolation:
    """One divergence between the suite and the manifest.

    Attributes:
        kind: Closed `ViolationKind` vocabulary.
        detail: Human-readable statement naming the tests or the entry.
    """

    kind: ViolationKind
    detail: str


def _failed_ids(junit: Path) -> List[str]:
    """Return ``classname::name`` for each testcase with a failure or error."""
    try:
        root = ET.parse(junit).getroot()
    except ET.ParseError as error:
        raise MalformedJUnitError(f"{junit}: invalid XML: {error}") from error
    if root.tag not in ("testsuite", "testsuites"):
        raise MalformedJUnitError(f"{junit}: root must be <testsuite> or <testsuites>, got <{root.tag}>")
    out: List[str] = []
    for case in root.iter("testcase"):
        classname = case.get("classname")
        name = case.get("name")
        if not classname or not name:
            context = ET.tostring(case, encoding="unicode")
            raise MalformedJUnitError(f"{junit}: testcase requires non-empty classname and name: {context}")
        if any(child.tag in ("failure", "error") for child in case):
            out.append(f"{classname}::{name}")
    return out


def _load_manifest(manifest: Path) -> List[ManifestEntry]:
    try:
        payload: object = json.loads(manifest.read_text())
    except json.JSONDecodeError as error:
        raise MalformedManifestError(f"{manifest}: invalid JSON: {error}") from error
    if not isinstance(payload, dict):
        raise MalformedManifestError(f"{manifest}: root must be an object, got {type(payload).__name__}")
    if set(payload) != {"expected_red"}:
        raise MalformedManifestError(f"{manifest}: root must contain only the expected_red key")
    raw_entries = payload["expected_red"]
    if not isinstance(raw_entries, list):
        raise MalformedManifestError(f"{manifest}: expected_red must be a list, got {type(raw_entries).__name__}")
    return [ManifestEntry.build(raw, index) for index, raw in enumerate(raw_entries)]


def check(junit: Path, manifest: Path) -> List[ManifestViolation]:
    """Diff the suite's reds against the manifest, both directions.

    Args:
        junit: Path to a pytest ``--junitxml`` report.
        manifest: Path to the red-manifest JSON.

    Returns:
        Violations; empty means the invariant holds.

    Raises:
        MalformedManifestError: The manifest violates its schema.
        MalformedJUnitError: The JUnit report is malformed.
    """
    entries = sorted(
        _load_manifest(manifest),
        key=lambda entry: (
            entry.match,
            entry.reason,
            entry.closes_with,
            entry.count,
        ),
    )
    reds = _failed_ids(junit)
    violations: List[ManifestViolation] = []
    owners: List[List[ManifestEntry]] = [[] for _ in reds]
    for entry in entries:
        matched_indices = [index for index, red in enumerate(reds) if entry.pattern.search(red)]
        for index in matched_indices:
            owners[index].append(entry)
        if not matched_indices:
            violations.append(
                ManifestViolation(
                    "expected-red-went-green",
                    f"{entry.match} (reason: {entry.reason}) — a finding, not a pass",
                )
            )
        elif len(matched_indices) != entry.count:
            violations.append(
                ManifestViolation(
                    "count-mismatch",
                    f"{entry.match}: expected {entry.count} red, saw {len(matched_indices)}",
                )
            )
    for red, red_owners in zip(reds, owners):
        if not red_owners:
            violations.append(ManifestViolation("unexpected-red", red))
        elif len(red_owners) > 1:
            patterns = sorted(owner.match for owner in red_owners)
            violations.append(
                ManifestViolation(
                    "overlapping-red",
                    f"{red}: matched {len(patterns)} manifest entries {patterns!r}",
                )
            )
    return sorted(
        violations,
        key=lambda violation: (
            VIOLATION_KIND_ORDER[violation.kind],
            violation.detail,
        ),
    )


def main(argv: List[str] | None = None) -> int:
    """Run the command-line red-manifest check."""
    parser = argparse.ArgumentParser(description="Diff suite reds against docs/red_manifest.json.")
    parser.add_argument("junit", type=Path)
    parser.add_argument("--manifest", type=Path, default=Path("docs/red_manifest.json"))
    args = parser.parse_args(argv)
    try:
        violations = check(args.junit, args.manifest)
    except MalformedManifestError as error:
        print(f"malformed-manifest: {error}", file=sys.stderr)
        return 2
    except MalformedJUnitError as error:
        print(f"malformed-junit: {error}", file=sys.stderr)
        return 2
    for violation in violations:
        print(f"{violation.kind}: {violation.detail}")
    if not violations:
        print("red set == manifest, both directions")
    return 1 if violations else 0


if __name__ == "__main__":
    sys.exit(main())
