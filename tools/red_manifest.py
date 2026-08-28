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

REQUIRED_KEYS = ("match", "count", "reason", "closes_with")


class MalformedManifestError(Exception):
    """The manifest is missing a required key or is not valid JSON."""


@dataclass(frozen=True)
class ManifestViolation:
    """One divergence between the suite and the manifest.

    Attributes:
        kind: ``unexpected-red`` | ``expected-red-went-green`` |
            ``count-mismatch``.
        detail: Human-readable statement naming the tests or the entry.
    """

    kind: str
    detail: str


def _failed_ids(junit: Path) -> List[str]:
    """Return ``classname::name`` for each testcase with a failure or error."""
    root = ET.parse(junit).getroot()
    out: List[str] = []
    for case in root.iter("testcase"):
        if any(child.tag in ("failure", "error") for child in case):
            out.append(f"{case.get('classname', '')}::{case.get('name', '')}")
    return out


def check(junit: Path, manifest: Path) -> List[ManifestViolation]:
    """Diff the suite's reds against the manifest, both directions.

    Args:
        junit: Path to a pytest ``--junitxml`` report.
        manifest: Path to the red-manifest JSON.

    Returns:
        Violations; empty means the invariant holds.

    Raises:
        MalformedManifestError: An entry lacks one of `REQUIRED_KEYS`.
    """
    entries = json.loads(manifest.read_text()).get("expected_red", [])
    for entry in entries:
        missing = [key for key in REQUIRED_KEYS if key not in entry]
        if missing:
            raise MalformedManifestError(f"manifest entry {entry!r} lacks {missing}")
    reds = _failed_ids(junit)
    violations: List[ManifestViolation] = []
    unclaimed = list(reds)
    for entry in entries:
        pattern = re.compile(entry["match"])
        matched = [red for red in unclaimed if pattern.search(red)]
        for red in matched:
            unclaimed.remove(red)
        if not matched and entry["count"] > 0:
            violations.append(
                ManifestViolation(
                    "expected-red-went-green",
                    f"{entry['match']} (reason: {entry['reason']}) — a finding, not a pass",
                )
            )
        elif matched and len(matched) != entry["count"]:
            violations.append(
                ManifestViolation(
                    "count-mismatch",
                    f"{entry['match']}: expected {entry['count']} red, saw {len(matched)}",
                )
            )
    violations.extend(ManifestViolation("unexpected-red", red) for red in unclaimed)
    return violations


def main(argv: List[str] | None = None) -> int:
    """Run the command-line red-manifest check."""
    parser = argparse.ArgumentParser(description="Diff suite reds against docs/red_manifest.json.")
    parser.add_argument("junit", type=Path)
    parser.add_argument("--manifest", type=Path, default=Path("docs/red_manifest.json"))
    args = parser.parse_args(argv)
    violations = check(args.junit, args.manifest)
    for violation in violations:
        print(f"{violation.kind}: {violation.detail}")
    if not violations:
        print("red set == manifest, both directions")
    return 1 if violations else 0


if __name__ == "__main__":
    sys.exit(main())
