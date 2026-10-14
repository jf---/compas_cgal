"""Census coverage certificate: every ledger file has exactly one record, read to the depth its slice requires.

Exit status 0 only when:
  * every slice output exists;
  * every file assigned to a slice has exactly one ``### <path> | lines N | read R`` record in that slice's output;
  * the recorded line count equals the ledger's;
  * full-depth slices (everything except the validation slice) cover every line 1..N with their read ranges;
  * the validation slice records either full ranges or ``structural``, and the files the brief names are read in full.
Prints one line per gap; the gaps are the re-dispatch list.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).parent
STRUCTURAL_SLICE = "v_validation_tests_tools_build"
MUST_READ_FULLY = {
    "benchmarks/quality.py",
    "benchmarks/quality_observations.py",
    "benchmarks/survey.py",
    "benchmarks/depletion.py",
    "benchmarks/gate.py",
    "benchmarks/coverage.py",
    "benchmarks/pathmetrics.py",
    "benchmarks/exceedance.py",
    "benchmarks/spec.py",
    "tests/benchmarks/test_quality.py",
}
HEADER = re.compile(r"^### (?P<path>\S+) \| lines (?P<lines>\d+) \| read (?P<read>.+?)\s*$")


def covered(read: str, total: int) -> bool:
    """Whether comma-separated ``a-b`` ranges cover every line ``1..total``."""
    lines: set[int] = set()
    for part in read.split(","):
        match = re.fullmatch(r"\s*(\d+)\s*-\s*(\d+)\s*", part)
        if match is None:
            return False
        lines.update(range(int(match.group(1)), int(match.group(2)) + 1))
    return total == 0 or set(range(1, total + 1)) <= lines


gaps: list[str] = []
for slice_file in sorted((HERE / "slices").glob("*.tsv")):
    name = slice_file.stem
    assigned = {path: int(count) for count, path in (row.split("\t") for row in slice_file.read_text().splitlines())}
    output = HERE / "out" / f"{name}.md"
    if not output.exists():
        gaps.append(f"{name}: output missing ({len(assigned)} files unrecorded)")
        continue
    records: dict[str, list[tuple[int, str]]] = {}
    for line in output.read_text().splitlines():
        match = HEADER.match(line)
        if match:
            records.setdefault(match["path"], []).append((int(match["lines"]), match["read"]))
    for path, total in assigned.items():
        found = records.get(path, [])
        if len(found) != 1:
            gaps.append(f"{name}: {path}: {len(found)} records")
            continue
        lines, read = found[0]
        if lines != total:
            gaps.append(f"{name}: {path}: recorded {lines} lines, ledger {total}")
        if name == STRUCTURAL_SLICE and read.strip() == "structural":
            if path in MUST_READ_FULLY:
                gaps.append(f"{name}: {path}: must be read fully, recorded structural")
            continue
        if not covered(read, total):
            gaps.append(f"{name}: {path}: read '{read}' does not cover 1-{total}")

commits = re.findall(r"^([0-9a-f]{7,40}) \d{4}-\d{2}-\d{2} ", (HERE / "slices" / "g7_unmerged_commits.txt").read_text(), re.MULTILINE)
g7 = HERE / "out" / "g7_unmerged_commits.md"
recorded = set(re.findall(r"^### commit ([0-9a-f]{7,40})", g7.read_text(), re.MULTILINE)) if g7.exists() else set()
for sha in commits:
    if not any(r.startswith(sha) or sha.startswith(r) for r in recorded):
        gaps.append(f"g7_unmerged_commits: commit {sha}: no record")

for gap in gaps:
    print(gap)
print(f"{len(gaps)} gap(s)")
sys.exit(1 if gaps else 0)
