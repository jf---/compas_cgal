"""Traceability certificate: every census finding has exactly one disposition in the plan's matrix.

Finding ids are derived from the census records, independently of the matrix:
  * register rows ``F-<slice>-<n> | ...`` in ``out/*.md``;
  * from each record's ``## Slice summary`` (up to the first ``### <...> | ...`` record header):
    ``<slice>-T<n>`` for each top-level numbered item under "Top findings",
    ``<slice>-U<n>`` for each top-level bullet (or table body row) under the "capabilities" heading,
    ``<slice>-P<n>`` for each top-level bullet under the "wrong or incomplete" heading,
    numbered in order of appearance. ``<slice>`` is the record stem up to its first underscore.

Exit status 0 only when every id appears in exactly one matrix row, no matrix row names an unknown id,
no row is ``UNMAPPED``, and the ``## Unmapped`` section lists exactly the ``UNMAPPED`` rows. Prints the counts by disposition. ``--dump`` prints the derived items instead.
"""

from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).parent
OUT = HERE / "out"
MATRIX = HERE.parents[1] / "plans" / "2026-10-13-product-pipeline-traceability.md"
DISPOSITIONS = {"task", "constraint", "frozen", "parked", "closed", "superseded", "UNMAPPED"}
REGISTER = re.compile(r"^(F-[A-Za-z0-9_]+-\d+) \|")
RECORD_HEADER = re.compile(r"^###+ .* \| ")
HEADING = re.compile(r"^(#+) ")
NUMBERED = re.compile(r"^\d+\. ")
BULLET = re.compile(r"^[-*] ")
TABLE_ROW = re.compile(r"^\|")
TABLE_RULE = re.compile(r"^\|[\s:|-]+\|\s*$")
MATRIX_ROW = re.compile(r"^\|\s*(?P<id>[A-Za-z0-9_-]+)\s*\|.*\|\s*(?P<disp>[A-Za-z]+)\s*\|[^|]*\|\s*$")


def section_kind(heading: str) -> str | None:
    """Classify a slice-summary sub-heading as T, U, P or irrelevant."""
    text = heading.lower()
    if "top findings" in text:
        return "T"
    if "capabilit" in text and ("not use" in text or "off the live" in text):
        return "U"
    if "wrong or incomplete" in text:
        return "P"
    return None


def summary_items(path: Path) -> list[tuple[str, str]]:
    """Return ``(id, first line)`` for every T/U/P item of one record's slice summary."""
    slice_name = path.stem.split("_", 1)[0]
    lines = path.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("## Slice summary"))
    items: list[tuple[str, str]] = []
    counts: Counter[str] = Counter()
    kind: str | None = None
    table_header_seen = False
    for line in lines[start + 1 :]:
        if RECORD_HEADER.match(line):
            break
        heading = HEADING.match(line)
        if heading:
            if len(heading.group(1)) <= 2:
                break
            kind = section_kind(line)
            table_header_seen = False
            continue
        if kind is None:
            continue
        hit = NUMBERED.match(line) or BULLET.match(line)
        if not hit and TABLE_ROW.match(line):
            if TABLE_RULE.match(line):
                continue
            if not table_header_seen:
                table_header_seen = True
                continue
            hit = TABLE_ROW.match(line)
        if hit:
            counts[kind] += 1
            items.append((f"{slice_name}-{kind}{counts[kind]}", line.strip()[:160]))
    return items


def expected_ids() -> dict[str, str]:
    """Every finding id the census defines, mapped to a short source line."""
    found: dict[str, str] = {}
    for record in sorted(OUT.glob("*.md")):
        for line in record.read_text().splitlines():
            match = REGISTER.match(line)
            if match:
                found[match.group(1)] = line[:160]
        for item_id, text in summary_items(record):
            found[item_id] = text
    return found


UNMAPPED_ITEM = re.compile(r"^- `([A-Za-z0-9_-]+)`", re.MULTILINE)


def unmapped_section(text: str) -> str:
    """The body of the matrix's ``## Unmapped`` section."""
    start = text.index("\n## Unmapped\n")
    end = text.find("\n## ", start + 1)
    return text[start : end if end != -1 else len(text)]


def main() -> int:
    expected = expected_ids()
    if "--dump" in sys.argv:
        for item_id, text in expected.items():
            if not item_id.startswith("F-"):
                print(f"{item_id}\t{text}")
        return 0
    rows: Counter[str] = Counter()
    disposition: dict[str, str] = {}
    problems: list[str] = []
    for line in MATRIX.read_text().splitlines():
        match = MATRIX_ROW.match(line)
        if match is None or match["id"] in {"id"} or set(match["id"]) <= {"-"}:
            continue
        item_id, disp = match["id"], match["disp"]
        rows[item_id] += 1
        disposition[item_id] = disp
        if disp not in DISPOSITIONS:
            problems.append(f"bad disposition {disp!r} for {item_id}")
        if item_id not in expected:
            problems.append(f"unknown id in matrix: {item_id}")
    for item_id in expected:
        if rows[item_id] == 0:
            problems.append(f"missing: {item_id}")
        elif rows[item_id] > 1:
            problems.append(f"duplicate ({rows[item_id]}x): {item_id}")
    unmapped = sorted(i for i, d in disposition.items() if d == "UNMAPPED")
    problems.extend(f"UNMAPPED: {i}" for i in unmapped)
    listed = set(UNMAPPED_ITEM.findall(unmapped_section(MATRIX.read_text())))
    problems.extend(f"UNMAPPED row not listed in ## Unmapped: {i}" for i in sorted(set(unmapped) - listed))
    problems.extend(f"listed in ## Unmapped but row is not UNMAPPED: {i}" for i in sorted(listed - set(unmapped)))
    for problem in problems:
        print(problem)
    tally = Counter(disposition[i] for i in expected if i in disposition)
    register = sum(1 for i in expected if i.startswith("F-"))
    print(
        f"ids {len(expected)} (register {register}, summary {len(expected) - register}); "
        + ", ".join(f"{d} {tally[d]}" for d in sorted(DISPOSITIONS, key=str.lower))
        + f"; {'FAIL' if problems else 'OK'}"
    )
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
