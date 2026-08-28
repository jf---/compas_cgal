"""Validate the frozen measurement-claim ledger identity and state."""

from __future__ import annotations

import ast
import hashlib
import json
import pathlib
import re
import subprocess
from typing import Literal
from typing import Optional
from typing import Sequence
from typing import TypedDict
from typing import cast

Disposition = Literal["pending", "re-earned", "corrected", "historical", "deleted", "not-a-claim"]


class LedgerRow(TypedDict):
    """One parsed measurement-claim ledger row."""

    ordinal: str
    claim_id: str
    extracted_location: str
    stable_anchor: str
    anchor_match: str
    disposition: Disposition
    evidence: str


class _FrozenHit(TypedDict):
    path: str
    line: int
    text: str
    stable_anchor: str
    anchor_match: str


DISPOSITIONS = ("pending", "re-earned", "corrected", "historical", "deleted", "not-a-claim")
IN_AUDIT_STATUS = "> **status: in audit — {done}/14 Task-5 extractor rows dispositioned**"
COMPLETE_STATUS = "> **status: complete — 14/14 Task-5 extractor rows dispositioned**"
FROZEN_SOURCE_COMMIT = "eec665c1df1cd8d1e98dd9dd1001b5984e17a703"
IMMUTABLE_COLUMNS_SHA256 = "e0761be0cce4f4bb5dfcc4335ec3b54ff5751c7e1454870c2273f2147250883e"

_ROW_COUNT = 14
_FILE_COUNT = 5
_LEDGER_HEADER = "| ordinal | claim ID | extracted location | stable anchor | anchor match | disposition | evidence / change |"
_LEDGER_SEPARATOR = "| --- | --- | --- | --- | --- | --- | --- |"
_SECTION_HEADINGS = (
    "## Scope and limits",
    "## Frozen extraction",
    "## Disposition semantics",
    "## Ledger",
    "## Frozen matched lines",
)
_EXTRACTION_PATTERN = re.compile(r"# .*(MEASURED|[Mm]easured (on|at|against)|measures [0-9]|[0-9]+(\.[0-9]+)?[x×] (faster|slower))")


class InvalidMeasurementClaimLedgerError(RuntimeError):
    """The measurement-claim ledger violates its frozen structural contract."""


def _git(repository: pathlib.Path, *arguments: str) -> bytes:
    command = ["git", "-C", str(repository), *arguments]
    try:
        completed = subprocess.run(command, check=False, capture_output=True)
    except OSError as exc:
        raise InvalidMeasurementClaimLedgerError(f"Git invocation failed in frozen-ledger repository {repository}: {' '.join(arguments)}") from exc
    if completed.returncode != 0:
        detail = completed.stderr.decode("utf-8", errors="replace").strip()
        raise InvalidMeasurementClaimLedgerError(f"Git could not read frozen ledger source {FROZEN_SOURCE_COMMIT}: {' '.join(arguments)}: {detail}")
    return completed.stdout


def _repository_root() -> pathlib.Path:
    start = pathlib.Path(__file__).resolve().parent
    raw = _git(start, "rev-parse", "--show-toplevel")
    try:
        root = pathlib.Path(raw.decode("utf-8").strip()).resolve()
    except UnicodeDecodeError as exc:
        raise InvalidMeasurementClaimLedgerError(f"Git repository path is not UTF-8: {start}") from exc
    if not root.is_dir():
        raise InvalidMeasurementClaimLedgerError(f"Git repository path is not a directory: {root}")
    return root


def _decode_frozen_blob(repository: pathlib.Path, path: str) -> str:
    raw = _git(repository, "show", f"{FROZEN_SOURCE_COMMIT}:{path}")
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise InvalidMeasurementClaimLedgerError(f"frozen source blob is not UTF-8: {path}") from exc


def _frozen_python_paths(repository: pathlib.Path) -> tuple[str, ...]:
    raw = _git(
        repository,
        "ls-tree",
        "-r",
        "--name-only",
        FROZEN_SOURCE_COMMIT,
        "--",
        "src/compas_cgal",
        "benchmarks",
    )
    try:
        paths = raw.decode("utf-8").splitlines()
    except UnicodeDecodeError as exc:
        raise InvalidMeasurementClaimLedgerError("frozen Git path list is not UTF-8") from exc
    selected = [path for path in paths if path.endswith(".py")]
    return tuple(sorted(selected, key=lambda path: (0 if path.startswith("src/compas_cgal/") else 1, path)))


def _assignment_name(node: ast.stmt) -> Optional[str]:
    target: Optional[ast.expr] = None
    if isinstance(node, ast.Assign) and len(node.targets) == 1:
        target = node.targets[0]
    elif isinstance(node, ast.AnnAssign):
        target = node.target
    if isinstance(target, ast.Name):
        return target.id
    return None


def _module_name(path: str) -> str:
    without_suffix = path[:-3]
    if without_suffix.startswith("src/"):
        without_suffix = without_suffix[len("src/") :]
    return without_suffix.replace("/", ".")


def _anchor_for_hit(path: str, source: str, line: int) -> str:
    try:
        tree = ast.parse(source, filename=path)
    except SyntaxError as exc:
        raise InvalidMeasurementClaimLedgerError(f"frozen source cannot be parsed for stable anchors: {path}") from exc
    assignments = [(node.lineno, name) for node in tree.body if (name := _assignment_name(node)) is not None and node.lineno > line]
    if not assignments:
        raise InvalidMeasurementClaimLedgerError(f"frozen hit has no following module assignment anchor: {path}:{line}")
    assignment_line, name = min(assignments)
    intervening = source.splitlines()[line : assignment_line - 1]
    if any(text.strip() and not text.lstrip().startswith("#") for text in intervening):
        raise InvalidMeasurementClaimLedgerError(f"frozen hit is not a leading comment for {name}: {path}:{line}")
    return f"{_module_name(path)}.{name}"


def _frozen_hits(repository: pathlib.Path) -> tuple[_FrozenHit, ...]:
    hits: list[_FrozenHit] = []
    for path in _frozen_python_paths(repository):
        source = _decode_frozen_blob(repository, path)
        for line_number, text in enumerate(source.splitlines(), start=1):
            if "superseded" in text or _EXTRACTION_PATTERN.search(text) is None:
                continue
            hits.append(
                {
                    "path": path,
                    "line": line_number,
                    "text": text,
                    "stable_anchor": _anchor_for_hit(path, source, line_number),
                    "anchor_match": "",
                }
            )
    hits.sort(
        key=lambda hit: (
            0 if hit["path"].startswith("src/compas_cgal/") else 1,
            hit["path"],
            hit["line"],
            hit["text"],
        )
    )
    if len(hits) != _ROW_COUNT or len({hit["path"] for hit in hits}) != _FILE_COUNT:
        raise InvalidMeasurementClaimLedgerError(
            f"frozen extraction must produce {_ROW_COUNT} hits across {_FILE_COUNT} files, got {len(hits)} hits across {len({hit['path'] for hit in hits})} files"
        )
    totals = {anchor: sum(hit["stable_anchor"] == anchor for hit in hits) for anchor in {hit["stable_anchor"] for hit in hits}}
    seen: dict[str, int] = {}
    for hit in hits:
        anchor = hit["stable_anchor"]
        seen[anchor] = seen.get(anchor, 0) + 1
        hit["anchor_match"] = f"{seen[anchor]}/{totals[anchor]}"
    return tuple(hits)


def _read_ledger(ledger: pathlib.Path) -> str:
    try:
        return ledger.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise InvalidMeasurementClaimLedgerError(f"measurement claim ledger is not readable UTF-8: {ledger}") from exc


def _unique_line_index(lines: Sequence[str], value: str, *, field: str) -> int:
    indices = [index for index, line in enumerate(lines) if line == value]
    if len(indices) != 1:
        raise InvalidMeasurementClaimLedgerError(f"ledger {field} must occur exactly once: {value}")
    return indices[0]


def _validate_page_frame(text: str) -> tuple[str, ...]:
    lines = tuple(text.splitlines())
    if not lines or lines[0] != "# Measurement-claim ledger":
        raise InvalidMeasurementClaimLedgerError("ledger title heading is missing or malformed")
    section_indices = [_unique_line_index(lines, heading, field="section heading") for heading in _SECTION_HEADINGS]
    if section_indices != sorted(section_indices):
        raise InvalidMeasurementClaimLedgerError("ledger section headings are out of order")
    required_fragments = (
        f"Extraction source commit: `{FROZEN_SOURCE_COMMIT}`",
        "Population: 14 line hits across 5 files",
        "- Opened (UTC): `2026-08-28`",
        "- Programme: coherence Wave 1 backlog A",
    )
    for fragment in required_fragments:
        if lines.count(fragment) != 1:
            raise InvalidMeasurementClaimLedgerError(f"ledger frozen metadata is missing or duplicated: {fragment}")
    return lines


def _unwrap_code(value: str, *, row: int, field: str) -> str:
    if len(value) < 3 or not value.startswith("`") or not value.endswith("`") or "`" in value[1:-1]:
        raise InvalidMeasurementClaimLedgerError(f"ledger row {row:03d} {field} must be one exact inline-code span")
    return value[1:-1]


def _parse_table(lines: Sequence[str]) -> tuple[LedgerRow, ...]:
    header_index = _unique_line_index(lines, _LEDGER_HEADER, field="ledger table heading")
    if header_index + 1 >= len(lines) or lines[header_index + 1] != _LEDGER_SEPARATOR:
        raise InvalidMeasurementClaimLedgerError("ledger table separator is missing or malformed")
    body: list[str] = []
    for line in lines[header_index + 2 :]:
        if not line.startswith("|"):
            break
        body.append(line)
    if len(body) != _ROW_COUNT:
        raise InvalidMeasurementClaimLedgerError(f"ledger table must contain {_ROW_COUNT} body rows, got {len(body)}")
    rows: list[LedgerRow] = []
    for expected, line in enumerate(body, start=1):
        if not line.endswith("|"):
            raise InvalidMeasurementClaimLedgerError(f"ledger row {expected:03d} is not pipe-delimited")
        cells = [cell.strip() for cell in line[1:-1].split("|")]
        if len(cells) != 7:
            raise InvalidMeasurementClaimLedgerError(f"ledger row {expected:03d} column count must be 7, got {len(cells)}")
        ordinal, claim_id, location_cell, anchor_cell, anchor_match, disposition_text, evidence = cells
        expected_ordinal = f"{expected:03d}"
        expected_id = f"MC-{expected:03d}"
        if ordinal != expected_ordinal or claim_id != expected_id:
            raise InvalidMeasurementClaimLedgerError(f"ledger row {expected:03d} must use ordinal {expected_ordinal} and claim ID {expected_id}")
        location = _unwrap_code(location_cell, row=expected, field="extracted location")
        anchor_suffix = " / leading comment"
        if not anchor_cell.endswith(anchor_suffix):
            raise InvalidMeasurementClaimLedgerError(f"ledger row {expected:03d} stable anchor suffix is malformed")
        anchor = _unwrap_code(anchor_cell[: -len(anchor_suffix)], row=expected, field="stable anchor")
        if disposition_text not in DISPOSITIONS:
            raise InvalidMeasurementClaimLedgerError(f"ledger row {expected:03d} disposition is empty or unknown: {disposition_text!r}")
        disposition = cast(Disposition, disposition_text)
        if disposition == "pending" and evidence != "—":
            raise InvalidMeasurementClaimLedgerError(f"ledger row {expected:03d} pending disposition must have em-dash evidence")
        if disposition != "pending" and (not evidence or evidence == "—"):
            raise InvalidMeasurementClaimLedgerError(f"ledger row {expected:03d} terminal disposition requires evidence")
        rows.append(
            {
                "ordinal": ordinal,
                "claim_id": claim_id,
                "extracted_location": location,
                "stable_anchor": anchor,
                "anchor_match": anchor_match,
                "disposition": disposition,
                "evidence": evidence,
            }
        )
    return tuple(rows)


def _immutable_digest(rows: Sequence[LedgerRow]) -> str:
    immutable = [
        {
            "ordinal": row["ordinal"],
            "claim_id": row["claim_id"],
            "extracted_location": row["extracted_location"],
            "stable_anchor": row["stable_anchor"],
            "anchor_match": row["anchor_match"],
        }
        for row in rows
    ]
    raw = json.dumps(immutable, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _validate_rows(rows: Sequence[LedgerRow], hits: Sequence[_FrozenHit]) -> None:
    if len({row["claim_id"] for row in rows}) != _ROW_COUNT:
        raise InvalidMeasurementClaimLedgerError("ledger contains a duplicate claim ID")
    if len({row["extracted_location"] for row in rows}) != _ROW_COUNT:
        raise InvalidMeasurementClaimLedgerError("ledger contains a duplicate extracted location")
    for index, (row, hit) in enumerate(zip(rows, hits), start=1):
        expected = (
            f"{index:03d}",
            f"MC-{index:03d}",
            f"{hit['path']}:{hit['line']}",
            hit["stable_anchor"],
            hit["anchor_match"],
        )
        actual = (
            row["ordinal"],
            row["claim_id"],
            row["extracted_location"],
            row["stable_anchor"],
            row["anchor_match"],
        )
        if actual != expected:
            raise InvalidMeasurementClaimLedgerError(f"ledger row {index:03d} immutable identity differs from frozen source: expected {expected!r}, got {actual!r}")
    digest = _immutable_digest(rows)
    if digest != IMMUTABLE_COLUMNS_SHA256:
        raise InvalidMeasurementClaimLedgerError(f"ledger immutable-column digest differs: expected {IMMUTABLE_COLUMNS_SHA256}, got {digest}")


def _validate_frozen_blocks(text: str, hits: Sequence[_FrozenHit]) -> None:
    pattern = re.compile(r"(?m)^### (MC-[0-9]{3})\n\n```text\n([^\n]*)\n```\n")
    blocks = pattern.findall(text)
    if len(blocks) != _ROW_COUNT:
        raise InvalidMeasurementClaimLedgerError(f"ledger must contain {_ROW_COUNT} exact frozen text blocks, got {len(blocks)}")
    for index, ((claim_id, frozen_text), hit) in enumerate(zip(blocks, hits), start=1):
        expected_id = f"MC-{index:03d}"
        if claim_id != expected_id or frozen_text != hit["text"]:
            raise InvalidMeasurementClaimLedgerError(f"ledger frozen block {expected_id} differs from {FROZEN_SOURCE_COMMIT}: got heading {claim_id!r} and text {frozen_text!r}")


def _validate_status(lines: Sequence[str], rows: Sequence[LedgerRow]) -> None:
    status_lines = [line for line in lines if line.startswith("> **status:")]
    done = sum(row["disposition"] != "pending" for row in rows)
    expected = COMPLETE_STATUS if done == _ROW_COUNT else IN_AUDIT_STATUS.format(done=done)
    if status_lines != [expected]:
        raise InvalidMeasurementClaimLedgerError(f"ledger status must equal observed disposition count byte-for-byte: expected {expected!r}, got {status_lines!r}")


def validate_ledger_structure(ledger: pathlib.Path) -> tuple[LedgerRow, ...]:
    """Validate the sole Markdown claim ledger against its frozen Git source.

    Args:
        ledger: Markdown ledger to parse.

    Returns:
        The ordered, validated ledger rows.

    Raises:
        InvalidMeasurementClaimLedgerError: The page, Git source, or frozen
            identity violates the contract.
    """
    text = _read_ledger(ledger)
    lines = _validate_page_frame(text)
    rows = _parse_table(lines)
    hits = _frozen_hits(_repository_root())
    _validate_rows(rows, hits)
    _validate_frozen_blocks(text, hits)
    _validate_status(lines, rows)
    return rows
