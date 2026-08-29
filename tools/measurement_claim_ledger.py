"""Validate the frozen measurement-claim ledger identity and state."""

from __future__ import annotations

import ast
import datetime
import hashlib
import json
import pathlib
import re
import subprocess
from typing import Dict
from typing import Literal
from typing import Optional
from typing import Sequence
from typing import TypedDict
from typing import cast

from tools.measurement_artifact import ValidatedEnvelope
from tools.measurement_claim_result import GeneratorCasePayload
from tools.measurement_claim_result import GeneratorClaimPayload
from tools.measurement_claim_result import GeneratorClaimRecord
from tools.measurement_claim_result import ValidatedArtifactDirectory
from tools.measurement_claim_result import ValidatedArtifactStartedUtc
from tools.measurement_claim_result import validate_generator_payload

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
_SCOPE_REQUIRED = (
    "This ledger freezes the Task-5 Python-comment regex population. It is not a",
    "claim of repository-wide invariant-I2 closure. Tasks 6 and 7 may change only",
    "the `disposition` and `evidence / change` cells; all identity columns and",
    "frozen matched lines remain immutable.",
)
_FROZEN_METADATA_REQUIRED = (
    f"Extraction source commit: `{FROZEN_SOURCE_COMMIT}`",
    "",
    "Population: 14 line hits across 5 files",
)
_EXTRACTOR_REQUIRED = (
    "```bash",
    'grep -rnE "# .*(MEASURED|[Mm]easured (on|at|against)|measures [0-9]|[0-9]+(\\.[0-9]+)?[x×] (faster|slower))" \\',
    "  src/compas_cgal benchmarks --include='*.py' | grep -v superseded",
    "```",
)
_ORDER_REQUIRED = (
    "Canonical order is `src/compas_cgal` before `benchmarks`, then `LC_ALL=C`",
    "relative-path order, numeric source line, and raw hit as the final tie-breaker.",
)
_DRIFT_REQUIRED = (
    "Execution scan (2026-08-28 UTC): exact match — 14/14 frozen hits present;",
    "0 changed, 0 missing, 0 new. Normalized stream SHA-256:",
    "`b78ac690bceda37fdebcdcdeb9731a8b213e7d946a897be45181cb95b1100f66`.",
    "This scan is drift evidence only and does not redefine the frozen population.",
)
_DISPOSITION_SEMANTICS_REQUIRED = (
    "- `pending` — final adjudication has not happened.",
    ("- `re-earned` — an authenticated rerun confirmed every material assertion and records the exact command, configuration, artifact, result digest, and full input commit."),
    ("- `corrected` — the assertion was wrong or incomplete; authenticated evidence and the source correction commit are recorded."),
    ("- `historical` — the original configuration cannot be reconstructed; the source is explicitly labelled and names the missing identity or configuration."),
    "- `deleted` — the assertion was removed while its frozen identity remains here.",
    ("- `not-a-claim` — semantic review found a regex false positive and records an explicit rationale."),
    "",
    "A case-level result does not automatically disposition every mapped row. Partial",
    "reproduction cannot become `re-earned`; every material assertion is adjudicated",
    "row by row. `reproduced` is not a ledger disposition.",
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


def _section(lines: Sequence[str], heading: str, next_heading: Optional[str]) -> tuple[str, ...]:
    start = _unique_line_index(lines, heading, field="section heading") + 1
    end = len(lines) if next_heading is None else _unique_line_index(lines, next_heading, field="section heading")
    return tuple(lines[start:end])


def _block_index(section: Sequence[str], block: Sequence[str], *, field: str) -> int:
    width = len(block)
    indices = [index for index in range(len(section) - width + 1) if tuple(section[index : index + width]) == tuple(block)]
    if len(indices) != 1:
        raise InvalidMeasurementClaimLedgerError(f"ledger {field} must occur exactly once in its owning section")
    return indices[0]


def _validate_page_frame(text: str) -> tuple[str, ...]:
    lines = tuple(text.splitlines())
    opening = (
        "# Measurement-claim ledger",
        "",
        lines[2] if len(lines) > 2 else "",
        "",
        "- Opened (UTC): `2026-08-28`",
        "- Programme: coherence Wave 1 backlog A",
        "",
        "## Scope and limits",
    )
    if len(lines) < len(opening) or tuple(lines[: len(opening)]) != opening or not opening[2].startswith("> **status:"):
        raise InvalidMeasurementClaimLedgerError("ledger opening status and metadata layout is malformed")

    actual_headings = tuple(line for line in lines if line.startswith("## "))
    if actual_headings != _SECTION_HEADINGS:
        raise InvalidMeasurementClaimLedgerError(f"ledger section headings must be exact and ordered: expected {_SECTION_HEADINGS!r}, got {actual_headings!r}")

    scope = _section(lines, "## Scope and limits", "## Frozen extraction")
    _block_index(scope, _SCOPE_REQUIRED, field="scope and I2 limitation")

    extraction = _section(lines, "## Frozen extraction", "## Disposition semantics")
    extraction_blocks = (
        ("frozen extraction metadata", _FROZEN_METADATA_REQUIRED),
        ("canonical extractor", _EXTRACTOR_REQUIRED),
        ("canonical order", _ORDER_REQUIRED),
        ("separate drift report", _DRIFT_REQUIRED),
    )
    positions = [_block_index(extraction, block, field=field) for field, block in extraction_blocks]
    if positions != sorted(positions):
        raise InvalidMeasurementClaimLedgerError("ledger extractor, canonical order, and drift report are out of order")

    semantics = _section(lines, "## Disposition semantics", "## Ledger")
    _block_index(semantics, _DISPOSITION_SEMANTICS_REQUIRED, field="disposition semantics")
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


def _validate_frozen_blocks(lines: Sequence[str], hits: Sequence[_FrozenHit]) -> None:
    expected_headings = tuple(f"### MC-{index:03d}" for index in range(1, _ROW_COUNT + 1))
    actual_headings = tuple(line for line in lines if line.startswith("### MC-"))
    if actual_headings != expected_headings:
        raise InvalidMeasurementClaimLedgerError(
            f"ledger claim headings must be exact and confined to the frozen section: expected {expected_headings!r}, got {actual_headings!r}"
        )
    expected_section: list[str] = [""]
    for index, hit in enumerate(hits, start=1):
        expected_section.extend((f"### MC-{index:03d}", "", "```text", hit["text"], "```"))
        if index != _ROW_COUNT:
            expected_section.append("")
    actual_section = _section(lines, "## Frozen matched lines", None)
    if actual_section != tuple(expected_section):
        raise InvalidMeasurementClaimLedgerError("ledger frozen blocks must be exact and wholly contained by the Frozen matched lines section")


def _validate_status(lines: Sequence[str], rows: Sequence[LedgerRow]) -> None:
    status_lines = [line for line in lines if line.startswith("> **status:")]
    done = sum(row["disposition"] != "pending" for row in rows)
    expected = COMPLETE_STATUS if done == _ROW_COUNT else IN_AUDIT_STATUS.format(done=done)
    if status_lines != [expected]:
        raise InvalidMeasurementClaimLedgerError(f"ledger status must equal observed disposition count byte-for-byte: expected {expected!r}, got {status_lines!r}")


def _artifact_text(
    payload: GeneratorClaimPayload,
    envelope: ValidatedEnvelope,
    *,
    started: ValidatedArtifactStartedUtc,
    artifact_directory: ValidatedArtifactDirectory,
) -> str:
    if type(started) is not datetime.datetime or started.tzinfo is None or started.utcoffset() != datetime.timedelta(0):
        raise InvalidMeasurementClaimLedgerError("Task-6 artifact started value must be one timezone-aware UTC datetime")
    if started > envelope.finished:
        raise InvalidMeasurementClaimLedgerError("Task-6 artifact started value must not follow the authenticated finish")
    if payload["source_commit"] != envelope.commit:
        raise InvalidMeasurementClaimLedgerError("Task-6 payload source commit differs from the authenticated artifact commit")
    expected = pathlib.PurePosixPath(
        "benchmarks",
        "measurement_claim_results",
        f"{started.date().isoformat()}-{str(envelope.commit)[:12]}-generator-{str(envelope.input_sha256)[:12]}",
    )
    if type(artifact_directory) is not pathlib.PurePosixPath or artifact_directory != expected:
        raise InvalidMeasurementClaimLedgerError(f"Task-6 artifact directory must equal the authenticated repository-relative path: {expected}")
    text = artifact_directory.as_posix()
    if any(token in text for token in ("|", "\r", "\n", "\u2028", "\u2029")):
        raise InvalidMeasurementClaimLedgerError("Task-6 artifact directory is not Markdown-table-safe")
    return text


def _canonical_json(value: object, *, field: str) -> str:
    try:
        return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)
    except (TypeError, ValueError) as exc:
        raise InvalidMeasurementClaimLedgerError(f"Task-6 {field} is not canonical finite JSON") from exc


def _render_validated_ledger_evidence(
    payload: GeneratorClaimPayload,
    envelope: ValidatedEnvelope,
    *,
    started: ValidatedArtifactStartedUtc,
    artifact_directory: ValidatedArtifactDirectory,
) -> Dict[str, str]:
    artifact = _artifact_text(payload, envelope, started=started, artifact_directory=artifact_directory)
    cases: Dict[str, GeneratorCasePayload] = {case["case"]: case for case in payload["cases"]}
    rendered: Dict[str, str] = {}
    for claim in payload["claims"]:
        claim_id = claim["claim_id"]
        case = cases[claim["case"]]
        cell = (
            f"artifact={artifact}; claim={claim_id}; input={envelope.input_sha256}; result={envelope.result_sha256}; "
            f"source={payload['source_commit']}; case={claim['case']}; disposition={claim['disposition']}; reason={claim['reason']}; "
            f"config={_canonical_json(case['config'], field=f'{claim_id} config')}; "
            f"evidence={_canonical_json(claim['evidence'], field=f'{claim_id} evidence')}; "
            f"selection={_canonical_json(claim['selection_decision_provenance'], field=f'{claim_id} selection')}; "
            "continuous_certificate=null"
        )
        if any(token in cell for token in ("|", "\r", "\n", "\u2028", "\u2029")):
            raise InvalidMeasurementClaimLedgerError(f"Task-6 {claim_id} evidence must be Markdown-safe and one physical line")
        rendered[claim_id] = cell
    return rendered


def render_ledger_evidence(
    payload: GeneratorClaimPayload,
    envelope: ValidatedEnvelope,
    *,
    started: ValidatedArtifactStartedUtc,
    artifact_directory: ValidatedArtifactDirectory,
) -> Dict[str, str]:
    """Render the ten authenticated generator-claim ledger evidence cells."""
    validated = validate_generator_payload(payload)
    return _render_validated_ledger_evidence(
        validated,
        envelope,
        started=started,
        artifact_directory=artifact_directory,
    )


def validate_ledger_evidence(
    rows: Sequence[LedgerRow],
    payload: GeneratorClaimPayload,
    envelope: ValidatedEnvelope,
    *,
    started: ValidatedArtifactStartedUtc,
    artifact_directory: ValidatedArtifactDirectory,
) -> tuple[LedgerRow, ...]:
    """Require the exact Task-6 ten-of-fourteen ledger acceptance state."""
    if len(rows) != _ROW_COUNT:
        raise InvalidMeasurementClaimLedgerError(f"Task-6 ledger acceptance requires exactly {_ROW_COUNT} rows")
    validated = validate_generator_payload(payload)
    evidence = _render_validated_ledger_evidence(
        validated,
        envelope,
        started=started,
        artifact_directory=artifact_directory,
    )
    claims: Dict[str, GeneratorClaimRecord] = {claim["claim_id"]: claim for claim in validated["claims"]}
    for index, row in enumerate(rows, start=1):
        claim_id = f"MC-{index:03d}"
        if row["ordinal"] != f"{index:03d}" or row["claim_id"] != claim_id:
            raise InvalidMeasurementClaimLedgerError(f"Task-6 ledger row {index:03d} identity is not canonical")
        if index <= 10:
            claim = claims[claim_id]
            if row["disposition"] != claim["disposition"]:
                raise InvalidMeasurementClaimLedgerError(f"Task-6 {claim_id} disposition differs from the authenticated payload")
            if row["evidence"] != evidence[claim_id]:
                raise InvalidMeasurementClaimLedgerError(f"Task-6 {claim_id} evidence is not byte-equal to the authenticated rendering")
        elif row["disposition"] != "pending" or row["evidence"] != "—":
            raise InvalidMeasurementClaimLedgerError(f"Task-6 {claim_id} must remain pending with em-dash evidence")
    return tuple(rows)


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
    _validate_frozen_blocks(lines, hits)
    _validate_status(lines, rows)
    return rows
