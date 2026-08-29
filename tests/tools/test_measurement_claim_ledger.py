from __future__ import annotations

import datetime
import hashlib
import importlib
import importlib.util
import json
import pathlib
import subprocess
from typing import Any

import pytest


PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[2]
LEDGER = PROJECT_ROOT / "docs" / "measurement_claims.md"
FROZEN_COMMIT = "eec665c1df1cd8d1e98dd9dd1001b5984e17a703"
FROZEN_SCAN_SHA256 = "b78ac690bceda37fdebcdcdeb9731a8b213e7d946a897be45181cb95b1100f66"

OPENING_BLOCK = (
    "# Measurement-claim ledger\n\n> **status: in audit — 0/14 Task-5 extractor rows dispositioned**\n\n- Opened (UTC): `2026-08-28`\n- Programme: coherence Wave 1 backlog A\n\n"
)
SCOPE_BLOCK = (
    "This ledger freezes the Task-5 Python-comment regex population. It is not a\n"
    "claim of repository-wide invariant-I2 closure. Tasks 6 and 7 may change only\n"
    "the `disposition` and `evidence / change` cells; all identity columns and\n"
    "frozen matched lines remain immutable.\n"
)
EXTRACTOR_BLOCK = (
    "```bash\n"
    'grep -rnE "# .*(MEASURED|[Mm]easured (on|at|against)|measures [0-9]|[0-9]+(\\.[0-9]+)?[x×] (faster|slower))" \\\n'
    "  src/compas_cgal benchmarks --include='*.py' | grep -v superseded\n"
    "```\n"
)
ORDER_BLOCK = "Canonical order is `src/compas_cgal` before `benchmarks`, then `LC_ALL=C`\nrelative-path order, numeric source line, and raw hit as the final tie-breaker.\n"
DRIFT_BLOCK = (
    "Execution scan (2026-08-28 UTC): exact match — 14/14 frozen hits present;\n"
    "0 changed, 0 missing, 0 new. Normalized stream SHA-256:\n"
    "`b78ac690bceda37fdebcdcdeb9731a8b213e7d946a897be45181cb95b1100f66`.\n"
    "This scan is drift evidence only and does not redefine the frozen population.\n"
)
DISPOSITION_BLOCK = (
    "- `pending` — final adjudication has not happened.\n"
    "- `re-earned` — an authenticated rerun confirmed every material assertion and records the exact command, configuration, artifact, result digest, and full input commit.\n"
    "- `corrected` — the assertion was wrong or incomplete; authenticated evidence and the source correction commit are recorded.\n"
    "- `historical` — the original configuration cannot be reconstructed; the source is explicitly labelled and names the missing identity or configuration.\n"
    "- `deleted` — the assertion was removed while its frozen identity remains here.\n"
    "- `not-a-claim` — semantic review found a regex false positive and records an explicit rationale.\n\n"
    "A case-level result does not automatically disposition every mapped row. Partial\n"
    "reproduction cannot become `re-earned`; every material assertion is adjudicated\n"
    "row by row. `reproduced` is not a ledger disposition.\n"
)

EXPECTED_ROWS = (
    ("001", "MC-001", "src/compas_cgal/engagement_radial_toolpath.py:193", "compas_cgal.engagement_radial_toolpath.RADIUS_LADDER_SUBDIVISIONS", "1/6"),
    ("002", "MC-002", "src/compas_cgal/engagement_radial_toolpath.py:194", "compas_cgal.engagement_radial_toolpath.RADIUS_LADDER_SUBDIVISIONS", "2/6"),
    ("003", "MC-003", "src/compas_cgal/engagement_radial_toolpath.py:195", "compas_cgal.engagement_radial_toolpath.RADIUS_LADDER_SUBDIVISIONS", "3/6"),
    ("004", "MC-004", "src/compas_cgal/engagement_radial_toolpath.py:201", "compas_cgal.engagement_radial_toolpath.RADIUS_LADDER_SUBDIVISIONS", "4/6"),
    ("005", "MC-005", "src/compas_cgal/engagement_radial_toolpath.py:207", "compas_cgal.engagement_radial_toolpath.RADIUS_LADDER_SUBDIVISIONS", "5/6"),
    ("006", "MC-006", "src/compas_cgal/engagement_radial_toolpath.py:225", "compas_cgal.engagement_radial_toolpath.RADIUS_LADDER_SUBDIVISIONS", "6/6"),
    ("007", "MC-007", "src/compas_cgal/engagement_radial_toolpath.py:257", "compas_cgal.engagement_radial_toolpath.RADIUS_LADDER_FLOOR_STEPS", "1/1"),
    ("008", "MC-008", "src/compas_cgal/engagement_radial_toolpath.py:284", "compas_cgal.engagement_radial_toolpath.RADIUS_LADDER_REFINEMENT_MARGIN", "1/1"),
    ("009", "MC-009", "src/compas_cgal/engagement_toolpath.py:120", "compas_cgal.engagement_toolpath.LOOP_PROBE_COUNT", "1/2"),
    ("010", "MC-010", "src/compas_cgal/engagement_toolpath.py:141", "compas_cgal.engagement_toolpath.LOOP_PROBE_COUNT", "2/2"),
    ("011", "MC-011", "benchmarks/gate.py:58", "benchmarks.gate.GATE_LARGE_RECT_WIDTH", "1/1"),
    ("012", "MC-012", "benchmarks/gate.py:66", "benchmarks.gate.L_ARM_TOOL_DIAMETERS", "1/1"),
    ("013", "MC-013", "benchmarks/mathsm.py:47", "benchmarks.mathsm.SPACING_SWEEP_TOOL_DIAMETERS", "1/1"),
    ("014", "MC-014", "benchmarks/quality.py:150", "benchmarks.quality.IMMERSION_STEADY_BAND_FRACTION", "1/1"),
)

EXPECTED_FROZEN_TEXT = (
    "# Measured on the 20x12 pocket at a 60 deg cap, station (18.482, 10.482), maximal",
    "# radius 0.5156, coarse step 0.05: rung 6 (radius 0.2156) measures 61.3 deg and",
    "# still cuts, rung 7 (radius 0.1656) measures 5.9 deg and cuts NOTHING -- one step",
    "# WHY THIS COUNT: MEASURED, NOT DERIVED -- the same footing as `LOOP_PROBE_COUNT`,",
    "# Measured on the 20x12 pocket at a 60 deg cap, tool diameter 2.0, over the 244",
    "# AND DO NOT KEEP THIS WHILE DROPPING `_least_bad_rung`. Measured on 6x4 at a",
    "# MEASURED, on the 20x12 pocket at a 60 deg cap with the gate below in place:",
    "# MEASURED, on the 6x4 pocket at a 40 deg cap -- the hard case, a pocket three tool",
    "# is FALSE. Measured on a 20x12 pocket, 2 mm tool, by walking every machining",
    "# WHY THIS COUNT: MEASURED CONVERGENCE, NOT A DERIVATION. There is a geometric",
    "# The pocket the corner defect was measured on. Ten by six tool diameters.",
    "# W > 4r, hence an arm STRICTLY WIDER THAN TWO TOOL DIAMETERS. Measured at",
    "# spacing stops helping -- on the reference pocket 0.025 measures 131.14 degrees",
    "# measured at the same cap, which is also what makes two generators comparable.",
)


def _module() -> Any:
    return importlib.import_module("tools.measurement_claim_ledger")


def _write_page(tmp_path: pathlib.Path, text: str) -> pathlib.Path:
    page = tmp_path / "measurement_claims.md"
    page.write_text(text, encoding="utf-8")
    return page


def _page_text() -> str:
    return LEDGER.read_text(encoding="utf-8")


def _replace_row(text: str, ordinal: int, column: int, value: str) -> str:
    prefix = f"| {ordinal:03d} | MC-{ordinal:03d} |"
    lines = text.splitlines()
    matches = [index for index, line in enumerate(lines) if line.startswith(prefix)]
    assert len(matches) == 1
    index = matches[0]
    cells = [cell.strip() for cell in lines[index][1:-1].split("|")]
    assert len(cells) == 7
    cells[column] = value
    lines[index] = "| " + " | ".join(cells) + " |"
    return "\n".join(lines) + "\n"


def _terminalize(text: str, count: int, *, complete: bool = False) -> str:
    module = _module()
    for ordinal in range(1, count + 1):
        text = _replace_row(text, ordinal, 5, "historical")
        text = _replace_row(text, ordinal, 6, f"reviewed=MC-{ordinal:03d}")
    wanted = module.COMPLETE_STATUS if complete else module.IN_AUDIT_STATUS.format(done=count)
    return text.replace(module.IN_AUDIT_STATUS.format(done=0), wanted)


def _immutable_digest(rows: tuple[dict[str, str], ...]) -> str:
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


def _git_show(path: str) -> str:
    completed = subprocess.run(
        ["git", "-C", str(PROJECT_ROOT), "show", f"{FROZEN_COMMIT}:{path}"],
        check=True,
        capture_output=True,
    )
    return completed.stdout.decode("utf-8")


SOURCE_COMMIT = "a" * 40
INPUT_SHA256 = "b" * 64
RESULT_SHA256 = "c" * 64
PAYLOAD_SHA256 = "d" * 64
STARTED = datetime.datetime(2026, 8, 29, 9, 15, 30, 123456, tzinfo=datetime.timezone.utc)
ARTIFACT_DIRECTORY = pathlib.PurePosixPath("benchmarks/measurement_claim_results/2026-08-29-aaaaaaaaaaaa-generator-bbbbbbbbbbbb")
GENERATOR_CASES = (
    "radial-station",
    "radial-subdivisions",
    "radial-floor",
    "radial-margin",
    "advance-placement",
    "advance-probe-count",
)
CLAIM_CASES = (
    "radial-station",
    "radial-station",
    "radial-station",
    "radial-subdivisions",
    "radial-subdivisions",
    "radial-subdivisions",
    "radial-floor",
    "radial-margin",
    "advance-placement",
    "advance-probe-count",
)


def _claim_payload() -> dict[str, Any]:
    cases = [
        {
            "case": case,
            "config": {"z": index, "a": [index, True]},
            "selection_decision_provenance": {"policy": f"policy-{index}"},
            "continuous_certificate": None,
        }
        for index, case in enumerate(GENERATOR_CASES, start=1)
    ]
    case_by_name = {case["case"]: case for case in cases}
    claims = [
        {
            "claim_id": f"MC-{index:03d}",
            "case": case,
            "disposition": "re-earned" if index == 1 else "corrected",
            "reason": f"reason {index}",
            "selection_decision_provenance": case_by_name[case]["selection_decision_provenance"],
            "evidence": {"value": index + 0.25, "ordinal": index},
        }
        for index, case in enumerate(CLAIM_CASES, start=1)
    ]
    return {
        "schema_version": "measurement-claim-payload/v1",
        "batch": "generator",
        "extraction_commit": FROZEN_COMMIT,
        "source_commit": SOURCE_COMMIT,
        "case_order": list(GENERATOR_CASES),
        "cases": cases,
        "claims": claims,
    }


def _envelope() -> Any:
    from tools.measurement_artifact import ValidatedEnvelope

    return ValidatedEnvelope.build(
        finished=STARTED + datetime.timedelta(minutes=1),
        commit=SOURCE_COMMIT,
        input_sha256=INPUT_SHA256,
        result_sha256=RESULT_SHA256,
        payload_sha256={"generator-claims.json": PAYLOAD_SHA256},
    )


def _trust_payload_validator(monkeypatch: pytest.MonkeyPatch, module: Any) -> None:
    monkeypatch.setattr(module, "validate_generator_payload", lambda payload: payload)


def _task6_rows(payload: dict[str, Any], evidence: dict[str, str]) -> tuple[dict[str, Any], ...]:
    rows: list[dict[str, Any]] = []
    claims = {claim["claim_id"]: claim for claim in payload["claims"]}
    for ordinal in range(1, 15):
        claim_id = f"MC-{ordinal:03d}"
        claim = claims.get(claim_id)
        rows.append(
            {
                "ordinal": f"{ordinal:03d}",
                "claim_id": claim_id,
                "extracted_location": f"path:{ordinal}",
                "stable_anchor": f"anchor.{ordinal}",
                "anchor_match": "1/1",
                "disposition": "pending" if claim is None else claim["disposition"],
                "evidence": "—" if claim is None else evidence[claim_id],
            }
        )
    return tuple(rows)


def test_task6_evidence_renderer_emits_exact_authenticated_one_line_cells(monkeypatch: pytest.MonkeyPatch) -> None:
    module = _module()
    payload = _claim_payload()
    _trust_payload_validator(monkeypatch, module)

    rendered = module.render_ledger_evidence(
        payload,
        _envelope(),
        started=STARTED,
        artifact_directory=ARTIFACT_DIRECTORY,
    )

    assert tuple(rendered) == tuple(f"MC-{index:03d}" for index in range(1, 11))
    assert rendered["MC-001"] == (
        f"artifact={ARTIFACT_DIRECTORY}; claim=MC-001; input={INPUT_SHA256}; result={RESULT_SHA256}; "
        f'source={SOURCE_COMMIT}; case=radial-station; disposition=re-earned; reason=reason 1; config={{"a":[1,true],"z":1}}; '
        'evidence={"ordinal":1,"value":1.25}; selection={"policy":"policy-1"}; continuous_certificate=null'
    )
    assert all("\n" not in cell and "\r" not in cell and "|" not in cell for cell in rendered.values())


def test_task6_ledger_acceptance_is_exactly_ten_of_fourteen(monkeypatch: pytest.MonkeyPatch) -> None:
    module = _module()
    payload = _claim_payload()
    _trust_payload_validator(monkeypatch, module)
    rendered = module.render_ledger_evidence(payload, _envelope(), started=STARTED, artifact_directory=ARTIFACT_DIRECTORY)
    rows = _task6_rows(payload, rendered)

    assert module.validate_ledger_evidence(rows, payload, _envelope(), started=STARTED, artifact_directory=ARTIFACT_DIRECTORY) == rows
    assert sum(row["disposition"] != "pending" for row in rows) == 10
    assert tuple((row["disposition"], row["evidence"]) for row in rows[10:]) == (("pending", "—"),) * 4


@pytest.mark.parametrize("field", ["disposition", "evidence", "pending"])
def test_task6_ledger_rejects_any_non_byte_equal_or_post_task6_row(monkeypatch: pytest.MonkeyPatch, field: str) -> None:
    module = _module()
    payload = _claim_payload()
    _trust_payload_validator(monkeypatch, module)
    rendered = module.render_ledger_evidence(payload, _envelope(), started=STARTED, artifact_directory=ARTIFACT_DIRECTORY)
    rows = list(_task6_rows(payload, rendered))
    target = dict(rows[10] if field == "pending" else rows[0])
    if field == "disposition":
        target["disposition"] = "historical"
    elif field == "evidence":
        target["evidence"] += " changed"
    else:
        target["disposition"] = "historical"
        target["evidence"] = "premature"
    rows[10 if field == "pending" else 0] = target

    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="MC-001|MC-011|byte|pending|disposition|evidence"):
        module.validate_ledger_evidence(tuple(rows), payload, _envelope(), started=STARTED, artifact_directory=ARTIFACT_DIRECTORY)


@pytest.mark.parametrize("unsafe", ["pipe | reason", "two\nlines", "carriage\rreturn"])
def test_task6_renderer_rejects_unsafe_reason_even_after_payload_validation(monkeypatch: pytest.MonkeyPatch, unsafe: str) -> None:
    module = _module()
    payload = _claim_payload()
    payload["claims"][0]["reason"] = unsafe
    _trust_payload_validator(monkeypatch, module)
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="one physical line|unsafe|Markdown"):
        module.render_ledger_evidence(payload, _envelope(), started=STARTED, artifact_directory=ARTIFACT_DIRECTORY)


@pytest.mark.parametrize(
    ("started", "directory"),
    [
        (STARTED.replace(tzinfo=None), ARTIFACT_DIRECTORY),
        (STARTED.astimezone(datetime.timezone(datetime.timedelta(hours=1))), ARTIFACT_DIRECTORY),
        (STARTED + datetime.timedelta(minutes=2), ARTIFACT_DIRECTORY),
        (STARTED, pathlib.PurePosixPath("benchmarks/measurement_claim_results/../unsafe|artifact")),
        (STARTED, pathlib.PurePosixPath("/benchmarks/measurement_claim_results/2026-08-29-aaaaaaaaaaaa-generator-bbbbbbbbbbbb")),
        (STARTED, pathlib.PurePosixPath("benchmarks/measurement_claim_results/2026-08-28-aaaaaaaaaaaa-generator-bbbbbbbbbbbb")),
    ],
)
def test_task6_renderer_rejects_unsafe_started_or_artifact_path(
    monkeypatch: pytest.MonkeyPatch,
    started: datetime.datetime,
    directory: pathlib.PurePosixPath,
) -> None:
    module = _module()
    _trust_payload_validator(monkeypatch, module)
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="started|artifact"):
        module.render_ledger_evidence(_claim_payload(), _envelope(), started=started, artifact_directory=directory)


def test_task6_renderer_cross_binds_payload_source_to_authenticated_commit(monkeypatch: pytest.MonkeyPatch) -> None:
    module = _module()
    payload = _claim_payload()
    payload["source_commit"] = "e" * 40
    _trust_payload_validator(monkeypatch, module)
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="source commit|artifact commit"):
        module.render_ledger_evidence(payload, _envelope(), started=STARTED, artifact_directory=ARTIFACT_DIRECTORY)


def test_task6_renderer_invokes_the_shared_payload_validator() -> None:
    module = _module()
    result_module = importlib.import_module("tools.measurement_claim_result")
    with pytest.raises(result_module.InvalidMeasurementClaimPayloadError):
        module.render_ledger_evidence(_claim_payload(), _envelope(), started=STARTED, artifact_directory=ARTIFACT_DIRECTORY)


def test_measurement_claim_ledger_module_exists() -> None:
    assert importlib.util.find_spec("tools.measurement_claim_ledger") is not None


def test_real_ledger_has_exact_initial_contract() -> None:
    module = _module()
    text = _page_text()
    rows = module.validate_ledger_structure(LEDGER)
    assert text.startswith(OPENING_BLOCK)
    assert SCOPE_BLOCK in text
    assert EXTRACTOR_BLOCK in text
    assert ORDER_BLOCK in text
    assert DRIFT_BLOCK in text
    assert DISPOSITION_BLOCK in text
    assert module.FROZEN_SOURCE_COMMIT == FROZEN_COMMIT
    assert tuple((row["ordinal"], row["claim_id"], row["extracted_location"], row["stable_anchor"], row["anchor_match"]) for row in rows) == EXPECTED_ROWS
    assert tuple(row["disposition"] for row in rows) == ("pending",) * 14
    assert tuple(row["evidence"] for row in rows) == ("—",) * 14
    assert _immutable_digest(rows) == module.IMMUTABLE_COLUMNS_SHA256


def test_frozen_rows_match_the_full_git_source_object() -> None:
    frozen_lines = {location: _git_show(location.rsplit(":", 1)[0]).splitlines()[int(location.rsplit(":", 1)[1]) - 1] for _, _, location, _, _ in EXPECTED_ROWS}
    assert tuple(frozen_lines[location] for _, _, location, _, _ in EXPECTED_ROWS) == EXPECTED_FROZEN_TEXT
    normalized = "".join(f"{location}:{frozen_lines[location]}\n" for _, _, location, _, _ in EXPECTED_ROWS).encode("utf-8")
    assert hashlib.sha256(normalized).hexdigest() == FROZEN_SCAN_SHA256


@pytest.mark.parametrize("field", ["status", "metadata"])
def test_opening_status_and_metadata_cannot_be_displaced(tmp_path: pathlib.Path, field: str) -> None:
    module = _module()
    text = _page_text()
    status = "> **status: in audit — 0/14 Task-5 extractor rows dispositioned**\n\n"
    metadata = "- Opened (UTC): `2026-08-28`\n- Programme: coherence Wave 1 backlog A\n\n"
    if field == "status":
        text = text.replace(status, "", 1).replace("## Scope and limits\n", "## Scope and limits\n\n" + status, 1)
    else:
        text = text.replace(metadata, "", 1).replace("## Scope and limits\n", "## Scope and limits\n\n" + metadata, 1)
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="opening"):
        module.validate_ledger_structure(_write_page(tmp_path, text))


@pytest.mark.parametrize("damage", ["missing", "moved"])
def test_scope_and_i2_limitation_must_remain_in_scope_section(tmp_path: pathlib.Path, damage: str) -> None:
    module = _module()
    text = _page_text().replace(SCOPE_BLOCK, "", 1)
    if damage == "moved":
        text = text.replace("## Frozen extraction\n", "## Frozen extraction\n\n" + SCOPE_BLOCK, 1)
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="scope|I2"):
        module.validate_ledger_structure(_write_page(tmp_path, text))


@pytest.mark.parametrize(
    ("block", "message"),
    [(EXTRACTOR_BLOCK, "extractor"), (ORDER_BLOCK, "order"), (DRIFT_BLOCK, "drift")],
)
@pytest.mark.parametrize("damage", ["missing", "moved"])
def test_extraction_contract_must_remain_in_frozen_extraction_section(tmp_path: pathlib.Path, block: str, message: str, damage: str) -> None:
    module = _module()
    text = _page_text().replace(block, "", 1)
    if damage == "moved":
        text = text.replace("## Scope and limits\n", "## Scope and limits\n\n" + block, 1)
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match=message):
        module.validate_ledger_structure(_write_page(tmp_path, text))


@pytest.mark.parametrize("damage", ["missing", "moved", "altered"])
def test_disposition_semantics_are_exact_and_section_scoped(tmp_path: pathlib.Path, damage: str) -> None:
    module = _module()
    text = _page_text().replace(DISPOSITION_BLOCK, "", 1)
    if damage == "moved":
        text = text.replace("## Scope and limits\n", "## Scope and limits\n\n" + DISPOSITION_BLOCK, 1)
    elif damage == "altered":
        weakened = DISPOSITION_BLOCK.replace("full input commit", "input commit")
        text = text.replace("## Ledger\n", weakened + "\n## Ledger\n", 1)
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="disposition semantics"):
        module.validate_ledger_structure(_write_page(tmp_path, text))


def test_frozen_block_cannot_move_outside_frozen_section(tmp_path: pathlib.Path) -> None:
    module = _module()
    text = _page_text()
    block = "### MC-001\n\n```text\n" + EXPECTED_FROZEN_TEXT[0] + "\n```\n\n"
    text = text.replace(block, "", 1).replace("## Frozen matched lines\n", block + "## Frozen matched lines\n", 1)
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="frozen|section|heading"):
        module.validate_ledger_structure(_write_page(tmp_path, text))


def test_extra_claim_heading_is_rejected(tmp_path: pathlib.Path) -> None:
    module = _module()
    text = _page_text().replace("## Frozen matched lines\n", "## Frozen matched lines\n\n### MC-999\n", 1)
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="claim heading|frozen"):
        module.validate_ledger_structure(_write_page(tmp_path, text))


@pytest.mark.parametrize(
    ("column", "replacement"),
    [
        (0, "999"),
        (1, "MC-999"),
        (2, "`benchmarks/gate.py:999`"),
        (3, "`benchmarks.gate.WRONG` / leading comment"),
        (4, "9/9"),
    ],
)
def test_each_immutable_column_mutation_is_rejected(tmp_path: pathlib.Path, column: int, replacement: str) -> None:
    module = _module()
    page = _write_page(tmp_path, _replace_row(_page_text(), 1, column, replacement))
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="row|immutable"):
        module.validate_ledger_structure(page)


@pytest.mark.parametrize("damage", ["missing", "duplicate", "out-of-order"])
def test_missing_duplicate_and_out_of_order_rows_are_rejected(tmp_path: pathlib.Path, damage: str) -> None:
    module = _module()
    text = _page_text()
    row_1 = next(line for line in text.splitlines() if line.startswith("| 001 | MC-001 |"))
    row_2 = next(line for line in text.splitlines() if line.startswith("| 002 | MC-002 |"))
    if damage == "missing":
        text = text.replace(f"{row_1}\n", "", 1)
    elif damage == "duplicate":
        text = text.replace(f"{row_2}\n", f"{row_1}\n", 1)
    else:
        text = text.replace(f"{row_1}\n{row_2}\n", f"{row_2}\n{row_1}\n", 1)
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="row|ordinal|duplicate"):
        module.validate_ledger_structure(_write_page(tmp_path, text))


@pytest.mark.parametrize("damage", ["missing", "duplicate", "altered"])
def test_missing_duplicate_and_altered_frozen_blocks_are_rejected(tmp_path: pathlib.Path, damage: str) -> None:
    module = _module()
    text = _page_text()
    block = "### MC-001\n\n```text\n" + EXPECTED_FROZEN_TEXT[0] + "\n```\n"
    if damage == "missing":
        text = text.replace(block, "", 1)
    elif damage == "duplicate":
        text = text.replace(block, block + "\n" + block, 1)
    else:
        text = text.replace(EXPECTED_FROZEN_TEXT[0], EXPECTED_FROZEN_TEXT[0] + " changed", 1)
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="frozen|MC-001"):
        module.validate_ledger_structure(_write_page(tmp_path, text))


def test_all_six_dispositions_are_accepted_with_matching_status(tmp_path: pathlib.Path) -> None:
    module = _module()
    text = _page_text()
    for ordinal, disposition in enumerate(module.DISPOSITIONS, start=1):
        text = _replace_row(text, ordinal, 5, disposition)
        if disposition != "pending":
            text = _replace_row(text, ordinal, 6, f"reviewed=MC-{ordinal:03d}")
    text = text.replace(module.IN_AUDIT_STATUS.format(done=0), module.IN_AUDIT_STATUS.format(done=5))
    rows = module.validate_ledger_structure(_write_page(tmp_path, text))
    assert tuple(row["disposition"] for row in rows[:6]) == module.DISPOSITIONS


@pytest.mark.parametrize("disposition", ["reproduced", "unknown", ""])
def test_unknown_or_empty_disposition_is_rejected(tmp_path: pathlib.Path, disposition: str) -> None:
    module = _module()
    page = _write_page(tmp_path, _replace_row(_page_text(), 1, 5, disposition))
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="disposition"):
        module.validate_ledger_structure(page)


@pytest.mark.parametrize("done", [0, 10])
def test_in_audit_status_is_exact_for_observed_terminal_count(tmp_path: pathlib.Path, done: int) -> None:
    module = _module()
    text = _page_text() if done == 0 else _terminalize(_page_text(), done)
    rows = module.validate_ledger_structure(_write_page(tmp_path, text))
    assert sum(row["disposition"] != "pending" for row in rows) == done


def test_complete_status_is_required_at_fourteen_of_fourteen(tmp_path: pathlib.Path) -> None:
    module = _module()
    text = _terminalize(_page_text(), 14, complete=True)
    assert len(module.validate_ledger_structure(_write_page(tmp_path, text))) == 14
    wrong = text.replace(module.COMPLETE_STATUS, module.IN_AUDIT_STATUS.format(done=14))
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="status"):
        module.validate_ledger_structure(_write_page(tmp_path, wrong))


@pytest.mark.parametrize("suffix", [" extra", "; opened 2026-08-28", "\n> **status: in audit — 0/14 Task-5 extractor rows dispositioned**"])
def test_status_prefix_suffix_and_duplicate_lines_are_rejected(tmp_path: pathlib.Path, suffix: str) -> None:
    module = _module()
    exact = module.IN_AUDIT_STATUS.format(done=0)
    page = _write_page(tmp_path, _page_text().replace(exact, exact + suffix, 1))
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="status"):
        module.validate_ledger_structure(page)


def test_status_count_disagreement_is_rejected(tmp_path: pathlib.Path) -> None:
    module = _module()
    text = _replace_row(_page_text(), 1, 5, "historical")
    text = _replace_row(text, 1, 6, "reviewed=MC-001")
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="status"):
        module.validate_ledger_structure(_write_page(tmp_path, text))


def test_pending_row_cannot_carry_evidence(tmp_path: pathlib.Path) -> None:
    module = _module()
    page = _write_page(tmp_path, _replace_row(_page_text(), 1, 6, "unsupported evidence"))
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="pending|evidence"):
        module.validate_ledger_structure(page)


@pytest.mark.parametrize(
    ("old", "new", "message"),
    [
        ("## Ledger", "## Claim table", "heading"),
        ("| --- | --- | --- | --- | --- | --- | --- |", "| --- | --- |", "separator"),
        ("| 001 | MC-001 |", "| 001 | extra | MC-001 |", "column"),
    ],
)
def test_malformed_heading_separator_and_column_count_raise_named_error(tmp_path: pathlib.Path, old: str, new: str, message: str) -> None:
    module = _module()
    page = _write_page(tmp_path, _page_text().replace(old, new, 1))
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match=message):
        module.validate_ledger_structure(page)


def test_missing_frozen_blobs_raise_named_error(tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch) -> None:
    module = _module()
    monkeypatch.setattr(module, "FROZEN_SOURCE_COMMIT", "0" * 40)
    page = _write_page(tmp_path, _page_text())
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="Git|frozen"):
        module.validate_ledger_structure(page)
