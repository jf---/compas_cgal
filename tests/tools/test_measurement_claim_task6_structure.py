from __future__ import annotations

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
RADIAL_SOURCE = "src/compas_cgal/engagement_radial_toolpath.py"
ADVANCE_SOURCE = "src/compas_cgal/engagement_toolpath.py"
SOURCE_FILE_SHA256 = {
    RADIAL_SOURCE: "4a12b0d8a404a7355271eafb10baa864e33d64411a60ca3195ccf85bd1a19af0",
    ADVANCE_SOURCE: "6bc5095fd7853949c4c6f32bcfbe7b2d85ffd038b82e6ce8cf6d259e224117cf",
}
SOURCE_REGION_BASELINES = (
    ("radial-module-which-circle", RADIAL_SOURCE, 70, 76, "1ab67dae4cbc8dea240b7f7a572be932eed99022427e7d773bb2442918ec7ccc"),
    ("radius-ladder-subdivisions", RADIAL_SOURCE, 186, 246, "e3898ee926c10c0dcac29a92d4811aa2f98ef1dc0fa9591b54f5b1e99a769a34"),
    ("radius-ladder-refinement-margin", RADIAL_SOURCE, 264, 305, "ba810f64945c4da67a1e994b2bab241dbba0bbb3d24241c9430f36178b1d0461"),
    ("gentlest-rung-peak", RADIAL_SOURCE, 393, 398, "19357c059c7b7f67247485b1613c8b043077ac5be332e16e4b9c8e0bde89b4b5"),
    ("least-bad-rung-double", RADIAL_SOURCE, 607, 614, "4aeee76e5df523410d1f0adc00d26fd6207d3c2f6b52937f9eceae80d6dfc311"),
    ("regulation-cap-angle", ADVANCE_SOURCE, 248, 254, "5ecd9cc249e9138e4dc9a7339ae7fd6cf33bc76432fab2c6a756c5347f45ca4a"),
    ("measured-peak-reporting", ADVANCE_SOURCE, 580, 587, "39fa70015f93d09c420d337aa80dbca94dbfe2dad4b86b72060c9698c5b96680"),
    ("loop-probe-count", ADVANCE_SOURCE, 109, 161, "a913fa0d7ce1c67304c72c1751e445c702798c6a6a887464f37d93394b6ad2ad"),
    ("radius-ladder-floor-steps", RADIAL_SOURCE, 249, 261, "a3b0954dde586e1abd2d3399c81c05d70e2eb748c1db247f538b945217ea8128"),
)
MANDATORY_SOURCE_REGIONS = tuple(region[0] for region in SOURCE_REGION_BASELINES[:8])

SOURCE_CORRECTIONS = {
    "radial-module-which-circle": b"WHICH circle is selected by the current reporting-driven policy; it is not an exact cap decision.\n",
    "radius-ladder-subdivisions": b"# Corrected subdivision evidence is authenticated by the Task-6 artifact.\n",
    "radius-ladder-refinement-margin": b"# This reporting comparison controls whether the refined scan executes.\n",
    "gentlest-rung-peak": b"    The peak travels with the rung because current forced-radius selection consumes this reported value.\n",
    "least-bad-rung-double": b"    Reported engagement ranks refused radii and therefore determines the forced circle emitted.\n",
    "regulation-cap-angle": b"        cap_angle: Reported cap angle used by the current refinement-control comparison.\n",
    "measured-peak-reporting": b"    This reporting value participates in forced-radius selection and refined-scan control.\n",
    "loop-probe-count": b"# Corrected finite-sweep evidence is authenticated by the Task-6 artifact.\n",
    "radius-ladder-floor-steps": b"# Corrected floor evidence is authenticated by the Task-6 artifact.\n",
}

INITIAL_STATUS = "> **status: in audit — 0/14 Task-5 extractor rows dispositioned**"
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
    return importlib.import_module("tools.measurement_claim_task6_ledger")


def _write_page(tmp_path: pathlib.Path, text: str) -> pathlib.Path:
    page = tmp_path / "measurement_claims.md"
    page.write_text(text, encoding="utf-8")
    return page


def _live_page_text() -> str:
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


def _initial_page_text() -> str:
    text = _live_page_text()
    for ordinal in range(1, 15):
        text = _replace_row(text, ordinal, 5, "pending")
        text = _replace_row(text, ordinal, 6, "—")
    lines = text.splitlines()
    status_indexes = [index for index, line in enumerate(lines) if line.startswith("> **status:")]
    assert len(status_indexes) == 1
    lines[status_indexes[0]] = INITIAL_STATUS
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


def _git_at(repository: pathlib.Path, *arguments: str, stdin: bytes | None = None) -> bytes:
    return subprocess.run(
        ["git", "-C", str(repository), *arguments],
        check=True,
        capture_output=True,
        input=stdin,
    ).stdout


def _commit(repository: pathlib.Path, message: str) -> str:
    _git_at(repository, "add", "-A")
    _git_at(
        repository,
        "-c",
        "user.name=Jelle Feringa",
        "-c",
        "user.email=jelleferinga@gmail.com",
        "commit",
        "-qm",
        message,
    )
    return _git_at(repository, "rev-parse", "HEAD^{commit}").decode("ascii").strip()


def _corrected_source_bytes(path: str, baseline: bytes, *, omitted: frozenset[str], correct_floor: bool) -> bytes:
    lines = baseline.splitlines(keepends=True)
    regions = [region for region in SOURCE_REGION_BASELINES if region[1] == path]
    for name, _, start, end, _ in sorted(regions, key=lambda region: region[2], reverse=True):
        if name in omitted or (name == "radius-ladder-floor-steps" and not correct_floor):
            continue
        lines[start - 1 : end] = [SOURCE_CORRECTIONS[name]]
    return b"".join(lines)


def _source_repository(
    tmp_path: pathlib.Path,
    *,
    omitted: frozenset[str] = frozenset(),
    correct_floor: bool = False,
    baseline_suffix: bytes = b"",
    candidate_mutation: str | None = None,
    extra_path: bool = False,
    executable_path: bool = False,
    object_format: str = "sha1",
) -> tuple[pathlib.Path, str, str]:
    repository = tmp_path / "source-repository"
    repository.mkdir()
    _git_at(repository, "init", "-q", f"--object-format={object_format}")
    baseline: dict[str, bytes] = {}
    for path in SOURCE_FILE_SHA256:
        raw = _git_at(PROJECT_ROOT, "show", f"{FROZEN_COMMIT}:{path}")
        if path == RADIAL_SOURCE:
            raw += baseline_suffix
        baseline[path] = raw
        target = repository / path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(raw)
    parent = _commit(repository, "baseline")

    candidate = {path: _corrected_source_bytes(path, raw, omitted=omitted, correct_floor=correct_floor) for path, raw in baseline.items()}
    if candidate_mutation == "unowned-docstring":
        candidate[RADIAL_SOURCE] = candidate[RADIAL_SOURCE].replace(
            b"READ THAT NUMBER PRECISELY:",
            b"READ THIS NUMBER PRECISELY:",
            1,
        )
    elif candidate_mutation == "executable":
        candidate[RADIAL_SOURCE] = candidate[RADIAL_SOURCE].replace(
            b"FULL_RADIUS_RUNG = 0\n",
            b"FULL_RADIUS_RUNG = 1\n",
            1,
        )
    elif candidate_mutation == "whitespace":
        candidate[RADIAL_SOURCE] = candidate[RADIAL_SOURCE].replace(
            b"FULL_RADIUS_RUNG = 0\n",
            b"FULL_RADIUS_RUNG = 0  \n",
            1,
        )
    for path, raw in candidate.items():
        (repository / path).write_bytes(raw)
    if extra_path:
        (repository / "unexpected.txt").write_text("not source correction\n", encoding="utf-8")
    _git_at(repository, "add", "-A")
    if executable_path:
        _git_at(repository, "update-index", "--chmod=+x", RADIAL_SOURCE)
    _git_at(
        repository,
        "-c",
        "user.name=Jelle Feringa",
        "-c",
        "user.email=jelleferinga@gmail.com",
        "commit",
        "-qm",
        "correction",
    )
    correction = _git_at(repository, "rev-parse", "HEAD^{commit}").decode("ascii").strip()
    return repository, parent, correction


def _merge_commit(repository: pathlib.Path, correction: str, parent: str) -> str:
    tree = _git_at(repository, "show", "-s", "--format=%T", correction).decode("ascii").strip()
    return (
        _git_at(
            repository,
            "-c",
            "user.name=Jelle Feringa",
            "-c",
            "user.email=jelleferinga@gmail.com",
            "commit-tree",
            tree,
            "-p",
            correction,
            "-p",
            parent,
            stdin=b"merge\n",
        )
        .decode("ascii")
        .strip()
    )


def _commit_tree(repository: pathlib.Path, tree_commit: str, parent: str, message: bytes) -> str:
    tree = _git_at(repository, "show", "-s", "--format=%T", tree_commit).decode("ascii").strip()
    return (
        _git_at(
            repository,
            "-c",
            "user.name=Jelle Feringa",
            "-c",
            "user.email=jelleferinga@gmail.com",
            "commit-tree",
            tree,
            "-p",
            parent,
            stdin=message,
        )
        .decode("ascii")
        .strip()
    )


def _merge_with_tree(repository: pathlib.Path, tree_commit: str, first_parent: str, second_parent: str) -> str:
    tree = _git_at(repository, "show", "-s", "--format=%T", tree_commit).decode("ascii").strip()
    return (
        _git_at(
            repository,
            "-c",
            "user.name=Jelle Feringa",
            "-c",
            "user.email=jelleferinga@gmail.com",
            "commit-tree",
            tree,
            "-p",
            first_parent,
            "-p",
            second_parent,
            stdin=b"merge execution\n",
        )
        .decode("ascii")
        .strip()
    )


def _test_source_region(lines: list[bytes], line_start: int, line_end: int) -> dict[str, object]:
    byte_start = sum(len(line) for line in lines[: line_start - 1])
    byte_end = sum(len(line) for line in lines[:line_end])
    return {
        "raw": b"".join(lines[line_start - 1 : line_end]),
        "byte_start": byte_start,
        "byte_end": byte_end,
        "line_start": line_start,
        "line_end": line_end,
    }


def _zero_count_hunk_repository(
    tmp_path: pathlib.Path,
    *,
    operation: str,
    zero_position: int,
    delete_count: int = 1,
) -> tuple[pathlib.Path, str, str, dict[str, dict[str, dict[str, object]]], dict[str, dict[str, dict[str, object]]]]:
    repository = tmp_path / "hunk-repository"
    repository.mkdir()
    _git_at(repository, "init", "-q")
    baseline_lines = [f"line-{line}\n".encode("ascii") for line in range(1, 15)]
    target = repository / RADIAL_SOURCE
    target.parent.mkdir(parents=True)
    target.write_bytes(b"".join(baseline_lines))
    parent = _commit(repository, "baseline")

    candidate_lines = baseline_lines.copy()
    if operation == "insertion":
        candidate_lines.insert(zero_position, b"inserted\n")
    else:
        del candidate_lines[zero_position : zero_position + delete_count]
    target.write_bytes(b"".join(candidate_lines))
    correction = _commit(repository, operation)

    region_name = "radial-module-which-circle"
    line_start = 5
    line_end = 8
    candidate_end = line_end
    if operation == "insertion" and line_start - 1 <= zero_position <= line_end:
        candidate_end += 1
    elif operation == "deletion":
        first_deleted = zero_position + 1
        last_deleted = zero_position + delete_count
        candidate_end -= max(0, min(line_end, last_deleted) - max(line_start, first_deleted) + 1)
    baseline_regions = {RADIAL_SOURCE: {region_name: _test_source_region(baseline_lines, line_start, line_end)}}
    candidate_regions = {RADIAL_SOURCE: {region_name: _test_source_region(candidate_lines, line_start, candidate_end)}}
    return repository, parent, correction, baseline_regions, candidate_regions


def test_frozen_rows_match_the_full_git_source_object() -> None:
    frozen_lines = {location: _git_show(location.rsplit(":", 1)[0]).splitlines()[int(location.rsplit(":", 1)[1]) - 1] for _, _, location, _, _ in EXPECTED_ROWS}
    assert tuple(frozen_lines[location] for _, _, location, _, _ in EXPECTED_ROWS) == EXPECTED_FROZEN_TEXT
    normalized = "".join(f"{location}:{frozen_lines[location]}\n" for _, _, location, _, _ in EXPECTED_ROWS).encode("utf-8")
    assert hashlib.sha256(normalized).hexdigest() == FROZEN_SCAN_SHA256


@pytest.mark.parametrize("field", ["status", "metadata"])
def test_opening_status_and_metadata_cannot_be_displaced(tmp_path: pathlib.Path, field: str) -> None:
    module = _module()
    text = _initial_page_text()
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
    text = _initial_page_text().replace(SCOPE_BLOCK, "", 1)
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
    text = _initial_page_text().replace(block, "", 1)
    if damage == "moved":
        text = text.replace("## Scope and limits\n", "## Scope and limits\n\n" + block, 1)
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match=message):
        module.validate_ledger_structure(_write_page(tmp_path, text))


@pytest.mark.parametrize("damage", ["missing", "moved", "altered"])
def test_disposition_semantics_are_exact_and_section_scoped(tmp_path: pathlib.Path, damage: str) -> None:
    module = _module()
    text = _initial_page_text().replace(DISPOSITION_BLOCK, "", 1)
    if damage == "moved":
        text = text.replace("## Scope and limits\n", "## Scope and limits\n\n" + DISPOSITION_BLOCK, 1)
    elif damage == "altered":
        weakened = DISPOSITION_BLOCK.replace("full input commit", "input commit")
        text = text.replace("## Ledger\n", weakened + "\n## Ledger\n", 1)
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="disposition semantics"):
        module.validate_ledger_structure(_write_page(tmp_path, text))


def test_frozen_block_cannot_move_outside_frozen_section(tmp_path: pathlib.Path) -> None:
    module = _module()
    text = _initial_page_text()
    block = "### MC-001\n\n```text\n" + EXPECTED_FROZEN_TEXT[0] + "\n```\n\n"
    text = text.replace(block, "", 1).replace("## Frozen matched lines\n", block + "## Frozen matched lines\n", 1)
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="frozen|section|heading"):
        module.validate_ledger_structure(_write_page(tmp_path, text))


def test_extra_claim_heading_is_rejected(tmp_path: pathlib.Path) -> None:
    module = _module()
    text = _initial_page_text().replace("## Frozen matched lines\n", "## Frozen matched lines\n\n### MC-999\n", 1)
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
    page = _write_page(tmp_path, _replace_row(_initial_page_text(), 1, column, replacement))
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="row|immutable"):
        module.validate_ledger_structure(page)


@pytest.mark.parametrize("damage", ["missing", "duplicate", "out-of-order"])
def test_missing_duplicate_and_out_of_order_rows_are_rejected(tmp_path: pathlib.Path, damage: str) -> None:
    module = _module()
    text = _initial_page_text()
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
    text = _initial_page_text()
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
    text = _initial_page_text()
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
    page = _write_page(tmp_path, _replace_row(_initial_page_text(), 1, 5, disposition))
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="disposition"):
        module.validate_ledger_structure(page)


@pytest.mark.parametrize("done", [0, 10])
def test_in_audit_status_is_exact_for_observed_terminal_count(tmp_path: pathlib.Path, done: int) -> None:
    module = _module()
    text = _initial_page_text() if done == 0 else _terminalize(_initial_page_text(), done)
    rows = module.validate_ledger_structure(_write_page(tmp_path, text))
    assert sum(row["disposition"] != "pending" for row in rows) == done


def test_complete_status_is_required_at_fourteen_of_fourteen(tmp_path: pathlib.Path) -> None:
    module = _module()
    text = _terminalize(_initial_page_text(), 14, complete=True)
    assert len(module.validate_ledger_structure(_write_page(tmp_path, text))) == 14
    wrong = text.replace(module.COMPLETE_STATUS, module.IN_AUDIT_STATUS.format(done=14))
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="status"):
        module.validate_ledger_structure(_write_page(tmp_path, wrong))


@pytest.mark.parametrize("suffix", [" extra", "; opened 2026-08-28", "\n> **status: in audit — 0/14 Task-5 extractor rows dispositioned**"])
def test_status_prefix_suffix_and_duplicate_lines_are_rejected(tmp_path: pathlib.Path, suffix: str) -> None:
    module = _module()
    exact = module.IN_AUDIT_STATUS.format(done=0)
    page = _write_page(tmp_path, _initial_page_text().replace(exact, exact + suffix, 1))
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="status"):
        module.validate_ledger_structure(page)


def test_status_count_disagreement_is_rejected(tmp_path: pathlib.Path) -> None:
    module = _module()
    text = _replace_row(_initial_page_text(), 1, 5, "historical")
    text = _replace_row(text, 1, 6, "reviewed=MC-001")
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="status"):
        module.validate_ledger_structure(_write_page(tmp_path, text))


def test_pending_row_cannot_carry_evidence(tmp_path: pathlib.Path) -> None:
    module = _module()
    page = _write_page(tmp_path, _replace_row(_initial_page_text(), 1, 6, "unsupported evidence"))
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
    page = _write_page(tmp_path, _initial_page_text().replace(old, new, 1))
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match=message):
        module.validate_ledger_structure(page)


def test_missing_frozen_blobs_raise_named_error(tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch) -> None:
    module = _module()
    monkeypatch.setattr(module, "FROZEN_SOURCE_COMMIT", "0" * 40)
    page = _write_page(tmp_path, _initial_page_text())
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="Git|frozen"):
        module.validate_ledger_structure(page)
