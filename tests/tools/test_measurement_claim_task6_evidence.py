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


SOURCE_COMMIT = "a" * 40
SOURCE_CORRECTION_COMMIT = "53135e04390e84bf69aa74dc4d0c1ce6ca308eb4"
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
        "schema_version": "measurement-claim-payload/v2",
        "batch": "generator",
        "extraction_commit": FROZEN_COMMIT,
        "source_commit": SOURCE_COMMIT,
        "source_correction_commit": SOURCE_CORRECTION_COMMIT,
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
    payload: object = _claim_payload()
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
        f'source={SOURCE_COMMIT}; source_correction={SOURCE_CORRECTION_COMMIT}; case=radial-station; disposition=re-earned; reason=reason 1; config={{"a":[1,true],"z":1}}; '
        'evidence={"ordinal":1,"value":1.25}; selection={"policy":"policy-1"}; continuous_certificate=null'
    )
    assert all("\n" not in cell and "\r" not in cell and "|" not in cell for cell in rendered.values())


def test_semantic_repair_renderer_carries_execution_and_correction_identities(monkeypatch: pytest.MonkeyPatch) -> None:
    module = _module()
    payload = _claim_payload()
    payload["source_correction_commit"] = SOURCE_CORRECTION_COMMIT
    _trust_payload_validator(monkeypatch, module)

    rendered = module.render_ledger_evidence(payload, _envelope(), started=STARTED, artifact_directory=ARTIFACT_DIRECTORY)

    assert all(f"source={SOURCE_COMMIT}; source_correction={SOURCE_CORRECTION_COMMIT};" in cell for cell in rendered.values())


def test_task6_private_row_comparator_accepts_exactly_ten_of_fourteen(monkeypatch: pytest.MonkeyPatch) -> None:
    module = _module()
    payload = _claim_payload()
    _trust_payload_validator(monkeypatch, module)
    rendered = module.render_ledger_evidence(payload, _envelope(), started=STARTED, artifact_directory=ARTIFACT_DIRECTORY)
    rows = _task6_rows(payload, rendered)

    assert module._validate_task6_ledger_rows(rows, payload, _envelope(), started=STARTED, artifact_directory=ARTIFACT_DIRECTORY) is None
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
        module._validate_task6_ledger_rows(tuple(rows), payload, _envelope(), started=STARTED, artifact_directory=ARTIFACT_DIRECTORY)


def test_public_task6_ledger_consumer_owns_structure_and_artifact_authentication(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: pathlib.Path,
) -> None:
    module = _module()
    payload = _claim_payload()
    _trust_payload_validator(monkeypatch, module)
    envelope = _envelope()
    rendered = module.render_ledger_evidence(payload, envelope, started=STARTED, artifact_directory=ARTIFACT_DIRECTORY)
    rows = _task6_rows(payload, rendered)
    ledger = tmp_path / "measurement_claims.md"
    artifact = tmp_path / "benchmarks" / "measurement_claim_results" / ARTIFACT_DIRECTORY.name
    calls: list[tuple[str, pathlib.Path]] = []

    def validate_structure(candidate: pathlib.Path) -> tuple[Any, ...]:
        calls.append(("ledger", candidate))
        return rows

    def validate_artifact(candidate: pathlib.Path) -> tuple[Any, ...]:
        calls.append(("artifact", candidate))
        return payload, envelope, STARTED, ARTIFACT_DIRECTORY

    def validate_source(
        repository: pathlib.Path,
        correction_commit: str,
        *,
        mc007_disposition: str,
    ) -> None:
        assert repository == tmp_path
        assert correction_commit == SOURCE_CORRECTION_COMMIT
        assert mc007_disposition == "corrected"
        calls.append(("source", repository))

    def validate_lineage(
        repository: pathlib.Path,
        correction_commit: str,
        execution_commit: str,
        *,
        mc007_disposition: str,
    ) -> None:
        assert repository == tmp_path
        assert correction_commit == SOURCE_CORRECTION_COMMIT
        assert execution_commit == SOURCE_COMMIT
        assert mc007_disposition == "corrected"
        calls.append(("lineage", repository))

    def compare_rows(*values: Any, **named: Any) -> None:
        del values, named
        calls.append(("rows", ledger))

    monkeypatch.setattr(module, "validate_ledger_structure", validate_structure)
    monkeypatch.setattr(module, "validate_claim_artifact", validate_artifact)
    monkeypatch.setattr(module, "validate_task6_source_correction", validate_source)
    monkeypatch.setattr(module, "validate_task6_source_lineage", validate_lineage)
    monkeypatch.setattr(module, "_validate_task6_ledger_rows", compare_rows)

    assert module.validate_ledger_evidence(ledger, [artifact]) is None
    assert calls == [("ledger", ledger), ("artifact", artifact), ("source", tmp_path), ("lineage", tmp_path), ("rows", ledger)]


@pytest.mark.parametrize("artifact_count", [0, 2])
def test_public_task6_ledger_consumer_rejects_noncanonical_artifact_cardinality_before_authentication(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: pathlib.Path,
    artifact_count: int,
) -> None:
    module = _module()
    calls: list[pathlib.Path] = []

    def validate_artifact(candidate: pathlib.Path) -> tuple[Any, ...]:
        calls.append(candidate)
        raise AssertionError("artifact authentication must not run")

    monkeypatch.setattr(module, "validate_claim_artifact", validate_artifact)
    artifacts = [tmp_path / f"artifact-{index}" for index in range(artifact_count)]

    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="exactly one.*artifact"):
        module.validate_ledger_evidence(tmp_path / "measurement_claims.md", artifacts)
    assert calls == []


def test_public_task6_ledger_consumer_does_not_authenticate_artifact_after_structure_failure(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: pathlib.Path,
) -> None:
    module = _module()
    ledger_error = module.InvalidMeasurementClaimLedgerError("invalid ledger")
    artifact_calls: list[pathlib.Path] = []

    def reject_structure(candidate: pathlib.Path) -> tuple[Any, ...]:
        assert candidate == tmp_path / "measurement_claims.md"
        raise ledger_error

    def validate_artifact(candidate: pathlib.Path) -> tuple[Any, ...]:
        artifact_calls.append(candidate)
        raise AssertionError("artifact authentication must not run")

    monkeypatch.setattr(module, "validate_ledger_structure", reject_structure)
    monkeypatch.setattr(module, "validate_claim_artifact", validate_artifact)

    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="invalid ledger"):
        module.validate_ledger_evidence(tmp_path / "measurement_claims.md", [tmp_path / "artifact"])
    assert artifact_calls == []


def test_public_task6_ledger_consumer_does_not_compare_unauthenticated_artifact(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: pathlib.Path,
) -> None:
    module = _module()
    result_module = importlib.import_module("tools.measurement_claim_result")
    rows: tuple[Any, ...] = ()
    artifact_calls: list[pathlib.Path] = []
    comparator_calls: list[tuple[Any, ...]] = []

    monkeypatch.setattr(module, "validate_ledger_structure", lambda candidate: rows)

    def reject_artifact(candidate: pathlib.Path) -> tuple[Any, ...]:
        artifact_calls.append(candidate)
        raise result_module.InvalidMeasurementClaimPayloadError("unauthenticated")

    def compare(*values: Any, **named: Any) -> None:
        comparator_calls.append((*values, named))

    monkeypatch.setattr(module, "validate_claim_artifact", reject_artifact)
    monkeypatch.setattr(module, "_validate_task6_ledger_rows", compare)
    artifact = tmp_path / "artifact"

    with pytest.raises(result_module.InvalidMeasurementClaimPayloadError, match="unauthenticated"):
        module.validate_ledger_evidence(tmp_path / "measurement_claims.md", [artifact])
    assert artifact_calls == [artifact]
    assert comparator_calls == []


def test_public_task6_ledger_consumer_does_not_compare_source_unverified_artifact(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: pathlib.Path,
) -> None:
    module = _module()
    payload = _claim_payload()
    _trust_payload_validator(monkeypatch, module)
    source_calls: list[tuple[pathlib.Path, str, str]] = []
    comparator_calls: list[tuple[Any, ...]] = []
    ledger = tmp_path / "docs" / "measurement_claims.md"
    artifact = tmp_path / "benchmarks" / "measurement_claim_results" / ARTIFACT_DIRECTORY.name
    monkeypatch.setattr(module, "validate_ledger_structure", lambda candidate: ())
    monkeypatch.setattr(
        module,
        "validate_claim_artifact",
        lambda candidate: (payload, _envelope(), STARTED, ARTIFACT_DIRECTORY),
    )

    def reject_source(repository: pathlib.Path, correction_commit: str, *, mc007_disposition: str) -> None:
        source_calls.append((repository, correction_commit, mc007_disposition))
        raise module.InvalidMeasurementClaimLedgerError("source unverified")

    def compare(*values: Any, **named: Any) -> None:
        comparator_calls.append((*values, named))

    monkeypatch.setattr(module, "validate_task6_source_correction", reject_source)
    monkeypatch.setattr(module, "_validate_task6_ledger_rows", compare)

    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="source unverified"):
        module.validate_ledger_evidence(ledger, [artifact])

    assert source_calls == [(tmp_path, SOURCE_CORRECTION_COMMIT, "corrected")]
    assert comparator_calls == []


def test_semantic_repair_public_consumer_stops_before_rows_when_lineage_fails(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: pathlib.Path,
) -> None:
    module = _module()
    payload = _claim_payload()
    ledger = tmp_path / "docs" / "measurement_claims.md"
    artifact = tmp_path / "benchmarks" / "measurement_claim_results" / ARTIFACT_DIRECTORY.name
    events: list[str] = []
    monkeypatch.setattr(module, "validate_ledger_structure", lambda candidate: ())
    monkeypatch.setattr(module, "validate_claim_artifact", lambda candidate: (payload, _envelope(), STARTED, ARTIFACT_DIRECTORY))
    monkeypatch.setattr(module, "validate_task6_source_correction", lambda *args, **kwargs: events.append("correction"))

    def reject_lineage(*args: Any, **kwargs: Any) -> None:
        del args, kwargs
        events.append("lineage")
        raise module.InvalidMeasurementClaimLedgerError("lineage unverified")

    monkeypatch.setattr(module, "validate_task6_source_lineage", reject_lineage)
    monkeypatch.setattr(module, "_validate_task6_ledger_rows", lambda *args, **kwargs: events.append("rows"))

    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="lineage unverified"):
        module.validate_ledger_evidence(ledger, [artifact])
    assert events == ["correction", "lineage"]


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


def test_measurement_claim_task6_ledger_module_exists() -> None:
    assert importlib.util.find_spec("tools.measurement_claim_task6_ledger") is not None


def test_committed_task6_ledger_has_authenticated_acceptance(tmp_path: pathlib.Path) -> None:
    module = _module()
    raw = subprocess.run(
        [
            "git",
            "--no-replace-objects",
            "show",
            "a8b6adbb40d2bea630068990828119011ece156c:docs/measurement_claims.md",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        stdout=subprocess.PIPE,
    ).stdout
    ledger = tmp_path / "measurement_claims.md"
    ledger.write_bytes(raw)
    text = raw.decode("utf-8")
    artifact_root = PROJECT_ROOT / "benchmarks" / "measurement_claim_results"
    artifacts = tuple(sorted(path for path in artifact_root.glob("*-generator-*") if path.is_dir()))
    assert len(artifacts) == 1
    assert module.validate_ledger_evidence(ledger, artifacts) is None
    rows = module.validate_ledger_structure(ledger)
    assert text.startswith("# Measurement-claim ledger\n\n")
    assert SCOPE_BLOCK in text
    assert EXTRACTOR_BLOCK in text
    assert ORDER_BLOCK in text
    assert DRIFT_BLOCK in text
    assert DISPOSITION_BLOCK in text
    assert module.FROZEN_SOURCE_COMMIT == FROZEN_COMMIT
    assert tuple((row["ordinal"], row["claim_id"], row["extracted_location"], row["stable_anchor"], row["anchor_match"]) for row in rows) == EXPECTED_ROWS
    assert sum(row["disposition"] != "pending" for row in rows) == 10
    assert tuple((row["disposition"], row["evidence"]) for row in rows[10:]) == (("pending", "—"),) * 4
    assert _immutable_digest(rows) == module.IMMUTABLE_COLUMNS_SHA256
