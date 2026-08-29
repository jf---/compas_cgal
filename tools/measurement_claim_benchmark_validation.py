"""Strict Figure-6 and derived benchmark-claim validation."""

from __future__ import annotations

import json
import math
from typing import Dict
from typing import List
from typing import Mapping
from typing import Optional
from typing import Sequence
from typing import Tuple
from typing import cast

from tools.measurement_artifact import GitObjectId
from tools.measurement_artifact import ValidatedEnvelope
from tools.measurement_claim_benchmark_identity import EXTRACTION_COMMIT
from tools.measurement_claim_benchmark_identity import HISTORY_COMMIT
from tools.measurement_claim_benchmark_identity import MC013_HISTORY_COMMIT
from tools.measurement_claim_benchmark_schema import BenchmarkClaimPayload
from tools.measurement_claim_benchmark_schema import BenchmarkClaimRecord
from tools.measurement_claim_benchmark_schema import MC011ClaimPayload
from tools.measurement_claim_benchmark_schema import MC012ClaimPayload
from tools.measurement_claim_benchmark_schema import MC013ClaimPayload
from tools.measurement_claim_benchmark_schema import MC013SelectedValuesPayload
from tools.measurement_claim_benchmark_schema import MC014ClaimPayload
from tools.measurement_claim_benchmark_semantic_input import FIGURE6_CAPS
from tools.measurement_claim_benchmark_semantic_input import FIGURE6_SEMANTIC_COMMAND
from tools.measurement_claim_benchmark_semantic_input import FIGURE6_SPACINGS
from tools.measurement_claim_benchmark_semantic_input import figure6_config
from tools.measurement_claim_case_validation import fail_payload
from tools.measurement_claim_case_validation import validate_array
from tools.measurement_claim_case_validation import validate_literal
from tools.measurement_claim_case_validation import validate_object
from tools.measurement_claim_case_validation import validate_same
from tools.measurement_claim_errors import InvalidMeasurementClaimPayloadError
from tools.measurement_claim_schema import Degrees
from tools.measurement_claim_schema import ToolDiameters
from tools.measurement_claim_schema import ValidatedArtifactDirectory
from tools.measurement_claim_schema import ValidatedArtifactStartedUtc

_ROOT_KEYS = ("pocket", "engagement_measured_at_cap_deg", "points", "spacing_trials")
_POCKET_KEYS = ("name", "family", "tool_diameter", "params", "holes")
_POINT_KEYS = (
    "cap_deg",
    "controlled_length",
    "controlled_cut_motions",
    "controlled_entry_cuts",
    "controlled_max_tea_deg",
    "controlled_max_tea_after_entry_deg",
    "controlled_exceedances_after_entry",
    "controlled_meets_cap",
    "mathsm_spacing_tool_diameters",
    "mathsm_length",
    "mathsm_max_tea_after_entry_deg",
    "length_ratio",
)
_TRIAL_KEYS = ("spacing_tool_diameters", "length", "cut_motions", "entry_cuts", "max_tea_deg", "max_tea_after_entry_deg")
_PAYLOAD_KEYS = ("schema_version", "batch", "extraction_commit", "source_commit", "semantic_command", "config", "claims")
_COMPARISON_HEADER = (
    "| cap (deg) | controlled length | controlled max TEA after entry (deg) | controlled meets cap | "
    "controlled exceedances after entry | spacing (tool diam.) | constant-spacing length | "
    "constant-spacing max TEA after entry (deg) | length ratio |"
)
_COMPARISON_RULE = "| ---: | ---: | ---: | :--- | ---: | ---: | ---: | ---: | ---: |"
_TRIAL_HEADER = "| spacing (tool diam.) | length | cut motions | max TEA after entry (deg) | raw max TEA (deg) |"
_TRIAL_RULE = "| ---: | ---: | ---: | ---: | ---: |"


def _float(value: object, field: str) -> float:
    if type(value) is not float or not math.isfinite(value):
        fail_payload(field, "must be a finite exact float")
    return cast(float, value)


def _integer(value: object, field: str) -> int:
    if type(value) is not int or value < 0:
        fail_payload(field, "must be a non-negative exact integer")
    return cast(int, value)


def _optional_float(value: object, field: str) -> None:
    if value is not None:
        _float(value, field)


def _validate_holes(value: object) -> List[object]:
    holes = validate_array(value, "figure6.json.pocket.holes")
    for ring_index, ring_value in enumerate(holes):
        ring = validate_array(ring_value, f"figure6.json.pocket.holes[{ring_index}]")
        for point_index, point_value in enumerate(ring):
            field = f"figure6.json.pocket.holes[{ring_index}][{point_index}]"
            point = validate_array(point_value, field)
            if len(point) != 3:
                fail_payload(field, "must contain exactly three coordinates")
            for coordinate_index, coordinate in enumerate(point):
                if type(coordinate) not in (int, float) or not math.isfinite(cast(float, coordinate)):
                    fail_payload(f"{field}[{coordinate_index}]", "must be finite numeric")
    return holes


def _table_rows(lines: Sequence[str], header: str, rule: str, width: int, field: str) -> List[List[str]]:
    if lines.count(header) != 1 or lines.count(rule) != 1:
        fail_payload(field, "requires its exact header and rule once")
    index = lines.index(header)
    if index + 1 >= len(lines) or lines[index + 1] != rule:
        fail_payload(field, "rule must immediately follow its header")
    rows: List[List[str]] = []
    for line in lines[index + 2 :]:
        if not line.startswith("|"):
            break
        if not line.endswith("|"):
            fail_payload(field, "row must start and end with an outer pipe")
        cells = [cell.strip() for cell in line.strip()[1:-1].split("|")]
        if len(cells) != width:
            fail_payload(field, f"row must contain exactly {width} cells")
        rows.append(cells)
    return rows


def _section(lines: Sequence[str], heading: str, next_heading: Optional[str]) -> Sequence[str]:
    start = lines.index(heading) + 1
    end = lines.index(next_heading) if next_heading is not None else len(lines)
    return lines[start:end]


def _validate_figure6_markdown(markdown_bytes: bytes, *, pocket_name: str, point_count: int, trial_count: int) -> None:
    try:
        text = markdown_bytes.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise InvalidMeasurementClaimPayloadError("figure6.md: invalid UTF-8") from exc
    if not text.strip():
        fail_payload("figure6.md", "must not be empty")
    lines = text.splitlines()
    h1 = f"# Figure 6 reproduction — path length against engagement cap ({pocket_name})"
    if [line for line in lines if line.startswith("# ")] != [h1]:
        fail_payload("figure6.md", "requires the exact Figure-6 H1 once")
    headings = [line for line in lines if line.startswith("## ")]
    if headings != ["## Per-cap comparison", "## Constant-spacing trials"]:
        fail_payload("figure6.md", "requires exact ordered H2 headings")
    if lines.index(h1) > lines.index(headings[0]) or lines.index(headings[0]) > lines.index(headings[1]):
        fail_payload("figure6.md", "headings are reordered")
    comparison_section = _section(lines, headings[0], headings[1])
    trial_section = _section(lines, headings[1], None)
    comparison = _table_rows(comparison_section, _COMPARISON_HEADER, _COMPARISON_RULE, 9, "figure6.md comparison")
    trials = _table_rows(trial_section, _TRIAL_HEADER, _TRIAL_RULE, 5, "figure6.md trials")
    if point_count != 8 or len(comparison) != 8:
        fail_payload("figure6.md comparison", "requires exactly eight rows")
    if trial_count != 12 or len(trials) != 12:
        fail_payload("figure6.md trials", "requires exactly twelve rows")


def validate_figure6_payload(figure6_payload: object, figure6_markdown: bytes) -> Tuple[str, Degrees, Degrees]:
    root = validate_object(figure6_payload, _ROOT_KEYS, "figure6.json")
    pocket = validate_object(root["pocket"], _POCKET_KEYS, "figure6.json.pocket")
    validate_literal(pocket["name"], "rect_20x12", "figure6.json.pocket.name")
    validate_literal(pocket["family"], "analytic", "figure6.json.pocket.family")
    validate_literal(pocket["tool_diameter"], 2.0, "figure6.json.pocket.tool_diameter")
    validate_same(pocket["params"], {"width": 20.0, "height": 12.0}, "figure6.json.pocket.params")
    holes = _validate_holes(pocket["holes"])
    if holes:
        fail_payload("figure6.json.pocket.holes", "fixed Figure-6 config requires exactly []")
    validate_literal(root["engagement_measured_at_cap_deg"], 180.0, "figure6.json.engagement_measured_at_cap_deg")
    points = validate_array(root["points"], "figure6.json.points")
    if len(points) != len(FIGURE6_CAPS):
        fail_payload("figure6.json.points", "requires the exact cap count")
    for index, (item, cap) in enumerate(zip(points, FIGURE6_CAPS)):
        point = validate_object(item, _POINT_KEYS, f"figure6.json.points[{index}]")
        validate_literal(point["cap_deg"], cap, f"figure6.json.points[{index}].cap_deg")
        for name in ("controlled_length", "controlled_max_tea_deg", "controlled_max_tea_after_entry_deg"):
            _float(point[name], f"figure6.json.points[{index}].{name}")
        for name in ("controlled_cut_motions", "controlled_entry_cuts", "controlled_exceedances_after_entry"):
            _integer(point[name], f"figure6.json.points[{index}].{name}")
        if type(point["controlled_meets_cap"]) is not bool:
            fail_payload(f"figure6.json.points[{index}].controlled_meets_cap", "must be an exact boolean")
        for name in ("mathsm_spacing_tool_diameters", "mathsm_length", "mathsm_max_tea_after_entry_deg", "length_ratio"):
            _optional_float(point[name], f"figure6.json.points[{index}].{name}")
    trials = validate_array(root["spacing_trials"], "figure6.json.spacing_trials")
    if len(trials) != len(FIGURE6_SPACINGS):
        fail_payload("figure6.json.spacing_trials", "requires the exact spacing count")
    selected: Dict[float, float] = {}
    for index, (item, spacing) in enumerate(zip(trials, FIGURE6_SPACINGS)):
        trial = validate_object(item, _TRIAL_KEYS, f"figure6.json.spacing_trials[{index}]")
        validate_literal(trial["spacing_tool_diameters"], spacing, f"figure6.json.spacing_trials[{index}].spacing_tool_diameters")
        for name in ("length", "max_tea_deg", "max_tea_after_entry_deg"):
            _float(trial[name], f"figure6.json.spacing_trials[{index}].{name}")
        _integer(trial["cut_motions"], f"figure6.json.spacing_trials[{index}].cut_motions")
        _integer(trial["entry_cuts"], f"figure6.json.spacing_trials[{index}].entry_cuts")
        if spacing in (0.025, 0.1):
            selected[spacing] = cast(float, trial["max_tea_after_entry_deg"])
    if set(selected) != {0.025, 0.1} or selected[0.025] <= selected[0.1]:
        fail_payload("figure6.json.spacing_trials", "requires unique 0.025 > 0.1 selected measurements")
    _validate_figure6_markdown(figure6_markdown, pocket_name=cast(str, pocket["name"]), point_count=len(points), trial_count=len(trials))
    return cast(str, pocket["name"]), Degrees(selected[0.025]), Degrees(selected[0.1])


def compose_benchmark_payload(*, source_commit: GitObjectId, figure6_payload: object, figure6_markdown: bytes) -> BenchmarkClaimPayload:
    _, fine, comparison = validate_figure6_payload(figure6_payload, figure6_markdown)
    selected = MC013SelectedValuesPayload(
        fine_spacing=ToolDiameters(0.025),
        fine_max_tea_after_entry=fine,
        comparison_spacing=ToolDiameters(0.1),
        comparison_max_tea_after_entry=comparison,
        angle_unit="degree",
        spacing_unit="tool-diameter",
    )
    claims: List[BenchmarkClaimRecord] = [
        MC011ClaimPayload(
            claim_id="MC-011",
            source="benchmarks/gate.py:58",
            disposition="not-a-claim",
            source_commit=source_commit,
            history_commit=GitObjectId(HISTORY_COMMIT),
            reason="Pocket dimensions identify where the defect was first observed; they are not a measured result.",
            missing_inputs=[],
        ),
        MC012ClaimPayload(
            claim_id="MC-012",
            source="benchmarks/gate.py:66",
            disposition="deleted",
            source_commit=source_commit,
            history_commit=GitObjectId(HISTORY_COMMIT),
            reason="The anecdote is under-specified and cannot be reconstructed without guessing inputs.",
            missing_inputs=["generator", "circle-selection", "entry-treatment", "operation-enumeration"],
        ),
        MC013ClaimPayload(
            claim_id="MC-013",
            source="benchmarks/mathsm.py:47",
            disposition="corrected",
            source_commit=source_commit,
            history_commit=GitObjectId(MC013_HISTORY_COMMIT),
            reason="Authenticated Figure-6 JSON establishes the non-monotone spacing observation.",
            missing_inputs=[],
            selected_values=selected,
        ),
        MC014ClaimPayload(
            claim_id="MC-014",
            source="benchmarks/quality.py:150",
            disposition="not-a-claim",
            source_commit=source_commit,
            history_commit=GitObjectId(HISTORY_COMMIT),
            reason="The shared cap defines metric comparability and reports no empirical value.",
            missing_inputs=[],
        ),
    ]
    return BenchmarkClaimPayload(
        schema_version="measurement-claim-payload/v1",
        batch="benchmark",
        extraction_commit=GitObjectId(EXTRACTION_COMMIT),
        source_commit=source_commit,
        semantic_command=list(FIGURE6_SEMANTIC_COMMAND),
        config=figure6_config(),
        claims=claims,
    )


def validate_benchmark_payload(payload: object, *, figure6_payload: object, figure6_markdown: bytes) -> BenchmarkClaimPayload:
    root = validate_object(payload, _PAYLOAD_KEYS, "benchmark-claims.json")
    source = root["source_commit"]
    if type(source) is not str or len(source) not in (40, 64) or any(character not in "0123456789abcdef" for character in source):
        fail_payload("benchmark-claims.json.source_commit", "must be a full lowercase Git object ID")
    expected = compose_benchmark_payload(source_commit=GitObjectId(cast(str, source)), figure6_payload=figure6_payload, figure6_markdown=figure6_markdown)
    validate_same(root, expected, "benchmark-claims.json")
    return cast(BenchmarkClaimPayload, root)


def render_benchmark_ledger_evidence(
    payload: BenchmarkClaimPayload,
    envelope: ValidatedEnvelope,
    *,
    started: ValidatedArtifactStartedUtc,
    artifact_directory: ValidatedArtifactDirectory,
) -> Mapping[str, str]:
    if payload["source_commit"] != envelope.commit:
        fail_payload("benchmark-claims.json.source_commit", "must equal authenticated envelope commit")
    artifact = f"{artifact_directory}/benchmark-claims.json@sha256:{envelope.payload_sha256['benchmark-claims.json']}"
    semantic_command = json.dumps(payload["semantic_command"], sort_keys=True, separators=(",", ":"), allow_nan=False)
    config = json.dumps(payload["config"], sort_keys=True, separators=(",", ":"), allow_nan=False)
    values: Dict[str, str] = {}
    for claim in payload["claims"]:
        record = cast(Dict[str, object], claim)
        selected = record.get("selected_values")
        cell = (
            f"artifact={artifact}; started={started.isoformat(timespec='microseconds')}; claim={record['claim_id']}; "
            f"input={envelope.input_sha256}; result={envelope.result_sha256}; source={record['source_commit']}; "
            f"history={record['history_commit']}; disposition={record['disposition']}; reason={record['reason']}; "
            f"semantic_command={semantic_command}; config={config}; "
            f"missing_inputs={json.dumps(record['missing_inputs'], separators=(',', ':'))}"
        )
        if selected is not None:
            cell += f"; selected_values={json.dumps(selected, sort_keys=True, separators=(',', ':'))}"
        if any(token in cell for token in ("|", "\r", "\n", "\u2028", "\u2029")):
            fail_payload(str(record["claim_id"]), "rendered evidence is not Markdown-safe")
        values[cast(str, record["claim_id"])] = cell
    return values
