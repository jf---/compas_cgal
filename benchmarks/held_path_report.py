"""Pure deterministic Markdown projection of Held Figure 5 evidence."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from datetime import timedelta
from html import escape
from typing import Callable
from typing import Union

from typing_extensions import Self
from typing_extensions import TypeAlias

from benchmarks.errors import InvalidHeldPathEvidenceError
from benchmarks.errors import InvalidHeldPathReportContextError
from benchmarks.held_path_evidence import HeldFigure5Characterization
from benchmarks.held_post_qualification import post_qualification_failures
from benchmarks.quality_observations import CRITERION_NAMES
from benchmarks.quality_observations import CountCriterion
from benchmarks.quality_observations import DegreesCriterion
from benchmarks.quality_observations import FractionCriterion
from benchmarks.quality_observations import PathQualityAssessment
from benchmarks.quality_observations import ToolRadiusMultipleCriterion

Criterion: TypeAlias = Union[FractionCriterion, CountCriterion, DegreesCriterion, ToolRadiusMultipleCriterion]
_CRITERION_ACCESSORS: tuple[tuple[str, Callable[[PathQualityAssessment], Criterion]], ...] = (
    ("uncut fraction", lambda assessment: assessment.uncut_fraction),
    ("gouging motions", lambda assessment: assessment.gouging_motions),
    ("unsafe rapids", lambda assessment: assessment.unsafe_rapids),
    ("continuity breaks", lambda assessment: assessment.continuity_breaks),
    ("zero-length motions", lambda assessment: assessment.zero_length_motions),
    ("degenerate loops", lambda assessment: assessment.degenerate_loops),
    ("redundant operations", lambda assessment: assessment.redundant_operations),
    ("cap exceedances", lambda assessment: assessment.cap_exceedances),
    ("slotting motions", lambda assessment: assessment.slotting_motions),
    ("max engagement step (deg)", lambda assessment: assessment.max_engagement_step),
    ("max loop radius step (tool radii)", lambda assessment: assessment.max_loop_radius_step),
    ("tangent breaks", lambda assessment: assessment.tangent_breaks),
)

_REPORT_ONLY_FIELDS: tuple[tuple[str, tuple[tuple[str, str], ...]], ...] = (
    (
        "elementary",
        (
            ("gouge_free", "boolean"),
            ("rapid_safety", "boolean"),
            ("marginal_loops", "count"),
            ("recut_fraction", "fraction"),
        ),
    ),
    (
        "cut",
        (
            ("max_engagement_deg", "degrees"),
            ("engagement_p95_deg", "degrees"),
            ("engagement_variance_deg2", "degrees^2"),
            ("max_chip_thickness_ratio", "ratio"),
            ("low_chip_thickness_ratio", "ratio"),
            ("max_engagement_gradient_deg_per_length", "degrees per benchmark-normalized mm"),
            ("immersion_steady_fraction", "fraction"),
            ("immersion_at_design_fraction", "fraction"),
            ("immersion_excursions", "count"),
            ("mean_radial_depth", "benchmark-normalized mm"),
            ("radial_depth_variance", "benchmark-normalized mm^2"),
            ("wall_scallop_height", "benchmark-normalized mm"),
            ("loop_radius_cv", "ratio"),
        ),
    ),
    (
        "speed",
        (
            ("cutting_length", "benchmark-normalized mm"),
            ("air_length", "benchmark-normalized mm"),
            ("air_fraction", "fraction"),
            ("max_curvature", "inverse benchmark-normalized mm"),
            ("curvature_breaks", "count"),
            ("direction_reversals", "count"),
            ("retract_count", "count"),
            ("reentry_count", "count"),
        ),
    ),
    (
        "longevity",
        (
            ("material_entries", "count"),
            ("cut_air_alternations", "count"),
            ("alternations_per_length", "inverse benchmark-normalized mm"),
        ),
    ),
    (
        "program",
        (
            ("block_count", "count"),
            ("cut_blocks", "count"),
            ("min_block_length", "benchmark-normalized mm"),
            ("median_block_length", "benchmark-normalized mm"),
            ("short_block_length", "benchmark-normalized mm"),
            ("block_length_cv", "ratio"),
            ("blocks_per_unit_length", "inverse benchmark-normalized mm"),
            ("arc_length_fraction", "fraction"),
        ),
    ),
)


@dataclass(frozen=True, init=False)
class HeldPathReportContext:
    """Validated invocation metadata for one deterministic report."""

    generated_at_utc: datetime
    pixi_command: str
    generator_policy_name: str

    def __init__(self) -> None:
        raise TypeError("HeldPathReportContext must be created with HeldPathReportContext.build().")

    @classmethod
    def build(cls, *, generated_at_utc: datetime, pixi_command: str, generator_policy_name: str) -> Self:
        if type(generated_at_utc) is not datetime or generated_at_utc.utcoffset() != timedelta(0):
            raise InvalidHeldPathReportContextError("report generation time must be one aware UTC datetime.")
        if type(pixi_command) is not str or not pixi_command.strip():
            raise InvalidHeldPathReportContextError("report Pixi command must be one non-empty string.")
        if type(generator_policy_name) is not str or not generator_policy_name.strip():
            raise InvalidHeldPathReportContextError("report generator policy must be one non-empty string.")
        context = object.__new__(cls)
        object.__setattr__(context, "generated_at_utc", generated_at_utc)
        object.__setattr__(context, "pixi_command", pixi_command)
        object.__setattr__(context, "generator_policy_name", generator_policy_name)
        return context


def _markdown_table_cell(value: object) -> str:
    text = str(value).replace("\r\n", "\n").replace("\r", "\n")
    escaped = escape(text, quote=True)
    return escaped.replace("\\", "\\\\").replace("|", "\\|").replace("\n", "<br>")


def _number(value: object) -> str:
    if type(value) is bool:
        return "true" if value else "false"
    if type(value) is int:
        return str(value)
    if not isinstance(value, float):
        raise TypeError("report numeric fields require validated integers or binary64 values.")
    return f"{float(value):.6f}"


def _row(label: object, value: object) -> str:
    return f"| {_markdown_table_cell(label)} | {_markdown_table_cell(value)} |"


def render_held_figure5_2d_path(characterization: HeldFigure5Characterization, context: HeldPathReportContext) -> str:
    """Render validated evidence without geometry, clocks, writes, or new decisions."""
    if type(characterization) is not HeldFigure5Characterization:
        raise InvalidHeldPathEvidenceError("Held report rendering requires one validated HeldFigure5Characterization.")
    if type(context) is not HeldPathReportContext:
        raise InvalidHeldPathReportContextError("Held report rendering requires one validated HeldPathReportContext.")
    failures = post_qualification_failures(characterization)
    lines = [
        "# Held Figure 5 2D path characterization",
        "",
        "**HISTORICAL 2D BENCHMARK CHARACTERIZATION - NOT A MANUFACTURING RELEASE**",
        "",
        "## Run context",
        "",
        "| Field | Value |",
        "| --- | --- |",
        _row("Generated at UTC", context.generated_at_utc.isoformat(timespec="seconds")),
        _row("Pixi invocation", context.pixi_command),
        _row("Generator policy", context.generator_policy_name),
        "",
        "## Claim boundaries",
        "",
        "Benchmark coordinates and lengths are normalized-scale millimetres; they are not a machine setup.",
        "Sampled negative observations are bounded evidence, not global proofs.",
        "The exercised guarded replay is not evidence of native audit compatibility.",
        "No controller, machine, setup, or manufacturing-release claim is made.",
        "",
        "## Case and source",
        "",
        "| Field | Value |",
        "| --- | --- |",
        _row("Case", characterization.case_name),
        _row("Tool diameter (benchmark-normalized mm)", _number(characterization.tool_diameter)),
        _row("TEA cap (degrees)", _number(characterization.tea_cap)),
        _row("Reference primitives", _number(characterization.reference_primitive_count)),
        _row("Projection vertices", _number(characterization.projection_vertex_count)),
        "",
        "## Path and timing",
        "",
        "| Field | Value |",
        "| --- | --- |",
        _row("Total source operations", _number(characterization.source_operation_count)),
        _row("TEA-audited lateral operations", _number(characterization.tea_audited_operation_count)),
        _row("TEA-audit-excluded operations", _number(characterization.excluded_operation_count)),
        _row("Sampled material-contact operations", _number(characterization.sampled_material_contact_operations)),
        _row("Generation (seconds)", _number(characterization.generation_seconds)),
        _row("Guarded audit (seconds)", _number(characterization.audit_seconds)),
        _row("Survey (seconds)", _number(characterization.survey_seconds)),
        _row("Coverage plus reduction (seconds)", _number(characterization.reduction_seconds)),
        "",
        "## Engagement dispositions",
        "",
        "| Disposition | Count |",
        "| --- | --- |",
        _row("certified", _number(characterization.engagement.certified_count)),
        _row("demonstrated_exceeded", _number(characterization.engagement.demonstrated_exceeded_count)),
        _row("unresolved", _number(characterization.engagement.unresolved_count)),
        "",
        "## Exceedance witnesses",
        "",
        "| Operation index | World X (benchmark-normalized mm) | World Y (benchmark-normalized mm) |",
        "| --- | --- | --- |",
    ]
    if characterization.witnesses:
        for witness in characterization.witnesses:
            lines.append(
                f"| {_markdown_table_cell(_number(witness.operation_index))} | "
                f"{_markdown_table_cell(_number(witness.position.x))} | {_markdown_table_cell(_number(witness.position.y))} |"
            )
    else:
        lines.append("| none | no sampled exact-predicate exceedance witness | no sampled exact-predicate exceedance witness |")
    lines.extend(
        [
            "",
            "## Quality criteria",
            "",
            "| Criterion | Measured | Required | Evidence kind | Outcome |",
            "| --- | --- | --- | --- | --- |",
        ]
    )
    assert tuple(name for name, _ in _CRITERION_ACCESSORS) == CRITERION_NAMES
    for name, accessor in _CRITERION_ACCESSORS:
        criterion = accessor(characterization.assessment)
        lines.append(
            f"| {_markdown_table_cell(name)} | {_markdown_table_cell(_number(criterion.measured))} | "
            f"{_markdown_table_cell(_number(criterion.required))} | {_markdown_table_cell(criterion.evidence)} | "
            f"{_markdown_table_cell(criterion.outcome)} |"
        )
    lines.extend(["", "## Report-only PathQuality", "", "| Group | Field | Unit | Value |", "| --- | --- | --- | --- |"])
    quality = characterization.path_quality
    for group_name, fields in _REPORT_ONLY_FIELDS:
        group = getattr(quality, group_name)
        for field_name, unit in fields:
            lines.append(
                f"| {_markdown_table_cell(group_name)} | {_markdown_table_cell(field_name)} | "
                f"{_markdown_table_cell(unit)} | {_markdown_table_cell(_number(getattr(group, field_name)))} |"
            )
        if group_name == "longevity":
            for index, (low, high, length) in enumerate(group.engagement_length_histogram):
                value = f"{_number(low)} to {_number(high)} degrees: {_number(length)}"
                lines.append(f"| longevity | engagement_length_histogram[{index}] | band: degrees; length: benchmark-normalized mm | {_markdown_table_cell(value)} |")
    lines.append(f"| top_level | cut_operations | count | {_markdown_table_cell(_number(quality.cut_operations))} |")
    lines.append(f"| top_level | path_length | benchmark-normalized mm | {_markdown_table_cell(_number(quality.path_length))} |")
    lines.extend(["", "## Post-qualification entry verdict", ""])
    if failures:
        lines.append("**not eligible for postprocessor qualification**")
        lines.extend(["", "| Ordered refusal |", "| --- |"])
        for failure in failures:
            lines.append(f"| {_markdown_table_cell(failure)} |")
    else:
        lines.append("**eligible for postprocessor qualification evidence entry**")
    lines.extend(["", "This verdict is not a manufacturing release."])
    return "\n".join(lines) + "\n"
