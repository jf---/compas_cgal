from __future__ import annotations

import time
from datetime import datetime
from datetime import timedelta
from datetime import timezone
from pathlib import Path

import pytest

from benchmarks.errors import HeldPathNotEligibleForPostQualificationError
from benchmarks.held_path_characterize import CharacterizationPhase
from benchmarks.held_post_qualification import post_qualification_failures
from benchmarks.held_post_qualification import require_post_qualification_candidate
from benchmarks.quality_observations import CRITERION_NAMES
from benchmarks.quality_observations import EVIDENCE_BY_CRITERION
from tools.held_path_characterization import DEFAULT_REPORT_PATH
from tools.held_path_characterization import write_held_figure5_report

LIVE_COMMAND = "pixi run held-figure5-characterize-live"
REPORT_PATH = Path("docs/benchmarks/held_figure5_2d_path.md")
LIVE_RUN_LEDGER_PATH = Path("docs/superpowers/state/held-figure5-live-run.md")
EXPECTED_SOURCE_OPERATIONS = 2_289
EXPECTED_PHASES: tuple[CharacterizationPhase, ...] = (
    "generation",
    "guarded_audit",
    "survey",
    "quality_reduction",
)
EXIT_ZERO_BOUNDARY = "Exit zero means evidence completion, not path eligibility or machine release."
ALLOWED_OUTCOMES = frozenset(
    {
        "no_failure_observed",
        "failure_observed",
        "criterion_satisfied",
        "criterion_violated",
        "within_declared_tolerance",
        "outside_declared_tolerance",
    }
)
EXPECTED_EVIDENCE_BY_CRITERION = {
    "uncut fraction": "sampled_diagnostic",
    "gouging motions": "sampled_diagnostic",
    "unsafe rapids": "tolerance_diagnostic",
    "continuity breaks": "tolerance_diagnostic",
    "zero-length motions": "derived_geometry",
    "degenerate loops": "derived_geometry",
    "redundant operations": "exact_depletion_replay",
    "cap exceedances": "sampled_exact_predicate",
    "slotting motions": "sampled_exact_predicate",
    "max engagement step (deg)": "sampled_diagnostic",
    "max loop radius step (tool radii)": "derived_geometry",
    "tangent breaks": "tolerance_diagnostic",
}


def _write_live_run_ledger(
    *,
    run_started_at: datetime,
    outcome: str,
    observed_phases: tuple[CharacterizationPhase, ...],
    elapsed_seconds: float,
    prior_report_generated_at: datetime | None,
    report_generated_at: datetime | None,
    command_result: str,
    qualification_verdict: str,
    qualification_failures: tuple[str, ...] | None,
) -> None:
    phase = "none" if not observed_phases else observed_phases[-1]
    phase_history = "none" if not observed_phases else " -> ".join(observed_phases)
    prior_report = "No report existed before this run." if prior_report_generated_at is None else f"Prior/stale evidence generated at {prior_report_generated_at.isoformat()}."
    report_provenance = "not available" if report_generated_at is None else f"{report_generated_at.isoformat()}; current evidence"
    rendered_failures = "pending" if qualification_failures is None else "none" if not qualification_failures else " | ".join(qualification_failures)
    markdown = "\n".join(
        (
            "# Held Figure 5 live run",
            "",
            f"- Run start UTC: {run_started_at.isoformat()}",
            f"- Outcome: {outcome}",
            f"- Active phase: {phase}",
            f"- Observed phases: {phase_history}",
            f"- Elapsed seconds: {elapsed_seconds:.6f}",
            f"- Command: `{LIVE_COMMAND}`",
            f"- Command result: {command_result}",
            f"- Report path: `{REPORT_PATH}`",
            f"- Pre-run report state: {prior_report}",
            f"- Report generated UTC: {report_provenance}",
            f"- Qualification verdict: {qualification_verdict}",
            f"- Qualification failures: {rendered_failures}",
            "",
            "This ledger records evidence execution only; geometric claims remain bounded by the report.",
            "",
        )
    )
    LIVE_RUN_LEDGER_PATH.write_text(markdown, encoding="utf-8")


def _rendered_number(value: object) -> str:
    if type(value) is bool:
        return "true" if value else "false"
    if type(value) is int:
        return str(value)
    assert isinstance(value, float)
    return f"{value:.6f}"


def _report_utc(markdown: str) -> datetime:
    prefix = "| Generated at UTC | "
    rows = [line for line in markdown.splitlines() if line.startswith(prefix)]
    assert len(rows) == 1
    generated_at = datetime.fromisoformat(rows[0].removeprefix(prefix).removesuffix(" |"))
    assert generated_at.utcoffset() == timedelta(0)
    return generated_at


def test_held_figure5_path_live_oracle() -> None:
    assert REPORT_PATH == DEFAULT_REPORT_PATH
    run_started_at = datetime.now(timezone.utc).replace(microsecond=0)
    operator_started = time.monotonic()
    prior_report_generated_at = _report_utc(REPORT_PATH.read_text(encoding="utf-8")) if REPORT_PATH.exists() else None
    observed_phases: list[CharacterizationPhase] = []
    _write_live_run_ledger(
        run_started_at=run_started_at,
        outcome="pending",
        observed_phases=(),
        elapsed_seconds=0.0,
        prior_report_generated_at=prior_report_generated_at,
        report_generated_at=None,
        command_result="pending",
        qualification_verdict="pending",
        qualification_failures=None,
    )

    def observe_phase(phase: CharacterizationPhase) -> None:
        assert len(observed_phases) < len(EXPECTED_PHASES)
        expected_phase = EXPECTED_PHASES[len(observed_phases)]
        assert phase == expected_phase
        observed_phases.append(phase)
        assert tuple(observed_phases) == EXPECTED_PHASES[: len(observed_phases)]
        _write_live_run_ledger(
            run_started_at=run_started_at,
            outcome="running",
            observed_phases=tuple(observed_phases),
            elapsed_seconds=time.monotonic() - operator_started,
            prior_report_generated_at=prior_report_generated_at,
            report_generated_at=None,
            command_result="running",
            qualification_verdict="pending",
            qualification_failures=None,
        )

    characterization = write_held_figure5_report(
        pixi_command=LIVE_COMMAND,
        phase_observer=observe_phase,
        path=REPORT_PATH,
    )
    assert tuple(observed_phases) == EXPECTED_PHASES

    engagement = characterization.engagement
    assert characterization.source_operation_count == EXPECTED_SOURCE_OPERATIONS
    assert len(characterization.snapshot) == EXPECTED_SOURCE_OPERATIONS
    assert tuple(operation.ordinal for operation in characterization.snapshot) == tuple(range(EXPECTED_SOURCE_OPERATIONS))
    assert not set(engagement.tea_audited) & set(engagement.excluded)
    assert tuple(sorted((*engagement.tea_audited, *engagement.excluded))) == tuple(range(EXPECTED_SOURCE_OPERATIONS))
    disposition_sets = (
        set(engagement.certified),
        set(engagement.demonstrated_exceeded),
        set(engagement.unresolved),
    )
    assert not disposition_sets[0] & disposition_sets[1]
    assert not disposition_sets[0] & disposition_sets[2]
    assert not disposition_sets[1] & disposition_sets[2]
    assert set().union(*disposition_sets) == set(engagement.tea_audited)
    assert engagement.tea_audited_count == len(engagement.tea_audited)
    assert engagement.excluded_count == len(engagement.excluded)
    assert engagement.certified_count == len(engagement.certified)
    assert engagement.demonstrated_exceeded_count == len(engagement.demonstrated_exceeded)
    assert engagement.unresolved_count == len(engagement.unresolved)
    assert engagement.tea_audited_count + engagement.excluded_count == EXPECTED_SOURCE_OPERATIONS
    assert engagement.certified_count + engagement.demonstrated_exceeded_count + engagement.unresolved_count == engagement.tea_audited_count
    assert characterization.tea_audited_operation_count == engagement.tea_audited_count
    assert characterization.excluded_operation_count == engagement.excluded_count
    assert 0 <= characterization.sampled_material_contact_operations <= engagement.tea_audited_count
    for witness in characterization.witnesses:
        assert 0 <= witness.operation_index < EXPECTED_SOURCE_OPERATIONS
        assert witness.operation_index in engagement.demonstrated_exceeded

    assessment = characterization.assessment
    criteria = (
        assessment.uncut_fraction,
        assessment.gouging_motions,
        assessment.unsafe_rapids,
        assessment.continuity_breaks,
        assessment.zero_length_motions,
        assessment.degenerate_loops,
        assessment.redundant_operations,
        assessment.cap_exceedances,
        assessment.slotting_motions,
        assessment.max_engagement_step,
        assessment.max_loop_radius_step,
        assessment.tangent_breaks,
    )
    assert len(criteria) == 12
    assert tuple(criterion.name for criterion in criteria) == CRITERION_NAMES
    assert len({criterion.name for criterion in criteria}) == 12
    assert EVIDENCE_BY_CRITERION == EXPECTED_EVIDENCE_BY_CRITERION
    assert tuple(criterion.evidence for criterion in criteria) == tuple(EVIDENCE_BY_CRITERION[name] for name in CRITERION_NAMES)
    assert all(criterion.outcome in ALLOWED_OUTCOMES for criterion in criteria)

    quality = characterization.path_quality
    path_quality_projections = (
        quality.elementary.uncut_fraction,
        quality.elementary.gouging_motions,
        quality.elementary.unsafe_rapids,
        quality.elementary.continuity_breaks,
        quality.elementary.zero_length_motions,
        quality.elementary.degenerate_loops,
        quality.elementary.redundant_operations,
        quality.cut.cap_exceedances,
        quality.cut.slotting_motions,
        quality.cut.max_engagement_step_deg,
        quality.cut.max_loop_radius_step,
        quality.speed.tangent_breaks,
    )
    assert tuple(criterion.measured for criterion in criteria) == path_quality_projections

    markdown = REPORT_PATH.read_text(encoding="utf-8")
    assert markdown.startswith("# Held Figure 5 2D path characterization\n\n**HISTORICAL 2D BENCHMARK CHARACTERIZATION - NOT A MANUFACTURING RELEASE**")
    assert f"| Pixi invocation | {LIVE_COMMAND} |" in markdown
    assert "| Generator policy | benchmarks.runner.generate_toolpath |" in markdown
    assert f"| Case | {characterization.case_name} |" in markdown
    assert f"| Reference primitives | {characterization.reference_primitive_count} |" in markdown
    assert f"| Projection vertices | {characterization.projection_vertex_count} |" in markdown
    for normalized_unit_label in (
        "| Tool diameter (benchmark-normalized mm) |",
        "| Operation index | World X (benchmark-normalized mm) | World Y (benchmark-normalized mm) |",
        "| cut | engagement_variance_deg2 | degrees^2 |",
        "| cut | max_engagement_gradient_deg_per_length | degrees per benchmark-normalized mm |",
        "| cut | radial_depth_variance | benchmark-normalized mm^2 |",
        "| speed | max_curvature | inverse benchmark-normalized mm |",
    ):
        assert normalized_unit_label in markdown
    assert "Sampled negative observations are bounded evidence, not global proofs." in markdown
    assert "The exercised guarded replay is not evidence of native audit compatibility." in markdown
    for label, value in (
        ("Total source operations", characterization.source_operation_count),
        ("TEA-audited lateral operations", characterization.tea_audited_operation_count),
        ("TEA-audit-excluded operations", characterization.excluded_operation_count),
        ("Sampled material-contact operations", characterization.sampled_material_contact_operations),
        ("Generation (seconds)", characterization.generation_seconds),
        ("Guarded audit (seconds)", characterization.audit_seconds),
        ("Survey (seconds)", characterization.survey_seconds),
        ("Coverage plus reduction (seconds)", characterization.reduction_seconds),
        ("certified", engagement.certified_count),
        ("demonstrated_exceeded", engagement.demonstrated_exceeded_count),
        ("unresolved", engagement.unresolved_count),
    ):
        assert f"| {label} | {_rendered_number(value)} |" in markdown

    witness_cursor = 0
    for witness in characterization.witnesses:
        row = f"| {witness.operation_index} | {_rendered_number(witness.position.x)} | {_rendered_number(witness.position.y)} |"
        witness_cursor = markdown.index(row, witness_cursor) + len(row)
    for criterion in criteria:
        row = f"| {criterion.name} | {_rendered_number(criterion.measured)} | {_rendered_number(criterion.required)} | {criterion.evidence} | {criterion.outcome} |"
        assert row in markdown

    report_generated_at = _report_utc(markdown)
    assert report_generated_at >= run_started_at
    failures = post_qualification_failures(characterization)
    if failures:
        verdict = "not eligible for postprocessor qualification: " + "; ".join(failures)
        assert "**not eligible for postprocessor qualification**" in markdown
        failure_cursor = 0
        for failure in failures:
            failure_cursor = markdown.index(failure, failure_cursor) + len(failure)
        with pytest.raises(HeldPathNotEligibleForPostQualificationError):
            require_post_qualification_candidate(characterization)
    else:
        verdict = "eligible for postprocessor qualification evidence entry"
        assert "**eligible for postprocessor qualification evidence entry**" in markdown
        candidate = require_post_qualification_candidate(characterization)
        assert candidate.characterization is characterization
        assert candidate.snapshot is characterization.snapshot
    assert "This verdict is not a manufacturing release." in markdown
    assert EXIT_ZERO_BOUNDARY in markdown

    _write_live_run_ledger(
        run_started_at=run_started_at,
        outcome="completed",
        observed_phases=tuple(observed_phases),
        elapsed_seconds=time.monotonic() - operator_started,
        prior_report_generated_at=prior_report_generated_at,
        report_generated_at=report_generated_at,
        command_result="exit 0; evidence completion only",
        qualification_verdict=verdict,
        qualification_failures=failures,
    )
    print(f"REPORT: {REPORT_PATH}")
    print("CHARACTERIZATION COMPLETED")
    print(verdict)
    print(EXIT_ZERO_BOUNDARY)
