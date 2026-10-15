"""Measure the twelve gate cases into a `GateVectors` document.

The criterion table below is the single definition of what the results file
carries: each entry names its group, whether the gate decides on it, the
precision it is stored at, the gate's threshold, and where its value is read
from the canonical quality reduction. The gate test
(`tests.benchmarks.test_quality`) and this measurement share that reduction, so
the file and the test cannot disagree about what a criterion means.

Run as ``pixi run gate-vectors``: twelve cases in parallel, bounded by the
slowest (radius_regulated on L_shape at the 40-degree cap).
"""

from __future__ import annotations

import argparse
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Callable
from typing import Sequence

from benchmarks.coverage import COVERAGE_GRID_SAMPLES
from benchmarks.gate import GATE_CAP_CASES
from benchmarks.gate import GATE_GENERATOR_NAMES
from benchmarks.gate import GATE_GENERATORS
from benchmarks.gate import GATE_POCKET_NAMES
from benchmarks.gate import GATE_TOOL_DIAMETER
from benchmarks.gate import gate_pocket
from benchmarks.held_path_snapshot import snapshot_toolpath
from benchmarks.quality import ENGAGEMENT_BANDS_DEG
from benchmarks.quality_observations import QualityEvidence
from benchmarks.quality_observations import reduce_quality_evidence
from benchmarks.survey import QUALITY_SAMPLES_PER_MOTION
from benchmarks.survey import PathSurvey
from benchmarks.survey import survey_path
from tools.gate_vectors import CASE_KEY_SEPARATOR
from tools.gate_vectors import GATE_TABLE_PATH
from tools.gate_vectors import GATE_VECTORS_PATH
from tools.gate_vectors import REQUIRED_CAP
from tools.gate_vectors import REQUIRED_COMPONENTS
from tools.gate_vectors import RULES_VERSION
from tools.gate_vectors import CriterionSpec
from tools.gate_vectors import Direction
from tools.gate_vectors import GateVectors
from tools.gate_vectors import RebaselineRecord
from tools.gate_vectors import load_gate_vectors

Reader = Callable[[QualityEvidence, PathSurvey], "float | int"]


@dataclass(frozen=True)
class MeasuredCriterion:
    """One criterion's storage specification and its reader.

    Attributes:
        spec: How it is stored and compared.
        read: Its value from one reduced case; None for a wall-clock cell, which
            the measurement fills itself.
    """

    spec: CriterionSpec
    read: Reader | None


def _criterion(name: str, group: str, direction: Direction, decimals: int, required: str | None, read: Reader | None) -> MeasuredCriterion:
    return MeasuredCriterion(spec=CriterionSpec(name=name, group=group, direction=direction, decimals=decimals, required=required), read=read)


def _band_name(low: float, high: float) -> str:
    return f"engagement_length_{low:03.0f}_{high:03.0f}_deg"


def _band_reader(position: int) -> Reader:
    return lambda evidence, _survey: evidence.path_quality.longevity.engagement_length_histogram[position][2]


# Counts are stored with no decimals; ratios and fractions at four or six; the
# precisions of the gated criteria are the gate's own report precisions
# (`tests.benchmarks.test_quality._violations`).
CRITERIA: tuple[MeasuredCriterion, ...] = (
    _criterion("uncut_fraction", "elementary", "lower", 6, "0", lambda e, _s: float(e.assessment.uncut_fraction.measured)),
    _criterion("gouging_motions", "elementary", "lower", 0, "0", lambda e, _s: int(e.assessment.gouging_motions.measured)),
    _criterion("unsafe_rapids", "elementary", "lower", 0, "0", lambda e, _s: int(e.assessment.unsafe_rapids.measured)),
    _criterion("continuity_breaks", "elementary", "lower", 0, "0", lambda e, _s: int(e.assessment.continuity_breaks.measured)),
    _criterion("zero_length_motions", "elementary", "lower", 0, "0", lambda e, _s: int(e.assessment.zero_length_motions.measured)),
    _criterion("redundant_operations", "elementary", "lower", 0, "0", lambda e, _s: int(e.assessment.redundant_operations.measured)),
    _criterion("material_plunges", "elementary", "lower", 0, REQUIRED_COMPONENTS, lambda e, _s: int(e.assessment.material_plunges.measured)),
    _criterion("misplaced_plunges", "elementary", "lower", 0, "0", lambda e, _s: int(e.assessment.misplaced_plunges.measured)),
    _criterion("degenerate_loops", "elementary", "reported", 0, None, lambda e, _s: e.path_quality.elementary.degenerate_loops),
    _criterion("marginal_loops", "elementary", "reported", 0, None, lambda e, _s: e.path_quality.elementary.marginal_loops),
    _criterion("recut_fraction", "elementary", "reported", 4, None, lambda e, _s: e.path_quality.elementary.recut_fraction),
    _criterion("within_band_motions", "elementary", "reported", 0, None, lambda _e, s: s.within_band_motions),
    _criterion("centre_domain_components", "elementary", "reported", 0, None, lambda e, _s: e.coverage.centre_domain_components),
    _criterion("cap_exceedances", "cut", "lower", 0, "0", lambda e, _s: int(e.assessment.cap_exceedances.measured)),
    _criterion("slotting_motions", "cut", "lower", 0, "0", lambda e, _s: int(e.assessment.slotting_motions.measured)),
    _criterion("max_engagement_step_deg", "cut", "lower", 2, REQUIRED_CAP, lambda e, _s: float(e.assessment.max_engagement_step.measured)),
    _criterion("max_loop_radius_step", "cut", "lower", 3, "2", lambda e, _s: float(e.assessment.max_loop_radius_step.measured)),
    _criterion("max_engagement_deg", "cut", "reported", 2, None, lambda e, _s: e.path_quality.cut.max_engagement_deg),
    _criterion("engagement_p95_deg", "cut", "reported", 2, None, lambda e, _s: e.path_quality.cut.engagement_p95_deg),
    _criterion("engagement_variance_deg2", "cut", "reported", 1, None, lambda e, _s: e.path_quality.cut.engagement_variance_deg2),
    _criterion("max_chip_thickness_ratio", "cut", "reported", 4, None, lambda e, _s: e.path_quality.cut.max_chip_thickness_ratio),
    _criterion("low_chip_thickness_ratio", "cut", "reported", 4, None, lambda e, _s: e.path_quality.cut.low_chip_thickness_ratio),
    _criterion("max_engagement_gradient_deg_per_length", "cut", "reported", 1, None, lambda e, _s: e.path_quality.cut.max_engagement_gradient_deg_per_length),
    _criterion("immersion_steady_fraction", "cut", "reported", 4, None, lambda e, _s: e.path_quality.cut.immersion_steady_fraction),
    _criterion("immersion_at_design_fraction", "cut", "reported", 4, None, lambda e, _s: e.path_quality.cut.immersion_at_design_fraction),
    _criterion("immersion_excursions", "cut", "reported", 0, None, lambda e, _s: e.path_quality.cut.immersion_excursions),
    _criterion("mean_radial_depth", "cut", "reported", 4, None, lambda e, _s: e.path_quality.cut.mean_radial_depth),
    _criterion("radial_depth_variance", "cut", "reported", 4, None, lambda e, _s: e.path_quality.cut.radial_depth_variance),
    _criterion("wall_scallop_height", "cut", "reported", 4, None, lambda e, _s: e.path_quality.cut.wall_scallop_height),
    _criterion("loop_radius_cv", "cut", "reported", 4, None, lambda e, _s: e.path_quality.cut.loop_radius_cv),
    _criterion("tangent_breaks", "speed", "lower", 0, "0", lambda e, _s: int(e.assessment.tangent_breaks.measured)),
    _criterion("path_length", "speed", "reported", 3, None, lambda e, _s: e.path_quality.path_length),
    _criterion("cutting_length", "speed", "reported", 3, None, lambda e, _s: e.path_quality.speed.cutting_length),
    _criterion("air_length", "speed", "reported", 3, None, lambda e, _s: e.path_quality.speed.air_length),
    _criterion("air_at_depth_length", "speed", "reported", 3, None, lambda e, _s: e.path_quality.speed.air_at_depth_length),
    _criterion("air_fraction", "speed", "reported", 4, None, lambda e, _s: e.path_quality.speed.air_fraction),
    _criterion("max_curvature", "speed", "reported", 3, None, lambda e, _s: e.path_quality.speed.max_curvature),
    _criterion("curvature_breaks", "speed", "reported", 0, None, lambda e, _s: e.path_quality.speed.curvature_breaks),
    _criterion("direction_reversals", "speed", "reported", 0, None, lambda e, _s: e.path_quality.speed.direction_reversals),
    _criterion("retract_count", "speed", "reported", 0, None, lambda e, _s: e.path_quality.speed.retract_count),
    _criterion("reentry_count", "speed", "reported", 0, None, lambda e, _s: e.path_quality.speed.reentry_count),
    _criterion("material_entries", "longevity", "reported", 0, None, lambda e, _s: e.path_quality.longevity.material_entries),
    _criterion("cut_air_alternations", "longevity", "reported", 0, None, lambda e, _s: e.path_quality.longevity.cut_air_alternations),
    _criterion("alternations_per_length", "longevity", "reported", 5, None, lambda e, _s: e.path_quality.longevity.alternations_per_length),
    *(
        _criterion(_band_name(low, high), "longevity", "reported", 3, None, _band_reader(position))
        for position, (low, high) in enumerate(zip(ENGAGEMENT_BANDS_DEG, ENGAGEMENT_BANDS_DEG[1:]))
    ),
    _criterion("generation_s", "cost", "wallclock", 2, None, None),
    _criterion("gate_s", "cost", "wallclock", 2, None, None),
)

GateCase = tuple[str, float, str, str]


def gate_cases() -> tuple[GateCase, ...]:
    """Every gate case as ``(cap id, cap degrees, generator, pocket)``, in table order."""
    return tuple((cap_id, float(cap), generator, pocket) for cap_id, cap in GATE_CAP_CASES for generator in GATE_GENERATOR_NAMES for pocket in GATE_POCKET_NAMES)


def case_key(cap_id: str, generator: str, pocket: str) -> str:
    """The results-file key of one case."""
    return CASE_KEY_SEPARATOR.join((cap_id, generator, pocket))


def measure_case(case: GateCase) -> tuple[str, dict[str, str]]:
    """Generate, survey and reduce one case into its stored cells.

    Args:
        case: ``(cap id, cap degrees, generator, pocket)``.

    Returns:
        The case key and its cells, in criterion order.
    """
    cap_id, cap, generator, pocket = case
    spec = gate_pocket(pocket, tea_cap_deg=cap)  # type: ignore[arg-type]
    started = time.perf_counter()
    result = GATE_GENERATORS[generator](spec)
    generated = time.perf_counter()
    survey = survey_path(spec, result)
    evidence = reduce_quality_evidence(spec, snapshot_toolpath(result), survey)
    finished = time.perf_counter()
    wallclock = {"generation_s": generated - started, "gate_s": finished - generated}
    cells: dict[str, str] = {}
    for criterion in CRITERIA:
        value = wallclock[criterion.spec.name] if criterion.read is None else criterion.read(evidence, survey)
        cells[criterion.spec.name] = criterion.spec.format(value)
    return case_key(cap_id, generator, pocket), cells


class MissingRebaselineApprovalError(Exception):
    """Re-measuring under a new rules version without a cited re-baseline approval."""


class MissingRebaselineCriteriaError(Exception):
    """Re-measuring under a new rules version without naming a criterion the bump changed or retired."""


def measure_gate_vectors(cases: Sequence[GateCase] | None = None, *, rebaseline: RebaselineRecord | None = None) -> GateVectors:
    """Measure *cases* (default: every gate case) in parallel.

    Args:
        cases: The cases to measure.
        rebaseline: The re-baseline record the written document carries; see
            `rebaseline_for`. None for a measurement that is only compared.

    Returns:
        The measured document at the current `RULES_VERSION`.
    """
    selected = gate_cases() if cases is None else tuple(cases)
    with ProcessPoolExecutor(max_workers=len(selected)) as pool:
        measured = dict(pool.map(measure_case, selected))
    return GateVectors(
        rules_version=RULES_VERSION,
        instrument={"tool_diameter": GATE_TOOL_DIAMETER, "samples_per_motion": QUALITY_SAMPLES_PER_MOTION, "coverage_grid": COVERAGE_GRID_SAMPLES},
        criteria=tuple(criterion.spec for criterion in CRITERIA),
        cases={case_key(cap_id, generator, pocket): measured[case_key(cap_id, generator, pocket)] for cap_id, _cap, generator, pocket in selected},
        rebaseline=rebaseline,
    )


def rebaseline_for(existing: GateVectors | None, approval: str | None, criteria: Sequence[str], retired: Sequence[str]) -> RebaselineRecord | None:
    """The re-baseline record a re-measurement must carry, never dropping one.

    Under the same rules version the existing record is kept as it is. Under a
    new one, the record names the old version, the approval (which must cite a
    ledger ruling), the gated criteria whose definition changed and the
    criteria the old version gated that the new one no longer does; the lists
    are validated against the measured criteria by `GateVectors.build` and
    against the parent's by the ratchet.

    Args:
        existing: The vectors currently in the file, or None for the first.
        approval: The user-approved ledger line, for a rules-version bump.
        criteria: The gated criteria the bump changed, for a rules-version bump.
        retired: The criteria the bump stops gating, for a rules-version bump.

    Returns:
        The record to write, or None when there is none.

    Raises:
        MissingRebaselineApprovalError: The rules version changed and no
            approval was given.
        MissingRebaselineCriteriaError: The rules version changed and no
            changed or retired criterion was named.
    """
    if existing is None:
        return None
    if existing.rules_version == RULES_VERSION:
        return existing.rebaseline
    if approval is None:
        raise MissingRebaselineApprovalError(
            f"rules {existing.rules_version} -> {RULES_VERSION} needs a recorded re-baseline; pass --rebaseline-approval '<Ruling id: user-approved ledger line>'."
        )
    if not criteria and not retired:
        raise MissingRebaselineCriteriaError(
            f"rules {existing.rules_version} -> {RULES_VERSION} re-baselines only the criteria it changed or retires; "
            "pass --rebaseline-criteria NAME [NAME ...] and/or --rebaseline-retired NAME [NAME ...]."
        )
    return RebaselineRecord(from_version=existing.rules_version, approval=approval, criteria=tuple(criteria), retired=tuple(retired))


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line entry: measure and write the vectors and their table.

    Args:
        argv: Arguments, defaulting to ``sys.argv[1:]``.

    Returns:
        The process exit status.
    """
    from tools.gate_table import render_gate_table

    parser = argparse.ArgumentParser(prog="python -m tools.gate_measure", description="Measure the gate cases into the pinned results file.")
    parser.add_argument("--write", type=Path, default=GATE_VECTORS_PATH)
    parser.add_argument("--table", type=Path, default=GATE_TABLE_PATH)
    parser.add_argument("--rebaseline-approval", default=None, help="the ledger line approving a rules-version bump")
    parser.add_argument("--rebaseline-criteria", nargs="+", default=(), metavar="NAME", help="the gated criteria whose definition a rules-version bump changed")
    parser.add_argument("--rebaseline-retired", nargs="+", default=(), metavar="NAME", help="the criteria the previous rules gated that a bump no longer gates")
    args = parser.parse_args(argv)
    existing = load_gate_vectors(args.write) if args.write.exists() else None
    record = rebaseline_for(existing, args.rebaseline_approval, args.rebaseline_criteria, args.rebaseline_retired)
    vectors = GateVectors.build(measure_gate_vectors(rebaseline=record).to_json())
    args.write.write_text(vectors.dumps())
    args.table.write_text(render_gate_table(vectors))
    print(f"wrote {args.write} and {args.table}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
