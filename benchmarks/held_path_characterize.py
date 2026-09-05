"""Generate Figure 5 once and reduce its immutable evidence transactionally."""

from __future__ import annotations

import time
from collections.abc import Callable
from typing import Literal
from typing import Protocol

from typing_extensions import TypeAlias

from benchmarks.held_path_evidence import HeldFigure5Characterization
from benchmarks.held_path_snapshot import HeldOperationSnapshot
from benchmarks.held_path_snapshot import assert_toolpath_matches_snapshot
from benchmarks.held_path_snapshot import snapshot_toolpath
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.quality_observations import QualityEvidence
from benchmarks.spec import PocketSpec
from benchmarks.survey import PathSurvey
from benchmarks.units import Seconds
from benchmarks.units import seconds_value
from compas_cgal.engagement import EngagementReport
from compas_cgal.toolpath import ToolpathResult

CharacterizationPhase: TypeAlias = Literal[
    "generation",
    "guarded_audit",
    "survey",
    "quality_reduction",
]
Clock: TypeAlias = Callable[[], Seconds]


class Generator(Protocol):
    def __call__(self, spec: PocketSpec) -> ToolpathResult: ...


class EngagementAuditor(Protocol):
    def __call__(self, spec: PocketSpec, result: ToolpathResult) -> EngagementReport: ...


class PathSurveyor(Protocol):
    def __call__(self, spec: PocketSpec, result: ToolpathResult) -> PathSurvey: ...


class QualityEvidenceReducer(Protocol):
    def __call__(
        self,
        spec: PocketSpec,
        snapshot: tuple[HeldOperationSnapshot, ...],
        survey: PathSurvey,
    ) -> QualityEvidence: ...


def _monotonic_seconds() -> Seconds:
    return seconds_value(time.perf_counter(), name="monotonic clock")


def characterize_figure5(
    generator: Generator,
    engagement_auditor: EngagementAuditor,
    path_surveyor: PathSurveyor,
    quality_evidence_reducer: QualityEvidenceReducer,
    *,
    phase_observer: Callable[[CharacterizationPhase], None],
    clock: Clock = _monotonic_seconds,
) -> HeldFigure5Characterization:
    """Build one immutable characterization from one generated operation stream."""
    case = load_held_reference_case("figure5")
    spec = case.pocket_spec()

    phase_observer("generation")
    started = clock()
    result = generator(spec)
    generation_seconds = seconds_value(clock() - started, name="generation duration")

    snapshot = snapshot_toolpath(result)
    assert_toolpath_matches_snapshot(result, snapshot)

    phase_observer("guarded_audit")
    started = clock()
    audit = engagement_auditor(spec, result)
    audit_seconds = seconds_value(clock() - started, name="guarded audit duration")
    assert_toolpath_matches_snapshot(result, snapshot)

    phase_observer("survey")
    started = clock()
    survey = path_surveyor(spec, result)
    survey_seconds = seconds_value(clock() - started, name="survey duration")
    assert_toolpath_matches_snapshot(result, snapshot)

    phase_observer("quality_reduction")
    started = clock()
    quality = quality_evidence_reducer(spec, snapshot, survey)
    reduction_seconds = seconds_value(clock() - started, name="quality reduction duration")
    assert_toolpath_matches_snapshot(result, snapshot)

    return HeldFigure5Characterization.build(
        case=case,
        snapshot=snapshot,
        audit=audit,
        survey=survey,
        quality=quality,
        generation_seconds=generation_seconds,
        audit_seconds=audit_seconds,
        survey_seconds=survey_seconds,
        reduction_seconds=reduction_seconds,
    )
