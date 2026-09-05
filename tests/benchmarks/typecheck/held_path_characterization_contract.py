from collections.abc import Callable
from typing import Optional

from typing_extensions import assert_type

from benchmarks.held_path_snapshot import HeldArcSnapshot
from benchmarks.held_path_snapshot import HeldCircleSnapshot
from benchmarks.held_path_snapshot import HeldLineSnapshot
from benchmarks.held_path_snapshot import HeldOperationSnapshot
from benchmarks.held_path_snapshot import assert_toolpath_matches_snapshot
from benchmarks.held_path_snapshot import snapshot_toolpath
from benchmarks.held_path_evidence import EngagementDispositionCounts
from benchmarks.held_path_evidence import EngagementExceedanceWitness
from benchmarks.held_path_evidence import HeldFigure5Characterization
from benchmarks.held_consumer_adapters import audit_figure5_engagement
from benchmarks.held_path_characterize import CharacterizationPhase
from benchmarks.held_path_characterize import Clock
from benchmarks.held_path_characterize import EngagementAuditor
from benchmarks.held_path_characterize import Generator
from benchmarks.held_path_characterize import PathSurveyor
from benchmarks.held_path_characterize import QualityEvidenceReducer
from benchmarks.held_path_characterize import _monotonic_seconds
from benchmarks.held_path_characterize import characterize_figure5
from benchmarks.held_post_qualification import HeldPostQualificationCandidate
from benchmarks.held_post_qualification import post_qualification_failures
from benchmarks.held_post_qualification import require_post_qualification_candidate
from benchmarks.survey import EngagementSample
from benchmarks.survey import MotionQuality
from benchmarks.units import Degrees
from benchmarks.units import MotionCount
from benchmarks.units import OperationIndex
from benchmarks.units import Seconds
from benchmarks.units import SquareMillimetre
from benchmarks.units import ToolRadiusMultiple
from benchmarks.units import UnitFraction
from benchmarks.units import closed_unit_fraction
from benchmarks.units import degrees_value
from benchmarks.units import motion_count
from benchmarks.units import operation_index
from benchmarks.units import seconds_value
from benchmarks.units import tool_radius_multiple
from compas_cgal.adaptive.units import Direction3
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Point3
from compas_cgal.adaptive.units import Radian
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.adaptive.units import WorldXYZ
from compas_cgal.toolpath import ToolpathResult
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.quality import PathQuality
from benchmarks.quality_observations import PathQualityAssessment
from benchmarks.quality_observations import QualityEvidence
from benchmarks.quality_observations import reduce_quality_evidence
from benchmarks.runner import generate_toolpath
from benchmarks.spec import PocketSpec
from benchmarks.survey import survey_path
from benchmarks.survey import PathSurvey
from compas_cgal.engagement import EngagementReport

seconds = assert_type(seconds_value(1.0, name="audit"), Seconds)
degrees = assert_type(degrees_value(80.0, name="cap"), Degrees)
fraction = assert_type(closed_unit_fraction(0.25, name="uncut fraction"), UnitFraction)
count = assert_type(motion_count(2, name="gouging motions"), MotionCount)
multiple = assert_type(tool_radius_multiple(2.0, name="step length"), ToolRadiusMultiple)
index = assert_type(operation_index(3, operation_count=4), OperationIndex)

seconds = degrees  # type: ignore[assignment]
degrees = fraction  # type: ignore[assignment]
fraction = count  # type: ignore[assignment]
count = multiple  # type: ignore[assignment]
multiple = index  # type: ignore[assignment]
index = seconds  # type: ignore[assignment]

point = assert_type(Point2[WorldXY].build(1.0, 2.0), Point2[WorldXY])
assert_type(Point3[WorldXYZ].build(1.0, 2.0, 0.0), Point3[WorldXYZ])


def _survey_contract(sample: EngagementSample, motion: MotionQuality) -> None:
    assert_type(sample.position, Point2[WorldXY])
    assert_type(sample.cap_exceeded, bool)
    assert_type(motion.cap_exceeded, bool)


def _snapshot_contract(
    result: ToolpathResult,
    line: HeldLineSnapshot,
    arc: HeldArcSnapshot,
    circle: HeldCircleSnapshot,
) -> None:
    snapshots = assert_type(snapshot_toolpath(result), tuple[HeldOperationSnapshot, ...])
    assert_toolpath_matches_snapshot(result, snapshots)
    assert_type(line.ordinal, OperationIndex)
    assert_type(line.start, Point3[WorldXYZ])
    assert_type(line.start_tangent, Optional[Direction3[WorldXYZ]])
    assert_type(arc.centre, Point3[WorldXYZ])
    assert_type(arc.xaxis, Direction3[WorldXYZ])
    assert_type(arc.radius, Millimetre)
    assert_type(arc.start_angle, Radian)
    assert_type(circle.centre, Point3[WorldXYZ])


def _characterization_contract(
    case: HeldReferenceCase,
    snapshots: tuple[HeldOperationSnapshot, ...],
    audit: EngagementReport,
    survey: PathSurvey,
    quality: QualityEvidence,
) -> None:
    characterization = assert_type(
        HeldFigure5Characterization.build(
            case=case,
            snapshot=snapshots,
            audit=audit,
            survey=survey,
            quality=quality,
            generation_seconds=Seconds(1.0),
            audit_seconds=Seconds(2.0),
            survey_seconds=Seconds(3.0),
            reduction_seconds=Seconds(4.0),
        ),
        HeldFigure5Characterization,
    )
    assert_type(characterization.snapshot, tuple[HeldOperationSnapshot, ...])
    assert_type(characterization.reference_primitive_count, int)
    assert_type(characterization.projection_vertex_count, int)
    assert_type(characterization.engagement, EngagementDispositionCounts)
    assert_type(characterization.witnesses, tuple[EngagementExceedanceWitness, ...])
    assert_type(characterization.witnesses[0].operation_index, OperationIndex)
    assert_type(characterization.witnesses[0].position, Point2[WorldXY])
    assert_type(characterization.generation_seconds, Seconds)
    assert_type(characterization.coverage_cell_area, SquareMillimetre)
    assert_type(characterization.coverage_wall_scallop_height, Millimetre)
    assert_type(characterization.coverage_uncut_fraction, UnitFraction)
    assert_type(characterization.coverage_remaining_area, SquareMillimetre)
    assert_type(characterization.path_quality, PathQuality)
    assert_type(characterization.assessment, PathQualityAssessment)
    assert_type(post_qualification_failures(characterization), tuple[str, ...])
    candidate = assert_type(
        HeldPostQualificationCandidate.build(characterization),
        HeldPostQualificationCandidate,
    )
    assert_type(
        require_post_qualification_candidate(characterization),
        HeldPostQualificationCandidate,
    )
    assert_type(candidate.characterization, HeldFigure5Characterization)
    assert_type(candidate.snapshot, tuple[HeldOperationSnapshot, ...])


def _orchestration_contract(spec: PocketSpec, result: ToolpathResult) -> None:
    generation_phase: CharacterizationPhase = "generation"
    audit_phase: CharacterizationPhase = "guarded_audit"
    survey_phase: CharacterizationPhase = "survey"
    reduction_phase: CharacterizationPhase = "quality_reduction"
    invalid_phase: CharacterizationPhase = "render"  # type: ignore[assignment]
    assert_type(generation_phase, CharacterizationPhase)
    assert_type(audit_phase, CharacterizationPhase)
    assert_type(survey_phase, CharacterizationPhase)
    assert_type(reduction_phase, CharacterizationPhase)
    assert_type(invalid_phase, CharacterizationPhase)

    def observer(phase: CharacterizationPhase) -> None:
        pass

    typed_observer: Callable[[CharacterizationPhase], None] = observer
    clock: Clock = _monotonic_seconds
    generator: Generator = generate_toolpath
    auditor: EngagementAuditor = audit_figure5_engagement
    surveyor: PathSurveyor = survey_path
    reducer: QualityEvidenceReducer = reduce_quality_evidence

    assert_type(_monotonic_seconds(), Seconds)
    assert_type(clock(), Seconds)
    assert_type(audit_figure5_engagement(spec, result), EngagementReport)
    assert_type(
        characterize_figure5(
            generator,
            auditor,
            surveyor,
            reducer,
            phase_observer=typed_observer,
            clock=clock,
        ),
        HeldFigure5Characterization,
    )
