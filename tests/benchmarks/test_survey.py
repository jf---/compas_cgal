"""Canonical path-survey evidence retained at every engagement station."""

import math

import numpy as np
import pytest
from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Frame
from compas.geometry import Line
from compas.geometry import Polygon

from compas_cgal import _coverage_2

from benchmarks.errors import UnreplayableOperationError
from benchmarks.spec import PocketSpec
from benchmarks.survey import QUALITY_SAMPLES_PER_MOTION
from benchmarks.survey import MotionQuality
from benchmarks.survey import survey_path
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult

SPEC = PocketSpec.build(
    name="survey_evidence",
    family="analytic",
    polygon=Polygon([[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]]),
    tool_diameter=2.0,
    tea_cap_deg=120.0,
)


def _survey_motion(geometry: object) -> MotionQuality:
    operations = [
        ToolpathOperation(
            geometry=Line([0.0, 0.0, 4.0], [0.0, 0.0, 0.0]),
            operation=OperationType.PLUNGE,
            path_index=0,
        ),
        ToolpathOperation(geometry=geometry, operation=OperationType.CUT, path_index=0),
    ]
    result = ToolpathResult(operations=operations, polyline=np.zeros((0, 3), dtype=float))
    return survey_path(SPEC, result).motions[0]


def test_survey_rejects_tilted_cut_circle() -> None:
    tilted = Circle(2.0, frame=Frame([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]))

    with pytest.raises(UnreplayableOperationError):
        _survey_motion(tilted)


def test_survey_rejects_curve_outside_the_inferred_cut_plane() -> None:
    displaced = Arc(radius=2.0, start_angle=0.0, end_angle=0.5 * math.pi, frame=Frame([0.0, 0.0, 1.0]))

    with pytest.raises(UnreplayableOperationError):
        _survey_motion(displaced)


def test_survey_rejects_curve_below_the_line_inferred_cut_plane() -> None:
    displaced = Arc(radius=2.0, start_angle=0.0, end_angle=0.5 * math.pi, frame=Frame([0.0, 0.0, -1.0]))

    with pytest.raises(UnreplayableOperationError):
        _survey_motion(displaced)


def test_survey_accepts_a_coherent_circle_only_cut_plane() -> None:
    operations = [ToolpathOperation(geometry=Circle(radius, frame=Frame([0.0, 0.0, 2.0])), operation=OperationType.CUT, path_index=0) for radius in (2.0, 3.0)]
    result = ToolpathResult(operations=operations, polyline=np.zeros((0, 3), dtype=float))

    assert len(survey_path(SPEC, result).motions) == 2


def test_survey_rejects_mixed_planes_in_a_curve_only_stream() -> None:
    operations = [ToolpathOperation(geometry=Circle(2.0, frame=Frame([0.0, 0.0, z])), operation=OperationType.CUT, path_index=0) for z in (0.0, 1.0)]
    result = ToolpathResult(operations=operations, polyline=np.zeros((0, 3), dtype=float))

    with pytest.raises(UnreplayableOperationError):
        survey_path(SPEC, result)


def test_survey_accepts_cut_circle_with_horizontal_clearance_link() -> None:
    operations = [
        ToolpathOperation(geometry=Circle(2.0, frame=Frame([0.0, 0.0, 0.0])), operation=OperationType.CUT, path_index=0),
        ToolpathOperation(geometry=Line([2.0, 0.0, 2.0], [-2.0, 0.0, 2.0]), operation=OperationType.LINK, path_index=0),
    ]
    result = ToolpathResult(operations=operations, polyline=np.zeros((0, 3), dtype=float))

    survey = survey_path(SPEC, result)

    assert len(survey.motions) == 1
    assert len(survey.rapids) == 1


def test_survey_accepts_cut_circle_with_curved_clearance_link() -> None:
    operations = [
        ToolpathOperation(geometry=Circle(2.0, frame=Frame([0.0, 0.0, 0.0])), operation=OperationType.CUT, path_index=0),
        ToolpathOperation(geometry=Circle(3.0, frame=Frame([0.0, 0.0, 2.0])), operation=OperationType.LINK, path_index=0),
    ]
    result = ToolpathResult(operations=operations, polyline=np.zeros((0, 3), dtype=float))

    survey = survey_path(SPEC, result)

    assert len(survey.motions) == 1
    assert len(survey.rapids) == 1


def test_survey_preserves_retract_precedence_for_tilted_curve() -> None:
    operations = [
        ToolpathOperation(geometry=Circle(2.0, frame=Frame([0.0, 0.0, 0.0])), operation=OperationType.CUT, path_index=0),
        ToolpathOperation(
            geometry=Circle(1.0, frame=Frame([0.0, 0.0, 2.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0])),
            operation=OperationType.RETRACT,
            path_index=0,
        ),
    ]
    result = ToolpathResult(operations=operations, polyline=np.zeros((0, 3), dtype=float))

    survey = survey_path(SPEC, result)

    assert survey.retract_indices == (1,)
    assert len(survey.rapids) == 1


def test_standard_circle_samples_retain_one_typed_seam() -> None:
    motion = _survey_motion(Circle(2.0, frame=Frame([0.0, 0.0, 0.0])))

    assert len(motion.samples) == QUALITY_SAMPLES_PER_MOTION == 45
    assert isinstance(motion.samples[0].position, Point2)
    assert motion.samples[0].position == Point2[WorldXY].build(2.0, 0.0)
    assert motion.samples[-1].position.x == pytest.approx(1.9805361374831407)
    assert motion.samples[-1].position.y == pytest.approx(-0.27834620192013124)
    assert motion.samples[-1].position != motion.samples[0].position
    assert motion.cap_exceeded is any(sample.cap_exceeded for sample in motion.samples)


@pytest.mark.parametrize(
    ("geometry", "expected_start", "expected_end"),
    [
        (Line([0.0, 0.0, 0.0], [3.0, 1.0, 0.0]), Point2[WorldXY].build(0.0, 0.0), Point2[WorldXY].build(3.0, 1.0)),
        (
            Arc(radius=2.0, start_angle=0.0, end_angle=0.5 * math.pi, frame=Frame([0.0, 0.0, 0.0])),
            Point2[WorldXY].build(2.0, 0.0),
            Point2[WorldXY].build(0.0, 2.0),
        ),
    ],
    ids=("line", "arc"),
)
def test_standard_open_motion_samples_retain_both_typed_endpoints(
    geometry: object,
    expected_start: Point2[WorldXY],
    expected_end: Point2[WorldXY],
) -> None:
    motion = _survey_motion(geometry)

    assert len(motion.samples) == QUALITY_SAMPLES_PER_MOTION + 1 == 46
    assert isinstance(motion.samples[0].position, Point2)
    assert motion.samples[0].position.x == pytest.approx(expected_start.x)
    assert motion.samples[0].position.y == pytest.approx(expected_start.y)
    assert isinstance(motion.samples[-1].position, Point2)
    assert motion.samples[-1].position.x == pytest.approx(expected_end.x)
    assert motion.samples[-1].position.y == pytest.approx(expected_end.y)
    assert motion.cap_exceeded is any(sample.cap_exceeded for sample in motion.samples)


def test_survey_uses_lightweight_exact_cutter_centre_domain(monkeypatch: pytest.MonkeyPatch) -> None:
    native_domain = _coverage_2.CutterCentreDomain2
    constructions: list[tuple[np.ndarray, list[np.ndarray], float]] = []
    membership_queries: list[tuple[float, float]] = []

    class RecordingCutterCentreDomain:
        def __init__(self, delegate: object) -> None:
            self._delegate = delegate

        @classmethod
        def build(
            cls,
            boundary: np.ndarray,
            holes: list[np.ndarray],
            tool_radius: float,
        ) -> "RecordingCutterCentreDomain":
            constructions.append((boundary.copy(), [hole.copy() for hole in holes], tool_radius))
            return cls(native_domain.build(boundary, holes, tool_radius))

        def contains(self, x: float, y: float) -> bool:
            membership_queries.append((x, y))
            return self._delegate.contains(x, y)

    def reject_eager_domain(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("survey constructed the eager reachable domain")

    monkeypatch.setattr("benchmarks.survey._coverage_2.CutterCentreDomain2", RecordingCutterCentreDomain)
    monkeypatch.setattr("benchmarks.survey._coverage_2.ReachableDomain2", reject_eager_domain)

    motion = _survey_motion(Line([-6.0, 0.0, 0.0], [-5.0, 0.0, 0.0]))

    assert len(constructions) == 1
    boundary, holes, tool_radius = constructions[0]
    np.testing.assert_array_equal(
        boundary,
        np.asarray(((-6.0, -4.0, 0.0), (6.0, -4.0, 0.0), (6.0, 4.0, 0.0), (-6.0, 4.0, 0.0))),
    )
    assert holes == []
    assert tool_radius == 1.0
    assert len(membership_queries) == QUALITY_SAMPLES_PER_MOTION + 1
    assert membership_queries[0] == (-6.0, 0.0)
    assert membership_queries[-1] == (-5.0, 0.0)
    assert not motion.samples[0].inside_centre_domain
    assert motion.samples[-1].inside_centre_domain


def test_cap_verdict_is_not_reconstructed_from_reporting_degrees(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("benchmarks.survey.math.degrees", lambda _radians: 0.0)

    motion = _survey_motion(Line([0.0, 0.0, 0.0], [3.0, 0.0, 0.0]))

    assert all(sample.engagement_deg == 0.0 for sample in motion.samples)
    assert any(sample.cap_exceeded for sample in motion.samples)
    assert motion.cap_exceeded is any(sample.cap_exceeded for sample in motion.samples)
