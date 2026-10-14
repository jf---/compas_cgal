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
from compas_cgal import _stock_2

from benchmarks.errors import UnreplayableOperationError
from benchmarks.spec import PocketSpec
from benchmarks.survey import QUALITY_SAMPLES_PER_MOTION
from benchmarks.survey import MotionQuality
from benchmarks.survey import survey_path
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY
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


# Exact surrogates for the in-process verdicts (boundary doctrine): 4*sin^2(cap/2).
HALF_TURN_RATIO = 4.0
JUST_UNDER_HALF_TURN_RATIO = 4.0 * math.sin(math.radians(179.0) / 2.0) ** 2
# degrees; the gate's angular sample spacing (QUALITY_SAMPLES_PER_MOTION stations per
# full turn, 8 deg). The acceptance's reporting allowance "180 deg +- one sample
# spacing"; the verdicts themselves are exact and asserted separately.
SLOT_SAMPLE_SPACING_DEG = 360.0 / QUALITY_SAMPLES_PER_MOTION


def test_a_straight_slot_into_solid_stock_reads_a_half_turn() -> None:
    # verify_direction.py read [0.0, 192.8, 205.7] ... [360.0, 360.0, 360.0] before Task 0.1.
    motion = _survey_motion(Line([0.0, 0.0, 0.0], [5.0, 0.0, 0.0]))

    assert motion.samples[0].engagement_deg == 0.0  # the plunge disk is already cleared
    for sample in motion.samples[1:]:
        assert abs(sample.engagement_deg - 180.0) <= SLOT_SAMPLE_SPACING_DEG
    assert motion.slot_exceeded  # slot threshold 60 deg < 180 deg, exact verdict
    assert motion.cap_exceeded  # cap 120 deg < 180 deg, exact verdict


def test_the_slot_half_turn_is_an_exact_cap_verdict() -> None:
    stock = _stock_2.Stock2(np.array([[-6.0, -4.0], [6.0, -4.0], [6.0, 4.0], [-6.0, 4.0]], dtype=float), [])
    stock.subtract_disk(0.0, 0.0, 1.0)
    parameters = [i / QUALITY_SAMPLES_PER_MOTION for i in range(1, QUALITY_SAMPLES_PER_MOTION + 1)]

    within = _stock_2.engagement_line_prefixes(stock, 0.0, 0.0, 5.0, 0.0, parameters, 1.0, HALF_TURN_RATIO, 0.0)
    over = _stock_2.engagement_line_prefixes(stock, 0.0, 0.0, 5.0, 0.0, parameters, 1.0, JUST_UNDER_HALF_TURN_RATIO, 0.0)

    assert not any(exceeded for _total, _run, exceeded in within)
    assert all(exceeded for _total, _run, exceeded in over)


def test_a_clockwise_loop_is_sampled_in_travel_order() -> None:
    loop = ToolpathOperation(geometry=Circle(2.0, frame=Frame([0.0, 0.0, 0.0])), operation=OperationType.CUT, path_index=0, clockwise=True)
    plunge = ToolpathOperation(geometry=Line([0.0, 0.0, 4.0], [0.0, 0.0, 0.0]), operation=OperationType.PLUNGE, path_index=0)
    motion = survey_path(SPEC, ToolpathResult(operations=[plunge, loop], polyline=np.zeros((0, 3), dtype=float))).motions[0]

    assert motion.samples[0].position == Point2[WorldXY].build(2.0, 0.0)
    assert motion.samples[1].position.y < 0.0  # clockwise: the first step goes below the seam
    assert motion.samples[-1].position.y > 0.0


# Plunge disks off every symmetry line of the loop below: a clockwise and a
# counter-clockwise prefix meet different material.
ASYMMETRIC_PLUNGES = ((0.0, 0.0), (-1.2, 2.4), (2.4, -1.6))
SURVEY_LOOP_RADIUS = 2.0
# radians; resolves the tangent direction in binary64 for the slow oracle only.
SLOW_ORACLE_APPROACH_ANGLE = 1e-6
# degrees; the slow sweep under-covers by at most 1e-4 r per its construction
# contract, which 0.1 degree comfortably isolates from a reversed prefix.
SLOW_ORACLE_REPORTING_TOLERANCE_DEG = 0.1


def _asymmetric_plunges() -> list[ToolpathOperation]:
    return [ToolpathOperation(geometry=Line([x, y, 4.0], [x, y, 0.0]), operation=OperationType.PLUNGE, path_index=0) for x, y in ASYMMETRIC_PLUNGES]


def _slow_oracle_deg(signed_start: float, signed_angle: float, clockwise: bool) -> float:
    """Static engagement at the station after sweeping the stock to just before it."""
    stock = _stock_2.Stock2(np.array([[-6.0, -4.0], [6.0, -4.0], [6.0, 4.0], [-6.0, 4.0]], dtype=float), [])
    for x, y in ASYMMETRIC_PLUNGES:
        stock.subtract_disk(x, y, 1.0)
    direction = -1.0 if clockwise else 1.0
    prior = signed_angle - direction * SLOW_ORACLE_APPROACH_ANGLE
    rho = SURVEY_LOOP_RADIUS
    stock.subtract_arc_sweep(0.0, 0.0, rho * math.cos(signed_start), rho * math.sin(signed_start), rho * math.cos(prior), rho * math.sin(prior), clockwise, 1.0)
    run = _stock_2.engagement_at(stock, rho * math.cos(signed_angle), rho * math.sin(signed_angle), 1.0, 4.0, 0.0)[1]
    return math.degrees(run)


@pytest.mark.parametrize("clockwise", [False, True])
def test_a_loop_reads_its_own_prefix_in_its_travel_direction(clockwise: bool) -> None:
    loop = ToolpathOperation(geometry=Circle(SURVEY_LOOP_RADIUS, frame=Frame([0.0, 0.0, 0.0])), operation=OperationType.CUT, path_index=0, clockwise=clockwise)
    motion = survey_path(SPEC, ToolpathResult(operations=[*_asymmetric_plunges(), loop], polyline=np.zeros((0, 3), dtype=float))).motions[0]

    direction = -1.0 if clockwise else 1.0
    for k, sample in enumerate(motion.samples[1:], start=1):
        angle = direction * math.tau * k / QUALITY_SAMPLES_PER_MOTION
        assert sample.position.x == pytest.approx(SURVEY_LOOP_RADIUS * math.cos(angle))
        assert sample.position.y == pytest.approx(SURVEY_LOOP_RADIUS * math.sin(angle))
        expected = _slow_oracle_deg(0.0, angle, clockwise)
        assert sample.engagement_deg == pytest.approx(expected, abs=SLOW_ORACLE_REPORTING_TOLERANCE_DEG), f"station {k}"


def _assert_arc_follows_its_travel(geometry: Arc, clockwise: bool, signed_start: float, travel_sweep: float) -> None:
    arc = ToolpathOperation(geometry=geometry, operation=OperationType.CUT, path_index=0, clockwise=clockwise)
    motion = survey_path(SPEC, ToolpathResult(operations=[*_asymmetric_plunges(), arc], polyline=np.zeros((0, 3), dtype=float))).motions[0]

    direction = -1.0 if clockwise else 1.0
    for k, sample in enumerate(motion.samples[1:], start=1):
        angle = signed_start + direction * travel_sweep * k / QUALITY_SAMPLES_PER_MOTION
        assert sample.position.x == pytest.approx(SURVEY_LOOP_RADIUS * math.cos(angle)), f"station {k}"
        assert sample.position.y == pytest.approx(SURVEY_LOOP_RADIUS * math.sin(angle)), f"station {k}"
        expected = _slow_oracle_deg(signed_start, angle, clockwise)
        assert sample.engagement_deg == pytest.approx(expected, abs=SLOW_ORACLE_REPORTING_TOLERANCE_DEG), f"station {k}"


@pytest.mark.parametrize(("clockwise", "start_angle", "end_angle"), [(False, 0.0, 1.5 * math.pi), (True, 1.5 * math.pi, 2.0 * math.pi)])
def test_an_arc_reads_its_own_prefix_in_its_travel_direction(clockwise: bool, start_angle: float, end_angle: float) -> None:
    # Three quarters of a turn in the travel direction. Arcs follow the snapshot
    # convention: start_angle <= end_angle, the flag names the turn direction from
    # point_at(0) to point_at(1). Clockwise from 270 degrees to 360 is the long way.
    geometry = Arc(radius=SURVEY_LOOP_RADIUS, start_angle=start_angle, end_angle=end_angle, frame=Frame([0.0, 0.0, 0.0]))

    _assert_arc_follows_its_travel(geometry, clockwise, start_angle, 1.5 * math.pi)


def test_a_clockwise_flag_decides_the_arc_travel_not_its_parametrisation() -> None:
    # point_at runs counter-clockwise a quarter turn; flagged clockwise, the replay
    # sweeps the three quarters the other way round between the same endpoints, and
    # the survey must sample that travel, not the parametrisation.
    geometry = Arc(radius=SURVEY_LOOP_RADIUS, start_angle=0.0, end_angle=0.5 * math.pi, frame=Frame([0.0, 0.0, 0.0]))

    _assert_arc_follows_its_travel(geometry, True, 0.0, 1.5 * math.pi)
