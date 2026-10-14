"""The survey reads the replay's exact sweeps: authored loop radii, three-valued removal."""

import math
from dataclasses import replace

import numpy as np
import pytest
from compas.geometry import Circle
from compas.geometry import Frame
from compas.geometry import Line
from compas.geometry import Polygon

from benchmarks.depletion import CutMotion
from benchmarks.errors import InconsistentRemovalVerdictError
from benchmarks.spec import PocketSpec
from benchmarks.survey import _motion_prefix_samples
from benchmarks.survey import survey_path
from compas_cgal import _stock_2
from compas_cgal.engagement import _cap_chord_ratio
from compas_cgal.engagement import _subtract_operation
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.stock import Stock
from compas_cgal.stock import SweepRemoval
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import WorldXY

TOOL_RADIUS = 1.0
GUIDE_RADIUS = 2.5
TRIALS = 12
PREFIX_SAMPLES = 9
_POCKET = Polygon([[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]])
SPEC = PocketSpec.build(name="synthetic_12x8", family="analytic", polygon=_POCKET, tool_diameter=2.0, tea_cap_deg=120.0)


def _loop(phase: float) -> ToolpathOperation:
    ux, uy = math.cos(phase), math.sin(phase)
    frame = Frame([0.5, 0.25, 0.0], [ux, uy, 0.0], [-uy, ux, 0.0])
    return ToolpathOperation(geometry=Circle(GUIDE_RADIUS, frame=frame), operation=OperationType.CUT, path_index=0, clockwise=True)


def _op(geometry: object, operation: OperationType) -> ToolpathOperation:
    return ToolpathOperation(geometry=geometry, operation=operation, path_index=0)


@pytest.mark.parametrize("k", range(TRIALS))
def test_a_loop_through_its_own_cut_annulus_engages_nothing_past_its_entry(k: int) -> None:
    """Stations p > 0 lie on the AUTHORED radius: every cutter disk is in the annulus.

    On the circle through the double ``point_at(0)`` -- ulps off the authored
    radius for a non-axis phase -- the rims poke out of the annulus and read a
    sliver of engagement. Station 0 IS that double entry (where the tool is),
    so it may read the sliver; an exactly incident entry reads nothing (below).
    """
    stock = Stock(_POCKET)
    _subtract_operation(stock, _loop(2 * math.pi * k / TRIALS), TOOL_RADIUS)
    motion = CutMotion(index=0, operation=_loop(2 * math.pi * k / TRIALS + 1.0), stock=stock)
    samples = _motion_prefix_samples(motion, PREFIX_SAMPLES, TOOL_RADIUS, _cap_chord_ratio(SPEC.tea_cap_rad))
    assert [total for total, _run, _over in samples[1:]] == [0.0] * (PREFIX_SAMPLES - 1)


@pytest.mark.parametrize("k", range(TRIALS))
def test_an_incident_loop_through_its_own_cut_annulus_engages_nothing_at_any_station(k: int) -> None:
    """With the entry exactly on the authored circle, station 0 reads nothing too."""
    authored = IncidentLoop.build(Point2[WorldXY].build(0.5, 0.25), GuideRadius.build(GUIDE_RADIUS), Radian(2 * math.pi * k / TRIALS))
    stock = Stock(_POCKET)
    stock.subtract_circle_sweep(authored.centre.x, authored.centre.y, authored.radius.value, TOOL_RADIUS)
    parameters = [index / (PREFIX_SAMPLES - 1) for index in range(PREFIX_SAMPLES)]
    samples = _stock_2.engagement_circle_prefixes(
        stock.raw,
        authored.centre.x,
        authored.centre.y,
        authored.radius.value,
        authored.entry.x,
        authored.entry.y,
        True,
        parameters,
        TOOL_RADIUS,
        _cap_chord_ratio(SPEC.tea_cap_rad),
    )
    assert [total for total, _run, _over in samples] == [0.0] * PREFIX_SAMPLES


def test_the_survey_carries_the_verdict_and_counts_within_band_motions() -> None:
    oblique = _op(Line([-4.0, -3.0, 0.0], [3.0, 2.0, 0.0]), OperationType.CUT)
    axis = _op(Line([-4.0, 3.0, 0.0], [4.0, 3.0, 0.0]), OperationType.CUT)
    operations = [
        _op(Line([-4.0, -3.0, 4.0], [-4.0, -3.0, 0.0]), OperationType.PLUNGE),
        oblique,
        oblique,
        _op(Line([3.0, 2.0, 0.0], [-4.0, 3.0, 0.0]), OperationType.LINK),
        axis,
        axis,
    ]
    survey = survey_path(SPEC, ToolpathResult(operations=operations, polyline=np.zeros((0, 3))), samples_per_motion=4)
    verdicts = [motion.removal for motion in survey.motions]
    assert verdicts == [SweepRemoval.REMOVES, SweepRemoval.WITHIN_BAND, SweepRemoval.REMOVES, SweepRemoval.REMOVES, SweepRemoval.NOTHING]
    assert [motion.removes_material for motion in survey.motions] == [True, True, True, True, False]
    assert survey.within_band_motions == 1


def test_a_record_whose_flag_contradicts_its_verdict_is_refused() -> None:
    oblique = _op(Line([-4.0, -3.0, 0.0], [3.0, 2.0, 0.0]), OperationType.CUT)
    operations = [_op(Line([-4.0, -3.0, 4.0], [-4.0, -3.0, 0.0]), OperationType.PLUNGE), oblique]
    motion = survey_path(SPEC, ToolpathResult(operations=operations, polyline=np.zeros((0, 3))), samples_per_motion=4).motions[0]
    with pytest.raises(InconsistentRemovalVerdictError):
        replace(motion, removes_material=False)


def _bridge_ending_at(ex: float, ey: float) -> ToolpathOperation:
    """A non-Pythagorean bridge with (ex, ey) as an exact endpoint.

    compas recomputes `Line.end` from start + direction * length, so only the
    START is stored bit for bit; the swept region of a segment does not depend
    on its travel direction, so the bridge is authored from the entry outwards.
    """
    return _op(Line([ex, ey, 0.0], [-4.3, -3.1, 0.0]), OperationType.CUT)


@pytest.mark.parametrize("phase", [0.3, 1.1, 2.0, 4.4])
def test_a_loop_entered_from_a_bridge_reads_no_engagement_at_its_entry(phase: float) -> None:
    """Station 0 is the authored entry, which is where the bridge ended.

    The bridge is non-Pythagorean, so its sides are a certified band; its end
    disks are exact, so a cutter at the bridge's end meets nothing.
    """
    loop = _loop(phase)
    assert isinstance(loop.geometry, Circle)
    entry = loop.geometry.point_at(0.0)
    ex, ey = float(entry[0]), float(entry[1])
    bridge = _bridge_ending_at(ex, ey)
    stock = Stock(_POCKET)
    _subtract_operation(stock, bridge, TOOL_RADIUS)
    motion = CutMotion(index=1, operation=loop, stock=stock)
    total, _run, _over = _motion_prefix_samples(motion, PREFIX_SAMPLES, TOOL_RADIUS, _cap_chord_ratio(SPEC.tea_cap_rad))[0]
    assert total == 0.0


def _zero_length(x: float, y: float, operation: OperationType) -> ToolpathOperation:
    return _op(Line([x, y, 0.0], [x, y, 0.0]), operation)


def test_a_zero_length_link_between_two_loops_is_one_station_and_reads_nothing() -> None:
    """Ruling 0.2-zerolen: a degenerate motion is measured, not refused."""
    loop = _loop(0.7)
    assert isinstance(loop.geometry, Circle)
    entry = loop.geometry.point_at(0.0)
    ex, ey = float(entry[0]), float(entry[1])
    operations = [
        _op(Line([ex, ey, 4.0], [ex, ey, 0.0]), OperationType.PLUNGE),
        loop,
        _zero_length(ex, ey, OperationType.LINK),
        loop,
    ]
    survey = survey_path(SPEC, ToolpathResult(operations=operations, polyline=np.zeros((0, 3))), samples_per_motion=4)
    link = survey.motions[1]
    assert link.kind.value == "line"
    assert len(link.samples) == 1
    assert link.samples[0].engagement_deg == 0.0
    assert (link.samples[0].position.x, link.samples[0].position.y) == (ex, ey)


def test_a_zero_length_cut_inside_material_reads_the_disks_engagement() -> None:
    plunge = (-4.0, -3.0)
    point = (3.0, 2.0)
    operations = [
        _op(Line([plunge[0], plunge[1], 4.0], [plunge[0], plunge[1], 0.0]), OperationType.PLUNGE),
        _op(Line([plunge[0], plunge[1], 0.0], [-3.0, -3.0, 0.0]), OperationType.CUT),
        _zero_length(point[0], point[1], OperationType.CUT),
    ]
    survey = survey_path(SPEC, ToolpathResult(operations=operations, polyline=np.zeros((0, 3))), samples_per_motion=4)
    stationary = survey.motions[1]
    assert len(stationary.samples) == 1
    reference = Stock(_POCKET)
    reference.subtract_disk(plunge[0], plunge[1], TOOL_RADIUS)
    _subtract_operation(reference, operations[1], TOOL_RADIUS)
    _total, run, exceeded = _stock_2.engagement_at(reference.raw, point[0], point[1], TOOL_RADIUS, _cap_chord_ratio(SPEC.tea_cap_rad))
    assert stationary.samples[0].engagement_deg == math.degrees(run)
    assert stationary.samples[0].cap_exceeded is exceeded
    assert stationary.cap_exceeded is exceeded
    assert stationary.peak_engagement_deg == 360.0
