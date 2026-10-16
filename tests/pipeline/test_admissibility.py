"""In-process admissibility: exact cap verdicts on prefixes, exact centre-domain containment, no mutation."""

from __future__ import annotations

import math
from collections.abc import Callable

import pytest
from compas.geometry import Circle
from compas.geometry import Frame
from compas.geometry import Polygon

from benchmarks import survey
from compas_cgal import _stock_2
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.pipeline import admissibility as admissibility_module
from compas_cgal.pipeline.admissibility import LINK_PROBE_SPACING_TOOL_RADII
from compas_cgal.pipeline.admissibility import MIN_LINK_PROBES
from compas_cgal.pipeline.admissibility import MAX_CAP_CHORD_RATIO
from compas_cgal.pipeline.admissibility import Admissibility
from compas_cgal.pipeline.admissibility import CapChordRatio
from compas_cgal.pipeline.admissibility import cap_ratio
from compas_cgal.pipeline.admissibility import centre_domain
from compas_cgal.pipeline.admissibility import straight_probe_positions
from compas_cgal.pipeline.direction import Direction2
from compas_cgal.pipeline.errors import InvalidCapChordRatioError
from compas_cgal.pipeline.errors import InvalidEngagementCapError
from compas_cgal.stock import Stock
from compas_cgal.stock import SweepRemoval
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

TOOL = ToolRadius.build(1.0)
RECT_12X8 = Polygon([[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]])


def _admissibility(cap_deg: float = 120.0) -> Admissibility:
    return Admissibility.build(TOOL, cap_ratio(Radian(math.radians(cap_deg))), centre_domain(RECT_12X8, [], TOOL))


def _loop(cx: float, cy: float, rho: float, phase: float = 0.0) -> IncidentLoop:
    return IncidentLoop.build(Point2[WorldXY].build(cx, cy), GuideRadius.build(rho), Radian(phase))


def _point(x: float, y: float) -> Point2[WorldXY]:
    return Point2[WorldXY].build(x, y)


def _cleared_start() -> Stock:
    # The tool already sits at the start (a plunge or the previous loop cleared it),
    # as every straight the pipeline emits does.
    stock = Stock(RECT_12X8)
    stock.subtract_disk(-4.0, 0.0, TOOL.value)
    return stock


def test_a_straight_slot_into_solid_stock_reads_a_half_turn_and_is_refused_at_120() -> None:
    # In-process a slot reads 180 deg: the 120 cap refuses it, a cap of 180 admits it (0.1 ruling C5).
    assert not _admissibility(120.0).straight_within_cap(_cleared_start(), _point(-4.0, 0.0), _point(4.0, 0.0))
    assert _admissibility(180.0).straight_within_cap(_cleared_start(), _point(-4.0, 0.0), _point(4.0, 0.0))


def test_a_straight_starting_inside_solid_stock_is_full_immersion_at_its_first_station() -> None:
    # Station 0 has no swept prefix behind it: a tool centred in virgin stock reads
    # a full turn, which no cap admits.
    assert not _admissibility(180.0).straight_within_cap(Stock(RECT_12X8), _point(-4.0, 0.0), _point(4.0, 0.0))


def test_a_link_that_leaves_the_centre_domain_is_refused() -> None:
    # The pocket spans x in [-6, 6]; a tool of radius 1 may not be centred past x = 5.
    adm = _admissibility()
    assert adm.straight_in_domain(_point(0.0, 0.0), _point(4.9, 0.0))
    assert not adm.straight_in_domain(_point(0.0, 0.0), _point(5.9, 0.0))


def test_a_loop_in_its_own_cleared_annulus_is_within_any_cap_and_removes_nothing() -> None:
    loop = _loop(0.0, 0.0, 2.0, 0.9)  # non-axis phase on purpose
    stock = Stock(RECT_12X8)
    stock.subtract_circle_sweep_local(loop.centre.x, loop.centre.y, loop.radius.value, TOOL.value)
    adm = _admissibility(1.0)
    assert adm.loop_within_cap(stock, loop, clockwise=True)
    assert adm.removal(stock, _circle_op(loop)) is SweepRemoval.NOTHING


def test_a_loop_into_solid_stock_exceeds_the_cap() -> None:
    assert not _admissibility(120.0).loop_within_cap(Stock(RECT_12X8), _loop(0.0, 0.0, 2.0, 0.9), clockwise=True)


def test_admissibility_never_mutates_the_stock() -> None:
    stock = Stock(RECT_12X8)
    before = stock.clone()
    adm = _admissibility()
    adm.straight_within_cap(stock, _point(-3.0, 0.0), _point(3.0, 1.0))
    adm.loop_within_cap(stock, _loop(0.0, 0.0, 2.0, 0.3), clockwise=True)
    adm.removal(stock, _circle_op(_loop(1.0, 0.0, 1.5, 0.3)))
    assert stock.exactly_equals(before)


def test_straight_probes_overlap_by_a_tool_diameter_and_keep_both_ends() -> None:
    positions = straight_probe_positions(_point(0.0, 0.0), _point(10.0, 0.0), TOOL)
    assert positions[0] == (0.0, 0.0) and positions[-1] == (10.0, 0.0)
    gaps = [later[0] - earlier[0] for earlier, later in zip(positions, positions[1:])]
    assert max(gaps) <= LINK_PROBE_SPACING_TOOL_RADII * TOOL.value
    assert len(straight_probe_positions(_point(0.0, 0.0), _point(0.1, 0.0), TOOL)) == MIN_LINK_PROBES


@pytest.mark.parametrize("cap", [0.0, -0.1, math.pi + 1e-9, math.nan])
def test_a_cap_outside_the_half_turn_is_refused(cap: float) -> None:
    with pytest.raises(InvalidEngagementCapError):
        cap_ratio(Radian(cap))


def _circle_op(loop: IncidentLoop) -> ToolpathOperation:
    circle = loop.circle()
    frame = Frame([loop.centre.x, loop.centre.y, 0.0], circle.frame.xaxis, circle.frame.yaxis)
    return ToolpathOperation(geometry=Circle(loop.radius.value, frame=frame), operation=OperationType.CUT, path_index=0, clockwise=True)


@pytest.mark.parametrize("ratio", [0.0, -1.0, 4.5, math.nan])
def test_a_chord_ratio_no_cap_can_have_is_refused_at_build(ratio: float) -> None:
    """B-minors: `Admissibility.build` validates the surrogate; a bare `CapChordRatio(...)` cannot bypass it."""
    with pytest.raises(InvalidCapChordRatioError):
        Admissibility.build(TOOL, CapChordRatio(ratio), centre_domain(RECT_12X8, [], TOOL))


def test_the_half_turn_cap_has_the_largest_chord_ratio() -> None:
    assert Admissibility.build(TOOL, cap_ratio(Radian(math.pi)), centre_domain(RECT_12X8, [], TOOL)).cap_ratio == MAX_CAP_CHORD_RATIO


def test_the_decisions_are_taken_at_the_gates_own_samples() -> None:
    """I4: loops and straights are decided at exactly the parameters the gate samples (no other station set)."""
    assert admissibility_module.GATE_SAMPLES_PER_MOTION == survey.QUALITY_SAMPLES_PER_MOTION
    assert admissibility_module.loop_parameters() == survey._sample_parameters(Circle(1.0, frame=Frame.worldXY()), survey.QUALITY_SAMPLES_PER_MOTION)
    assert admissibility_module.motion_parameters() == survey._sample_parameters(None, survey.QUALITY_SAMPLES_PER_MOTION)


def test_an_approach_arc_into_virgin_stock_is_refused_and_one_on_a_cut_annulus_is_not() -> None:
    """I2: an arc piece is decided against the cap in process, like every other emitted motion."""
    adm = _admissibility(120.0)
    loop = _loop(0.0, 0.0, 2.0)
    end = _point(-2.0, 0.0)
    virgin = Stock(RECT_12X8)
    assert not adm.arc_within_cap(virgin, loop.centre, loop.entry, end, True)
    virgin.subtract_circle_sweep_local(0.0, 0.0, loop.radius.value, TOOL.value)
    assert adm.arc_within_cap(virgin, loop.centre, loop.entry, end, True)


def _refuse_reporting_calls(monkeypatch: pytest.MonkeyPatch) -> None:
    def reporting_call(*_args: object, **_kwargs: object) -> object:
        raise AssertionError("a decision asked for the reporting prefix samples")

    for name in ("engagement_line_prefixes", "engagement_arc_prefixes", "engagement_circle_prefixes"):
        monkeypatch.setattr(_stock_2, name, reporting_call)


def test_the_cap_decisions_use_the_decision_only_native_queries(monkeypatch: pytest.MonkeyPatch) -> None:
    """Perf round 2: a decision reads only whether some station exceeds, so it never pays for reporting samples."""
    adm = _admissibility(120.0)
    loop = _loop(0.0, 0.0, 2.0, 0.9)
    stock = Stock(RECT_12X8)
    expected = (
        adm.straight_within_cap(_cleared_start(), _point(-4.0, 0.0), _point(4.0, 0.0)),
        adm.loop_within_cap(stock, loop, clockwise=True),
        adm.arc_within_cap(stock, loop.centre, loop.entry, _point(-2.0, 0.0), True),
    )
    _refuse_reporting_calls(monkeypatch)
    assert (
        adm.straight_within_cap(_cleared_start(), _point(-4.0, 0.0), _point(4.0, 0.0)),
        adm.loop_within_cap(stock, loop, clockwise=True),
        adm.arc_within_cap(stock, loop.centre, loop.entry, _point(-2.0, 0.0), True),
    ) == expected


# Generic doubles (no dyadic structure): the witness that the identities below are
# not artefacts of integer fixtures.
GENERIC_LOOPS = [(0.3170000000000002, -0.7310000000000001, 1.9130000000000003, 0.4130000000000001), (-1.1, 0.7, 0.8, 2.9), (0.0, 0.0, 2.0, 0.9)]


def _generic_stock() -> Stock:
    stock = Stock(RECT_12X8)
    stock.subtract_disk(-1.7310000000000001, 0.3170000000000002, 1.0730000000000002)
    stock.subtract_line_sweep_local(-1.7310000000000001, 0.3170000000000002, 0.9130000000000001, -0.2170000000000001, TOOL.value)
    return stock


@pytest.mark.parametrize("clockwise", [True, False])
@pytest.mark.parametrize("loop_args", GENERIC_LOOPS)
@pytest.mark.parametrize("make_stock", [lambda: Stock(RECT_12X8), _generic_stock], ids=["virgin", "generic"])
def test_a_loops_station_zero_is_the_at_rest_entry_probe(make_stock: Callable[[], Stock], loop_args: tuple[float, float, float, float], clockwise: bool) -> None:
    """The exact identity the loop's station-0 reuse rests on: the empty prefix is `engagement_at` at the entry, same stock."""
    stock = make_stock()
    loop = _loop(*loop_args)
    ratio = float(cap_ratio(Radian(math.radians(40.0))))
    zero = _stock_2.engagement_circle_prefixes(stock.raw, loop.centre.x, loop.centre.y, loop.radius.value, loop.entry.x, loop.entry.y, clockwise, [0.0], TOOL.value, ratio, 0.0)[0]
    assert zero == _stock_2.engagement_at(stock.raw, loop.entry.x, loop.entry.y, TOOL.value, ratio, 0.0)


def _partly_cut_annulus(loop: IncidentLoop) -> Stock:
    """The loop's own annulus already cut: the ring reads clear, so the in-process loop is consulted."""
    stock = Stock(RECT_12X8)
    stock.subtract_circle_sweep_local(loop.centre.x, loop.centre.y, loop.radius.value, TOOL.value)
    return stock


CONJUNCTION_STOCKS: dict[str, Callable[[IncidentLoop], Stock]] = {
    "virgin": lambda _loop: Stock(RECT_12X8),
    "generic": lambda _loop: _generic_stock(),
    "own_annulus": _partly_cut_annulus,
}


def _conjunctions(stock_name: str, loop_args: tuple[float, float, float, float], clockwise: bool, cap_deg: float) -> tuple[bool, bool, bool]:
    adm = _admissibility(cap_deg)
    loop = _loop(*loop_args)
    stock = CONJUNCTION_STOCKS[stock_name](loop)
    ring = adm.loop_within_cap_at_rest(stock, loop, Direction2[WorldXY].build(1.0, 0.0))
    return ring, ring and adm.loop_within_cap_past_entry(stock, loop, clockwise), ring and adm.loop_within_cap(stock, loop, clockwise)


@pytest.mark.parametrize("cap_deg", [40.0, 120.0])
@pytest.mark.parametrize("clockwise", [True, False])
@pytest.mark.parametrize("loop_args", GENERIC_LOOPS)
@pytest.mark.parametrize("stock_name", CONJUNCTION_STOCKS)
def test_the_loop_past_its_entry_and_the_at_rest_ring_decide_what_the_full_loop_and_the_ring_decide(
    stock_name: str, loop_args: tuple[float, float, float, float], clockwise: bool, cap_deg: float
) -> None:
    """Conjunction identity: ring AND loop(k >= 1) == ring AND loop(k >= 0), because the ring's entry probe IS station 0."""
    _ring, past_entry, full = _conjunctions(stock_name, loop_args, clockwise, cap_deg)
    assert past_entry == full


def test_the_conjunction_fixtures_are_not_vacuous() -> None:
    """The ring admits on some fixtures, so the past-entry loop is actually consulted there."""
    rings = {_conjunctions("own_annulus", loop_args, True, cap_deg)[0] for loop_args in GENERIC_LOOPS for cap_deg in (40.0, 120.0)}
    assert True in rings


# The gate's L pocket: its one reflex vertex is (6, 6).
L_SHAPE = Polygon([[0.0, 0.0, 0.0], [12.0, 0.0, 0.0], [12.0, 6.0, 0.0], [6.0, 6.0, 0.0], [6.0, 10.0, 0.0], [0.0, 10.0, 0.0]])
# A link chord whose two ends lie 1.0001 r from the reflex vertex, symmetric about
# 225 deg and 0.985 r apart (closer than the old station spacing): its centre path
# enters the vertex's excluded r-disk by about 0.13 r. Generic doubles, not integers.
REFLEX_CHORD_START = (5.036273183746056, 5.732734900084135)
REFLEX_CHORD_END = (5.732734900084135, 5.036273183746056)


def _l_admissibility() -> Admissibility:
    return Admissibility.build(TOOL, cap_ratio(Radian(math.radians(120.0))), centre_domain(L_SHAPE, [], TOOL))


def test_a_link_chord_through_a_reflex_corner_disk_is_refused() -> None:
    """A-I2: both ends are legal centres, yet the chord between them gouges the reflex corner."""
    adm = _l_admissibility()
    assert adm.domain.contains(*REFLEX_CHORD_START) and adm.domain.contains(*REFLEX_CHORD_END)
    assert not adm.straight_in_domain(_point(*REFLEX_CHORD_START), _point(*REFLEX_CHORD_END))


def test_a_straight_at_exactly_tool_radius_from_the_reflex_corner_is_admitted() -> None:
    """The domain is closed: y = 5 passes the reflex vertex (6, 6) and the wall y = 6 at exactly r (integer fixture)."""
    adm = _l_admissibility()
    assert adm.straight_in_domain(_point(4.0, 5.0), _point(8.0, 5.0))
    assert not adm.straight_in_domain(_point(4.0, math.nextafter(5.0, 6.0)), _point(8.0, math.nextafter(5.0, 6.0)))


def test_a_generic_double_straight_clear_of_every_edge_is_admitted() -> None:
    """Generic-double witness: an admitted straight far from the boundary is decided promptly and admitted."""
    adm = _l_admissibility()
    assert adm.straight_in_domain(_point(0.1 + 0.2 + 2.0, math.e), _point(math.pi + 1.0, 3.3))
