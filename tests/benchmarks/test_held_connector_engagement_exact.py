"""Exact sweep types of the connector-engagement certificate: rational / one-root lines, arcs, circle junctions.

Every refusal fixture is checked against exact in-process station verdicts on a
large square minus ONE annulus (``engagement_line_prefixes`` at prefix
parameter 0: the full rim against the stock, no own prefix), so a refusal at
the cap is shown to be a true exceedance of the single-sweep bound, and every
certificate in the property test is checked for exceeded stations.
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
import pytest
from hypothesis import HealthCheck
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.held_connector_engagement import ConnectorEngagementAssessment
from benchmarks.held_connector_engagement import ConnectorTier
from benchmarks.held_connector_engagement import InvalidConnectorEngagementInputError
from benchmarks.held_connector_engagement_arcs import ArcSweep
from benchmarks.held_connector_engagement_exact import ExactConnector
from benchmarks.held_connector_engagement_exact import ExactConnectorAssessment
from benchmarks.held_connector_engagement_exact import ExactConnectorCertificate
from benchmarks.held_connector_engagement_exact import ExactConnectorPath
from benchmarks.held_connector_engagement_exact import RefusalReason
from benchmarks.held_connector_engagement_exact import UncertifiedExactConnectorError
from benchmarks.held_connector_engagement_lines import LineSweep
from benchmarks.held_connector_engagement_lines import RationalLine
from benchmarks.held_connector_engagement_range import DegenerateSweepError
from benchmarks.held_connector_engagement_range import DiscontinuousConnectorError
from benchmarks.held_connector_engagement_range import MixedRootLineSweepError
from benchmarks.held_connector_engagement_range import PointOffSweepSupportError
from benchmarks.held_connector_engagement_range import PriorCircle
from benchmarks.held_connector_engagement_range import RangeEngagementLemma
from benchmarks.held_connector_engagement_range import SweepKind
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_boundary import Sense
from benchmarks.held_omega_connector_engagement import circle_junction_point
from benchmarks.held_one_root import OneRootNumber
from benchmarks.held_one_root import OneRootPoint2
from benchmarks.held_one_root import compare_on_arc
from compas_cgal import _stock_2
from compas_cgal.motion import EngagementCap
from compas_cgal.units import WorldXY

R1 = Fraction(1)
# Half-width of the square stock the single-annulus oracle carves, mm: far
# beyond every fixture so the outer boundary never touches a tool disk.
ORACLE_HALF_WIDTH_MM = 50.0
# Direction length of the dummy segment a parameter-0 station is read on, mm;
# at parameter 0 the own prefix is empty, so the heading does not matter.
STATION_HEADING_MM = 1e-3
# Arc stations per property example (binary64 samples of the exact arc).
ARC_STATIONS = 96


def _p(x: Fraction | int, y: Fraction | int) -> OneRootPoint2:
    return OneRootPoint2.rational(Fraction(x), Fraction(y))


def _e(x: Fraction | int, y: Fraction | int) -> ExactPoint2[WorldXY]:
    point: ExactPoint2[WorldXY] = ExactPoint2[WorldXY].build(Fraction(x), Fraction(y))
    return point


def _on_circle(centre: tuple[Fraction, Fraction], radius: Fraction, x_offset: Fraction, side: int) -> OneRootPoint2:
    """``centre + (x_offset, side·√(radius² − x_offset²))``: a one-root point exactly on the circle."""
    return OneRootPoint2.build(OneRootNumber.of(centre[0] + x_offset), OneRootNumber.build(centre[1], side, radius * radius - x_offset * x_offset))


def _chord_ratio(cos_sq_half_cap: Fraction) -> float:
    ratio = 4 * (1 - cos_sq_half_cap)
    assert Fraction(float(ratio)) == ratio
    return float(ratio)


def _lemma(chord_ratio: float) -> RangeEngagementLemma:
    return RangeEngagementLemma.build(R1, chord_ratio)


def _cap_ratio(degrees: float) -> float:
    return EngagementCap.build(math.radians(degrees)).chord_ratio


def _stations(circle: tuple[float, float, float], points: list[tuple[float, float]], chord_ratio: float) -> list[tuple[float, bool]]:
    """Exact station verdicts (radians, exceeded) at ``points`` on a large square minus ONE annulus."""
    cx, cy, rho = circle
    square = np.array(
        [
            (-ORACLE_HALF_WIDTH_MM, -ORACLE_HALF_WIDTH_MM, 0.0),
            (ORACLE_HALF_WIDTH_MM, -ORACLE_HALF_WIDTH_MM, 0.0),
            (ORACLE_HALF_WIDTH_MM, ORACLE_HALF_WIDTH_MM, 0.0),
            (-ORACLE_HALF_WIDTH_MM, ORACLE_HALF_WIDTH_MM, 0.0),
        ],
        dtype=np.float64,
    )
    stock = _stock_2.Stock2(square, [])
    stock.subtract_circle_sweep(cx, cy, rho, float(R1))
    rows = []
    for x, y in points:
        (sample,) = _stock_2.engagement_line_prefixes(stock, x, y, x + STATION_HEADING_MM, y, [0.0], float(R1), chord_ratio, 0.0)
        rows.append((sample[1], sample[2]))
    return rows


def _arc_samples(arc: ArcSweep, count: int) -> list[tuple[float, float]]:
    cx, cy = float(arc.centre.x), float(arc.centre.y)
    sx, sy = arc.start.as_floats()
    ex, ey = arc.end.as_floats()
    start = math.atan2(sy - cy, sx - cx)
    sweep = math.tau if arc.full else (arc.sense * (math.atan2(ey - cy, ex - cx) - start)) % math.tau
    radius = float(arc.radius)
    return [(cx + radius * math.cos(start + arc.sense * sweep * k / count), cy + radius * math.sin(start + arc.sense * sweep * k / count)) for k in range(count + 1)]


# --------------------------------------------------------------------------
# Lines with exact rational endpoints (not binary64), one-root and cross-root.
# --------------------------------------------------------------------------

# c = (0, 0), rho = 1, r = 1 (R = 2, no hole arc); the far end (7/5, 1/5) has
# D = 2, where the outer test is an equality at C = 1/8, s = 7/2.
UNIT = PriorCircle(_e(0, 0), Fraction(1))
EQUALITY_LINE_END = _p(Fraction(7, 5), Fraction(1, 5))


def test_rational_line_endpoints_are_not_binary64_and_equality_certifies() -> None:
    assert Fraction(float(Fraction(7, 5))) != Fraction(7, 5)
    sweep = LineSweep.build(_p(1, 0), EQUALITY_LINE_END)
    ratio = _chord_ratio(Fraction(1, 8))
    assert _lemma(ratio).certifies(sweep, UNIT)
    assert not _lemma(math.nextafter(ratio, 0.0)).certifies(sweep, UNIT)


# Cross-root line on the rational support y = 1/2 from (√2, 1/2) to (√3, 1/2);
# circle c = (3/2, 2), rho = 1.6: the foot (3/2, 1/2) lies strictly inside and
# decides the hole arc (31.13 deg there; both ends read less).
CROSS_LINE_CIRCLE = PriorCircle(_e(Fraction(3, 2), 2), Fraction(16, 10))
CROSS_LINE_SUPPORT = RationalLine.build(_e(0, Fraction(1, 2)), (Fraction(1), Fraction(0)))
CROSS_LINE = LineSweep.build(
    OneRootPoint2.build(OneRootNumber.build(0, 1, 2), OneRootNumber.of(Fraction(1, 2))),
    OneRootPoint2.build(OneRootNumber.build(0, 1, 3), OneRootNumber.of(Fraction(1, 2))),
    CROSS_LINE_SUPPORT,
)


def test_cross_root_line_is_decided_by_its_interior_foot_and_the_refusal_is_true() -> None:
    assert _lemma(_cap_ratio(31.2)).certifies(CROSS_LINE, CROSS_LINE_CIRCLE)
    assert not _lemma(_cap_ratio(31.0)).certifies(CROSS_LINE, CROSS_LINE_CIRCLE)
    rows = _stations((1.5, 2.0, 1.6), [(1.5, 0.5)], _cap_ratio(31.0))
    assert rows[0][1], "the foot station must exceed the 31.0 deg cap on the single-annulus stock"
    assert not _stations((1.5, 2.0, 1.6), [(1.5, 0.5)], _cap_ratio(31.2))[0][1]


def test_cross_root_line_without_a_rational_support_is_refused_by_name() -> None:
    with pytest.raises(MixedRootLineSweepError):
        LineSweep.build(CROSS_LINE.start, CROSS_LINE.end)


def test_line_endpoint_off_its_support_is_named() -> None:
    with pytest.raises(PointOffSweepSupportError):
        LineSweep.build(_p(0, 0), _p(1, 1), CROSS_LINE_SUPPORT)


def test_zero_length_line_is_named() -> None:
    with pytest.raises(DegenerateSweepError):
        LineSweep.build(_p(1, 0), _p(1, 0))


def test_common_extension_line_splits_at_its_exact_midpoint() -> None:
    start = OneRootPoint2.build(OneRootNumber.build(0, 1, 2), OneRootNumber.of(0))
    end = OneRootPoint2.build(OneRootNumber.build(1, 1, 2), OneRootNumber.build(1, -1, 2))
    halves = LineSweep.build(start, end).split()
    assert halves is not None
    first, second = halves
    assert first.start == start and second.end == end and first.end == second.start
    assert first.end == OneRootPoint2.build(OneRootNumber.build(Fraction(1, 2), 1, 2), OneRootNumber.build(Fraction(1, 2), Fraction(-1, 2), 2))


# --------------------------------------------------------------------------
# Arcs: the tool centre on a circle of ∂Ω*; extremes at the ends or at the
# interior near / far points of the prior circle's centre.
# --------------------------------------------------------------------------

# Arc about q = (3/10, 4/10), radius 1, CCW from (13/10, 4/10) to (3/10, 14/10).
# Against c = (0, 0), rho = 1 (R = 2) the far point (9/10, 6/5) is interior,
# D = 9/4, and the outer test is an equality at C = 1/16 (s = 15/4).
ARC_CENTRE = (Fraction(3, 10), Fraction(4, 10))
RATIONAL_ARC = ArcSweep.build(_e(*ARC_CENTRE), Fraction(1), 1, _p(Fraction(13, 10), Fraction(4, 10)), _p(Fraction(3, 10), Fraction(14, 10)))
# The same circle up to (11/10, 1): its sweep stops before the far direction.
SHORT_ARC = ArcSweep.build(_e(*ARC_CENTRE), Fraction(1), 1, _p(Fraction(13, 10), Fraction(4, 10)), _p(Fraction(11, 10), 1))


def test_arc_interior_far_point_equality_certifies_and_next_binary64_refuses() -> None:
    ratio = _chord_ratio(Fraction(1, 16))
    assert _lemma(ratio).certifies(RATIONAL_ARC, UNIT)
    below = math.nextafter(ratio, 0.0)
    assert not _lemma(below).certifies(RATIONAL_ARC, UNIT)
    assert _lemma(below).certifies(SHORT_ARC, UNIT), "the ends alone pass: only the interior far point refuses"


def test_arc_refused_at_the_cap_is_a_true_exceedance() -> None:
    ratio = _chord_ratio(Fraction(17, 256))
    assert not _lemma(ratio).certifies(RATIONAL_ARC, UNIT)
    assert _lemma(ratio).certifies(SHORT_ARC, UNIT)
    ((_, exceeded),) = _stations((0.0, 0.0, 1.0), [(0.9, 1.2)], ratio)
    assert exceeded


# Cross-root arc about q = (1/2, 1/2) (|q - c|² = 1/2 is not a square), radius
# 1, CCW from q + (9/10, √19/10) to q + (1/2, √3/2): ends in Q(√19) and Q(√3),
# the far point (L + 1)² in Q(√2). C = 1/1024: the ends pass, the far point fails.
CROSS_ARC = ArcSweep.build(
    _e(Fraction(1, 2), Fraction(1, 2)),
    Fraction(1),
    1,
    _on_circle((Fraction(1, 2), Fraction(1, 2)), Fraction(1), Fraction(9, 10), 1),
    _on_circle((Fraction(1, 2), Fraction(1, 2)), Fraction(1), Fraction(1, 2), 1),
)


def test_cross_root_arc_refuses_at_its_far_point_and_the_station_exceeds() -> None:
    assert CROSS_ARC.start.radicand != CROSS_ARC.end.radicand
    ratio = _chord_ratio(Fraction(1, 1024))
    assert not _lemma(ratio).certifies(CROSS_ARC, UNIT)
    far = 0.5 + math.sqrt(0.5)
    ((_, exceeded),) = _stations((0.0, 0.0, 1.0), [(far, far)], ratio)
    assert exceeded
    assert _lemma(_chord_ratio(Fraction(1, 2048))).certifies(CROSS_ARC, UNIT)


# Hole branch: c = (0, 0), rho = 3/2 (h = 1/2, d*² = 3/4); arc about (0, 1),
# radius 1/2, CCW from (-1/2, 1) through the near point (0, 1/2) (D = 1/4) to
# (1/2, 1) (D = 5/4 at both ends). Only the near point puts d* inside the range;
# the d* test is an equality at C = 3/4 (s = 1, a 60 deg cap).
HOLE_CIRCLE = PriorCircle(_e(0, 0), Fraction(3, 2))
HOLE_ARC = ArcSweep.build(_e(0, 1), Fraction(1, 2), 1, _p(Fraction(-1, 2), 1), _p(Fraction(1, 2), 1))


def test_arc_near_point_brings_the_interior_hole_minimum_into_range() -> None:
    assert _lemma(1.0).certifies(HOLE_ARC, HOLE_CIRCLE)
    assert not _lemma(math.nextafter(1.0, 0.0)).certifies(HOLE_ARC, HOLE_CIRCLE)
    ratio = _cap_ratio(56.0)
    assert not _lemma(ratio).certifies(HOLE_ARC, HOLE_CIRCLE)
    ((_, exceeded),) = _stations((0.0, 0.0, 1.5), [(math.sqrt(3.0) / 4, 0.75)], ratio)
    assert exceeded, "the station at d = d* on the arc reads the 60 deg hole arc"


def test_arc_endpoint_off_the_circle_is_named() -> None:
    with pytest.raises(PointOffSweepSupportError):
        ArcSweep.build(_e(0, 0), Fraction(1), 1, _p(1, 0), _p(1, 1))


def test_arc_split_point_is_rational_on_the_circle_and_strictly_inside() -> None:
    halves = CROSS_ARC.split()
    assert halves is not None
    first, second = halves
    middle = first.end
    assert middle == second.start and middle.radicand == 0
    dx, dy = middle.minus(CROSS_ARC.centre)
    assert dx * dx + dy * dy == 1
    assert first.start == CROSS_ARC.start and second.end == CROSS_ARC.end
    assert not first.full and not second.full


def test_full_circle_arc_covers_both_extremes() -> None:
    start = _p(Fraction(13, 10), Fraction(4, 10))
    full = ArcSweep.build(_e(*ARC_CENTRE), Fraction(1), -1, start, start)
    assert full.full
    # the far point (9/10, 6/5) is on the full turn, so the equality cap still certifies and below refuses
    ratio = _chord_ratio(Fraction(1, 16))
    assert _lemma(ratio).certifies(full, UNIT)
    assert not _lemma(math.nextafter(ratio, 0.0)).certifies(full, UNIT)


@st.composite
def one_root_arcs(draw: st.DrawFn) -> tuple[ArcSweep, PriorCircle, float]:
    """Arcs with ends in two different extensions near a prior circle, and a cap."""
    eighth = st.integers(min_value=-16, max_value=16).map(lambda k: Fraction(k, 8))
    centre = (draw(eighth), draw(eighth))
    radius = Fraction(draw(st.integers(min_value=2, max_value=24)), 8)
    offsets = st.integers(min_value=-99, max_value=99).map(lambda k: radius * Fraction(k, 100))
    start = _on_circle(centre, radius, draw(offsets), draw(st.sampled_from((1, -1))))
    end = _on_circle(centre, radius, draw(offsets), draw(st.sampled_from((1, -1))))
    sense: Sense = draw(st.sampled_from((1, -1)))
    arc = ArcSweep.build(_e(*centre), radius, sense, start, end)
    angle = draw(st.floats(min_value=0.0, max_value=2 * math.pi))
    rho = draw(st.floats(min_value=0.1, max_value=4.0))
    sx, sy = start.as_floats()
    reach = rho + draw(st.floats(min_value=-0.5, max_value=0.5))
    cx, cy = sx - reach * math.cos(angle), sy - reach * math.sin(angle)
    circle = PriorCircle(ExactPoint2[WorldXY].from_binary64(cx, cy), Fraction(rho))
    return arc, circle, _cap_ratio(draw(st.floats(min_value=20.0, max_value=180.0)))


@settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(one_root_arcs())
def test_certified_one_root_arc_has_no_exceeded_station(case: tuple[ArcSweep, PriorCircle, float]) -> None:
    arc, circle, ratio = case
    if not _lemma(ratio).certifies(arc, circle):
        return
    rows = _stations((float(circle.centre.x), float(circle.centre.y), float(circle.guide_radius)), _arc_samples(arc, ARC_STATIONS), ratio)
    assert not any(exceeded for _, exceeded in rows)


def test_the_property_draws_reach_certificates() -> None:
    """The property above is not vacuous: a short arc next to its own circle is certified."""
    centre = (Fraction(0), Fraction(0))
    arc = ArcSweep.build(_e(*centre), Fraction(2), 1, _on_circle(centre, Fraction(2), Fraction(199, 100), 1), _on_circle(centre, Fraction(2), Fraction(197, 100), 1))
    sx, sy = arc.start.as_floats()
    circle = PriorCircle(ExactPoint2[WorldXY].from_binary64(sx - 1.0, sy), Fraction(1))
    assert _lemma(_cap_ratio(80.0)).certifies(arc, circle)


# --------------------------------------------------------------------------
# Circle junctions: binary64 circle -> exact site as a short exact line.
# --------------------------------------------------------------------------


def test_circle_junction_point_is_rational_and_exactly_on_the_binary64_circle() -> None:
    site = _on_circle((Fraction(0), Fraction(0)), Fraction(3), Fraction(1), 1)  # (1, √8)
    centre, radius = (0.25, 0.5), 2.0
    junction = circle_junction_point(centre, radius, site)
    assert junction.radicand == 0
    dx, dy = junction.minus(ExactPoint2[WorldXY].from_binary64(*centre))
    assert dx * dx + dy * dy == Fraction(radius) ** 2
    jx, jy = junction.as_floats()
    sx, sy = site.as_floats()
    expected = math.hypot(sx - centre[0], sy - centre[1]) - radius  # the gap, radial
    assert math.isclose(math.hypot(jx - sx, jy - sy), abs(expected), rel_tol=1e-9)


def _junction_path(gap: float, *, cross_root: bool) -> ExactConnectorPath:
    """Circle 0 tangent from inside to the line y = 0 near the site, then circle 1 further along.

    The site of circle 0 is (0, 0) (or one-root on a unit arc when
    ``cross_root``); circle 0 has centre (0, rho + gap), so its junction line
    has length ``gap``.
    """
    rho = 1.0
    if cross_root:
        site0 = _on_circle((Fraction(0), Fraction(1)), Fraction(1), Fraction(1, 3), -1)
        site1 = _on_circle((Fraction(0), Fraction(1)), Fraction(1), Fraction(1, 2), -1)
        motion: LineSweep | ArcSweep = ArcSweep.build(_e(0, 1), Fraction(1), 1, site0, site1)
    else:
        site0, site1 = _p(0, 0), _p(Fraction(1, 5), 0)
        motion = LineSweep.build(site0, site1)
    x0, y0 = site0.as_floats()
    x1, y1 = site1.as_floats()
    circles = (
        PriorCircle(ExactPoint2[WorldXY].from_binary64(x0, y0 + rho + gap), Fraction(rho)),
        PriorCircle(ExactPoint2[WorldXY].from_binary64(x1, y1 + rho), Fraction(rho)),
    )
    exit_point = circle_junction_point((x0, y0 + rho + gap), rho, site0)
    entry_point = circle_junction_point((x1, y1 + rho), rho, site1)
    sweeps = (
        LineSweep.build(exit_point, site0, kind=SweepKind.CIRCLE_JUNCTION),
        motion,
        LineSweep.build(site1, entry_point, kind=SweepKind.CIRCLE_JUNCTION),
    )
    return ExactConnectorPath.build(circles, (ExactConnector.build(sweeps),), (exit_point, entry_point))


def test_short_junction_connector_is_certified() -> None:
    assessment = ExactConnectorAssessment.build(_junction_path(2.0**-40, cross_root=False), R1, _cap_ratio(80.0))
    assert assessment.refused == ()
    assert assessment.verdicts[0].tier is ConnectorTier.TIER_1


def test_cross_root_junction_connector_is_certified() -> None:
    path = _junction_path(2.0**-40, cross_root=True)
    assert path.connectors[0].sweeps[0].end.radicand != path.connectors[0].sweeps[2].start.radicand
    assessment = ExactConnectorAssessment.build(path, R1, _cap_ratio(80.0))
    assert assessment.refused == ()


def test_long_junction_is_refused_by_name_and_is_a_true_exceedance() -> None:
    path = _junction_path(0.5, cross_root=False)
    assessment = ExactConnectorAssessment.build(path, R1, _cap_ratio(80.0))
    refusal = assessment.refusals[0]
    assert (refusal.connector, refusal.sweep, refusal.kind, refusal.reason) == (0, 0, SweepKind.CIRCLE_JUNCTION, RefusalReason.NOT_PROVED_AT_DEPTH)
    assert [r.sweep for r in assessment.refusals] == [0, 1, 2], "every refused sweep is named, not only the first"
    ((_, exceeded),) = _stations((0.0, 1.5, 1.0), [(0.0, 0.0)], _cap_ratio(80.0))
    assert exceeded
    with pytest.raises(UncertifiedExactConnectorError) as raised:
        ExactConnectorCertificate.build(path, R1, _cap_ratio(80.0))
    assert raised.value.connectors == (0,)
    assert raised.value.refusals == assessment.refusals


# --------------------------------------------------------------------------
# Path contract.
# --------------------------------------------------------------------------


def test_discontinuous_connector_is_named() -> None:
    with pytest.raises(DiscontinuousConnectorError):
        ExactConnector.build((LineSweep.build(_p(0, 0), _p(1, 0)), LineSweep.build(_p(1, 1), _p(2, 1))))


def test_connector_count_must_match_the_seam() -> None:
    circles = (UNIT, UNIT)
    connector = ExactConnector.build((LineSweep.build(_p(1, 0), _p(2, 0)),))
    with pytest.raises(InvalidConnectorEngagementInputError):
        ExactConnectorPath.build(circles, (connector, connector, connector), (_p(1, 0), _p(2, 0)))


def test_connector_must_start_and_end_at_its_circles_junctions() -> None:
    connector = ExactConnector.build((LineSweep.build(_p(1, 0), _p(2, 0)),))
    with pytest.raises(InvalidConnectorEngagementInputError):
        ExactConnectorPath.build((UNIT, UNIT), (connector,), (_p(1, 0), _p(3, 0)))


def test_a_future_circle_never_certifies() -> None:
    # connector 0 runs far from circle 0 and right next to circle 1: only circle 1 (cut later) would certify it
    far = PriorCircle(_e(-10, 0), Fraction(1))
    near = PriorCircle(_e(5, 1), Fraction(1))
    sweeps = (LineSweep.build(_p(-10, -1), _p(5, 0)),)
    path = ExactConnectorPath.build((far, near), (ExactConnector.build(sweeps),), (_p(-10, -1), _p(5, 0)))
    assessment = ExactConnectorAssessment.build(path, R1, _cap_ratio(80.0))
    assert assessment.refused == (0,)


def test_binary64_injection_agrees_with_the_integer_grid_on_figure5() -> None:
    from benchmarks.held_toolpath_benchmark import generate

    path = generate("figure5_refined")
    cap = EngagementCap.build(math.radians(float(path.case.tea_cap)))
    grid = ConnectorEngagementAssessment.build(path.circles, path.connectors, path.case.tool_radius, cap)
    exact = ExactConnectorAssessment.build(ExactConnectorPath.from_figure5(path.circles, path.connectors), Fraction(float(path.case.tool_radius.value)), cap.chord_ratio)
    assert [verdict.tier for verdict in exact.verdicts] == [verdict.tier for verdict in grid.verdicts]


@settings(max_examples=200, deadline=None)
@given(one_root_arcs(), st.integers(min_value=-20, max_value=20), st.integers(min_value=-20, max_value=20))
def test_filtered_arc_order_equals_compare_on_arc(case: tuple[ArcSweep, PriorCircle, float], dx: int, dy: int) -> None:
    """The sign-filtered order rule decides exactly as the repo's exact ``compare_on_arc``."""
    arc = case[0]
    if dx == 0 and dy == 0:
        return
    for point in (_p(arc.centre.x + Fraction(dx, 7), arc.centre.y + Fraction(dy, 7)), arc.start, arc.end):
        assert arc.order(point) == compare_on_arc(arc.centre, arc.sense, arc.start, point, arc.end)
