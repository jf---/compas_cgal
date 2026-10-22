"""Exact continuous connector-engagement certificate (tiers 1 and 1b) against exact station oracles.

The numbered sections follow the failing-test list in
``docs/held_connector_engagement.md``. Every "no false certificate" claim is
checked against exact in-process station verdicts (``engagement_line_prefixes``)
on a replayed stock: a single annulus for the lemma, the tier-2 replay
(``Stock_OVER``, so an exceeded station is a true exceedance) for whole paths.
"""

from __future__ import annotations

import math
import time
from fractions import Fraction

import numpy as np
import pytest
from hypothesis import HealthCheck
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.held_connector_engagement import MAX_BISECTION_DEPTH
from benchmarks.held_connector_engagement import TIERS_BUDGET_S
from benchmarks.held_connector_engagement import ConnectorEngagementAssessment
from benchmarks.held_connector_engagement import ConnectorEngagementCertificate
from benchmarks.held_connector_engagement import ConnectorEngagementLemma
from benchmarks.held_connector_engagement import ConnectorTier
from benchmarks.held_connector_engagement import InvalidConnectorEngagementInputError
from benchmarks.held_connector_engagement import LemmaCircle
from benchmarks.held_connector_engagement import UncertifiedConnectorEngagementError
from benchmarks.held_connector_engagement import moving_pieces
from benchmarks.held_figure5_boundary_path import Figure5BoundaryConnectors
from benchmarks.held_figure5_boundary_path import Figure5BoundaryTransition
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_boundary_path import Figure5SeamClosure
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from compas_cgal import _stock_2
from compas_cgal.motion import EngagementCap
from compas_cgal.pipeline.admissibility import motion_parameters
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY
from tools.held_connector_engagement_tier2 import StationVerdict
from tools.held_connector_engagement_tier2 import replay_full
from tools.held_connector_engagement_tier2 import replay_windowed

R1 = ToolRadius.build(1.0)
CAP_80 = EngagementCap.build(math.radians(80.0))
# The design's dense station count for the single-annulus property: 1,000
# stations (1,001 prefix parameters), twenty times the gate's 46.
DENSE_STATIONS = 1000
# Reporting slack between two reported angles (a station's native arc angle
# and the lemma's acos), degrees; never a decision (the exact check is the
# exceeded flag). Near an internal tangency of the rim, a reported run angle
# carries sqrt-of-ulp error: 2 sqrt(2 delta / r) with a distance error delta
# of 1e-13 mm (100x the coordinate ulp at 10 mm), r = 1, is 9e-7 rad = 5e-5 deg.
REPORTED_ANGLE_SLACK_DEG = 1e-4
# Half-width of the square stock the single-annulus oracle carves, mm: far
# beyond every fixture so the outer boundary never touches a tool disk.
ORACLE_HALF_WIDTH_MM = 50.0


def _p(x: float, y: float) -> Point2[WorldXY]:
    point: Point2[WorldXY] = Point2[WorldXY].build(x, y)
    return point


def _circle(index: int, centre: tuple[float, float], contact: tuple[float, float], rho: float) -> Figure5CounterclockwiseCircle:
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), Fraction(index, 10_000), 4)
    return Figure5CounterclockwiseCircle(GuideRunId(0), GuideRunStationOrdinal(index), site, site, _p(*centre), _p(*contact), GuideRadius.build(rho))


def _transition(points: tuple[tuple[float, float], ...], *, within_fan: bool = False) -> Figure5BoundaryTransition:
    length = sum(math.dist(a, b) for a, b in zip(points, points[1:]))
    return Figure5BoundaryTransition(tuple(_p(*point) for point in points), Millimetre(length), within_fan)


def _connectors(transitions: tuple[Figure5BoundaryTransition, ...], seam: Figure5SeamClosure) -> Figure5BoundaryConnectors:
    return Figure5BoundaryConnectors(transitions, seam, Millimetre(0.0))


def _square(half: float) -> np.ndarray:
    return np.array([(-half, -half, 0.0), (half, -half, 0.0), (half, half, 0.0), (-half, half, 0.0)], dtype=np.float64)


def _single_annulus_stations(
    circle: tuple[tuple[float, float], float], a: tuple[float, float], b: tuple[float, float], cap: EngagementCap, stations: int
) -> list[tuple[float, float, bool]]:
    """Exact station verdicts of piece ``a -> b`` on a large square minus ONE annulus."""
    (cx, cy), rho = circle
    stock = _stock_2.Stock2(_square(ORACLE_HALF_WIDTH_MM), [])
    stock.subtract_circle_sweep(cx, cy, rho, float(R1.value))
    parameters = [k / stations for k in range(stations + 1)]
    samples: list[tuple[float, float, bool]] = _stock_2.engagement_line_prefixes(stock, a[0], a[1], b[0], b[1], parameters, float(R1.value), cap.chord_ratio, 0.0)
    return samples


def _lemma(cap: EngagementCap = CAP_80) -> ConnectorEngagementLemma:
    return ConnectorEngagementLemma.build(R1, cap)


def _rational(point: tuple[float, float]) -> tuple[Fraction, Fraction]:
    return Fraction(point[0]), Fraction(point[1])


def _lemma_circle(centre: tuple[float, float], rho: float) -> LemmaCircle:
    return LemmaCircle(_rational(centre), Fraction(rho))


# --------------------------------------------------------------------------
# 1. Lemma against the oracle (Hypothesis): no false certificate.
# --------------------------------------------------------------------------

caps_deg = st.floats(min_value=20.0, max_value=180.0)


@st.composite
def wall_connectors(draw: st.DrawFn) -> tuple[tuple[tuple[float, float], float], tuple[float, float], tuple[float, float], EngagementCap]:
    """A circle tangent to the wall y = 0 (contact on the offset y = 1) and a connector along that offset."""
    rho = draw(st.floats(min_value=0.05, max_value=6.0))
    x0 = draw(st.floats(min_value=-3.0, max_value=3.0))
    length = draw(st.floats(min_value=0.001, max_value=3.0))
    cap = EngagementCap.build(math.radians(draw(caps_deg)))
    return ((x0, 1.0 + rho), rho), (x0, 1.0), (x0 + length, 1.0), cap


@st.composite
def near_circle_segments(draw: st.DrawFn) -> tuple[tuple[tuple[float, float], float], tuple[float, float], tuple[float, float], EngagementCap]:
    """A random straight piece starting on or near a circle (as a Held connector does), any heading."""
    rho = draw(st.one_of(st.floats(min_value=0.05, max_value=0.9), st.floats(min_value=1.0, max_value=6.0)))
    cx, cy = draw(st.floats(min_value=-3.0, max_value=3.0)), draw(st.floats(min_value=-3.0, max_value=3.0))
    phi = draw(st.floats(min_value=0.0, max_value=2 * math.pi))
    offset = draw(st.floats(min_value=-0.3, max_value=0.3))
    a = (cx + (rho + offset) * math.cos(phi), cy + (rho + offset) * math.sin(phi))
    length, heading = draw(st.floats(min_value=0.001, max_value=1.5)), draw(st.floats(min_value=0.0, max_value=2 * math.pi))
    b = (a[0] + length * math.cos(heading), a[1] + length * math.sin(heading))
    cap = EngagementCap.build(math.radians(draw(caps_deg)))
    return ((cx, cy), rho), a, b, cap


def _check_single_annulus(circle: tuple[tuple[float, float], float], a: tuple[float, float], b: tuple[float, float], cap: EngagementCap, stations: int) -> None:
    if a == b:
        return
    lemma = _lemma(cap)
    exact_circle = _lemma_circle(*circle)
    samples = _single_annulus_stations(circle, a, b, cap, stations)
    measured = max(math.degrees(sample[1]) for sample in samples)
    bound = lemma.bound_deg(_rational(a), _rational(b), exact_circle)
    assert measured <= bound + REPORTED_ANGLE_SLACK_DEG
    if lemma.piece_certified_by(_rational(a), _rational(b), exact_circle):
        assert not any(sample[2] for sample in samples)


@settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(wall_connectors())
def test_1_wall_connector_certificate_has_no_exceeded_station(case: tuple[tuple[tuple[float, float], float], tuple[float, float], tuple[float, float], EngagementCap]) -> None:
    _check_single_annulus(*case, stations=len(motion_parameters()) - 1)
    _check_single_annulus(*case, stations=DENSE_STATIONS)


@settings(max_examples=150, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(near_circle_segments())
def test_1_random_piece_certificate_has_no_exceeded_station(case: tuple[tuple[tuple[float, float], float], tuple[float, float], tuple[float, float], EngagementCap]) -> None:
    _check_single_annulus(*case, stations=len(motion_parameters()) - 1)


def _station_exceeded(circle: tuple[tuple[float, float], float], a: tuple[float, float], b: tuple[float, float], cap: EngagementCap) -> bool:
    return any(sample[2] for sample in _single_annulus_stations(circle, a, b, cap, DENSE_STATIONS))


def test_1_hole_arc_is_checked_and_its_removal_is_caught(monkeypatch: pytest.MonkeyPatch) -> None:
    """The piece dips into a big circle's hole: the outer arc passes, the hole arc does not."""
    circle, a, b = ((0.0, 0.0), 3.0), (3.0, 0.0), (0.0, 1.5)
    assert _station_exceeded(circle, a, b, CAP_80)
    assert not _lemma().piece_certified_by(_rational(a), _rational(b), _lemma_circle(*circle))
    monkeypatch.setattr(ConnectorEngagementLemma, "hole_interval_within", lambda *_args: True)
    assert _lemma().piece_certified_by(_rational(a), _rational(b), _lemma_circle(*circle)), "defect injection must create the false certificate"


def test_1_loose_outer_test_is_caught(monkeypatch: pytest.MonkeyPatch) -> None:
    """N_O <= 0 alone is not enough: the outer arc at d ~ 1.58 is about 162 degrees."""
    circle, a, b = ((0.0, 0.0), 1.0), (1.0, 0.0), (1.5, 0.5)
    assert _station_exceeded(circle, a, b, CAP_80)
    assert not _lemma().piece_certified_by(_rational(a), _rational(b), _lemma_circle(*circle))
    monkeypatch.setattr(ConnectorEngagementLemma, "outer_within", lambda self, squared_distance, outer_radius: squared_distance + self.tool_radius**2 - outer_radius**2 <= 0)
    assert _lemma().piece_certified_by(_rational(a), _rational(b), _lemma_circle(*circle)), "defect injection must create the false certificate"


@st.composite
def wall_paths(draw: st.DrawFn) -> tuple[tuple[Figure5CounterclockwiseCircle, ...], tuple[Figure5BoundaryTransition, ...], EngagementCap]:
    """An open Held-like path along the wall y = 0: circles tangent to it, connectors along the offset y = 1.

    A drawn kink lifts a connector off the offset so non-Pythagorean pieces
    (UNDER/OVER line sweeps) are exercised too.
    """
    count = draw(st.integers(min_value=2, max_value=5))
    xs = [0.0]
    for _ in range(count - 1):
        xs.append(xs[-1] + draw(st.floats(min_value=0.05, max_value=2.5)))
    rhos = [draw(st.floats(min_value=0.3, max_value=3.0)) for _ in range(count)]
    circles = tuple(_circle(k, (x, 1.0 + rho), (x, 1.0), rho) for k, (x, rho) in enumerate(zip(xs, rhos)))
    transitions = []
    for x0, x1 in zip(xs, xs[1:]):
        lift = draw(st.sampled_from((0.0, 0.0, 0.25, 0.7)))
        middle = ((x0 + x1) / 2, 1.0 + lift)
        transitions.append(_transition(((x0, 1.0), middle, (x1, 1.0)) if lift else ((x0, 1.0), (x1, 1.0))))
    cap = EngagementCap.build(math.radians(draw(st.floats(min_value=30.0, max_value=150.0))))
    return circles, tuple(transitions), cap


# Rectangle pocket for the path fixtures: the wall y = 0, wide and tall
# enough that no tool disk reaches its other sides, mm.
POCKET = tuple(_p(x, y) for x, y in ((-10.0, 0.0), (25.0, 0.0), (25.0, 20.0), (-10.0, 20.0)))


@settings(max_examples=25, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(wall_paths())
def test_1_no_false_certificate_on_the_replayed_path(case: tuple[tuple[Figure5CounterclockwiseCircle, ...], tuple[Figure5BoundaryTransition, ...], EngagementCap]) -> None:
    circles, transitions, cap = case
    assessment = ConnectorEngagementAssessment.build(circles, _connectors(transitions, Figure5SeamClosure.COINCIDENT), R1, cap)
    stations = {row.connector: row for row in replay_full(POCKET, circles, transitions, R1, cap, motion_parameters())}
    for verdict in assessment.verdicts:
        if verdict.certified:
            assert stations[verdict.connector].verdict is not StationVerdict.VIOLATED
            assert stations[verdict.connector].exceeded_over == 0


@settings(max_examples=40, deadline=None)
@given(wall_paths())
def test_1_integer_grid_decides_as_the_rationals_do(case: tuple[tuple[Figure5CounterclockwiseCircle, ...], tuple[Figure5BoundaryTransition, ...], EngagementCap]) -> None:
    """The assessment runs the lemma on one integer dyadic grid; its tier-1 verdicts equal the Fraction lemma's."""
    circles, transitions, cap = case
    assessment = ConnectorEngagementAssessment.build(circles, _connectors(transitions, Figure5SeamClosure.COINCIDENT), R1, cap)
    lemma = _lemma(cap)
    for verdict, circle, transition in zip(assessment.verdicts, circles[: len(transitions)], transitions, strict=True):
        own = LemmaCircle.from_circle(circle)
        assert (verdict.tier is ConnectorTier.TIER_1) == all(lemma.piece_certified_by(a, b, own) for a, b in moving_pieces(transition))


# --------------------------------------------------------------------------
# 2. Far endpoint decides on a straight tangent wall.
# --------------------------------------------------------------------------


@settings(max_examples=200, deadline=None)
@given(wall_connectors())
def test_2_far_endpoint_decides_on_a_tangent_wall(case: tuple[tuple[tuple[float, float], float], tuple[float, float], tuple[float, float], EngagementCap]) -> None:
    (centre, rho), a, b, cap = case
    lemma = _lemma(cap)
    circle = _lemma_circle(centre, rho)
    far = _rational(b)
    squared = (far[0] - circle.centre[0]) ** 2 + (far[1] - circle.centre[1]) ** 2
    expected = lemma.outer_within(squared, circle.guide_radius + lemma.tool_radius)
    assert lemma.piece_certified_by(_rational(a), far, circle) is expected


# --------------------------------------------------------------------------
# 3. The interior hole minimum d* is checked (constructed fixture).
# --------------------------------------------------------------------------

DSTAR_CAP = EngagementCap.build(math.radians(22.0))
DSTAR_CIRCLE = ((0.0, 0.0), 1.2)
DSTAR_PIECE = ((1.1, 0.0), (0.9, 0.0))


def test_3_interior_hole_minimum_refuses(monkeypatch: pytest.MonkeyPatch) -> None:
    lemma = _lemma(DSTAR_CAP)
    circle = _lemma_circle(*DSTAR_CIRCLE)
    a, b = (_rational(point) for point in DSTAR_PIECE)
    h = circle.guide_radius - lemma.tool_radius
    assert lemma.hole_within(Fraction(0.9) ** 2, h) and lemma.hole_within(Fraction(1.1) ** 2, h), "both endpoint hole arcs pass"
    assert lemma.outer_within(Fraction(1.1) ** 2, circle.guide_radius + lemma.tool_radius)
    assert not lemma.piece_certified_by(a, b, circle)
    assert _station_exceeded(DSTAR_CIRCLE, *DSTAR_PIECE, DSTAR_CAP), "the oracle confirms the d* exceedance"

    def endpoints_only(self: ConnectorEngagementLemma, d_min_sq: Fraction, d_max_sq: Fraction, inner: Fraction) -> bool:
        return self.hole_within(d_min_sq, inner) and self.hole_within(d_max_sq, inner)

    monkeypatch.setattr(ConnectorEngagementLemma, "hole_interval_within", endpoints_only)
    assert _lemma(DSTAR_CAP).piece_certified_by(a, b, circle), "deleting the d* branch must create the false certificate"


# --------------------------------------------------------------------------
# 4. Equality certifies; the next binary64 above refuses. No epsilon.
# --------------------------------------------------------------------------


def test_4_equality_certifies_in_the_squared_half() -> None:
    """s = 2 (C = 1/2), r = 1, rho = 4 (R = 5), point (3, 3): N = -6, N^2 = 36 = 4 C D r^2."""
    lemma = ConnectorEngagementLemma.from_chord_ratio(R1, 2.0)
    circle = _lemma_circle((0.0, 0.0), 4.0)
    assert lemma.cos_sq_half_cap == Fraction(1, 2)
    assert lemma.outer_within(Fraction(18), Fraction(5))
    assert lemma.piece_certified_by(_rational((4.0, 0.0)), _rational((3.0, 3.0)), circle)
    above = math.nextafter(3.0, math.inf)
    assert not lemma.piece_certified_by(_rational((4.0, 0.0)), _rational((above, 3.0)), circle)


def test_4_equality_certifies_in_the_sign_half() -> None:
    """s = 4 (C = 0), r = 3, rho = 2 (R = 5), point (4, 0): N = 16 + 9 - 25 = 0."""
    lemma = ConnectorEngagementLemma.from_chord_ratio(ToolRadius.build(3.0), 4.0)
    circle = _lemma_circle((0.0, 0.0), 2.0)
    assert lemma.cos_sq_half_cap == 0
    assert lemma.piece_certified_by(_rational((0.0, 2.0)), _rational((4.0, 0.0)), circle)
    assert not lemma.piece_certified_by(_rational((0.0, 2.0)), _rational((math.nextafter(4.0, math.inf), 0.0)), circle)


@pytest.mark.parametrize("ratio", (0.0, -1.0, 4.000000000000001, math.inf, math.nan))
def test_4_chord_surrogate_outside_its_domain_is_named(ratio: float) -> None:
    with pytest.raises(InvalidConnectorEngagementInputError):
        ConnectorEngagementLemma.from_chord_ratio(R1, ratio)


# --------------------------------------------------------------------------
# 5. Future circles never certify.
# --------------------------------------------------------------------------


def _future_only_path() -> tuple[tuple[Figure5CounterclockwiseCircle, ...], tuple[Figure5BoundaryTransition, ...]]:
    """Circle 0 is small and left behind; circle 1 is huge and tangent to the connector's far end.

    Only circle 1's annulus bounds the connector, and circle 1 is cut after it.
    """
    circles = (_circle(0, (0.0, 0.0), (1.0, 0.0), 1.0), _circle(1, (3.0, 100.0), (3.0, 0.0), 100.0))
    return circles, (_transition(((1.0, 0.0), (3.0, 0.0))),)


def test_5_future_circle_never_certifies(monkeypatch: pytest.MonkeyPatch) -> None:
    circles, transitions = _future_only_path()
    connectors = _connectors(transitions, Figure5SeamClosure.COINCIDENT)
    later = _lemma().piece_certified_by(_rational((1.0, 0.0)), _rational((3.0, 0.0)), LemmaCircle.from_circle(circles[1]))
    assert later, "the fixture is meaningful only if the future circle would certify"
    verdict = ConnectorEngagementAssessment.build(circles, connectors, R1, CAP_80).verdicts[0]
    assert verdict.tier is ConnectorTier.REFUSED
    with pytest.raises(UncertifiedConnectorEngagementError) as caught:
        ConnectorEngagementCertificate.build(circles, connectors, R1, CAP_80)
    assert caught.value.connectors == (0,)

    import benchmarks.held_connector_engagement as module

    monkeypatch.setattr(module, "prior_circle_indices", lambda connector, circle_count: range(circle_count - 1, -1, -1))
    injected = ConnectorEngagementAssessment.build(circles, connectors, R1, CAP_80).verdicts[0]
    assert injected.tier is ConnectorTier.TIER_1B, "allowing k > i must create the false certificate"


# --------------------------------------------------------------------------
# 6. Seam, COINCIDENT seam, zero-motion fan connectors.
# --------------------------------------------------------------------------


def _closed_seam_path() -> tuple[tuple[Figure5CounterclockwiseCircle, ...], tuple[Figure5BoundaryTransition, ...]]:
    """Three circles along y = 1; the seam returns from the far circle to circle 0 along the offset.

    The seam leaves circle 2 (small) far behind, so tier 1 refuses it; circle 0
    is huge and tangent at the seam's end, so only an earlier circle certifies it.
    """
    circles = (
        _circle(0, (0.0, 101.0), (0.0, 1.0), 100.0),
        _circle(1, (1.0, 1.5), (1.0, 1.0), 0.5),
        _circle(2, (2.0, 1.5), (2.0, 1.0), 0.5),
    )
    transitions = (
        _transition(((0.0, 1.0), (1.0, 1.0))),
        _transition(((1.0, 1.0), (2.0, 1.0))),
        _transition(((2.0, 1.0), (0.0, 1.0))),
    )
    return circles, transitions


def test_6_closed_path_gets_n_verdicts_and_the_seam_uses_every_circle() -> None:
    circles, transitions = _closed_seam_path()
    certificate = ConnectorEngagementCertificate.build(circles, _connectors(transitions, Figure5SeamClosure.CONNECTOR), R1, CAP_80)
    assert len(certificate.verdicts) == len(circles)
    seam = certificate.verdicts[-1]
    assert seam.seam and seam.connector == len(circles) - 1
    assert seam.tier is ConnectorTier.TIER_1B
    assert 0 in {sub.circle for sub in seam.sub_pieces}, "the seam is certified by circle 0, cut before it"
    assert all(sub.circle <= verdict.connector for verdict in certificate.verdicts[:-1] for sub in verdict.sub_pieces)


def test_6_coincident_seam_gets_n_minus_one_verdicts() -> None:
    circles, transitions = _closed_seam_path()
    certificate = ConnectorEngagementCertificate.build(circles, _connectors(transitions[:-1], Figure5SeamClosure.COINCIDENT), R1, CAP_80)
    assert len(certificate.verdicts) == len(circles) - 1
    assert not any(verdict.seam for verdict in certificate.verdicts)


def test_6_seam_kind_must_match_the_connector_count() -> None:
    circles, transitions = _closed_seam_path()
    with pytest.raises(InvalidConnectorEngagementInputError):
        ConnectorEngagementAssessment.build(circles, _connectors(transitions[:-1], Figure5SeamClosure.CONNECTOR), R1, CAP_80)


def test_6_zero_motion_fan_connector_is_certified_without_a_query(monkeypatch: pytest.MonkeyPatch) -> None:
    circles = (_circle(0, (0.0, 2.0), (0.0, 1.0), 1.0), _circle(1, (0.0, 3.0), (0.0, 1.0), 2.0))
    transitions = (_transition(((0.0, 1.0), (0.0, 1.0)), within_fan=True),)

    def no_query(*_args: object) -> bool:
        raise AssertionError("a zero-motion connector must not ask the lemma")

    monkeypatch.setattr(ConnectorEngagementLemma, "piece_certified_by", no_query)
    certificate = ConnectorEngagementCertificate.build(circles, _connectors(transitions, Figure5SeamClosure.COINCIDENT), R1, CAP_80)
    assert certificate.verdicts[0].tier is ConnectorTier.STATIONARY


# --------------------------------------------------------------------------
# 7. Repair (certificate side): a midpoint circle splits a refused connector into certified halves.
# 8. Exhaustion is named: a refusal after MAX_BISECTION_DEPTH names the connector.
# --------------------------------------------------------------------------


def _wall_run(xs: tuple[float, ...], rho: float) -> tuple[tuple[Figure5CounterclockwiseCircle, ...], Figure5BoundaryConnectors]:
    circles = tuple(_circle(k, (x, 1.0 + rho), (x, 1.0), rho) for k, x in enumerate(xs))
    transitions = tuple(_transition(((x0, 1.0), (x1, 1.0))) for x0, x1 in zip(xs, xs[1:]))
    return circles, _connectors(transitions, Figure5SeamClosure.COINCIDENT)


def test_7_midpoint_circle_insertion_repairs_a_refused_connector() -> None:
    """Refinement's repair, on the certificate side: halving the connector drives the outer arc to 0 (q_i lies on circle i)."""
    xs: tuple[float, ...] = (0.0, 1.0)
    circles, connectors = _wall_run(xs, 0.5)
    with pytest.raises(UncertifiedConnectorEngagementError):
        ConnectorEngagementCertificate.build(circles, connectors, R1, CAP_80)
    for _round in range(4):
        xs = tuple(sorted(set(xs) | {(a + b) / 2 for a, b in zip(xs, xs[1:])}))
        circles, connectors = _wall_run(xs, 0.5)
        try:
            certificate = ConnectorEngagementCertificate.build(circles, connectors, R1, CAP_80)
        except UncertifiedConnectorEngagementError:
            continue
        assert all(verdict.certified for verdict in certificate.verdicts)
        return
    pytest.fail("midpoint insertion did not converge in 4 rounds")


def test_8_exhausted_bisection_names_the_connector() -> None:
    circles, transitions = _future_only_path()
    with pytest.raises(UncertifiedConnectorEngagementError, match=r"\b0\b") as caught:
        ConnectorEngagementCertificate.build(circles, _connectors(transitions, Figure5SeamClosure.COINCIDENT), R1, CAP_80)
    assert caught.value.connectors == (0,)
    assert MAX_BISECTION_DEPTH == 12


# --------------------------------------------------------------------------
# 9. Windowing contract (tier 2 oracle): windowed replay == full replay.
# --------------------------------------------------------------------------


def _slanted_path() -> tuple[tuple[Figure5CounterclockwiseCircle, ...], tuple[Figure5BoundaryTransition, ...]]:
    xs = (0.0, 1.3, 2.1, 3.7, 4.4, 6.0)
    rhos = (1.0, 0.6, 1.8, 0.4, 1.1, 0.9)
    circles = tuple(_circle(k, (x, 1.0 + rho), (x, 1.0), rho) for k, (x, rho) in enumerate(zip(xs, rhos)))
    transitions = tuple(_transition(((x0, 1.0), ((x0 + x1) / 2, 1.4), (x1, 1.0))) for x0, x1 in zip(xs, xs[1:]))
    return circles, transitions


def test_9_windowed_replay_equals_the_full_replay() -> None:
    circles, transitions = _slanted_path()
    parameters = motion_parameters()
    full = replay_full(POCKET, circles, transitions, R1, CAP_80, parameters)
    windowed = replay_windowed(POCKET, circles, transitions, R1, CAP_80, parameters, range(len(transitions)))
    assert [row.connector for row in windowed] == [row.connector for row in full]
    for full_row, window_row in zip(full, windowed):
        assert window_row.stations == full_row.stations
        assert window_row.max_run_deg_under == full_row.max_run_deg_under
        assert window_row.exceeded_under == full_row.exceeded_under
        assert window_row.max_run_deg_over == full_row.max_run_deg_over
        assert window_row.exceeded_over == full_row.exceeded_over


# --------------------------------------------------------------------------
# 10. Bracket: an UNDER-band sliver is not a violation.
# --------------------------------------------------------------------------

# A cap below the rim's crossing of the <= 1e-4 r UNDER band of a
# non-Pythagorean line sweep (about sqrt(2e-4) rad, ~0.8 degrees per side).
SLIVER_CAP = EngagementCap.build(math.radians(0.2))


def test_10_under_band_sliver_is_clear_on_stock_over() -> None:
    """The second connector re-traverses the first's non-Pythagorean line: only the UNDER sliver is left."""
    circles = (_circle(0, (0.0, 0.2), (0.0, 0.0), 0.2), _circle(1, (3.0, 2.2), (3.0, 2.0), 0.2), _circle(2, (0.0, 0.2), (0.0, 0.0), 0.2))
    transitions = (_transition(((0.0, 0.0), (3.0, 2.0))), _transition(((3.0, 2.0), (0.0, 0.0))))
    pocket = tuple(_p(x, y) for x, y in ((-5.0, -5.0), (10.0, -5.0), (10.0, 10.0), (-5.0, 10.0)))
    assert _stock_2.line_sweep_model(0.0, 0.0, 3.0, 2.0, float(R1.value)) is not _stock_2.LineSweepModel.EXACT_CAPSULE
    rows = replay_full(pocket, circles, transitions, R1, SLIVER_CAP, motion_parameters())
    second = rows[1]
    assert second.exceeded_under > 0, "the fixture must leave an UNDER sliver the stations engage"
    assert second.exceeded_over == 0
    assert second.verdict is StationVerdict.UNRESOLVED, "inside the bracket: not a violation"


# --------------------------------------------------------------------------
# 11. Workloads: the post-gate on the generated paths, tiers 1 + 1b within budget.
# --------------------------------------------------------------------------


def test_11_figure5_certificate_builds_within_budget() -> None:
    from benchmarks.held_toolpath_benchmark import generate

    path = generate("figure5_refined")
    cap = EngagementCap.build(math.radians(float(path.case.tea_cap)))
    started = time.perf_counter()
    certificate = ConnectorEngagementCertificate.build(path.circles, path.connectors, path.case.tool_radius, cap)
    seconds = time.perf_counter() - started
    assert len(certificate.verdicts) == len(path.transitions)
    assert certificate.verdicts[-1].seam and certificate.verdicts[-1].tier is ConnectorTier.TIER_1B
    assert seconds <= TIERS_BUDGET_S


def test_11_crossed_skis_connectors_are_all_certified_within_budget() -> None:
    """Flipped by the refinement wiring: its 2 connectors without a certificate are repaired by midpoint circles."""
    from benchmarks.held_toolpath_benchmark import generate

    path = generate("figure8_crossed_skis_boundary")
    cap = EngagementCap.build(math.radians(float(path.case.tea_cap)))
    assessment = ConnectorEngagementAssessment.build(path.circles, path.connectors, path.case.tool_radius, cap)
    assert assessment.seconds <= TIERS_BUDGET_S
    assert not assessment.refused
    certificate = ConnectorEngagementCertificate.build(path.circles, path.connectors, path.case.tool_radius, cap)
    assert len(certificate.verdicts) == len(path.transitions)
