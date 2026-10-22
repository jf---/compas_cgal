"""The native declared-boundary clearance oracle against the Python oracle: identical verdicts and reported primitives.

Fixtures are exactly representable (the exact-centre-domain pockets: Pythagorean
directions, dyadic or small rational data), so exact touches are reachable by
construction: ∂Ω* pieces touch the boundary at exactly r, circles at (10, 15)
of radius 4 touch two primitives at once. Every comparison is exact equality of
the verdict and of the reported primitive; reported depths are binary64 views
and are not compared here (the differential tool reports them).
"""

from __future__ import annotations

from fractions import Fraction

import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.held_centre_domain_supports import OffsetLineSupport
from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_clearance import ArcElement
from benchmarks.held_declared_clearance import CircleElement
from benchmarks.held_declared_clearance import ClearanceVerdict
from benchmarks.held_declared_clearance import DeclaredClearance
from benchmarks.held_declared_clearance import Element
from benchmarks.held_declared_clearance import SegmentElement
from benchmarks.held_declared_clearance_native import NativeDeclaredClearance
from benchmarks.held_declared_clearance_native import last_counters
from benchmarks.held_declared_clearance_wire import offset_owners
from benchmarks.held_exact_centre_domain import CentreDomain
from benchmarks.held_exact_centre_domain import CentreDomainPiece
from benchmarks.held_gouge_audit import Clearance
from benchmarks.held_gouge_audit import HeldGougeAuditReport
from benchmarks.held_omega_connectors import OmegaConnectors
from benchmarks.held_clearance_backend import ClearanceBackend
from benchmarks.held_omega_gouge_audit import audit_omega_path
from benchmarks.held_omega_gouge_audit import audit_polygon_path
from benchmarks.held_omega_gouge_audit_native import audit_omega_path_native
from benchmarks.held_omega_gouge_audit_native import audit_polygon_path_native
from benchmarks.held_omega_sites import single_loop_pieces
from benchmarks.held_one_root import OneRootNumber
from benchmarks.held_one_root import OneRootPoint2
from compas_cgal import _held_declared_clearance_2
from compas_cgal.units import WorldXY
from tests.benchmarks.test_held_exact_centre_domain import concave_arc_pocket
from tests.benchmarks.test_held_exact_centre_domain import convex_arc_pocket
from tests.benchmarks.test_held_exact_centre_domain import l_shape_pocket
from tests.benchmarks.test_held_omega_path import _arc_contacts
from tests.benchmarks.test_held_omega_path import _omega_circles
from tests.benchmarks.test_held_omega_path import _polygon_path

F = Fraction
R = F(1)
_FIXTURES = {"concave": concave_arc_pocket, "convex": convex_arc_pocket, "l_shape": l_shape_pocket}
_ORACLES: dict[str, tuple[DeclaredBoundary, DeclaredClearance, NativeDeclaredClearance]] = {}


def _oracles(name: str) -> tuple[DeclaredBoundary, DeclaredClearance, NativeDeclaredClearance]:
    if name not in _ORACLES:
        boundary = _FIXTURES[name]()
        _ORACLES[name] = (boundary, DeclaredClearance.build(boundary, R), NativeDeclaredClearance.build(boundary, R))
    return _ORACLES[name]


def _decision(verdict: ClearanceVerdict) -> tuple[Clearance, int | None]:
    return verdict.clearance, verdict.primitive


def _point(x: Fraction, y: Fraction) -> ExactPoint2[WorldXY]:
    return ExactPoint2[WorldXY].build(x, y)


def _piece_element(piece: CentreDomainPiece) -> Element:
    support = piece.support
    if isinstance(support, OffsetLineSupport):
        return SegmentElement(piece.source, piece.target, support.direction, support.normal, support.level)
    return ArcElement(support.centre, support.radius, support.sense, piece.source, piece.target)


def _decisions(report: HeldGougeAuditReport) -> list[tuple[object, ...]]:
    return [(item.element, item.index, item.clearance, item.primitive, item.primitive_kind, item.segment) for item in report.elements]


# ---------------------------------------------------------------------------
# Exact touches and gouges
# ---------------------------------------------------------------------------


def test_native_decides_touch_and_gouge_exactly() -> None:
    """Convex arc about (10, 7.5), ρ = 12.5, top y = 20: the circle about (10, 15) of radius 4 touches both; 2^-40 more gouges."""
    _, _, native = _oracles("convex")
    centre = _point(F(10), F(15))
    touch, gouge = native.circles([CircleElement(centre, F(4)), CircleElement(centre, F(4) + F(1, 2**40))])
    assert touch.clearance is Clearance.TOUCH
    assert gouge.clearance is Clearance.GOUGE
    left, right, lower = _point(F(5), F(19)), _point(F(15), F(19)), _point(F(15), F(19) + F(1, 2**40))
    flush, dipping = native.segments_binary64([(left, right), (left, lower)])
    assert flush.clearance is Clearance.TOUCH
    assert dipping.clearance is Clearance.GOUGE


@pytest.mark.parametrize("fixture", sorted(_FIXTURES))
def test_every_omega_piece_touches_through_its_owner_lemma_and_reports_like_python(fixture: str) -> None:
    boundary, python, native = _oracles(fixture)
    pieces = single_loop_pieces(CentreDomain.build(boundary, R))
    elements = [_piece_element(piece) for piece in pieces]
    verdicts = native.elements(elements, [offset_owners(piece.owners) for piece in pieces])
    counters = last_counters()
    assert [_decision(verdict) for verdict in verdicts] == [_decision(python.element(element)) for element in elements]
    assert all(verdict.clearance is Clearance.TOUCH for verdict in verdicts)
    assert counters.owner_touches > 0


@pytest.mark.parametrize("fixture", sorted(_FIXTURES))
def test_a_false_owner_claim_never_changes_a_verdict(fixture: str) -> None:
    """Every primitive named as an owner of every piece: the lemma's hypotheses reject the false ones."""
    boundary, python, native = _oracles(fixture)
    pieces = single_loop_pieces(CentreDomain.build(boundary, R))
    elements = [_piece_element(piece) for piece in pieces]
    everyone = tuple(range(len(boundary.primitives)))
    verdicts = native.elements(elements, [everyone] * len(elements))
    assert [_decision(verdict) for verdict in verdicts] == [_decision(python.element(element)) for element in elements]


def test_an_offset_segment_beyond_its_owner_is_decided_in_full() -> None:
    """The L-shape's bottom line y = 0, x in [0, 20]: its offset y = 1 extended to x = 25 leaves the owner's support region."""
    boundary, python, native = _oracles("l_shape")
    start, end = OneRootPoint2.rational(F(2), F(1)), OneRootPoint2.rational(F(25), F(1))
    element = SegmentElement(start, end, (F(1), F(0)), (F(0), F(1)), F(1))
    (verdict,) = native.elements([element], [(0,)])
    assert last_counters().owner_touches == 0
    assert _decision(verdict) == _decision(python.element(element))
    assert verdict.clearance is Clearance.GOUGE


# ---------------------------------------------------------------------------
# Differential: random elements on every fixture
# ---------------------------------------------------------------------------

_COORDINATE = st.floats(min_value=-2, max_value=22, allow_nan=False)


@settings(max_examples=300, deadline=None)
@given(st.sampled_from(tuple(_FIXTURES)), _COORDINATE, _COORDINATE, st.floats(min_value=2**-10, max_value=12), _COORDINATE, _COORDINATE)
def test_binary64_circles_and_segments_decide_like_python(fixture: str, x: float, y: float, radius: float, u: float, v: float) -> None:
    _, python, native = _oracles(fixture)
    centre, other = ExactPoint2[WorldXY].from_binary64(x, y), ExactPoint2[WorldXY].from_binary64(u, v)
    (circle,) = native.circles([CircleElement(centre, F(radius))])
    assert _decision(circle) == _decision(python.circle(CircleElement(centre, F(radius))))
    if centre != other:
        (segment,) = native.segments_binary64([(centre, other)])
        assert _decision(segment) == _decision(python.element(SegmentElement.between(centre, other)))


def _on_circle(centre: ExactPoint2[WorldXY], radius: Fraction, t: Fraction) -> OneRootPoint2:
    """The rational point of the circle at half-angle tangent t."""
    scale = 1 + t * t
    return OneRootPoint2.rational(centre.x + radius * (1 - t * t) / scale, centre.y + radius * 2 * t / scale)


def _on_circle_at_x(centre: ExactPoint2[WorldXY], radius: Fraction, x: Fraction, upper: bool) -> OneRootPoint2:
    """The one-root point of the circle above or below x (y = cy ± √(ρ² - (x - cx)²))."""
    return OneRootPoint2.build(OneRootNumber.of(x), OneRootNumber.build(centre.y, 1 if upper else -1, radius * radius - (x - centre.x) ** 2))


_SMALL = st.fractions(min_value=-2, max_value=22, max_denominator=64)
_TANGENT = st.fractions(min_value=-3, max_value=3, max_denominator=16)


@settings(max_examples=300, deadline=None)
@given(st.sampled_from(tuple(_FIXTURES)), _SMALL, _SMALL, st.fractions(min_value=F(1, 8), max_value=14, max_denominator=32), _TANGENT, _TANGENT, st.sampled_from((1, -1)))
def test_rational_arcs_decide_like_python(fixture: str, cx: Fraction, cy: Fraction, radius: Fraction, t0: Fraction, t1: Fraction, sense: int) -> None:
    _, python, native = _oracles(fixture)
    centre = _point(cx, cy)
    element = ArcElement(centre, radius, 1 if sense == 1 else -1, _on_circle(centre, radius, t0), _on_circle(centre, radius, t1))
    (verdict,) = native.elements([element])
    assert _decision(verdict) == _decision(python.element(element))


@settings(max_examples=200, deadline=None)
@given(
    st.sampled_from(tuple(_FIXTURES)),
    _SMALL,
    _SMALL,
    st.fractions(min_value=F(1, 2), max_value=14, max_denominator=32),
    st.fractions(min_value=F(-15, 16), max_value=F(15, 16), max_denominator=64),
    st.booleans(),
    _TANGENT,
    st.sampled_from((1, -1)),
)
def test_one_root_arcs_and_segments_decide_like_python(fixture: str, cx: Fraction, cy: Fraction, radius: Fraction, offset: Fraction, upper: bool, t: Fraction, sense: int) -> None:
    """Endpoints in a quadratic extension: an arc from a one-root point to a rational one, and a vertical segment from it."""
    _, python, native = _oracles(fixture)
    centre = _point(cx, cy)
    x = cx + offset * radius
    irrational = _on_circle_at_x(centre, radius, x, upper)
    arc = ArcElement(centre, radius, 1 if sense == 1 else -1, irrational, _on_circle(centre, radius, t))
    above = OneRootPoint2.build(irrational.x, irrational.y + 1)
    segment = SegmentElement(irrational, above, (F(0), F(1)), (F(-1), F(0)), -x)
    elements: list[Element] = [arc, segment]
    verdicts = native.elements(elements)
    assert [_decision(verdict) for verdict in verdicts] == [_decision(python.element(element)) for element in elements]


def test_one_batch_equals_singletons() -> None:
    _, _, native = _oracles("concave")
    elements = [CircleElement(_point(F(x), F(y)), F(r)) for x, y, r in ((10, 10, 3), (2, 2, 1), (10, 6, F(1, 2)), (18, 18, 1))]
    assert native.circles(elements) == [native.circles([element])[0] for element in elements]


# ---------------------------------------------------------------------------
# Whole-path audits on the concave-arc fixture
# ---------------------------------------------------------------------------


def test_native_audits_report_like_the_python_audits_on_the_concave_fixture() -> None:
    boundary, python, native = _oracles("concave")
    domain = CentreDomain.build(boundary, R)
    circles = _omega_circles(domain, python, _arc_contacts())
    connectors = OmegaConnectors.build(single_loop_pieces(domain), circles).connectors
    assert _decisions(audit_omega_path_native(native, circles, connectors)) == _decisions(audit_omega_path(python, circles, connectors, ClearanceBackend.PYTHON))
    polygon_circles, polylines = _polygon_path(boundary, _arc_contacts())
    expected = audit_polygon_path(python, polygon_circles, polylines, ClearanceBackend.PYTHON)
    found = audit_polygon_path_native(native, polygon_circles, polylines)
    assert _decisions(found) == _decisions(expected)
    assert found.gouging_circles


def test_the_audit_switch_runs_the_native_oracle_with_identical_decisions() -> None:
    boundary, python, _ = _oracles("concave")
    domain = CentreDomain.build(boundary, R)
    circles = _omega_circles(domain, python, _arc_contacts())
    connectors = OmegaConnectors.build(single_loop_pieces(domain), circles).connectors
    native = audit_omega_path(python, circles, connectors, ClearanceBackend.NATIVE)
    assert _decisions(native) == _decisions(audit_omega_path(python, circles, connectors, ClearanceBackend.PYTHON))
    polygon_circles, polylines = _polygon_path(boundary, _arc_contacts())
    assert _decisions(audit_polygon_path(python, polygon_circles, polylines, ClearanceBackend.NATIVE)) == _decisions(
        audit_polygon_path(python, polygon_circles, polylines, ClearanceBackend.PYTHON)
    )


# ---------------------------------------------------------------------------
# Named failure modes
# ---------------------------------------------------------------------------


def test_named_errors() -> None:
    boundary, _, native = _oracles("l_shape")
    with pytest.raises(_held_declared_clearance_2.NonPositiveClearanceRadiusError):
        NativeDeclaredClearance.build(boundary, F(0))
    with pytest.raises(_held_declared_clearance_2.NonPositiveClearanceRadiusError):
        NativeDeclaredClearance.build(boundary, F(1, 3))
    with pytest.raises(_held_declared_clearance_2.NonPositiveClearanceRadiusError):
        native.circles([CircleElement(_point(F(5), F(5)), F(0))])
    with pytest.raises(_held_declared_clearance_2.DegenerateClearanceElementError):
        native.segments_binary64([(_point(F(5), F(5)), _point(F(5), F(5)))])
    mixed = OneRootPoint2(OneRootNumber.build(F(1), F(1), F(2)), OneRootNumber.build(F(1), F(1), F(3)))
    with pytest.raises(_held_declared_clearance_2.MixedRootPointError):
        native.elements([SegmentElement(mixed, OneRootPoint2.rational(F(5), F(5)), (F(1), F(0)), (F(0), F(1)), F(5))])
    with pytest.raises(_held_declared_clearance_2.InvalidClearanceWireError):
        native.handle.elements(
            [(0, 0, [99], ((("1", "1"), ("0", "1"), ("0", "1")),) * 2, ((("2", "1"), ("0", "1"), ("0", "1")),) * 2, [("1", "1"), ("0", "1"), ("0", "1"), ("1", "1"), ("0", "1")])]
        )
    with pytest.raises(_held_declared_clearance_2.InvalidDeclaredPrimitiveError):
        _held_declared_clearance_2.DeclaredClearance2.build([(0, 0, [("0", "1"), ("0", "1"), ("1", "1"), ("0", "1"), ("1", "1"), ("1", "1"), ("1", "1")])], ("1", "1"))
