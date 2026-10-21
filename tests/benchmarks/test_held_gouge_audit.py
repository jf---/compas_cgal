"""Exact gouge audit on tiny analytic pockets with closed-form expected depths."""

from __future__ import annotations

import math
from fractions import Fraction

import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.held_gouge_audit import AuditCircle
from benchmarks.held_gouge_audit import AuditConnector
from benchmarks.held_gouge_audit import Clearance
from benchmarks.held_gouge_audit import ElementClearance
from benchmarks.held_gouge_audit import HeldGougeAudit
from benchmarks.held_gouge_audit import HeldPathGougesBoundaryError
from benchmarks.held_gouge_audit import require_gouge_free
from benchmarks.held_gouge_boundary import Box
from benchmarks.held_gouge_boundary import DyadicGrid
from benchmarks.held_gouge_boundary import exact_box
from benchmarks.held_gouge_boundary import exact_primitive
from benchmarks.held_gouge_boundary import OffGridValueError
from benchmarks.held_gouge_surd import InvalidSurdError
from benchmarks.held_gouge_surd import Surd
from benchmarks.held_reference_geometry import ReferenceArc
from benchmarks.held_reference_geometry import ReferenceLine
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import WorldXY

R1 = Millimetre(1.0)


def _p(x: float, y: float) -> Point2[WorldXY]:
    return Point2(Millimetre(x), Millimetre(y))


def _line(a: tuple[float, float], b: tuple[float, float]) -> ReferenceLine:
    return ReferenceLine(_p(*a), _p(*b))


def _arc(a: tuple[float, float], b: tuple[float, float], c: tuple[float, float], sweep: float) -> ReferenceArc:
    return ReferenceArc(_p(*a), _p(*b), _p(*c), Radian(sweep))


def _dip_pocket() -> list[ReferenceLine | ReferenceArc]:
    """10 x 10 CCW square whose top edge dips through a concave arc.

    The arc is centred (5, 12), radius 5/2, from (6.5, 10) to (3.5, 10)
    clockwise; its chord is y = 10 and its sagitta is exactly 1/2.
    """
    sweep = -2 * math.atan2(1.5, 2.0)
    return [
        _line((0, 0), (10, 0)),
        _line((10, 0), (10, 10)),
        _line((10, 10), (6.5, 10)),
        _arc((6.5, 10), (3.5, 10), (5, 12), sweep),
        _line((3.5, 10), (0, 10)),
        _line((0, 10), (0, 0)),
    ]


def _square() -> list[ReferenceLine | ReferenceArc]:
    return [_line((0, 0), (10, 0)), _line((10, 0), (10, 10)), _line((10, 10), (0, 10)), _line((0, 10), (0, 0))]


def _disc() -> list[ReferenceLine | ReferenceArc]:
    """Radius-5 disc pocket as two convex half arcs about the origin."""
    return [_arc((5, 0), (-5, 0), (0, 0), math.pi), _arc((-5, 0), (5, 0), (0, 0), math.pi)]


def _circle(x: float, y: float, radius: float) -> AuditCircle:
    return AuditCircle(_p(x, y), Millimetre(radius))


def _connector(*points: tuple[float, float]) -> AuditConnector:
    return AuditConnector(tuple(_p(*point) for point in points))


def _one(audit: HeldGougeAudit, element: AuditCircle | AuditConnector) -> ElementClearance | None:
    """The verdict on a single element; ``None`` when it is clear."""
    report = audit.audit([element], []) if isinstance(element, AuditCircle) else audit.audit([], [element])
    return report.elements[0] if report.elements else None


def test_circle_tangent_to_the_chord_offset_gouges_by_exactly_the_sagitta() -> None:
    audit = HeldGougeAudit.build(_dip_pocket(), R1)
    # Top of the circle is on y = 9, the chord's inward r-offset.
    verdict = _one(audit, _circle(5, 7, 2))
    assert verdict is not None
    assert verdict.clearance is Clearance.GOUGE
    assert verdict.depth_mm == 0.5
    assert (verdict.primitive, verdict.primitive_kind) == (3, "concave-arc")


def test_circle_below_the_dip_by_the_sagitta_touches() -> None:
    audit = HeldGougeAudit.build(_dip_pocket(), R1)
    verdict = _one(audit, _circle(5, 6.5, 2))
    assert verdict is not None
    assert verdict.clearance is Clearance.TOUCH
    assert verdict.depth_mm == 0.0


def test_convex_arc_tie_is_touching_not_a_gouge() -> None:
    audit = HeldGougeAudit.build(_disc(), R1)
    touch = _one(audit, _circle(0, 0, 4))
    assert touch is not None and touch.clearance is Clearance.TOUCH
    assert touch.primitive_kind == "convex-arc"
    assert _one(audit, _circle(0, 0, 3.5)) is None
    deep = _one(audit, _circle(0, 0, 4.25))
    assert deep is not None and deep.clearance is Clearance.GOUGE and deep.depth_mm == 0.25


def test_straight_line_ties_are_touching() -> None:
    audit = HeldGougeAudit.build(_square(), R1)
    circle = _one(audit, _circle(5, 5, 4))
    assert circle is not None and circle.clearance is Clearance.TOUCH and circle.primitive_kind == "line"
    connector = _one(audit, _connector((2, 1), (8, 1)))
    assert connector is not None and connector.clearance is Clearance.TOUCH
    assert _one(audit, _connector((2, 1.5), (8, 1.5), (8, 8))) is None


def test_connector_under_the_dip_gouges_at_the_interior_critical_point() -> None:
    audit = HeldGougeAudit.build(_dip_pocket(), R1)
    # Arc endpoints are exactly r above the piece (a tie); the arc's bottom (5, 9.5) is r/2 above.
    verdict = _one(audit, _connector((5, 5), (3, 9), (7, 9)))
    assert verdict is not None
    assert verdict.clearance is Clearance.GOUGE
    assert (verdict.depth_mm, verdict.segment, verdict.primitive) == (0.5, 1, 3)


def test_connector_crossing_the_arc_gouges_by_the_full_radius() -> None:
    audit = HeldGougeAudit.build(_dip_pocket(), R1)
    verdict = _one(audit, _connector((5, 8), (5, 9.75)))
    assert verdict is not None and verdict.clearance is Clearance.GOUGE and verdict.depth_mm == 1.0


def test_connector_near_a_line_gouges_by_its_exact_depth() -> None:
    audit = HeldGougeAudit.build(_square(), R1)
    verdict = _one(audit, _connector((2, 0.75), (8, 0.75)))
    assert verdict is not None and verdict.clearance is Clearance.GOUGE and verdict.depth_mm == 0.25


def test_irrational_arc_end_tie_is_decided_exactly() -> None:
    # The stored end is one ulp off the circle; the true end is (0, 5) exactly, at
    # √k·(end - centre) with k not a square. (-1.5, 7) is 5/2 from it.
    end = (0.0, math.nextafter(5.0, 6.0))
    audit = HeldGougeAudit.build([_arc((5, 0), end, (0, 0), math.pi / 2)], Millimetre(2.5))
    verdict = _one(audit, _connector((-1.5, 7)))
    assert verdict is not None and verdict.clearance is Clearance.TOUCH


def test_enforcement_lists_every_gouge_with_its_depth() -> None:
    audit = HeldGougeAudit.build(_dip_pocket(), R1)
    report = audit.audit([_circle(5, 7, 2), _circle(5, 5, 1)], [_connector((2, 0.75), (8, 0.75))])
    assert (len(report.gouging_circles), len(report.gouging_connectors)) == (1, 1)
    with pytest.raises(HeldPathGougesBoundaryError, match=r"circle 0: depth 0\.5 mm against concave-arc 3") as caught:
        require_gouge_free(report)
    assert "connector 0 piece 0: depth 0.25 mm against line 0" in str(caught.value)
    require_gouge_free(audit.audit([_circle(5, 5, 1)], []))


def test_invalid_surds_and_off_grid_doubles_are_refused() -> None:
    with pytest.raises(InvalidSurdError):
        Surd.build(1, 1, -2, 1)
    with pytest.raises(OffGridValueError):
        DyadicGrid(3).integer(1 / 16)
    assert DyadicGrid(4).integer(1 / 16) == 1


def test_mixed_grid_inputs_are_injected_on_the_finest_grid() -> None:
    # A circle centre far finer than the boundary's grid still decides exactly.
    audit = HeldGougeAudit.build(_square(), R1)
    tiny = 2.0**-40  # 5 + tiny is a double; the boundary grid is integers
    verdict = _one(audit, _circle(5, 5 + tiny, 4))
    assert verdict is not None and verdict.clearance is Clearance.GOUGE and verdict.depth_mm == tiny


DYADIC = st.integers(min_value=1, max_value=1023).map(lambda n: n / 128)


@settings(max_examples=150, deadline=None)
@given(x=DYADIC, y=DYADIC, radius=st.integers(min_value=1, max_value=600).map(lambda n: n / 128))
def test_square_circle_verdict_matches_the_rational_oracle(x: float, y: float, radius: float) -> None:
    clearance_to_wall = min(x, 10 - x, y, 10 - y) - radius
    if clearance_to_wall < 0:
        return  # outside the stated precondition: the circle must lie in the pocket
    verdict = _one(HeldGougeAudit.build(_square(), R1), _circle(x, y, radius))
    expected = (clearance_to_wall > 1) - (clearance_to_wall < 1)
    assert (Clearance.CLEAR if verdict is None else verdict.clearance) == expected


@settings(max_examples=150, deadline=None)
@given(x=st.integers(-400, 400).map(lambda n: n / 128), y=st.integers(-400, 400).map(lambda n: n / 128), radius=st.integers(1, 512).map(lambda n: n / 128))
def test_disc_circle_verdict_matches_the_rational_oracle(x: float, y: float, radius: float) -> None:
    centre2 = Fraction(x) ** 2 + Fraction(y) ** 2
    reach = Fraction(4) - Fraction(radius)  # dist = 5 - |c| - R compared with r = 1
    if Fraction(5) - Fraction(radius) < 0 or centre2 > (Fraction(5) - Fraction(radius)) ** 2:
        return  # outside the precondition
    expected = -1 if reach < 0 else (centre2 < reach * reach) - (centre2 > reach * reach)
    verdict = _one(HeldGougeAudit.build(_disc(), R1), _circle(x, y, radius))
    assert (Clearance.CLEAR if verdict is None else verdict.clearance) == expected


@pytest.mark.parametrize(
    ("arc", "expected"),
    [
        (_arc((5, 0), (0, 5), (0, 0), math.pi / 2), Box(0.0, 0.0, 5.0, 5.0)),
        (_dip_pocket()[3], Box(3.5, 9.5, 6.5, 10.0)),
        (_disc()[0], Box(-5.0, 0.0, 5.0, 5.0)),
    ],
)
def test_arc_pruning_box_is_the_exact_span_box(arc: ReferenceArc, expected: Box) -> None:
    grid = DyadicGrid(8)
    assert exact_box(exact_primitive(0, arc, grid), grid) == expected
