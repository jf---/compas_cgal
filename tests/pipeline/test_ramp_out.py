"""The tangent pencil and its furthest-first ramp-out ladder (Task 1.2, N1 / P-C / L1).

Members of the pencil at the root's entry ``P`` are internally tangent at
``P``. The ramp-out walks from the helix footprint to the root by one
furthest-first scan per member, never a bisection (engagement is not
monotone), with the structural coverage bound ``rho_next - rho_prev <= r``
prepended by the instance and the injected conditions after it in order.
"""

from __future__ import annotations

import math
from fractions import Fraction

import pytest
from compas.geometry import Polygon
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from compas_cgal.helix_units import EntryRampLength
from compas_cgal.incident_loop import INCIDENCE_GRID_EXPONENT
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.pipeline.admissibility import Admissibility
from compas_cgal.pipeline.admissibility import cap_ratio
from compas_cgal.pipeline.admissibility import centre_domain
from compas_cgal.pipeline.admissibility import entry_rise_ratio
from compas_cgal.pipeline.admissibility import loop_parameters
from compas_cgal.pipeline.errors import RampOutStuckError
from compas_cgal.pipeline.pencil import PENCIL_RUNG_TOOL_RADII
from compas_cgal.pipeline.pencil import PencilConditions
from compas_cgal.pipeline.pencil import PencilStepCondition
from compas_cgal.pipeline.pencil import RampOut
from compas_cgal.pipeline.pencil import TangentPencil
from compas_cgal.pipeline.pencil import entry_containment
from compas_cgal.stock import Stock
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

GRID = 2**INCIDENCE_GRID_EXPONENT
RECT_12X8 = Polygon([[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]])
TOOL = ToolRadius.build(1.0)
ROOT_RADIUS = 3.0
# A pocket wide enough to contain every root disk the Hypothesis test draws (radius <= 40 + r at the origin).
FREE_POCKET = entry_containment(Polygon([[-50.0, -50.0, 0.0], [50.0, -50.0, 0.0], [50.0, 50.0, 0.0], [-50.0, 50.0, 0.0]]), [], TOOL)


def _pencil(root: IncidentLoop) -> TangentPencil:
    return TangentPencil.build(root, TOOL, FREE_POCKET)


def _root(radius: float = ROOT_RADIUS, phase: float = 0.0) -> IncidentLoop:
    return IncidentLoop.build(Point2[WorldXY].build(0.0, 0.0), GuideRadius.build(radius), Radian(phase))


def _accept(_stock: Stock, _origin: IncidentLoop, _candidate: IncidentLoop, /) -> bool:
    return True


def _refuse(_stock: Stock, _origin: IncidentLoop, _candidate: IncidentLoop, /) -> bool:
    return False


def _run(ramp: RampOut, stock: Stock) -> list[IncidentLoop]:
    """Every member the ramp admits, root last; the stock is not depleted (the conditions here ignore it)."""
    members = [ramp.footprint]
    while not ramp.is_root(members[-1]):
        members.append(ramp.next_member(stock, members[-1]))
    return members


@settings(max_examples=100, deadline=None)
@given(radius=st.floats(min_value=1.05, max_value=40.0), phase=st.floats(min_value=-math.pi, max_value=math.pi), fraction=st.floats(min_value=0.0, max_value=1.0))
def test_every_member_is_internally_tangent_to_the_root_at_its_entry(radius: float, phase: float, fraction: float) -> None:
    root = _root(radius, phase)
    pencil = _pencil(root)
    steps = max(1, round(fraction * pencil.root_steps))
    member = pencil.member(steps)
    assert member.certifies_incidence()
    assert member.entry == root.entry
    # c_k - P is parallel to c_root - P with the same sense, and |c_k - c_root| == rho_root - rho_k, exactly.
    px, py = Fraction(root.entry.x), Fraction(root.entry.y)
    kx, ky = Fraction(member.centre.x) - px, Fraction(member.centre.y) - py
    rx, ry = Fraction(root.centre.x) - px, Fraction(root.centre.y) - py
    assert kx * ry - ky * rx == 0 and kx * rx + ky * ry > 0
    gap = Fraction(root.radius.value) - Fraction(member.radius.value)
    dx, dy = Fraction(member.centre.x) - Fraction(root.centre.x), Fraction(member.centre.y) - Fraction(root.centre.y)
    assert gap >= 0 and dx * dx + dy * dy == gap * gap
    assert pencil.steps_of(member) == steps
    assert pencil.member(pencil.root_steps) == root


def test_the_pencil_step_bound_is_exact_at_r() -> None:
    """P-C at its boundary: a rise of exactly r is admitted, one grid step more is refused."""
    pencil = _pencil(_root())
    bound = PencilStepCondition(TOOL)
    reach = Fraction(TOOL.value) * GRID / pencil.hypotenuse
    assert reach.denominator == 1
    origin = pencil.member(1)
    stock = Stock(RECT_12X8)
    assert bound(stock, origin, pencil.member(1 + int(reach)))
    assert not bound(stock, origin, pencil.member(2 + int(reach)))


def test_the_ramp_out_prepends_the_pencil_step_bound() -> None:
    """Injected `(accept,)` only: every rise is still <= r, the last member is the root."""
    root = _root()
    ramp = RampOut.build(_pencil(root), TOOL, (_accept,))
    members = _run(ramp, Stock(RECT_12X8))
    assert members[-1] == root
    rises = [Fraction(b.radius.value) - Fraction(a.radius.value) for a, b in zip(members, members[1:])]
    assert rises and all(0 < rise <= Fraction(TOOL.value) for rise in rises)


def test_the_ramp_out_is_furthest_first_and_never_bisects() -> None:
    """One scan: candidates are asked from the furthest P-C-admissible rung down, one rung apart, and the first admitted is returned."""
    root = _root()
    pencil = _pencil(root)
    asked: list[int] = []
    admit_at = 5

    def scripted(_stock: Stock, _origin: IncidentLoop, candidate: IncidentLoop, /) -> bool:
        asked.append(pencil.steps_of(candidate))
        return len(asked) == admit_at

    ramp = RampOut.build(pencil, TOOL, (scripted,))
    got = ramp.next_member(Stock(RECT_12X8), ramp.footprint)
    footprint_steps = pencil.steps_of(ramp.footprint)
    reach = int(Fraction(TOOL.value) * GRID // pencil.hypotenuse)
    first_admissible = max(step for step in range(pencil.root_steps, footprint_steps, -ramp.rung_steps) if step <= footprint_steps + reach)
    assert asked[0] == first_admissible
    assert all(a - b == ramp.rung_steps for a, b in zip(asked, asked[1:]))
    assert pencil.steps_of(got) == asked[-1] and len(asked) == admit_at


def test_the_ramp_rung_is_the_declared_fraction_of_the_coverage_window_on_the_root_grid() -> None:
    pencil = _pencil(_root())
    ramp = RampOut.build(pencil, TOOL, ())
    assert ramp.rung_steps == math.ceil(PENCIL_RUNG_TOOL_RADII * Fraction(TOOL.value) * GRID / pencil.hypotenuse)


def test_no_admissible_ramp_member_raises_by_name() -> None:
    ramp = RampOut.build(_pencil(_root()), TOOL, (_refuse,))
    with pytest.raises(RampOutStuckError):
        ramp.next_member(Stock(RECT_12X8), ramp.footprint)


def test_the_ramp_out_conditions_run_in_order_cheapest_first() -> None:
    """N1 order: a candidate the bound refuses never reaches the junction; one the junction refuses never reaches the cap."""
    pencil = _pencil(_root())
    seen_junction: list[IncidentLoop] = []
    passed_junction: list[IncidentLoop] = []
    seen_cap: list[IncidentLoop] = []

    def junction(_stock: Stock, origin: IncidentLoop, candidate: IncidentLoop, /) -> bool:
        seen_junction.append(candidate)
        assert Fraction(candidate.radius.value) - Fraction(origin.radius.value) <= Fraction(TOOL.value)
        admitted = pencil.steps_of(candidate) % 3 == 0
        if admitted:
            passed_junction.append(candidate)
        return admitted

    def cap(_stock: Stock, _origin: IncidentLoop, candidate: IncidentLoop, /) -> bool:
        seen_cap.append(candidate)
        return True

    ramp = RampOut.build(pencil, TOOL, (junction, cap))
    _run(ramp, Stock(RECT_12X8))
    assert seen_junction and seen_cap == passed_junction


def _admissibility(cap_deg: float) -> Admissibility:
    cap = Radian(math.radians(cap_deg))
    return Admissibility.build(TOOL, cap_ratio(cap), centre_domain(RECT_12X8, [], TOOL), entry_rise_ratio=entry_rise_ratio(cap, EntryRampLength.default()))


def test_the_pencil_conditions_are_the_junction_then_the_full_cap() -> None:
    admissibility = _admissibility(120.0)
    conditions = PencilConditions(admissibility, clockwise=True)
    assert conditions.ordered == (conditions.junction_admissible, conditions.loop_admissible)


def test_the_ramp_out_asks_every_loop_station_and_never_the_at_rest_ring(monkeypatch: pytest.MonkeyPatch) -> None:
    """The full cap asks k/45 for k = 0..44 (station 0 included: no ring supplies it); the at-rest ring is never asked (AR scope)."""
    from compas_cgal import _stock_2

    admissibility = _admissibility(120.0)
    asked: list[list[float]] = []
    native = _stock_2.engagement_circle_exceeded

    def recorded(*args: object) -> bool:
        asked.append(list(args[7]))  # type: ignore[call-overload]
        return bool(native(*args))  # type: ignore[arg-type]

    def ring(*_args: object) -> bool:
        raise AssertionError("the ramp-out asked the at-rest ring")

    monkeypatch.setattr(_stock_2, "engagement_circle_exceeded", recorded)
    monkeypatch.setattr(Admissibility, "loop_within_cap_at_rest", ring)
    root = _root()
    ramp = RampOut.build(_pencil(root), TOOL, PencilConditions(admissibility, clockwise=True).ordered)
    stock = Stock(RECT_12X8)
    stock.subtract_circle_sweep_local(ramp.footprint.centre.x, ramp.footprint.centre.y, ramp.footprint.radius.value, TOOL.value)
    ramp.next_member(stock, ramp.footprint)
    assert loop_parameters() in asked
    assert all(parameters in ([loop_parameters()[1]], loop_parameters()) for parameters in asked)
