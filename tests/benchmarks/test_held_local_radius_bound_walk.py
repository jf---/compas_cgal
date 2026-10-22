"""The ulp walk of :func:`benchmarks.held_omega_circles.local_radius_bound` starts next to its answer.

The radius bound against an offset circle is ``√D - ρ'`` (Ω* outside) or
``ρ' - √D`` (inside): a difference of nearly equal numbers when the centre is
near the support. Its estimate must be the cancellation-free
``±(D - ρ'²)/(√D + ρ')``, so the exact walk takes a bounded number of steps;
the walk alone decides the value (the largest binary64 meeting the exact
predicate), checked here against that predicate.
"""

from __future__ import annotations

import math
from fractions import Fraction
from types import SimpleNamespace

import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

import benchmarks.held_omega_circles as circles
from benchmarks.held_centre_domain_supports import OffsetCircleSupport
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_exact_centre_domain import CentreDomain
from compas_cgal.units import WorldXY
from tests.benchmarks.test_held_exact_centre_domain import concave_arc_pocket
from tests.benchmarks.test_held_exact_centre_domain import convex_arc_pocket

# Tool radius (mm) of the fixtures.
TOOL_RADIUS = Fraction(1)
# Exact predicate evaluations a walk may take: the estimate is within a few
# ulps of the answer (one rounded numerator, one rounded sqrt, one sum and one
# division: <= 4 ulps), and a walk from k ulps above costs k + 2 calls, from k
# below k + 2; 4 ulps of slack on top makes 10. The cancelled estimate
# ρ' - √D costs ~ulp(ρ')/ulp(answer) calls (2^10 and more below).
MAX_HOLDS_CALLS = 10
# Clearances (mm) of the centre from the support circle are ±2^-k, k in this
# range: answers 2^-4 to 2^-8 mm against ρ' = 11.5 / 13.5 mm, where the
# cancelled estimate is off by ~ulp(ρ') = 2^-49 mm, i.e. 2^8 to 2^12 ulps of
# the answer (kept small so the cancelled walk still terminates quickly).
CLEARANCE_EXPONENTS = (4, 8)
# The guide radius (mm) asked for: larger than every clearance drawn.
WANTED = 1.0


def _arc_support(pocket: str) -> OffsetCircleSupport:
    boundary = concave_arc_pocket() if pocket == "concave" else convex_arc_pocket()
    (support,) = CentreDomain.build(boundary, TOOL_RADIUS).supports.offset_arcs
    return support


SUPPORTS = {name: _arc_support(name) for name in ("concave", "convex")}


def _holds(support: OffsetCircleSupport, centre: ExactPoint2[WorldXY], radius: float) -> bool:
    """The exact containment predicate, stated independently of the module."""
    dx, dy = centre.minus(support.centre)
    squared, rho, value = dx * dx + dy * dy, support.radius, Fraction(radius)
    if support.sense == -1:
        return value + rho <= 0 or (value + rho) ** 2 <= squared
    return rho - value >= 0 and squared <= (rho - value) ** 2


def _centre(support: OffsetCircleSupport, angle: float, clearance: float) -> ExactPoint2[WorldXY]:
    """The binary64 point at ``ρ' - σ·clearance`` from the support centre (Ω* side for clearance > 0)."""
    distance = float(support.radius) - support.sense * clearance
    cx, cy = float(support.centre.x), float(support.centre.y)
    return ExactPoint2[WorldXY].from_binary64(cx + distance * math.cos(angle), cy + distance * math.sin(angle))


def _counted_bound(monkeypatch: pytest.MonkeyPatch, support: OffsetCircleSupport, centre: ExactPoint2[WorldXY]) -> tuple[float, int]:
    calls = [0]
    walk = circles._largest_double

    def spy(estimate: float, holds: object) -> float:
        def counted(radius: Fraction) -> bool:
            calls[0] += 1
            return holds(radius)  # type: ignore[operator]

        return walk(estimate, counted)

    monkeypatch.setattr(circles, "_largest_double", spy)
    site = SimpleNamespace(piece=SimpleNamespace(support=support))
    return circles.local_radius_bound(site, centre, WANTED), calls[0]  # type: ignore[arg-type]


@settings(max_examples=200, deadline=None)
@given(
    pocket=st.sampled_from(sorted(SUPPORTS)),
    angle=st.floats(math.pi / 8, 7 * math.pi / 8),
    exponent=st.integers(*CLEARANCE_EXPONENTS),
    side=st.sampled_from((1, -1)),
)
def test_walk_is_bounded_and_lands_on_the_largest_double(pocket: str, angle: float, exponent: int, side: int) -> None:
    support = SUPPORTS[pocket]
    centre = _centre(support, angle if support.sense == -1 else -angle, side * 2.0**-exponent)
    with pytest.MonkeyPatch.context() as monkeypatch:
        radius, calls = _counted_bound(monkeypatch, support, centre)
    assert calls <= MAX_HOLDS_CALLS
    assert _holds(support, centre, radius)
    assert not _holds(support, centre, math.nextafter(radius, math.inf))


def test_walk_near_the_support_is_bounded(monkeypatch: pytest.MonkeyPatch) -> None:
    """A clearance 2^-8 mm on a 13.5 mm support: the cancelled estimate walks thousands of ulps."""
    support = SUPPORTS["concave"]
    centre = _centre(support, 1.1, 2.0**-8)
    radius, calls = _counted_bound(monkeypatch, support, centre)
    assert calls <= MAX_HOLDS_CALLS
    assert _holds(support, centre, radius) and not _holds(support, centre, math.nextafter(radius, math.inf))
