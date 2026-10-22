"""Native site and placement batch on Ω* (``CentreDomain2.place_binary64``) against the Python reference.

The batch must build the circle :func:`benchmarks.held_omega_circles.place_local`
builds on :meth:`benchmarks.held_omega_sites.OmegaSite.decided`, field for
field, and its centre must be the correctly rounded binary64 of the exact
centre ``q + t·n``: an exact midpoint between two doubles rounds to the even
one, a value a sub-ulp fraction off a midpoint rounds to its own side. The
oracle for one-root coordinates is independent of the native code: an exact
sign of ``a - m + b√r`` over Fractions.
"""

from __future__ import annotations

import math
import struct
from fractions import Fraction

import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.held_centre_domain_native import OUTSIDE_PIECE
from benchmarks.held_centre_domain_native import PLACED
from benchmarks.held_centre_domain_native import PROBE_AT_CIRCLE_CENTRE
from benchmarks.held_centre_domain_native import native_place_binary64
from benchmarks.held_centre_domain_supports import OffsetCircleSupport
from benchmarks.held_centre_domain_supports import OffsetLineSupport
from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_exact_centre_domain import AmbiguousContactProjectionError
from benchmarks.held_exact_centre_domain import CentreDomain
from benchmarks.held_exact_centre_domain import DecidedContact
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_omega_circles import NoContainedCircleError
from benchmarks.held_omega_circles import place_local
from benchmarks.held_omega_circles import place_local_batch
from benchmarks.held_omega_sites import OmegaLoop
from benchmarks.held_omega_sites import OmegaSite
from compas_cgal import _held_centre_domain_2
from compas_cgal.units import ExactMillimetre
from tests.benchmarks.test_held_exact_centre_domain import _Pocket
from tests.benchmarks.test_held_exact_centre_domain import _point
from tests.benchmarks.test_held_exact_centre_domain import concave_arc_pocket
from tests.benchmarks.test_held_exact_centre_domain import convex_arc_pocket
from tests.benchmarks.test_held_exact_centre_domain import l_shape_pocket

F = Fraction
R = F(1)
# Generous admissible gap (mm): every fixture probe below is within it.
GAP = ExactMillimetre(F(10))
# Guide radii (mm) the property tests draw: below, at and above the fixtures'
# convex arc radius 11.5 and the tool radius.
GUIDE_RADII = (0.25, 0.5, 1.0, 2.0, 3.75, 11.5, 12.0)
# Fixture coordinates are within (−5, 25) mm.
PROBE_RANGE = (-2.0, 22.0)
# Near-tie search: guide radii 1 + k·2^-40 mm, k < 4096, keep the one whose exact
# centre x is nearest a rounding midpoint; it must be within 10^-3 ulp of one.
NEAR_TIE_SEARCH_STEPS = 4096
NEAR_TIE_SEARCH_STRIDE = 2.0**-40
NEAR_TIE_ULPS = F(1, 1000)
# Precision (bits) of the rational square root used to select near-tie inputs.
APPROXIMATION_BITS = 200


def _provenance(index: int) -> tuple[GuideRunId, GuideRunStationOrdinal, ProjectionBoundarySite]:
    return GuideRunId(0), GuideRunStationOrdinal(index), ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), F(0), 1)


def _square() -> DeclaredBoundary:
    """10 x 10 square; with r = 1 its Ω* is [1, 9]²."""
    return _Pocket(_point(0, 0)).line(0, 10).line(1, 10).line(2, 10).line(3, 10).close("square_10")


def _domain(boundary: DeclaredBoundary) -> tuple[CentreDomain, OmegaLoop]:
    domain = CentreDomain.build(boundary, R)
    return domain, OmegaLoop.build(domain)


def _decided(domain: CentreDomain, query: tuple[float, float]) -> DecidedContact:
    (decision,) = domain.decide_contacts([query], GAP).contacts
    assert decision is not None
    return decision


def _reference(loop: OmegaLoop, contact: DecidedContact, radius: float) -> object:
    """``place_local`` on the Python site, or the exception class it raises."""
    try:
        return place_local(OmegaSite.decided(loop, contact), radius, _provenance(0))
    except NoContainedCircleError:
        return NoContainedCircleError


def _batched(loop: OmegaLoop, contact: DecidedContact, radius: float) -> object:
    try:
        (circle,) = place_local_batch(loop, [contact], [radius], [_provenance(0)])
    except NoContainedCircleError:
        return NoContainedCircleError
    return circle


def _sign_one_root(a: Fraction, b: Fraction, r: Fraction) -> int:
    """Exact sign of ``a + b√r`` over Fractions (the oracle; independent of the native code)."""
    first = (a > 0) - (a < 0)
    second = 0 if r == 0 else (b > 0) - (b < 0)
    if second == 0:
        return first
    if first in (0, second):
        return second
    difference = a * a - b * b * r
    return first * ((difference > 0) - (difference < 0))


def _is_nearest(value: tuple[Fraction, Fraction, Fraction], candidate: float) -> bool:
    """Whether ``candidate`` is the correctly rounded binary64 (ties to even) of ``a + b√r``."""
    a, b, r = value
    up, down = math.nextafter(candidate, math.inf), math.nextafter(candidate, -math.inf)
    above = _sign_one_root(a - (F(candidate) + F(up)) / 2, b, r)
    below = _sign_one_root(a - (F(down) + F(candidate)) / 2, b, r)
    even = _even(candidate)
    return (above < 0 or (above == 0 and even)) and (below > 0 or (below == 0 and even))


def _even(value: float) -> bool:
    return struct.unpack("<q", struct.pack("<d", value))[0] & 1 == 0


def _approximate(a: Fraction, b: Fraction, r: Fraction) -> Fraction:
    """``a + b√r`` to ``2^-APPROXIMATION_BITS`` (selects test inputs only; never an oracle verdict)."""
    root = F(math.isqrt(r.numerator * r.denominator * 4**APPROXIMATION_BITS), r.denominator * 2**APPROXIMATION_BITS)
    return a + b * root


def test_an_exact_midpoint_centre_rounds_to_even() -> None:
    domain, loop = _domain(_square())
    # Bottom piece y = 1, normal (0, 1): c* = (5, 1 + t) with t = 2^-53 is the
    # exact midpoint of 1 and 1 + 2^-52; ties to even give 1.
    contact = _decided(domain, (5.0, 1.5))
    (placement,) = native_place_binary64(domain.native, [contact.query], [contact.piece.index], [contact.is_vertex], [2.0**-53])
    assert placement.status == PLACED and not placement.at_source
    assert placement.centre == (5.0, 1.0)
    # The rounded centre is on the support: no positive radius fits, as in Python.
    assert placement.radius == 0.0
    assert _reference(loop, contact, 2.0**-53) is NoContainedCircleError
    assert _batched(loop, contact, 2.0**-53) is NoContainedCircleError


@pytest.mark.parametrize(("offset", "expected"), [(2.0**-105, 1.0 + 2.0**-52), (-(2.0**-106), 1.0)])
def test_a_sub_ulp_offset_from_a_midpoint_rounds_to_its_side(offset: float, expected: float) -> None:
    domain, loop = _domain(_square())
    contact = _decided(domain, (5.0, 1.5))
    distance = 2.0**-53 + offset
    assert F(distance) == F(1, 2**53) + F(offset)
    (placement,) = native_place_binary64(domain.native, [contact.query], [contact.piece.index], [contact.is_vertex], [distance])
    assert placement.centre == (5.0, expected)
    assert repr(_batched(loop, contact, distance)) == repr(_reference(loop, contact, distance))


def test_a_one_root_centre_is_correctly_rounded_on_a_near_tie() -> None:
    """Concave arc about (10, -7.5), ρ' = 13.5 (Ω* outside): c* = C + (ρ' + t) w/√D, irrational."""
    domain, loop = _domain(concave_arc_pocket())
    contact = _decided(domain, (7.0, -7.5 + math.sqrt(13.5**2 - 9.0) + 0.25))
    support = contact.piece.support
    assert isinstance(support, OffsetCircleSupport) and support.sense == -1
    probe = (F(contact.query[0]), F(contact.query[1]))
    wx, wy = probe[0] - support.centre.x, probe[1] - support.centre.y
    radial = wx * wx + wy * wy
    # Search the guide radius whose exact centre x is nearest a midpoint of two
    # doubles: a near tie no enclosure a few ulps wide decides.
    best: tuple[Fraction, float] | None = None
    for step in range(NEAR_TIE_SEARCH_STEPS):
        distance = 1.0 + step * NEAR_TIE_SEARCH_STRIDE
        scale = (support.radius - support.sense * F(distance)) / radial
        value = _approximate(support.centre.x, scale * wx, radial)
        nearest = float(value)
        ulp = F(math.ulp(nearest))
        gap = min(abs(value - (F(nearest) + side * ulp / 2)) for side in (-1, 1)) / ulp
        if best is None or gap < best[0]:
            best = (gap, distance)
    assert best is not None and best[0] < NEAR_TIE_ULPS
    distance = best[1]
    (placement,) = native_place_binary64(domain.native, [contact.query], [contact.piece.index], [contact.is_vertex], [distance])
    scale = (support.radius - support.sense * F(distance)) / radial
    assert _is_nearest((support.centre.x, scale * wx, radial), placement.centre[0])
    assert _is_nearest((support.centre.y, scale * wy, radial), placement.centre[1])
    assert repr(_batched(loop, contact, distance)) == repr(_reference(loop, contact, distance))


_CACHE: dict[str, tuple[CentreDomain, OmegaLoop]] = {}


@pytest.mark.parametrize("pocket", [l_shape_pocket, concave_arc_pocket, convex_arc_pocket, _square])
@settings(max_examples=150, deadline=None)
@given(
    x=st.floats(*PROBE_RANGE, allow_nan=False, allow_infinity=False),
    y=st.floats(*PROBE_RANGE, allow_nan=False, allow_infinity=False),
    radius=st.sampled_from(GUIDE_RADII),
)
def test_the_batch_builds_place_locals_circle_field_for_field(pocket: object, x: float, y: float, radius: float) -> None:
    domain, loop = _CACHE.setdefault(pocket.__name__, _domain(pocket()))  # type: ignore[attr-defined, operator]
    try:
        (decision,) = domain.decide_contacts([(x, y)], GAP).contacts
    except AmbiguousContactProjectionError:
        # Equidistant from two ∂Ω* points: no decided contact to place.
        return
    if decision is None:
        return
    assert repr(_batched(loop, decision, radius)) == repr(_reference(loop, decision, radius))


@pytest.mark.parametrize("pocket", [l_shape_pocket, concave_arc_pocket, convex_arc_pocket])
def test_vertex_contacts_build_place_locals_circle(pocket: object) -> None:
    domain, loop = _domain(pocket())  # type: ignore[operator]
    for index, piece in enumerate(loop.pieces):
        for radius in GUIDE_RADII:
            contact = DecidedContact((0.0, 0.0), piece, True)
            assert repr(_batched(loop, contact, radius)) == repr(_reference(loop, contact, radius)), (index, radius)


def test_a_probe_projecting_beyond_its_piece_is_refused() -> None:
    domain, loop = _domain(_square())
    bottom = next(piece for piece in loop.pieces if isinstance(piece.support, OffsetLineSupport) and piece.support.normal == (F(0), F(1)))
    # The bottom piece runs from (1, 1) to (9, 1): x = 9.5 projects beyond its target, x = 1 onto its source.
    (beyond, at_source) = native_place_binary64(domain.native, [(9.5, 1.5), (1.0, 1.5)], [bottom.index] * 2, [False, False], [1.0, 1.0])
    assert beyond.status == OUTSIDE_PIECE
    assert at_source.status == PLACED and at_source.at_source


def test_a_probe_at_a_circle_centre_is_refused() -> None:
    domain, loop = _domain(convex_arc_pocket())
    arc = next(piece for piece in loop.pieces if isinstance(piece.support, OffsetCircleSupport) and piece.support.centre == _point(10, F(15, 2)))
    (placement,) = native_place_binary64(domain.native, [(10.0, 7.5)], [arc.index], [False], [1.0])
    assert placement.status == PROBE_AT_CIRCLE_CENTRE


def test_malformed_batches_fail_loud() -> None:
    domain, loop = _domain(_square())
    with pytest.raises(_held_centre_domain_2.InvalidCentreDomainInputError):
        native_place_binary64(domain.native, [(5.0, 1.5)], [len(loop.pieces)], [False], [1.0])
    with pytest.raises(_held_centre_domain_2.InvalidCentreDomainInputError):
        native_place_binary64(domain.native, [(5.0, 1.5)], [0], [False], [math.inf])
    with pytest.raises(_held_centre_domain_2.InvalidCentreDomainInputError):
        domain.native.place_binary64([5.0], [1.5, 2.0], [0], [0], [1.0])
