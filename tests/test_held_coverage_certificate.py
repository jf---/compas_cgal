"""Contract tests for the exact Held coverage certificate (docs/held_coverage_certificate.md).

The certificate decides whether the union of the swept regions covers the
machinable target T = centre domain (+) tool disk. These fixtures are tiny and
use only binary64-exact (dyadic) coordinates, so every expected outcome is an
exact statement, re-checked here with `Fraction` arithmetic as an independent
oracle.

Fixture: centre domain = CCW triangle (0,0), (8,0), (4,7); tool radius r = 1.
At each vertex a concentric stack of full-circle cuts sweeps a solid disk of
reach rho (outer annulus, middle annulus, inner disk with R <= r). The stacks
share boundary circles exactly (outer edge of one = inner edge of the next),
which exercises coincident supports with opposite sides.

- Three stacks of reach 17/4 leave a curvilinear triangle around the
  circumcentre (4, 33/14), whose distance 4.643 exceeds 17/4 from all three.
- A fourth cut (R = 3/4 at (4, 19/8), solid disk of radius 7/4) closes it.
- Shrinking the stack at (8,0) to reach 31/8 opens a gap on the target's own
  boundary, the outer flank y = -1 of edge (0,0)-(8,0), for x in
  (sqrt(273)/4, 8 - sqrt(897)/8) = (4.1307, 4.2563).
- A boundary connector along that edge, parameters 1/4 to 3/4, closes the
  boundary gap with a capsule whose outer flank lies exactly on the target's.
"""

from __future__ import annotations

from collections.abc import Sequence
from fractions import Fraction

import pytest

from compas_cgal.coverage_certificate.certificate import CoverageCertificate
from compas_cgal.coverage_certificate.errors import UncoveredPocketError
from compas_cgal.coverage_certificate.sweeps import BoundaryConnectorSweep
from compas_cgal.coverage_certificate.sweeps import FullCircleSweep
from compas_cgal.coverage_certificate.target import CentreBoundarySite
from compas_cgal.coverage_certificate.target import CentreDomain
from compas_cgal.coverage_certificate.witness import ExposedPieceKind
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

TOOL_RADIUS = Fraction(1)
TRIANGLE = ((Fraction(0), Fraction(0)), (Fraction(8), Fraction(0)), (Fraction(4), Fraction(7)))

# (centre, path radius R) for every full-circle cut, all dyadic so binary64-exact.
GAP_STACKS = (
    ((0.0, 0.0), (3.25, 1.25, 0.25)),  # reach 17/4
    ((8.0, 0.0), (3.25, 1.25, 0.25)),  # reach 17/4
    ((4.0, 7.0), (3.25, 1.25, 0.25)),  # reach 17/4
)
CLOSING_CUT = ((4.0, 2.375), (0.75,))  # solid disk of radius 7/4
EDGE_GAP_STACKS = (
    ((0.0, 0.0), (3.25, 1.25, 0.25)),  # reach 17/4
    ((8.0, 0.0), (2.875, 0.875)),  # reach 31/8: opens the flank gap
    ((4.0, 7.0), (3.5, 1.5, 0.5)),  # reach 9/2: keeps edge (8,0)-(4,7) covered
    CLOSING_CUT,
)

# Exact bounds of the expected gaps (squared distances, so no square roots).
CIRCUMCENTRE = (Fraction(4), Fraction(33, 14))
EDGE_GAP_X = (Fraction(4), Fraction(43, 10))  # brackets (4.1307, 4.2563)
EDGE_GAP_Y = (Fraction(-1), Fraction(-7, 10))  # top vertex of the gap at y = -0.709


def _cuts(stacks: Sequence[tuple[tuple[float, float], tuple[float, ...]]]) -> list[FullCircleSweep]:
    return [FullCircleSweep.build(Point2[WorldXY].build(*centre), GuideRadius.build(radius)) for centre, radii in stacks for radius in radii]


def _domain() -> CentreDomain:
    return CentreDomain.build([Point2[WorldXY].build(float(x), float(y)) for x, y in TRIANGLE])


def _certify(stacks: Sequence[tuple[tuple[float, float], tuple[float, ...]]], connectors: Sequence[BoundaryConnectorSweep] = ()) -> CoverageCertificate:
    return CoverageCertificate.build(_domain(), ToolRadius.build(float(TOOL_RADIUS)), _cuts(stacks), connectors)


def _squared_distance_to_segment(p: tuple[Fraction, Fraction], a: tuple[Fraction, Fraction], b: tuple[Fraction, Fraction]) -> Fraction:
    dx, dy = b[0] - a[0], b[1] - a[1]
    t = min(max(((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / (dx * dx + dy * dy), Fraction(0)), Fraction(1))
    ex, ey = p[0] - a[0] - t * dx, p[1] - a[1] - t * dy
    return ex * ex + ey * ey


def _in_target(p: tuple[Fraction, Fraction]) -> bool:
    edges = [(TRIANGLE[k], TRIANGLE[(k + 1) % 3]) for k in range(3)]
    inside = all((b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]) >= 0 for a, b in edges)
    return inside or min(_squared_distance_to_segment(p, a, b) for a, b in edges) <= TOOL_RADIUS**2


def _swept_by_cut(p: tuple[Fraction, Fraction], centre: tuple[float, float], radius: float) -> bool:
    big_r = Fraction(radius)
    d2 = (p[0] - Fraction(centre[0])) ** 2 + (p[1] - Fraction(centre[1])) ** 2
    inner = max(big_r - TOOL_RADIUS, Fraction(0))
    return inner**2 <= d2 <= (big_r + TOOL_RADIUS) ** 2


def _assert_exactly_uncovered(error: UncoveredPocketError, stacks: Sequence[tuple[tuple[float, float], tuple[float, ...]]]) -> None:
    assert error.witnesses, "an uncovered certificate must name at least one exposed edge"
    for witness in error.witnesses:
        p = witness.uncovered_point
        assert isinstance(p[0], Fraction) and isinstance(p[1], Fraction), "the uncovered point is an exact rational"
        assert _in_target(p), f"witness {witness} lies outside the machinable target"
        assert not any(_swept_by_cut(p, centre, radius) for centre, radii in stacks for radius in radii), f"witness {witness} is swept"


def test_three_annulus_gap_raises_with_exact_witnesses_inside_the_gap() -> None:
    with pytest.raises(UncoveredPocketError) as caught:
        _certify(GAP_STACKS)
    _assert_exactly_uncovered(caught.value, GAP_STACKS)
    for witness in caught.value.witnesses:
        x, y = witness.uncovered_point
        assert (x - CIRCUMCENTRE[0]) ** 2 + (y - CIRCUMCENTRE[1]) ** 2 < 1, f"{witness} is not in the circumcentre gap"
        assert witness.piece_kind is not ExposedPieceKind.TARGET_BOUNDARY


def test_fourth_cut_closes_the_gap_and_the_certificate_holds() -> None:
    certificate = _certify((*GAP_STACKS, CLOSING_CUT))
    assert isinstance(certificate, CoverageCertificate)
    assert certificate.exposed_edge_count == 0


def test_gap_on_the_target_boundary_is_reported_as_an_exposed_target_edge() -> None:
    with pytest.raises(UncoveredPocketError) as caught:
        _certify(EDGE_GAP_STACKS)
    _assert_exactly_uncovered(caught.value, EDGE_GAP_STACKS)
    for witness in caught.value.witnesses:
        x, y = witness.uncovered_point
        assert EDGE_GAP_X[0] < x < EDGE_GAP_X[1] and EDGE_GAP_Y[0] <= y < EDGE_GAP_Y[1], f"{witness} is not in the flank gap"
    assert any(witness.piece_kind is ExposedPieceKind.TARGET_BOUNDARY for witness in caught.value.witnesses)


def test_boundary_connector_closes_the_target_edge_gap_on_a_coincident_flank() -> None:
    domain = _domain()
    connector = BoundaryConnectorSweep.build(
        domain,
        CentreBoundarySite.build(loop=0, segment=0, parameter=Fraction(1, 4)),
        CentreBoundarySite.build(loop=0, segment=0, parameter=Fraction(3, 4)),
    )
    certificate = CoverageCertificate.build(domain, ToolRadius.build(float(TOOL_RADIUS)), _cuts(EDGE_GAP_STACKS), [connector])
    assert certificate.exposed_edge_count == 0
