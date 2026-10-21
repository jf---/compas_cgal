"""Contract tests for the capsule and target pieces of the coverage certificate (phase 2).

Connector capsules enter as two end disks and a rectangle whose flanks lie at
the irrational offset r sqrt(L); target pieces are the edge-rectangle flanks and
vertex circles. These tests hold the lazy front and the tower back end to their
contract:

- forced EXACT_ONLY is byte-identical to LAZY with capsules, the audit oracle on;
- generic doubles: every witness re-checks exactly with ``Fraction`` here,
  against cuts AND connector capsules, independently of the native code, and the
  cross-base comparison (two flanks of different irrational edges) is exercised;
- a circle nearly tangent to an irrational flank reaches the exact fallback and
  is decided transversal (not tangent);
- connector caps tangent to their own flanks are exact tangencies.
"""

from __future__ import annotations

import math
import random
from collections.abc import Sequence
from fractions import Fraction

import pytest

from compas_cgal.coverage_certificate.certificate import CoverageCertificate
from compas_cgal.coverage_certificate.circle_core import FilterMode
from compas_cgal.coverage_certificate.errors import UncoveredPocketError
from compas_cgal.coverage_certificate.native_certificate import CertificateReport
from compas_cgal.coverage_certificate.native_certificate import run_certificate
from compas_cgal.coverage_certificate.sweeps import BoundaryConnectorSweep
from compas_cgal.coverage_certificate.sweeps import FullCircleSweep
from compas_cgal.coverage_certificate.target import CentreBoundarySite
from compas_cgal.coverage_certificate.target import CentreDomain
from compas_cgal.coverage_certificate.witness import ExposedPieceKind
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

Cut = tuple[float, float, float]
XY = tuple[Fraction, Fraction]


def _domain(points: Sequence[tuple[float, float]]) -> CentreDomain:
    return CentreDomain.build([Point2[WorldXY].build(x, y) for x, y in points])


def _sweeps(cuts: Sequence[Cut]) -> list[FullCircleSweep]:
    return [FullCircleSweep.build(Point2[WorldXY].build(x, y), GuideRadius.build(big_r)) for x, y, big_r in cuts]


def _both_modes(domain: CentreDomain, r: float, cuts: Sequence[Cut], connectors: Sequence[BoundaryConnectorSweep]) -> tuple[CertificateReport, CertificateReport]:
    tool = ToolRadius.build(r)
    lazy = run_certificate(domain, tool, _sweeps(cuts), connectors, filter_mode=FilterMode.LAZY, audit=True)
    exact = run_certificate(domain, tool, _sweeps(cuts), connectors, filter_mode=FilterMode.EXACT_ONLY, audit=True)
    assert lazy.canonical_text() == exact.canonical_text()
    return lazy, exact


# --- independent Fraction oracle -------------------------------------------------


def _dist2_to_segment(p: XY, a: XY, b: XY) -> Fraction:
    dx, dy = b[0] - a[0], b[1] - a[1]
    t = min(max(((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / (dx * dx + dy * dy), Fraction(0)), Fraction(1))
    ex, ey = p[0] - a[0] - t * dx, p[1] - a[1] - t * dy
    return ex * ex + ey * ey


def _in_polygon(p: XY, loop: Sequence[XY]) -> bool:
    inside = False
    for i in range(len(loop)):
        a, b = loop[i], loop[(i + 1) % len(loop)]
        if (a[1] > p[1]) != (b[1] > p[1]):
            x = a[0] + (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
            if p[0] < x:
                inside = not inside
    return inside


def _connector_pieces(loop: Sequence[XY], c: BoundaryConnectorSweep) -> list[tuple[XY, XY]]:
    n = len(loop)

    def point(segment: int, tau: Fraction) -> XY:
        a, b = loop[segment], loop[(segment + 1) % n]
        return (a[0] + tau * (b[0] - a[0]), a[1] + tau * (b[1] - a[1]))

    s0, t0, s1, t1 = c.start.segment, c.start.parameter, c.end.segment, c.end.parameter
    if s0 == s1 and t0 < t1:
        return [(point(s0, t0), point(s1, t1))]
    pieces = [(point(s0, t0), loop[(s0 + 1) % n])]
    s = (s0 + 1) % n
    while s != s1:
        pieces.append((loop[s], loop[(s + 1) % n]))
        s = (s + 1) % n
    pieces.append((loop[s1], point(s1, t1)))
    return pieces


def _assert_witnesses_exact(report: CertificateReport, loop: Sequence[XY], r: float, cuts: Sequence[Cut], connectors: Sequence[BoundaryConnectorSweep]) -> None:
    rr = Fraction(r) ** 2
    capsules = [piece for c in connectors for piece in _connector_pieces(loop, c)]
    for e in report.exposed:
        u = (Fraction(e.uncovered[0]), Fraction(e.uncovered[1]))
        in_target = _in_polygon(u, loop) or any(_dist2_to_segment(u, loop[i], loop[(i + 1) % len(loop)]) <= rr for i in range(len(loop)))
        assert in_target, f"{e} uncovered point is outside T"
        for x, y, big_r in cuts:
            d2 = (u[0] - Fraction(x)) ** 2 + (u[1] - Fraction(y)) ** 2
            inner = max(Fraction(big_r) - Fraction(r), Fraction(0))
            assert not inner**2 <= d2 <= (Fraction(big_r) + Fraction(r)) ** 2, f"{e} uncovered point is in a cut"
        for a, b in capsules:
            assert _dist2_to_segment(u, a, b) > rr, f"{e} uncovered point is in a connector capsule"


# --- fixtures --------------------------------------------------------------------

TRIANGLE = [(0.0, 0.0), (8.0, 0.0), (4.0, 7.0)]  # edges 1 and 2 have L = 65: irrational flanks
EDGE_GAP_CUTS = [(0.0, 0.0, 3.25), (0.0, 0.0, 1.25), (0.0, 0.0, 0.25), (8.0, 0.0, 2.875), (8.0, 0.0, 0.875), (4.0, 7.0, 3.5), (4.0, 7.0, 1.5), (4.0, 7.0, 0.5), (4.0, 2.375, 0.75)]


def test_connector_flank_coincides_exactly_and_lazy_equals_exact() -> None:
    domain = _domain(TRIANGLE)
    connector = BoundaryConnectorSweep.build(domain, CentreBoundarySite.build(0, 0, Fraction(1, 4)), CentreBoundarySite.build(0, 0, Fraction(3, 4)))
    lazy, _ = _both_modes(domain, 1.0, EDGE_GAP_CUTS, [connector])
    assert not lazy.exposed
    assert lazy.counters["coincident_line_supports"] >= 1, "the connector flank shares the target flank's support"
    # The connector's end disks are tangent to its own flanks: exact double roots.
    assert lazy.counters["line_tangencies"] >= 4


def test_connectors_on_irrational_edges_lazy_equals_exact() -> None:
    # Connectors along the two L = 65 edges and across vertex (4, 7): their flanks
    # are towers over Q(sqrt(65)); the gap fixture stays open, so witnesses exist.
    domain = _domain(TRIANGLE)
    connectors = [
        BoundaryConnectorSweep.build(domain, CentreBoundarySite.build(0, 1, Fraction(1, 3)), CentreBoundarySite.build(0, 2, Fraction(1, 5))),
        BoundaryConnectorSweep.build(domain, CentreBoundarySite.build(0, 2, Fraction(1, 2)), CentreBoundarySite.build(0, 0, Fraction(1, 8))),
    ]
    cuts = EDGE_GAP_CUTS[:-1]  # without the closing cut: the centre gap is open
    lazy, _ = _both_modes(domain, 1.0, cuts, connectors)
    assert lazy.exposed
    loop = [(Fraction(x), Fraction(y)) for x, y in TRIANGLE]
    _assert_witnesses_exact(lazy, loop, 1.0, cuts, connectors)


def _star_polygon(rng: random.Random, n: int) -> list[tuple[float, float]]:
    angles = sorted(rng.uniform(0.0, 2 * math.pi) for _ in range(n))
    return [(10.0 + rng.uniform(6.0, 10.0) * math.cos(a), 10.0 + rng.uniform(6.0, 10.0) * math.sin(a)) for a in angles]


def test_generic_double_capsule_witnesses_are_exactly_uncovered() -> None:
    rng = random.Random(20261021)
    outer = _star_polygon(rng, 9)
    domain = _domain(outer)
    r = rng.uniform(0.6, 0.9)
    cuts: list[Cut] = [(rng.uniform(2.0, 18.0), rng.uniform(2.0, 18.0), rng.uniform(0.3, 2.5)) for _ in range(40)]
    connectors = []
    for _ in range(8):
        s0 = rng.randrange(9)
        t0 = Fraction(rng.random())
        s1 = (s0 + rng.randrange(2)) % 9
        t1 = Fraction(rng.random())
        if s1 == s0 and t1 <= t0:
            t0, t1 = t1, t0
        if (s0, t0) == (s1, t1):
            continue
        connectors.append(BoundaryConnectorSweep.build(domain, CentreBoundarySite.build(0, s0, t0), CentreBoundarySite.build(0, s1, t1)))
    lazy, exact = _both_modes(domain, r, cuts, connectors)
    assert lazy.exposed
    assert lazy.counters["line_pieces"] > 0
    assert exact.counters["cross_base_separations"] > 0, "two flanks over different irrational bases were compared"
    loop = [(Fraction(x), Fraction(y)) for x, y in outer]
    _assert_witnesses_exact(lazy, loop, r, cuts, connectors)
    assert any(e.piece_kind is ExposedPieceKind.TARGET_BOUNDARY for e in lazy.exposed)


# Triangle with edge (0,0)-(3,1): L = 10, inner flank -x + 3y = sqrt(10). A cut at
# (0, cy) with outer radius 3/2 and cy = sqrt(10) (1 + 3/2) / 3 rounded to binary64
# is within one rounding of tangency to that flank: the certified discriminant
# enclosure straddles zero, the exact tower decides two transversal crossings.
NEAR_FLANK_DOMAIN = [(0.0, 0.0), (3.0, 1.0), (-1.0, 4.0)]
NEAR_FLANK_CUT = (0.0, math.sqrt(10.0) * 2.5 / 3.0, 0.5)


def test_near_tangent_flank_circle_pair_hits_the_exact_fallback() -> None:
    domain = _domain(NEAR_FLANK_DOMAIN)
    lazy, _ = _both_modes(domain, 1.0, [NEAR_FLANK_CUT], [])
    assert lazy.counters["circle_line_exact"] >= 1, "the filter could not decide the near-tangent pair"
    assert lazy.counters["line_tangencies"] == 0, "binary64 cannot be exactly tangent to an irrational flank"
    loop = [(Fraction(x), Fraction(y)) for x, y in NEAR_FLANK_DOMAIN]
    _assert_witnesses_exact(lazy, loop, 1.0, [NEAR_FLANK_CUT], [])


def test_uncovered_target_raises_with_target_boundary_witnesses() -> None:
    domain = _domain(NEAR_FLANK_DOMAIN)
    with pytest.raises(UncoveredPocketError) as caught:
        CoverageCertificate.build(domain, ToolRadius.build(1.0), _sweeps([NEAR_FLANK_CUT]))
    assert any(w.piece_kind is ExposedPieceKind.TARGET_BOUNDARY for w in caught.value.witnesses)
