"""Fixtures that reach the certificate's tangency rules an injected defect used to survive (RV7 I3.3, RV10 Important 2).

- **Fast path 3** (line_sweep_2.cpp, two tangent vertex roots on one line): a
  polygon capsule whose flank x = 2 is collinear with the left offset line of
  a pocket with concave arcs top and bottom, and spans both of its one-root
  kinks (y = −7.5 + √(13.5² − 9²) ≈ 2.56 and 27.5 − √(13.5² − 9²) ≈ 17.44), so the two target
  vertex circles touch the flank inside its range. The sort compares two events
  only when they end up adjacent: with concave arcs the declared circles cross
  x = 2 outside the kinks (y ≈ 2.1 and 17.9), so nothing lies between the two
  tangencies, and EXACT_ONLY orders the pair exactly (counter
  ``vertex_tangent_pairs``). The convex lens does not reach it: its declared
  circles cross x = 2 between the kinks.
- **Support-only tangencies** (lazy_crossings_2.cpp): an arc–arc G1 junction at
  a rational kink. The bottom of P runs along the arc about (10, 7.5), ρ = 12.5,
  from (0, 0) to (10, −5) and on, tangentially, along the arc about (10, 9.5),
  ρ = 14.5, to (20, −1). The Ω* kink is w = (10, −4); the connectors' inner
  supports (radii 10.5 and 12.5) touch at p = (10, −3). With the cuts below
  y = 1 removed, the exposed edge on the first inner support ends at p.

Both configurations are pinned by exposure digest and checked LAZY == EXACT_ONLY,
and every witness is re-checked independently of the certificate.
"""

from __future__ import annotations

from collections.abc import Iterator
from fractions import Fraction
from functools import cache

import pytest

from benchmarks.held_coverage_target import exact_coverage_target
from benchmarks.held_exact_centre_domain import CentreDomain as ExactCentreDomain
from compas_cgal.coverage_certificate.circle_core import FilterMode
from compas_cgal.coverage_certificate.native_certificate import CertificateReport
from compas_cgal.coverage_certificate.native_certificate import run_certificate
from compas_cgal.coverage_certificate.sweeps import BoundaryConnectorSweep
from compas_cgal.coverage_certificate.target import CentreBoundarySite
from compas_cgal.coverage_certificate.target import CentreDomain
from compas_cgal.coverage_certificate.witness import ExposedPieceKind
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY
from tests.benchmarks.test_coverage_certificate_exact_connectors import LENS_GRID
from tests.benchmarks.test_coverage_certificate_exact_connectors import R
from tests.benchmarks.test_coverage_certificate_exact_connectors import TOOL
from tests.benchmarks.test_coverage_certificate_exact_connectors import UNUSED_PATH
from tests.benchmarks.test_coverage_certificate_exact_connectors import Fixture
from tests.benchmarks.test_coverage_certificate_exact_connectors import _check_witnesses
from tests.benchmarks.test_coverage_certificate_exact_connectors import _digest
from tests.benchmarks.test_coverage_certificate_exact_connectors import _grid
from tests.benchmarks.test_coverage_certificate_exact_connectors import _point
from tests.benchmarks.test_coverage_certificate_exact_connectors import _Pocket
from tests.benchmarks.test_coverage_certificate_exact_connectors import _sweeps
from tests.benchmarks.test_coverage_certificate_exact_connectors import _whole_boundary

F = Fraction
# The path polygon: its left edge x = 1 is the pocket's left offset line,
# extended past both kinks.
FLANK_PATH = CentreDomain.build([Point2[WorldXY].build(x, y) for x, y in ((1.0, -4.0), (19.0, -4.0), (19.0, 24.0), (1.0, 24.0))])
# Interior cuts kept right of x = 3.5, so the capsule flank x = 2 is exposed.
FLANK_CUTS = [c for c in LENS_GRID if c[0] > 3.5 and 3.0 < c[1] < 17.0]
G1_GRID = _grid(2.0, -2.0, 2.5, 7, 11)
# Cuts kept above y = 1, so the inner supports meet the gap at p = (10, -3).
G1_CUTS = [c for c in G1_GRID if c[1] > 1.0]
DIGESTS = {
    "concave_lens_polygon_flank": "5f136c5ad9c0c9b8",
    "g1_arcs_bottom_gap": "a0151c9602e9e144",
}


@cache
def concave_lens() -> Fixture:
    """Concave arcs top and bottom (about (10, -7.5) and (10, 27.5), ρ = 12.5): one-root kinks at both ends of each line piece."""
    boundary = _Pocket(_point(0, 0)).arc(_point(20, 0), F(1, 2), 0, -1).line(1, 20).arc(_point(0, 20), F(1, 2), 2, -1).line(3, 20).close("concave_lens")
    domain = ExactCentreDomain.build(boundary, R)
    return boundary, domain, exact_coverage_target(boundary, domain)


@cache
def g1_arcs() -> Fixture:
    """Bottom of P: two CCW arcs meeting G1 at the rational point (10, -5); walls and a top close it."""
    boundary = _Pocket(_point(0, 0)).arc(_point(10, -5), F(-1, 2), 0, 1).arc(_point(20, -1), F(0), 0, 1).line(1, 21).line(2, 20).line(3, 20).close("g1_arcs")
    domain = ExactCentreDomain.build(boundary, R)
    return boundary, domain, exact_coverage_target(boundary, domain)


@pytest.fixture(autouse=True, scope="module")
def _release_native_domains() -> Iterator[None]:
    yield
    g1_arcs.cache_clear()
    concave_lens.cache_clear()


def _both(
    target: object, cuts: list[tuple[float, float, float]], path: CentreDomain, connectors: list[BoundaryConnectorSweep], exact_connectors: list[object]
) -> tuple[CertificateReport, CertificateReport]:
    runs = tuple(
        run_certificate(path, TOOL, _sweeps(cuts), connectors, exact_target=target, exact_connectors=exact_connectors, filter_mode=mode, audit=True)  # type: ignore[arg-type]
        for mode in (FilterMode.LAZY, FilterMode.EXACT_ONLY)
    )
    assert runs[0].canonical_text() == runs[1].canonical_text()
    return runs


def test_two_tangent_vertex_roots_on_a_polygon_flank_are_ordered_exactly() -> None:
    # No exact connectors: their arc supports would cross x = 2 between the kinks.
    _, domain, target = concave_lens()
    left_edge = BoundaryConnectorSweep.build(FLANK_PATH, CentreBoundarySite.build(0, 3, F(0)), CentreBoundarySite.build(0, 0, F(0)))
    lazy, exact = _both(target, FLANK_CUTS, FLANK_PATH, [left_edge], [])
    assert exact.counters["vertex_tangent_pairs"] > 0
    assert any(p.piece_kind is ExposedPieceKind.CAPSULE_FLANK for p in lazy.exposed)
    assert _digest(lazy) == DIGESTS["concave_lens_polygon_flank"]


def test_an_arc_arc_g1_junction_at_a_rational_kink_cuts_the_inner_supports_at_the_kink() -> None:
    _, domain, target = g1_arcs()
    lazy, _ = _both(target, G1_CUTS, UNUSED_PATH, [], _whole_boundary(target))
    _check_witnesses(domain, lazy, G1_CUTS)
    assert lazy.exposed
    assert _digest(lazy) == DIGESTS["g1_arcs_bottom_gap"]
