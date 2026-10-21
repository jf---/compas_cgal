"""Contract tests: the coverage certificate on the exact target T = Omega* (+) D_r.

The target is cut in completely (declared boundary supports, one-root vertex
circles), so ``CoverageCertificate.build`` decides both ways on Omega* targets
and never raises ``UnmodelledCertificatePieceError``. Every witness is
re-checked here independently of the certificate's own target code:

- in T: the exact squared distance to Omega* from the native centre-domain
  module (``CentreDomain.squared_distance_to``, a different implementation) is
  at most r^2;
- outside every sweep: Fraction arithmetic against every cut annulus and every
  connector capsule.

Fixtures are Pythagorean declared boundaries with r = 1, so Omega* has one-root
vertices wherever an offset arc meets an offset line.
"""

from __future__ import annotations

import hashlib
import random
from collections.abc import Iterator
from collections.abc import Sequence
from fractions import Fraction
from functools import cache

import pytest

from compas_cgal import _coverage_certificate_2 as native

from benchmarks.held_coverage_target import exact_coverage_target
from benchmarks.held_declared_boundary import DeclaredArc
from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_declared_boundary import DeclaredLine
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_boundary import PythagoreanDirection
from benchmarks.held_exact_centre_domain import CentreDomain as ExactCentreDomain
from compas_cgal.coverage_certificate.certificate import CoverageCertificate
from compas_cgal.coverage_certificate.circle_core import FilterMode
from compas_cgal.coverage_certificate.errors import CoincidentTargetSupportError
from compas_cgal.coverage_certificate.errors import UncoveredPocketError
from compas_cgal.coverage_certificate.exact_target import ExactCentreDomainTarget
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

F = Fraction
R = F(1)
TOOL = ToolRadius.build(1.0)
Cut = tuple[float, float, float]
XY = tuple[Fraction, Fraction]


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _point(x: Fraction | int, y: Fraction | int) -> ExactPoint2[WorldXY]:
    return ExactPoint2[WorldXY].build(F(x), F(y))


class _Pocket:
    """Sequential builder of a declared CCW cycle."""

    def __init__(self, start: ExactPoint2[WorldXY]) -> None:
        self.cursor = start
        self.primitives: list[DeclaredLine | DeclaredArc] = []

    def line(self, quarter: int, length: Fraction | int) -> _Pocket:
        line = DeclaredLine.build(self.cursor, PythagoreanDirection.build(F(0), quarter), F(length), len(self.primitives))
        self.primitives.append(line)
        self.cursor = line.end
        return self

    def arc(self, end: ExactPoint2[WorldXY], t: Fraction, quarter: int, sense: int) -> _Pocket:
        arc = DeclaredArc.build(self.cursor, end, PythagoreanDirection.build(t, quarter), 1 if sense == 1 else -1, len(self.primitives))
        self.primitives.append(arc)
        self.cursor = arc.end
        return self

    def close(self, name: str) -> DeclaredBoundary:
        return DeclaredBoundary.build(name, tuple(self.primitives), hashlib.sha256(name.encode()).hexdigest())


@cache
def square() -> tuple[DeclaredBoundary, ExactCentreDomain, ExactCentreDomainTarget]:
    boundary = _Pocket(_point(0, 0)).line(0, 20).line(1, 20).line(2, 20).line(3, 20).close("square")
    domain = ExactCentreDomain.build(boundary, R)
    return boundary, domain, exact_coverage_target(boundary, domain)


@cache
def convex_arc() -> tuple[DeclaredBoundary, ExactCentreDomain, ExactCentreDomainTarget]:
    """20 x 20 square whose bottom edge is a CCW arc about (10, 7.5), rho = 12.5, bulging to (10, -5)."""
    boundary = _Pocket(_point(0, 0)).arc(_point(20, 0), F(-1, 2), 0, 1).line(1, 20).line(2, 20).line(3, 20).close("convex_arc")
    domain = ExactCentreDomain.build(boundary, R)
    return boundary, domain, exact_coverage_target(boundary, domain)


@cache
def concave_arc() -> tuple[DeclaredBoundary, ExactCentreDomain, ExactCentreDomainTarget]:
    """20 x 20 square whose bottom edge is a CW arc about (10, -7.5), rho = 12.5, bulging into P to (10, 5)."""
    boundary = _Pocket(_point(0, 0)).arc(_point(20, 0), F(1, 2), 0, -1).line(1, 20).line(2, 20).line(3, 20).close("concave_arc")
    domain = ExactCentreDomain.build(boundary, R)
    return boundary, domain, exact_coverage_target(boundary, domain)


@cache
def l_shape() -> tuple[DeclaredBoundary, ExactCentreDomain, ExactCentreDomainTarget]:
    """L-shape with a reflex corner at (10, 10): a junction r-circle piece of dOmega*."""
    boundary = _Pocket(_point(0, 0)).line(0, 20).line(1, 10).line(2, 10).line(1, 10).line(2, 10).line(3, 20).close("l_shape")
    domain = ExactCentreDomain.build(boundary, R)
    return boundary, domain, exact_coverage_target(boundary, domain)


@cache
def tight_convex_arc() -> tuple[DeclaredBoundary, ExactCentreDomain, ExactCentreDomainTarget]:
    """20 x 20 square whose (20, 20) corner is a convex arc of rho = 1/2 < r (unreachable material, no piece)."""
    boundary = _Pocket(_point(0, 0)).line(0, 20).line(1, F(39, 2)).arc(_point(F(39, 2), 20), F(0), 1, 1).line(2, F(39, 2)).line(3, 20).close("tight_convex_arc")
    domain = ExactCentreDomain.build(boundary, R)
    return boundary, domain, exact_coverage_target(boundary, domain)


@pytest.fixture(autouse=True, scope="module")
def _release_native_domains() -> Iterator[None]:
    """Drop the cached native Omega* handles before interpreter shutdown (nanobind reports live instances at exit)."""
    yield
    for fixture in (square, convex_arc, concave_arc, l_shape, tight_convex_arc):
        fixture.cache_clear()


def _path_domain(points: Sequence[tuple[float, float]]) -> CentreDomain:
    return CentreDomain.build([Point2[WorldXY].build(x, y) for x, y in points])


# A path domain for runs without connectors: only carries sites, never the target.
UNUSED_PATH = _path_domain([(5.0, 5.0), (15.0, 5.0), (15.0, 15.0), (5.0, 15.0)])


def _sweeps(cuts: Sequence[Cut]) -> list[FullCircleSweep]:
    return [FullCircleSweep.build(Point2[WorldXY].build(x, y), GuideRadius.build(big_r)) for x, y, big_r in cuts]


def _grid(x0: float, y0: float, step: float, nx: int, ny: int, big_r: float = 1.0) -> list[Cut]:
    return [(x0 + step * i, y0 + step * j, big_r) for i in range(nx) for j in range(ny)]


def _both_modes(target: ExactCentreDomainTarget, cuts: Sequence[Cut], path: CentreDomain = UNUSED_PATH, connectors: Sequence[BoundaryConnectorSweep] = ()) -> CertificateReport:
    lazy = run_certificate(path, TOOL, _sweeps(cuts), connectors, exact_target=target, filter_mode=FilterMode.LAZY, audit=True)
    exact = run_certificate(path, TOOL, _sweeps(cuts), connectors, exact_target=target, filter_mode=FilterMode.EXACT_ONLY, audit=True)
    assert lazy.canonical_text() == exact.canonical_text()
    return lazy


# ---------------------------------------------------------------------------
# Independent oracles
# ---------------------------------------------------------------------------


def _in_target(domain: ExactCentreDomain, point: XY) -> bool:
    return domain.squared_distance_to(_point(*point)).compare(R * R) <= 0


def _segment_distance_sq(p: XY, a: XY, b: XY) -> Fraction:
    dx, dy = b[0] - a[0], b[1] - a[1]
    t = ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / (dx * dx + dy * dy)
    t = min(max(t, F(0)), F(1))
    qx, qy = a[0] + t * dx, a[1] + t * dy
    return (p[0] - qx) ** 2 + (p[1] - qy) ** 2


def _swept(point: XY, cuts: Sequence[Cut], capsules: Sequence[tuple[XY, XY]] = ()) -> bool:
    for x, y, big_r in cuts:
        d2 = (point[0] - F(x)) ** 2 + (point[1] - F(y)) ** 2
        outer, inner = (F(big_r) + R) ** 2, (max(F(big_r) - R, F(0))) ** 2
        if inner <= d2 <= outer:
            return True
    return any(_segment_distance_sq(point, a, b) <= R * R for a, b in capsules)


def _check_witnesses(domain: ExactCentreDomain, report: CertificateReport, cuts: Sequence[Cut], capsules: Sequence[tuple[XY, XY]] = ()) -> None:
    assert report.exposed
    for piece in report.exposed:
        point = (F(piece.uncovered[0]), F(piece.uncovered[1]))
        assert _in_target(domain, point), f"witness {piece} is not in T"
        assert not _swept(point, cuts, capsules), f"witness {piece} is swept"


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("fixture", [square, convex_arc, concave_arc, l_shape, tight_convex_arc], ids=lambda f: f.__name__)
def test_exact_target_membership_matches_the_centre_domain_distance_oracle(fixture: object) -> None:
    # T membership (exact_target_2.cpp) against dist(x, Omega*) <= r from the
    # native centre-domain module: random rationals over the pocket's box, and
    # a dense lattice around every kink and corner where the members meet.
    assert callable(fixture)
    _, domain, target = fixture()
    rng = random.Random(7)
    points: list[XY] = [(F(rng.randint(-120, 1320), 60), F(rng.randint(-420, 1320), 60)) for _ in range(600)]
    for cx, cy in ((0, 0), (20, 0), (20, 20), (0, 20), (10, 10), (10, -5), (10, 5), (F(39, 2), F(39, 2))):
        points += [(F(cx) + F(i, 8), F(cy) + F(j, 8)) for i in range(-12, 13) for j in range(-12, 13)]
    boundary_rows, piece_rows, disk_rows = target.wire()
    wire_points = [((str(x.numerator), str(x.denominator)), (str(y.numerator), str(y.denominator))) for x, y in points]
    lazy = native.contract_probe_exact_target_contains(1.0, boundary_rows, piece_rows, disk_rows, wire_points, False)
    exact = native.contract_probe_exact_target_contains(1.0, boundary_rows, piece_rows, disk_rows, wire_points, True)
    assert lazy == exact
    expected = [_in_target(domain, p) for p in points]
    mismatches = [(p, got) for p, got, want in zip(points, lazy, expected) if got != want]
    assert not mismatches, mismatches[:5]
    assert any(expected) and not all(expected)


def test_dense_disk_grid_certifies_the_square_exact_target() -> None:
    _, _, target = square()
    cuts = _grid(-1.0, -1.0, 2.5, 10, 10)
    report = _both_modes(target, cuts)
    assert not report.exposed
    assert report.undecided == 0
    certificate = CoverageCertificate.build(UNUSED_PATH, TOOL, _sweeps(cuts), exact_target=target)
    assert certificate.exposed_edge_count == 0


def test_a_missing_cut_is_proven_with_witnesses_in_t() -> None:
    _, domain, target = square()
    cuts = [c for c in _grid(-1.0, -1.0, 2.5, 10, 10) if (c[0], c[1]) != (9.0, 9.0)]
    report = _both_modes(target, cuts)
    _check_witnesses(domain, report, cuts)
    for piece in report.exposed:
        assert (F(piece.uncovered[0]) - 9) ** 2 + (F(piece.uncovered[1]) - 9) ** 2 < F(25, 4)
    with pytest.raises(UncoveredPocketError):
        CoverageCertificate.build(UNUSED_PATH, TOOL, _sweeps(cuts), exact_target=target)


def test_one_root_vertex_circles_cut_the_sweeps_and_a_dense_grid_passes() -> None:
    _, _, target = convex_arc()
    assert any(d.radicand != 0 for d in target.vertex_disks)
    cuts = _grid(-1.0, -6.0, 2.5, 10, 12)
    report = _both_modes(target, cuts)
    assert not report.exposed
    assert report.counters["vertex_circle_pairs"] > 0


def test_cusp_material_outside_t_is_never_a_witness() -> None:
    # Cuts only at tool centres inside Omega*: the cusps at the convex kinks are
    # unreachable material, outside T, and must not be reported.
    _, domain, target = convex_arc()
    rng = random.Random(20261021)
    cuts: list[Cut] = []
    for i in range(9):
        for j in range(10):
            x, y = 1.25 + 2.2 * i + rng.uniform(-0.2, 0.2), -2.75 + 2.2 * j + rng.uniform(-0.2, 0.2)
            if domain.locate(_point(F(x), F(y))).value != "exterior":
                cuts.append((x, y, rng.choice([0.6, 1.0, 1.3])))
    report = _both_modes(target, cuts)
    _check_witnesses(domain, report, cuts)
    assert report.counters["vertex_circle_pairs"] > 0


def test_a_hole_reaching_t_only_through_a_vertex_disk_is_found() -> None:
    # Found by a raster search: an annulus's inner hole that meets T only inside
    # D_r(w) at the convex kink (0, 0), where w is one-root. A certificate that
    # drops vertex disks from T misses it; the exact one must not.
    _, domain, target = convex_arc()
    corner = (0.3, 0.0)
    cuts: list[Cut] = []
    for i in range(10):
        for j in range(11):
            x, y = 0.25 + 0.7 * i, -4.5 + 0.7 * j
            if (x - corner[0]) ** 2 + (y - corner[1]) ** 2 > 2.6**2:
                cuts.append((x, y, 0.5))
    cuts.append((float.fromhex("-0x1.ecabe0ef293d8p-5"), float.fromhex("-0x1.bcde9f54d3658p-3"), float.fromhex("0x1.3780c43a4ca0ep+0")))
    report = _both_modes(target, cuts)
    _check_witnesses(domain, report, cuts)
    near_kink = [p for p in report.exposed if p.uncovered[0] < F(1, 2) and -F(1, 2) < p.uncovered[1] < F(1, 4)]
    assert near_kink, [piece.uncovered for piece in report.exposed]


def test_cuts_far_from_the_pocket_leave_the_interior_anchor_uncovered() -> None:
    _, domain, target = square()
    cuts: list[Cut] = [(100.0, 100.0, 1.0)]
    report = _both_modes(target, cuts)
    assert [piece.piece_kind for piece in report.exposed] == [ExposedPieceKind.TARGET_INTERIOR]
    _check_witnesses(domain, report, cuts)


def test_connector_flanks_inside_omega_star_expose_the_wall_strip() -> None:
    _, domain, target = square()
    path = _path_domain([(1.5, 1.5), (18.5, 1.5), (18.5, 18.5), (1.5, 18.5)])
    connectors = [
        BoundaryConnectorSweep.build(path, CentreBoundarySite.build(0, 0, F(0)), CentreBoundarySite.build(0, 2, F(0))),
        BoundaryConnectorSweep.build(path, CentreBoundarySite.build(0, 2, F(0)), CentreBoundarySite.build(0, 0, F(0))),
    ]
    cuts = _grid(4.0, 4.0, 2.5, 6, 6)
    report = _both_modes(target, cuts, path, connectors)
    corners = [(F(1.5), F(1.5)), (F(18.5), F(1.5)), (F(18.5), F(18.5)), (F(1.5), F(18.5))]
    capsules = [(corners[k], corners[(k + 1) % 4]) for k in range(4)]
    _check_witnesses(domain, report, cuts, capsules)
    assert ExposedPieceKind.CAPSULE_FLANK in {piece.piece_kind for piece in report.exposed}


def test_connector_flanks_cross_one_root_vertex_circles_exactly() -> None:
    # A polygonal path domain inside Omega* of the convex-arc pocket: its bottom
    # chain follows the offset arc at radius 11.4 < 11.5, so the outer flanks
    # pass 0.1 mm off the walls and cross the vertex circles at both kinks.
    _, domain, target = convex_arc()
    bottom = [(x, 7.5 - (129.96 - (x - 10.0) ** 2) ** 0.5) for x in (1.1, 3.0, 6.0, 10.0, 14.0, 17.0, 18.9)]
    outline = [*bottom, (18.9, 18.9), (1.1, 18.9)]
    path = _path_domain(outline)
    n = len(outline)
    connectors = [
        BoundaryConnectorSweep.build(path, CentreBoundarySite.build(0, 0, F(0)), CentreBoundarySite.build(0, n // 2, F(0))),
        BoundaryConnectorSweep.build(path, CentreBoundarySite.build(0, n // 2, F(0)), CentreBoundarySite.build(0, 0, F(0))),
    ]
    cuts = _grid(4.0, 4.0, 2.5, 6, 6)
    report = _both_modes(target, cuts, path, connectors)
    corners = [(F(x), F(y)) for x, y in outline]
    capsules = [(corners[k], corners[(k + 1) % n]) for k in range(n)]
    _check_witnesses(domain, report, cuts, capsules)
    assert report.counters["vertex_line_pairs"] > 0
    assert ExposedPieceKind.CAPSULE_FLANK in {piece.piece_kind for piece in report.exposed}


def test_connector_on_the_exact_boundary_of_omega_star_is_a_named_unmodelled_coincidence() -> None:
    _, _, target = square()
    path = _path_domain([(1.0, 1.0), (19.0, 1.0), (19.0, 19.0), (1.0, 19.0)])
    connectors = [BoundaryConnectorSweep.build(path, CentreBoundarySite.build(0, 0, F(0)), CentreBoundarySite.build(0, 2, F(0)))]
    with pytest.raises(CoincidentTargetSupportError):
        CoverageCertificate.build(path, TOOL, [], connectors, exact_target=target)
