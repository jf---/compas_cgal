"""Contract tests for the circle-circle core of the coverage certificate (phase 1).

The core runs behind a certified interval front with Circular_kernel_2 as the
exact back end. These tests hold the front to its contract:

- forced EXACT_ONLY mode gives byte-identical decisions and witnesses to LAZY
  on the spike's degenerate integer grid (the 4367 Monstera circles: slow tier);
- generic doubles: every witness re-checks exactly with ``Fraction`` here,
  independently of the native code;
- a near-tangent pair whose two certified branch boxes overlap is labelled by
  the exact rational branch matcher (the case the spike could only abort on);
- the audit oracle (every covered edge re-witnessed by pure exact membership)
  agrees with the counters;
- VERTEX_DISK target members (one-root centres) decide witness membership.
"""

from __future__ import annotations

import random
from collections.abc import Sequence
from fractions import Fraction
from pathlib import Path

import pytest

from compas_cgal.coverage_certificate.certificate import CoverageCertificate
from compas_cgal.coverage_certificate.circle_core import CircleCoreReport
from compas_cgal.coverage_certificate.circle_core import FilterMode
from compas_cgal.coverage_certificate.circle_core import run_circle_core
from compas_cgal.coverage_certificate.errors import DegenerateConnectorError
from compas_cgal.coverage_certificate.errors import InvalidCentreBoundarySiteError
from compas_cgal.coverage_certificate.errors import InvalidCentreDomainError
from compas_cgal.coverage_certificate.errors import InvalidVertexDiskError
from compas_cgal.coverage_certificate.errors import UncoveredPocketError
from compas_cgal.coverage_certificate.errors import UnmodelledCertificatePieceError
from compas_cgal.coverage_certificate.sweeps import BoundaryConnectorSweep
from compas_cgal.coverage_certificate.sweeps import FullCircleSweep
from compas_cgal.coverage_certificate.target import CentreBoundarySite
from compas_cgal.coverage_certificate.target import CentreDomain
from compas_cgal.coverage_certificate.target import OneRootMillimetre
from compas_cgal.coverage_certificate.target import VertexDisk
from compas_cgal.coverage_certificate.witness import ExposedPieceKind
from compas_cgal.coverage_certificate.witness import UncoveredWitness
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

FIXTURES = Path(__file__).parent / "fixtures" / "coverage_certificate"

Cut = tuple[float, float, float]  # centre x, centre y, path radius R


def _read_circles(name: str) -> tuple[float, list[Cut]]:
    lines = (FIXTURES / name).read_text().split("\n")
    count, radius = lines[0].split()
    cuts = []
    for line in lines[1 : 1 + int(count)]:
        x, y, big_r = (float.fromhex(token) for token in line.split())
        cuts.append((x, y, big_r))
    return float.fromhex(radius), cuts


def _domain(points: Sequence[tuple[float, float]]) -> CentreDomain:
    return CentreDomain.build([Point2[WorldXY].build(x, y) for x, y in points])


def _sweeps(cuts: Sequence[Cut]) -> list[FullCircleSweep]:
    return [FullCircleSweep.build(Point2[WorldXY].build(x, y), GuideRadius.build(big_r)) for x, y, big_r in cuts]


def _run(cuts: Sequence[Cut], tool_radius: float, domain: CentreDomain, mode: FilterMode, audit: bool = False) -> CircleCoreReport:
    return run_circle_core(domain, ToolRadius.build(tool_radius), _sweeps(cuts), filter_mode=mode, audit=audit)


def _assert_witnesses_exact(report: CircleCoreReport, cuts: Sequence[Cut], tool_radius: float) -> None:
    """Independent Fraction oracle for every exposed edge's witnesses."""
    r = Fraction(tool_radius)
    regions = [(Fraction(x), Fraction(y), Fraction(big_r) + r, max(Fraction(big_r) - r, Fraction(0))) for x, y, big_r in cuts]

    def in_closed(p: tuple[Fraction, Fraction], region: tuple[Fraction, Fraction, Fraction, Fraction]) -> bool:
        cx, cy, outer, inner = region
        d2 = (p[0] - cx) ** 2 + (p[1] - cy) ** 2
        return inner * inner <= d2 <= outer * outer

    for edge in report.exposed:
        cx, cy, outer, inner = regions[edge.region_index]
        rho = inner if edge.piece_kind is ExposedPieceKind.ANNULUS_INNER else outer
        w = edge.on_edge
        assert (w[0] - cx) ** 2 + (w[1] - cy) ** 2 == rho * rho, f"{edge} witness is not on its circle"
        assert not any(in_closed(w, region) for i, region in enumerate(regions) if i != edge.region_index), f"{edge} witness is covered"
        assert not any(in_closed(edge.uncovered, region) for region in regions), f"{edge} uncovered point is swept"


GRID_DOMAIN = [(-1.0, -1.0), (11.0, -1.0), (11.0, 11.0), (-1.0, 11.0)]


def test_forced_exact_mode_is_byte_identical_on_the_degenerate_grid() -> None:
    tool_radius, cuts = _read_circles("degenerate_grid.txt")
    domain = _domain(GRID_DOMAIN)
    lazy = _run(cuts, tool_radius, domain, FilterMode.LAZY, audit=True)
    exact = _run(cuts, tool_radius, domain, FilterMode.EXACT_ONLY, audit=True)
    assert lazy.canonical_text() == exact.canonical_text()
    # The fixture really is degenerate, and the lazy front really fell back.
    assert lazy.counters["tangencies"] > 0
    assert lazy.counters["coincident_pairs"] >= 2
    assert lazy.counters["start_retries"] > 0
    assert lazy.counters["exact_compare_x"] > 0
    assert lazy.counters["exact_half"] > 0
    assert lazy.exposed, "the grid's outer envelope is exposed"
    _assert_witnesses_exact(lazy, cuts, tool_radius)


def test_monstera_circles_match_the_spike_counts() -> None:
    # LAZY only: the forced EXACT_ONLY byte-identity on these 4367 circles takes
    # 12-90 s and lives in the slow tier
    # (benchmarks/acceptance/test_coverage_certificate_monstera_exact_slow.py).
    tool_radius, cuts = _read_circles("monstera_circles.txt")
    domain = _monstera_domain()
    lazy = _run(cuts, tool_radius, domain, FilterMode.LAZY)
    # The spike's independent count: 1,181,504 sub-arcs, 4194 exposed.
    assert lazy.counters["sub_arcs"] == 1_181_504
    assert len(lazy.exposed) == 4194
    assert lazy.counters["pairs_exact"] <= 1, "generic doubles must not sit on a filter boundary"


def _monstera_domain() -> CentreDomain:
    lines = (FIXTURES / "monstera_centre_domain.txt").read_text().split("\n")
    return _domain([(float.fromhex(x), float.fromhex(y)) for x, y in (line.split() for line in lines[1 : 1 + int(lines[0])])])


def test_generic_double_witnesses_are_exactly_uncovered() -> None:
    rng = random.Random(20261017)
    tool_radius = rng.uniform(0.5, 1.0)
    cuts = [(rng.uniform(0.0, 20.0), rng.uniform(0.0, 20.0), rng.uniform(0.2, 3.0)) for _ in range(80)]
    domain = _domain([(0.0, 0.0), (20.0, 0.0), (20.0, 20.0), (0.0, 20.0)])
    lazy = _run(cuts, tool_radius, domain, FilterMode.LAZY, audit=True)
    exact = _run(cuts, tool_radius, domain, FilterMode.EXACT_ONLY)
    assert lazy.canonical_text() == exact.canonical_text()
    assert lazy.exposed
    assert lazy.counters["pairs_exact"] == 0 and lazy.counters["exact_compare_x"] == 0, "generic input is decided by the filter"
    _assert_witnesses_exact(lazy, cuts, tool_radius)


# Two disks of radius 1 (R = r = 1/2) whose centres are 2 - 2^-26 apart, at height
# 2^40: the crossings are 2 s = 1.2e-4 apart in y, below one ulp (2.4e-4) of
# 2^40, so the certified + and - branch boxes overlap although the interval
# discriminant certifies two transversal crossings.
NEAR_TANGENT = [(0.0, 2.0**40, 0.5), (2.0 - 2.0**-26, 2.0**40, 0.5)]
NEAR_TANGENT_DOMAIN = [(-2.0, 2.0**40 - 2.0), (4.0, 2.0**40 - 2.0), (1.0, 2.0**40 + 3.0)]


def test_near_tangent_pair_hits_the_exact_branch_matcher() -> None:
    domain = _domain(NEAR_TANGENT_DOMAIN)
    lazy = _run(NEAR_TANGENT, 0.5, domain, FilterMode.LAZY, audit=True)
    exact = _run(NEAR_TANGENT, 0.5, domain, FilterMode.EXACT_ONLY, audit=True)
    assert lazy.counters["pairs_filtered_transversal"] == 1, "the filter certified the pair transversal"
    assert lazy.counters["exact_branch_matches"] >= 1, "the overlapping branch boxes went to the exact matcher"
    assert lazy.canonical_text() == exact.canonical_text()
    # Each disk boundary is cut at the two crossings; the arc outside the other
    # disk is exposed, the thin lens arc is covered.
    assert sorted((edge.region_index, edge.piece_kind) for edge in lazy.exposed) == [(0, ExposedPieceKind.DISK), (1, ExposedPieceKind.DISK)]
    _assert_witnesses_exact(lazy, NEAR_TANGENT, 0.5)


# The three-stack gap of tests/test_held_coverage_certificate.py, but with a tiny
# centre domain far from the gap, so that only a VERTEX_DISK member puts the gap
# in T. Its centre (4, 33/14 + sqrt(2)/100) is a one-root point.
GAP_CUTS = [(cx, cy, big_r) for cx, cy in ((0.0, 0.0), (8.0, 0.0), (4.0, 7.0)) for big_r in (3.25, 1.25, 0.25)]
FAR_DOMAIN = [(-30.0, -30.0), (-29.0, -30.0), (-29.5, -29.0)]
GAP_VERTEX = VertexDisk.build(OneRootMillimetre(Fraction(4), Fraction(0)), OneRootMillimetre(Fraction(33, 14), Fraction(1, 100)), Fraction(2))


def _sign_one_root(a: Fraction, b: Fraction, d: Fraction) -> int:
    """sign(a + b sqrt(d)) for rational a, b, d >= 0, by squaring."""
    sa, sb = (a > 0) - (a < 0), (b > 0) - (b < 0)
    if sb == 0 or d == 0:
        return sa
    if sa == 0 or sa == sb:
        return sb
    return sa if a * a > b * b * d else (0 if a * a == b * b * d else sb)


def _near_gap(witness: UncoveredWitness) -> bool:
    x, y = witness.uncovered_point
    return (x - 4) ** 2 + (y - Fraction(33, 14)) ** 2 < 1


def test_vertex_disk_member_alone_puts_the_gap_in_the_target() -> None:
    domain = _domain(FAR_DOMAIN)
    sweeps = _sweeps(GAP_CUTS)
    # Without the vertex disk T is the far triangle alone: no cut reaches it, so
    # it is proven uncovered there, and nothing near the gap is target material.
    with pytest.raises(UncoveredPocketError) as bare:
        CoverageCertificate.build(domain, ToolRadius.build(1.0), sweeps)
    assert bare.value.witnesses
    assert not any(_near_gap(witness) for witness in bare.value.witnesses)
    with pytest.raises(UncoveredPocketError) as caught:
        CoverageCertificate.build(domain, ToolRadius.build(1.0), sweeps, vertex_disks=[GAP_VERTEX])
    gap = [witness for witness in caught.value.witnesses if _near_gap(witness)]
    assert gap, "the vertex disk puts the gap in T"
    for witness in gap:
        x, y = witness.uncovered_point
        ex, ey = x - Fraction(4), y - Fraction(33, 14)
        # |p - w|^2 - 1 = (ex^2 + ey^2 + d/10^4 - 1) - (2 ey / 100) sqrt(d), d = 2.
        assert _sign_one_root(ex * ex + ey * ey + Fraction(2, 10_000) - 1, -2 * ey / 100, Fraction(2)) <= 0, f"{witness} is outside the vertex disk"


def test_covered_target_with_a_vertex_disk_is_unmodelled_not_passed() -> None:
    # One disk cut (radius 2) covers the far triangle's T and a VERTEX_DISK inside
    # it. No uncovered point exists, but the vertex disk's one-root boundary is not
    # cut into the arrangement yet, so the certificate refuses to pass.
    domain = _domain(FAR_DOMAIN)
    cover = _sweeps([(-29.5, -29.625, 1.0)])
    vertex = VertexDisk.build(OneRootMillimetre(Fraction(-59, 2), Fraction(0)), OneRootMillimetre(Fraction(-237, 8), Fraction(1, 100)), Fraction(2))
    assert CoverageCertificate.build(domain, ToolRadius.build(1.0), cover).exposed_edge_count == 0
    with pytest.raises(UnmodelledCertificatePieceError):
        CoverageCertificate.build(domain, ToolRadius.build(1.0), cover, vertex_disks=[vertex])


def test_isolated_cut_is_exposed_on_its_whole_circle() -> None:
    domain = _domain([(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)])
    report = _run([(0.25, 0.25, 2.0)], 1.0, domain, FilterMode.LAZY)
    kinds = sorted(edge.piece_kind.name for edge in report.exposed)
    assert kinds == ["ANNULUS_INNER", "ANNULUS_OUTER"]
    assert all(edge.whole_circle for edge in report.exposed)
    inner = next(edge for edge in report.exposed if edge.piece_kind is ExposedPieceKind.ANNULUS_INNER)
    assert inner.uncovered_in_target, "the unswept hole of the annulus contains target material"
    with pytest.raises(UncoveredPocketError):
        CoverageCertificate.build(domain, ToolRadius.build(1.0), _sweeps([(0.25, 0.25, 2.0)]))


def test_opposite_side_coincident_support_covers_the_shared_circle() -> None:
    # Annulus R = 3 (radii 2..4) and disk R = 1 <= r (radius 2) share the circle of
    # radius 2 with opposite sides: each covers the other's outer side there (P6),
    # so only the annulus's outer circle is exposed.
    domain = _domain([(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)])
    report = _run([(0.0, 0.0, 3.0), (0.0, 0.0, 1.0)], 1.0, domain, FilterMode.LAZY)
    assert report.counters["coincident_pairs"] == 1
    assert report.counters["covered_by_coincident_support"] == 2
    assert [(edge.region_index, edge.piece_kind) for edge in report.exposed] == [(0, ExposedPieceKind.ANNULUS_OUTER)]


@pytest.mark.parametrize(
    "outer",
    [
        [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0)],  # clockwise
        [(0.0, 0.0), (1.0, 1.0), (1.0, 0.0), (0.0, 1.0)],  # self-intersecting
        [(0.0, 0.0), (1.0, 0.0), (1.0, 0.0), (0.0, 1.0)],  # zero-length edge
        [(0.0, 0.0), (1.0, 0.0)],  # too few vertices
    ],
)
def test_invalid_centre_domain_is_rejected(outer: list[tuple[float, float]]) -> None:
    with pytest.raises(InvalidCentreDomainError):
        _domain(outer)


def test_hole_must_be_clockwise_and_strictly_inside() -> None:
    outer = [Point2[WorldXY].build(x, y) for x, y in [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]]
    hole_cw = [Point2[WorldXY].build(x, y) for x, y in [(1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0)]]
    assert len(CentreDomain.build(outer, [hole_cw]).loops) == 2
    with pytest.raises(InvalidCentreDomainError):
        CentreDomain.build(outer, [list(reversed(hole_cw))])
    outside = [Point2[WorldXY].build(x, y) for x, y in [(5.0, 1.0), (5.0, 2.0), (6.0, 2.0), (6.0, 1.0)]]
    with pytest.raises(InvalidCentreDomainError):
        CentreDomain.build(outer, [outside])


def test_connector_sites_are_validated() -> None:
    domain = _domain([(0.0, 0.0), (8.0, 0.0), (4.0, 7.0)])
    site = CentreBoundarySite.build(loop=0, segment=0, parameter=Fraction(1, 4))
    with pytest.raises(DegenerateConnectorError):
        BoundaryConnectorSweep.build(domain, site, site)
    with pytest.raises(InvalidCentreBoundarySiteError):
        BoundaryConnectorSweep.build(domain, site, CentreBoundarySite.build(loop=0, segment=3, parameter=Fraction(0)))
    with pytest.raises(InvalidCentreBoundarySiteError):
        CentreBoundarySite.build(loop=0, segment=0, parameter=Fraction(1))


def test_vertex_disk_rejects_a_negative_radicand() -> None:
    with pytest.raises(InvalidVertexDiskError):
        VertexDisk.build(OneRootMillimetre(Fraction(0), Fraction(1)), OneRootMillimetre(Fraction(0), Fraction(0)), Fraction(-1))
