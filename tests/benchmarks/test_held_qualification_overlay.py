"""Contract tests of the Held qualification overlay: it renders, and its caption says exactly what the certificate decided.

Fixtures use dyadic coordinates only, so every expected exact outcome is exact:

- pocket: the square [0, 12]^2, r = 1, centre domain Omega = [1, 11]^2, so the
  machinable target T = Omega (+) D_1 is the pocket with rounded corners;
- ``gap``: one cut R = 5 about (6, 6) sweeps the annulus 4 <= |x - c| <= 6; the
  disk |x - c| < 4 is in T and uncut, so the certificate rejects the path;
- ``covered``: pocket [0, 4]^2, Omega = [1, 3]^2, one cut R = 1/2 about (2, 2)
  and a closed boundary traversal of two connectors, (1, 1) -> (3, 3) and the
  seam (3, 3) -> (1, 1). Their capsules hold every point within 1 of the
  boundary of Omega, which is all of T, so the certificate passes; without the
  seam connector the far half of the band is uncut and it fails;
- gouge: a cut R = 5.25 about (6, 6) reaches 6.25 > 6, past every wall by 0.25.
"""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path

import pytest

from benchmarks.held_reference_geometry import ReferenceLine
from benchmarks.tools.held_qualification_overlay import CoverageQualification
from benchmarks.tools.held_qualification_overlay import CoverageStatus
from benchmarks.tools.held_qualification_overlay import OverlayCircle
from benchmarks.tools.held_qualification_overlay import OverlayInput
from benchmarks.tools.held_qualification_overlay import caption
from benchmarks.tools.held_qualification_overlay import choose_insets
from benchmarks.tools.held_qualification_overlay import gouge_audit
from benchmarks.tools.held_qualification_overlay import qualify_and_render
from benchmarks.tools.held_qualification_overlay import qualify_coverage
from benchmarks.tools.held_qualification_overlay import witness_clusters
from compas_cgal.coverage_certificate.sweeps import BoundaryConnectorSweep
from compas_cgal.coverage_certificate.target import CentreBoundarySite
from compas_cgal.coverage_certificate.target import CentreDomain
from compas_cgal.coverage_certificate.witness import ExposedPieceKind
from compas_cgal.coverage_certificate.witness import UncoveredWitness
from compas_cgal.units import ExactMillimetre
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

# Coarse raster for the fixtures, mm: they are a few mm across.
FIXTURE_PIXEL_MM = 0.05


def _square(low: float, high: float) -> tuple[Point2[WorldXY], ...]:
    return tuple(Point2[WorldXY].build(x, y) for x, y in ((low, low), (high, low), (high, high), (low, high)))


def _pocket(low: float, high: float) -> tuple[ReferenceLine, ...]:
    corners = _square(low, high)
    return tuple(ReferenceLine.build(a, b) for a, b in zip(corners, corners[1:] + corners[:1]))


def _corner(index: int) -> CentreBoundarySite:
    return CentreBoundarySite.build(loop=0, segment=index, parameter=Fraction(0))


def _overlay(
    pocket: tuple[float, float],
    omega: tuple[float, float],
    circles: list[tuple[float, float, float]],
    corner_hops: tuple[tuple[int, int], ...] = (),
) -> OverlayInput:
    """A fixture whose connectors run CCW along Omega between its corners (polyline and exact sweep agree)."""
    corners = _square(*omega)
    domain = CentreDomain.build(list(corners))
    polylines = tuple(tuple(corners[k % 4] for k in range(start, start + (end - start) % 4 + 1)) for start, end in corner_hops)
    return OverlayInput(
        name="fixture",
        title="Fixture",
        pocket=_pocket(*pocket),
        domain=domain,
        tool_radius=ToolRadius.build(1.0),
        circles=tuple(OverlayCircle(Point2[WorldXY].build(x, y), GuideRadius.build(r)) for x, y, r in circles),
        connectors=polylines,
        connector_sweeps=tuple(BoundaryConnectorSweep.build(domain, _corner(start), _corner(end)) for start, end in corner_hops),
    )


def test_an_annulus_with_an_uncut_core_is_not_covered_and_the_caption_counts_the_witnesses() -> None:
    overlay = _overlay((0.0, 12.0), (1.0, 11.0), [(6.0, 6.0, 5.0)])
    coverage = qualify_coverage(overlay)
    assert coverage.status is CoverageStatus.NOT_COVERED
    assert coverage.witnesses and coverage.undecided_edge_count == 0
    assert any(abs(x - 6.0) < 4.0 and abs(y - 6.0) < 4.0 for x, y in coverage.witnesses_mm)
    certified, diagnostic = caption(coverage, gouge_audit(overlay))
    assert certified.startswith("EXACT - coverage certificate: NOT COVERED")
    assert f"{len(coverage.witnesses)} exact witness" in certified
    assert "not yet complete" not in certified
    assert diagnostic.startswith("DIAGNOSTIC") and "never a decision" in diagnostic
    assert "raster" not in certified


def test_a_closed_boundary_traversal_is_certified_covered() -> None:
    overlay = _overlay((0.0, 4.0), (1.0, 3.0), [(2.0, 2.0, 0.5)], ((0, 2), (2, 0)))
    coverage = qualify_coverage(overlay)
    assert coverage.status is CoverageStatus.CERTIFIED_COVERED
    assert coverage.witnesses == ()
    assert coverage.regions == 3
    certified, _ = caption(coverage, gouge_audit(overlay))
    assert certified.startswith("EXACT - coverage certificate: CERTIFIED COVERED")
    assert "NOT COVERED" not in certified


def test_dropping_the_seam_connector_is_decided_not_covered() -> None:
    overlay = _overlay((0.0, 4.0), (1.0, 3.0), [(2.0, 2.0, 0.5)], ((0, 2),))
    coverage = qualify_coverage(overlay)
    assert coverage.status is CoverageStatus.NOT_COVERED
    # The open side is the band beyond the edges (3,3)-(1,3)-(1,1): x < 1 or y > 3.
    assert all(x < 1.0 or y > 3.0 for x, y in coverage.witnesses_mm)


def test_a_cut_past_the_walls_is_an_exact_gouge_in_the_caption() -> None:
    overlay = _overlay((0.0, 12.0), (1.0, 11.0), [(6.0, 6.0, 5.25)])
    report = gouge_audit(overlay)
    assert len(report.gouging_circles) == 1
    assert report.max_depth_mm("circle") == pytest.approx(0.25)
    certified, _ = caption(qualify_coverage(overlay), report)
    assert "gouge audit: 1 circles + 0 connectors cut outside the analytic pocket" in certified


def test_no_gouge_is_stated_as_an_exact_clearance() -> None:
    overlay = _overlay((0.0, 12.0), (1.0, 11.0), [(6.0, 6.0, 5.0)])
    certified, _ = caption(CoverageQualification(CoverageStatus.CERTIFIED_COVERED, (), 0, 1), gouge_audit(overlay))
    assert "no element cuts outside the analytic pocket" in certified


def test_witness_clusters_link_at_one_tool_radius_largest_first() -> None:
    clusters = witness_clusters([(0.0, 0.0), (10.0, 0.0), (0.5, 0.0), (1.4, 0.0), (10.5, 0.0), (20.0, 0.0), (0.0, 0.9)])
    assert [len(group) for group in clusters] == [4, 2, 1]


def _witness(x: Fraction, y: Fraction) -> UncoveredWitness:
    return UncoveredWitness(ExposedPieceKind.ANNULUS_INNER, 0, (ExactMillimetre(x), ExactMillimetre(y)))


def test_inset_labels_rank_equal_clusters_by_uncut_area_not_by_discovery_order() -> None:
    # Cut R = 5 about (6, 6): the window about the centre is wholly uncut core,
    # the window about (6, 2.5) straddles the annulus's inner edge |x - c| = 4.
    # Both clusters hold two witnesses; the edge cluster sorts first by position.
    overlay = _overlay((0.0, 12.0), (1.0, 11.0), [(6.0, 6.0, 5.0)])
    edge = (_witness(Fraction(6), Fraction(5, 2)), _witness(Fraction(25, 4), Fraction(5, 2)))
    core = (_witness(Fraction(6), Fraction(6)), _witness(Fraction(25, 4), Fraction(6)))
    coverage = CoverageQualification(CoverageStatus.NOT_COVERED, edge + core, 0, 1)
    insets = choose_insets(overlay, coverage, ())
    assert [(inset.label, inset.centre) for inset in insets] == [("A", (6.125, 6.0)), ("B", (6.125, 2.5))]


def test_the_fixture_renders_with_witness_insets_and_a_gouge_inset(tmp_path: Path) -> None:
    overlay = _overlay((0.0, 12.0), (1.0, 11.0), [(6.0, 6.0, 5.25)])
    output = tmp_path / "held_fixture_qualification.png"
    coverage, report, insets = qualify_and_render(overlay, output, main_pixel_mm=FIXTURE_PIXEL_MM)
    assert output.stat().st_size > 0
    assert coverage.status is CoverageStatus.NOT_COVERED
    assert [inset.note for inset in insets] == ["uncut", "uncut", "gouge"]
    assert report.gouges
