"""Slow-tier acceptance: the published Monstera path is NOT covered, and the certificate says where.

Not collected by pytest (benchmarks/acceptance is outside `testpaths`); this is
the repo's slow tier, alongside test_held_toolpath_performance_slow.py. Run with
`pixi run pytest benchmarks/acceptance/test_held_coverage_certificate_monstera_slow.py -q`.
Path generation alone takes about 14 minutes (824.6 s measured 2026-10-17).

Expected outcome, from the diagnostic raster in
.superpowers/sdd/2026-10-13-product-pipeline/held-monstera-uncut-report.md
(reporting only, never a certificate): two interior uncut slivers,
A = 0.036 mm^2 at (45.45, 61.09) between circles 4018 / 1008 / 2706 and
B = 0.029 mm^2 at (33.97, 71.41) between circles 159 / 4039 / 570. The exact
certificate must reject the path and name an exposed edge at each.
"""

from __future__ import annotations

import math

import pytest

from benchmarks.held_boundary_order_spacing import refine_boundary_ordered_path
from benchmarks.held_figure5_boundary_path import Figure5BoundaryConnectors
from benchmarks.held_figure5_boundary_path import Figure5SeamClosure
from benchmarks.held_figure5_path import build_held_reference_path
from benchmarks.held_figure5_raw_guide import build_held_reference_raw_guide
from benchmarks.held_material_domain import PartPolygon
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import figure7_inward_offset
from compas_cgal.coverage_certificate.certificate import CoverageCertificate
from compas_cgal.coverage_certificate.errors import UncoveredPocketError
from compas_cgal.coverage_certificate.sweeps import BoundaryConnectorSweep
from compas_cgal.coverage_certificate.sweeps import FullCircleSweep
from compas_cgal.coverage_certificate.target import CentreBoundarySite
from compas_cgal.coverage_certificate.target import CentreDomain
from compas_cgal.motion import EngagementCap

# Raster-located uncut components, mm (reporting values from the diagnostic).
SLIVER_A_MM = (45.45, 61.09)
SLIVER_B_MM = (33.97, 71.41)
# Reporting-only proximity for matching an exact witness to a raster centroid,
# mm: each sliver is ~0.2 mm across at 0.02 mm/px, so 0.5 mm is unambiguous.
SLIVER_MATCH_MM = 0.5


def test_published_monstera_path_is_rejected_with_witnesses_at_both_slivers() -> None:
    case = load_held_reference_case("figure8_monstera")
    guide = build_held_reference_raw_guide(case)
    components = figure7_inward_offset(case).components
    assert len(components) == 1, "the Monstera centre domain is one loop without holes"
    result = build_held_reference_path(guide, components, qualify_contacts=True)
    cap = EngagementCap.build(math.radians(float(case.tea_cap)))
    _, refined = refine_boundary_ordered_path(result, components[0], case.tool_radius, cap, part=PartPolygon.from_case(case), corner_approaches=True)
    circles = refined.circles

    domain = CentreDomain.build(list(components[0]))
    cuts = [FullCircleSweep.build(circle.center, circle.radius) for circle in circles]
    # Boundary-order connectors run CCW along the centre loop between consecutive
    # contact sites, then close CCW from the last site back to the first unless
    # the seam contacts coincide; a fan (same site) is stationary and sweeps
    # nothing new.
    seam = Figure5BoundaryConnectors.build(components[0], circles).seam
    successors = circles[1:] + (circles[:1] if seam is Figure5SeamClosure.CONNECTOR else ())
    connectors = [
        BoundaryConnectorSweep.build(
            domain,
            CentreBoundarySite.build(loop=0, segment=int(a.boundary_site.segment_id), parameter=a.boundary_site.parameter),
            CentreBoundarySite.build(loop=0, segment=int(b.boundary_site.segment_id), parameter=b.boundary_site.parameter),
        )
        for a, b in zip(circles[: len(successors)], successors, strict=True)
        if a.boundary_site != b.boundary_site
    ]

    with pytest.raises(UncoveredPocketError) as caught:
        CoverageCertificate.build(domain, case.tool_radius, cuts, connectors)

    witnesses = caught.value.witnesses
    report = "\n".join(f"{w.piece_kind.name} {w.uncovered_point_mm}" for w in witnesses)
    for name, (sx, sy) in (("A", SLIVER_A_MM), ("B", SLIVER_B_MM)):
        assert any(math.dist(w.uncovered_point_mm, (sx, sy)) < SLIVER_MATCH_MM for w in witnesses), f"no witness at sliver {name}:\n{report}"
