"""Slow-tier contracts of the production certificate on the generated Monstera path.

Not collected by pytest (``testpaths = ["tests"]``); run deliberately with the
editable-skip pixi invocation and
``pytest benchmarks/acceptance/test_coverage_certificate_monstera_path_slow.py -n auto --dist=loadscope -q``
(``loadscope`` keeps the module-scoped paths on one worker). Generating the
path takes about 14 minutes (824.6 s measured 2026-10-17); the unrepaired and
the repaired path are each built once per module.

The input is the generator's own closed path: ``generate`` (closing seam
connector included) and ``path_sweeps``, exactly what
``tools.coverage_certificate_acceptance`` records in
``benchmarks/results/coverage_certificate_acceptance.md``.
"""

from __future__ import annotations

import math

import pytest

from benchmarks.held_coverage_repair import path_sweeps
from benchmarks.held_coverage_repair import refine_boundary_ordered_path_certified
from benchmarks.held_figure5_path import build_held_reference_path
from benchmarks.held_figure5_raw_guide import build_held_reference_raw_guide
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import figure7_inward_offset
from benchmarks.held_toolpath_benchmark import generate
from compas_cgal.coverage_certificate.certificate import CoverageCertificate
from compas_cgal.coverage_certificate.circle_core import FilterMode
from compas_cgal.coverage_certificate.errors import UncoveredPocketError
from compas_cgal.coverage_certificate.native_certificate import run_certificate
from compas_cgal.coverage_certificate.sweeps import BoundaryConnectorSweep
from compas_cgal.coverage_certificate.sweeps import FullCircleSweep
from compas_cgal.coverage_certificate.target import CentreDomain
from compas_cgal.motion import EngagementCap
from compas_cgal.units import ToolRadius
from tools.coverage_certificate_clusters import cluster_witnesses

MONSTERA = "figure8_monstera_boundary_qualified_corners"
# Raster-located uncut slivers of the unrepaired path, mm (reporting values,
# .superpowers/sdd/2026-10-13-product-pipeline/held-monstera-uncut-report.md).
SLIVERS_MM = {"A": (45.45, 61.08), "B": (33.97, 71.41)}
# Reporting-only proximity for matching an exact witness to a sliver, mm: each
# sliver is ~0.2 mm across and they are >15 mm apart, so 0.5 mm is unambiguous.
SLIVER_MATCH_MM = 0.5

CertificateInput = tuple[CentreDomain, ToolRadius, tuple[FullCircleSweep, ...], tuple[BoundaryConnectorSweep, ...]]


@pytest.fixture(scope="module")
def unrepaired() -> CertificateInput:
    path = generate(MONSTERA)
    domain = CentreDomain.build(list(path.boundary))
    cuts, connectors = path_sweeps(domain, path.boundary, path.circles)
    return domain, path.case.tool_radius, cuts, connectors


@pytest.fixture(scope="module")
def repaired() -> CertificateInput:
    case = load_held_reference_case("figure8_monstera")
    components = figure7_inward_offset(case).components
    assert len(components) == 1, "the Monstera centre domain is one loop without holes"
    source = build_held_reference_path(build_held_reference_raw_guide(case), components, qualify_contacts=True)
    cap = EngagementCap.build(math.radians(float(case.tea_cap)))
    certified = refine_boundary_ordered_path_certified(source, components[0], case.tool_radius, cap, corner_approaches=True)
    domain = CentreDomain.build(list(components[0]))
    cuts, connectors = path_sweeps(domain, components[0], certified.circles)
    return domain, case.tool_radius, cuts, connectors


def test_full_path_forced_exact_mode_is_byte_identical_to_lazy(unrepaired: CertificateInput) -> None:
    domain, tool, cuts, connectors = unrepaired
    lazy = run_certificate(domain, tool, cuts, connectors, filter_mode=FilterMode.LAZY)
    exact = run_certificate(domain, tool, cuts, connectors, filter_mode=FilterMode.EXACT_ONLY)
    assert connectors, "the full path has connectors"
    assert lazy.exposed, "the unrepaired path is uncovered, so the identity covers exact witnesses"
    assert lazy.canonical_text() == exact.canonical_text()


def test_unrepaired_path_witnesses_lie_at_the_two_slivers_and_nowhere_else(unrepaired: CertificateInput) -> None:
    with pytest.raises(UncoveredPocketError) as caught:
        CoverageCertificate.build(*unrepaired)
    witnesses = caught.value.witnesses
    report = "\n".join(f"{w.piece_kind.name} {w.uncovered_point_mm}" for w in witnesses)
    for w in witnesses:
        assert any(math.dist(w.uncovered_point_mm, s) < SLIVER_MATCH_MM for s in SLIVERS_MM.values()), f"witness off both slivers:\n{report}"
    for name, sliver in SLIVERS_MM.items():
        assert any(math.dist(w.uncovered_point_mm, sliver) < SLIVER_MATCH_MM for w in witnesses), f"no witness at sliver {name}:\n{report}"
    assert len(cluster_witnesses(witnesses)) == len(SLIVERS_MM), report


def test_repaired_path_is_certified_covered(repaired: CertificateInput) -> None:
    certificate = CoverageCertificate.build(*repaired)
    assert certificate.exposed_edge_count == 0
