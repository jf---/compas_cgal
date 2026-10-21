"""Slow-tier contract: forced EXACT_ONLY is byte-identical to LAZY on the 4367 Monstera circles.

Not collected by pytest (``testpaths = ["tests"]``); run deliberately with
``pixi run pytest benchmarks/acceptance/test_coverage_certificate_monstera_exact_slow.py -q``.
EXACT_ONLY sends every one of the 454k pair decisions to the exact back end:
12-90 s depending on load, too slow for the default suite, where the LAZY counts
are checked by ``test_monstera_circles_match_the_spike_counts``.
"""

from __future__ import annotations

from pathlib import Path

from compas_cgal.coverage_certificate.circle_core import FilterMode
from compas_cgal.coverage_certificate.circle_core import run_circle_core
from compas_cgal.coverage_certificate.sweeps import FullCircleSweep
from compas_cgal.coverage_certificate.target import CentreDomain
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

FIXTURES = Path(__file__).parents[2] / "tests" / "fixtures" / "coverage_certificate"


def test_forced_exact_mode_is_byte_identical_on_monstera_circles() -> None:
    lines = (FIXTURES / "monstera_circles.txt").read_text().split("\n")
    count, radius = lines[0].split()
    cuts = []
    for line in lines[1 : 1 + int(count)]:
        x, y, big_r = (float.fromhex(token) for token in line.split())
        cuts.append(FullCircleSweep.build(Point2[WorldXY].build(x, y), GuideRadius.build(big_r)))
    rows = (FIXTURES / "monstera_centre_domain.txt").read_text().split("\n")
    domain = CentreDomain.build([Point2[WorldXY].build(float.fromhex(x), float.fromhex(y)) for x, y in (row.split() for row in rows[1 : 1 + int(rows[0])])])
    tool = ToolRadius.build(float.fromhex(radius))
    lazy = run_circle_core(domain, tool, cuts, filter_mode=FilterMode.LAZY)
    exact = run_circle_core(domain, tool, cuts, filter_mode=FilterMode.EXACT_ONLY)
    assert lazy.canonical_text() == exact.canonical_text()
    assert lazy.counters["sub_arcs"] == 1_181_504
    assert len(lazy.exposed) == 4194
