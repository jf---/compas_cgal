"""The published Figure 5 workload generated on Ω*: zero exact gouges, exact arc connectors, closed seam."""

from __future__ import annotations

from functools import cache

import pytest

from benchmarks.held_figure5_boundary_path import Figure5SeamClosure
from benchmarks.held_omega_toolpath import OmegaGeneratedToolpath
from benchmarks.held_omega_toolpath import PathDomain
from benchmarks.held_omega_toolpath import generate_omega
from benchmarks.held_omega_toolpath import omega_metrics
from benchmarks.held_toolpath_benchmark import generate_on

pytestmark = pytest.mark.xdist_group("held_omega_figure5")


@cache
def _figure5() -> OmegaGeneratedToolpath:
    return generate_omega("figure5_refined")


def test_figure5_on_omega_has_no_exact_gouge() -> None:
    report = _figure5().audit()
    assert not report.gouging_circles
    assert not report.gouging_connectors
    assert report.circles == report.connectors == len(_figure5().circles)


def test_figure5_on_omega_follows_arcs_and_closes_its_seam() -> None:
    path = _figure5()
    metrics = omega_metrics(path)
    assert metrics.arc_connector_pieces > 0
    assert path.connectors.seam is Figure5SeamClosure.CONNECTOR
    for previous, following in zip(path.connectors.connectors, path.connectors.connectors[1:]):
        assert previous.end == following.start
    assert path.connectors.connectors[-1].end == path.connectors.connectors[0].start


def test_figure5_on_omega_reaches_every_guide_run() -> None:
    path = _figure5()
    assert path.initial.reached_run_ids == frozenset(run for runs in path.initial.circle_sources for run in runs)


def test_the_benchmark_switch_keeps_the_polygon_path_by_default() -> None:
    assert type(generate_on("figure5_refined", PathDomain.POLYGON)).__name__ == "GeneratedToolpath"
