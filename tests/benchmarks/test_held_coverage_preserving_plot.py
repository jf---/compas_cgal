"""Coverage-preserving comparison plot accepts open and closed (seam) paths."""

import json
from fractions import Fraction
from pathlib import Path

import pytest

from benchmarks.held_figure5_boundary_path import Figure5BoundaryPath
from benchmarks.held_figure5_boundary_path import Figure5BoundaryTransition
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.tools.held_coverage_preserving_plot import InvalidCoveragePreservingPlotError
from benchmarks.tools.held_coverage_preserving_plot import render_coverage_preserving_comparison
from benchmarks.tools.held_sequence_plot import CoverageVerdict
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import WorldXY


def _circles(xs: tuple[float, ...]) -> tuple[Figure5CounterclockwiseCircle, ...]:
    circles = []
    for index, x in enumerate(xs):
        site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), Fraction(x / 10), 4)
        circles.append(
            Figure5CounterclockwiseCircle(GuideRunId(0), GuideRunStationOrdinal(index), site, site, Point2[WorldXY].build(x, 1), Point2[WorldXY].build(x, 0), GuideRadius.build(1))
        )
    return tuple(circles)


def _connect(a: Figure5CounterclockwiseCircle, b: Figure5CounterclockwiseCircle) -> Figure5BoundaryTransition:
    return Figure5BoundaryTransition((a.contact_point, b.contact_point), Millimetre(abs(float(b.contact_point.x) - float(a.contact_point.x))), False)


def _path(xs: tuple[float, ...], *, closed: bool) -> Figure5BoundaryPath:
    circles = _circles(xs)
    transitions = tuple(_connect(a, b) for a, b in zip(circles, circles[1:]))
    if closed:
        transitions += (_connect(circles[-1], circles[0]),)
    return Figure5BoundaryPath(circles, transitions, frozenset({GuideRunId(0)}))


@pytest.mark.parametrize("closed", (True, False))
def test_renders_open_and_closed_paths(closed: bool, tmp_path: Path) -> None:
    before, after = _path((2, 3, 4, 5), closed=closed), _path((2, 4, 5), closed=closed)
    output = tmp_path / "plot.png"
    render_coverage_preserving_comparison(load_held_reference_case("figure5"), before, after, (Radian(0.5), Radian(0.6)), output)
    report = json.loads(output.with_suffix(".json").read_text())
    assert output.stat().st_size > 0
    assert report["circles"] == 3
    assert report["connectors"] == (3 if closed else 2)


def test_rejects_connector_count_neither_open_nor_closed(tmp_path: Path) -> None:
    before = _path((2, 3, 4), closed=True)
    after = _path((2, 4, 5), closed=True)
    broken = Figure5BoundaryPath(after.circles, after.transitions[:1], after.reached_run_ids)
    with pytest.raises(InvalidCoveragePreservingPlotError, match="connectors"):
        render_coverage_preserving_comparison(load_held_reference_case("figure5"), before, broken, (Radian(0.5), Radian(0.6)), tmp_path / "plot.png")


def test_rejects_bound_count_other_than_successor_count(tmp_path: Path) -> None:
    path = _path((2, 4, 5), closed=True)
    with pytest.raises(InvalidCoveragePreservingPlotError, match="bound"):
        render_coverage_preserving_comparison(load_held_reference_case("figure5"), path, path, (Radian(0.5), Radian(0.6), Radian(0.7)), tmp_path / "plot.png")


def test_report_records_the_dense_count_and_both_certificate_verdicts(tmp_path: Path) -> None:
    before, after = _path((2, 3, 4, 5), closed=True), _path((2, 4, 5), closed=True)
    output = tmp_path / "plot.png"
    uncovered = CoverageVerdict.build(covered=False, witnesses=((3.0, 1.0),), piece_kinds=("annulus_outer",))
    covered = CoverageVerdict.build(covered=True, witnesses=(), piece_kinds=())
    render_coverage_preserving_comparison(
        load_held_reference_case("figure5"), before, after, (Radian(0.5), Radian(0.6)), output, before_coverage=uncovered, after_coverage=covered, dense_circles=7
    )
    report = json.loads(output.with_suffix(".json").read_text())
    assert report["dense_circles"] == 7
    assert report["before_coverage"] == {"covered": False, "witnesses": 1, "sites": 1}
    assert report["after_coverage"] == {"covered": True, "witnesses": 0, "sites": 0}
