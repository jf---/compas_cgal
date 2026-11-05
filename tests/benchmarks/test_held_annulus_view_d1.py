"""The D1 harness (engine == Stock2 on the circle view) of the actual-stock refinement, on Figure 5's first circles, plus its report contract."""

from __future__ import annotations

import json

import pytest

from benchmarks.held_actual_stock_last16 import CIRCLE
from benchmarks.held_actual_stock_last16 import SweepView
from benchmarks.held_annulus_view_d1 import D1Report
from benchmarks.held_annulus_view_d1 import D1Row
from benchmarks.held_annulus_view_d1 import EngineStockMismatchError
from benchmarks.held_annulus_view_d1 import d1_report
from benchmarks.held_annulus_view_d1 import engine_view
from benchmarks.held_annulus_view_d1 import part_rows
from benchmarks.held_annulus_view_d1 import view_range
from benchmarks.held_actual_stock_last16 import Last16Oracle
from benchmarks.held_actual_stock_oracle import differential_workload
from benchmarks.held_material_domain import PartPolygon
from benchmarks.held_material_domain import require_part_material


@pytest.fixture(scope="module")
def figure5():  # noqa: ANN201
    return differential_workload("figure5_refined")


def test_view_range_is_every_prior_circle_and_never_the_circle_itself() -> None:
    assert view_range(0, None) == (0, 0)
    assert view_range(7, None) == (0, 7)
    assert view_range(7, 1) == (6, 7)
    assert view_range(7, 16) == (0, 7)
    assert view_range(40, 16) == (24, 40)


def test_engine_view_stores_every_circle(figure5) -> None:  # noqa: ANN001
    assert engine_view(figure5).size == len(figure5.circles)


def test_first_figure5_circles_equal_the_arrangement(figure5) -> None:  # noqa: ANN001
    report = d1_report(figure5, range(25), last=None)
    report.require()
    assert report.stations == 25 * 45
    assert report.verdict_mismatches == 0
    assert report.exact_run_mismatches == 0


def test_circle_zero_is_refused_on_the_part(figure5) -> None:  # noqa: ANN001
    report = d1_report(figure5, [0], last=None)
    assert report.rows[0].exceeded


def test_engine_reads_the_part_polygon(figure5) -> None:  # noqa: ANN001
    require_part_material(part_rows(figure5), PartPolygon.from_stock_rows(figure5.polygon), figure5.name)


def test_d1_reference_is_the_circle_view_on_the_part(figure5) -> None:  # noqa: ANN001
    oracle = Last16Oracle(figure5, sweeps=SweepView.CIRCLES)
    assert all(kind == CIRCLE for _owner, kind, _item in oracle.candidates(3, None))
    assert [item for _owner, _kind, item in oracle.candidates(3, None)] == [0, 1, 2]


def _row(verdict: int = 0, runs: int = 0, reported: float = 0.0) -> D1Row:
    return D1Row(circle=3, stations=45, verdict_mismatches=verdict, exact_run_mismatches=runs, max_reported_delta_rad=reported, exceeded=False, engine_ms=1.0, reference_ms=2.0)


def _report(*rows: D1Row) -> D1Report:
    return D1Report(workload="w", last=None, load_before=(0.0, 0.0, 0.0), load_after=(0.0, 0.0, 0.0), rows=rows)


def test_report_require_passes_when_equal() -> None:
    _report(_row()).require()


@pytest.mark.parametrize("row", [_row(verdict=1), _row(runs=1), _row(reported=1.0)])
def test_report_require_raises_on_any_mismatch(row: D1Row) -> None:
    with pytest.raises(EngineStockMismatchError):
        _report(row).require()


def test_report_json_round_trips() -> None:
    record = json.loads(json.dumps(_report(_row(), _row()).to_json()))
    assert record["schema"] == "held-annulus-view-d1/1"
    assert record["circles"] == 2
    assert record["stations"] == 90
    assert record["verdict_mismatches"] == 0


# The stations where a cell between two near-coincident events first broke the witness construction (circle, station).
THIN_CELL_STATIONS = ((1159, 17), (1189, 28), (1191, 23), (1272, 23), (1277, 25), (1365, 24), (1449, 24))


@pytest.mark.parametrize(("circle", "station"), THIN_CELL_STATIONS)
def test_thin_cell_stations_resolve_and_equal_the_arrangement(figure5, circle: int, station: int) -> None:  # noqa: ANN001
    report = d1_report(figure5, [circle], last=None)
    report.require()
    assert report.stations == 45
    assert 0 <= station < 45
