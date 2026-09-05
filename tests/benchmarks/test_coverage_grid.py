from dataclasses import replace

import pytest
from compas.geometry import Polygon

import benchmarks.coverage as coverage_module
from benchmarks.errors import CoarseCoverageGridError
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.spec import PocketSpec


def test_figure5_minimum_coverage_grid_is_the_true_rounded_minimum() -> None:
    spec = load_held_reference_case("figure5").pocket_spec()

    assert coverage_module.minimum_coverage_grid(spec) == 668
    with pytest.raises(CoarseCoverageGridError):
        coverage_module._grid_axes(spec, 667)
    xs, ys, cell = coverage_module._grid_axes(spec, 668)
    assert (len(xs), len(ys)) == (668, 458)
    assert cell == pytest.approx(0.09993886008290782)
    assert cell <= coverage_module.MAX_CELL_TOOL_RADIUS_FRACTION * spec.tool_radius


def test_minimum_grid_accounts_for_short_axis_rounding() -> None:
    spec = PocketSpec.build(
        name="rounding-counterexample",
        family="analytic",
        polygon=Polygon(((0.0, 0.0), (100.0, 0.0), (100.0, 94.49), (0.0, 94.49))),
        tool_diameter=20.0,
        tea_cap_deg=120.0,
    )

    simple_ceiling = 100
    with pytest.raises(CoarseCoverageGridError):
        coverage_module._grid_axes(spec, simple_ceiling)
    assert coverage_module.minimum_coverage_grid(spec) == 101
    xs, ys, cell = coverage_module._grid_axes(spec, 101)
    assert (len(xs), len(ys)) == (101, 95)
    assert cell <= 1.0


def test_minimum_grid_preserves_the_default_cost_floor_for_small_cases() -> None:
    figure5 = load_held_reference_case("figure5").pocket_spec()
    small = replace(figure5, polygon=Polygon(((-2.0, -2.0), (2.0, -2.0), (2.0, 2.0), (-2.0, 2.0))))

    assert coverage_module.minimum_coverage_grid(small) == 40


def test_measure_coverage_builds_one_material_predicate_and_queries_each_cell(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    figure5 = load_held_reference_case("figure5").pocket_spec()
    spec = replace(
        figure5,
        polygon=Polygon(((-2.0, -2.0), (2.0, -2.0), (2.0, 2.0), (-2.0, 2.0))),
    )
    calls = {"build": 0, "contains": 0}

    class Predicate:
        @classmethod
        def build(cls, boundary: object, holes: object, radius: float) -> "Predicate":
            calls["build"] += 1
            return cls()

        def contains(self, x: float, y: float) -> bool:
            calls["contains"] += 1
            return True

    class Stock:
        def contains(self, x: float, y: float) -> bool:
            return True

    monkeypatch.setattr(
        coverage_module._coverage_2,
        "ReachableMaterialPredicate2",
        Predicate,
        raising=False,
    )

    estimate = coverage_module.measure_coverage(spec, Stock(), grid=40)  # type: ignore[arg-type]

    assert calls == {"build": 1, "contains": 1600}
    assert estimate.reachable_samples == 1600
    assert estimate.uncut_reachable_samples == 1600
    assert estimate.remaining_samples == 1600
