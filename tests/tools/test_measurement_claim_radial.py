from __future__ import annotations

import dataclasses
import importlib
import importlib.util
import math
from types import SimpleNamespace
from typing import Any

import pytest
from compas.geometry import Circle
from compas.geometry import Line

from compas_cgal import engagement_radial_toolpath
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation


def _module() -> Any:
    return importlib.import_module("tools.measurement_claim_radial")


def test_measurement_claim_radial_module_exists() -> None:
    assert importlib.util.find_spec("tools.measurement_claim_radial") is not None


def test_radial_case_order_is_fixed() -> None:
    assert _module().RADIAL_CASE_ORDER == (
        "radial-station",
        "radial-subdivisions",
        "radial-floor",
        "radial-margin",
    )


def test_radial_runner_rejects_unknown_case_before_generator(monkeypatch: pytest.MonkeyPatch) -> None:
    module = _module()
    calls: list[object] = []
    monkeypatch.setattr(module, "radius_regulated_toolpath", lambda **kwargs: calls.append(kwargs))
    with pytest.raises(module.UnknownMeasurementClaimCaseError):
        module.run_radial_case("unknown")
    assert calls == []


def test_radial_public_call_records_every_effective_argument_without_defaults() -> None:
    module = _module()

    public_call = module._public_call(module.Millimetres(20.0), module.Millimetres(12.0), module.Degrees(60.0))

    assert public_call == {
        "polygon": ((0.0, 0.0, 0.0), (20.0, 0.0, 0.0), (20.0, 12.0, 0.0), (0.0, 12.0, 0.0)),
        "holes": [],
        "tool_diameter": 2.0,
        "tea_cap": 60.0,
        "guide_step": 0.025,
        "max_advance": 1.0,
        "radial_clearance": 0.002,
        "climb": True,
        "cut_z": 0.0,
        "clearance_z": 2.0,
        "max_passes": 1000,
        "samples_per_radian": 10.0,
    }


def test_radial_provenance_resolves_real_decision_symbols_and_peak_carrier() -> None:
    module = _module()

    provenance = module._provenance()

    assert provenance["native_cap_decision_site"] == "compas_cgal._stock_2.engagement_at[2]"
    assert provenance["engagement_reporting_value_site"] == "compas_cgal._stock_2.engagement_at[1]"
    assert provenance["reporting_driven_decision_sites"] == [
        {
            "symbol": "compas_cgal.engagement_radial_toolpath._least_bad_rung",
            "value_source": "compas_cgal._stock_2.engagement_at[1]",
            "effect": "forced-radius-selection",
        },
        {
            "symbol": "compas_cgal.engagement_radial_toolpath._largest_admissible_radius",
            "value_source": "compas_cgal.engagement_radial_toolpath._GentlestRung.peak",
            "effect": "refined-scan-control",
        },
    ]
    assert callable(engagement_radial_toolpath._least_bad_rung)
    assert callable(engagement_radial_toolpath._largest_admissible_radius)
    assert callable(module._stock_2.engagement_at)
    assert [field.name for field in dataclasses.fields(engagement_radial_toolpath._GentlestRung)] == ["rung", "peak"]


def test_radial_audit_uses_exactly_sixteen_entry_phased_positions_and_separate_native_and_reporting_indices(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    module = _module()
    station = engagement_radial_toolpath._GuideStation(cx=10.0, cy=20.0, radius=2.0, clockwise=True, tx=1.0, ty=0.0)
    regulation = SimpleNamespace(tool_radius=1.0, cap_ratio=0.75)
    stock = SimpleNamespace(raw=object())
    positions: list[tuple[float, float, float]] = []

    def engagement_at(raw: object, x: float, y: float, tool_radius: float, cap_ratio: float, z: float) -> tuple[bool, float, bool]:
        assert raw is stock.raw
        assert tool_radius == 1.0
        assert cap_ratio == 0.75
        positions.append((x, y, z))
        return True, math.pi / 6.0, False

    monkeypatch.setattr(module._stock_2, "engagement_at", engagement_at)
    monkeypatch.setattr(engagement_radial_toolpath, "_loop_reaches_material", lambda *args: True)

    observed = module._audit_circle(stock, station, 2.0, (1.0, 0.0), regulation)

    assert len(positions) == 16
    assert positions[0] == pytest.approx((*station.entry, 0.0))
    assert len({(round(x, 12), round(y, 12)) for x, y, _ in positions}) == 16
    assert observed["reporting"] == pytest.approx([30.0] * 16)
    assert observed["verdicts"] == [False] * 16
    assert observed["cuts_material"] is True


def test_radial_wrapper_observes_once_before_bridge_mutation_and_excludes_each_chain_entry(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    module = _module()
    stock = SimpleNamespace(raw=object())
    regulation = SimpleNamespace(tool_radius=1.0, cap_ratio=0.75)
    station = engagement_radial_toolpath._GuideStation(cx=1.0, cy=2.0, radius=0.5, clockwise=True, tx=1.0, ty=0.0)
    bridge_removed = False
    original_calls: list[tuple[object, ...]] = []
    audited_pre_bridge: list[bool] = []

    def original(*args: object) -> Any:
        original_calls.append(args)
        return engagement_radial_toolpath._RadiusChoice.of(station.radius, station, forced=False)

    def audit(*args: object) -> dict[str, object]:
        del args
        audited_pre_bridge.append(not bridge_removed)
        return {
            "centre": (1.0, 2.0),
            "maximal_radius": 0.5,
            "selected_radius": 0.5,
            "reporting": [30.0] * 16,
            "verdicts": [False] * 16,
            "cuts_material": True,
        }

    def run_generator(*args: object) -> Any:
        nonlocal bridge_removed
        del args
        wrapped = engagement_radial_toolpath._largest_admissible_radius
        for _ in range(3):
            wrapped(stock, station, 0.0, (1.0, 0.0), regulation)
            bridge_removed = True
            bridge_removed = False
        return SimpleNamespace(
            operations=[
                ToolpathOperation(Circle(0.5), OperationType.CUT, 4),
                ToolpathOperation(Circle(0.5), OperationType.CUT, 4),
                ToolpathOperation(Circle(0.5), OperationType.CUT, 9),
            ]
        )

    monkeypatch.setattr(module, "_audit_circle", audit)
    monkeypatch.setattr(module, "_run_generator", run_generator)

    captured = module._capture_run(module.Millimetres(20.0), module.Millimetres(12.0), module.Degrees(60.0), original)

    assert len(original_calls) == 3
    assert audited_pre_bridge == [True, True, True]
    assert len(captured["all_circles"]) == 3
    assert captured["non_entry_circles"] == [captured["all_circles"][1]]


def test_radial_station_repeats_exact_serialized_refined_radii_in_native_rows(monkeypatch: pytest.MonkeyPatch) -> None:
    module = _module()

    def observation(radius: float, peak: float, *, exceeded: bool = False) -> dict[str, object]:
        return {
            "centre": (18.482, 10.482),
            "maximal_radius": 0.5156,
            "selected_radius": radius,
            "reporting": [peak] * 16,
            "verdicts": [exceeded] * 16,
            "cuts_material": True,
        }

    selected = {**observation(0.212345, 59.0), "forced": False, "finishes": False}
    station = {
        "centre": (18.482, 10.482),
        "maximal_radius": 0.515555,
        "coarse_step": 0.05,
        "rung_6": observation(0.2156, 61.3, exceeded=True),
        "rung_7": observation(0.1656, 5.9),
        "refined": [observation(0.515555, 107.0, exceeded=True), observation(0.212345, 59.0)],
        "selected": selected,
    }
    monkeypatch.setattr(
        module,
        "_capture_run",
        lambda *args: {"result": SimpleNamespace(operations=[]), "all_circles": [], "non_entry_circles": [], "stations": [station]},
    )

    payload = module._station_payload(lambda *args: None)

    reconstructed = payload["reconstruction"]["refined_radius_sequence"]
    native = [row["radius"] for row in payload["native_sampled_decisions"]["refined_candidates"]]
    assert native == reconstructed == [0.5156, 0.2123]


def test_radial_cutting_length_is_cut_circle_circumference_plus_cut_line_xy_length() -> None:
    module = _module()
    result = SimpleNamespace(
        operations=[
            ToolpathOperation(Circle(1.0), OperationType.CUT, 0),
            ToolpathOperation(Line([0.0, 0.0, 0.0], [3.0, 4.0, 0.0]), OperationType.CUT, 0),
            ToolpathOperation(Line([0.0, 0.0, 0.0], [100.0, 0.0, 0.0]), OperationType.LINK, 0),
        ]
    )

    assert module._cutting_length(result) == round(2.0 * math.pi + 5.0)


def test_radial_sweep_reconstruction_preserves_each_configuration_count() -> None:
    module = _module()
    circle = {
        "centre": (1.0, 2.0),
        "maximal_radius": 0.5,
        "selected_radius": 0.4,
        "reporting": [30.0] * 16,
        "verdicts": [False] * 16,
        "cuts_material": True,
        "forced": False,
        "finishes": False,
    }

    reconstruction = module._sweep_reconstruction([[circle], [circle, circle], []])

    assert reconstruction["non_entry_circle_counts"] == [1, 2, 0]


@pytest.mark.parametrize(
    ("binding", "value"),
    (
        ("LOOP_PROBE_COUNT", 31),
        ("LOOP_PROBE_ANGLES_DEG", (0.0,)),
    ),
)
def test_radial_runner_rejects_drifted_shared_probe_configuration_before_generation(
    monkeypatch: pytest.MonkeyPatch,
    binding: str,
    value: object,
) -> None:
    module = _module()
    calls: list[object] = []
    monkeypatch.setattr(module.engagement_toolpath, binding, value)
    monkeypatch.setattr(module, "radius_regulated_toolpath", lambda **kwargs: calls.append(kwargs))

    with pytest.raises(module.InvalidMeasurementClaimConfigError, match="probe ring"):
        module.run_radial_case("radial-station")

    assert calls == []


@pytest.mark.parametrize(
    ("binding", "value"),
    (
        ("RADIUS_LADDER_RUNGS", 39),
        ("MAX_RADIAL_SWEEPS_PER_CHAIN", 39),
    ),
)
def test_radial_runner_rejects_drifted_ladder_budget_before_generation(
    monkeypatch: pytest.MonkeyPatch,
    binding: str,
    value: int,
) -> None:
    module = _module()
    calls: list[object] = []
    monkeypatch.setattr(engagement_radial_toolpath, binding, value)
    monkeypatch.setattr(module, "radius_regulated_toolpath", lambda **kwargs: calls.append(kwargs))

    with pytest.raises(module.InvalidMeasurementClaimConfigError, match="rung and sweep budgets"):
        module.run_radial_case("radial-station")

    assert calls == []


@pytest.mark.parametrize(
    "case",
    ("radial-station", "radial-subdivisions", "radial-floor", "radial-margin"),
)
def test_radial_runner_restores_owned_bindings_after_failure(
    monkeypatch: pytest.MonkeyPatch,
    case: str,
) -> None:
    module = _module()
    original_subdivisions = engagement_radial_toolpath.RADIUS_LADDER_SUBDIVISIONS
    original_floor = engagement_radial_toolpath.RADIUS_LADDER_FLOOR_STEPS
    original_margin = engagement_radial_toolpath.RADIUS_LADDER_REFINEMENT_MARGIN
    original_radius = engagement_radial_toolpath._largest_admissible_radius

    def fail(**kwargs: object) -> None:
        del kwargs
        raise RuntimeError("generator failed")

    monkeypatch.setattr(module, "radius_regulated_toolpath", fail)
    with pytest.raises(RuntimeError, match="generator failed"):
        module.run_radial_case(case)

    assert engagement_radial_toolpath.RADIUS_LADDER_SUBDIVISIONS == original_subdivisions
    assert engagement_radial_toolpath.RADIUS_LADDER_FLOOR_STEPS == original_floor
    assert engagement_radial_toolpath.RADIUS_LADDER_REFINEMENT_MARGIN == original_margin
    assert engagement_radial_toolpath._largest_admissible_radius is original_radius


def test_radial_runner_restores_owned_bindings_after_success(monkeypatch: pytest.MonkeyPatch) -> None:
    module = _module()
    original_subdivisions = engagement_radial_toolpath.RADIUS_LADDER_SUBDIVISIONS
    original_floor = engagement_radial_toolpath.RADIUS_LADDER_FLOOR_STEPS
    original_margin = engagement_radial_toolpath.RADIUS_LADDER_REFINEMENT_MARGIN
    original_radius = engagement_radial_toolpath._largest_admissible_radius
    sentinel = object()

    def succeed(original_search: object) -> object:
        del original_search
        engagement_radial_toolpath.RADIUS_LADDER_SUBDIVISIONS = 99
        engagement_radial_toolpath.RADIUS_LADDER_FLOOR_STEPS = 99.0
        engagement_radial_toolpath.RADIUS_LADDER_REFINEMENT_MARGIN = 99.0
        engagement_radial_toolpath._largest_admissible_radius = lambda *args: None
        return sentinel

    monkeypatch.setattr(module, "_station_payload", succeed)

    assert module.run_radial_case("radial-station") is sentinel
    assert engagement_radial_toolpath.RADIUS_LADDER_SUBDIVISIONS == original_subdivisions
    assert engagement_radial_toolpath.RADIUS_LADDER_FLOOR_STEPS == original_floor
    assert engagement_radial_toolpath.RADIUS_LADDER_REFINEMENT_MARGIN == original_margin
    assert engagement_radial_toolpath._largest_admissible_radius is original_radius
