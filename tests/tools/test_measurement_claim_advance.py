from __future__ import annotations

import importlib
import math
from typing import Any

import pytest

from compas.geometry import Polygon
from compas_cgal import engagement_toolpath


PLACEMENT_CAPS = (80.0, 100.0)
COUNT_CAPS = (40.0, 80.0)
COUNT_CONFIGURATIONS = (
    ("3-old", 3, (-60.0, 0.0, 60.0)),
    ("8", 8, tuple(45.0 * index for index in range(8))),
    ("12", 12, tuple(30.0 * index for index in range(12))),
    ("16", 16, tuple(22.5 * index for index in range(16))),
    ("24", 24, tuple(15.0 * index for index in range(24))),
    ("32", 32, tuple(11.25 * index for index in range(32))),
    ("40", 40, tuple(9.0 * index for index in range(40))),
    ("48", 48, tuple(7.5 * index for index in range(48))),
)
PLACEMENT_AUDIT_OFFSETS = tuple(11.25 * index for index in range(32))
COUNT_AUDIT_OFFSETS = tuple(3.0 + 6.0 * index for index in range(60))
PLACEMENT_SIGNED_AUDIT_OFFSETS = tuple(offset if offset <= 180.0 else offset - 360.0 for offset in PLACEMENT_AUDIT_OFFSETS)
POLYGON = (
    (0.0, 0.0, 0.0),
    (20.0, 0.0, 0.0),
    (20.0, 12.0, 0.0),
    (0.0, 12.0, 0.0),
)


def _module() -> Any:
    return importlib.import_module("tools.measurement_claim_advance")


def _polygon_points(polygon: Polygon) -> tuple[tuple[float, float, float], ...]:
    return tuple((float(point.x), float(point.y), float(point.z)) for point in polygon.points)


class _FakeStock:
    def __init__(self) -> None:
        self.raw = self
        self.bridge_removed = False


class _AdvanceRuntime:
    def __init__(self) -> None:
        self.generator_calls: list[dict[str, object]] = []
        self.original_calls: list[tuple[object, ...]] = []
        self.engagement_calls: list[tuple[_FakeStock, float, float]] = []
        self.raise_from_generator = False
        self.force_selection = False

    def original_search(
        self,
        stock: object,
        stations: object,
        origin_index: int,
        window_end: int,
        tool_radius: float,
        cap_ratio: float,
    ) -> tuple[int, bool]:
        self.original_calls.append((stock, stations, origin_index, window_end, tool_radius, cap_ratio))
        return 1, self.force_selection

    def generator(self, **kwargs: object) -> object:
        self.generator_calls.append(kwargs)
        if self.raise_from_generator:
            raise RuntimeError("generator failed")
        stock = _FakeStock()
        stations = [
            engagement_toolpath._GuideStation(cx=0.0, cy=0.0, radius=1.0, clockwise=bool(kwargs["climb"]), tx=1.0, ty=0.0),
            engagement_toolpath._GuideStation(cx=2.0, cy=0.0, radius=1.0, clockwise=bool(kwargs["climb"]), tx=1.0, ty=0.0),
        ]
        result = engagement_toolpath._largest_admissible_advance(stock, stations, 0, 1, 1.0, 2.0)
        stock.bridge_removed = True
        return result

    def engagement_at(
        self,
        stock: _FakeStock,
        x: float,
        y: float,
        tool_radius: float,
        cap_ratio: float,
        gap_close_ratio: float,
    ) -> tuple[float, float, bool]:
        assert stock.bridge_removed is False
        assert tool_radius == 1.0
        assert cap_ratio == 2.0
        assert gap_close_ratio == 0.0
        self.engagement_calls.append((stock, x, y))
        offset = round(math.degrees(math.atan2(y, x - 2.0)) % 360.0, 9)
        # Deliberately contradict angle and verdict so a reporting-angle
        # comparison cannot accidentally masquerade as the native decision.
        return 0.0, math.radians(179.0 if offset == 0.0 else 1.0), offset in (90.0, 93.0)


def _install_runtime(monkeypatch: pytest.MonkeyPatch, runtime: _AdvanceRuntime) -> None:
    module = _module()
    monkeypatch.setattr(module, "engagement_controlled_toolpath", runtime.generator)
    monkeypatch.setattr(module._stock_2, "engagement_at", runtime.engagement_at)
    monkeypatch.setattr(engagement_toolpath, "_largest_admissible_advance", runtime.original_search)


def _assert_public_call(call: dict[str, object], *, cap: float, climb: bool) -> None:
    assert set(call) == {
        "polygon",
        "holes",
        "tool_diameter",
        "tea_cap_deg",
        "guide_step_tool_diameters",
        "max_advance_tool_diameters",
        "radial_clearance",
        "climb",
        "cut_z",
        "clearance_z",
        "max_passes",
        "samples_per_radian",
    }
    assert _polygon_points(call["polygon"]) == POLYGON
    assert call["holes"] == []
    assert call["tool_diameter"] == 2.0
    assert call["tea_cap_deg"] == cap
    assert call["guide_step_tool_diameters"] == 0.025
    assert call["max_advance_tool_diameters"] == 1.0
    assert call["radial_clearance"] == 0.002
    assert call["climb"] is climb
    assert call["cut_z"] == 0.0
    assert call["clearance_z"] == 2.0
    assert call["max_passes"] == 1000
    assert call["samples_per_radian"] == 10.0


def _assert_advance_provenance(payload: dict[str, object]) -> None:
    expected = {
        "policy": "advance-native-cap-reporting-observation/v1",
        "native_cap_decision_site": "compas_cgal._stock_2.engagement_at[2]",
        "engagement_reporting_value_site": "compas_cgal._stock_2.engagement_at[1]",
        "reporting_driven_decision_sites": [],
    }
    assert payload["selection_decision_provenance"] == expected


def test_advance_case_order_is_fixed() -> None:
    assert _module().ADVANCE_CASE_ORDER == ("advance-placement", "advance-probe-count")


def test_advance_placement_runs_complete_calls_and_audits_pre_bridge_stock(monkeypatch: pytest.MonkeyPatch) -> None:
    module = _module()
    runtime = _AdvanceRuntime()
    _install_runtime(monkeypatch, runtime)
    payload = module.run_advance_case("advance-placement")

    assert payload["case"] == "advance-placement"
    assert payload["source_claim_ids"] == ["MC-009"]
    assert payload["continuous_certificate"] is None
    _assert_advance_provenance(payload)
    assert [(call["tea_cap_deg"], call["climb"]) for call in runtime.generator_calls] == [
        (80.0, True),
        (80.0, False),
        (100.0, True),
        (100.0, False),
    ]
    for call, (cap, climb) in zip(runtime.generator_calls, ((80.0, True), (80.0, False), (100.0, True), (100.0, False))):
        _assert_public_call(call, cap=cap, climb=climb)

    config = payload["config"]
    assert config["public_calls"] == [
        {
            "polygon": POLYGON,
            "holes": [],
            "tool_diameter": 2.0,
            "tea_cap": cap,
            "guide_step": 0.025,
            "max_advance": 1.0,
            "radial_clearance": 0.002,
            "climb": climb,
            "cut_z": 0.0,
            "clearance_z": 2.0,
            "max_passes": 1000,
            "samples_per_radian": 10.0,
        }
        for cap, climb in ((80.0, True), (80.0, False), (100.0, True), (100.0, False))
    ]
    assert len(config["public_calls"]) == 4
    assert config["generator_probe_angles"] == [-60.0, 0.0, 60.0]
    assert config["generator_prepends_entry_probe"] is True
    assert config["audit"] == {
        "phase": "advance-direction",
        "probe_offsets": list(PLACEMENT_AUDIT_OFFSETS),
        "excludes_entry_probe": True,
    }
    assert payload["reconstruction"] == {
        "stock_model": "generator-pre-bridge/v1",
        "selected_circle_policy": "accepted-non-forced/v1",
        "original_calls_per_wrapper": 1,
        "generator_prepends_entry_probe": True,
        "audit_excludes_entry_probe": True,
        "selected_circle_count": 4,
    }
    assert len(runtime.original_calls) == 4
    assert all(call[2:] == (0, 1, 1.0, 2.0) for call in runtime.original_calls)

    native = payload["native_sampled_decisions"]
    assert len(native) == 4
    assert native[0] == {
        "cap": 80.0,
        "climb": True,
        "selected_circles": 1,
        "observations": {"observations": 32, "accepted": 31, "exceeded": 1},
    }
    reporting = payload["reporting_values"]
    assert reporting[0]["positions_over_cap"] == 1
    assert reporting[0]["worst_peak"] == 179.0
    assert reporting[0]["worst_peak_offset"] == 0.0
    assert reporting[0]["old_probe_peak"] == 179.0
    assert reporting[0]["angle_unit"] == "degree"
    raw_offset_rows = reporting[0]["offset_bin_counts"]
    # Preserve all 32 native samples: the historical 30-degree boundary policy
    # is unrecoverable, so the corrected artifact must not infer sector counts.
    assert [row["offset"] for row in raw_offset_rows] == list(PLACEMENT_SIGNED_AUDIT_OFFSETS)
    assert raw_offset_rows[8] == {
        "offset": 90.0,
        "positions_over_cap": 1,
        "angle_unit": "degree",
    }
    assert sum(row["positions_over_cap"] for row in raw_offset_rows) == reporting[0]["positions_over_cap"]
    assert reporting[0]["positions_over_cap"] == native[0]["observations"]["exceeded"]
    assert native[0]["observations"]["observations"] == 32
    # 32 independent audit positions plus the separately reported old triple.
    assert len(runtime.engagement_calls) == 4 * 35


def test_advance_probe_count_uses_exact_generator_and_half_step_audit_angles(monkeypatch: pytest.MonkeyPatch) -> None:
    module = _module()
    runtime = _AdvanceRuntime()
    _install_runtime(monkeypatch, runtime)
    payload = module.run_advance_case("advance-probe-count")

    assert payload["case"] == "advance-probe-count"
    assert payload["source_claim_ids"] == ["MC-010"]
    assert payload["continuous_certificate"] is None
    _assert_advance_provenance(payload)
    config = payload["config"]
    assert config["public_calls"] == [
        {
            "polygon": POLYGON,
            "holes": [],
            "tool_diameter": 2.0,
            "tea_cap": cap,
            "guide_step": 0.025,
            "max_advance": 1.0,
            "radial_clearance": 0.002,
            "climb": True,
            "cut_z": 0.0,
            "clearance_z": 2.0,
            "max_passes": 1000,
            "samples_per_radian": 10.0,
        }
        for cap in COUNT_CAPS
    ]
    assert config["generator_configurations"] == [{"label": label, "probe_count": count, "probe_angles": list(angles)} for label, count, angles in COUNT_CONFIGURATIONS]
    assert config["generator_prepends_entry_probe"] is True
    assert config["audit"] == {
        "phase": "advance-half-step",
        "probe_offsets": list(COUNT_AUDIT_OFFSETS),
        "excludes_entry_probe": True,
    }
    assert len(runtime.generator_calls) == 16
    assert [(call["tea_cap_deg"], call["climb"]) for call in runtime.generator_calls] == [(cap, True) for _label, _count, _angles in COUNT_CONFIGURATIONS for cap in COUNT_CAPS]
    for call in runtime.generator_calls:
        cap = call["tea_cap_deg"]
        assert type(cap) is float
        _assert_public_call(call, cap=cap, climb=True)
    assert len(runtime.original_calls) == 16
    assert len(runtime.engagement_calls) == 16 * 60
    assert payload["reconstruction"]["selected_circle_count"] == 16
    assert payload["native_sampled_decisions"][0] == {
        "label": "3-old",
        "probe_count": 3,
        "cap": 40.0,
        "selected_circles": 1,
        "observations": {"observations": 60, "accepted": 59, "exceeded": 1},
    }
    assert payload["reporting_values"][0] == {
        "label": "3-old",
        "probe_count": 3,
        "cap": 40.0,
        "worst_peak": 1.0,
        "angle_unit": "degree",
    }


@pytest.mark.parametrize("case", ["advance-placement", "advance-probe-count"])
def test_advance_runner_restores_owned_bindings_after_success(monkeypatch: pytest.MonkeyPatch, case: str) -> None:
    module = _module()
    runtime = _AdvanceRuntime()
    _install_runtime(monkeypatch, runtime)
    original_count = engagement_toolpath.LOOP_PROBE_COUNT
    original_angles = engagement_toolpath.LOOP_PROBE_ANGLES_DEG
    original_search = engagement_toolpath._largest_admissible_advance

    module.run_advance_case(case)

    assert engagement_toolpath.LOOP_PROBE_COUNT is original_count
    assert engagement_toolpath.LOOP_PROBE_ANGLES_DEG is original_angles
    assert engagement_toolpath._largest_admissible_advance is original_search


@pytest.mark.parametrize("case", ["advance-placement", "advance-probe-count"])
def test_advance_runner_restores_owned_bindings_after_exception(monkeypatch: pytest.MonkeyPatch, case: str) -> None:
    module = _module()
    runtime = _AdvanceRuntime()
    runtime.raise_from_generator = True
    _install_runtime(monkeypatch, runtime)
    original_count = engagement_toolpath.LOOP_PROBE_COUNT
    original_angles = engagement_toolpath.LOOP_PROBE_ANGLES_DEG
    original_search = engagement_toolpath._largest_admissible_advance

    with pytest.raises(RuntimeError, match="generator failed"):
        module.run_advance_case(case)

    assert engagement_toolpath.LOOP_PROBE_COUNT is original_count
    assert engagement_toolpath.LOOP_PROBE_ANGLES_DEG is original_angles
    assert engagement_toolpath._largest_admissible_advance is original_search


def test_advance_runner_rejects_unknown_case_before_generator(monkeypatch: pytest.MonkeyPatch) -> None:
    module = _module()
    runtime = _AdvanceRuntime()
    _install_runtime(monkeypatch, runtime)
    with pytest.raises(module.UnknownMeasurementClaimCaseError):
        module.run_advance_case("unknown")
    assert runtime.generator_calls == []


def test_advance_runner_fails_loud_when_a_configuration_selects_no_nonforced_circle(monkeypatch: pytest.MonkeyPatch) -> None:
    module = _module()
    runtime = _AdvanceRuntime()
    runtime.force_selection = True
    _install_runtime(monkeypatch, runtime)
    with pytest.raises(module.ProbeInstrumentationContractError, match="accepted non-forced"):
        module.run_advance_case("advance-placement")
