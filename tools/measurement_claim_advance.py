"""Reconstruct the two advance-generator measurement-claim cases.

The generator remains authoritative. This module observes accepted non-forced
circles at the read-only search seam, before the caller subtracts the bridge.
"""

from __future__ import annotations

import math
from typing import Callable
from typing import List
from typing import Literal
from typing import Optional
from typing import TypedDict

from compas.geometry import Polygon

from compas_cgal import _stock_2
from compas_cgal import engagement_toolpath
from compas_cgal.engagement_toolpath import engagement_controlled_toolpath
from compas_cgal.stock import Stock
from tools.measurement_claim_errors import ProbeInstrumentationContractError
from tools.measurement_claim_errors import UnknownMeasurementClaimCaseError
from tools.measurement_claim_schema import AdvanceOffsetBinCountPayload
from tools.measurement_claim_schema import AdvancePlacementAuditConfigPayload
from tools.measurement_claim_schema import AdvancePlacementCasePayload
from tools.measurement_claim_schema import AdvancePlacementConfigPayload
from tools.measurement_claim_schema import AdvancePlacementNativeRowPayload
from tools.measurement_claim_schema import AdvancePlacementReportingRowPayload
from tools.measurement_claim_schema import AdvanceProbeConfigurationPayload
from tools.measurement_claim_schema import AdvanceProbeCountAuditConfigPayload
from tools.measurement_claim_schema import AdvanceProbeCountCasePayload
from tools.measurement_claim_schema import AdvanceProbeCountConfigPayload
from tools.measurement_claim_schema import AdvanceProbeCountNativeRowPayload
from tools.measurement_claim_schema import AdvanceProbeCountReportingRowPayload
from tools.measurement_claim_schema import AdvanceReconstructionPayload
from tools.measurement_claim_schema import AdvanceSelectionDecisionProvenancePayload
from tools.measurement_claim_schema import Degrees
from tools.measurement_claim_schema import GeneratorCasePayload
from tools.measurement_claim_schema import GeneratorPublicCallPayload
from tools.measurement_claim_schema import Millimetres
from tools.measurement_claim_schema import NativeObservationCountsPayload
from tools.measurement_claim_schema import ProbeConfigurationLabel
from tools.measurement_claim_schema import SamplesPerRadian
from tools.measurement_claim_schema import ToolDiameters
from tools.measurement_claim_schema import WorldMillimetres
from tools.measurement_claim_schema import WorldPointMillimetres
from tools.measurement_claim_schema import WorldRectangleMillimetres

AdvanceCase = Literal["advance-placement", "advance-probe-count"]
ADVANCE_CASE_ORDER: tuple[AdvanceCase, ...] = ("advance-placement", "advance-probe-count")

_POCKET_WIDTH = Millimetres(20.0)
_POCKET_HEIGHT = Millimetres(12.0)
_TOOL_DIAMETER = Millimetres(2.0)
_GUIDE_STEP = ToolDiameters(0.025)
_MAX_ADVANCE = ToolDiameters(1.0)
_RADIAL_CLEARANCE = Millimetres(0.002)
_CUT_Z = WorldMillimetres(0.0)
_CLEARANCE_Z = WorldMillimetres(2.0)
_MAX_PASSES = 1000
_SAMPLES_PER_RADIAN = SamplesPerRadian(10.0)

_PLACEMENT_CAPS: tuple[Degrees, ...] = (Degrees(80.0), Degrees(100.0))
_COUNT_CAPS: tuple[Degrees, ...] = (Degrees(40.0), Degrees(80.0))
_OLD_PROBE_ANGLES: tuple[Degrees, ...] = (Degrees(-60.0), Degrees(0.0), Degrees(60.0))
_PLACEMENT_AUDIT_OFFSETS: tuple[Degrees, ...] = tuple(Degrees(360.0 * index / 32) for index in range(32))
_COUNT_AUDIT_OFFSETS: tuple[Degrees, ...] = tuple(Degrees(360.0 * (index + 0.5) / 60) for index in range(60))
_COUNT_CONFIGURATIONS: tuple[tuple[ProbeConfigurationLabel, int, tuple[Degrees, ...]], ...] = (
    ("3-old", 3, _OLD_PROBE_ANGLES),
    ("8", 8, tuple(Degrees(360.0 * index / 8) for index in range(8))),
    ("12", 12, tuple(Degrees(360.0 * index / 12) for index in range(12))),
    ("16", 16, tuple(Degrees(360.0 * index / 16) for index in range(16))),
    ("24", 24, tuple(Degrees(360.0 * index / 24) for index in range(24))),
    ("32", 32, tuple(Degrees(360.0 * index / 32) for index in range(32))),
    ("40", 40, tuple(Degrees(360.0 * index / 40) for index in range(40))),
    ("48", 48, tuple(Degrees(360.0 * index / 48) for index in range(48))),
)

_Search = Callable[[Stock, list[engagement_toolpath._GuideStation], int, int, float, float], tuple[int, bool]]


class _SelectedCircleObservation(TypedDict):
    reporting: list[Degrees]
    verdicts: list[bool]
    worst_peak: Degrees
    worst_peak_offset: Degrees
    old_probe_peak: Optional[Degrees]


def _polygon_payload() -> WorldRectangleMillimetres:
    lower_left: WorldPointMillimetres = (WorldMillimetres(0.0), WorldMillimetres(0.0), WorldMillimetres(0.0))
    width = WorldMillimetres(float(_POCKET_WIDTH))
    height = WorldMillimetres(float(_POCKET_HEIGHT))
    lower_right: WorldPointMillimetres = (width, WorldMillimetres(0.0), WorldMillimetres(0.0))
    upper_right: WorldPointMillimetres = (width, height, WorldMillimetres(0.0))
    upper_left: WorldPointMillimetres = (WorldMillimetres(0.0), height, WorldMillimetres(0.0))
    return lower_left, lower_right, upper_right, upper_left


def _polygon() -> Polygon:
    return Polygon(_polygon_payload())


def _public_call(cap: Degrees, *, climb: bool) -> GeneratorPublicCallPayload:
    payload: GeneratorPublicCallPayload = {
        "polygon": _polygon_payload(),
        "holes": [],
        "tool_diameter": _TOOL_DIAMETER,
        "tea_cap": cap,
        "guide_step": _GUIDE_STEP,
        "max_advance": _MAX_ADVANCE,
        "radial_clearance": _RADIAL_CLEARANCE,
        "climb": climb,
        "cut_z": _CUT_Z,
        "clearance_z": _CLEARANCE_Z,
        "max_passes": _MAX_PASSES,
        "samples_per_radian": _SAMPLES_PER_RADIAN,
    }
    return payload


def _run_generator(cap: Degrees, *, climb: bool) -> None:
    engagement_controlled_toolpath(
        polygon=_polygon(),
        holes=[],
        tool_diameter=float(_TOOL_DIAMETER),
        tea_cap_deg=float(cap),
        guide_step_tool_diameters=float(_GUIDE_STEP),
        max_advance_tool_diameters=float(_MAX_ADVANCE),
        radial_clearance=float(_RADIAL_CLEARANCE),
        climb=climb,
        cut_z=float(_CUT_Z),
        clearance_z=float(_CLEARANCE_Z),
        max_passes=_MAX_PASSES,
        samples_per_radian=_SAMPLES_PER_RADIAN,
    )


def _provenance() -> AdvanceSelectionDecisionProvenancePayload:
    payload: AdvanceSelectionDecisionProvenancePayload = {
        "policy": "advance-native-cap-reporting-observation/v1",
        "native_cap_decision_site": "compas_cgal._stock_2.engagement_at[2]",
        "engagement_reporting_value_site": "compas_cgal._stock_2.engagement_at[1]",
        "reporting_driven_decision_sites": [],
    }
    return payload


def _signed_offset(offset: Degrees) -> Degrees:
    value = float(offset)
    return Degrees(value - 360.0 if value > 180.0 else value)


def _position(
    station: engagement_toolpath._GuideStation,
    advance: tuple[float, float],
    offset: Degrees,
) -> tuple[float, float]:
    angle = math.radians(float(offset))
    dx, dy = advance
    cos_a = math.cos(angle)
    sin_a = math.sin(angle)
    ux = dx * cos_a - dy * sin_a
    uy = dx * sin_a + dy * cos_a
    return station.cx + station.radius * ux, station.cy + station.radius * uy


def _observe(
    stock: Stock,
    station: engagement_toolpath._GuideStation,
    advance: tuple[float, float],
    tool_radius: float,
    cap_ratio: float,
    offsets: tuple[Degrees, ...],
) -> tuple[list[Degrees], list[bool]]:
    reporting: list[Degrees] = []
    verdicts: list[bool] = []
    for offset in offsets:
        x, y = _position(station, advance, offset)
        observation = _stock_2.engagement_at(stock.raw, x, y, tool_radius, cap_ratio, 0.0)
        reported = observation[1]
        exceeded = observation[2]
        if type(reported) is not float or not math.isfinite(reported):
            raise ProbeInstrumentationContractError("advance audit reporting value must be one finite native float at engagement_at index 1")
        if type(exceeded) is not bool:
            raise ProbeInstrumentationContractError("advance audit cap verdict must be one native bool at engagement_at index 2")
        reporting.append(Degrees(math.degrees(reported)))
        verdicts.append(exceeded)
    return reporting, verdicts


def _observe_selected_circle(
    stock: Stock,
    stations: list[engagement_toolpath._GuideStation],
    origin_index: int,
    selected_index: int,
    tool_radius: float,
    cap_ratio: float,
    audit_offsets: tuple[Degrees, ...],
    old_probe_offsets: Optional[tuple[Degrees, ...]],
) -> _SelectedCircleObservation:
    origin = stations[origin_index]
    selected = stations[selected_index]
    advance = engagement_toolpath._unit_tangent(origin.cx, origin.cy, selected.cx, selected.cy)
    if advance == (0.0, 0.0):
        raise ProbeInstrumentationContractError("accepted advance has coincident origin and selected circle centres")
    reporting, verdicts = _observe(stock, selected, advance, tool_radius, cap_ratio, audit_offsets)
    worst_index = max(range(len(reporting)), key=reporting.__getitem__)
    old_peak: Optional[Degrees] = None
    if old_probe_offsets is not None:
        old_reporting, _old_verdicts = _observe(stock, selected, advance, tool_radius, cap_ratio, old_probe_offsets)
        old_peak = max(old_reporting)
    return {
        "reporting": reporting,
        "verdicts": verdicts,
        "worst_peak": reporting[worst_index],
        "worst_peak_offset": _signed_offset(audit_offsets[worst_index]),
        "old_probe_peak": old_peak,
    }


def _counts(observations: list[_SelectedCircleObservation]) -> NativeObservationCountsPayload:
    verdicts = [verdict for observation in observations for verdict in observation["verdicts"]]
    exceeded = sum(verdicts)
    payload: NativeObservationCountsPayload = {
        "observations": len(verdicts),
        "accepted": len(verdicts) - exceeded,
        "exceeded": exceeded,
    }
    return payload


def _worst_observation(observations: list[_SelectedCircleObservation]) -> _SelectedCircleObservation:
    return max(observations, key=lambda observation: float(observation["worst_peak"]))


def _reconstruction(selected_circle_count: int) -> AdvanceReconstructionPayload:
    payload: AdvanceReconstructionPayload = {
        "stock_model": "generator-pre-bridge/v1",
        "selected_circle_policy": "accepted-non-forced/v1",
        "original_calls_per_wrapper": 1,
        "generator_prepends_entry_probe": True,
        "audit_excludes_entry_probe": True,
        "selected_circle_count": selected_circle_count,
    }
    return payload


def _placement_payload(original_search: _Search) -> AdvancePlacementCasePayload:
    observed: list[_SelectedCircleObservation] = []
    active: list[_SelectedCircleObservation] = []

    def observe_after_search(
        stock: Stock,
        stations: list[engagement_toolpath._GuideStation],
        origin_index: int,
        window_end: int,
        tool_radius: float,
        cap_ratio: float,
    ) -> tuple[int, bool]:
        selected_index, forced = original_search(stock, stations, origin_index, window_end, tool_radius, cap_ratio)
        if not forced:
            active.append(
                _observe_selected_circle(
                    stock,
                    stations,
                    origin_index,
                    selected_index,
                    tool_radius,
                    cap_ratio,
                    _PLACEMENT_AUDIT_OFFSETS,
                    _OLD_PROBE_ANGLES,
                )
            )
        return selected_index, forced

    engagement_toolpath.LOOP_PROBE_COUNT = len(_OLD_PROBE_ANGLES)
    engagement_toolpath.LOOP_PROBE_ANGLES_DEG = tuple(float(angle) for angle in _OLD_PROBE_ANGLES)
    engagement_toolpath._largest_admissible_advance = observe_after_search
    native_rows: List[AdvancePlacementNativeRowPayload] = []
    reporting_rows: List[AdvancePlacementReportingRowPayload] = []
    for cap in _PLACEMENT_CAPS:
        for climb in (True, False):
            active.clear()
            _run_generator(cap, climb=climb)
            if not active:
                raise ProbeInstrumentationContractError(f"advance-placement cap={float(cap)!r}, climb={climb!r} selected no accepted non-forced circle")
            observed.extend(active)
            counts = _counts(active)
            worst = _worst_observation(active)
            offset_counts: List[AdvanceOffsetBinCountPayload] = []
            for index, offset in enumerate(_PLACEMENT_AUDIT_OFFSETS):
                exceeded = sum(observation["verdicts"][index] for observation in active)
                offset_row: AdvanceOffsetBinCountPayload = {
                    "offset": _signed_offset(offset),
                    "positions_over_cap": exceeded,
                    "angle_unit": "degree",
                }
                offset_counts.append(offset_row)
            native_row: AdvancePlacementNativeRowPayload = {
                "cap": cap,
                "climb": climb,
                "selected_circles": len(active),
                "observations": counts,
            }
            native_rows.append(native_row)
            old_probe_peak = worst["old_probe_peak"]
            if old_probe_peak is None:
                raise ProbeInstrumentationContractError("advance-placement observation omitted the old-probe reporting peak")
            reporting_row: AdvancePlacementReportingRowPayload = {
                "cap": cap,
                "climb": climb,
                "selected_circles": len(active),
                "positions_over_cap": counts["exceeded"],
                "worst_peak": worst["worst_peak"],
                "worst_peak_offset": worst["worst_peak_offset"],
                "old_probe_peak": old_probe_peak,
                "offset_bin_counts": offset_counts,
                "angle_unit": "degree",
            }
            reporting_rows.append(reporting_row)
    audit: AdvancePlacementAuditConfigPayload = {
        "phase": "advance-direction",
        "probe_offsets": [Degrees(offset) for offset in _PLACEMENT_AUDIT_OFFSETS],
        "excludes_entry_probe": True,
    }
    config: AdvancePlacementConfigPayload = {
        "public_calls": [_public_call(cap, climb=climb) for cap in _PLACEMENT_CAPS for climb in (True, False)],
        "generator_probe_angles": [Degrees(angle) for angle in _OLD_PROBE_ANGLES],
        "generator_prepends_entry_probe": True,
        "audit": audit,
    }
    payload: AdvancePlacementCasePayload = {
        "case": "advance-placement",
        "source_claim_ids": ["MC-009"],
        "config": config,
        "reconstruction": _reconstruction(len(observed)),
        "native_sampled_decisions": native_rows,
        "reporting_values": reporting_rows,
        "selection_decision_provenance": _provenance(),
        "continuous_certificate": None,
    }
    return payload


def _probe_count_payload(original_search: _Search) -> AdvanceProbeCountCasePayload:
    observed: list[_SelectedCircleObservation] = []
    active: list[_SelectedCircleObservation] = []

    def observe_after_search(
        stock: Stock,
        stations: list[engagement_toolpath._GuideStation],
        origin_index: int,
        window_end: int,
        tool_radius: float,
        cap_ratio: float,
    ) -> tuple[int, bool]:
        selected_index, forced = original_search(stock, stations, origin_index, window_end, tool_radius, cap_ratio)
        if not forced:
            active.append(
                _observe_selected_circle(
                    stock,
                    stations,
                    origin_index,
                    selected_index,
                    tool_radius,
                    cap_ratio,
                    _COUNT_AUDIT_OFFSETS,
                    None,
                )
            )
        return selected_index, forced

    engagement_toolpath._largest_admissible_advance = observe_after_search
    native_rows: List[AdvanceProbeCountNativeRowPayload] = []
    reporting_rows: List[AdvanceProbeCountReportingRowPayload] = []
    for label, probe_count, probe_angles in _COUNT_CONFIGURATIONS:
        engagement_toolpath.LOOP_PROBE_COUNT = probe_count
        engagement_toolpath.LOOP_PROBE_ANGLES_DEG = tuple(float(angle) for angle in probe_angles)
        for cap in _COUNT_CAPS:
            active.clear()
            _run_generator(cap, climb=True)
            if not active:
                raise ProbeInstrumentationContractError(f"advance-probe-count label={label!r}, cap={float(cap)!r} selected no accepted non-forced circle")
            observed.extend(active)
            counts = _counts(active)
            worst = _worst_observation(active)
            native_row: AdvanceProbeCountNativeRowPayload = {
                "label": label,
                "probe_count": probe_count,
                "cap": cap,
                "selected_circles": len(active),
                "observations": counts,
            }
            native_rows.append(native_row)
            reporting_row: AdvanceProbeCountReportingRowPayload = {
                "label": label,
                "probe_count": probe_count,
                "cap": cap,
                "worst_peak": worst["worst_peak"],
                "angle_unit": "degree",
            }
            reporting_rows.append(reporting_row)
    generator_configurations: List[AdvanceProbeConfigurationPayload] = []
    for label, probe_count, probe_angles in _COUNT_CONFIGURATIONS:
        configuration: AdvanceProbeConfigurationPayload = {
            "label": label,
            "probe_count": probe_count,
            "probe_angles": [Degrees(angle) for angle in probe_angles],
        }
        generator_configurations.append(configuration)
    audit: AdvanceProbeCountAuditConfigPayload = {
        "phase": "advance-half-step",
        "probe_offsets": [Degrees(offset) for offset in _COUNT_AUDIT_OFFSETS],
        "excludes_entry_probe": True,
    }
    config: AdvanceProbeCountConfigPayload = {
        "public_calls": [_public_call(cap, climb=True) for cap in _COUNT_CAPS],
        "generator_configurations": generator_configurations,
        "generator_prepends_entry_probe": True,
        "audit": audit,
    }
    payload: AdvanceProbeCountCasePayload = {
        "case": "advance-probe-count",
        "source_claim_ids": ["MC-010"],
        "config": config,
        "reconstruction": _reconstruction(len(observed)),
        "native_sampled_decisions": native_rows,
        "reporting_values": reporting_rows,
        "selection_decision_provenance": _provenance(),
        "continuous_certificate": None,
    }
    return payload


def run_advance_case(case: AdvanceCase) -> GeneratorCasePayload:
    """Run one fixed advance claim case and restore its instrumentation seam."""
    if case not in ADVANCE_CASE_ORDER:
        raise UnknownMeasurementClaimCaseError(f"unknown advance measurement-claim case: {case!r}")
    original_count = engagement_toolpath.LOOP_PROBE_COUNT
    original_angles = engagement_toolpath.LOOP_PROBE_ANGLES_DEG
    original_search = engagement_toolpath._largest_admissible_advance
    try:
        if case == "advance-placement":
            return _placement_payload(original_search)
        return _probe_count_payload(original_search)
    finally:
        engagement_toolpath.LOOP_PROBE_COUNT = original_count
        engagement_toolpath.LOOP_PROBE_ANGLES_DEG = original_angles
        engagement_toolpath._largest_admissible_advance = original_search
