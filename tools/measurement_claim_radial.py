"""Reconstruct the four radial-generator measurement-claim cases."""

from __future__ import annotations

import math
from dataclasses import replace
from typing import Callable
from typing import List
from typing import Literal
from typing import Optional
from typing import Sequence
from typing import Tuple
from typing import TypedDict

from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Line
from compas.geometry import Polygon

from compas_cgal import _stock_2
from compas_cgal import engagement_radial_toolpath
from compas_cgal import engagement_toolpath
from compas_cgal.engagement_radial_toolpath import radius_regulated_toolpath
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathResult
from tools.measurement_claim_errors import InvalidMeasurementClaimConfigError
from tools.measurement_claim_errors import ProbeInstrumentationContractError
from tools.measurement_claim_errors import UnknownMeasurementClaimCaseError
from tools.measurement_claim_schema import CoarseSteps
from tools.measurement_claim_schema import Degrees
from tools.measurement_claim_schema import DimensionlessRatio
from tools.measurement_claim_schema import GeneratorCasePayload
from tools.measurement_claim_schema import GeneratorPublicCallPayload
from tools.measurement_claim_schema import Millimetres
from tools.measurement_claim_schema import NativeObservationCountsPayload
from tools.measurement_claim_schema import RadialAuditConfigPayload
from tools.measurement_claim_schema import RadialFloorCasePayload
from tools.measurement_claim_schema import RadialFloorConfigPayload
from tools.measurement_claim_schema import RadialFloorNativeRowPayload
from tools.measurement_claim_schema import RadialFloorReportingRowPayload
from tools.measurement_claim_schema import RadialMarginCasePayload
from tools.measurement_claim_schema import RadialMarginConfigPayload
from tools.measurement_claim_schema import RadialMarginNativeRowPayload
from tools.measurement_claim_schema import RadialMarginReportingRowPayload
from tools.measurement_claim_schema import RadialSelectionDecisionProvenancePayload
from tools.measurement_claim_schema import RadialStationCasePayload
from tools.measurement_claim_schema import RadialStationConfigPayload
from tools.measurement_claim_schema import RadialStationNativeDecisionsPayload
from tools.measurement_claim_schema import RadialStationReconstructionPayload
from tools.measurement_claim_schema import RadialStationRefinedCandidateNativePayload
from tools.measurement_claim_schema import RadialStationReportingPayload
from tools.measurement_claim_schema import RadialSubdivisionNativeRowPayload
from tools.measurement_claim_schema import RadialSubdivisionReportingRowPayload
from tools.measurement_claim_schema import RadialSubdivisionsCasePayload
from tools.measurement_claim_schema import RadialSubdivisionsConfigPayload
from tools.measurement_claim_schema import RadialSweepReconstructionPayload
from tools.measurement_claim_schema import SamplesPerRadian
from tools.measurement_claim_schema import ToolDiameters
from tools.measurement_claim_schema import WorldMillimetres
from tools.measurement_claim_schema import WorldPointMillimetres
from tools.measurement_claim_schema import WorldRectangleMillimetres
from tools.measurement_claim_schema import WorldXYMillimetres

RadialCase = Literal["radial-station", "radial-subdivisions", "radial-floor", "radial-margin"]
RADIAL_CASE_ORDER: Tuple[RadialCase, ...] = ("radial-station", "radial-subdivisions", "radial-floor", "radial-margin")

_TOOL_DIAMETER = Millimetres(2.0)
_GUIDE_STEP = ToolDiameters(0.025)
_MAX_ADVANCE = ToolDiameters(1.0)
_RADIAL_CLEARANCE = Millimetres(0.002)
_CUT_Z = WorldMillimetres(0.0)
_CLEARANCE_Z = WorldMillimetres(2.0)
_MAX_PASSES = 1000
_SAMPLES_PER_RADIAN = SamplesPerRadian(10.0)
_GENERATOR_PROBE_ANGLES: Tuple[Degrees, ...] = tuple(Degrees(360.0 * index / 32) for index in range(32))
_AUDIT_OFFSETS: Tuple[Degrees, ...] = tuple(Degrees(360.0 * index / 16) for index in range(16))
_SUBDIVISIONS = (1, 2, 4, 8, 16)
_FLOORS: Tuple[CoarseSteps, ...] = (CoarseSteps(0.25), CoarseSteps(0.5), CoarseSteps(1.0))
_MARGINS: Tuple[DimensionlessRatio, ...] = tuple(DimensionlessRatio(value) for value in (1.25, 1.4, 1.5, 1.75, 2.0))
_TARGET_CENTRE: WorldXYMillimetres = (WorldMillimetres(18.482), WorldMillimetres(10.482))
_TARGET_MAXIMAL_RADIUS = Millimetres(0.5156)
_TARGET_CENTRE_DECIMALS = 3
_TARGET_RADIUS_DECIMALS = 4

_RadiusSearch = Callable[
    [
        Stock,
        engagement_radial_toolpath._GuideStation,
        float,
        Tuple[float, float],
        engagement_radial_toolpath._Regulation,
    ],
    Optional[engagement_radial_toolpath._RadiusChoice],
]


class _CircleAudit(TypedDict):
    centre: WorldXYMillimetres
    maximal_radius: Millimetres
    selected_radius: Millimetres
    reporting: List[Degrees]
    verdicts: List[bool]
    cuts_material: bool


class _SelectedCircleObservation(_CircleAudit):
    forced: bool
    finishes: bool


class _StationObservation(TypedDict):
    centre: WorldXYMillimetres
    maximal_radius: Millimetres
    coarse_step: Millimetres
    rung_6: _CircleAudit
    rung_7: _CircleAudit
    refined: List[_CircleAudit]
    selected: _SelectedCircleObservation


class _CapturedRun(TypedDict):
    result: ToolpathResult
    all_circles: List[_SelectedCircleObservation]
    non_entry_circles: List[_SelectedCircleObservation]
    stations: List[_StationObservation]


def _polygon_payload(width: Millimetres, height: Millimetres) -> WorldRectangleMillimetres:
    zero = WorldMillimetres(0.0)
    world_width = WorldMillimetres(float(width))
    world_height = WorldMillimetres(float(height))
    lower_left: WorldPointMillimetres = (zero, zero, zero)
    lower_right: WorldPointMillimetres = (world_width, zero, zero)
    upper_right: WorldPointMillimetres = (world_width, world_height, zero)
    upper_left: WorldPointMillimetres = (zero, world_height, zero)
    return lower_left, lower_right, upper_right, upper_left


def _public_call(width: Millimetres, height: Millimetres, cap: Degrees) -> GeneratorPublicCallPayload:
    return {
        "polygon": _polygon_payload(width, height),
        "holes": [],
        "tool_diameter": _TOOL_DIAMETER,
        "tea_cap": cap,
        "guide_step": _GUIDE_STEP,
        "max_advance": _MAX_ADVANCE,
        "radial_clearance": _RADIAL_CLEARANCE,
        "climb": True,
        "cut_z": _CUT_Z,
        "clearance_z": _CLEARANCE_Z,
        "max_passes": _MAX_PASSES,
        "samples_per_radian": _SAMPLES_PER_RADIAN,
    }


def _run_generator(width: Millimetres, height: Millimetres, cap: Degrees) -> ToolpathResult:
    return radius_regulated_toolpath(
        polygon=Polygon(_polygon_payload(width, height)),
        holes=[],
        tool_diameter=float(_TOOL_DIAMETER),
        tea_cap_deg=float(cap),
        guide_step_tool_diameters=float(_GUIDE_STEP),
        max_advance_tool_diameters=float(_MAX_ADVANCE),
        radial_clearance=float(_RADIAL_CLEARANCE),
        climb=True,
        cut_z=float(_CUT_Z),
        clearance_z=float(_CLEARANCE_Z),
        max_passes=_MAX_PASSES,
        samples_per_radian=float(_SAMPLES_PER_RADIAN),
    )


def _provenance() -> RadialSelectionDecisionProvenancePayload:
    return {
        "policy": "radial-known-reporting-driven-selection/v1",
        "native_cap_decision_site": "compas_cgal._stock_2.engagement_at[2]",
        "engagement_reporting_value_site": "compas_cgal._stock_2.engagement_at[1]",
        "reporting_driven_decision_sites": [
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
        ],
    }


def _audit_circle(
    stock: Stock,
    station: engagement_radial_toolpath._GuideStation,
    maximal_radius: float,
    advance: Tuple[float, float],
    regulation: engagement_radial_toolpath._Regulation,
) -> _CircleAudit:
    entry_x, entry_y = station.entry
    radial_x = (entry_x - station.cx) / station.radius
    radial_y = (entry_y - station.cy) / station.radius
    reporting: List[Degrees] = []
    verdicts: List[bool] = []
    for offset in _AUDIT_OFFSETS:
        angle = math.radians(float(offset))
        cos_angle = math.cos(angle)
        sin_angle = math.sin(angle)
        direction_x = radial_x * cos_angle - radial_y * sin_angle
        direction_y = radial_x * sin_angle + radial_y * cos_angle
        x = station.cx + station.radius * direction_x
        y = station.cy + station.radius * direction_y
        observation = _stock_2.engagement_at(stock.raw, x, y, regulation.tool_radius, regulation.cap_ratio, 0.0)
        reported = observation[1]
        exceeded = observation[2]
        if type(reported) is not float or not math.isfinite(reported):
            raise ProbeInstrumentationContractError("radial audit reporting value must be one finite native float at engagement_at index 1")
        if type(exceeded) is not bool:
            raise ProbeInstrumentationContractError("radial audit cap verdict must be one native bool at engagement_at index 2")
        reporting.append(Degrees(math.degrees(reported)))
        verdicts.append(exceeded)
    return {
        "centre": (WorldMillimetres(station.cx), WorldMillimetres(station.cy)),
        "maximal_radius": Millimetres(maximal_radius),
        "selected_radius": Millimetres(station.radius),
        "reporting": reporting,
        "verdicts": verdicts,
        "cuts_material": engagement_radial_toolpath._loop_reaches_material(stock, station, advance, regulation.tool_radius),
    }


def _is_target(station: engagement_radial_toolpath._GuideStation, emitted_radius: float) -> bool:
    return (
        emitted_radius == 0.0
        and round(station.cx, _TARGET_CENTRE_DECIMALS) == float(_TARGET_CENTRE[0])
        and round(station.cy, _TARGET_CENTRE_DECIMALS) == float(_TARGET_CENTRE[1])
        and round(station.radius, _TARGET_RADIUS_DECIMALS) == float(_TARGET_MAXIMAL_RADIUS)
    )


def _station_observation(
    stock: Stock,
    station: engagement_radial_toolpath._GuideStation,
    emitted_radius: float,
    advance: Tuple[float, float],
    regulation: engagement_radial_toolpath._Regulation,
    selected: _SelectedCircleObservation,
) -> _StationObservation:
    coarse = engagement_radial_toolpath._radius_ladder(
        station.radius,
        emitted_radius,
        regulation.guide_step,
        rungs=engagement_radial_toolpath.RADIUS_LADDER_RUNGS,
        floor_radius=engagement_radial_toolpath.NO_RADIUS_FLOOR,
    )
    if len(coarse) <= 7:
        raise ProbeInstrumentationContractError("target station does not expose coarse rungs 6 and 7")
    coarse_observations = [_audit_circle(stock, replace(station, radius=coarse[index]), station.radius, advance, regulation) for index in (6, 7)]
    refined_radii = engagement_radial_toolpath._refined_radius_ladder(station.radius, emitted_radius, regulation.guide_step)
    refined = [_audit_circle(stock, replace(station, radius=radius), station.radius, advance, regulation) for radius in refined_radii]
    return {
        "centre": (WorldMillimetres(station.cx), WorldMillimetres(station.cy)),
        "maximal_radius": Millimetres(station.radius),
        "coarse_step": Millimetres(regulation.guide_step),
        "rung_6": coarse_observations[0],
        "rung_7": coarse_observations[1],
        "refined": refined,
        "selected": selected,
    }


def _capture_run(
    width: Millimetres,
    height: Millimetres,
    cap: Degrees,
    original_search: _RadiusSearch,
) -> _CapturedRun:
    observations: List[_SelectedCircleObservation] = []
    stations: List[_StationObservation] = []

    def observe_after_search(
        stock: Stock,
        station: engagement_radial_toolpath._GuideStation,
        emitted_radius: float,
        advance: Tuple[float, float],
        regulation: engagement_radial_toolpath._Regulation,
    ) -> Optional[engagement_radial_toolpath._RadiusChoice]:
        choice = original_search(stock, station, emitted_radius, advance, regulation)
        if choice is None:
            return None
        candidate = replace(station, radius=choice.radius)
        audit = _audit_circle(stock, candidate, station.radius, advance, regulation)
        observed: _SelectedCircleObservation = {
            "centre": audit["centre"],
            "maximal_radius": audit["maximal_radius"],
            "selected_radius": audit["selected_radius"],
            "reporting": audit["reporting"],
            "verdicts": audit["verdicts"],
            "cuts_material": audit["cuts_material"],
            "forced": choice.forced,
            "finishes": choice.finishes,
        }
        observations.append(observed)
        if _is_target(station, emitted_radius):
            stations.append(_station_observation(stock, station, emitted_radius, advance, regulation, observed))
        return choice

    engagement_radial_toolpath._largest_admissible_radius = observe_after_search
    result = _run_generator(width, height, cap)
    circles = [operation for operation in result.operations if operation.operation is OperationType.CUT and isinstance(operation.geometry, Circle)]
    if len(circles) != len(observations):
        raise ProbeInstrumentationContractError("radial selection observations do not align one-for-one with emitted CUT circles")
    seen_paths: set[int] = set()
    non_entry: List[_SelectedCircleObservation] = []
    for operation, observation in zip(circles, observations):
        if operation.path_index in seen_paths:
            non_entry.append(observation)
        else:
            seen_paths.add(operation.path_index)
    return {"result": result, "all_circles": observations, "non_entry_circles": non_entry, "stations": stations}


def _counts(circles: List[_SelectedCircleObservation]) -> NativeObservationCountsPayload:
    verdicts = [verdict for circle in circles for verdict in circle["verdicts"]]
    exceeded = sum(verdicts)
    return {"observations": len(verdicts), "accepted": len(verdicts) - exceeded, "exceeded": exceeded}


def _worst_peak(circles: List[_SelectedCircleObservation]) -> Degrees:
    if not circles:
        raise ProbeInstrumentationContractError("radial audit contains no non-entry machining circle")
    return Degrees(max(float(value) for circle in circles for value in circle["reporting"]))


def _circles_over_cap(circles: List[_SelectedCircleObservation]) -> int:
    return sum(any(circle["verdicts"]) for circle in circles)


def _cutting_length(result: ToolpathResult) -> Millimetres:
    total = 0.0
    for operation in result.operations:
        if operation.operation is not OperationType.CUT:
            continue
        geometry = operation.geometry
        if isinstance(geometry, Circle):
            total += 2.0 * math.pi * float(geometry.radius)
        elif isinstance(geometry, Line):
            total += math.hypot(float(geometry.end.x - geometry.start.x), float(geometry.end.y - geometry.start.y))
        elif isinstance(geometry, Arc):
            raise ProbeInstrumentationContractError("radial generator emitted an unexpected CUT arc")
        else:
            raise ProbeInstrumentationContractError("radial generator emitted an unknown CUT geometry")
    return Millimetres(float(round(total)))


def _audit_config() -> RadialAuditConfigPayload:
    return {
        "phase": "entry-angle",
        "probe_offsets": list(_AUDIT_OFFSETS),
        "includes_entry_phase": True,
        "adds_separate_entry_probe": False,
        "excludes_chain_entry_circles": True,
    }


def _sweep_reconstruction(circles_by_configuration: Sequence[List[_SelectedCircleObservation]]) -> RadialSweepReconstructionPayload:
    return {
        "stock_model": "generator-faithful-radial-pre-bridge/v1",
        "audit_position_count": len(_AUDIT_OFFSETS),
        "audit_includes_entry_phase": True,
        "audit_adds_separate_entry_probe": False,
        "excludes_chain_entry_circles": True,
        "non_entry_circle_counts": [len(circles) for circles in circles_by_configuration],
    }


def _station_payload(original_search: _RadiusSearch) -> RadialStationCasePayload:
    run = _capture_run(Millimetres(20.0), Millimetres(12.0), Degrees(60.0), original_search)
    if len(run["stations"]) != 1:
        raise ProbeInstrumentationContractError(f"radial-station requires one target stock occurrence, observed {len(run['stations'])}")
    station = run["stations"][0]
    refined_radii = [Millimetres(round(float(row["selected_radius"]), _TARGET_RADIUS_DECIMALS)) for row in station["refined"]]
    refined_native: List[RadialStationRefinedCandidateNativePayload] = [
        {"radius": radius, "cap_exceeded": any(row["verdicts"]), "cuts_material": row["cuts_material"]} for radius, row in zip(refined_radii, station["refined"])
    ]
    compliant = [row for row in station["refined"] if not any(row["verdicts"]) and row["cuts_material"]]
    if not compliant:
        raise ProbeInstrumentationContractError("target station has no compliant material-cutting refined band")
    rung_6 = station["rung_6"]
    rung_7 = station["rung_7"]
    selected = station["selected"]
    maximal = station["refined"][0]
    if maximal["selected_radius"] != station["maximal_radius"] or not any(maximal["verdicts"]) or not maximal["cuts_material"]:
        raise ProbeInstrumentationContractError("target station maximal refined candidate is not the over-cap material-cutting forced counterfactual")
    selected_candidates = [row for row in compliant if row["selected_radius"] == selected["selected_radius"]]
    if (
        selected["forced"]
        or selected["finishes"]
        or any(selected["verdicts"])
        or not selected["cuts_material"]
        or len(selected_candidates) != 1
        or selected_candidates[0]["reporting"] != selected["reporting"]
        or selected_candidates[0]["verdicts"] != selected["verdicts"]
    ):
        raise ProbeInstrumentationContractError("target station selected circle is not a non-forced, non-finishing material-cutting refined rescue")
    reconstruction: RadialStationReconstructionPayload = {
        "stock_model": "generator-faithful-radial-pre-bridge/v1",
        "target_centre": _TARGET_CENTRE,
        "centre_decimal_places": _TARGET_CENTRE_DECIMALS,
        "target_maximal_radius": _TARGET_MAXIMAL_RADIUS,
        "radius_decimal_places": _TARGET_RADIUS_DECIMALS,
        "occurrence_count": len(run["stations"]),
        "coarse_step": Millimetres(round(float(station["coarse_step"]), 4)),
        "coarse_rungs": [6, 7],
        "refined_radius_sequence": refined_radii,
    }
    native: RadialStationNativeDecisionsPayload = {
        "rung_6_cap_exceeded": any(rung_6["verdicts"]),
        "rung_6_cuts_material": rung_6["cuts_material"],
        "rung_7_cap_exceeded": any(rung_7["verdicts"]),
        "rung_7_cuts_material": rung_7["cuts_material"],
        "refined_candidates": refined_native,
    }
    reporting: RadialStationReportingPayload = {
        "maximal_radius": Millimetres(round(float(station["maximal_radius"]), 4)),
        "rung_6_radius": Millimetres(round(float(rung_6["selected_radius"]), 4)),
        "rung_6_peak": Degrees(round(max(float(value) for value in rung_6["reporting"]), 1)),
        "rung_7_radius": Millimetres(round(float(rung_7["selected_radius"]), 4)),
        "rung_7_peak": Degrees(round(max(float(value) for value in rung_7["reporting"]), 1)),
        "refined_band_min_radius": Millimetres(round(min(float(row["selected_radius"]) for row in compliant), 4)),
        "refined_band_max_radius": Millimetres(round(max(float(row["selected_radius"]) for row in compliant), 4)),
        "forced_peak": Degrees(round(max(float(value) for value in maximal["reporting"]), 1)),
        "rescued_peak": Degrees(round(max(float(value) for value in selected["reporting"]), 1)),
        "angle_unit": "degree",
        "length_unit": "mm",
    }
    config: RadialStationConfigPayload = {
        "public_call": _public_call(Millimetres(20.0), Millimetres(12.0), Degrees(60.0)),
        "subdivisions": 8,
        "floor_steps": CoarseSteps(0.5),
        "refinement_margin": DimensionlessRatio(1.4),
        "generator_probe_angles": list(_GENERATOR_PROBE_ANGLES),
        "ladder_step": ToolDiameters(0.025),
        "ladder_span": ToolDiameters(1.0),
        "ladder_rungs": 40,
        "max_radial_sweeps": 40,
    }
    return {
        "case": "radial-station",
        "source_claim_ids": ["MC-001", "MC-002", "MC-003"],
        "config": config,
        "reconstruction": reconstruction,
        "native_sampled_decisions": native,
        "reporting_values": reporting,
        "selection_decision_provenance": _provenance(),
        "continuous_certificate": None,
    }


def _subdivisions_payload(original_search: _RadiusSearch) -> RadialSubdivisionsCasePayload:
    rows: List[Tuple[int, _CapturedRun]] = []
    for value in _SUBDIVISIONS:
        engagement_radial_toolpath.RADIUS_LADDER_SUBDIVISIONS = value
        rows.append((value, _capture_run(Millimetres(20.0), Millimetres(12.0), Degrees(60.0), original_search)))
    native: List[RadialSubdivisionNativeRowPayload] = []
    reporting: List[RadialSubdivisionReportingRowPayload] = []
    for value, run in rows:
        circles = run["non_entry_circles"]
        native.append({"subdivisions": value, "selected_circles": len(circles), "observations": _counts(circles)})
        reporting.append(
            {"subdivisions": value, "worst_peak": Degrees(round(float(_worst_peak(circles)), 1)), "circles_over_cap": _circles_over_cap(circles), "angle_unit": "degree"}
        )
    config: RadialSubdivisionsConfigPayload = {
        "public_call": _public_call(Millimetres(20.0), Millimetres(12.0), Degrees(60.0)),
        "subdivision_values": list(_SUBDIVISIONS),
        "floor_steps": CoarseSteps(0.5),
        "refinement_margin": DimensionlessRatio(1.4),
        "generator_probe_angles": list(_GENERATOR_PROBE_ANGLES),
        "audit": _audit_config(),
    }
    return {
        "case": "radial-subdivisions",
        "source_claim_ids": ["MC-004", "MC-005", "MC-006"],
        "config": config,
        "reconstruction": _sweep_reconstruction([run["non_entry_circles"] for _, run in rows]),
        "native_sampled_decisions": native,
        "reporting_values": reporting,
        "selection_decision_provenance": _provenance(),
        "continuous_certificate": None,
    }


def _floor_payload(original_search: _RadiusSearch) -> RadialFloorCasePayload:
    rows: List[Tuple[CoarseSteps, _CapturedRun]] = []
    for value in _FLOORS:
        engagement_radial_toolpath.RADIUS_LADDER_FLOOR_STEPS = float(value)
        rows.append((value, _capture_run(Millimetres(20.0), Millimetres(12.0), Degrees(60.0), original_search)))
    native: List[RadialFloorNativeRowPayload] = []
    reporting: List[RadialFloorReportingRowPayload] = []
    for value, run in rows:
        circles = run["non_entry_circles"]
        native.append({"floor_steps": value, "selected_circles": len(circles), "observations": _counts(circles)})
        reporting.append(
            {
                "floor_steps": value,
                "worst_peak": Degrees(round(float(_worst_peak(circles)), 1)),
                "circles_over_cap": _circles_over_cap(circles),
                "angle_unit": "degree",
            }
        )
    config: RadialFloorConfigPayload = {
        "public_call": _public_call(Millimetres(20.0), Millimetres(12.0), Degrees(60.0)),
        "floor_values": list(_FLOORS),
        "subdivisions": 8,
        "refinement_margin": DimensionlessRatio(1.4),
        "generator_probe_angles": list(_GENERATOR_PROBE_ANGLES),
        "audit": _audit_config(),
    }
    return {
        "case": "radial-floor",
        "source_claim_ids": ["MC-007"],
        "config": config,
        "reconstruction": _sweep_reconstruction([run["non_entry_circles"] for _, run in rows]),
        "native_sampled_decisions": native,
        "reporting_values": reporting,
        "selection_decision_provenance": _provenance(),
        "continuous_certificate": None,
    }


def _margin_payload(original_search: _RadiusSearch) -> RadialMarginCasePayload:
    rows: List[Tuple[DimensionlessRatio, _CapturedRun]] = []
    for value in _MARGINS:
        engagement_radial_toolpath.RADIUS_LADDER_REFINEMENT_MARGIN = float(value)
        rows.append((value, _capture_run(Millimetres(6.0), Millimetres(4.0), Degrees(40.0), original_search)))
    native: List[RadialMarginNativeRowPayload] = []
    reporting: List[RadialMarginReportingRowPayload] = []
    for value, run in rows:
        circles = run["non_entry_circles"]
        native.append({"refinement_margin": value, "selected_circles": len(circles), "observations": _counts(circles)})
        reporting.append(
            {
                "refinement_margin": value,
                "worst_peak": Degrees(round(float(_worst_peak(circles)), 1)),
                "circles_over_cap": _circles_over_cap(circles),
                "cutting_length": _cutting_length(run["result"]),
                "angle_unit": "degree",
                "length_unit": "mm",
            }
        )
    config: RadialMarginConfigPayload = {
        "public_call": _public_call(Millimetres(6.0), Millimetres(4.0), Degrees(40.0)),
        "margin_values": list(_MARGINS),
        "subdivisions": 8,
        "floor_steps": CoarseSteps(0.5),
        "generator_probe_angles": list(_GENERATOR_PROBE_ANGLES),
        "audit": _audit_config(),
    }
    return {
        "case": "radial-margin",
        "source_claim_ids": ["MC-008"],
        "config": config,
        "reconstruction": _sweep_reconstruction([run["non_entry_circles"] for _, run in rows]),
        "native_sampled_decisions": native,
        "reporting_values": reporting,
        "selection_decision_provenance": _provenance(),
        "continuous_certificate": None,
    }


def run_radial_case(case: RadialCase) -> GeneratorCasePayload:
    """Run one fixed radial claim case and restore every owned mutation seam."""
    if case not in RADIAL_CASE_ORDER:
        raise UnknownMeasurementClaimCaseError(f"unknown radial measurement-claim case: {case!r}")
    if engagement_toolpath.LOOP_PROBE_COUNT != len(_GENERATOR_PROBE_ANGLES) or engagement_toolpath.LOOP_PROBE_ANGLES_DEG != tuple(
        float(angle) for angle in _GENERATOR_PROBE_ANGLES
    ):
        raise InvalidMeasurementClaimConfigError("radial generator probe ring does not equal the fixed 32-position measurement configuration")
    if engagement_radial_toolpath.RADIUS_LADDER_RUNGS != 40 or engagement_radial_toolpath.MAX_RADIAL_SWEEPS_PER_CHAIN != 40:
        raise InvalidMeasurementClaimConfigError("radial ladder rung and sweep budgets do not equal the fixed measurement configuration")
    original_subdivisions = engagement_radial_toolpath.RADIUS_LADDER_SUBDIVISIONS
    original_floor = engagement_radial_toolpath.RADIUS_LADDER_FLOOR_STEPS
    original_margin = engagement_radial_toolpath.RADIUS_LADDER_REFINEMENT_MARGIN
    original_search = engagement_radial_toolpath._largest_admissible_radius
    try:
        engagement_radial_toolpath.RADIUS_LADDER_SUBDIVISIONS = 8
        engagement_radial_toolpath.RADIUS_LADDER_FLOOR_STEPS = 0.5
        engagement_radial_toolpath.RADIUS_LADDER_REFINEMENT_MARGIN = 1.4
        if case == "radial-station":
            return _station_payload(original_search)
        if case == "radial-subdivisions":
            return _subdivisions_payload(original_search)
        if case == "radial-floor":
            return _floor_payload(original_search)
        return _margin_payload(original_search)
    finally:
        engagement_radial_toolpath.RADIUS_LADDER_SUBDIVISIONS = original_subdivisions
        engagement_radial_toolpath.RADIUS_LADDER_FLOOR_STEPS = original_floor
        engagement_radial_toolpath.RADIUS_LADDER_REFINEMENT_MARGIN = original_margin
        engagement_radial_toolpath._largest_admissible_radius = original_search
