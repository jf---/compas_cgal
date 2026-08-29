"""Typed validation for authenticated generator-claim artifacts."""

from __future__ import annotations

import datetime
import math
import pathlib
import re
from typing import Dict
from typing import List
from typing import Literal
from typing import NewType
from typing import Optional
from typing import Sequence
from typing import Tuple
from typing import TypedDict
from typing import Union
from typing import cast

from tools import measurement_artifact
from tools import measurement_claim_json
from tools.measurement_artifact import GitObjectId
from tools.measurement_artifact import MeasurementArtifactError
from tools.measurement_artifact import ValidatedEnvelope
from tools.measurement_claim_adjudication import adjudicate_station_claims
from tools.measurement_claim_units import Degrees as Degrees
from tools.measurement_claim_units import Millimetres as Millimetres
from tools.measurement_claim_units import WorldMillimetres as WorldMillimetres

ToolDiameters = NewType("ToolDiameters", float)
CoarseSteps = NewType("CoarseSteps", float)
DimensionlessRatio = NewType("DimensionlessRatio", float)
SamplesPerRadian = NewType("SamplesPerRadian", float)
ValidatedArtifactStartedUtc = NewType("ValidatedArtifactStartedUtc", datetime.datetime)
ValidatedArtifactDirectory = NewType("ValidatedArtifactDirectory", pathlib.PurePosixPath)
AngleUnit = Literal["degree"]
LengthUnit = Literal["mm"]
SpacingUnit = Literal["tool-diameter"]
WorldXYMillimetres = Tuple[WorldMillimetres, WorldMillimetres]
WorldPointMillimetres = Tuple[WorldMillimetres, WorldMillimetres, WorldMillimetres]
WorldRectangleMillimetres = Tuple[WorldPointMillimetres, WorldPointMillimetres, WorldPointMillimetres, WorldPointMillimetres]

ClaimDisposition = Literal["re-earned", "corrected", "historical", "deleted"]
RadialCase = Literal["radial-station", "radial-subdivisions", "radial-floor", "radial-margin"]
AdvanceCase = Literal["advance-placement", "advance-probe-count"]
GeneratorCase = Union[RadialCase, AdvanceCase]
ProbeConfigurationLabel = Literal["3-old", "8", "12", "16", "24", "32", "40", "48"]
RadialStockModel = Literal["generator-faithful-radial-pre-bridge/v1"]
AdvanceStockModel = Literal["generator-pre-bridge/v1"]
SelectedCirclePolicy = Literal["accepted-non-forced/v1"]
RadialAuditPhase = Literal["entry-angle"]
AdvancePlacementAuditPhase = Literal["advance-direction"]
AdvanceProbeCountAuditPhase = Literal["advance-half-step"]
MissingConfiguration = Literal["refinement-without-reporting-ranking"]


class LeastBadRungDecisionSitePayload(TypedDict):
    symbol: Literal["compas_cgal.engagement_radial_toolpath._least_bad_rung"]
    value_source: Literal["compas_cgal._stock_2.engagement_at[1]"]
    effect: Literal["forced-radius-selection"]


class LargestAdmissibleRadiusDecisionSitePayload(TypedDict):
    symbol: Literal["compas_cgal.engagement_radial_toolpath._largest_admissible_radius"]
    value_source: Literal["compas_cgal.engagement_radial_toolpath._GentlestRung.peak"]
    effect: Literal["refined-scan-control"]


SelectionDecisionSitePayload = Union[LeastBadRungDecisionSitePayload, LargestAdmissibleRadiusDecisionSitePayload]


class RadialSelectionDecisionProvenancePayload(TypedDict):
    policy: Literal["radial-known-reporting-driven-selection/v1"]
    native_cap_decision_site: Literal["compas_cgal._stock_2.engagement_at[2]"]
    engagement_reporting_value_site: Literal["compas_cgal._stock_2.engagement_at[1]"]
    reporting_driven_decision_sites: List[SelectionDecisionSitePayload]


class AdvanceSelectionDecisionProvenancePayload(TypedDict):
    policy: Literal["advance-native-cap-reporting-observation/v1"]
    native_cap_decision_site: Literal["compas_cgal._stock_2.engagement_at[2]"]
    engagement_reporting_value_site: Literal["compas_cgal._stock_2.engagement_at[1]"]
    reporting_driven_decision_sites: List[SelectionDecisionSitePayload]


SelectionDecisionProvenancePayload = Union[RadialSelectionDecisionProvenancePayload, AdvanceSelectionDecisionProvenancePayload]


class GeneratorPublicCallPayload(TypedDict):
    polygon: WorldRectangleMillimetres
    holes: List[List[WorldPointMillimetres]]
    tool_diameter: Millimetres
    tea_cap: Degrees
    guide_step: ToolDiameters
    max_advance: ToolDiameters
    radial_clearance: Millimetres
    climb: bool
    cut_z: WorldMillimetres
    clearance_z: WorldMillimetres
    max_passes: int
    samples_per_radian: SamplesPerRadian


class RadialAuditConfigPayload(TypedDict):
    phase: RadialAuditPhase
    probe_offsets: List[Degrees]
    includes_entry_phase: Literal[True]
    adds_separate_entry_probe: Literal[False]
    excludes_chain_entry_circles: Literal[True]


class AdvancePlacementAuditConfigPayload(TypedDict):
    phase: AdvancePlacementAuditPhase
    probe_offsets: List[Degrees]
    excludes_entry_probe: Literal[True]


class AdvanceProbeCountAuditConfigPayload(TypedDict):
    phase: AdvanceProbeCountAuditPhase
    probe_offsets: List[Degrees]
    excludes_entry_probe: Literal[True]


class AdvanceProbeConfigurationPayload(TypedDict):
    label: ProbeConfigurationLabel
    probe_count: int
    probe_angles: List[Degrees]


class RadialStationConfigPayload(TypedDict):
    public_call: GeneratorPublicCallPayload
    subdivisions: int
    floor_steps: CoarseSteps
    refinement_margin: DimensionlessRatio
    generator_probe_angles: List[Degrees]
    ladder_step: ToolDiameters
    ladder_span: ToolDiameters
    ladder_rungs: int
    max_radial_sweeps: int


class RadialSubdivisionsConfigPayload(TypedDict):
    public_call: GeneratorPublicCallPayload
    subdivision_values: List[int]
    floor_steps: CoarseSteps
    refinement_margin: DimensionlessRatio
    generator_probe_angles: List[Degrees]
    audit: RadialAuditConfigPayload


class RadialFloorConfigPayload(TypedDict):
    public_call: GeneratorPublicCallPayload
    floor_values: List[CoarseSteps]
    subdivisions: int
    refinement_margin: DimensionlessRatio
    generator_probe_angles: List[Degrees]
    audit: RadialAuditConfigPayload


class RadialMarginConfigPayload(TypedDict):
    public_call: GeneratorPublicCallPayload
    margin_values: List[DimensionlessRatio]
    subdivisions: int
    floor_steps: CoarseSteps
    generator_probe_angles: List[Degrees]
    audit: RadialAuditConfigPayload


class AdvancePlacementConfigPayload(TypedDict):
    public_calls: List[GeneratorPublicCallPayload]
    generator_probe_angles: List[Degrees]
    generator_prepends_entry_probe: Literal[True]
    audit: AdvancePlacementAuditConfigPayload


class AdvanceProbeCountConfigPayload(TypedDict):
    public_calls: List[GeneratorPublicCallPayload]
    generator_configurations: List[AdvanceProbeConfigurationPayload]
    generator_prepends_entry_probe: Literal[True]
    audit: AdvanceProbeCountAuditConfigPayload


class RadialStationReconstructionPayload(TypedDict):
    stock_model: RadialStockModel
    target_centre: WorldXYMillimetres
    centre_decimal_places: int
    target_maximal_radius: Millimetres
    radius_decimal_places: int
    occurrence_count: int
    coarse_step: Millimetres
    coarse_rungs: List[int]
    refined_radius_sequence: List[Millimetres]


class RadialSweepReconstructionPayload(TypedDict):
    stock_model: RadialStockModel
    audit_position_count: int
    audit_includes_entry_phase: Literal[True]
    audit_adds_separate_entry_probe: Literal[False]
    excludes_chain_entry_circles: Literal[True]
    non_entry_circle_counts: List[int]


class AdvanceReconstructionPayload(TypedDict):
    stock_model: AdvanceStockModel
    selected_circle_policy: SelectedCirclePolicy
    original_calls_per_wrapper: Literal[1]
    generator_prepends_entry_probe: Literal[True]
    audit_excludes_entry_probe: Literal[True]
    selected_circle_count: int


class NativeObservationCountsPayload(TypedDict):
    observations: int
    accepted: int
    exceeded: int


class RadialStationRefinedCandidateNativePayload(TypedDict):
    radius: Millimetres
    cap_exceeded: bool
    cuts_material: bool


class RadialStationNativeDecisionsPayload(TypedDict):
    rung_6_cap_exceeded: Optional[bool]
    rung_6_cuts_material: Optional[bool]
    rung_7_cap_exceeded: Optional[bool]
    rung_7_cuts_material: Optional[bool]
    refined_candidates: List[RadialStationRefinedCandidateNativePayload]


class RadialSubdivisionNativeRowPayload(TypedDict):
    subdivisions: int
    selected_circles: int
    observations: NativeObservationCountsPayload


class RadialFloorNativeRowPayload(TypedDict):
    floor_steps: CoarseSteps
    selected_circles: int
    observations: NativeObservationCountsPayload


class RadialMarginNativeRowPayload(TypedDict):
    refinement_margin: DimensionlessRatio
    selected_circles: int
    observations: NativeObservationCountsPayload


class AdvancePlacementNativeRowPayload(TypedDict):
    cap: Degrees
    climb: bool
    selected_circles: int
    observations: NativeObservationCountsPayload


class AdvanceProbeCountNativeRowPayload(TypedDict):
    label: ProbeConfigurationLabel
    probe_count: int
    cap: Degrees
    selected_circles: int
    observations: NativeObservationCountsPayload


class AdvanceOffsetBinCountPayload(TypedDict):
    offset: Degrees
    positions_over_cap: int
    angle_unit: AngleUnit


class RadialStationReportingPayload(TypedDict):
    maximal_radius: Optional[Millimetres]
    rung_6_radius: Optional[Millimetres]
    rung_6_peak: Optional[Degrees]
    rung_7_radius: Optional[Millimetres]
    rung_7_peak: Optional[Degrees]
    refined_band_min_radius: Optional[Millimetres]
    refined_band_max_radius: Optional[Millimetres]
    forced_peak: Optional[Degrees]
    rescued_peak: Optional[Degrees]
    angle_unit: AngleUnit
    length_unit: LengthUnit


class RadialSubdivisionReportingRowPayload(TypedDict):
    subdivisions: int
    worst_peak: Degrees
    circles_over_cap: int
    angle_unit: AngleUnit


class RadialFloorReportingRowPayload(TypedDict):
    floor_steps: CoarseSteps
    worst_peak: Degrees
    circles_over_cap: int
    angle_unit: AngleUnit


class RadialMarginReportingRowPayload(TypedDict):
    refinement_margin: DimensionlessRatio
    worst_peak: Degrees
    circles_over_cap: int
    cutting_length: Millimetres
    angle_unit: AngleUnit
    length_unit: LengthUnit


class AdvancePlacementReportingRowPayload(TypedDict):
    cap: Degrees
    climb: bool
    selected_circles: int
    positions_over_cap: int
    worst_peak: Degrees
    worst_peak_offset: Degrees
    old_probe_peak: Degrees
    offset_bin_counts: List[AdvanceOffsetBinCountPayload]
    angle_unit: AngleUnit


class AdvanceProbeCountReportingRowPayload(TypedDict):
    label: ProbeConfigurationLabel
    probe_count: int
    cap: Degrees
    worst_peak: Degrees
    angle_unit: AngleUnit


class RadialStationCasePayload(TypedDict):
    case: Literal["radial-station"]
    source_claim_ids: List[Literal["MC-001", "MC-002", "MC-003"]]
    config: RadialStationConfigPayload
    reconstruction: RadialStationReconstructionPayload
    native_sampled_decisions: RadialStationNativeDecisionsPayload
    reporting_values: RadialStationReportingPayload
    selection_decision_provenance: RadialSelectionDecisionProvenancePayload
    continuous_certificate: None


class RadialSubdivisionsCasePayload(TypedDict):
    case: Literal["radial-subdivisions"]
    source_claim_ids: List[Literal["MC-004", "MC-005", "MC-006"]]
    config: RadialSubdivisionsConfigPayload
    reconstruction: RadialSweepReconstructionPayload
    native_sampled_decisions: List[RadialSubdivisionNativeRowPayload]
    reporting_values: List[RadialSubdivisionReportingRowPayload]
    selection_decision_provenance: RadialSelectionDecisionProvenancePayload
    continuous_certificate: None


class RadialFloorCasePayload(TypedDict):
    case: Literal["radial-floor"]
    source_claim_ids: List[Literal["MC-007"]]
    config: RadialFloorConfigPayload
    reconstruction: RadialSweepReconstructionPayload
    native_sampled_decisions: List[RadialFloorNativeRowPayload]
    reporting_values: List[RadialFloorReportingRowPayload]
    selection_decision_provenance: RadialSelectionDecisionProvenancePayload
    continuous_certificate: None


class RadialMarginCasePayload(TypedDict):
    case: Literal["radial-margin"]
    source_claim_ids: List[Literal["MC-008"]]
    config: RadialMarginConfigPayload
    reconstruction: RadialSweepReconstructionPayload
    native_sampled_decisions: List[RadialMarginNativeRowPayload]
    reporting_values: List[RadialMarginReportingRowPayload]
    selection_decision_provenance: RadialSelectionDecisionProvenancePayload
    continuous_certificate: None


class AdvancePlacementCasePayload(TypedDict):
    case: Literal["advance-placement"]
    source_claim_ids: List[Literal["MC-009"]]
    config: AdvancePlacementConfigPayload
    reconstruction: AdvanceReconstructionPayload
    native_sampled_decisions: List[AdvancePlacementNativeRowPayload]
    reporting_values: List[AdvancePlacementReportingRowPayload]
    selection_decision_provenance: AdvanceSelectionDecisionProvenancePayload
    continuous_certificate: None


class AdvanceProbeCountCasePayload(TypedDict):
    case: Literal["advance-probe-count"]
    source_claim_ids: List[Literal["MC-010"]]
    config: AdvanceProbeCountConfigPayload
    reconstruction: AdvanceReconstructionPayload
    native_sampled_decisions: List[AdvanceProbeCountNativeRowPayload]
    reporting_values: List[AdvanceProbeCountReportingRowPayload]
    selection_decision_provenance: AdvanceSelectionDecisionProvenancePayload
    continuous_certificate: None


GeneratorCasePayload = Union[
    RadialStationCasePayload,
    RadialSubdivisionsCasePayload,
    RadialFloorCasePayload,
    RadialMarginCasePayload,
    AdvancePlacementCasePayload,
    AdvanceProbeCountCasePayload,
]


class MC001EvidencePayload(TypedDict):
    occurrence_count: int
    station_centre: Optional[WorldXYMillimetres]
    maximal_radius: Optional[Millimetres]
    length_unit: LengthUnit


class MC002EvidencePayload(TypedDict):
    coarse_step: Millimetres
    rung_6_radius: Optional[Millimetres]
    rung_6_peak: Optional[Degrees]
    rung_6_cuts_material: Optional[bool]
    rung_7_radius: Optional[Millimetres]
    rung_7_peak: Optional[Degrees]
    rung_7_cuts_material: Optional[bool]
    angle_unit: AngleUnit
    length_unit: LengthUnit


class MC003EvidencePayload(TypedDict):
    rung_7_radius: Optional[Millimetres]
    rung_7_peak: Optional[Degrees]
    rung_7_cuts_material: Optional[bool]
    refined_band_min_radius: Optional[Millimetres]
    refined_band_max_radius: Optional[Millimetres]
    forced_peak: Optional[Degrees]
    rescued_peak: Optional[Degrees]
    angle_unit: AngleUnit
    length_unit: LengthUnit


class MC004EvidencePayload(TypedDict):
    audit_position_count: int
    audit_phase: RadialAuditPhase
    includes_entry_phase: Literal[True]
    adds_separate_entry_probe: Literal[False]


class MC005EvidencePayload(TypedDict):
    non_entry_circle_count: int
    rows: List[RadialSubdivisionReportingRowPayload]


class MC006EvidencePayload(TypedDict):
    history_commit: GitObjectId
    missing_configuration: MissingConfiguration


class MC007EvidencePayload(TypedDict):
    audit_position_count: int
    rows: List[RadialFloorReportingRowPayload]


class MC008EvidencePayload(TypedDict):
    rows: List[RadialMarginReportingRowPayload]
    baseline_available: Literal[False]
    open_ended_gate_claim_removed: Literal[True]
    reporting_selection_disclosed: Literal[True]


class MC009EvidencePayload(TypedDict):
    selected_circle_count: int
    audit_position_count: int
    rows: List[AdvancePlacementReportingRowPayload]


class MC010EvidencePayload(TypedDict):
    audit_position_count: int
    rows: List[AdvanceProbeCountReportingRowPayload]
    timing_claim_removed: Literal[True]
    relative_cost_claim_removed: Literal[True]


class MC001ClaimPayload(TypedDict):
    claim_id: Literal["MC-001"]
    case: Literal["radial-station"]
    disposition: Literal["re-earned", "corrected"]
    reason: str
    selection_decision_provenance: RadialSelectionDecisionProvenancePayload
    evidence: MC001EvidencePayload


class MC002ClaimPayload(TypedDict):
    claim_id: Literal["MC-002"]
    case: Literal["radial-station"]
    disposition: Literal["re-earned", "corrected"]
    reason: str
    selection_decision_provenance: RadialSelectionDecisionProvenancePayload
    evidence: MC002EvidencePayload


class MC003ClaimPayload(TypedDict):
    claim_id: Literal["MC-003"]
    case: Literal["radial-station"]
    disposition: Literal["re-earned", "corrected"]
    reason: str
    selection_decision_provenance: RadialSelectionDecisionProvenancePayload
    evidence: MC003EvidencePayload


class MC004ClaimPayload(TypedDict):
    claim_id: Literal["MC-004"]
    case: Literal["radial-subdivisions"]
    disposition: Literal["corrected", "historical"]
    reason: str
    selection_decision_provenance: RadialSelectionDecisionProvenancePayload
    evidence: MC004EvidencePayload


class MC005ClaimPayload(TypedDict):
    claim_id: Literal["MC-005"]
    case: Literal["radial-subdivisions"]
    disposition: Literal["re-earned", "corrected", "historical"]
    reason: str
    selection_decision_provenance: RadialSelectionDecisionProvenancePayload
    evidence: MC005EvidencePayload


class MC006ClaimPayload(TypedDict):
    claim_id: Literal["MC-006"]
    case: Literal["radial-subdivisions"]
    disposition: Literal["historical", "deleted"]
    reason: str
    selection_decision_provenance: RadialSelectionDecisionProvenancePayload
    evidence: MC006EvidencePayload


class MC007ClaimPayload(TypedDict):
    claim_id: Literal["MC-007"]
    case: Literal["radial-floor"]
    disposition: Literal["re-earned", "corrected", "historical"]
    reason: str
    selection_decision_provenance: RadialSelectionDecisionProvenancePayload
    evidence: MC007EvidencePayload


class MC008ClaimPayload(TypedDict):
    claim_id: Literal["MC-008"]
    case: Literal["radial-margin"]
    disposition: Literal["corrected", "historical"]
    reason: str
    selection_decision_provenance: RadialSelectionDecisionProvenancePayload
    evidence: MC008EvidencePayload


class MC009ClaimPayload(TypedDict):
    claim_id: Literal["MC-009"]
    case: Literal["advance-placement"]
    disposition: Literal["re-earned", "corrected", "historical"]
    reason: str
    selection_decision_provenance: AdvanceSelectionDecisionProvenancePayload
    evidence: MC009EvidencePayload


class MC010ClaimPayload(TypedDict):
    claim_id: Literal["MC-010"]
    case: Literal["advance-probe-count"]
    disposition: Literal["corrected", "historical", "deleted"]
    reason: str
    selection_decision_provenance: AdvanceSelectionDecisionProvenancePayload
    evidence: MC010EvidencePayload


GeneratorClaimRecord = Union[
    MC001ClaimPayload,
    MC002ClaimPayload,
    MC003ClaimPayload,
    MC004ClaimPayload,
    MC005ClaimPayload,
    MC006ClaimPayload,
    MC007ClaimPayload,
    MC008ClaimPayload,
    MC009ClaimPayload,
    MC010ClaimPayload,
]


class RadialStationCaseInputPayload(TypedDict):
    case: Literal["radial-station"]
    source_claim_ids: List[Literal["MC-001", "MC-002", "MC-003"]]
    config: RadialStationConfigPayload
    selection_decision_provenance: RadialSelectionDecisionProvenancePayload


class RadialSubdivisionsCaseInputPayload(TypedDict):
    case: Literal["radial-subdivisions"]
    source_claim_ids: List[Literal["MC-004", "MC-005", "MC-006"]]
    config: RadialSubdivisionsConfigPayload
    selection_decision_provenance: RadialSelectionDecisionProvenancePayload


class RadialFloorCaseInputPayload(TypedDict):
    case: Literal["radial-floor"]
    source_claim_ids: List[Literal["MC-007"]]
    config: RadialFloorConfigPayload
    selection_decision_provenance: RadialSelectionDecisionProvenancePayload


class RadialMarginCaseInputPayload(TypedDict):
    case: Literal["radial-margin"]
    source_claim_ids: List[Literal["MC-008"]]
    config: RadialMarginConfigPayload
    selection_decision_provenance: RadialSelectionDecisionProvenancePayload


class AdvancePlacementCaseInputPayload(TypedDict):
    case: Literal["advance-placement"]
    source_claim_ids: List[Literal["MC-009"]]
    config: AdvancePlacementConfigPayload
    selection_decision_provenance: AdvanceSelectionDecisionProvenancePayload


class AdvanceProbeCountCaseInputPayload(TypedDict):
    case: Literal["advance-probe-count"]
    source_claim_ids: List[Literal["MC-010"]]
    config: AdvanceProbeCountConfigPayload
    selection_decision_provenance: AdvanceSelectionDecisionProvenancePayload


GeneratorCaseInputPayload = Union[
    RadialStationCaseInputPayload,
    RadialSubdivisionsCaseInputPayload,
    RadialFloorCaseInputPayload,
    RadialMarginCaseInputPayload,
    AdvancePlacementCaseInputPayload,
    AdvanceProbeCountCaseInputPayload,
]


class GeneratorClaimInputPayload(TypedDict):
    extraction_commit: GitObjectId
    source_correction_commit: GitObjectId
    case_order: List[GeneratorCase]
    case_inputs: List[GeneratorCaseInputPayload]


class GeneratorClaimPayload(TypedDict):
    schema_version: Literal["measurement-claim-payload/v2"]
    batch: Literal["generator"]
    extraction_commit: GitObjectId
    source_commit: GitObjectId
    source_correction_commit: GitObjectId
    case_order: List[GeneratorCase]
    cases: List[GeneratorCasePayload]
    claims: List[GeneratorClaimRecord]


class MeasurementClaimError(MeasurementArtifactError): ...


class UnknownMeasurementClaimCaseError(MeasurementClaimError): ...


class InvalidMeasurementClaimConfigError(MeasurementClaimError): ...


class ProbeInstrumentationContractError(MeasurementClaimError): ...


class InvalidMeasurementClaimPayloadError(MeasurementClaimError): ...


class MeasurementClaimChildError(MeasurementClaimError): ...


ARTIFACT_KIND = "generator-measurement-claims/v2"
INPUT_VERSION = "generator-measurement-claim-input/v2"
RESULT_VERSION = "generator-measurement-claim-result/v2"
PAYLOAD_NAME = "generator-claims.json"
EXTRACTION_COMMIT = "eec665c1df1cd8d1e98dd9dd1001b5984e17a703"
SOURCE_CORRECTION_COMMIT = "53135e04390e84bf69aa74dc4d0c1ce6ca308eb4"
_HISTORY_COMMIT = "29050b01e656ea7bf577b18f7bb50a04ff9a23c9"
_CASE_ORDER: List[GeneratorCase] = ["radial-station", "radial-subdivisions", "radial-floor", "radial-margin", "advance-placement", "advance-probe-count"]
_CLAIM_IDS = [f"MC-{index:03d}" for index in range(1, 11)]
_CASE_CLAIMS = {
    "radial-station": ["MC-001", "MC-002", "MC-003"],
    "radial-subdivisions": ["MC-004", "MC-005", "MC-006"],
    "radial-floor": ["MC-007"],
    "radial-margin": ["MC-008"],
    "advance-placement": ["MC-009"],
    "advance-probe-count": ["MC-010"],
}
_OBJECT_ID = re.compile(r"(?:[0-9a-f]{40}|[0-9a-f]{64})\Z")
_FINAL_NAME = re.compile(r"\d{4}-\d{2}-\d{2}-[0-9a-f]{12}-generator-[0-9a-f]{12}\Z")
_PLACEMENT_BIN_CENTRES = [angle if angle <= 180.0 else angle - 360.0 for angle in (360.0 * index / 32 for index in range(32))]
_COMMON_KEYS = (
    "case",
    "source_claim_ids",
    "config",
    "reconstruction",
    "native_sampled_decisions",
    "reporting_values",
    "selection_decision_provenance",
    "continuous_certificate",
)


def _fail(field: str, detail: str) -> None:
    raise InvalidMeasurementClaimPayloadError(f"{field}: {detail}")


def _object(value: object, keys: Sequence[str], field: str) -> Dict[str, object]:
    if type(value) is not dict:
        _fail(field, "must be an exact JSON object")
    result = cast(Dict[str, object], value)
    if set(result) != set(keys):
        _fail(field, f"keys are {tuple(result)}, expected {tuple(keys)}")
    return result


def _array(value: object, field: str) -> List[object]:
    if type(value) is not list:
        _fail(field, "must be an exact JSON array")
    return cast(List[object], value)


def _integer(value: object, field: str) -> int:
    if type(value) is not int or value < 0:
        _fail(field, "must be a non-negative exact integer")
    return cast(int, value)


def _number(value: object, field: str, *, optional: bool = False) -> Optional[float]:
    if optional and value is None:
        return None
    if type(value) is not float or not math.isfinite(value):
        _fail(field, "must be a finite exact float")
    return cast(float, value)


def _boolean(value: object, field: str, *, optional: bool = False) -> Optional[bool]:
    if optional and value is None:
        return None
    if type(value) is not bool:
        _fail(field, "must be an exact boolean")
    return cast(bool, value)


def _literal(value: object, expected: object, field: str) -> None:
    if type(value) is not type(expected) or value != expected:
        _fail(field, f"must equal {expected!r}")


def _same(value: object, expected: object, field: str) -> None:
    if type(value) is not type(expected):
        _fail(field, f"type differs from canonical {type(expected).__name__}")
    if type(value) is dict:
        actual_object = cast(Dict[str, object], value)
        expected_object = cast(Dict[str, object], expected)
        if set(actual_object) != set(expected_object):
            _fail(field, "keys differ from canonical configuration")
        for key in expected_object:
            _same(actual_object[key], expected_object[key], f"{field}.{key}")
    elif type(value) is list:
        actual_list = cast(List[object], value)
        expected_list = cast(List[object], expected)
        if len(actual_list) != len(expected_list):
            _fail(field, "length differs from canonical configuration")
        for index, item in enumerate(expected_list):
            _same(actual_list[index], item, f"{field}[{index}]")
    elif value != expected:
        _fail(field, f"must equal {expected!r}")


def _finite_tree(value: object, field: str) -> None:
    measurement_claim_json.validate_finite(value, field, InvalidMeasurementClaimPayloadError)


def _public_call(width: float, height: float, cap: float, climb: bool = True) -> Dict[str, object]:
    return {
        "polygon": [[0.0, 0.0, 0.0], [width, 0.0, 0.0], [width, height, 0.0], [0.0, height, 0.0]],
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


def _radial_audit() -> Dict[str, object]:
    return {
        "phase": "entry-angle",
        "probe_offsets": [360.0 * index / 16 for index in range(16)],
        "includes_entry_phase": True,
        "adds_separate_entry_probe": False,
        "excludes_chain_entry_circles": True,
    }


def _expected_config(case: str) -> Dict[str, object]:
    generator_angles = [360.0 * index / 32 for index in range(32)]
    if case == "radial-station":
        return {
            "public_call": _public_call(20.0, 12.0, 60.0),
            "subdivisions": 8,
            "floor_steps": 0.5,
            "refinement_margin": 1.4,
            "generator_probe_angles": generator_angles,
            "ladder_step": 0.025,
            "ladder_span": 1.0,
            "ladder_rungs": 40,
            "max_radial_sweeps": 40,
        }
    if case == "radial-subdivisions":
        return {
            "public_call": _public_call(20.0, 12.0, 60.0),
            "subdivision_values": [1, 2, 4, 8, 16],
            "floor_steps": 0.5,
            "refinement_margin": 1.4,
            "generator_probe_angles": generator_angles,
            "audit": _radial_audit(),
        }
    if case == "radial-floor":
        return {
            "public_call": _public_call(20.0, 12.0, 60.0),
            "floor_values": [0.25, 0.5, 1.0],
            "subdivisions": 8,
            "refinement_margin": 1.4,
            "generator_probe_angles": generator_angles,
            "audit": _radial_audit(),
        }
    if case == "radial-margin":
        return {
            "public_call": _public_call(6.0, 4.0, 40.0),
            "margin_values": [1.25, 1.4, 1.5, 1.75, 2.0],
            "subdivisions": 8,
            "floor_steps": 0.5,
            "generator_probe_angles": generator_angles,
            "audit": _radial_audit(),
        }
    if case == "advance-placement":
        return {
            "public_calls": [_public_call(20.0, 12.0, cap, climb) for cap in (80.0, 100.0) for climb in (True, False)],
            "generator_probe_angles": [-60.0, 0.0, 60.0],
            "generator_prepends_entry_probe": True,
            "audit": {"phase": "advance-direction", "probe_offsets": generator_angles, "excludes_entry_probe": True},
        }
    if case == "advance-probe-count":
        labels = [("3-old", 3), ("8", 8), ("12", 12), ("16", 16), ("24", 24), ("32", 32), ("40", 40), ("48", 48)]
        return {
            "public_calls": [_public_call(20.0, 12.0, cap) for cap in (40.0, 80.0)],
            "generator_configurations": [
                {
                    "label": label,
                    "probe_count": count,
                    "probe_angles": [-60.0, 0.0, 60.0] if label == "3-old" else [360.0 * index / count for index in range(count)],
                }
                for label, count in labels
            ],
            "generator_prepends_entry_probe": True,
            "audit": {
                "phase": "advance-half-step",
                "probe_offsets": [360.0 * (index + 0.5) / 60 for index in range(60)],
                "excludes_entry_probe": True,
            },
        }
    raise UnknownMeasurementClaimCaseError(case)


def _expected_provenance(radial: bool) -> Dict[str, object]:
    sites: List[object] = []
    policy = "advance-native-cap-reporting-observation/v1"
    if radial:
        policy = "radial-known-reporting-driven-selection/v1"
        sites = [
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
    return {
        "policy": policy,
        "native_cap_decision_site": "compas_cgal._stock_2.engagement_at[2]",
        "engagement_reporting_value_site": "compas_cgal._stock_2.engagement_at[1]",
        "reporting_driven_decision_sites": sites,
    }


def _counts(value: object, field: str) -> Dict[str, object]:
    counts = _object(value, ("observations", "accepted", "exceeded"), field)
    observations = _integer(counts["observations"], f"{field}.observations")
    accepted = _integer(counts["accepted"], f"{field}.accepted")
    exceeded = _integer(counts["exceeded"], f"{field}.exceeded")
    if observations != accepted + exceeded:
        _fail(field, "observations must equal accepted + exceeded")
    return counts


def _station_case(case: Dict[str, object], field: str) -> None:
    reconstruction = _object(
        case["reconstruction"],
        (
            "stock_model",
            "target_centre",
            "centre_decimal_places",
            "target_maximal_radius",
            "radius_decimal_places",
            "occurrence_count",
            "coarse_step",
            "coarse_rungs",
            "refined_radius_sequence",
        ),
        f"{field}.reconstruction",
    )
    _literal(reconstruction["stock_model"], "generator-faithful-radial-pre-bridge/v1", f"{field}.reconstruction.stock_model")
    centre = _array(reconstruction["target_centre"], f"{field}.reconstruction.target_centre")
    if len(centre) != 2:
        _fail(f"{field}.reconstruction.target_centre", "must contain two world coordinates")
    for index, value in enumerate(centre):
        _number(value, f"{field}.reconstruction.target_centre[{index}]")
    for name in ("centre_decimal_places", "radius_decimal_places", "occurrence_count"):
        _integer(reconstruction[name], f"{field}.reconstruction.{name}")
    for name in ("target_maximal_radius", "coarse_step"):
        _number(reconstruction[name], f"{field}.reconstruction.{name}")
    coarse_rungs = _array(reconstruction["coarse_rungs"], f"{field}.reconstruction.coarse_rungs")
    for index, value in enumerate(coarse_rungs):
        _integer(value, f"{field}.reconstruction.coarse_rungs[{index}]")
    sequence = _array(reconstruction["refined_radius_sequence"], f"{field}.reconstruction.refined_radius_sequence")
    for index, value in enumerate(sequence):
        _number(value, f"{field}.reconstruction.refined_radius_sequence[{index}]")

    native = _object(
        case["native_sampled_decisions"],
        ("rung_6_cap_exceeded", "rung_6_cuts_material", "rung_7_cap_exceeded", "rung_7_cuts_material", "refined_candidates"),
        f"{field}.native_sampled_decisions",
    )
    for name in ("rung_6_cap_exceeded", "rung_6_cuts_material", "rung_7_cap_exceeded", "rung_7_cuts_material"):
        _boolean(native[name], f"{field}.native_sampled_decisions.{name}", optional=True)
    refined = _array(native["refined_candidates"], f"{field}.native_sampled_decisions.refined_candidates")
    if len(refined) != len(sequence):
        _fail(field, "refined candidates must align with refined radius sequence")
    for index, value in enumerate(refined):
        row = _object(value, ("radius", "cap_exceeded", "cuts_material"), f"{field}.native_sampled_decisions.refined_candidates[{index}]")
        _number(row["radius"], f"{field}.native_sampled_decisions.refined_candidates[{index}].radius")
        _literal(row["radius"], sequence[index], f"{field}.native_sampled_decisions.refined_candidates[{index}].radius")
        _boolean(row["cap_exceeded"], f"{field}.native_sampled_decisions.refined_candidates[{index}].cap_exceeded")
        _boolean(row["cuts_material"], f"{field}.native_sampled_decisions.refined_candidates[{index}].cuts_material")

    reporting = _object(
        case["reporting_values"],
        (
            "maximal_radius",
            "rung_6_radius",
            "rung_6_peak",
            "rung_7_radius",
            "rung_7_peak",
            "refined_band_min_radius",
            "refined_band_max_radius",
            "forced_peak",
            "rescued_peak",
            "angle_unit",
            "length_unit",
        ),
        f"{field}.reporting_values",
    )
    for name in (
        "maximal_radius",
        "rung_6_radius",
        "rung_6_peak",
        "rung_7_radius",
        "rung_7_peak",
        "refined_band_min_radius",
        "refined_band_max_radius",
        "forced_peak",
        "rescued_peak",
    ):
        _number(reporting[name], f"{field}.reporting_values.{name}", optional=True)
    _literal(reporting["angle_unit"], "degree", f"{field}.reporting_values.angle_unit")
    _literal(reporting["length_unit"], "mm", f"{field}.reporting_values.length_unit")


def _sweep_case(case: Dict[str, object], field: str, discriminator: str, report_keys: Sequence[str]) -> None:
    reconstruction = _object(
        case["reconstruction"],
        ("stock_model", "audit_position_count", "audit_includes_entry_phase", "audit_adds_separate_entry_probe", "excludes_chain_entry_circles", "non_entry_circle_counts"),
        f"{field}.reconstruction",
    )
    _literal(reconstruction["stock_model"], "generator-faithful-radial-pre-bridge/v1", f"{field}.reconstruction.stock_model")
    audit_count = _integer(reconstruction["audit_position_count"], f"{field}.reconstruction.audit_position_count")
    selected_counts = [
        _integer(value, f"{field}.reconstruction.non_entry_circle_counts[{index}]")
        for index, value in enumerate(_array(reconstruction["non_entry_circle_counts"], f"{field}.reconstruction.non_entry_circle_counts"))
    ]
    _literal(reconstruction["audit_includes_entry_phase"], True, f"{field}.reconstruction.audit_includes_entry_phase")
    _literal(reconstruction["audit_adds_separate_entry_probe"], False, f"{field}.reconstruction.audit_adds_separate_entry_probe")
    _literal(reconstruction["excludes_chain_entry_circles"], True, f"{field}.reconstruction.excludes_chain_entry_circles")
    config = cast(Dict[str, object], case["config"])
    _literal(audit_count, len(cast(List[object], cast(Dict[str, object], config["audit"])["probe_offsets"])), f"{field}.reconstruction.audit_position_count")
    native_rows = _array(case["native_sampled_decisions"], f"{field}.native_sampled_decisions")
    report_rows = _array(case["reporting_values"], f"{field}.reporting_values")
    expected_values = cast(List[object], config[{"subdivisions": "subdivision_values", "floor_steps": "floor_values", "refinement_margin": "margin_values"}[discriminator]])
    if len(selected_counts) != len(expected_values) or len(native_rows) != len(expected_values) or len(report_rows) != len(expected_values):
        _fail(field, "reconstruction/native/reporting rows must align with configuration")
    for index, expected in enumerate(expected_values):
        selected_count = selected_counts[index]
        native = _object(native_rows[index], (discriminator, "selected_circles", "observations"), f"{field}.native_sampled_decisions[{index}]")
        _literal(native[discriminator], expected, f"{field}.native_sampled_decisions[{index}].{discriminator}")
        _literal(
            _integer(native["selected_circles"], f"{field}.native_sampled_decisions[{index}].selected_circles"),
            selected_count,
            f"{field}.native_sampled_decisions[{index}].selected_circles",
        )
        observations = _counts(native["observations"], f"{field}.native_sampled_decisions[{index}].observations")
        _literal(observations["observations"], selected_count * audit_count, f"{field}.native_sampled_decisions[{index}].observations.observations")
        report = _object(report_rows[index], report_keys, f"{field}.reporting_values[{index}]")
        _literal(report[discriminator], expected, f"{field}.reporting_values[{index}].{discriminator}")
        _number(report["worst_peak"], f"{field}.reporting_values[{index}].worst_peak")
        circles_over_cap = _integer(report["circles_over_cap"], f"{field}.reporting_values[{index}].circles_over_cap")
        exceeded_positions = cast(int, observations["exceeded"])
        if circles_over_cap > selected_count or circles_over_cap > exceeded_positions or (circles_over_cap == 0) != (exceeded_positions == 0):
            _fail(f"{field}.reporting_values[{index}].circles_over_cap", "must be sound for selected circles and exceeded audit positions")
        _literal(report["angle_unit"], "degree", f"{field}.reporting_values[{index}].angle_unit")
        if discriminator == "refinement_margin":
            _number(report["cutting_length"], f"{field}.reporting_values[{index}].cutting_length")
            _literal(report["length_unit"], "mm", f"{field}.reporting_values[{index}].length_unit")


def _advance_reconstruction(case: Dict[str, object], field: str, audit_count: int) -> Tuple[Dict[str, object], int]:
    reconstruction = _object(
        case["reconstruction"],
        ("stock_model", "selected_circle_policy", "original_calls_per_wrapper", "generator_prepends_entry_probe", "audit_excludes_entry_probe", "selected_circle_count"),
        f"{field}.reconstruction",
    )
    _literal(reconstruction["stock_model"], "generator-pre-bridge/v1", f"{field}.reconstruction.stock_model")
    _literal(reconstruction["selected_circle_policy"], "accepted-non-forced/v1", f"{field}.reconstruction.selected_circle_policy")
    _literal(reconstruction["original_calls_per_wrapper"], 1, f"{field}.reconstruction.original_calls_per_wrapper")
    _literal(reconstruction["generator_prepends_entry_probe"], True, f"{field}.reconstruction.generator_prepends_entry_probe")
    _literal(reconstruction["audit_excludes_entry_probe"], True, f"{field}.reconstruction.audit_excludes_entry_probe")
    selected = _integer(reconstruction["selected_circle_count"], f"{field}.reconstruction.selected_circle_count")
    config = cast(Dict[str, object], case["config"])
    _literal(audit_count, len(cast(List[object], cast(Dict[str, object], config["audit"])["probe_offsets"])), f"{field}.audit count")
    return reconstruction, selected


def _advance_placement(case: Dict[str, object], field: str) -> None:
    _, total_selected = _advance_reconstruction(case, field, 32)
    config = cast(Dict[str, object], case["config"])
    calls = cast(List[object], config["public_calls"])
    native_rows = _array(case["native_sampled_decisions"], f"{field}.native_sampled_decisions")
    reporting_rows = _array(case["reporting_values"], f"{field}.reporting_values")
    if len(native_rows) != len(calls) or len(reporting_rows) != len(calls):
        _fail(field, "placement rows must align with public calls")
    selected_sum = 0
    for index, call_value in enumerate(calls):
        call = cast(Dict[str, object], call_value)
        native = _object(native_rows[index], ("cap", "climb", "selected_circles", "observations"), f"{field}.native_sampled_decisions[{index}]")
        report = _object(
            reporting_rows[index],
            ("cap", "climb", "selected_circles", "positions_over_cap", "worst_peak", "worst_peak_offset", "old_probe_peak", "offset_bin_counts", "angle_unit"),
            f"{field}.reporting_values[{index}]",
        )
        for name in ("cap", "climb"):
            _literal(native[name], call["tea_cap" if name == "cap" else name], f"{field}.native_sampled_decisions[{index}].{name}")
            _literal(report[name], native[name], f"{field}.reporting_values[{index}].{name}")
        selected = _integer(native["selected_circles"], f"{field}.native_sampled_decisions[{index}].selected_circles")
        selected_sum += selected
        _literal(_integer(report["selected_circles"], f"{field}.reporting_values[{index}].selected_circles"), selected, f"{field}.reporting_values[{index}].selected_circles")
        observations = _counts(native["observations"], f"{field}.native_sampled_decisions[{index}].observations")
        _literal(observations["observations"], selected * 32, f"{field}.native_sampled_decisions[{index}].observations.observations")
        over_cap = _integer(report["positions_over_cap"], f"{field}.reporting_values[{index}].positions_over_cap")
        _literal(over_cap, observations["exceeded"], f"{field}.reporting_values[{index}].positions_over_cap")
        for name in ("worst_peak", "worst_peak_offset", "old_probe_peak"):
            _number(report[name], f"{field}.reporting_values[{index}].{name}")
        _literal(report["angle_unit"], "degree", f"{field}.reporting_values[{index}].angle_unit")
        bins = _array(report["offset_bin_counts"], f"{field}.reporting_values[{index}].offset_bin_counts")
        offsets: List[object] = []
        bin_sum = 0
        for bin_index, bin_value in enumerate(bins):
            bin_row = _object(bin_value, ("offset", "positions_over_cap", "angle_unit"), f"{field}.reporting_values[{index}].offset_bin_counts[{bin_index}]")
            _number(bin_row["offset"], f"{field}.reporting_values[{index}].offset_bin_counts[{bin_index}].offset")
            offsets.append(bin_row["offset"])
            bin_sum += _integer(bin_row["positions_over_cap"], f"{field}.reporting_values[{index}].offset_bin_counts[{bin_index}].positions_over_cap")
            _literal(bin_row["angle_unit"], "degree", f"{field}.reporting_values[{index}].offset_bin_counts[{bin_index}].angle_unit")
        _same(offsets, _PLACEMENT_BIN_CENTRES, f"{field}.reporting_values[{index}].offset_bin_counts offsets")
        _literal(bin_sum, over_cap, f"{field}.reporting_values[{index}].offset_bin_counts")
    _literal(selected_sum, total_selected, f"{field}.reconstruction.selected_circle_count")


def _advance_probe_count(case: Dict[str, object], field: str) -> None:
    _, total_selected = _advance_reconstruction(case, field, 60)
    config = cast(Dict[str, object], case["config"])
    configurations = cast(List[object], config["generator_configurations"])
    calls = cast(List[object], config["public_calls"])
    expected = [
        (cast(Dict[str, object], configuration)["label"], cast(Dict[str, object], configuration)["probe_count"], cast(Dict[str, object], call)["tea_cap"])
        for configuration in configurations
        for call in calls
    ]
    native_rows = _array(case["native_sampled_decisions"], f"{field}.native_sampled_decisions")
    reporting_rows = _array(case["reporting_values"], f"{field}.reporting_values")
    if len(native_rows) != len(expected) or len(reporting_rows) != len(expected):
        _fail(field, "probe-count rows must align with configurations and caps")
    selected_sum = 0
    for index, discriminants in enumerate(expected):
        native = _object(native_rows[index], ("label", "probe_count", "cap", "selected_circles", "observations"), f"{field}.native_sampled_decisions[{index}]")
        report = _object(reporting_rows[index], ("label", "probe_count", "cap", "worst_peak", "angle_unit"), f"{field}.reporting_values[{index}]")
        for name, expected_value in zip(("label", "probe_count", "cap"), discriminants):
            _literal(native[name], expected_value, f"{field}.native_sampled_decisions[{index}].{name}")
            _literal(report[name], expected_value, f"{field}.reporting_values[{index}].{name}")
        selected = _integer(native["selected_circles"], f"{field}.native_sampled_decisions[{index}].selected_circles")
        selected_sum += selected
        observations = _counts(native["observations"], f"{field}.native_sampled_decisions[{index}].observations")
        _literal(observations["observations"], selected * 60, f"{field}.native_sampled_decisions[{index}].observations.observations")
        _number(report["worst_peak"], f"{field}.reporting_values[{index}].worst_peak")
        _literal(report["angle_unit"], "degree", f"{field}.reporting_values[{index}].angle_unit")
    _literal(selected_sum, total_selected, f"{field}.reconstruction.selected_circle_count")


def _validate_case(value: object, expected_case: str, index: int) -> Dict[str, object]:
    field = f"cases[{index}]"
    case = _object(value, _COMMON_KEYS, field)
    _literal(case["case"], expected_case, f"{field}.case")
    _same(case["source_claim_ids"], _CASE_CLAIMS[expected_case], f"{field}.source_claim_ids")
    try:
        _same(case["config"], _expected_config(expected_case), f"{field}.config")
    except InvalidMeasurementClaimPayloadError as exc:
        raise InvalidMeasurementClaimConfigError(str(exc)) from exc
    radial = expected_case.startswith("radial-")
    _same(case["selection_decision_provenance"], _expected_provenance(radial), f"{field}.selection_decision_provenance")
    _literal(case["continuous_certificate"], None, f"{field}.continuous_certificate")
    if expected_case == "radial-station":
        _station_case(case, field)
    elif expected_case == "radial-subdivisions":
        _sweep_case(case, field, "subdivisions", ("subdivisions", "worst_peak", "circles_over_cap", "angle_unit"))
    elif expected_case == "radial-floor":
        _sweep_case(case, field, "floor_steps", ("floor_steps", "worst_peak", "circles_over_cap", "angle_unit"))
    elif expected_case == "radial-margin":
        _sweep_case(case, field, "refinement_margin", ("refinement_margin", "worst_peak", "circles_over_cap", "cutting_length", "angle_unit", "length_unit"))
    elif expected_case == "advance-placement":
        _advance_placement(case, field)
    else:
        _advance_probe_count(case, field)
    return case


ClaimEvidencePayloads = Tuple[
    MC001EvidencePayload,
    MC002EvidencePayload,
    MC003EvidencePayload,
    MC004EvidencePayload,
    MC005EvidencePayload,
    MC006EvidencePayload,
    MC007EvidencePayload,
    MC008EvidencePayload,
    MC009EvidencePayload,
    MC010EvidencePayload,
]


def _claim_evidence(cases: Sequence[GeneratorCasePayload]) -> ClaimEvidencePayloads:
    station = cast(RadialStationCasePayload, cases[0])
    subdivisions = cast(RadialSubdivisionsCasePayload, cases[1])
    floor = cast(RadialFloorCasePayload, cases[2])
    margin = cast(RadialMarginCasePayload, cases[3])
    placement = cast(AdvancePlacementCasePayload, cases[4])
    probe_count = cast(AdvanceProbeCountCasePayload, cases[5])
    reconstruction = station["reconstruction"]
    native = station["native_sampled_decisions"]
    reporting = station["reporting_values"]
    occurrence = reconstruction["occurrence_count"]
    evidence_1: MC001EvidencePayload = {
        "occurrence_count": occurrence,
        "station_centre": reconstruction["target_centre"] if occurrence else None,
        "maximal_radius": reporting["maximal_radius"] if occurrence else None,
        "length_unit": "mm",
    }
    evidence_2: MC002EvidencePayload = {
        "coarse_step": reconstruction["coarse_step"],
        "rung_6_radius": reporting["rung_6_radius"],
        "rung_6_peak": reporting["rung_6_peak"],
        "rung_6_cuts_material": native["rung_6_cuts_material"],
        "rung_7_radius": reporting["rung_7_radius"],
        "rung_7_peak": reporting["rung_7_peak"],
        "rung_7_cuts_material": native["rung_7_cuts_material"],
        "angle_unit": "degree",
        "length_unit": "mm",
    }
    evidence_3: MC003EvidencePayload = {
        "rung_7_radius": reporting["rung_7_radius"],
        "rung_7_peak": reporting["rung_7_peak"],
        "rung_7_cuts_material": native["rung_7_cuts_material"],
        "refined_band_min_radius": reporting["refined_band_min_radius"],
        "refined_band_max_radius": reporting["refined_band_max_radius"],
        "forced_peak": reporting["forced_peak"],
        "rescued_peak": reporting["rescued_peak"],
        "angle_unit": "degree",
        "length_unit": "mm",
    }
    evidence_4: MC004EvidencePayload = {
        "audit_position_count": subdivisions["reconstruction"]["audit_position_count"],
        "audit_phase": "entry-angle",
        "includes_entry_phase": True,
        "adds_separate_entry_probe": False,
    }
    subdivision_counts = subdivisions["reconstruction"]["non_entry_circle_counts"]
    if not subdivision_counts or any(count != subdivision_counts[0] for count in subdivision_counts[1:]):
        _fail("cases[1].reconstruction.non_entry_circle_counts", "must be identical to support MC-005")
    evidence_5: MC005EvidencePayload = {
        "non_entry_circle_count": subdivision_counts[0],
        "rows": subdivisions["reporting_values"],
    }
    evidence_6: MC006EvidencePayload = {
        "history_commit": GitObjectId(_HISTORY_COMMIT),
        "missing_configuration": "refinement-without-reporting-ranking",
    }
    evidence_7: MC007EvidencePayload = {
        "audit_position_count": floor["reconstruction"]["audit_position_count"],
        "rows": floor["reporting_values"],
    }
    evidence_8: MC008EvidencePayload = {
        "rows": margin["reporting_values"],
        "baseline_available": False,
        "open_ended_gate_claim_removed": True,
        "reporting_selection_disclosed": True,
    }
    evidence_9: MC009EvidencePayload = {
        "selected_circle_count": placement["reconstruction"]["selected_circle_count"],
        "audit_position_count": len(placement["config"]["audit"]["probe_offsets"]),
        "rows": placement["reporting_values"],
    }
    evidence_10: MC010EvidencePayload = {
        "audit_position_count": len(probe_count["config"]["audit"]["probe_offsets"]),
        "rows": probe_count["reporting_values"],
        "timing_claim_removed": True,
        "relative_cost_claim_removed": True,
    }
    return evidence_1, evidence_2, evidence_3, evidence_4, evidence_5, evidence_6, evidence_7, evidence_8, evidence_9, evidence_10


def compose_generator_payload(source_commit: GitObjectId, cases: Sequence[GeneratorCasePayload]) -> GeneratorClaimPayload:
    """Compose and validate the sole ten-record claim adjudication."""
    if _OBJECT_ID.fullmatch(str(source_commit)) is None:
        _fail("source_commit", "must be one full lowercase Git object ID")
    if str(source_commit) == SOURCE_CORRECTION_COMMIT:
        _fail("source_commit", "must be later than the distinct source correction commit")
    case_list: List[GeneratorCasePayload] = list(cases)
    if len(case_list) != len(_CASE_ORDER):
        _fail("cases", "must contain exactly six cases")
    for index, case_name in enumerate(_CASE_ORDER):
        _literal(case_list[index]["case"], case_name, f"cases[{index}].case")
        _same(case_list[index]["source_claim_ids"], _CASE_CLAIMS[case_name], f"cases[{index}].source_claim_ids")
    station = cast(RadialStationCasePayload, case_list[0])
    subdivisions = cast(RadialSubdivisionsCasePayload, case_list[1])
    floor = cast(RadialFloorCasePayload, case_list[2])
    margin = cast(RadialMarginCasePayload, case_list[3])
    placement = cast(AdvancePlacementCasePayload, case_list[4])
    probe_count = cast(AdvanceProbeCountCasePayload, case_list[5])
    e1, e2, e3, e4, e5, e6, e7, e8, e9, e10 = _claim_evidence(case_list)
    (mc001_disposition, mc001_reason), (mc002_disposition, mc002_reason), (mc003_disposition, mc003_reason) = adjudicate_station_claims(e1, e2, e3)

    subdivision_rows = [(row["subdivisions"], row["worst_peak"], row["circles_over_cap"]) for row in e5["rows"]]
    subdivisions_match = e5["non_entry_circle_count"] == 244 and subdivision_rows == [(value, Degrees(88.6), count) for value, count in zip([1, 2, 4, 8, 16], [12, 12, 8, 8, 8])]
    subdivision_d: Literal["re-earned", "corrected"] = "re-earned" if subdivisions_match else "corrected"
    subdivision_r = "frozen subdivision table matches" if subdivisions_match else "authenticated subdivision table differs"
    floor_matches = e7["audit_position_count"] == 16 and [row["circles_over_cap"] for row in e7["rows"]] == [8, 8, 12]
    floor_d: Literal["re-earned", "corrected"] = "re-earned" if floor_matches else "corrected"
    floor_r = "frozen floor counts match" if floor_matches else "authenticated floor counts differ"
    station_p = station["selection_decision_provenance"]
    subdivision_p = subdivisions["selection_decision_provenance"]
    floor_p = floor["selection_decision_provenance"]
    margin_p = margin["selection_decision_provenance"]
    placement_p = placement["selection_decision_provenance"]
    count_p = probe_count["selection_decision_provenance"]
    reason_4 = "entry-phased audit statement corrected"
    reason_6 = "configuration absent in history"
    reason_8 = "finite sweep corrects open-ended claim"
    reason_9 = "sector tie policy and omitted forward-peak histogram prevent whole-row re-earning"
    reason_10 = "unsupported cost claims removed"
    claim_1 = MC001ClaimPayload(claim_id="MC-001", case="radial-station", disposition=mc001_disposition, reason=mc001_reason, selection_decision_provenance=station_p, evidence=e1)
    claim_2 = MC002ClaimPayload(claim_id="MC-002", case="radial-station", disposition=mc002_disposition, reason=mc002_reason, selection_decision_provenance=station_p, evidence=e2)
    claim_3 = MC003ClaimPayload(claim_id="MC-003", case="radial-station", disposition=mc003_disposition, reason=mc003_reason, selection_decision_provenance=station_p, evidence=e3)
    claim_4 = MC004ClaimPayload(claim_id="MC-004", case="radial-subdivisions", disposition="corrected", reason=reason_4, selection_decision_provenance=subdivision_p, evidence=e4)
    claim_5 = MC005ClaimPayload(
        claim_id="MC-005", case="radial-subdivisions", disposition=subdivision_d, reason=subdivision_r, selection_decision_provenance=subdivision_p, evidence=e5
    )
    claim_6 = MC006ClaimPayload(claim_id="MC-006", case="radial-subdivisions", disposition="historical", reason=reason_6, selection_decision_provenance=subdivision_p, evidence=e6)
    claim_7 = MC007ClaimPayload(claim_id="MC-007", case="radial-floor", disposition=floor_d, reason=floor_r, selection_decision_provenance=floor_p, evidence=e7)
    claim_8 = MC008ClaimPayload(claim_id="MC-008", case="radial-margin", disposition="corrected", reason=reason_8, selection_decision_provenance=margin_p, evidence=e8)
    claim_9 = MC009ClaimPayload(claim_id="MC-009", case="advance-placement", disposition="corrected", reason=reason_9, selection_decision_provenance=placement_p, evidence=e9)
    claim_10 = MC010ClaimPayload(claim_id="MC-010", case="advance-probe-count", disposition="corrected", reason=reason_10, selection_decision_provenance=count_p, evidence=e10)
    payload: GeneratorClaimPayload = {
        "schema_version": "measurement-claim-payload/v2",
        "batch": "generator",
        "extraction_commit": GitObjectId(EXTRACTION_COMMIT),
        "source_commit": source_commit,
        "source_correction_commit": GitObjectId(SOURCE_CORRECTION_COMMIT),
        "case_order": list(_CASE_ORDER),
        "cases": case_list,
        "claims": [claim_1, claim_2, claim_3, claim_4, claim_5, claim_6, claim_7, claim_8, claim_9, claim_10],
    }
    return payload


def _validate_claim(
    value: object,
    claim_id: str,
    case: Dict[str, object],
    expected_evidence: object,
    index: int,
) -> Dict[str, object]:
    field = f"claims[{index}]"
    claim = _object(value, ("claim_id", "case", "disposition", "reason", "selection_decision_provenance", "evidence"), field)
    expected_case = next(case for case, claim_ids in _CASE_CLAIMS.items() if claim_id in claim_ids)
    _literal(claim["claim_id"], claim_id, f"{field}.claim_id")
    _literal(claim["case"], expected_case, f"{field}.case")
    permitted = {
        "MC-001": ("re-earned", "corrected"),
        "MC-002": ("re-earned", "corrected"),
        "MC-003": ("re-earned", "corrected"),
        "MC-004": ("corrected", "historical"),
        "MC-005": ("re-earned", "corrected", "historical"),
        "MC-006": ("historical", "deleted"),
        "MC-007": ("re-earned", "corrected", "historical"),
        "MC-008": ("corrected", "historical"),
        "MC-009": ("re-earned", "corrected", "historical"),
        "MC-010": ("corrected", "historical", "deleted"),
    }[claim_id]
    if type(claim["disposition"]) is not str or claim["disposition"] not in permitted:
        _fail(f"{field}.disposition", f"not permitted for {claim_id}")
    reason = claim["reason"]
    if type(reason) is not str or not reason or any(token in reason for token in ("|", "\r", "\n", "\u2028", "\u2029")):
        _fail(f"{field}.reason", "must be non-empty and Markdown-table-safe")
    _same(claim["selection_decision_provenance"], case["selection_decision_provenance"], f"{field}.selection_decision_provenance")
    _same(claim["evidence"], expected_evidence, f"{field}.evidence")
    return claim


def validate_generator_payload(payload: object) -> GeneratorClaimPayload:
    """Validate the complete six-case, ten-claim payload without schema erasure."""
    _finite_tree(payload, "payload")
    root = _object(
        payload,
        ("schema_version", "batch", "extraction_commit", "source_commit", "source_correction_commit", "case_order", "cases", "claims"),
        "payload",
    )
    _literal(root["schema_version"], "measurement-claim-payload/v2", "schema_version")
    _literal(root["batch"], "generator", "batch")
    _literal(root["extraction_commit"], EXTRACTION_COMMIT, "extraction_commit")
    if type(root["source_commit"]) is not str or _OBJECT_ID.fullmatch(root["source_commit"]) is None:
        _fail("source_commit", "must be one full lowercase Git object ID")
    _literal(root["source_correction_commit"], SOURCE_CORRECTION_COMMIT, "source_correction_commit")
    _same(root["case_order"], _CASE_ORDER, "case_order")
    case_values = _array(root["cases"], "cases")
    if len(case_values) != len(_CASE_ORDER):
        _fail("cases", "must contain exactly six cases")
    cases: Dict[str, Dict[str, object]] = {case_name: _validate_case(case_values[index], case_name, index) for index, case_name in enumerate(_CASE_ORDER)}
    claim_values = _array(root["claims"], "claims")
    if len(claim_values) != len(_CLAIM_IDS):
        _fail("claims", "must contain exactly ten claims")
    typed_cases = cast(List[GeneratorCasePayload], case_values)
    evidence = _claim_evidence(typed_cases)
    for index, claim_id in enumerate(_CLAIM_IDS):
        expected_case = next(case for case, claim_ids in _CASE_CLAIMS.items() if claim_id in claim_ids)
        _validate_claim(claim_values[index], claim_id, cases[expected_case], evidence[index], index)
    canonical = compose_generator_payload(GitObjectId(cast(str, root["source_commit"])), typed_cases)
    _same(claim_values, canonical["claims"], "claims")
    return cast(GeneratorClaimPayload, payload)


def _decode(data: bytes, field: str) -> object:
    return measurement_claim_json.decode_strict(data, field, InvalidMeasurementClaimPayloadError)


def _read(path: pathlib.Path, field: str) -> bytes:
    if not path.is_file() or path.is_symlink():
        raise InvalidMeasurementClaimPayloadError(f"{field}: must be a regular non-symlink file")
    return path.read_bytes()


def _logical_name(actual_name: str) -> str:
    if actual_name.startswith("."):
        marker = ".stage-"
        if marker not in actual_name[1:]:
            raise InvalidMeasurementClaimPayloadError("artifact name is not an owned hidden stage")
        logical, suffix = actual_name[1:].split(marker, 1)
        if not suffix or marker in suffix or _FINAL_NAME.fullmatch(logical) is None:
            raise InvalidMeasurementClaimPayloadError("artifact stage wrapper is malformed")
        return logical
    if _FINAL_NAME.fullmatch(actual_name) is None:
        raise InvalidMeasurementClaimPayloadError("artifact final name is malformed")
    return actual_name


def generator_semantic_input(payload: GeneratorClaimPayload) -> GeneratorClaimInputPayload:
    cases = payload["cases"]
    s = cast(RadialStationCasePayload, cases[0])
    d = cast(RadialSubdivisionsCasePayload, cases[1])
    f = cast(RadialFloorCasePayload, cases[2])
    m = cast(RadialMarginCasePayload, cases[3])
    p = cast(AdvancePlacementCasePayload, cases[4])
    c = cast(AdvanceProbeCountCasePayload, cases[5])
    case_inputs: List[GeneratorCaseInputPayload] = [
        RadialStationCaseInputPayload(
            case=s["case"], source_claim_ids=s["source_claim_ids"], config=s["config"], selection_decision_provenance=s["selection_decision_provenance"]
        ),
        RadialSubdivisionsCaseInputPayload(
            case=d["case"], source_claim_ids=d["source_claim_ids"], config=d["config"], selection_decision_provenance=d["selection_decision_provenance"]
        ),
        RadialFloorCaseInputPayload(case=f["case"], source_claim_ids=f["source_claim_ids"], config=f["config"], selection_decision_provenance=f["selection_decision_provenance"]),
        RadialMarginCaseInputPayload(case=m["case"], source_claim_ids=m["source_claim_ids"], config=m["config"], selection_decision_provenance=m["selection_decision_provenance"]),
        AdvancePlacementCaseInputPayload(
            case=p["case"], source_claim_ids=p["source_claim_ids"], config=p["config"], selection_decision_provenance=p["selection_decision_provenance"]
        ),
        AdvanceProbeCountCaseInputPayload(
            case=c["case"], source_claim_ids=c["source_claim_ids"], config=c["config"], selection_decision_provenance=c["selection_decision_provenance"]
        ),
    ]
    return GeneratorClaimInputPayload(
        extraction_commit=payload["extraction_commit"],
        source_correction_commit=payload["source_correction_commit"],
        case_order=payload["case_order"],
        case_inputs=case_inputs,
    )


def validate_claim_artifact(
    result: pathlib.Path,
) -> Tuple[GeneratorClaimPayload, ValidatedEnvelope, ValidatedArtifactStartedUtc, ValidatedArtifactDirectory]:
    """Authenticate one canonical final or exact Task-4-owned hidden stage."""
    if result.is_symlink():
        raise InvalidMeasurementClaimPayloadError("artifact directory must not be a symlink")
    lexical = result.absolute()
    if lexical.parent.name != "measurement_claim_results" or lexical.parent.parent.name != "benchmarks":
        raise InvalidMeasurementClaimPayloadError("artifact must be directly below benchmarks/measurement_claim_results")
    logical_name = _logical_name(lexical.name)
    repository = lexical.parent.parent.parent
    stamp_path, payload_path = lexical / measurement_artifact.STAMP_NAME, lexical / PAYLOAD_NAME
    stamp_before = _read(stamp_path, "stamp")
    payload_before = _read(payload_path, PAYLOAD_NAME)
    envelope = measurement_artifact.validate_envelope(
        lexical,
        logical_name=logical_name,
        artifact_kind=measurement_artifact.ArtifactKind(ARTIFACT_KIND),
        repository=repository,
    )
    stamp_after = _read(stamp_path, "stamp")
    payload_after = _read(payload_path, PAYLOAD_NAME)
    if stamp_before != stamp_after or payload_before != payload_after:
        raise InvalidMeasurementClaimPayloadError("artifact bytes changed during validation")
    stamp = _object(_decode(stamp_before, "stamp"), measurement_artifact.ENVELOPE_KEYS, "stamp")
    payload = validate_generator_payload(_decode(payload_before, PAYLOAD_NAME))
    _literal(stamp["artifact_kind"], ARTIFACT_KIND, "stamp.artifact_kind")
    argv = _array(stamp["argv"], "stamp.argv")
    if len(argv) != 5 or type(argv[0]) is not str or not argv[0]:
        _fail("stamp.argv", "must contain one executable and the canonical command")
    _same(argv[1:], ["-m", "tools.measurement_claim_probes", "run-generator", "--all"], "stamp.argv[1:]")
    input_identity = _object(stamp["input_identity"], ("version", "payload", "sha256"), "stamp.input_identity")
    result_identity = _object(stamp["result_identity"], ("version", "payloads", "sha256"), "stamp.result_identity")
    _literal(input_identity["version"], INPUT_VERSION, "stamp.input_identity.version")
    _literal(result_identity["version"], RESULT_VERSION, "stamp.result_identity.version")
    payloads = _object(result_identity["payloads"], (PAYLOAD_NAME,), "stamp.result_identity.payloads")
    _object(payloads[PAYLOAD_NAME], ("sha256",), f"stamp.result_identity.payloads.{PAYLOAD_NAME}")
    _same(input_identity["payload"], generator_semantic_input(payload), "stamp.input_identity.payload")
    _literal(payload["source_commit"], str(envelope.commit), "source_commit")
    started_value = stamp["started"]
    if type(started_value) is not str:
        _fail("stamp.started", "must be a canonical UTC timestamp")
    try:
        started = datetime.datetime.fromisoformat(cast(str, started_value))
    except ValueError as exc:
        raise InvalidMeasurementClaimPayloadError("stamp.started: not ISO-8601") from exc
    if started.utcoffset() != datetime.timedelta(0) or started.isoformat(timespec="microseconds") != started_value:
        _fail("stamp.started", "must be canonical timezone-aware UTC with microseconds")
    expected_name = f"{started.date().isoformat()}-{str(envelope.commit)[:12]}-generator-{str(envelope.input_sha256)[:12]}"
    _literal(logical_name, expected_name, "artifact logical name")
    canonical = pathlib.PurePosixPath("benchmarks", "measurement_claim_results", expected_name)
    return payload, envelope, ValidatedArtifactStartedUtc(started), ValidatedArtifactDirectory(canonical)
