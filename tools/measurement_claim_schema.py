"""Typed schema for authenticated generator-claim artifacts."""

from __future__ import annotations

import datetime
import pathlib
from typing import List
from typing import Literal
from typing import NewType
from typing import Optional
from typing import Tuple
from typing import TypedDict
from typing import Union

from tools.measurement_artifact import GitObjectId
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
