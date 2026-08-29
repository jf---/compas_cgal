"""Sole composition and validation authority for v2 claim payloads."""

from __future__ import annotations

import re
from typing import Dict
from typing import List
from typing import Literal
from typing import Sequence
from typing import Tuple
from typing import cast

from tools.measurement_artifact import GitObjectId
from tools.measurement_claim_adjudication import adjudicate_station_claims
from tools.measurement_claim_case_validation import CASE_CLAIMS
from tools.measurement_claim_case_validation import CASE_ORDER
from tools.measurement_claim_case_validation import fail_payload
from tools.measurement_claim_case_validation import generator_cases_wire_view
from tools.measurement_claim_case_validation import validate_array
from tools.measurement_claim_case_validation import validate_finite_tree
from tools.measurement_claim_case_validation import validate_generator_cases
from tools.measurement_claim_case_validation import validate_literal
from tools.measurement_claim_case_validation import validate_object
from tools.measurement_claim_case_validation import validate_same
from tools.measurement_claim_identity import EXTRACTION_COMMIT
from tools.measurement_claim_identity import SOURCE_CORRECTION_COMMIT
from tools.measurement_claim_schema import AdvancePlacementCasePayload
from tools.measurement_claim_schema import AdvanceProbeCountCasePayload
from tools.measurement_claim_schema import Degrees
from tools.measurement_claim_schema import GeneratorCasePayload
from tools.measurement_claim_schema import GeneratorClaimPayload
from tools.measurement_claim_schema import MC001ClaimPayload
from tools.measurement_claim_schema import MC001EvidencePayload
from tools.measurement_claim_schema import MC002ClaimPayload
from tools.measurement_claim_schema import MC002EvidencePayload
from tools.measurement_claim_schema import MC003ClaimPayload
from tools.measurement_claim_schema import MC003EvidencePayload
from tools.measurement_claim_schema import MC004ClaimPayload
from tools.measurement_claim_schema import MC004EvidencePayload
from tools.measurement_claim_schema import MC005ClaimPayload
from tools.measurement_claim_schema import MC005EvidencePayload
from tools.measurement_claim_schema import MC006ClaimPayload
from tools.measurement_claim_schema import MC006EvidencePayload
from tools.measurement_claim_schema import MC007ClaimPayload
from tools.measurement_claim_schema import MC007EvidencePayload
from tools.measurement_claim_schema import MC008ClaimPayload
from tools.measurement_claim_schema import MC008EvidencePayload
from tools.measurement_claim_schema import MC009ClaimPayload
from tools.measurement_claim_schema import MC009EvidencePayload
from tools.measurement_claim_schema import MC010ClaimPayload
from tools.measurement_claim_schema import MC010EvidencePayload
from tools.measurement_claim_schema import RadialFloorCasePayload
from tools.measurement_claim_schema import RadialMarginCasePayload
from tools.measurement_claim_schema import RadialStationCasePayload
from tools.measurement_claim_schema import RadialSubdivisionsCasePayload

_HISTORY_COMMIT = "29050b01e656ea7bf577b18f7bb50a04ff9a23c9"
_CLAIM_IDS = [f"MC-{index:03d}" for index in range(1, 11)]
_OBJECT_ID = re.compile(r"(?:[0-9a-f]{40}|[0-9a-f]{64})\Z")


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


def claim_evidence(cases: Sequence[GeneratorCasePayload]) -> ClaimEvidencePayloads:
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
        fail_payload("cases[1].reconstruction.non_entry_circle_counts", "must be identical to support MC-005")
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


def _validated_source_commit(value: object) -> GitObjectId:
    if type(value) is not str or _OBJECT_ID.fullmatch(value) is None:
        fail_payload("source_commit", "must be one full lowercase Git object ID")
    if value == SOURCE_CORRECTION_COMMIT:
        fail_payload("source_commit", "must be later than the distinct source correction commit")
    return GitObjectId(cast(str, value))


def compose_generator_payload(source_commit: GitObjectId, cases: Sequence[GeneratorCasePayload]) -> GeneratorClaimPayload:
    """Compose and validate the sole ten-record claim adjudication."""
    validated_source = _validated_source_commit(source_commit)
    validate_generator_cases(generator_cases_wire_view(cases))
    case_list = list(cases)
    return _compose_validated_generator_payload(validated_source, case_list)


def _compose_validated_generator_payload(source_commit: GitObjectId, case_list: List[GeneratorCasePayload]) -> GeneratorClaimPayload:
    station = cast(RadialStationCasePayload, case_list[0])
    subdivisions = cast(RadialSubdivisionsCasePayload, case_list[1])
    floor = cast(RadialFloorCasePayload, case_list[2])
    margin = cast(RadialMarginCasePayload, case_list[3])
    placement = cast(AdvancePlacementCasePayload, case_list[4])
    probe_count = cast(AdvanceProbeCountCasePayload, case_list[5])
    e1, e2, e3, e4, e5, e6, e7, e8, e9, e10 = claim_evidence(case_list)
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
        "case_order": list(CASE_ORDER),
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
    claim = validate_object(value, ("claim_id", "case", "disposition", "reason", "selection_decision_provenance", "evidence"), field)
    expected_case = next(case for case, claim_ids in CASE_CLAIMS.items() if claim_id in claim_ids)
    validate_literal(claim["claim_id"], claim_id, f"{field}.claim_id")
    validate_literal(claim["case"], expected_case, f"{field}.case")
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
        fail_payload(f"{field}.disposition", f"not permitted for {claim_id}")
    reason = claim["reason"]
    if type(reason) is not str or not reason or any(token in reason for token in ("|", "\r", "\n", "\u2028", "\u2029")):
        fail_payload(f"{field}.reason", "must be non-empty and Markdown-table-safe")
    validate_same(claim["selection_decision_provenance"], case["selection_decision_provenance"], f"{field}.selection_decision_provenance")
    validate_same(claim["evidence"], expected_evidence, f"{field}.evidence")
    return claim


def validate_generator_payload(payload: object) -> GeneratorClaimPayload:
    """Validate the complete six-case, ten-claim payload without schema erasure."""
    validate_finite_tree(payload, "payload")
    root = validate_object(
        payload,
        ("schema_version", "batch", "extraction_commit", "source_commit", "source_correction_commit", "case_order", "cases", "claims"),
        "payload",
    )
    validate_literal(root["schema_version"], "measurement-claim-payload/v2", "schema_version")
    validate_literal(root["batch"], "generator", "batch")
    validate_literal(root["extraction_commit"], EXTRACTION_COMMIT, "extraction_commit")
    source_commit = _validated_source_commit(root["source_commit"])
    validate_literal(root["source_correction_commit"], SOURCE_CORRECTION_COMMIT, "source_correction_commit")
    validate_same(root["case_order"], CASE_ORDER, "case_order")
    typed_cases = validate_generator_cases(root["cases"])
    cases: Dict[str, Dict[str, object]] = {case_name: cast(Dict[str, object], typed_cases[index]) for index, case_name in enumerate(CASE_ORDER)}
    claim_values = validate_array(root["claims"], "claims")
    if len(claim_values) != len(_CLAIM_IDS):
        fail_payload("claims", "must contain exactly ten claims")
    evidence = claim_evidence(typed_cases)
    for index, claim_id in enumerate(_CLAIM_IDS):
        expected_case = next(case for case, claim_ids in CASE_CLAIMS.items() if claim_id in claim_ids)
        _validate_claim(claim_values[index], claim_id, cases[expected_case], evidence[index], index)
    canonical = _compose_validated_generator_payload(source_commit, typed_cases)
    validate_same(claim_values, canonical["claims"], "claims")
    return cast(GeneratorClaimPayload, payload)
