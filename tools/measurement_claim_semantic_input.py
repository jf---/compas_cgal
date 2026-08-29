"""Canonical semantic-input projection for v2 claim payloads."""

from typing import List
from typing import cast

from tools.measurement_claim_schema import AdvancePlacementCaseInputPayload
from tools.measurement_claim_schema import AdvancePlacementCasePayload
from tools.measurement_claim_schema import AdvanceProbeCountCaseInputPayload
from tools.measurement_claim_schema import AdvanceProbeCountCasePayload
from tools.measurement_claim_schema import GeneratorCaseInputPayload
from tools.measurement_claim_schema import GeneratorClaimInputPayload
from tools.measurement_claim_schema import GeneratorClaimPayload
from tools.measurement_claim_schema import RadialFloorCaseInputPayload
from tools.measurement_claim_schema import RadialFloorCasePayload
from tools.measurement_claim_schema import RadialMarginCaseInputPayload
from tools.measurement_claim_schema import RadialMarginCasePayload
from tools.measurement_claim_schema import RadialStationCaseInputPayload
from tools.measurement_claim_schema import RadialStationCasePayload
from tools.measurement_claim_schema import RadialSubdivisionsCaseInputPayload
from tools.measurement_claim_schema import RadialSubdivisionsCasePayload


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
