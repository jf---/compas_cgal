"""Typed Figure-6 benchmark claim boundary."""

from typing import List
from typing import Literal
from typing import TypedDict
from typing import Union

from tools.measurement_artifact import GitObjectId
from tools.measurement_claim_schema import Degrees
from tools.measurement_claim_schema import Millimetres
from tools.measurement_claim_schema import ToolDiameters
from tools.measurement_claim_schema import WorldPointMillimetres

MissingBenchmarkInput = Literal["generator", "circle-selection", "entry-treatment", "operation-enumeration"]


class Figure6ConfigPayload(TypedDict):
    width: Millimetres
    height: Millimetres
    tool_diameter: Millimetres
    holes: List[List[WorldPointMillimetres]]
    reporting_cap: Degrees
    caps: List[Degrees]
    spacings: List[ToolDiameters]
    length_unit: Literal["mm"]
    angle_unit: Literal["degree"]
    spacing_unit: Literal["tool-diameter"]


class MC013SelectedValuesPayload(TypedDict):
    fine_spacing: ToolDiameters
    fine_max_tea_after_entry: Degrees
    comparison_spacing: ToolDiameters
    comparison_max_tea_after_entry: Degrees
    angle_unit: Literal["degree"]
    spacing_unit: Literal["tool-diameter"]


class MC011ClaimPayload(TypedDict):
    claim_id: Literal["MC-011"]
    source: Literal["benchmarks/gate.py:58"]
    disposition: Literal["not-a-claim"]
    source_commit: GitObjectId
    history_commit: GitObjectId
    reason: str
    missing_inputs: List[MissingBenchmarkInput]


class MC012ClaimPayload(TypedDict):
    claim_id: Literal["MC-012"]
    source: Literal["benchmarks/gate.py:66"]
    disposition: Literal["deleted"]
    source_commit: GitObjectId
    history_commit: GitObjectId
    reason: str
    missing_inputs: List[MissingBenchmarkInput]


class MC013ClaimPayload(TypedDict):
    claim_id: Literal["MC-013"]
    source: Literal["benchmarks/mathsm.py:47"]
    disposition: Literal["corrected"]
    source_commit: GitObjectId
    history_commit: GitObjectId
    reason: str
    missing_inputs: List[MissingBenchmarkInput]
    selected_values: MC013SelectedValuesPayload


class MC014ClaimPayload(TypedDict):
    claim_id: Literal["MC-014"]
    source: Literal["benchmarks/quality.py:150"]
    disposition: Literal["not-a-claim"]
    source_commit: GitObjectId
    history_commit: GitObjectId
    reason: str
    missing_inputs: List[MissingBenchmarkInput]


BenchmarkClaimRecord = Union[MC011ClaimPayload, MC012ClaimPayload, MC013ClaimPayload, MC014ClaimPayload]


class BenchmarkClaimPayload(TypedDict):
    schema_version: Literal["measurement-claim-payload/v1"]
    batch: Literal["benchmark"]
    extraction_commit: GitObjectId
    source_commit: GitObjectId
    semantic_command: List[str]
    config: Figure6ConfigPayload
    claims: List[BenchmarkClaimRecord]
