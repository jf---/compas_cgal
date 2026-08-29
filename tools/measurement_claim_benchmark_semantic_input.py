"""Sole stage-free semantic input projection for Figure-6 claims."""

from typing import Dict
from typing import Tuple
from typing import cast

from tools.measurement_artifact import GitObjectId
from tools.measurement_claim_benchmark_identity import EXTRACTION_COMMIT
from tools.measurement_claim_benchmark_identity import HISTORY_COMMIT
from tools.measurement_claim_benchmark_identity import MC013_HISTORY_COMMIT
from tools.measurement_claim_benchmark_schema import BenchmarkClaimPayload
from tools.measurement_claim_benchmark_schema import Figure6ConfigPayload
from tools.measurement_claim_schema import Degrees
from tools.measurement_claim_schema import Millimetres
from tools.measurement_claim_schema import ToolDiameters

FIGURE6_CAPS: Tuple[float, ...] = (20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0)
FIGURE6_SPACINGS: Tuple[float, ...] = (0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6)
FIGURE6_SEMANTIC_COMMAND: Tuple[str, ...] = (
    "-m",
    "benchmarks.cli",
    "figure6",
    "--caps",
    "20",
    "40",
    "60",
    "80",
    "100",
    "120",
    "140",
    "160",
    "--spacings",
    "0.025",
    "0.05",
    "0.075",
    "0.1",
    "0.125",
    "0.15",
    "0.2",
    "0.25",
    "0.3",
    "0.4",
    "0.5",
    "0.6",
)


def figure6_config() -> Figure6ConfigPayload:
    return Figure6ConfigPayload(
        width=Millimetres(20.0),
        height=Millimetres(12.0),
        tool_diameter=Millimetres(2.0),
        holes=[],
        reporting_cap=Degrees(180.0),
        caps=[Degrees(value) for value in FIGURE6_CAPS],
        spacings=[ToolDiameters(value) for value in FIGURE6_SPACINGS],
        length_unit="mm",
        angle_unit="degree",
        spacing_unit="tool-diameter",
    )


def benchmark_semantic_input(source_commit: GitObjectId) -> Dict[str, object]:
    return {
        "extraction_commit": EXTRACTION_COMMIT,
        "source_commit": source_commit,
        "semantic_command": list(FIGURE6_SEMANTIC_COMMAND),
        "config": cast(Dict[str, object], figure6_config()),
        "claim_sources": [
            {
                "claim_id": "MC-011",
                "source": "benchmarks/gate.py:58",
                "disposition": "not-a-claim",
                "history_commit": HISTORY_COMMIT,
                "missing_inputs": [],
            },
            {
                "claim_id": "MC-012",
                "source": "benchmarks/gate.py:66",
                "disposition": "deleted",
                "history_commit": HISTORY_COMMIT,
                "missing_inputs": ["generator", "circle-selection", "entry-treatment", "operation-enumeration"],
            },
            {
                "claim_id": "MC-013",
                "source": "benchmarks/mathsm.py:47",
                "disposition": "corrected",
                "history_commit": MC013_HISTORY_COMMIT,
                "missing_inputs": [],
            },
            {
                "claim_id": "MC-014",
                "source": "benchmarks/quality.py:150",
                "disposition": "not-a-claim",
                "history_commit": HISTORY_COMMIT,
                "missing_inputs": [],
            },
        ],
    }


def benchmark_payload_semantic_input(payload: BenchmarkClaimPayload) -> Dict[str, object]:
    claim_sources = []
    for claim in payload["claims"]:
        record = cast(Dict[str, object], claim)
        claim_sources.append(
            {
                "claim_id": record["claim_id"],
                "source": record["source"],
                "disposition": record["disposition"],
                "history_commit": record["history_commit"],
                "missing_inputs": record["missing_inputs"],
            }
        )
    return {
        "extraction_commit": payload["extraction_commit"],
        "source_commit": payload["source_commit"],
        "semantic_command": payload["semantic_command"],
        "config": cast(Dict[str, object], payload["config"]),
        "claim_sources": claim_sources,
    }
