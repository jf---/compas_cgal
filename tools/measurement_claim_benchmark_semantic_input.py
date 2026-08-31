"""Sole stage-free semantic input projection for Figure-6 claims."""

from typing import Final
from typing import List
from typing import Tuple

from tools.measurement_artifact import GitObjectId
from tools.measurement_claim_benchmark_identity import EXTRACTION_COMMIT
from tools.measurement_claim_benchmark_identity import HISTORY_COMMIT
from tools.measurement_claim_benchmark_identity import MC013_HISTORY_COMMIT
from tools.measurement_claim_benchmark_schema import BenchmarkClaimInputPayload
from tools.measurement_claim_benchmark_schema import BenchmarkClaimPayload
from tools.measurement_claim_benchmark_schema import BenchmarkClaimSourcePayload
from tools.measurement_claim_benchmark_schema import Figure6ConfigPayload
from tools.measurement_claim_schema import Degrees
from tools.measurement_claim_schema import Millimetres
from tools.measurement_claim_schema import ToolDiameters

FIGURE6_CAPS: Tuple[Degrees, ...] = tuple(Degrees(value) for value in (20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0))
FIGURE6_FINE_SPACING: Final = ToolDiameters(0.025)
FIGURE6_COMPARISON_SPACING: Final = ToolDiameters(0.1)
FIGURE6_SPACINGS: Tuple[ToolDiameters, ...] = (
    FIGURE6_FINE_SPACING,
    ToolDiameters(0.05),
    ToolDiameters(0.075),
    FIGURE6_COMPARISON_SPACING,
    ToolDiameters(0.125),
    ToolDiameters(0.15),
    ToolDiameters(0.2),
    ToolDiameters(0.25),
    ToolDiameters(0.3),
    ToolDiameters(0.4),
    ToolDiameters(0.5),
    ToolDiameters(0.6),
)
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
    """Build the fixed unit-bearing Figure-6 configuration.

    Returns:
        Typed JSON-wire configuration for the benchmark child.
    """
    return Figure6ConfigPayload(
        width=Millimetres(20.0),
        height=Millimetres(12.0),
        tool_diameter=Millimetres(2.0),
        holes=[],
        reporting_cap=Degrees(180.0),
        caps=list(FIGURE6_CAPS),
        spacings=list(FIGURE6_SPACINGS),
        length_unit="mm",
        angle_unit="degree",
        spacing_unit="tool-diameter",
    )


def benchmark_semantic_input(source_commit: GitObjectId) -> BenchmarkClaimInputPayload:
    """Build the stage-free semantic input for one source commit.

    Args:
        source_commit: Full commit used to execute Figure 6.

    Returns:
        Typed semantic input for the authenticated envelope.
    """
    return BenchmarkClaimInputPayload(
        extraction_commit=GitObjectId(EXTRACTION_COMMIT),
        source_commit=source_commit,
        semantic_command=list(FIGURE6_SEMANTIC_COMMAND),
        config=figure6_config(),
        claim_sources=[
            BenchmarkClaimSourcePayload(
                claim_id="MC-011",
                source="benchmarks/gate.py:58",
                disposition="not-a-claim",
                history_commit=GitObjectId(HISTORY_COMMIT),
                missing_inputs=[],
            ),
            BenchmarkClaimSourcePayload(
                claim_id="MC-012",
                source="benchmarks/gate.py:66",
                disposition="deleted",
                history_commit=GitObjectId(HISTORY_COMMIT),
                missing_inputs=["generator", "circle-selection", "entry-treatment", "operation-enumeration"],
            ),
            BenchmarkClaimSourcePayload(
                claim_id="MC-013",
                source="benchmarks/mathsm.py:47",
                disposition="corrected",
                history_commit=GitObjectId(MC013_HISTORY_COMMIT),
                missing_inputs=[],
            ),
            BenchmarkClaimSourcePayload(
                claim_id="MC-014",
                source="benchmarks/quality.py:150",
                disposition="not-a-claim",
                history_commit=GitObjectId(HISTORY_COMMIT),
                missing_inputs=[],
            ),
        ],
    )


def benchmark_payload_semantic_input(payload: BenchmarkClaimPayload) -> BenchmarkClaimInputPayload:
    """Project one result payload onto its authenticated semantic input.

    Args:
        payload: Validated benchmark result payload.

    Returns:
        Typed semantic input reconstructed from result fields.
    """
    claim_sources: List[BenchmarkClaimSourcePayload] = []
    for claim in payload["claims"]:
        claim_sources.append(
            BenchmarkClaimSourcePayload(
                claim_id=claim["claim_id"],
                source=claim["source"],
                disposition=claim["disposition"],
                history_commit=claim["history_commit"],
                missing_inputs=claim["missing_inputs"],
            )
        )
    return BenchmarkClaimInputPayload(
        extraction_commit=payload["extraction_commit"],
        source_commit=payload["source_commit"],
        semantic_command=payload["semantic_command"],
        config=payload["config"],
        claim_sources=claim_sources,
    )
