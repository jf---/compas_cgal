"""Composition and validation of the four benchmark claim records."""

from __future__ import annotations

from typing import List
from typing import cast

from tools.measurement_artifact import GitObjectId
from tools.measurement_claim_benchmark_figure6 import validate_figure6_json
from tools.measurement_claim_benchmark_identity import EXTRACTION_COMMIT
from tools.measurement_claim_benchmark_identity import HISTORY_COMMIT
from tools.measurement_claim_benchmark_identity import MC013_HISTORY_COMMIT
from tools.measurement_claim_benchmark_markdown import validate_figure6_markdown
from tools.measurement_claim_benchmark_schema import BenchmarkClaimPayload
from tools.measurement_claim_benchmark_schema import BenchmarkClaimRecord
from tools.measurement_claim_benchmark_schema import MC011ClaimPayload
from tools.measurement_claim_benchmark_schema import MC012ClaimPayload
from tools.measurement_claim_benchmark_schema import MC013ClaimPayload
from tools.measurement_claim_benchmark_schema import MC013SelectedValuesPayload
from tools.measurement_claim_benchmark_schema import MC014ClaimPayload
from tools.measurement_claim_benchmark_semantic_input import FIGURE6_COMPARISON_SPACING
from tools.measurement_claim_benchmark_semantic_input import FIGURE6_FINE_SPACING
from tools.measurement_claim_benchmark_semantic_input import FIGURE6_SEMANTIC_COMMAND
from tools.measurement_claim_benchmark_semantic_input import figure6_config
from tools.measurement_claim_case_validation import fail_payload
from tools.measurement_claim_case_validation import validate_object
from tools.measurement_claim_case_validation import validate_same

_PAYLOAD_KEYS = ("schema_version", "batch", "extraction_commit", "source_commit", "semantic_command", "config", "claims")


def compose_benchmark_payload(
    *,
    source_commit: GitObjectId,
    figure6_payload: object,
    figure6_markdown: bytes,
) -> BenchmarkClaimPayload:
    """Derive the four benchmark claims from one validated Figure-6 result.

    Args:
        source_commit: Full commit used to execute the benchmark.
        figure6_payload: Decoded raw Figure-6 JSON.
        figure6_markdown: Raw Figure-6 Markdown bytes.

    Returns:
        Typed benchmark claim payload.

    Raises:
        InvalidMeasurementClaimPayloadError: Either Figure-6 representation is
            malformed or the selected measurements violate their contract.
    """
    pocket_name, fine, comparison, point_count, trial_count = validate_figure6_json(figure6_payload)
    validate_figure6_markdown(
        figure6_markdown,
        pocket_name=pocket_name,
        point_count=point_count,
        trial_count=trial_count,
    )
    selected = MC013SelectedValuesPayload(
        fine_spacing=FIGURE6_FINE_SPACING,
        fine_max_tea_after_entry=fine,
        comparison_spacing=FIGURE6_COMPARISON_SPACING,
        comparison_max_tea_after_entry=comparison,
        angle_unit="degree",
        spacing_unit="tool-diameter",
    )
    claims: List[BenchmarkClaimRecord] = [
        MC011ClaimPayload(
            claim_id="MC-011",
            source="benchmarks/gate.py:58",
            disposition="not-a-claim",
            source_commit=source_commit,
            history_commit=GitObjectId(HISTORY_COMMIT),
            reason="Pocket dimensions identify where the defect was first observed; they are not a measured result.",
            missing_inputs=[],
        ),
        MC012ClaimPayload(
            claim_id="MC-012",
            source="benchmarks/gate.py:66",
            disposition="deleted",
            source_commit=source_commit,
            history_commit=GitObjectId(HISTORY_COMMIT),
            reason="The anecdote is under-specified and cannot be reconstructed without guessing inputs.",
            missing_inputs=["generator", "circle-selection", "entry-treatment", "operation-enumeration"],
        ),
        MC013ClaimPayload(
            claim_id="MC-013",
            source="benchmarks/mathsm.py:47",
            disposition="corrected",
            source_commit=source_commit,
            history_commit=GitObjectId(MC013_HISTORY_COMMIT),
            reason="Authenticated Figure-6 JSON establishes the non-monotone spacing observation.",
            missing_inputs=[],
            selected_values=selected,
        ),
        MC014ClaimPayload(
            claim_id="MC-014",
            source="benchmarks/quality.py:150",
            disposition="not-a-claim",
            source_commit=source_commit,
            history_commit=GitObjectId(HISTORY_COMMIT),
            reason="The shared cap defines metric comparability and reports no empirical value.",
            missing_inputs=[],
        ),
    ]
    return BenchmarkClaimPayload(
        schema_version="measurement-claim-payload/v1",
        batch="benchmark",
        extraction_commit=GitObjectId(EXTRACTION_COMMIT),
        source_commit=source_commit,
        semantic_command=list(FIGURE6_SEMANTIC_COMMAND),
        config=figure6_config(),
        claims=claims,
    )


def validate_benchmark_payload(
    payload: object,
    *,
    figure6_payload: object,
    figure6_markdown: bytes,
) -> BenchmarkClaimPayload:
    """Require exact equality with a freshly composed benchmark payload.

    Args:
        payload: Candidate decoded claim payload.
        figure6_payload: Decoded raw Figure-6 JSON.
        figure6_markdown: Raw Figure-6 Markdown bytes.

    Returns:
        The validated typed payload.

    Raises:
        InvalidMeasurementClaimPayloadError: The source identity or any derived
            field differs from fresh composition.
    """
    root = validate_object(payload, _PAYLOAD_KEYS, "benchmark-claims.json")
    source = root["source_commit"]
    if type(source) is not str or len(source) not in (40, 64) or any(character not in "0123456789abcdef" for character in source):
        fail_payload("benchmark-claims.json.source_commit", "must be a full lowercase Git object ID")
    expected = compose_benchmark_payload(
        source_commit=GitObjectId(cast(str, source)),
        figure6_payload=figure6_payload,
        figure6_markdown=figure6_markdown,
    )
    validate_same(root, expected, "benchmark-claims.json")
    return cast(BenchmarkClaimPayload, root)
