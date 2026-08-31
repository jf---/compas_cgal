"""Strict Markdown grammar for the fixed Figure-6 benchmark result."""

from __future__ import annotations

from typing import Optional
from typing import Sequence

from tools.measurement_claim_benchmark_semantic_input import FIGURE6_CAPS
from tools.measurement_claim_benchmark_semantic_input import FIGURE6_SPACINGS
from tools.measurement_claim_errors import InvalidMeasurementClaimPayloadError
from tools.measurement_claim_markdown import parse_pipe_table

_COMPARISON_HEADER = (
    "| cap (deg) | controlled length | controlled max TEA after entry (deg) | controlled meets cap | "
    "controlled exceedances after entry | spacing (tool diam.) | constant-spacing length | "
    "constant-spacing max TEA after entry (deg) | length ratio |"
)
_COMPARISON_RULE = "| ---: | ---: | ---: | :--- | ---: | ---: | ---: | ---: | ---: |"
_COMPARISON_COLUMN_COUNT = 9
_TRIAL_HEADER = "| spacing (tool diam.) | length | cut motions | max TEA after entry (deg) | raw max TEA (deg) |"
_TRIAL_RULE = "| ---: | ---: | ---: | ---: | ---: |"
_TRIAL_COLUMN_COUNT = 5


def _section(lines: Sequence[str], heading: str, next_heading: Optional[str]) -> Sequence[str]:
    start = lines.index(heading) + 1
    end = lines.index(next_heading) if next_heading is not None else len(lines)
    return lines[start:end]


def validate_figure6_markdown(
    markdown_bytes: bytes,
    *,
    pocket_name: str,
    point_count: int,
    trial_count: int,
) -> None:
    """Validate the two owned Figure-6 tables and section boundaries.

    Args:
        markdown_bytes: Raw Figure-6 Markdown bytes.
        pocket_name: Validated pocket name required in the H1.
        point_count: Number of cap records in the JSON authority.
        trial_count: Number of spacing records in the JSON authority.

    Raises:
        InvalidMeasurementClaimPayloadError: UTF-8, headings, table grammar,
            row counts, or section remainder violate the contract.
    """
    try:
        text = markdown_bytes.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise InvalidMeasurementClaimPayloadError("figure6.md: invalid UTF-8") from exc
    if not text.strip():
        raise InvalidMeasurementClaimPayloadError("figure6.md: must not be empty")
    lines = text.splitlines()
    h1 = f"# Figure 6 reproduction — path length against engagement cap ({pocket_name})"
    if [line for line in lines if line.startswith("# ")] != [h1]:
        raise InvalidMeasurementClaimPayloadError("figure6.md: requires the exact Figure-6 H1 once")
    headings = [line for line in lines if line.startswith("## ")]
    if headings != ["## Per-cap comparison", "## Constant-spacing trials"]:
        raise InvalidMeasurementClaimPayloadError("figure6.md: requires exact ordered H2 headings")
    if lines.index(h1) > lines.index(headings[0]) or lines.index(headings[0]) > lines.index(headings[1]):
        raise InvalidMeasurementClaimPayloadError("figure6.md: headings are reordered")
    comparison = parse_pipe_table(
        _section(lines, headings[0], headings[1]),
        header=_COMPARISON_HEADER,
        rule=_COMPARISON_RULE,
        columns=_COMPARISON_COLUMN_COUNT,
        field="figure6.md comparison",
        error=InvalidMeasurementClaimPayloadError,
    )
    trials = parse_pipe_table(
        _section(lines, headings[1], None),
        header=_TRIAL_HEADER,
        rule=_TRIAL_RULE,
        columns=_TRIAL_COLUMN_COUNT,
        field="figure6.md trials",
        error=InvalidMeasurementClaimPayloadError,
    )
    if point_count != len(FIGURE6_CAPS) or len(comparison) != len(FIGURE6_CAPS):
        raise InvalidMeasurementClaimPayloadError(f"figure6.md comparison: requires exactly {len(FIGURE6_CAPS)} rows")
    if trial_count != len(FIGURE6_SPACINGS) or len(trials) != len(FIGURE6_SPACINGS):
        raise InvalidMeasurementClaimPayloadError(f"figure6.md trials: requires exactly {len(FIGURE6_SPACINGS)} rows")
