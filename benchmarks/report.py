"""Emit the corpus result as markdown and JSON.

The conclusion is emitted BEFORE any table, by construction, so the ordering
survives every re-run rather than depending on whoever last edited the output.

The cap appears as two columns, never one. `truly_exceeding` and `uncertified`
answer opposite questions and lean in opposite directions, and a report carrying
only one of them misstates the other -- which is exactly how a generator that had
improved came to read as a regression. The legend below the conclusion is emitted
unconditionally for the same reason: a reader who meets these numbers cold must
not be able to mistake one for the other.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import List
from typing import Tuple

from benchmarks.measurement import MeasurementRecord

MARKDOWN_NAME = "benchmark_report.md"
JSON_NAME = "benchmark_report.json"

# Emitted verbatim under the conclusion. The distinction it draws is the one this
# module exists to keep legible, so it is part of the artifact, not a comment.
CAP_COLUMN_LEGEND = (
    "**truly exceeding** counts cut motions where the exact engagement predicate "
    "fired at a sampled cutter position: a demonstrated LOWER BOUND on how many "
    "motions are over the cap, never a certificate. **uncertified** counts "
    "operations whose cap could not be *proved*, which includes every operation "
    "the certifier declined to measure because its guard would not close; it is "
    "sound and over-counts. Neither number estimates the other."
)

TABLE_HEADER = "| instance | family | generate (s) | certify (s) | stations | max TEA (deg) | truly exceeding | uncertified | unresolved | arr. vertices | max digits |"
TABLE_RULE = "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"


def render_markdown(records: List[MeasurementRecord]) -> str:
    """Render the corpus result, conclusion first.

    Args:
        records: Every measurement, in any order.

    Returns:
        The markdown document.
    """
    ok = [r for r in records if r.error is None]
    bad = [r for r in records if r.error is not None]
    lines: List[str] = ["# Benchmark corpus result", ""]
    lines += _conclusion(ok, bad)
    lines += [CAP_COLUMN_LEGEND, ""]
    lines += ["## Per-instance measurements", "", TABLE_HEADER, TABLE_RULE]
    for r in sorted(ok, key=lambda x: -x.certify_seconds):
        lines.append(
            f"| {r.name} | {r.family} | {r.generate_seconds:.3f} | {r.certify_seconds:.2f} | {r.stations} | "
            f"{r.max_tea_deg:.1f} | {r.truly_exceeding} | {r.uncertified} | {r.unresolved} | "
            f"{r.arrangement_vertices_final} | {r.max_coordinate_digits} |"
        )
    lines.append("")

    if bad:
        lines += ["## Instances that failed to measure", "", "| instance | family | error |", "| --- | --- | --- |"]
        for r in bad:
            lines.append(f"| {r.name} | {r.family} | {r.error} |")
        lines.append("")

    return "\n".join(lines)


def _conclusion(ok: List[MeasurementRecord], bad: List[MeasurementRecord]) -> List[str]:
    """The answer, stated before any evidence.

    Exceedance leads because it is the question about the toolpath; certifiability
    follows because it is the question about the certifier. Reversing them, or
    dropping either, is how the two get conflated.

    Args:
        ok: Measured instances.
        bad: Instances that failed to measure.

    Returns:
        The conclusion paragraph as markdown lines, blank-line terminated.
    """
    if not ok:
        return [f"Every one of the {len(bad)} instances failed to measure; no conclusion is available.", ""]

    slowest = max(ok, key=lambda r: r.certify_seconds)
    total_certify = sum(r.certify_seconds for r in ok)
    total_generate = sum(r.generate_seconds for r in ok)
    exceeding = [r for r in ok if r.truly_exceeding > 0]
    exceeding_motions = sum(r.truly_exceeding for r in ok)
    # Disjoint from `exceeding` on purpose: an instance with a demonstrated
    # exceedance is necessarily also uncertified, so counting both sets whole
    # would double-report it and reproduce the conflation this split removed.
    unproved_only = [r for r in ok if r.uncertified > 0 and r.truly_exceeding == 0]
    ratio = (total_certify / total_generate) if total_generate > 0.0 else float("inf")
    return [
        f"**{len(exceeding)} of {len(ok)}** measured instances put the tool over its engagement cap, "
        f"at **{exceeding_motions}** cut motions in total. A further **{len(unproved_only)}** could not be proved "
        f"under the cap without any exceedance being demonstrated. "
        f"The slowest instance is **{slowest.name}** at **{slowest.certify_seconds:.2f} s** to certify "
        f"({slowest.stations} stations); across all {len(ok)} instances certification costs **{ratio:.0f}x** generation "
        f"({total_certify:.1f} s against {total_generate:.2f} s). {len(bad)} instances failed to measure.",
        "",
    ]


def write_report(records: List[MeasurementRecord], out_dir: Path) -> Tuple[Path, Path]:
    """Write the markdown and JSON artifacts.

    Args:
        records: Every measurement.
        out_dir: Destination directory; created when missing.

    Returns:
        ``(markdown_path, json_path)``.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    md_path = out_dir / MARKDOWN_NAME
    json_path = out_dir / JSON_NAME
    md_path.write_text(render_markdown(records))
    json_path.write_text(json.dumps([r.to_dict() for r in records], indent=2))
    return md_path, json_path
