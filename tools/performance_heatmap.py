"""Where generation time goes: the workload x cost-bucket heatmap and its shares document.

Reads the profiles `tools.performance_profile` writes and splits each workload's
profiled time into disjoint buckets that sum to 100 %:

- Held lanes: cumulative-time subtrees (``load_held_reference_case``,
  ``_native_station_boundary_analyses`` (the wired clearance: the guide's one
  native AABB batch),
  ``build_guide_distance_admissible_hypotheses`` (every station's hypotheses
  from that batch) minus clearance,
  ``build_held_reference_path`` minus hypotheses, remainder).
- Stock lanes (pipeline, rho): self time per module of the native-call
  wrappers (``admissibility.py`` / ``_station_is_admissible`` = exact engagement
  predicates; ``stock.py`` and ``engagement.py`` ``_subtract*`` = stock
  depletion), remainder = Python search and orchestration.

Shares, not seconds: cProfile inflates Python roughly threefold. Writes
``docs/assets/images/performance_time_heatmap.png`` and
``benchmarks/results/performance_time_buckets.json``. Run as
``pixi run performance-figures``.
"""

from __future__ import annotations

import hashlib
import json
import sys
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any
from typing import Sequence

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from tools.performance_stats import MissingProfileFunctionError
from tools.performance_stats import ProfileStats
from tools.performance_stats import cumulative_s
from tools.performance_stats import load_stats
from tools.performance_stats import module_name
from tools.performance_style import BOLD_FROM_PERCENT
from tools.performance_style import CELL_INSET
from tools.performance_style import CELL_SIDE
from tools.performance_style import CMAP
from tools.performance_style import FIGURE_DPI
from tools.performance_style import FOOTER_FONT_SIZE
from tools.performance_style import INK
from tools.performance_style import MUTED
from tools.performance_style import NA
from tools.performance_style import NA_FONT_SIZE
from tools.performance_style import PERCENT
from tools.performance_style import SURFACE
from tools.performance_style import TITLE_FONT_SIZE
from tools.performance_style import VALUE_FONT_SIZE
from tools.performance_style import WHITE_TEXT_FROM_PERCENT
from tools.performance_workloads import PERFORMANCE_WORKLOADS
from tools.performance_workloads import ROOT
from tools.performance_workloads import Lane
from tools.performance_workloads import PerformanceWorkload
from tools.performance_workloads import read_profile_meta

HEATMAP_PATH = ROOT / "docs" / "assets" / "images" / "performance_time_heatmap.png"
BUCKETS_PATH = ROOT / "benchmarks" / "results" / "performance_time_buckets.json"
SCHEMA = "performance-time-buckets/1"
# A row's buckets partition the total exactly up to float summation; half a
# percentage point is far above that and far below any rendered difference.
ROW_SUM_TOLERANCE_PERCENT = 0.5


class Bucket(Enum):
    """The disjoint cost buckets, in column order; values are the JSON keys."""

    CLEARANCE = "boundary_clearance"
    CASE_LOAD = "case_load"
    HYPOTHESES = "hypotheses_selection"
    PREDICATES = "exact_engagement_predicates"
    DEPLETION = "stock_depletion"
    ORCHESTRATION = "search_orchestration"


COLUMN_LABELS: dict[Bucket, str] = {
    Bucket.CLEARANCE: "Boundary clearance\n(native AABB)",
    Bucket.CASE_LOAD: "Case load\n(biarc fit)",
    Bucket.HYPOTHESES: "Hypotheses &\nselection",
    Bucket.PREDICATES: "Exact engagement\npredicates",
    Bucket.DEPLETION: "Stock\ndepletion",
    Bucket.ORCHESTRATION: "Search &\norchestration",
}
HELD_BUCKETS = (Bucket.CLEARANCE, Bucket.CASE_LOAD, Bucket.HYPOTHESES, Bucket.ORCHESTRATION)
STOCK_BUCKETS = (Bucket.PREDICATES, Bucket.DEPLETION, Bucket.ORCHESTRATION)

# Held lane: the functions whose cumulative subtrees form the buckets.
HELD_CASE_LOAD = "load_held_reference_case"
HELD_CLEARANCE = "_native_station_boundary_analyses"
HELD_HYPOTHESES = "build_guide_distance_admissible_hypotheses"
HELD_PLACEMENT = "build_held_reference_path"
# Stock lane: the native-call wrappers whose self time forms the buckets.
PREDICATE_MODULE = "admissibility.py"
PREDICATE_FUNCTION = "_station_is_admissible"
DEPLETION_MODULE = "stock.py"
DEPLETION_SUBTRACT_MODULE = "engagement.py"
DEPLETION_SUBTRACT_PREFIX = "_subtract"

FIGURE_SIZE_IN = (11.0, 4.6)


class RowSumError(Exception):
    """A workload's buckets do not sum to 100 %."""


class NegativeBucketError(Exception):
    """A bucket came out negative: the subtrees it subtracts are not nested as assumed."""


@dataclass(frozen=True)
class BucketShares:
    """One workload's split of profiled time.

    Attributes:
        shares: Percent per applicable bucket; buckets not in the lane are absent.
        total_s: Total profiled seconds.
    """

    shares: dict[Bucket, float]
    total_s: float


def _checked(seconds: dict[Bucket, float], total: float) -> BucketShares:
    shares = {bucket: PERCENT * value / total for bucket, value in seconds.items()}
    negative = {bucket.value: share for bucket, share in shares.items() if share < 0.0}
    if negative:
        raise NegativeBucketError(f"negative buckets {negative}")
    row_sum = sum(shares.values())
    if abs(row_sum - PERCENT) > ROW_SUM_TOLERANCE_PERCENT:
        raise RowSumError(f"buckets sum to {row_sum:.3f} %, not {PERCENT:.0f} %")
    return BucketShares(shares=shares, total_s=total)


def held_buckets(stats: ProfileStats) -> BucketShares:
    """Split a Held workload's profile into cumulative-subtree buckets.

    Args:
        stats: The profile.

    Returns:
        The shares of `HELD_BUCKETS`.

    Raises:
        MissingProfileFunctionError: A bucket function is not in the profile.
        RowSumError: The buckets do not sum to 100 %.
        NegativeBucketError: A subtree difference came out negative.
    """
    total = stats.total_tt
    load, clearance = cumulative_s(stats, HELD_CASE_LOAD), cumulative_s(stats, HELD_CLEARANCE)
    hypotheses, placement = cumulative_s(stats, HELD_HYPOTHESES), cumulative_s(stats, HELD_PLACEMENT)
    seconds = {
        Bucket.CLEARANCE: clearance,
        Bucket.CASE_LOAD: load,
        Bucket.HYPOTHESES: (hypotheses - clearance) + (placement - hypotheses),
        Bucket.ORCHESTRATION: total - load - placement,
    }
    return _checked(seconds, total)


def stock_buckets(stats: ProfileStats) -> BucketShares:
    """Split a stock-lane workload's profile into wrapper self-time buckets.

    Args:
        stats: The profile.

    Returns:
        The shares of `STOCK_BUCKETS`.

    Raises:
        MissingProfileFunctionError: No profiled function falls in the predicate or the depletion bucket.
        RowSumError: The buckets do not sum to 100 %.
    """
    total = stats.total_tt
    predicates = depletion = 0.0
    predicate_hits = depletion_hits = 0
    for (filename, _, function), row in stats.stats.items():
        module = module_name(filename)
        if module == PREDICATE_MODULE or function == PREDICATE_FUNCTION:
            predicates += row[2]
            predicate_hits += 1
        elif module == DEPLETION_MODULE or (module == DEPLETION_SUBTRACT_MODULE and function.startswith(DEPLETION_SUBTRACT_PREFIX)):
            depletion += row[2]
            depletion_hits += 1
    if predicate_hits == 0:
        raise MissingProfileFunctionError(f"no {PREDICATE_MODULE} function and no {PREDICATE_FUNCTION} in the profile; was it renamed?")
    if depletion_hits == 0:
        raise MissingProfileFunctionError(f"no {DEPLETION_MODULE} or {DEPLETION_SUBTRACT_MODULE}:{DEPLETION_SUBTRACT_PREFIX}* function in the profile; was it renamed?")
    seconds = {Bucket.PREDICATES: predicates, Bucket.DEPLETION: depletion, Bucket.ORCHESTRATION: total - predicates - depletion}
    return _checked(seconds, total)


def workload_buckets(workload: PerformanceWorkload, stats: ProfileStats) -> BucketShares:
    """The buckets of *workload*'s lane applied to *stats*."""
    return held_buckets(stats) if workload.lane is Lane.HELD else stock_buckets(stats)


def render_heatmap(rows: Sequence[tuple[PerformanceWorkload, BucketShares]], profiled_on: str, path: Path) -> None:
    """Draw the workload x bucket heatmap.

    Args:
        rows: Workloads and their shares, in row order.
        profiled_on: The profile date(s) for the footer.
        path: The PNG to write.
    """
    columns = list(Bucket)
    fig, ax = plt.subplots(figsize=FIGURE_SIZE_IN, dpi=FIGURE_DPI)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    for i, (_, shares) in enumerate(rows):
        for j, bucket in enumerate(columns):
            corner = (j - 0.5 + CELL_INSET, i - 0.5 + CELL_INSET)
            if bucket not in shares.shares:
                ax.add_patch(Rectangle(corner, CELL_SIDE, CELL_SIDE, color=NA, lw=0))
                ax.text(j, i, "n/a", ha="center", va="center", fontsize=NA_FONT_SIZE, color=MUTED)
                continue
            value = shares.shares[bucket]
            ax.add_patch(Rectangle(corner, CELL_SIDE, CELL_SIDE, color=CMAP(value / PERCENT), lw=0))
            ax.text(
                j,
                i,
                f"{value:.0f} %",
                ha="center",
                va="center",
                fontsize=VALUE_FONT_SIZE,
                fontweight="bold" if value >= BOLD_FROM_PERCENT else "normal",
                color="white" if value >= WHITE_TEXT_FROM_PERCENT else INK,
            )
    ax.set_xlim(-0.5, len(columns) - 0.5)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_xticks(range(len(columns)), [COLUMN_LABELS[bucket] for bucket in columns], fontsize=9, color=INK)
    ax.set_yticks(range(len(rows)), [workload.label for workload, _ in rows], fontsize=10, color=INK)
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("Where generation time goes — share of each workload's profiled time", loc="left", fontsize=TITLE_FONT_SIZE, color=INK, pad=48)
    fig.text(
        0.01,
        0.015,
        f"cProfile shares, {profiled_on} (Python inflated ~3×) · rows sum to 100 % · n/a = stage does not exist in that lane · "
        "sources: build/performance/*.prof · regenerate: pixi run performance-figures",
        fontsize=FOOTER_FONT_SIZE,
        color=MUTED,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)


def buckets_document(rows: Sequence[tuple[PerformanceWorkload, BucketShares, dict[str, Any]]]) -> dict[str, Any]:
    """The machine-readable shares document.

    Args:
        rows: Workloads, their shares and their profile sidecars.

    Returns:
        The JSON-ready document; a bucket that does not apply to a lane is ``null``.
    """
    workloads: dict[str, Any] = {}
    for workload, shares, meta in rows:
        workloads[workload.key] = {
            "label": workload.label,
            "lane": workload.lane.value,
            "target": workload.target,
            "profile": str(workload.profile_path.relative_to(ROOT)),
            "profile_sha256": hashlib.sha256(workload.profile_path.read_bytes()).hexdigest(),
            "commit": meta["commit"],
            "dirty": meta["dirty"],
            "date": meta["date"],
            "profiled_total_s": round(shares.total_s, 3),
            "shares_percent": {bucket.value: (round(shares.shares[bucket], 2) if bucket in shares.shares else None) for bucket in Bucket},
        }
    return {"schema": SCHEMA, "generated_by": "tools.performance_heatmap", "buckets": [bucket.value for bucket in Bucket], "workloads": workloads}


def main() -> int:
    """Read every workload profile; write the heatmap and the shares document.

    Returns:
        The process exit status.
    """
    rows = []
    for workload in PERFORMANCE_WORKLOADS:
        meta = read_profile_meta(workload)
        rows.append((workload, workload_buckets(workload, load_stats(workload.profile_path)), meta))
    dates = sorted({str(meta["date"]) for _, _, meta in rows})
    render_heatmap([(workload, shares) for workload, shares, _ in rows], " / ".join(dates), HEATMAP_PATH)
    document = buckets_document(rows)
    BUCKETS_PATH.write_text(json.dumps(document, indent=2, ensure_ascii=False) + "\n")
    for key, entry in document["workloads"].items():
        print(key, {bucket: share for bucket, share in entry["shares_percent"].items() if share is not None})
    print(f"wrote {HEATMAP_PATH.relative_to(ROOT)} and {BUCKETS_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
