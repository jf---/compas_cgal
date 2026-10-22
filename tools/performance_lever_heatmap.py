"""The lever x workload heatmap: how much profiled time each performance lever can at most remove.

Reads the committed lever mapping (`tools.performance_levers`) and the profiles
`tools.performance_profile` writes. A cell is the share of the workload's
profiled time spent in the code the lever replaces, an upper bound on its
saving (for a ``wired`` lever: the code that replaced it, the time it still
costs); ``n/a`` where none of that code runs in the workload; hatched and
labelled ``est.`` where the profile cannot attribute it exactly (the mapping
states why). Writes ``docs/assets/images/performance_lever_heatmap.png`` and
``benchmarks/results/performance_lever_shares.json``. Run as
``pixi run performance-figures``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any
from typing import Sequence

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from tools.performance_levers import Lever
from tools.performance_levers import lever_share
from tools.performance_levers import load_levers
from tools.performance_stats import load_stats
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
from tools.performance_workloads import PerformanceWorkload
from tools.performance_workloads import read_profile_meta

LEVER_HEATMAP_PATH = ROOT / "docs" / "assets" / "images" / "performance_lever_heatmap.png"
LEVER_SHARES_PATH = ROOT / "benchmarks" / "results" / "performance_lever_shares.json"
SCHEMA = "performance-lever-shares/1"
FIGURE_SIZE_IN = (12.0, 5.8)
HATCH = "////"
HATCH_COLOR = "#ffffff"
# Hatch opacity: visible on the lightest ramp step, never hiding the label.
HATCH_ALPHA = 0.55
STORED_DECIMALS = 2

Cell = tuple[float | None, str | None]


def lever_cells(levers: Sequence[Lever], workloads: Sequence[PerformanceWorkload]) -> list[list[Cell]]:
    """Every lever's share and estimate reason per workload.

    Args:
        levers: The levers, in row order.
        workloads: The workloads, in column order.

    Returns:
        ``cells[i][j] = (share percent or None for n/a, estimate reason or None)``.
    """
    stats = {workload.key: load_stats(workload.profile_path) for workload in workloads}
    return [[(lever_share(lever, workload.key, stats[workload.key]), lever.estimated.get(workload.key)) for workload in workloads] for lever in levers]


def render_lever_heatmap(levers: Sequence[Lever], workloads: Sequence[PerformanceWorkload], cells: Sequence[Sequence[Cell]], profiled_on: str, path: Path) -> None:
    """Draw the lever x workload heatmap.

    Args:
        levers: Row levers.
        workloads: Column workloads.
        cells: `lever_cells` output.
        profiled_on: The profile date(s) for the footer.
        path: The PNG to write.
    """
    fig, ax = plt.subplots(figsize=FIGURE_SIZE_IN, dpi=FIGURE_DPI)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    for i, row in enumerate(cells):
        for j, (share, reason) in enumerate(row):
            corner = (j - 0.5 + CELL_INSET, i - 0.5 + CELL_INSET)
            if share is None:
                ax.add_patch(Rectangle(corner, CELL_SIDE, CELL_SIDE, color=NA, lw=0))
                ax.text(j, i, "n/a", ha="center", va="center", fontsize=NA_FONT_SIZE, color=MUTED)
                continue
            ax.add_patch(Rectangle(corner, CELL_SIDE, CELL_SIDE, color=CMAP(share / PERCENT), lw=0))
            label = f"{share:.0f} %"
            if reason is not None:
                ax.add_patch(Rectangle(corner, CELL_SIDE, CELL_SIDE, fill=False, hatch=HATCH, edgecolor=HATCH_COLOR, alpha=HATCH_ALPHA, lw=0))
                label = f"≤ {share:.0f} %\nest."
            ax.text(
                j,
                i,
                label,
                ha="center",
                va="center",
                fontsize=VALUE_FONT_SIZE,
                fontweight="bold" if share >= BOLD_FROM_PERCENT else "normal",
                color="white" if share >= WHITE_TEXT_FROM_PERCENT else INK,
                bbox=None
                if reason is None
                else {"boxstyle": "round,pad=0.15", "fc": SURFACE if share < WHITE_TEXT_FROM_PERCENT else CMAP(share / PERCENT), "ec": "none", "alpha": 0.85},
            )
    ax.set_xlim(-0.5, len(workloads) - 0.5)
    ax.set_ylim(len(levers) - 0.5, -0.5)
    ax.set_xticks(range(len(workloads)), [workload.label.replace(" ", "\n", 1) for workload in workloads], fontsize=9, color=INK)
    ax.set_yticks(range(len(levers)), [f"{lever.label}  [{lever.status.value}]" for lever in levers], fontsize=10, color=INK)
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("What each lever can at most remove — profiled share of the code it replaces", loc="left", fontsize=TITLE_FONT_SIZE, color=INK, pad=40)
    fig.text(
        0.01,
        0.015,
        f"cProfile shares, {profiled_on} (Python inflated ~3×) · upper bound on the saving; levers overlap, columns do not sum · "
        "n/a = code does not run there · wired = the replacement's remaining share\n"
        "est. = not attributable in the profile (reason in benchmarks/results/performance_levers.json) · regenerate: pixi run performance-figures",
        fontsize=FOOTER_FONT_SIZE,
        color=MUTED,
    )
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)


def lever_shares_document(levers: Sequence[Lever], workloads: Sequence[PerformanceWorkload], cells: Sequence[Sequence[Cell]], metas: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """The machine-readable lever shares.

    Args:
        levers: Row levers.
        workloads: Column workloads.
        cells: `lever_cells` output.
        metas: Profile sidecars, one per workload.

    Returns:
        The JSON-ready document; ``share_percent`` is ``null`` for n/a.
    """
    return {
        "schema": SCHEMA,
        "generated_by": "tools.performance_lever_heatmap",
        "mapping": "benchmarks/results/performance_levers.json",
        "profiles": {workload.key: {"commit": meta["commit"], "dirty": meta["dirty"], "date": meta["date"]} for workload, meta in zip(workloads, metas)},
        "levers": {
            lever.key: {
                "label": lever.label,
                "status": lever.status.value,
                "cells": {
                    workload.key: {"share_percent": None if share is None else round(share, STORED_DECIMALS), "estimated": reason is not None, "reason": reason}
                    for workload, (share, reason) in zip(workloads, row)
                },
            }
            for lever, row in zip(levers, cells)
        },
    }


def main() -> int:
    """Read the mapping and every profile; write the lever heatmap and its shares document.

    Returns:
        The process exit status.
    """
    levers = load_levers()
    workloads = PERFORMANCE_WORKLOADS
    metas = [read_profile_meta(workload) for workload in workloads]
    cells = lever_cells(levers, workloads)
    dates = " / ".join(sorted({str(meta["date"]) for meta in metas}))
    render_lever_heatmap(levers, workloads, cells, dates, LEVER_HEATMAP_PATH)
    LEVER_SHARES_PATH.write_text(json.dumps(lever_shares_document(levers, workloads, cells, metas), indent=2, ensure_ascii=False) + "\n")
    for lever, row in zip(levers, cells):
        print(lever.key, [None if share is None else round(share, 1) for share, _ in row])
    print(f"wrote {LEVER_HEATMAP_PATH.relative_to(ROOT)} and {LEVER_SHARES_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
