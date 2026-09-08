"""Render the curved-circle performance figure from a benchmark record and its baseline."""

from __future__ import annotations

import argparse
import math
import pathlib
from dataclasses import dataclass

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from benchmarks.held_curved_circle_benchmark import CurvedCircleBenchmark
from benchmarks.held_curved_circle_benchmark import read_json

CURRENT_COLOR = "#087e8b"
BASELINE_COLOR = "#9aa5a4"
STOP_COLOR = "#c0392b"
BACKGROUND = "#faf9f5"


@dataclass(frozen=True)
class FigureSummary:
    cases: tuple[str, ...]
    stops: dict[str, int]


def render_performance_figure(current: CurvedCircleBenchmark, baseline: CurvedCircleBenchmark | None, output: pathlib.Path) -> FigureSummary:
    """One panel per case: per-piece query time (ms, log) by piece index; baseline in grey; stops marked."""
    names = tuple(current.cases)
    columns = 2 if len(names) > 1 else 1
    rows = max(1, math.ceil(len(names) / columns))
    figure = Figure(figsize=(6.5 * columns, 3.6 * rows + 1.2), facecolor=BACKGROUND)
    FigureCanvasAgg(figure)
    axes = figure.subplots(rows, columns, squeeze=False)
    stops: dict[str, int] = {}
    for position, name in enumerate(names):
        axis = axes[position // columns][position % columns]
        case = current.cases[name]
        if case.pieces:
            axis.semilogy([p.piece_index for p in case.pieces], [p.seconds * 1000 for p in case.pieces], ".", color=CURRENT_COLOR, label="current, median of repeats")
        if baseline is not None and name in baseline.cases and baseline.cases[name].pieces:
            recorded = baseline.cases[name]
            axis.semilogy([p.piece_index for p in recorded.pieces], [p.seconds * 1000 for p in recorded.pieces], "_", color=BASELINE_COLOR, label=f"baseline {baseline.build.commit[:7]}")
        if case.stop_piece_index is not None:
            stops[name] = case.stop_piece_index
            axis.axvline(case.stop_piece_index, color=STOP_COLOR, linewidth=1, linestyle="--", label=f"stop: {case.stop_type}")
        if case.pieces:
            axis.set_title(f"{name.replace('_', ' ')} · {case.completed_queries}/{case.native_pieces} pieces · median {case.median_seconds * 1000:.2f} ms · max {case.max_seconds * 1000:.1f} ms", fontsize=10)
        else:
            axis.set_title(f"{name.replace('_', ' ')} · no completed query", fontsize=10)
        axis.set(xlabel="native piece index", ylabel="query time / ms (log)")
        axis.grid(True, which="both", linewidth=0.3, alpha=0.5)
        axis.legend(fontsize=7, loc="upper right")
    for position in range(len(names), rows * columns):
        axes[position // columns][position % columns].set_visible(False)
    dirty = " (dirty)" if current.build.dirty else ""
    figure.suptitle(f"Native curved-circle query · {current.identity.cpu_brand} · {current.build.commit[:10]}{dirty}", x=0.03, ha="left", fontsize=14)
    figure.text(0.03, 0.005, f"Wall time of one circle_on_piece query at parameter {current.parameter}, median of {current.repeats} repeats per piece; import, plotting and validation excluded. Recorded {current.build.recorded_at}.", fontsize=8)
    figure.subplots_adjust(left=0.075, right=0.99, top=0.9, bottom=0.12, hspace=0.45, wspace=0.2)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=160)
    return FigureSummary(cases=names, stops=stops)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--current", type=pathlib.Path, required=True)
    parser.add_argument("--baseline", type=pathlib.Path)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()
    baseline = read_json(args.baseline) if args.baseline is not None else None
    summary = render_performance_figure(read_json(args.current), baseline, args.output)
    print(f"PLOT {args.output}: {len(summary.cases)} cases, stops {summary.stops}")


if __name__ == "__main__":
    main()
