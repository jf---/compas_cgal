"""Tier 2 of the connector-engagement check: exact station verdicts on the replayed actual stock.

The oracle and the refutation tool for ``benchmarks/held_connector_engagement.py``,
never part of generation (``docs/held_connector_engagement.md``, Tier 2).
Stock is replayed in emission order, circle ``i`` then connector ``i``, with
the native exact stock. Each connector piece is decided at exact stations
(``engagement_line_prefixes``: stock minus the piece's own swept prefix)
before it is subtracted. Two stocks bracket the truth: ``Stock_UNDER``
(``subtract_line_sweep_local``, contains the true stock) and ``Stock_OVER``
(``subtract_line_sweep_over``, contained in it). Max-run engagement is
monotone in the stock, so a clear station on UNDER is truly clear and an
exceeded station on OVER is a true exceedance. Nothing between stations is
claimed.

Windowing is exact: engagement at p reads only S ∩ D(p, r). A window disk
strictly containing every station disk of the connector, minus every prior
sweep that can meet it, equals the full replay inside the window.

Two modes:

- ``--windowed`` (default; opt-in pixi task, minutes): every connector the
  certificate leaves without a tier-1/1b certificate, on every workload;
  writes ``benchmarks/results/held_connector_engagement_tier2.{json,md}``.
- ``--full WORKLOAD``: the whole sequential replay of one workload (every
  connector; about 20 minutes on Figure 5). One-off evidence; writes
  ``build/held-connector-engagement/<workload>.full.json``.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from collections.abc import Iterable
from collections.abc import Sequence
from dataclasses import asdict
from dataclasses import dataclass
from enum import Enum

import numpy as np
import numpy.typing as npt

from benchmarks.held_connector_engagement import ConnectorEngagementAssessment
from benchmarks.held_connector_engagement import ConnectorEngagementLemma
from benchmarks.held_connector_engagement import LemmaCircle
from benchmarks.held_connector_engagement import moving_pieces
from benchmarks.held_figure5_boundary_path import Figure5BoundaryTransition
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.performance_record import REPOSITORY
from benchmarks.performance_record import build_identity
from benchmarks.performance_record import machine_identity
from compas_cgal import _stock_2
from compas_cgal.motion import EngagementCap
from compas_cgal.pipeline.admissibility import motion_parameters
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

FloatArray = npt.NDArray[np.float64]
Segment = tuple[float, float, float, float]

# The pipeline's own gap-closure ratio for these verdicts: none.
NO_GAP_CLOSURE = 0.0
# Window margin beyond the connector's bounding box, tool radii: any positive
# margin keeps every station's closed tool disk inside the open window.
WINDOW_MARGIN_TOOL_RADII = 2.0
# Extra reach of the prior-sweep filter, tool radii, so an OVER line sweep
# (at most 1e-4 r beyond its capsule) is never dropped; extra sweeps outside
# the window change nothing.
FILTER_SLACK_TOOL_RADII = 1.0
# The window's outer clearing ring reaches this many region diameters (plus
# ten tool radii) out, so it removes everything outside the window, mm-free.
RING_REACH_DIAMETERS = 2.0
RING_REACH_TOOL_RADII = 10.0

OUTPUT_MD = REPOSITORY / "benchmarks/results/held_connector_engagement_tier2.md"
OUTPUT_JSON = REPOSITORY / "benchmarks/results/held_connector_engagement_tier2.json"
FULL_DIR = REPOSITORY / "build/held-connector-engagement"


class StationVerdict(Enum):
    """A connector's station verdict under the UNDER/OVER bracket."""

    CLEAR = "clear"  # no exceeded station on Stock_UNDER: truly clear at every station
    VIOLATED = "violated"  # an exceeded station on Stock_OVER: a true exceedance
    UNRESOLVED = "unresolved"  # exceeded on UNDER only: inside the bracket


@dataclass(frozen=True)
class ConnectorStations:
    """Exact station verdicts of one connector on both bracketing stocks."""

    connector: int
    pieces: int
    stations: int
    max_run_deg_under: float
    exceeded_under: int
    max_run_deg_over: float
    exceeded_over: int

    @property
    def verdict(self) -> StationVerdict:
        if self.exceeded_over:
            return StationVerdict.VIOLATED
        return StationVerdict.CLEAR if not self.exceeded_under else StationVerdict.UNRESOLVED


def _polygon(region: Sequence[Point2[WorldXY]]) -> FloatArray:
    return np.array([(float(p.x), float(p.y), 0.0) for p in region], dtype=np.float64)


def _segments(transition: Figure5BoundaryTransition) -> list[Segment]:
    return [(float(a[0]), float(a[1]), float(b[0]), float(b[1])) for a, b in moving_pieces(transition)]


class _Tally:
    def __init__(self) -> None:
        self.pieces = self.stations = self.exceeded_under = self.exceeded_over = 0
        self.max_under = self.max_over = 0.0

    def piece(self, under: _stock_2.Stock2, over: _stock_2.Stock2, segment: Segment, tool: float, ratio: float, parameters: Sequence[float]) -> None:
        x0, y0, x1, y1 = segment
        rows_under = _stock_2.engagement_line_prefixes(under, x0, y0, x1, y1, parameters, tool, ratio, NO_GAP_CLOSURE)
        rows_over = _stock_2.engagement_line_prefixes(over, x0, y0, x1, y1, parameters, tool, ratio, NO_GAP_CLOSURE)
        self.pieces += 1
        self.stations += len(rows_under)
        self.max_under = max(self.max_under, max(math.degrees(row[1]) for row in rows_under))
        self.max_over = max(self.max_over, max(math.degrees(row[1]) for row in rows_over))
        self.exceeded_under += sum(1 for row in rows_under if row[2])
        self.exceeded_over += sum(1 for row in rows_over if row[2])
        under.subtract_line_sweep_local(x0, y0, x1, y1, tool)
        over.subtract_line_sweep_over(x0, y0, x1, y1, tool)

    def row(self, connector: int) -> ConnectorStations:
        return ConnectorStations(connector, self.pieces, self.stations, self.max_under, self.exceeded_under, self.max_over, self.exceeded_over)


def replay_full(
    region: Sequence[Point2[WorldXY]],
    circles: tuple[Figure5CounterclockwiseCircle, ...],
    transitions: tuple[Figure5BoundaryTransition, ...],
    tool: ToolRadius,
    cap: EngagementCap,
    parameters: Sequence[float],
) -> tuple[ConnectorStations, ...]:
    """Sequential replay of the whole path; one row per connector."""
    r = float(tool.value)
    under = _stock_2.Stock2(_polygon(region), [])
    over = _stock_2.Stock2(_polygon(region), [])
    rows = []
    for index, transition in enumerate(transitions):
        circle = circles[index]
        cx, cy, rho = float(circle.center.x), float(circle.center.y), float(circle.radius.value)
        under.subtract_circle_sweep_local(cx, cy, rho, r)
        # Global, not local: local surgery after an OVER line sweep is a known
        # native defect (LocalDepletionEscapedError or a segfault; repro in
        # the implementation report). The removed region is the same annulus.
        over.subtract_circle_sweep(cx, cy, rho, r)
        tally = _Tally()
        for segment in _segments(transition):
            tally.piece(under, over, segment, r, cap.chord_ratio, parameters)
        rows.append(tally.row(index))
    return tuple(rows)


def replay_windowed(
    region: Sequence[Point2[WorldXY]],
    circles: tuple[Figure5CounterclockwiseCircle, ...],
    transitions: tuple[Figure5BoundaryTransition, ...],
    tool: ToolRadius,
    cap: EngagementCap,
    parameters: Sequence[float],
    connectors: Iterable[int],
) -> tuple[ConnectorStations, ...]:
    """Per-connector windowed replay; equal to :func:`replay_full` on every selected row."""
    r = float(tool.value)
    polygon = _polygon(region)
    xy = polygon[:, :2]
    ring_outer = RING_REACH_DIAMETERS * float(np.hypot(*(xy.max(axis=0) - xy.min(axis=0)))) + RING_REACH_TOOL_RADII * r
    centres = np.array([(float(c.center.x), float(c.center.y)) for c in circles], dtype=np.float64)
    outer = np.array([float(c.radius.value) + r for c in circles], dtype=np.float64)
    owned = [_segments(transition) for transition in transitions]
    rows = []
    for index in connectors:
        own = owned[index]
        if not own:
            rows.append(_Tally().row(index))
            continue
        points = np.array([(s[0], s[1]) for s in own] + [(own[-1][2], own[-1][3])], dtype=np.float64)
        low, high = points.min(axis=0), points.max(axis=0)
        cx, cy = (low + high) / 2.0
        radius = float(np.hypot(*(high - low))) / 2.0 + WINDOW_MARGIN_TOOL_RADII * r
        reach = radius + FILTER_SLACK_TOOL_RADII * r
        circle_hits = np.nonzero(np.hypot(centres[: index + 1, 0] - cx, centres[: index + 1, 1] - cy) <= reach + outer[: index + 1])[0]
        under = _stock_2.Stock2(polygon, [])
        under.subtract_annulus(float(cx), float(cy), radius, ring_outer)
        for k in circle_hits:
            c = circles[int(k)]
            under.subtract_circle_sweep_local(float(c.center.x), float(c.center.y), float(c.radius.value), r)
        # Every local circle surgery precedes the first OVER line sweep (see replay_full).
        over = under.clone()
        for prior in owned[:index]:
            for x0, y0, x1, y1 in prior:
                nx = min(max(float(cx), min(x0, x1) - r), max(x0, x1) + r)
                ny = min(max(float(cy), min(y0, y1) - r), max(y0, y1) + r)
                if math.hypot(nx - cx, ny - cy) <= reach:
                    under.subtract_line_sweep_local(x0, y0, x1, y1, r)
                    over.subtract_line_sweep_over(x0, y0, x1, y1, r)
        tally = _Tally()
        for segment in own:
            tally.piece(under, over, segment, r, cap.chord_ratio, parameters)
        rows.append(tally.row(index))
    return tuple(rows)


@dataclass(frozen=True)
class Tier2Row:
    """One uncertified connector of one workload, measured at its stations."""

    workload: str
    stations: ConnectorStations
    tier1_bound_deg: float


def _tier1_bound_deg(circle: Figure5CounterclockwiseCircle, transition: Figure5BoundaryTransition, lemma: ConnectorEngagementLemma) -> float:
    exact = LemmaCircle.from_circle(circle)
    return max(lemma.bound_deg(a, b, exact) for a, b in moving_pieces(transition))


def windowed_left(workloads: Sequence[str]) -> list[Tier2Row]:
    """Tier 2 on every connector tiers 1 and 1b leave without a certificate."""
    from benchmarks.held_toolpath_benchmark import generate

    rows = []
    for workload in workloads:
        path = generate(workload)
        cap = EngagementCap.build(math.radians(float(path.case.tea_cap)))
        tool = path.case.tool_radius
        assessment = ConnectorEngagementAssessment.build(path.circles, path.connectors, tool, cap)
        lemma = ConnectorEngagementLemma.build(tool, cap)
        stations = replay_windowed(path.case.projection.points, path.circles, path.transitions, tool, cap, motion_parameters(), assessment.refused)
        for row in stations:
            rows.append(Tier2Row(workload, row, _tier1_bound_deg(path.circles[row.connector], path.transitions[row.connector], lemma)))
    return rows


def _markdown(rows: list[Tier2Row], workloads: Sequence[str], seconds: float) -> str:
    build, machine = build_identity(), machine_identity()
    lines = [
        "<!-- Generated by `pixi run held-connector-engagement-tier2` (tools/held_connector_engagement_tier2.py); do not edit. -->",
        "",
        f"Commit `{build.commit[:10]}`{' (dirty)' if build.dirty else ''} · {machine.cpu_brand} · recorded {build.recorded_at} · {seconds:.0f} s.",
        "",
        "| workload | uncertified connectors | stations over the cap (UNDER / OVER) | station max-run (°) | measured / tier-1 bound |",
        "| --- | ---: | ---: | --- | --- |",
    ]
    for workload in workloads:
        mine = [row for row in rows if row.workload == workload]
        if not mine:
            lines.append(f"| `{workload}` | 0 | — | — | — |")
            continue
        runs = sorted(row.stations.max_run_deg_under for row in mine)
        ratios = [row.stations.max_run_deg_under / row.tier1_bound_deg for row in mine]
        under = sum(row.stations.exceeded_under for row in mine)
        over = sum(row.stations.exceeded_over for row in mine)
        spread = f"{runs[0]:.2f}" if len(runs) == 1 else f"{runs[0]:.2f}–{runs[-1]:.2f} (median {runs[len(runs) // 2]:.2f})"
        lines.append(f"| `{workload}` | {len(mine)} | {under} / {over} | {spread} | {f'{min(ratios):.2f}' if len(ratios) == 1 else f'{min(ratios):.2f}–{max(ratios):.2f}'} |")
    worst = max(rows, key=lambda row: row.stations.max_run_deg_under, default=None)
    if worst is not None:
        lines += ["", f"Closest call: `{worst.workload}` connector {worst.stations.connector}, {worst.stations.max_run_deg_under:.2f}° at a station."]
    return "\n".join(lines) + "\n"


def _full(workload: str) -> None:
    from benchmarks.held_toolpath_benchmark import generate

    path = generate(workload)
    cap = EngagementCap.build(math.radians(float(path.case.tea_cap)))
    tool = path.case.tool_radius
    lemma = ConnectorEngagementLemma.build(tool, cap)
    started = time.perf_counter()
    rows = replay_full(path.case.projection.points, path.circles, path.transitions, tool, cap, motion_parameters())
    seconds = time.perf_counter() - started
    moving = [row for row in rows if row.pieces]
    ratios = sorted(row.max_run_deg_under / bound for row in moving if (bound := _tier1_bound_deg(path.circles[row.connector], path.transitions[row.connector], lemma)) > 0)
    summary = {
        "workload": workload,
        "kind": "ONE-OFF EVIDENCE: exact station verdicts on the full sequential replay; no between-station claim",
        "connectors": len(moving),
        "stations": sum(row.stations for row in moving),
        "exceeded_under": sum(row.exceeded_under for row in moving),
        "exceeded_over": sum(row.exceeded_over for row in moving),
        "max_station_max_run_deg": max(row.max_run_deg_under for row in moving),
        "seam_max_run_deg": rows[-1].max_run_deg_under if len(rows) == len(path.circles) else None,
        "max_run_over_tier1_bound": {"median": ratios[len(ratios) // 2], "max": ratios[-1]},
        "replay_seconds": seconds,
    }
    FULL_DIR.mkdir(parents=True, exist_ok=True)
    out = FULL_DIR / f"{workload}.full.json"
    out.write_text(json.dumps({"summary": summary, "rows": [asdict(row) for row in rows]}, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def main(argv: Sequence[str] | None = None) -> None:
    from benchmarks.held_toolpath_benchmark import GENERATORS

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    parser.add_argument("--full", metavar="WORKLOAD", default=None, help="one-off: the full sequential replay of one workload")
    args = parser.parse_args(argv)
    if args.full is not None:
        _full(args.full)
        return
    workloads = tuple(GENERATORS)
    started = time.perf_counter()
    rows = windowed_left(workloads)
    seconds = time.perf_counter() - started
    OUTPUT_JSON.write_text(json.dumps([{"workload": row.workload, "tier1_bound_deg": row.tier1_bound_deg, **asdict(row.stations)} for row in rows], indent=2) + "\n")
    OUTPUT_MD.write_text(_markdown(rows, workloads, seconds))
    print(OUTPUT_MD.read_text())


if __name__ == "__main__":
    main()
