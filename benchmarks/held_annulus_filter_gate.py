"""The exit gate of unit U3 of the actual-stock refinement: the filtered decision on a real path, against the exact engine and the arrangement.

Per workload, every circle t >= 1, on the circle view [0, t) fed as a stream:

- **Soundness.** Each certified station is within on the exact engine on the view that
  certified it, with the filter's chord bound at or above the exact run. A station no
  view certifies is an exact station.
- **Decision.** The filtered decision equals the exact decision on C_t.
- **Fall-through.** Exact stations, station 0 and k > 0 apart.
- **Cost.** The whole decision through the binding, process CPU, minimum over passes, with
  the arrangement engine's decision in the same run for the load-independent ratio. A loaded
  run claims the ratio, never an absolute budget pass (ruling D2).

Off the production path: nothing here is called by a generator.

CLI (``build/pyrun.sh -m benchmarks.held_annulus_filter_gate ...``):

    WORKLOAD [--passes N] [--ratio-every K] [--no-soundness] [--out DIR]
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from collections.abc import Callable
from dataclasses import asdict
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from benchmarks.held_actual_stock_last16 import Last16CertificateError
from benchmarks.held_actual_stock_last16 import Last16Oracle
from benchmarks.held_actual_stock_last16 import SweepView
from benchmarks.held_actual_stock_oracle import differential_workload
from benchmarks.held_annulus_store_gate import stream_view
from benchmarks.performance_record import REPOSITORY
from compas_cgal import _stock_2
from compas_cgal.pipeline.admissibility import loop_parameters
from tools.actual_stock_profile import NO_GAP_CLOSURE
from tools.actual_stock_profile import ReplayWorkload

SCHEMA = "held-annulus-filter-gate/1"
DEFAULT_OUT = REPOSITORY / "build" / "u2"
US_PER_S = 1e6
UNCERTIFIED = 3
# Chord squared in tool-radius units of a reported span is to_double of exact endpoints, then a sine.
REPORTED_CHORD_AGREEMENT = 1e-9
# Added-cost ceilings of the design, section 4 (seconds per run), and the measured engine calls per run they divide.
CEILING_S = {"figure5_refined": 0.15, "figure8_crossed_skis_boundary": 0.15, "figure8_monstera_boundary_qualified_corners": 0.50}


class FilterUnsoundError(Last16CertificateError):
    """The filter certified a station the exact engine calls exceeded, or its chord bound is below an exact run."""


class FilteredDecisionMismatchError(Last16CertificateError):
    """The filtered decision and the exact decision differ on a circle view."""


@dataclass(frozen=True)
class FilterGateReport:
    """U3's gate on one workload."""

    workload: str
    decisions: int
    stations: int
    certified_last_one: int
    certified_last_sixteen: int
    certified_range: int
    exact_stations: int
    exact_station_zero: int
    refused: int
    soundness_compared: int
    unsound_stations: int
    chord_bound_breaks: int
    decision_mismatches: int
    filtered_decision_us: float
    arrangement_decision_us: float
    rim_exact_decision_us: float
    arrangement_over_filtered: float
    added_s: float
    ceiling_s: float
    load_before: tuple[float, float, float]
    load_after: tuple[float, float, float]

    def require(self) -> None:
        """Raises:
        FilterUnsoundError: a certified station is not within, or a chord bound is low.
        FilteredDecisionMismatchError: a filtered decision differs from the exact one.
        """
        if self.unsound_stations or self.chord_bound_breaks:
            raise FilterUnsoundError(f"{self.workload}: {self.unsound_stations} unsound stations, {self.chord_bound_breaks} chord bound breaks")
        if self.decision_mismatches:
            raise FilteredDecisionMismatchError(f"{self.workload}: {self.decision_mismatches} decisions differ")

    def to_json(self) -> dict[str, object]:
        return {"schema": SCHEMA, **asdict(self)}


def _args(workload: ReplayWorkload, index: int) -> tuple[float, float, float, float, float]:
    c = workload.circles[index]
    return (c.cx, c.cy, c.guide_radius, c.entry_x, c.entry_y)


def soundness(workload: ReplayWorkload, view: _stock_2.AnnulusView2) -> tuple[int, int, int, int]:
    """(stations compared, unsound, chord bound breaks, decision mismatches) over every circle t >= 1."""
    parameters = loop_parameters()
    tool_squared = workload.tool_radius * workload.tool_radius
    compared = unsound = breaks = mismatches = 0
    for index in range(1, len(workload.circles)):
        circle = _args(workload, index)
        args = (*circle, False, parameters, workload.cap_chord_ratio, NO_GAP_CLOSURE)
        report = view.circle_filter_report(0, index, *circle, parameters, workload.cap_chord_ratio)
        ranges = {0: (index - 1, index), 1: (max(0, index - 16), index), 2: (0, index)}
        for tier in {row[0] for row in report} - {UNCERTIFIED}:
            first, end = ranges[tier]
            exact = view.circle_prefixes(first, end, *args)
            for (station_tier, bound), (_total, span, exceeded) in zip(report, exact):
                if station_tier != tier:
                    continue
                compared += 1
                unsound += bool(exceeded)
                breaks += 4.0 * math.sin(span / 2.0) ** 2 * tool_squared > bound + REPORTED_CHORD_AGREEMENT
        mismatches += view.circle_exceeded_filtered(0, index, *args) != view.circle_exceeded(0, index, *args)
    return compared, unsound, breaks, mismatches


def decision_pass_s(workload: ReplayWorkload, view: _stock_2.AnnulusView2) -> float:
    """CPU of one pass of the production decision over every circle, the stream's views."""
    parameters = loop_parameters()
    cap = workload.cap_chord_ratio
    rows = [(index, *_args(workload, index)) for index in range(1, len(workload.circles))]
    decide = view.circle_exceeded_filtered
    start = time.process_time()
    for index, cx, cy, rho, ex, ey in rows:
        decide(0, index, cx, cy, rho, ex, ey, False, parameters, cap, NO_GAP_CLOSURE)
    return time.process_time() - start


def _warm_min(call: Callable[[], object], repeats: int) -> float:
    best = math.inf
    for _ in range(repeats):
        start = time.process_time()
        call()
        best = min(best, time.process_time() - start)
    return best


def reference_decision_us(workload: ReplayWorkload, view: _stock_2.AnnulusView2, every: int) -> tuple[float, float]:
    """The arrangement engine's and the exact rim's warm decision on the last-1 view, sampled, us: the same run's denominators."""
    parameters = loop_parameters()
    oracle = Last16Oracle(workload, sweeps=SweepView.CIRCLES)
    cap = workload.cap_chord_ratio
    arrangement = rim = 0.0
    count = 0
    for index in range(every, len(workload.circles), every):
        common = (*_args(workload, index), False)
        _window, stock = oracle.view(index, 1)
        arrangement += _warm_min(lambda: _stock_2.engagement_circle_exceeded(stock, *common, parameters, workload.tool_radius, cap, NO_GAP_CLOSURE), 3)
        rim += _warm_min(lambda: view.circle_exceeded(index - 1, index, *common, parameters, cap, NO_GAP_CLOSURE), 3)
        count += 1
    return US_PER_S * arrangement / max(count, 1), US_PER_S * rim / max(count, 1)


def filter_gate(workload: ReplayWorkload, passes: int, ratio_every: int, check_soundness: bool = True) -> FilterGateReport:
    load_before = os.getloadavg()
    view = stream_view(workload)
    compared = unsound = breaks = mismatches = 0
    if check_soundness:
        compared, unsound, breaks, mismatches = soundness(workload, view)
    arrangement_us, rim_us = reference_decision_us(workload, view, ratio_every)
    view.reset_filter_counters()
    best = math.inf
    decisions = max(len(workload.circles) - 1, 1)
    for _ in range(passes):
        view.reset_filter_counters()
        best = min(best, decision_pass_s(workload, view))
    counters = view.filter_counters
    filtered_us = US_PER_S * best / decisions
    return FilterGateReport(
        workload=workload.name,
        decisions=counters.decisions,
        stations=counters.stations,
        certified_last_one=counters.certified_last_one,
        certified_last_sixteen=counters.certified_last_sixteen,
        certified_range=counters.certified_range,
        exact_stations=counters.exact_stations,
        exact_station_zero=counters.exact_station_zero,
        refused=counters.refused,
        soundness_compared=compared,
        unsound_stations=unsound,
        chord_bound_breaks=breaks,
        decision_mismatches=mismatches,
        filtered_decision_us=filtered_us,
        arrangement_decision_us=arrangement_us,
        rim_exact_decision_us=rim_us,
        arrangement_over_filtered=arrangement_us / filtered_us if filtered_us > 0.0 else math.inf,
        added_s=best,
        ceiling_s=CEILING_S.get(workload.name, math.inf),
        load_before=load_before,
        load_after=os.getloadavg(),
    )


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("workload")
    parser.add_argument("--passes", type=int, default=9)
    parser.add_argument("--ratio-every", type=int, default=20)
    parser.add_argument("--no-soundness", action="store_true")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args(argv)
    report = filter_gate(differential_workload(args.workload), args.passes, args.ratio_every, not args.no_soundness)
    record = report.to_json()
    print(json.dumps(record, indent=1))
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / f"{args.workload}.filter_gate.json").write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
    report.require()
    return 0


if __name__ == "__main__":
    sys.exit(main())
