"""D1 of the actual-stock refinement: the rim engine ``AnnulusView2`` equals ``Stock2`` on the same circle view, station by station (design section 3).

For each compared circle t the reference side is the existing arrangement engine
on the oracle's circle-only window on the part polygon P
(``benchmarks.held_actual_stock_oracle.view_prefixes`` semantics). The engine reads
the same emission range of the same closed annuli. Three comparisons per station:

- the exact verdict;
- the exact endpoints of every maximal run (``circle_matches_stock``), the strong
  form of "reported max-run equal": it compares points, not their ``to_double``;
- the reported max-run doubles, within ``REPORTED_SPAN_AGREEMENT_RAD`` (reporting
  only: one exact point has several one-root representations, and ``to_double``
  of each differs by ulps).

Off the production path: nothing here is called by a generator.

CLI (``build/pyrun.sh -m benchmarks.held_annulus_view_d1 ...``):

    d1 WORKLOAD [--first N] [--every K] [--last M] [--out DIR]
    time WORKLOAD [--every K] [--repeats N] [--out DIR]
"""

from __future__ import annotations

import argparse
import json
import os
import statistics
import sys
import time
from collections.abc import Callable
from collections.abc import Iterable
from dataclasses import asdict
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

from benchmarks.held_actual_stock_last16 import Last16CertificateError
from benchmarks.held_actual_stock_last16 import Last16Oracle
from benchmarks.held_actual_stock_last16 import SweepView
from benchmarks.held_actual_stock_oracle import differential_workload
from benchmarks.held_material_domain import PartPolygon
from benchmarks.performance_record import REPOSITORY
from compas_cgal import _stock_2
from compas_cgal.pipeline.admissibility import loop_parameters
from tools.actual_stock_profile import NO_GAP_CLOSURE
from tools.actual_stock_profile import ReplayWorkload

SCHEMA = "held-annulus-view-d1/1"
TIME_SCHEMA = "held-annulus-view-time/1"
DEFAULT_OUT = REPOSITORY / "build" / "u1"
MS_PER_S = 1e3
# Runs equal at their exact endpoints report spans differing only by to_double of different one-root representations (radians; about 6e-8 deg).
REPORTED_SPAN_AGREEMENT_RAD = 1e-9
# Parameter of station 0 alone (the entry), for the station-0 decision cost c_x.
STATION_ZERO = [0.0]


class EngineStockMismatchError(Last16CertificateError):
    """D1: the rim engine and the arrangement engine disagree on a station of the same view."""


def view_range(index: int, last: Optional[int]) -> tuple[int, int]:
    """The half-open emission range of circle ``index``'s view: the ``last`` circles before it, or every prior circle. Never the circle itself."""
    first = 0 if last is None else max(0, index - last)
    return first, index


def part_rows(workload: ReplayWorkload) -> np.ndarray:
    """The engine's material domain: the workload's stock rows, through the typed part polygon P."""
    return PartPolygon.from_stock_rows(workload.polygon).stock_rows()


def annulus_rows(workload: ReplayWorkload) -> np.ndarray:
    """Every emitted circle as ``(cx, cy, guide_radius)``, in emission order."""
    return np.array([(c.cx, c.cy, c.guide_radius) for c in workload.circles], dtype=np.float64).reshape(-1, 3)


def engine_view(workload: ReplayWorkload) -> _stock_2.AnnulusView2:
    """The rim engine over the workload's part polygon and all its circles; a decision selects its view by emission range."""
    return _stock_2.AnnulusView2.build(part_rows(workload), annulus_rows(workload), workload.tool_radius)


@dataclass(frozen=True)
class D1Row:
    """One circle's comparison over its stations."""

    circle: int
    stations: int
    verdict_mismatches: int
    exact_run_mismatches: int
    max_reported_delta_rad: float
    exceeded: bool
    engine_ms: float
    reference_ms: float


@dataclass(frozen=True)
class D1Report:
    """D1 over the compared circles of one workload, in one view."""

    workload: str
    last: Optional[int]
    load_before: tuple[float, float, float]
    load_after: tuple[float, float, float]
    rows: tuple[D1Row, ...]

    @property
    def stations(self) -> int:
        return sum(row.stations for row in self.rows)

    @property
    def verdict_mismatches(self) -> int:
        return sum(row.verdict_mismatches for row in self.rows)

    @property
    def exact_run_mismatches(self) -> int:
        return sum(row.exact_run_mismatches for row in self.rows)

    @property
    def max_reported_delta_rad(self) -> float:
        return max((row.max_reported_delta_rad for row in self.rows), default=0.0)

    def require(self) -> None:
        """Raises:
        EngineStockMismatchError: a verdict, a run endpoint or a reported span differs.
        """
        bad = [row.circle for row in self.rows if row.verdict_mismatches or row.exact_run_mismatches or row.max_reported_delta_rad > REPORTED_SPAN_AGREEMENT_RAD]
        if bad:
            raise EngineStockMismatchError(
                f"{self.workload} (last={self.last}): {len(bad)} circles differ, verdicts {self.verdict_mismatches}, exact runs {self.exact_run_mismatches}, "
                f"max reported delta {self.max_reported_delta_rad:.3g} rad, first circles {bad[:20]}"
            )

    def to_json(self) -> dict[str, object]:
        return {
            "schema": SCHEMA,
            "workload": self.workload,
            "last": self.last,
            "circles": len(self.rows),
            "stations": self.stations,
            "verdict_mismatches": self.verdict_mismatches,
            "exact_run_mismatches": self.exact_run_mismatches,
            "max_reported_delta_rad": self.max_reported_delta_rad,
            "engine_ms_sum": sum(row.engine_ms for row in self.rows),
            "reference_ms_sum": sum(row.reference_ms for row in self.rows),
            "load_before": list(self.load_before),
            "load_after": list(self.load_after),
            "rows": [asdict(row) for row in self.rows],
        }


def d1_row(oracle: Last16Oracle, view: _stock_2.AnnulusView2, index: int, last: Optional[int]) -> D1Row:
    """One circle on one view: engine reports and exact agreement against ``Stock2``."""
    workload = oracle.workload
    circle = workload.circles[index]
    parameters = loop_parameters()
    args = (circle.cx, circle.cy, circle.guide_radius, circle.entry_x, circle.entry_y, False, parameters, workload.cap_chord_ratio, NO_GAP_CLOSURE)
    first, end = view_range(index, last)
    start = time.process_time()
    engine = view.circle_prefixes(first, end, *args)
    engine_s = time.process_time() - start
    start = time.process_time()
    _window, stock = oracle.view(index, last)
    reference = workload.prefixes(stock, index)
    reference_s = time.process_time() - start
    agreement = view.circle_matches_stock(stock, first, end, *args)
    return D1Row(
        circle=index,
        stations=len(parameters),
        verdict_mismatches=sum(ours[2] != theirs[2] for ours, theirs in zip(engine, reference)) + sum(not verdict for verdict, _runs in agreement),
        exact_run_mismatches=sum(not runs for _verdict, runs in agreement),
        max_reported_delta_rad=max(abs(ours[1] - theirs[1]) for ours, theirs in zip(engine, reference)),
        exceeded=any(sample[2] for sample in engine),
        engine_ms=MS_PER_S * engine_s,
        reference_ms=MS_PER_S * reference_s,
    )


def d1_report(workload: ReplayWorkload, indices: Iterable[int], last: Optional[int]) -> D1Report:
    """D1 over ``indices`` on the view ``last`` (``None``: every prior circle, the circle view C_t)."""
    load_before = os.getloadavg()
    oracle = Last16Oracle(workload, sweeps=SweepView.CIRCLES)
    view = engine_view(workload)
    rows = tuple(d1_row(oracle, view, index, last) for index in indices)
    return D1Report(workload.name, last, load_before, os.getloadavg(), rows)


# ---- the time gate -----------------------------------------------------------------


def _warm_min(call: Callable[[], bool], repeats: int) -> float:
    best = float("inf")
    for _ in range(repeats):
        start = time.process_time()
        call()
        best = min(best, time.process_time() - start)
    return best


def _distribution_ms(values: list[float]) -> dict[str, float]:
    ordered = sorted(values)
    return {
        "n": len(ordered),
        "p50_ms": MS_PER_S * statistics.median(ordered),
        "p90_ms": MS_PER_S * ordered[int(0.9 * (len(ordered) - 1))],
        "mean_ms": MS_PER_S * statistics.fmean(ordered),
        "max_ms": MS_PER_S * ordered[-1],
    }


def time_decisions(workload: ReplayWorkload, every: int, repeats: int) -> dict[str, object]:
    """Exact decision cost of the engine, process CPU, min of ``repeats`` per circle.

    Sampled circles: the station-0 decision and the 45-station decision on C_t
    (the U-1 cost definition), beside the arrangement engine's on the same circle
    at the same load, warm: its stock is already built and located (the first of
    its repeats pays the point-location build; the minimum discards it). Every circle, once: the 45-station decision on
    last-1 escalating to C_t, the tier the production decision will take, with the
    total against the ceilings of design section 4.
    """
    load_before = os.getloadavg()
    start = time.process_time()
    view = engine_view(workload)
    build_s = time.process_time() - start
    parameters = loop_parameters()
    oracle = Last16Oracle(workload, sweeps=SweepView.CIRCLES)
    station_zero: list[float] = []
    stations: list[float] = []
    reference_station_zero: list[float] = []
    reference_stations: list[float] = []
    for index in range(max(every, 1), len(workload.circles), every):
        c = workload.circles[index]
        first, end = view_range(index, None)
        common = (c.cx, c.cy, c.guide_radius, c.entry_x, c.entry_y, False)
        station_zero.append(_warm_min(lambda: view.circle_exceeded(first, end, *common, STATION_ZERO, workload.cap_chord_ratio, NO_GAP_CLOSURE), repeats))
        stations.append(_warm_min(lambda: view.circle_exceeded(first, end, *common, parameters, workload.cap_chord_ratio, NO_GAP_CLOSURE), repeats))
        _window, stock = oracle.view(index, None)
        reference_station_zero.append(
            _warm_min(lambda: _stock_2.engagement_circle_exceeded(stock, *common, STATION_ZERO, workload.tool_radius, workload.cap_chord_ratio, NO_GAP_CLOSURE), repeats)
        )
        reference_stations.append(
            _warm_min(lambda: _stock_2.engagement_circle_exceeded(stock, *common, parameters, workload.tool_radius, workload.cap_chord_ratio, NO_GAP_CLOSURE), repeats)
        )
    tiers = {"last1": 0, "full": 0, "refused": 0}
    decide_start = time.process_time()
    for index in range(1, len(workload.circles)):
        c = workload.circles[index]
        decision = (c.cx, c.cy, c.guide_radius, c.entry_x, c.entry_y, False, parameters, workload.cap_chord_ratio, NO_GAP_CLOSURE)
        if not view.circle_exceeded(index - 1, index, *decision):
            tiers["last1"] += 1
        elif not view.circle_exceeded(0, index, *decision):
            tiers["full"] += 1
        else:
            tiers["refused"] += 1
    decide_s = time.process_time() - decide_start
    return {
        "schema": TIME_SCHEMA,
        "workload": workload.name,
        "every": every,
        "repeats": repeats,
        "circles": len(workload.circles),
        "build_ms": MS_PER_S * build_s,
        "station_zero_on_full_view": _distribution_ms(station_zero),
        "stations45_on_full_view": _distribution_ms(stations),
        "reference_station_zero_warm": _distribution_ms(reference_station_zero),
        "reference_stations45_warm": _distribution_ms(reference_stations),
        "tiered_all_circles_s": decide_s,
        "tiered_per_decision_us": 1e6 * decide_s / max(len(workload.circles) - 1, 1),
        "tiers": tiers,
        "load_before": list(load_before),
        "load_after": list(os.getloadavg()),
    }


def _write(out: Path, stem: str, record: dict[str, object]) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{stem}.json"
    path.write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
    return path


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    d1 = commands.add_parser("d1")
    d1.add_argument("workload")
    d1.add_argument("--first", type=int, default=None, help="only the first N circles")
    d1.add_argument("--every", type=int, default=1)
    d1.add_argument("--last", type=int, default=None, help="view: the last M circles (default: every prior circle)")
    d1.add_argument("--out", type=Path, default=DEFAULT_OUT)
    timing = commands.add_parser("time")
    timing.add_argument("workload")
    timing.add_argument("--every", type=int, default=20)
    timing.add_argument("--repeats", type=int, default=5)
    timing.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args(argv)
    workload = differential_workload(args.workload)
    if args.command == "time":
        record = time_decisions(workload, args.every, args.repeats)
        print(json.dumps(record, indent=1))
        _write(args.out, f"{args.workload}.time", record)
        return 0
    count = len(workload.circles) if args.first is None else min(args.first, len(workload.circles))
    report = d1_report(workload, range(0, count, args.every), args.last)
    record = report.to_json()
    summary = {key: value for key, value in record.items() if key != "rows"}
    print(json.dumps(summary, indent=1))
    _write(args.out, f"{args.workload}.d1.first{args.first}.every{args.every}.last{args.last}", record)
    report.require()
    return 0


if __name__ == "__main__":
    sys.exit(main())
