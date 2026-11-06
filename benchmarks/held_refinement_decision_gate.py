"""The exit gate of unit U4 of the actual-stock refinement: production implies certificate, and the decision really ran.

Per workload, on the path its generator emits:

- **D2.** Every emitted circle t >= 1 is ``CERTIFIED_LAST16`` or ``CERTIFIED_FULL_WINDOW``
  on the actual stock (and, with ``circle_view``, on the circle view the decision read),
  and circle 0 is ``REFUSED``.
- **Call count.** The engine decided at least final circles - 1 targets; each refined
  circle t >= 1 has an accepted decision on its own doubles at its own position; the
  engine's own count equals the footprints; the first refinement's site reached the
  engine, and the seam site at least once per seam circle.

Off the production path: the harness wraps the generator's refinement entry points to
read what they return, and calls nothing the generator does not.

CLI (``build/pyrun.sh -m benchmarks.held_refinement_decision_gate ...``):

    WORKLOAD [--every K] [--processes N] [--circle-view] [--out DIR]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import asdict
from dataclasses import dataclass
from functools import wraps
from pathlib import Path
from typing import Any
from typing import Optional

from benchmarks import held_toolpath_benchmark as toolpath_benchmark
from benchmarks.errors import BenchmarkError
from benchmarks.held_actual_stock_decision import ENTRY_POSITION
from benchmarks.held_actual_stock_decision import read_disk
from benchmarks.held_actual_stock_last16 import SweepView
from benchmarks.held_actual_stock_last16 import Verdict
from benchmarks.held_actual_stock_oracle import D2Report
from benchmarks.held_actual_stock_oracle import production_implies_certificate
from benchmarks.held_actual_stock_oracle import workload_on_part
from benchmarks.held_figure5_engagement_refinement import TracedEngagementRefinement
from benchmarks.held_figure5_engagement_trace import DecisionSite
from benchmarks.performance_record import REPOSITORY

SCHEMA = "held-refinement-decision-gate/1"
DEFAULT_OUT = REPOSITORY / "build" / "u4"


class RefinementDecisionGateError(BenchmarkError):
    """Base class of the U4 gate's failures."""


class RefinementNotCapturedError(RefinementDecisionGateError):
    """The generator did not run exactly one refinement."""


class DecisionCallCountError(RefinementDecisionGateError):
    """The engine was not reached as often, or at the positions, the wiring requires."""


@contextmanager
def _captured_refinements() -> Iterator[list[TracedEngagementRefinement]]:
    """Record what the generators' refinement entry points return, unchanged."""
    captured: list[TracedEngagementRefinement] = []
    figure5 = getattr(toolpath_benchmark, "refine_figure5_engagement")
    boundary = getattr(toolpath_benchmark, "refine_boundary_ordered_path")

    @wraps(figure5)
    def on_figure5(*args: Any, **kwargs: Any) -> TracedEngagementRefinement:
        result: TracedEngagementRefinement = figure5(*args, **kwargs)
        captured.append(result)
        return result

    @wraps(boundary)
    def on_boundary(*args: Any, **kwargs: Any) -> Any:
        result = boundary(*args, **kwargs)
        captured.append(result[1])
        return result

    setattr(toolpath_benchmark, "refine_figure5_engagement", on_figure5)
    setattr(toolpath_benchmark, "refine_boundary_ordered_path", on_boundary)
    try:
        yield captured
    finally:
        setattr(toolpath_benchmark, "refine_figure5_engagement", figure5)
        setattr(toolpath_benchmark, "refine_boundary_ordered_path", boundary)


@dataclass(frozen=True)
class CallCount:
    """Where the stock decisions of one refinement were taken, and what the engine says it did."""

    refined_circles: int
    final_circles: int
    seam_circles: int
    decisions: int
    advance: int
    close_seam: int
    continuation: int
    undecided_circles: tuple[int, ...]
    engine_decisions: int
    engine_exact_decisions: int
    engine_stations: int
    engine_certified_last_one: int
    engine_certified_last_sixteen: int
    engine_certified_range: int
    engine_exact_stations: int
    engine_exact_station_zero: int
    engine_refused: int

    def require(self) -> None:
        """Raises:
        DecisionCallCountError: the wiring breaches any of the counts of the module docstring.
        """
        if self.decisions < self.final_circles - 1:
            raise DecisionCallCountError(f"{self.decisions} decisions for {self.final_circles} final circles")
        if self.engine_decisions != self.decisions:
            raise DecisionCallCountError(f"the engine counted {self.engine_decisions} decisions, the footprints {self.decisions}")
        if self.undecided_circles:
            raise DecisionCallCountError(f"{len(self.undecided_circles)} refined circles have no accepted decision of their own: {list(self.undecided_circles[:20])}")
        if self.advance == 0:
            raise DecisionCallCountError("no decision at the first refinement's site")
        if self.close_seam < self.seam_circles:
            raise DecisionCallCountError(f"{self.seam_circles} seam circles and {self.close_seam} seam decisions")


def call_count(refinement: TracedEngagementRefinement, final_circles: int) -> CallCount:
    """The refinement's counts, from its trace and the engine's own counters."""
    trace = refinement.trace
    footprints = [footprint for footprint in trace.footprints if footprint.read is not None]
    accepted = {(footprint.position, footprint.read) for footprint in footprints}
    undecided = tuple(t for t, circle in enumerate(refinement.circles) if t != ENTRY_POSITION and (t, read_disk(circle)) not in accepted)
    counters = trace.filter_counters
    if counters is None:
        raise RefinementNotCapturedError("the refinement's trace carries no engine counters")
    return CallCount(
        refined_circles=len(refinement.circles),
        final_circles=final_circles,
        seam_circles=len(refinement.circles) - 1 - refinement.original_indices[-1],
        decisions=len(footprints),
        advance=trace.decisions_at(DecisionSite.ADVANCE),
        close_seam=trace.decisions_at(DecisionSite.CLOSE_SEAM),
        continuation=trace.decisions_at(DecisionSite.CONTINUATION),
        undecided_circles=undecided,
        engine_decisions=counters.decisions,
        engine_exact_decisions=counters.exact_decisions,
        engine_stations=counters.stations,
        engine_certified_last_one=counters.certified_last_one,
        engine_certified_last_sixteen=counters.certified_last_sixteen,
        engine_certified_range=counters.certified_range,
        engine_exact_stations=counters.exact_stations,
        engine_exact_station_zero=counters.exact_station_zero,
        engine_refused=counters.refused,
    )


@dataclass(frozen=True)
class RefinementDecisionReport:
    """U4's gate on one workload."""

    workload: str
    every: int
    generation_cpu_s: float
    counts: CallCount
    d2_actual: D2Report
    d2_circle_view: Optional[D2Report]

    def require(self) -> None:
        """Raises:
        DecisionCallCountError: see ``CallCount.require``.
        EntryCircleNotRefusedError, CertificateViolationError: see ``D2Report.require``.
        """
        self.counts.require()
        self.d2_actual.require()
        if self.d2_circle_view is not None:
            self.d2_circle_view.require()

    def to_json(self) -> dict[str, object]:
        def plain(value: object) -> object:
            if isinstance(value, (Verdict, SweepView)):
                return value.value
            raise TypeError(f"not JSON serializable: {type(value).__name__}")

        record: dict[str, object] = json.loads(json.dumps({"schema": SCHEMA, **asdict(self)}, default=plain))
        return record


def generate_with_refinement(workload: str) -> tuple[toolpath_benchmark.GeneratedToolpath, TracedEngagementRefinement, float]:
    """Run the generator: its path, the one refinement behind it, its process CPU.

    Raises:
        RefinementNotCapturedError: the generator did not run exactly one refinement.
    """
    with _captured_refinements() as captured:
        start = time.process_time()
        path = toolpath_benchmark.generate(workload)
        cpu = time.process_time() - start
    if len(captured) != 1:
        raise RefinementNotCapturedError(f"{workload}: {len(captured)} refinements")
    return path, captured[0], cpu


def decision_gate(workload: str, every: int, processes: int, circle_view: bool = False) -> RefinementDecisionReport:
    """The whole gate on ``workload``: ``every``-th circle under D2, always circle 0."""
    path, refinement, cpu = generate_with_refinement(workload)
    replay = workload_on_part(path, cpu)
    indices = range(0, len(replay.circles), every)
    return RefinementDecisionReport(
        workload=workload,
        every=every,
        generation_cpu_s=cpu,
        counts=call_count(refinement, len(path.circles)),
        d2_actual=production_implies_certificate(replay, indices, processes, SweepView.ACTUAL),
        d2_circle_view=production_implies_certificate(replay, indices, processes, SweepView.CIRCLES) if circle_view else None,
    )


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workload")
    parser.add_argument("--every", type=int, default=1)
    parser.add_argument("--processes", type=int, default=2)
    parser.add_argument("--circle-view", action="store_true")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args(argv)
    report = decision_gate(args.workload, args.every, args.processes, args.circle_view)
    record = report.to_json()
    print(json.dumps(record, indent=1))
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / f"{args.workload}.decision_gate.json").write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
    report.require()
    return 0


if __name__ == "__main__":
    sys.exit(main())
