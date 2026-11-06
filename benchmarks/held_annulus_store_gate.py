"""The exit gate of unit U2 of the actual-stock refinement: ``AnnulusStore2`` behind ``AnnulusView2``, on a real emission stream.

Three statements, per workload:

- **Selection.** For every decision of the stream, the store returns every annulus of
  the view whose box meets the decision's read box (brute force over the same boxes).
  A box miss is an exact miss, so 0 misses is the superset the engine relies on.
- **Query cost.** One query per decision at full store size, process CPU, against
  ``QUERY_CEILING_US`` (design section 5). Stale references from a truncation are part
  of the measured cost. A loaded run states its load and claims nothing absolute.
- **Nested views.** On the sampled circles, last-1 >= last-16 >= the full view, station
  by station: verdicts ordered, max runs non-increasing (0 bound breaks).

Off the production path: nothing here is called by a generator.

CLI (``build/pyrun.sh -m benchmarks.held_annulus_store_gate ...``):

    WORKLOAD [--every K] [--repeats N] [--out DIR]
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from dataclasses import asdict
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
from numpy.typing import NDArray

from benchmarks.held_actual_stock_last16 import Last16CertificateError
from benchmarks.held_actual_stock_oracle import differential_workload
from benchmarks.held_annulus_view_d1 import part_rows
from benchmarks.performance_record import REPOSITORY
from compas_cgal import _stock_2
from compas_cgal.pipeline.admissibility import loop_parameters
from tools.actual_stock_profile import NO_GAP_CLOSURE
from tools.actual_stock_profile import CircleSweep
from tools.actual_stock_profile import ReplayWorkload

SCHEMA = "held-annulus-store-gate/1"
DEFAULT_OUT = REPOSITORY / "build" / "u2"
US_PER_S = 1e6
# Design section 5, U2's exit gate (us of CPU per decision's query, at Monstera size).
QUERY_CEILING_US = 5.0
# The second tier's size (design section 2.5).
TIER_LAST = 16
# Max-run reports of nested views may differ only by to_double of different one-root representations (radians).
REPORTED_SPAN_AGREEMENT_RAD = 1e-9


class StoreCandidateMissError(Last16CertificateError):
    """The store dropped an annulus whose box meets a decision's read box."""


class ViewBoundBreakError(Last16CertificateError):
    """A smaller view reported more engagement than a larger one: the nesting the tiers rest on is broken."""


def read_box(circle: CircleSweep, tool_radius: float) -> tuple[float, float, float, float]:
    """The box of the disk every rim of ``circle`` lies in: B(c, max(rho, |e - c|) + r), outward by one ulp per side."""
    reach = math.nextafter(max(circle.guide_radius, math.hypot(circle.entry_x - circle.cx, circle.entry_y - circle.cy)) + tool_radius, math.inf)
    return (
        math.nextafter(circle.cx - reach, -math.inf),
        math.nextafter(circle.cy - reach, -math.inf),
        math.nextafter(circle.cx + reach, math.inf),
        math.nextafter(circle.cy + reach, math.inf),
    )


def annulus_box(circle: CircleSweep, tool_radius: float) -> tuple[float, float, float, float]:
    """The annulus's box, wider than any rounding of ``c +- (rho + r)``: the brute-force side of the selection check."""
    outer = math.nextafter(circle.guide_radius + tool_radius, math.inf)
    return (
        math.nextafter(circle.cx - outer, -math.inf),
        math.nextafter(circle.cy - outer, -math.inf),
        math.nextafter(circle.cx + outer, math.inf),
        math.nextafter(circle.cy + outer, math.inf),
    )


def stream_view(workload: ReplayWorkload, count: Optional[int] = None) -> _stock_2.AnnulusView2:
    """The engine opened empty and fed in emission order, as refinement feeds it."""
    view = _stock_2.AnnulusView2.open(part_rows(workload), workload.tool_radius)
    for circle in workload.circles if count is None else workload.circles[:count]:
        view.append(circle.cx, circle.cy, circle.guide_radius)
    return view


def selection_misses(workload: ReplayWorkload) -> tuple[int, int]:
    """(misses, candidates) over every decision of the stream: the store against brute force on the same boxes."""
    view = _stock_2.AnnulusView2.open(part_rows(workload), workload.tool_radius)
    tool = workload.tool_radius
    boxes = np.array([annulus_box(circle, tool) for circle in workload.circles], dtype=np.float64)
    misses = 0
    candidates = 0
    for index, circle in enumerate(workload.circles):
        query = read_box(circle, tool)
        found = view.candidates(0, index, *query)
        candidates += len(found)
        prior = boxes[:index]
        meets = (prior[:, 2] >= query[0]) & (prior[:, 0] <= query[2]) & (prior[:, 3] >= query[1]) & (prior[:, 1] <= query[3])
        misses += int(np.count_nonzero(meets)) - len(set(found) & set(np.flatnonzero(meets).tolist()))
        view.append(circle.cx, circle.cy, circle.guide_radius)
    return misses, candidates


def _query_boxes(workload: ReplayWorkload) -> NDArray[np.float64]:
    return np.array([read_box(circle, workload.tool_radius) for circle in workload.circles], dtype=np.float64)


def query_us(view: _stock_2.AnnulusView2, boxes: NDArray[np.float64], end: int, repeats: int) -> float:
    """Minimum over ``repeats`` of the CPU per query, over all of ``boxes`` on the view [0, end)."""
    best = math.inf
    for _ in range(repeats):
        start = time.process_time()
        view.count_candidates(0, end, boxes)
        best = min(best, time.process_time() - start)
    return US_PER_S * best / max(len(boxes), 1)


def nested_breaks(workload: ReplayWorkload, every: int) -> tuple[int, int, list[int]]:
    """(stations compared, bound breaks, first circles with a break) on last-1, last-16 and the full view."""
    view = stream_view(workload)
    parameters = loop_parameters()
    stations = 0
    bad: list[int] = []
    for index in range(max(every, 1), len(workload.circles), every):
        c = workload.circles[index]
        args = (c.cx, c.cy, c.guide_radius, c.entry_x, c.entry_y, False, parameters, workload.cap_chord_ratio, NO_GAP_CLOSURE)
        one = view.circle_prefixes(index - 1, index, *args)
        sixteen = view.circle_prefixes(max(0, index - TIER_LAST), index, *args)
        full = view.circle_prefixes(0, index, *args)
        stations += len(full)
        broken = False
        for wide, middle, narrow in zip(one, sixteen, full):
            if narrow[2] and not middle[2] or middle[2] and not wide[2]:
                broken = True
            if middle[1] > wide[1] + REPORTED_SPAN_AGREEMENT_RAD or narrow[1] > middle[1] + REPORTED_SPAN_AGREEMENT_RAD:
                broken = True
        if broken:
            bad.append(index)
    return stations, len(bad), bad[:20]


@dataclass(frozen=True)
class StoreGateReport:
    """U2's gate on one workload."""

    workload: str
    circles: int
    candidates: int
    candidate_misses: int
    append_us: float
    query_us_full_store: float
    query_us_after_truncation: float
    truncate_us: float
    stale_references: int
    nested_stations: int
    nested_breaks: int
    nested_first_breaks: tuple[int, ...]
    load_before: tuple[float, float, float]
    load_after: tuple[float, float, float]

    def require(self) -> None:
        """Raises:
        StoreCandidateMissError: an annulus that meets a read box was not returned.
        ViewBoundBreakError: a nested view reported more engagement than a larger one.
        """
        if self.candidate_misses:
            raise StoreCandidateMissError(f"{self.workload}: {self.candidate_misses} missed candidates")
        if self.nested_breaks:
            raise ViewBoundBreakError(f"{self.workload}: {self.nested_breaks} circles break nesting, first {list(self.nested_first_breaks)}")

    def to_json(self) -> dict[str, object]:
        return {"schema": SCHEMA, "query_ceiling_us": QUERY_CEILING_US, **asdict(self)}


def store_gate(workload: ReplayWorkload, every: int, repeats: int) -> StoreGateReport:
    """U2's gate: selection on every decision, query cost at full size and through a truncation, nesting on every ``every``-th circle."""
    load_before = os.getloadavg()
    misses, candidates = selection_misses(workload)
    best_append = math.inf
    view = stream_view(workload)
    for _ in range(repeats):
        start = time.process_time()
        view = stream_view(workload)
        best_append = min(best_append, time.process_time() - start)
    circles = len(workload.circles)
    boxes = _query_boxes(workload)
    full_size = query_us(view, boxes, circles, repeats)
    # A splice drops a tail smaller than the live store, so the stale references stay in the grid.
    keep = circles - circles // 4
    start = time.process_time()
    view.truncate(keep)
    truncate_s = time.process_time() - start
    stale = view.stale_references
    after = query_us(view, boxes, keep, repeats)
    stations, breaks, first = nested_breaks(workload, every)
    return StoreGateReport(
        workload=workload.name,
        circles=circles,
        candidates=candidates,
        candidate_misses=misses,
        append_us=US_PER_S * best_append / max(circles, 1),
        query_us_full_store=full_size,
        query_us_after_truncation=after,
        truncate_us=US_PER_S * truncate_s,
        stale_references=stale,
        nested_stations=stations,
        nested_breaks=breaks,
        nested_first_breaks=tuple(first),
        load_before=load_before,
        load_after=os.getloadavg(),
    )


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("workload")
    parser.add_argument("--every", type=int, default=20)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args(argv)
    report = store_gate(differential_workload(args.workload), args.every, args.repeats)
    record = report.to_json()
    print(json.dumps(record, indent=1))
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / f"{args.workload}.store_gate.json").write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
    report.require()
    return 0


if __name__ == "__main__":
    sys.exit(main())
