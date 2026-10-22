"""Reference for the last-16 certificate: the full sequential exact replay of the actual stock.

Emission order, global overlay (``subtract_circle_sweep`` / ``subtract_line_sweep``,
connectors as ``Stock_UNDER``): circle ``i`` is decided on the stock as it
stands before it (every circle ``j < i`` and every connector they own
removed), then removed with its connector. No window, no clip, no selection:
the independent oracle the certificate is differenced against.
"""

from __future__ import annotations

import sys
import time
from collections.abc import Collection
from typing import Optional

from tools.actual_stock_profile import ReplayWorkload

# Progress line every this many circles on stderr (long Monstera runs).
PROGRESS_EVERY = 250


def replay_exceeded(workload: ReplayWorkload, query: Optional[Collection[int]] = None, progress: bool = False) -> dict[int, bool]:
    """Exact actual-stock decision of every circle in ``query`` (default: all), by one sequential replay."""
    owned = workload.segments_by_owner()
    wanted = set(range(len(workload.circles))) if query is None else set(query)
    last = max(wanted, default=-1)
    stock = workload.fresh_stock()
    out: dict[int, bool] = {}
    start = time.process_time()
    for index, circle in enumerate(workload.circles):
        if index > last:
            break
        if index in wanted:
            out[index] = bool(workload.exceeded(stock, index))
        stock.subtract_circle_sweep(circle.cx, circle.cy, circle.guide_radius, workload.tool_radius)
        for s in owned[index]:
            stock.subtract_line_sweep(s.x0, s.y0, s.x1, s.y1, workload.tool_radius)
        if progress and index % PROGRESS_EVERY == 0:
            print(f"replay {workload.name}: circle {index:,}/{len(workload.circles):,}, {time.process_time() - start:.1f} s CPU", file=sys.stderr, flush=True)
    return out
