"""Exact last-16 actual-stock circle certificate, one circle at a time (stage C1 of the actual-stock engine).

A machining circle ``i`` is CERTIFIED when no loop station exceeds the cap on
the stock ``Omega minus the last 16 circles (i-16 .. i-1) and their connectors``.
That stock contains the actual stock (Omega minus every prior sweep), and the
45-station decision is monotone under stock inclusion (the clip-monotonicity
lemma of ``src/swept_prefix_engagement_2.h``: a smaller stock has engaged arcs
inside the larger stock's runs), so a last-16 certificate certifies the actual
stock. On a last-16 refusal the circle ESCALATES (user ruling U1, 2026-10-22)
to the FULL window: every prior sweep that can reach the window. Its stock
equals the actual stock on the decision's read set, so its refusal is a true
violation at the stations (nothing is claimed between stations).

Both views are window stocks (``benchmarks.held_actual_stock_last16_window``):
the clip is admitted only by the exact containment proof, and a sweep is left
out only by an exact disjointness proof. Binary64 enters only as cost choices
(the clip radius, which candidates are worth an exact omission proof), never as
a decision. Leaving a sweep out only enlarges the stock, so an omission error
cannot forge a certificate: a certificate's soundness rests on the clip proof
and on never removing a sweep that is not prior; exact omission protects only
the claim that a refusal is a true violation. Sweeps are removed by the global
overlay (``subtract_circle_sweep``, ``subtract_line_sweep``: the connector's
``Stock_UNDER``), which AS2 measured bit-identical to and 7-9x cheaper than
local surgery on window-sized stocks.

Not in generation (ruling 2026-10-21: not until it is fast). Batch runner:
``benchmarks.held_actual_stock_last16_batch``.
"""

from __future__ import annotations

import enum
import time
from dataclasses import dataclass
from typing import NewType
from typing import Optional

import numpy as np

from benchmarks.held_actual_stock_last16_window import ExactCircleReach
from benchmarks.held_actual_stock_last16_window import ExactDisk
from benchmarks.held_actual_stock_last16_window import ExactSegmentReach
from benchmarks.held_actual_stock_last16_window import circle_misses_window
from benchmarks.held_actual_stock_last16_window import prove_clip_contains_stations
from benchmarks.held_actual_stock_last16_window import segment_misses_window
from compas_cgal import _stock_2
from tools.actual_stock_profile import ReplayWorkload

# The ruled certificate depth: circles i-16 .. i-1 and the connectors they own.
LAST_M = 16
# Clip radius beyond the circle's swept disk, in tool radii. A COST choice (the
# window's size): any value the exact clip proof admits is sound.
WINDOW_MARGIN_TOOL_RADII = 0.5
# The clip annulus' outer radius exceeds twice the polygon's bbox diagonal by
# this many tool radii. Cost only: the decision reads nothing outside the window.
CLIP_OUTER_TOOL_RADII = 10.0
# Binary64 pre-filter slack, in tool radii: a sweep whose estimated distance to
# the window is within it is included without an exact omission proof. Cost
# only: inclusion is always sound, omission always needs the exact proof.
PREFILTER_SLACK_TOOL_RADII = 1.0
# The decision runs with no gap closure, as the generator's gate does.
NO_GAP_CLOSURE = 0.0
MS_PER_S = 1e3

CpuMs = NewType("CpuMs", float)
CIRCLE, SEGMENT = 0, 1
# (owner circle, kind, item index): emission order sorts by owner, circle before its connector pieces.
SweepEvent = tuple[int, int, int]


class Last16CertificateError(Exception):
    """Base class of the last-16 certificate failures."""


class WindowStockMismatchError(Last16CertificateError):
    """A window stock differs from the reference stock built from every sweep of its view (an omission or clip was wrong)."""


class BoundBreakError(Last16CertificateError):
    """A circle is exceeded on a stock contained in one on which it is certified: monotonicity or exactness is broken."""


class Verdict(enum.Enum):
    """A circle's certificate outcome."""

    CERTIFIED_LAST16 = "certified_last16"
    CERTIFIED_FULL_WINDOW = "certified_full_window"
    REFUSED = "refused"


@dataclass(frozen=True)
class WindowView:
    """One window stock of one circle and its decision."""

    sweeps: int
    omitted: int
    halfedges: int
    exceeded: bool
    select_ms: CpuMs
    build_ms: CpuMs
    query_ms: CpuMs


@dataclass(frozen=True)
class CircleCertificate:
    """The certificate of one circle: its last-16 view, and the full window view when it escalated."""

    circle: int
    last16: WindowView
    full: Optional[WindowView]

    @property
    def verdict(self) -> Verdict:
        if not self.last16.exceeded:
            return Verdict.CERTIFIED_LAST16
        if self.full is not None and not self.full.exceeded:
            return Verdict.CERTIFIED_FULL_WINDOW
        return Verdict.REFUSED

    @property
    def cpu_ms(self) -> CpuMs:
        views = (self.last16,) if self.full is None else (self.last16, self.full)
        return CpuMs(sum(v.select_ms + v.build_ms + v.query_ms for v in views))


def _cpu_ms(start: float) -> CpuMs:
    return CpuMs(MS_PER_S * (time.process_time() - start))


class Last16Oracle:
    """The certificate over one workload: exact sweep reaches, window construction, and the per-circle decision."""

    def __init__(self, workload: ReplayWorkload, last: int = LAST_M) -> None:
        self.workload = workload
        self.last = last
        r = workload.tool_radius
        self.circle_reach = tuple(ExactCircleReach.build(c.cx, c.cy, c.guide_radius, r) for c in workload.circles)
        self.segment_reach = tuple(ExactSegmentReach.build(s.x0, s.y0, s.x1, s.y1, r) for s in workload.segments)
        self.centres = np.array([(c.cx, c.cy) for c in workload.circles], dtype=np.float64).reshape(-1, 2)
        self.outer = np.array([c.guide_radius for c in workload.circles], dtype=np.float64) + r
        self.owner = np.array([s.owner for s in workload.segments], dtype=np.int64)
        ends = np.array([(s.x0, s.y0, s.x1, s.y1) for s in workload.segments], dtype=np.float64).reshape(-1, 4)
        self.seg_a, self.seg_b = ends[:, :2], ends[:, 2:]
        xy = workload.polygon[:, :2]
        self.clip_outer = 2.0 * float(np.hypot(*(xy.max(axis=0) - xy.min(axis=0)))) + CLIP_OUTER_TOOL_RADII * r

    def window(self, index: int) -> ExactDisk:
        """Circle ``index``'s clip disk, admitted by the exact containment proof.

        Raises:
            WindowClipContainmentError: the binary64-chosen radius fails the proof.
        """
        w, c = self.workload, self.workload.circles[index]
        reach = max(c.guide_radius, float(np.hypot(c.entry_x - c.cx, c.entry_y - c.cy)))
        disk = ExactDisk.build(c.cx, c.cy, reach + (1.0 + WINDOW_MARGIN_TOOL_RADII) * w.tool_radius)
        prove_clip_contains_stations(disk, c.guide_radius, c.entry_x, c.entry_y, w.tool_radius)
        return disk

    def candidates(self, index: int, last: Optional[int]) -> list[SweepEvent]:
        """Every prior sweep of the view, emission-ordered: circles ``j < index`` and connectors owned by them (``j >= index - last``)."""
        first = 0 if last is None else max(0, index - last)
        circles = [(j, CIRCLE, j) for j in range(first, index)]
        segments = [(int(self.owner[k]), SEGMENT, int(k)) for k in np.nonzero((self.owner < index) & (self.owner >= first))[0]]
        return sorted(circles + segments)

    def select(self, index: int, window: ExactDisk, last: Optional[int]) -> tuple[list[SweepEvent], int]:
        """The view's sweeps minus those PROVEN (exactly) to miss the open window; returns (kept, omitted count)."""
        events = self.candidates(index, last)
        radius = float(window.radius)
        centre = self.centres[index]
        slack = PREFILTER_SLACK_TOOL_RADII * self.workload.tool_radius
        kept: list[SweepEvent] = []
        omitted = 0
        circle_ids = np.array([e[2] for e in events if e[1] == CIRCLE], dtype=np.int64)
        seg_ids = np.array([e[2] for e in events if e[1] == SEGMENT], dtype=np.int64)
        circle_far = dict(zip(circle_ids.tolist(), (np.hypot(*(self.centres[circle_ids] - centre).T) > radius + self.outer[circle_ids] + slack).tolist()))
        a, b = self.seg_a[seg_ids], self.seg_b[seg_ids]
        u = b - a
        length = np.einsum("ij,ij->i", u, u)
        t = np.clip(np.einsum("ij,ij->i", centre - a, u) / np.where(length > 0, length, 1.0), 0.0, 1.0)
        seg_far = dict(zip(seg_ids.tolist(), (np.hypot(*(a + t[:, None] * u - centre).T) > radius + self.workload.tool_radius + slack).tolist()))
        for event in events:
            _owner, kind, item = event
            if kind == CIRCLE:
                proven = circle_far[item] and circle_misses_window(window, self.circle_reach[item])
            else:
                proven = seg_far[item] and segment_misses_window(window, self.segment_reach[item])
            if proven:
                omitted += 1
            else:
                kept.append(event)
        return kept, omitted

    def build_stock(self, index: int, window: ExactDisk, events: list[SweepEvent], local: bool = False) -> _stock_2.Stock2:
        """Omega, clipped to ``window``, minus ``events`` (global overlay; ``local`` uses local surgery, for references)."""
        w = self.workload
        stock = w.fresh_stock()
        stock.subtract_annulus(w.circles[index].cx, w.circles[index].cy, float(window.radius), self.clip_outer)
        for _owner, kind, item in events:
            if kind == CIRCLE:
                c = w.circles[item]
                if local:
                    w.subtract_circle(stock, c)
                else:
                    stock.subtract_circle_sweep(c.cx, c.cy, c.guide_radius, w.tool_radius)
            else:
                s = w.segments[item]
                if local:
                    w.subtract_segment(stock, s)
                else:
                    stock.subtract_line_sweep(s.x0, s.y0, s.x1, s.y1, w.tool_radius)
        return stock

    def view(self, index: int, last: Optional[int]) -> tuple[WindowView, _stock_2.Stock2]:
        """Circle ``index`` decided on the window stock of view ``last`` (``None``: the full window)."""
        start = time.process_time()
        window = self.window(index)
        events, omitted = self.select(index, window, last)
        select_ms = _cpu_ms(start)
        start = time.process_time()
        stock = self.build_stock(index, window, events)
        build_ms = _cpu_ms(start)
        start = time.process_time()
        exceeded = self.workload.exceeded(stock, index)
        query_ms = _cpu_ms(start)
        view = WindowView(len(events), omitted, int(stock.arrangement_stats()[1]), bool(exceeded), select_ms, build_ms, query_ms)
        return view, stock

    def certify(self, index: int) -> CircleCertificate:
        """Last-16 decision; on refusal, escalation to the full window (U1)."""
        last16, _ = self.view(index, self.last)
        full = None if not last16.exceeded else self.view(index, None)[0]
        return CircleCertificate(index, last16, full)

    def audit_stock(self, index: int, stock: _stock_2.Stock2, last: Optional[int]) -> None:
        """Check ``stock`` against the reference of its view: clipped Omega minus EVERY sweep of the view, by local surgery.

        Raises:
            WindowStockMismatchError: the point sets differ.
        """
        reference = self.build_stock(index, self.window(index), self.candidates(index, last), local=True)
        if not stock.exactly_equals(reference):
            raise WindowStockMismatchError(f"{self.workload.name} circle {index}: window stock (last={last}) differs from its reference")

    def check_bound(self, index: int) -> bool:
        """Both views of circle ``index``: the full-window stock lies in the last-16 stock and its verdict never exceeds where last-16 certifies.

        Returns the full-window ``exceeded``.

        Raises:
            BoundBreakError: inclusion or verdict monotonicity fails.
        """
        last16, s16 = self.view(index, self.last)
        full, sfull = self.view(index, None)
        if not sfull.is_subset_of(s16):
            raise BoundBreakError(f"{self.workload.name} circle {index}: full-window stock is not inside the last-{self.last} stock")
        if full.exceeded and not last16.exceeded:
            raise BoundBreakError(f"{self.workload.name} circle {index}: exceeded on the full window, certified by last-{self.last}")
        return full.exceeded
