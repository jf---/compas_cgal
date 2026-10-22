"""Exact last-16 actual-stock certificate: window proofs, escalation (U1), injected false certificates, batch determinism."""

from __future__ import annotations

import math
import numpy as np
import pytest

from benchmarks.held_actual_stock_last16 import LAST_M
from benchmarks.held_actual_stock_last16 import BoundBreakError
from benchmarks.held_actual_stock_last16 import Last16Oracle
from benchmarks.held_actual_stock_last16 import Verdict
from benchmarks.held_actual_stock_last16 import WindowStockMismatchError
from benchmarks import held_actual_stock_last16_batch as batch
from benchmarks.held_actual_stock_last16_batch import WorkerNotInitializedError
from benchmarks.held_actual_stock_last16_batch import certify_batch
from benchmarks.held_actual_stock_last16_reference import replay_exceeded
from benchmarks.held_actual_stock_last16_window import ExactCircleReach
from benchmarks.held_actual_stock_last16_window import ExactDisk
from benchmarks.held_actual_stock_last16_window import ExactSegmentReach
from benchmarks.held_actual_stock_last16_window import WindowClipContainmentError
from benchmarks.held_actual_stock_last16_window import circle_misses_window
from benchmarks.held_actual_stock_last16_window import prove_clip_contains_stations
from benchmarks.held_actual_stock_last16_window import segment_misses_window
from compas_cgal import _stock_2
from tools.actual_stock_profile import CircleSweep
from tools.actual_stock_profile import ReplayWorkload
from tools.actual_stock_profile import SegmentSweep
from tools.actual_stock_profile import run_windows

TOOL = 1.0
CAP_DEG = 80.0
# A big circle at the origin: annulus [2, 4]. Its duplicate later reads only cut material.
BIG = (0.0, 0.0, 3.0)
FAR_X = 30.0


def _circle(cx: float, cy: float, rho: float) -> CircleSweep:
    return CircleSweep(cx, cy, rho, cx + rho, cy)


def _workload(circles: list[CircleSweep], segments: tuple[SegmentSweep, ...] = ()) -> ReplayWorkload:
    polygon = np.array([(-10.0, -10.0, 0.0), (50.0, -10.0, 0.0), (50.0, 10.0, 0.0), (-10.0, 10.0, 0.0)], dtype=np.float64)
    return ReplayWorkload(
        name="synthetic",
        polygon=polygon,
        tool_radius=TOOL,
        cap_chord_ratio=_stock_2.cap_chord_ratio(math.radians(CAP_DEG)),
        circles=tuple(circles),
        segments=segments,
        generation_s=0.0,
    )


def _escalation_workload() -> ReplayWorkload:
    """Circle 0 = BIG; LAST_M small circles far away; then BIG again: refused by last-16, clear on the actual stock."""
    fillers = [_circle(FAR_X + 2.0 * k / LAST_M, 0.0, 0.5) for k in range(LAST_M)]
    return _workload([_circle(*BIG), *fillers, _circle(*BIG)])


# ---- exact window geometry -------------------------------------------------


def test_clip_proof_admits_margin_and_refuses_tangent_and_shrunk() -> None:
    prove_clip_contains_stations(ExactDisk.build(0.0, 0.0, 4.5), 3.0, 3.0, 0.0, TOOL)
    with pytest.raises(WindowClipContainmentError):
        prove_clip_contains_stations(ExactDisk.build(0.0, 0.0, 4.0), 3.0, 3.0, 0.0, TOOL)
    with pytest.raises(WindowClipContainmentError):
        prove_clip_contains_stations(ExactDisk.build(0.0, 0.0, 2.2), 3.0, 3.0, 0.0, TOOL)
    # The entry station off the guide circle is checked on its own.
    with pytest.raises(WindowClipContainmentError):
        prove_clip_contains_stations(ExactDisk.build(0.0, 0.0, 4.5), 3.0, 3.6, 0.0, TOOL)


def test_omission_predicates_are_exact_at_tangency() -> None:
    window = ExactDisk.build(0.0, 0.0, 2.0)
    # Annulus outer circle externally tangent to the open window: disjoint.
    assert circle_misses_window(window, ExactCircleReach.build(5.0, 0.0, 2.0, TOOL))
    # One ulp closer: they meet.
    assert not circle_misses_window(window, ExactCircleReach.build(math.nextafter(5.0, 0.0), 0.0, 2.0, TOOL))
    # Window inside the annulus hole: hole radius 6, d + R = 4.
    assert circle_misses_window(window, ExactCircleReach.build(2.0, 0.0, 7.0, TOOL))
    assert not circle_misses_window(window, ExactCircleReach.build(2.0, 0.0, 4.5, TOOL))
    # Capsule at distance exactly R + r: disjoint; nearer: meets.
    assert segment_misses_window(window, ExactSegmentReach.build(-5.0, 3.0, 5.0, 3.0, TOOL))
    assert not segment_misses_window(window, ExactSegmentReach.build(-5.0, math.nextafter(3.0, 0.0), 5.0, 3.0, TOOL))
    assert segment_misses_window(window, ExactSegmentReach.build(3.0, 0.0, 9.0, 0.0, TOOL))


# ---- escalation (U1) ---------------------------------------------------------


def test_last16_refusal_escalates_to_full_window_certificate() -> None:
    w = _escalation_workload()
    oracle = Last16Oracle(w)
    index = len(w.circles) - 1
    certificate = oracle.certify(index)
    assert certificate.last16.exceeded
    assert certificate.full is not None and not certificate.full.exceeded
    assert certificate.verdict is Verdict.CERTIFIED_FULL_WINDOW
    assert replay_exceeded(w)[index] is False


def test_plunge_is_a_true_refusal_and_clear_circles_do_not_escalate() -> None:
    w = _escalation_workload()
    oracle = Last16Oracle(w)
    assert oracle.certify(0).verdict is Verdict.REFUSED
    assert replay_exceeded(w)[0] is True


# ---- the window stock equals its reference; injected defects are caught ----------


def test_window_stocks_equal_references_and_full_lies_in_last16() -> None:
    w = _escalation_workload()
    oracle = Last16Oracle(w)
    index = len(w.circles) - 1
    for last in (LAST_M, None):
        _view, stock = oracle.view(index, last)
        oracle.audit_stock(index, stock, last)
    assert oracle.check_bound(index) is False


def test_dropped_window_sweep_is_caught() -> None:
    w = _workload([_circle(0.0, 0.0, 3.0), _circle(0.5, 0.0, 3.0), _circle(1.0, 0.0, 3.0)])
    oracle = Last16Oracle(w)
    window = oracle.window(2)
    events, _ = oracle.select(2, window, LAST_M)
    assert len(events) == 2
    with pytest.raises(WindowStockMismatchError):
        oracle.audit_stock(2, oracle.build_stock(2, window, events[1:]), LAST_M)


def test_shrunk_clip_would_forge_a_certificate_and_is_refused_by_the_proof() -> None:
    w = _workload([_circle(*BIG)])
    oracle = Last16Oracle(w)
    assert oracle.certify(0).last16.exceeded
    shrunk = ExactDisk.build(0.0, 0.0, 2.2)
    # Bypassing the proof: the shrunk clip removes read material and forges a certificate.
    assert not w.exceeded(oracle.build_stock(0, shrunk, []), 0)
    with pytest.raises(WindowClipContainmentError):
        prove_clip_contains_stations(shrunk, 3.0, 3.0, 0.0, TOOL)


def test_clip_does_not_change_the_decision() -> None:
    w = _workload([_circle(0.0, 0.0, 3.0), _circle(0.8, 0.0, 3.0), _circle(1.6, 0.3, 2.5)])
    oracle = Last16Oracle(w)
    for index in range(len(w.circles)):
        unclipped = w.fresh_stock()
        for j in range(index):
            w.subtract_circle(unclipped, w.circles[j])
        assert oracle.view(index, None)[0].exceeded == w.exceeded(unclipped, index)


def test_bound_break_is_raised_on_an_inverted_pair(monkeypatch: pytest.MonkeyPatch) -> None:
    w = _escalation_workload()
    oracle = Last16Oracle(w)
    index = len(w.circles) - 1
    real = oracle.view

    def swapped(i: int, last: object) -> object:
        return real(i, None if last is not None else LAST_M)

    monkeypatch.setattr(oracle, "view", swapped)
    with pytest.raises(BoundBreakError):
        oracle.check_bound(index)


# ---- batch ---------------------------------------------------------------------


def test_batch_is_deterministic_across_process_counts() -> None:
    w = _escalation_workload()
    serial = certify_batch(w, range(len(w.circles)), processes=1)
    parallel = certify_batch(w, range(len(w.circles)), processes=2)
    assert [c.circle for c in parallel] == list(range(len(w.circles)))
    assert [(c.circle, c.verdict, c.last16.sweeps, c.last16.exceeded, c.full is None) for c in serial] == [
        (c.circle, c.verdict, c.last16.sweeps, c.last16.exceeded, c.full is None) for c in parallel
    ]


def test_unstarted_worker_raises_a_named_error(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(batch, "_ORACLE", None)
    with pytest.raises(WorkerNotInitializedError):
        batch._run([0])


# ---- connector segments ----------------------------------------------------------


def _connector_workload() -> ReplayWorkload:
    """Three overlapping circles; circle 0 owns an oblique connector (non-Pythagorean: removed as ``Stock_UNDER``)."""
    connector = SegmentSweep(0, 3.0, 0.0, 3.5, 0.7)
    return _workload([_circle(0.0, 0.0, 3.0), _circle(0.8, 0.0, 3.0), _circle(1.6, 0.3, 2.5)], (connector,))


def test_connector_segment_window_stocks_equal_references_and_replay() -> None:
    w = _connector_workload()
    oracle = Last16Oracle(w)
    reference = replay_exceeded(w)
    for index in range(1, len(w.circles)):
        events, _omitted = oracle.select(index, oracle.window(index), LAST_M)
        assert (0, 1, 0) in events, "the connector reaches every later window"
        for last in (LAST_M, None):
            _view, stock = oracle.view(index, last)
            oracle.audit_stock(index, stock, last)
        assert oracle.check_bound(index) == reference[index]
        assert (oracle.certify(index).verdict is Verdict.REFUSED) == reference[index]


# ---- the profile's window rows are the certificate's window stocks --------------


def test_profile_windows_are_the_certificate_views() -> None:
    """An entry station off the guide circle widens the certificate's window; the profile must use that same window and selection."""
    off_guide = CircleSweep(1.6, 0.3, 2.5, 1.6 + 2.5 + 0.6 * TOOL, 0.3)
    w = _workload([_circle(0.0, 0.0, 3.0), _circle(0.8, 0.0, 3.0), off_guide], (SegmentSweep(0, 3.0, 0.0, 3.5, 0.7),))
    oracle = Last16Oracle(w)
    rows = run_windows(w, 1)
    assert [r.circle for r in rows] == [1, 2]
    for row in rows:
        for stock, last in ((row.exact, None), (row.last16, LAST_M), (row.last1, 1)):
            view, _ = oracle.view(row.circle, last)
            assert (stock.sweeps, stock.halfedges, stock.query.exceeded) == (view.sweeps, view.halfedges, view.exceeded)
