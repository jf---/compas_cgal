"""Differential oracle of the actual-stock refinement decision (unit U0 of the 2026-11-04 design): D2 and D4, off the production path.

Two exact properties of an emitted polygon-lane path, both decided by the
last-16 certificate (``Last16Oracle``) on the part polygon P:

- **D4, circle view vs actual.** C_t (P minus every prior circle's swept
  annulus, connectors left out; ``SweepView.CIRCLES``) contains the actual
  stock S_t (``SweepView.ACTUAL``), so the circle view refuses at least what
  the actual stock refuses. D4 counts the circles it refuses BEYOND that:
  accepted on S_t, refused on C_t. Each such circle is a midpoint the future
  refinement decision would insert that the certificate does not need. The
  harness decides C_t on every circle and S_t only where C_t refuses (plus
  every ``audit_every``-th circle, where an S_t refusal under a C_t
  certificate is a ``BoundBreakError``: inclusion or exactness is broken).
  Circle 0 has no prior sweep, so its two views are one stock and it can
  never count.
- **D2, production implies certificate.** Every emitted circle t >= 1 is
  certified (last-16 or full window), and the entry circle 0 is refused (its
  engagement is owned by the entry move, ruling R3; a certified circle 0
  means the oracle is not reading P). ``sweeps=SweepView.CIRCLES`` checks the
  stronger circle-view statement the refinement decision will make.

Both are reports with a ``require()`` that raises a named error. Circles are
independent, so they run in a ``spawn`` pool; results come back in index
order, so output never depends on the process count. Nothing here is called
by a generator.

CLI (``build/pyrun.sh -m benchmarks.held_actual_stock_oracle ...``):

    d4 WORKLOAD [--every K] [--audit-every K] [--processes N] [--out DIR]
    d2 WORKLOAD [--every K] [--circles] [--processes N] [--out DIR]
"""

from __future__ import annotations

import argparse
import json
import multiprocessing
import os
import sys
import time
from collections.abc import Iterable
from collections.abc import Sequence
from dataclasses import asdict
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from benchmarks.held_actual_stock_last16 import BoundBreakError
from benchmarks.held_actual_stock_last16 import Last16CertificateError
from benchmarks.held_actual_stock_last16 import Last16Oracle
from benchmarks.held_actual_stock_last16 import SweepView
from benchmarks.held_actual_stock_last16 import Verdict
from benchmarks.held_actual_stock_last16_batch import workload_digest
from benchmarks.held_material_domain import PartPolygon
from benchmarks.held_material_domain import require_part_material
from benchmarks.held_toolpath_benchmark import GeneratedToolpath
from benchmarks.held_toolpath_benchmark import generate
from benchmarks.performance_record import REPOSITORY
from tools.actual_stock_profile import ReplayWorkload

SCHEMA = "held-actual-stock-oracle/1"
SKIS = "figure8_crossed_skis_boundary"
DEFAULT_OUT = REPOSITORY / "build" / "actual-stock-u0"
# The entry circle: its engagement is owned by the entry move (ruling R3).
ENTRY_CIRCLE = 0
# Circles per task handed to a worker: small enough to balance the uneven per-circle cost.
CHUNK = 8
CERTIFIED = (Verdict.CERTIFIED_LAST16, Verdict.CERTIFIED_FULL_WINDOW)

StationPrefix = tuple[float, float, bool]


class ActualStockOracleError(Last16CertificateError):
    """Base class of the differential-oracle failures."""


class CircleViewRefusalError(ActualStockOracleError):
    """D4: the circle view refuses circles the actual stock accepts."""


class CertificateViolationError(ActualStockOracleError):
    """D2: an emitted circle after the entry is refused by the certificate."""


class EntryCircleNotRefusedError(ActualStockOracleError):
    """D2: the entry circle is certified, so the oracle is not deciding the plunge on P."""


class OracleWorkerNotInitializedError(ActualStockOracleError):
    """An oracle task ran in a process whose workload was never installed (``_init_worker`` did not run)."""


def workload_on_part(path: GeneratedToolpath, generation_s: float) -> ReplayWorkload:
    """The replay workload of ``path``, checked exactly to decide on the part polygon P = ``path.case.projection.points``.

    Raises:
        MaterialDomainMismatchError: the workload's stock polygon is not P.
    """
    workload = ReplayWorkload.from_path(path, generation_s)
    require_part_material(workload.polygon, PartPolygon.from_case(path.case), path.workload)
    return workload


def differential_workload(name: str) -> ReplayWorkload:
    """Generate ``name`` (process CPU recorded) and return its workload on P.

    Raises:
        MaterialDomainMismatchError: the workload's stock polygon is not P.
    """
    start = time.process_time()
    path = generate(name)
    return workload_on_part(path, time.process_time() - start)


def view_prefixes(oracle: Last16Oracle, index: int, last: Optional[int]) -> tuple[StationPrefix, ...]:
    """Per station of circle ``index``: (start, max run, exceeded) from ``engagement_circle_prefixes`` on the oracle's window stock of view ``last``.

    The reference side of D1: with ``SweepView.CIRCLES`` this is the existing
    engine on the circle-only window on P.
    """
    _view, stock = oracle.view(index, last)
    return tuple(oracle.workload.prefixes(stock, index))


# ---- process pool ------------------------------------------------------------------

_WORKLOAD: Optional[ReplayWorkload] = None
_ORACLES: dict[SweepView, Last16Oracle] = {}


def _init_worker(workload: ReplayWorkload) -> None:
    global _WORKLOAD
    _WORKLOAD = workload
    _ORACLES.clear()


def _oracle(view: SweepView) -> Last16Oracle:
    """This process's oracle of ``view``, built on first use.

    Raises:
        OracleWorkerNotInitializedError: ``_init_worker`` has not run in this process.
    """
    if _WORKLOAD is None:
        raise OracleWorkerNotInitializedError("oracle task ran before _init_worker installed the workload")
    if view not in _ORACLES:
        _ORACLES[view] = Last16Oracle(_WORKLOAD, sweeps=view)
    return _ORACLES[view]


@dataclass(frozen=True)
class D4Row:
    """One circle's circle-view verdict, and its actual verdict where it was decided."""

    circle: int
    circle_view: Verdict
    actual: Optional[Verdict]
    cpu_ms: float

    @property
    def violation(self) -> bool:
        """Refused on C_t, accepted on S_t."""
        return self.circle_view is Verdict.REFUSED and self.actual is not None and self.actual is not Verdict.REFUSED


def _d4_row(index: int, audit_every: int) -> D4Row:
    """Decide circle ``index`` on C_t, and on S_t where C_t refuses or the circle is audited.

    Raises:
        BoundBreakError: S_t refuses where C_t (its superset) certifies.
    """
    start = time.process_time()
    circle_view = _oracle(SweepView.CIRCLES).certify(index).verdict
    actual: Optional[Verdict] = None
    if circle_view is Verdict.REFUSED or (audit_every > 0 and index % audit_every == 0):
        actual = _oracle(SweepView.ACTUAL).certify(index).verdict
        if actual is Verdict.REFUSED and circle_view is not Verdict.REFUSED:
            raise BoundBreakError(f"circle {index}: refused on the actual stock, {circle_view.value} on the circle view that contains it")
    return D4Row(index, circle_view, actual, 1e3 * (time.process_time() - start))


def _d4_task(task: tuple[Sequence[int], int]) -> list[D4Row]:
    indices, audit_every = task
    return [_d4_row(index, audit_every) for index in indices]


def _d2_task(task: tuple[Sequence[int], SweepView]) -> list[tuple[int, Verdict]]:
    indices, view = task
    oracle = _oracle(view)
    return [(index, oracle.certify(index).verdict) for index in indices]


def _chunks(indices: Iterable[int]) -> list[list[int]]:
    order = list(indices)
    return [order[k : k + CHUNK] for k in range(0, len(order), CHUNK)]


# ---- D4 ------------------------------------------------------------------------------


@dataclass(frozen=True)
class D4Report:
    """D4 over the compared circles of one workload."""

    workload: str
    digest: str
    compared: int
    audited: int  # circles whose actual verdict was decided (circle-view refusals and audit samples)
    rows: tuple[D4Row, ...]
    violations: tuple[int, ...]
    cpu_s: float
    wall_s: float
    load_before: tuple[float, float, float]
    load_after: tuple[float, float, float]

    def require(self) -> None:
        """Raises:
        CircleViewRefusalError: some circle is refused on C_t and accepted on S_t.
        """
        if self.violations:
            raise CircleViewRefusalError(f"{self.workload}: {len(self.violations)} circles refused on the circle view, accepted on the actual stock: {list(self.violations[:20])}")


def circle_view_vs_actual(workload: ReplayWorkload, indices: Iterable[int], processes: int, audit_every: int = 0) -> D4Report:
    """D4 on every circle of ``indices``; ``processes == 1`` runs in this process.

    Raises:
        BoundBreakError: an audited or refused circle breaks C_t ⊇ S_t.
    """
    load_before = os.getloadavg()
    wall = time.perf_counter()
    tasks = [(chunk, audit_every) for chunk in _chunks(indices)]
    if processes == 1:
        _init_worker(workload)
        rows = tuple(row for task in tasks for row in _d4_task(task))
    else:
        with multiprocessing.get_context("spawn").Pool(processes, initializer=_init_worker, initargs=(workload,)) as pool:
            rows = tuple(row for result in pool.imap(_d4_task, tasks) for row in result)
    return D4Report(
        workload=workload.name,
        digest=workload_digest(workload),
        compared=len(rows),
        audited=sum(r.actual is not None for r in rows),
        rows=rows,
        violations=tuple(r.circle for r in rows if r.violation),
        cpu_s=sum(r.cpu_ms for r in rows) / 1e3,
        wall_s=time.perf_counter() - wall,
        load_before=load_before,
        load_after=os.getloadavg(),
    )


def render_d4(report: D4Report) -> str:
    r = report
    refused = sum(row.circle_view is Verdict.REFUSED for row in r.rows)
    return (
        f"D4 `{r.workload}` (digest {r.digest[:16]}): {r.compared:,} circles, circle view refused {refused}, "
        f"actual decided on {r.audited}, violations {len(r.violations)} {list(r.violations[:20])}; "
        f"circle CPU {r.cpu_s:.1f} s, wall {r.wall_s:.1f} s, load {r.load_before[0]:.0f} -> {r.load_after[0]:.0f}\n"
    )


# ---- D2 ------------------------------------------------------------------------------


@dataclass(frozen=True)
class D2Report:
    """D2 over the compared circles of one workload, in one view."""

    workload: str
    digest: str
    sweeps: SweepView
    compared: int
    entry: Optional[Verdict]
    certified_last16: int
    certified_full_window: int
    violations: tuple[int, ...]
    wall_s: float
    load_before: tuple[float, float, float]

    def require(self) -> None:
        """Raises:
        EntryCircleNotRefusedError: circle 0 was compared and certified.
        CertificateViolationError: some circle t >= 1 is refused.
        """
        if self.entry is not None and self.entry is not Verdict.REFUSED:
            raise EntryCircleNotRefusedError(f"{self.workload}: entry circle {self.entry.value} on the {self.sweeps.value} view; it must be refused on P")
        if self.violations:
            raise CertificateViolationError(
                f"{self.workload}: {len(self.violations)} circles after the entry refused on the {self.sweeps.value} view: {list(self.violations[:20])}"
            )


def production_implies_certificate(workload: ReplayWorkload, indices: Iterable[int], processes: int, sweeps: SweepView = SweepView.ACTUAL) -> D2Report:
    """D2 on every circle of ``indices`` in view ``sweeps``; ``processes == 1`` runs in this process."""
    load_before = os.getloadavg()
    wall = time.perf_counter()
    tasks = [(chunk, sweeps) for chunk in _chunks(indices)]
    if processes == 1:
        _init_worker(workload)
        verdicts = [v for task in tasks for v in _d2_task(task)]
    else:
        with multiprocessing.get_context("spawn").Pool(processes, initializer=_init_worker, initargs=(workload,)) as pool:
            verdicts = [v for result in pool.imap(_d2_task, tasks) for v in result]
    entry = next((verdict for index, verdict in verdicts if index == ENTRY_CIRCLE), None)
    return D2Report(
        workload=workload.name,
        digest=workload_digest(workload),
        sweeps=sweeps,
        compared=len(verdicts),
        entry=entry,
        certified_last16=sum(v is Verdict.CERTIFIED_LAST16 for _i, v in verdicts),
        certified_full_window=sum(v is Verdict.CERTIFIED_FULL_WINDOW for _i, v in verdicts),
        violations=tuple(i for i, v in verdicts if i != ENTRY_CIRCLE and v not in CERTIFIED),
        wall_s=time.perf_counter() - wall,
        load_before=load_before,
    )


def render_d2(report: D2Report) -> str:
    r = report
    entry = "not compared" if r.entry is None else r.entry.value
    return (
        f"D2 `{r.workload}` (digest {r.digest[:16]}, {r.sweeps.value} view): {r.compared:,} circles, entry {entry}, "
        f"certified last-16 {r.certified_last16:,}, full window {r.certified_full_window}, violations {len(r.violations)} {list(r.violations[:20])}; "
        f"wall {r.wall_s:.1f} s, load {r.load_before[0]:.0f}\n"
    )


# ---- CLI ---------------------------------------------------------------------------


def _json_default(value: object) -> object:
    if isinstance(value, (Verdict, SweepView)):
        return value.value
    raise TypeError(f"not JSON serializable: {type(value).__name__}")


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    parser.add_argument("check", choices=("d4", "d2"))
    parser.add_argument("workload")
    parser.add_argument("--every", type=int, default=1)
    parser.add_argument("--audit-every", type=int, default=0)
    parser.add_argument("--circles", action="store_true", help="d2 on the circle view instead of the actual stock")
    parser.add_argument("--processes", type=int, default=2)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args(argv)
    workload = differential_workload(args.workload)
    indices = range(0, len(workload.circles), args.every)
    args.out.mkdir(parents=True, exist_ok=True)
    stem = f"{workload.name}.{args.check}.every{args.every}"
    if args.check == "d4":
        d4 = circle_view_vs_actual(workload, indices, args.processes, args.audit_every)
        record, text = asdict(d4), render_d4(d4)
    else:
        d2 = production_implies_certificate(workload, indices, args.processes, SweepView.CIRCLES if args.circles else SweepView.ACTUAL)
        record, text = asdict(d2), render_d2(d2)
        stem += ".circles" if args.circles else ""
    (args.out / f"{stem}.json").write_text(json.dumps({"schema": SCHEMA, "generation_s": workload.generation_s, **record}, default=_json_default) + "\n", encoding="utf-8")
    (args.out / f"{stem}.md").write_text(text, encoding="utf-8")
    sys.stdout.write(text)


if __name__ == "__main__":
    main()
