"""Batch runner of the exact last-16 certificate: every circle of a workload over processes, differenced against the full replay.

Circles are independent, so they run in a ``spawn`` process pool; each worker
builds one ``Last16Oracle`` from the parent's pickled workload, and results
are returned in circle order, so the output does not depend on the process
count (timings aside). The workload is identified by a SHA-256 of its sweeps,
so a reference replay recorded by another process is refused if the path
changed in between.

Usage (``build/pyrun.sh -m benchmarks.held_actual_stock_last16_batch ...``):

    WORKLOAD --processes 8 [--every K] [--audit-every K] [--reference REF.json] [--out DIR]
    WORKLOAD --write-reference REF.json        # full sequential exact replay verdicts

``--audit-every K`` additionally checks, on every K-th circle, both window
stocks against their local-surgery references (``exactly_equals``) and the
bound (full window inside last-16, verdict monotone). With ``--reference``,
every certified circle must be clear on the replay (else ``BoundBreakError``)
and every final verdict must equal the replay's (the full window equals the
actual stock on the read set).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing
import os
import resource
import sys
import time
from collections.abc import Iterable
from collections.abc import Sequence
from dataclasses import asdict
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from benchmarks.held_actual_stock_last16 import LAST_M
from benchmarks.held_actual_stock_last16 import BoundBreakError
from benchmarks.held_actual_stock_last16 import CircleCertificate
from benchmarks.held_actual_stock_last16 import Last16CertificateError
from benchmarks.held_actual_stock_last16 import Last16Oracle
from benchmarks.held_actual_stock_last16 import Verdict
from benchmarks.held_actual_stock_last16_reference import replay_exceeded
from benchmarks.performance_record import REPOSITORY
from tools.actual_stock_profile import ReplayWorkload

SCHEMA = "held-actual-stock-last16/1"
DEFAULT_OUT = REPOSITORY / "build" / "c1"
# Circles per task handed to a worker: small enough to balance the uneven per-circle cost.
CHUNK = 8
MS_PER_S = 1e3


class WorkloadMismatchError(Last16CertificateError):
    """A reference record belongs to a different workload (the path changed)."""


class VerdictMismatchError(Last16CertificateError):
    """A final certificate verdict differs from the full replay's exact decision."""


class WorkerNotInitializedError(Last16CertificateError):
    """A batch task ran in a process whose oracle was never built (``_init_worker`` did not run)."""


def workload_digest(workload: ReplayWorkload) -> str:
    """SHA-256 of the workload's sweeps and decision constants, in emission order (doubles by ``float.hex``)."""
    h = hashlib.sha256()
    h.update(f"{workload.name}|{workload.tool_radius.hex()}|{workload.cap_chord_ratio.hex()}".encode())
    for row in workload.polygon:
        h.update(("p" + ",".join(float(v).hex() for v in row)).encode())
    for c in workload.circles:
        h.update(("c" + ",".join(v.hex() for v in (c.cx, c.cy, c.guide_radius, c.entry_x, c.entry_y))).encode())
    for s in workload.segments:
        h.update((f"s{s.owner}," + ",".join(v.hex() for v in (s.x0, s.y0, s.x1, s.y1))).encode())
    return h.hexdigest()


_ORACLE: Optional[Last16Oracle] = None
_AUDIT_EVERY = 0


def _init_worker(workload: ReplayWorkload, audit_every: int) -> None:
    global _ORACLE, _AUDIT_EVERY
    _ORACLE = Last16Oracle(workload)
    _AUDIT_EVERY = audit_every


def _run(indices: Sequence[int]) -> list[CircleCertificate]:
    """Certify ``indices`` with this process's oracle.

    Raises:
        WorkerNotInitializedError: ``_init_worker`` has not run in this process.
    """
    if _ORACLE is None:
        raise WorkerNotInitializedError("last-16 batch task ran before _init_worker built the oracle")
    out = []
    for index in indices:
        certificate = _ORACLE.certify(index)
        if _AUDIT_EVERY and index % _AUDIT_EVERY == 0:
            for last in (LAST_M, None):
                _ORACLE.audit_stock(index, _ORACLE.view(index, last)[1], last)
            _ORACLE.check_bound(index)
        out.append(certificate)
    return out


def certify_batch(workload: ReplayWorkload, indices: Iterable[int], processes: int, audit_every: int = 0) -> tuple[CircleCertificate, ...]:
    """The certificate of every circle in ``indices``, in that order; ``processes == 1`` runs in this process."""
    order = list(indices)
    chunks = [order[k : k + CHUNK] for k in range(0, len(order), CHUNK)]
    if processes == 1:
        _init_worker(workload, audit_every)
        return tuple(c for chunk in chunks for c in _run(chunk))
    with multiprocessing.get_context("spawn").Pool(processes, initializer=_init_worker, initargs=(workload, audit_every)) as pool:
        return tuple(c for result in pool.imap(_run, chunks) for c in result)


@dataclass(frozen=True)
class BatchSummary:
    """Counts and costs of one batch run."""

    workload: str
    digest: str
    circles: int
    processes: int
    certified_last16: int
    certified_full_window: int
    refused: int
    refused_circles: tuple[int, ...]
    escalated_circles: tuple[int, ...]
    circle_cpu_s: float
    last16_cpu_s: float
    escalation_cpu_s: float
    per_circle_ms_mean: float
    per_circle_ms_p50: float
    per_circle_ms_p90: float
    per_circle_ms_max: float
    escalation_ms_mean: float
    children_cpu_s: float  # every worker's CPU including start-up, plus the parent's during the batch
    wall_s: float
    generation_s: float
    load_before: tuple[float, float, float]
    load_after: tuple[float, float, float]
    reference_compared: Optional[int]
    false_refusals_before_escalation: Optional[int]


def _quantile(values: Sequence[float], q: float) -> float:
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, int(q * len(ordered)))] if ordered else 0.0


def summarize(
    workload: ReplayWorkload,
    certificates: Sequence[CircleCertificate],
    processes: int,
    children_cpu_s: float,
    wall_s: float,
    load_before: tuple[float, float, float],
    reference: Optional[dict[int, bool]],
) -> BatchSummary:
    """Counts and costs; with ``reference``, also the exact differential.

    Raises:
        BoundBreakError: a certified circle is exceeded on the replay.
        VerdictMismatchError: a final verdict differs from the replay.
    """
    compared = false_refusals = None
    if reference is not None:
        compared = false_refusals = 0
        for c in certificates:
            if c.circle not in reference:
                continue
            compared += 1
            actual = reference[c.circle]
            if actual and c.verdict is not Verdict.REFUSED:
                raise BoundBreakError(f"{workload.name} circle {c.circle}: exceeded on the replay, certificate {c.verdict.value}")
            if not actual and c.verdict is Verdict.REFUSED:
                raise VerdictMismatchError(f"{workload.name} circle {c.circle}: clear on the replay, refused by the full window")
            false_refusals += int(c.last16.exceeded and not actual)
    per_circle = [c.cpu_ms for c in certificates]
    escalated = [c for c in certificates if c.full is not None]
    last16_ms = [c.last16.select_ms + c.last16.build_ms + c.last16.query_ms for c in certificates]
    escalation_ms = [c.full.select_ms + c.full.build_ms + c.full.query_ms for c in escalated if c.full is not None]
    return BatchSummary(
        workload=workload.name,
        digest=workload_digest(workload),
        circles=len(certificates),
        processes=processes,
        certified_last16=sum(c.verdict is Verdict.CERTIFIED_LAST16 for c in certificates),
        certified_full_window=sum(c.verdict is Verdict.CERTIFIED_FULL_WINDOW for c in certificates),
        refused=sum(c.verdict is Verdict.REFUSED for c in certificates),
        refused_circles=tuple(c.circle for c in certificates if c.verdict is Verdict.REFUSED),
        escalated_circles=tuple(c.circle for c in escalated),
        circle_cpu_s=sum(per_circle) / MS_PER_S,
        last16_cpu_s=sum(last16_ms) / MS_PER_S,
        escalation_cpu_s=sum(escalation_ms) / MS_PER_S,
        per_circle_ms_mean=sum(per_circle) / max(len(per_circle), 1),
        per_circle_ms_p50=_quantile(per_circle, 0.5),
        per_circle_ms_p90=_quantile(per_circle, 0.9),
        per_circle_ms_max=max(per_circle, default=0.0),
        escalation_ms_mean=sum(escalation_ms) / max(len(escalation_ms), 1),
        children_cpu_s=children_cpu_s,
        wall_s=wall_s,
        generation_s=workload.generation_s,
        load_before=load_before,
        load_after=os.getloadavg(),
        reference_compared=compared,
        false_refusals_before_escalation=false_refusals,
    )


def render_markdown(summary: BatchSummary) -> str:
    s = summary
    load = " / ".join(f"{v:.0f}" for v in s.load_before)
    lines = [
        f"<!-- Generated by benchmarks/held_actual_stock_last16_batch.py; workload digest {s.digest[:16]} -->",
        "",
        f"`{s.workload}`: {s.circles:,} circles, {s.processes} processes, load before {load}.",
        "",
        "| certified last-16 | escalated | certified full window | refused | circle CPU | of which escalation "
        "| children CPU | wall | per circle mean / p50 / p90 / max | escalation mean |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: |",
        f"| {s.certified_last16:,} | {len(s.escalated_circles)} | {s.certified_full_window} | {s.refused} | {s.circle_cpu_s:.1f} s | {s.escalation_cpu_s:.1f} s | "
        f"{s.children_cpu_s:.1f} s | {s.wall_s:.1f} s | "
        f"{s.per_circle_ms_mean:.1f} / {s.per_circle_ms_p50:.1f} / {s.per_circle_ms_p90:.1f} / {s.per_circle_ms_max:.1f} ms | {s.escalation_ms_mean:.0f} ms |",
        "",
        f"Refused circles: {list(s.refused_circles)}.",
    ]
    if s.reference_compared is not None:
        lines.append(
            f"Differential vs the full exact replay: {s.reference_compared:,} circles compared, 0 bound breaks, 0 verdict mismatches, "
            f"{s.false_refusals_before_escalation} last-16 false refusals before escalation."
        )
    return "\n".join(lines) + "\n"


def _children_cpu_s() -> float:
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    return usage.ru_utime + usage.ru_stime


def _load_reference(path: Path, digest: str) -> dict[int, bool]:
    record = json.loads(path.read_text(encoding="utf-8"))
    if record["digest"] != digest:
        raise WorkloadMismatchError(f"{path}: digest {record['digest'][:16]} != workload {digest[:16]}")
    return {int(k): bool(v) for k, v in record["exceeded"].items()}


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    parser.add_argument("workload")
    parser.add_argument("--processes", type=int, default=8)
    parser.add_argument("--every", type=int, default=1)
    parser.add_argument("--audit-every", type=int, default=0)
    parser.add_argument("--reference", type=Path, default=None)
    parser.add_argument("--write-reference", type=Path, default=None)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args(argv)
    workload = ReplayWorkload.build(args.workload)
    digest = workload_digest(workload)
    if args.write_reference is not None:
        start = time.process_time()
        exceeded = replay_exceeded(workload, progress=True)
        record = {
            "schema": SCHEMA + "/reference",
            "workload": workload.name,
            "digest": digest,
            "cpu_s": time.process_time() - start,
            "load": os.getloadavg(),
            "exceeded": exceeded,
        }
        args.write_reference.write_text(json.dumps(record) + "\n", encoding="utf-8")
        print(f"reference {workload.name}: {sum(exceeded.values())} of {len(exceeded)} exceeded, {record['cpu_s']:.1f} s CPU")
        return
    reference = None if args.reference is None else _load_reference(args.reference, digest)
    load_before = os.getloadavg()
    children_before = _children_cpu_s() + time.process_time()
    wall = time.perf_counter()
    certificates = certify_batch(workload, range(0, len(workload.circles), args.every), args.processes, args.audit_every)
    wall_s = time.perf_counter() - wall
    children = _children_cpu_s() + time.process_time() - children_before
    summary = summarize(workload, certificates, args.processes, children, wall_s, load_before, reference)
    args.out.mkdir(parents=True, exist_ok=True)
    stem = f"{workload.name}.last16.every{args.every}"
    rows = [{"circle": c.circle, "verdict": c.verdict.value, "last16": asdict(c.last16), "full": None if c.full is None else asdict(c.full)} for c in certificates]
    (args.out / f"{stem}.json").write_text(json.dumps({"schema": SCHEMA, "summary": asdict(summary), "rows": rows}) + "\n", encoding="utf-8")
    markdown = render_markdown(summary)
    (args.out / f"{stem}.md").write_text(markdown, encoding="utf-8")
    sys.stdout.write(markdown)


if __name__ == "__main__":
    main()
