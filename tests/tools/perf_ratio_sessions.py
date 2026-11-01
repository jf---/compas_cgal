"""Synthetic ``perf-ratio/1`` sessions for unit tests: exact CPUs, no timing.

A session is built from per-workload ``(ref_cpu, head_cpu)`` pairs through the
record's own checked builders (`PairRun.build`), with pair orders alternating as
the harness runs them (`tools.perf_ratio.pair_order`), so every derived value a
test reads is the record's own arithmetic on the given CPUs.
"""

from __future__ import annotations

from collections.abc import Mapping
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

from tools.perf_ratio import pair_order
from tools.perf_ratio_record import ArmIdentity
from tools.perf_ratio_record import ArmName
from tools.perf_ratio_record import CpuSeconds
from tools.perf_ratio_record import Injection
from tools.perf_ratio_record import LoadOneMinute
from tools.perf_ratio_record import PairRun
from tools.perf_ratio_record import PerfRatioSession
from tools.perf_ratio_record import WorkloadRatio

# The reference commit the landed A/A calibration records were measured at.
REF = "15fbed70ebb016bdd8ef49453e04606d03e0b1fe"
# A commit id other than REF, for an arm measured elsewhere.
OTHER_REV = "0123456789abcdef0123456789abcdef01234567"
# The CPU model of the calibrated machine.
CPU_MODEL = "Apple M1 Max"
# A harness digest (any 64-hex string).
HARNESS = "h" * 64
# 1-minute load written around every synthetic pair (any finite value; never decided on).
SYNTHETIC_LOAD = LoadOneMinute(3.0)


def synthetic_arm(rev: str, name: ArmName, *, native: str = "n" * 64) -> ArmIdentity:
    """A clean arm identity at *rev*, its paths named after the arm."""
    return ArmIdentity(
        rev=rev,
        dirty=False,
        patch_sha256=None,
        native_sha256=native,
        native_sources_sha256="s" * 64,
        environment="default",
        python_version="3.12.7",
        stock_2_sha256="k" * 64,
        pixi_lock_sha256="l" * 64,
        compas_cgal_file=f"/trees/{name.value}/src/compas_cgal/__init__.py",
        stock_2_origin=f"/trees/{name.value}/build/_stock_2.so",
        editable_build_dir=f"/trees/{name.value}/build",
        interpreter_prefix=f"/trees/{name.value}/.pixi/envs/default",
    )


def synthetic_session(
    ratios: Mapping[str, Sequence[tuple[float, float]]],
    *,
    ref_commit: str = REF,
    head_rev: str = REF,
    cpu_model: str = CPU_MODEL,
    harness: str = HARNESS,
    injected: Optional[Injection] = None,
    head_digest: str = "d",
    ref_digest: str = "d",
    native: str = "n" * 64,
) -> PerfRatioSession:
    """A session whose workload *w* runs the pairs ``ratios[w]``, each ``(ref_cpu, head_cpu)``.

    Args:
        ratios: Workload name -> its pairs' (reference CPU, HEAD CPU) seconds, in run order.
        ref_commit: The reference arm's rev (and the session's ``ref_commit``).
        head_rev: The HEAD arm's rev; REF by default, an A/A session.
        cpu_model: The machine's CPU model.
        harness: The session's ``harness_sha256``.
        injected: The audit's injected spin, None for an honest session.
        head_digest: Every workload's HEAD output digest.
        ref_digest: Every workload's reference output digest.
        native: Both arms' ``native_sha256``.

    Returns:
        The session, built from checked pairs.
    """
    workloads = tuple(
        WorkloadRatio(
            workload=name,
            family=name.split("/")[0],
            pairs=tuple(
                PairRun.build(
                    order=pair_order(index),
                    ref_cpu=CpuSeconds(ref_cpu),
                    head_cpu=CpuSeconds(head_cpu),
                    load_before=SYNTHETIC_LOAD,
                    load_after=SYNTHETIC_LOAD,
                )
                for index, (ref_cpu, head_cpu) in enumerate(pairs)
            ),
            ref_digest=ref_digest,
            head_digest=head_digest,
        )
        for name, pairs in ratios.items()
    )
    return PerfRatioSession(
        ref_commit=ref_commit,
        arms={ArmName.REF: synthetic_arm(ref_commit, ArmName.REF, native=native), ArmName.HEAD: synthetic_arm(head_rev, ArmName.HEAD, native=native)},
        workloads=workloads,
        host="bench-host",
        cpu_model=cpu_model,
        started_utc="2026-11-01T10:00:00+00:00",
        ended_utc="2026-11-01T10:16:00+00:00",
        harness_sha256=harness,
        injected=injected,
    )


def write_record(session: PerfRatioSession, path: Path) -> Path:
    """Write *session* as its record (`PerfRatioSession.dumps`) at *path*, creating parents; returns *path*."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(session.dumps(), encoding="utf-8")
    return path
