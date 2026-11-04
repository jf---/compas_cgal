"""Merging an extension session's pairs into a gated workload's run.

The gate (`tools.perf_ratio_gate_verdict.gate`) merges each extension of a red
workload as it measures it; the anchor chain
(`tools.perf_ratio_gate_chain.derive_anchors`) replays the same merges over a
promotion's committed records. Both read one function, so a link's ``rho_w`` is
the ratio the gate decided on. This module sits below both: it depends on the
record only.
"""

from __future__ import annotations

from collections.abc import Callable

from tools.perf_ratio_gate_errors import GateRefusal
from tools.perf_ratio_record import ArmName
from tools.perf_ratio_record import PerfRatioSession
from tools.perf_ratio_record import WorkloadRatio


class ExtensionMismatchError(GateRefusal):
    """An extension session did not measure what its base session and request measured."""


# Merges one extension session into a workload's run: (run so far, base session, extension) -> merged run.
Merge = Callable[[WorkloadRatio, PerfRatioSession, PerfRatioSession], WorkloadRatio]


def merge_extension(base: WorkloadRatio, base_session: PerfRatioSession, extension: PerfRatioSession) -> WorkloadRatio:
    """*base* with *extension*'s pairs of that workload appended, each pair as measured.

    Args:
        base: The workload's run so far.
        base_session: The session *base* comes from.
        extension: The extension session: *base*'s workload alone, measured as *base_session* was.

    Returns:
        The merged run: *base*'s pairs then the extension's, *base*'s name, family and digests.

    Raises:
        ExtensionMismatchError: The extension's arms, harness, CPU model, host, reference or
            injection differ from *base_session*'s, it measures other workloads than *base*'s,
            or its output digests (or family) differ from *base*'s.
    """
    for arm in ArmName:
        if extension.arms[arm] != base_session.arms[arm]:
            raise ExtensionMismatchError(f"extension's {arm.value} arm is {extension.arms[arm]}, not the base session's {base_session.arms[arm]}; it measured another build.")
    for label, measured, expected in (
        ("harness_sha256", extension.harness_sha256, base_session.harness_sha256),
        ("cpu_model", extension.cpu_model, base_session.cpu_model),
        ("host", extension.host, base_session.host),
        ("ref_commit", extension.ref_commit, base_session.ref_commit),
        ("injected", extension.injected, base_session.injected),
    ):
        if measured != expected:
            raise ExtensionMismatchError(f"extension's {label} is {measured!r}, not the base session's {expected!r}.")
    names = tuple(workload.workload for workload in extension.workloads)
    if names != (base.workload,):
        raise ExtensionMismatchError(f"extension measures {list(names)}, not exactly the requested {base.workload!r}.")
    (run,) = extension.workloads
    for label, measured, expected in (
        ("family", run.family, base.family),
        ("ref_digest", run.ref_digest, base.ref_digest),
        ("head_digest", run.head_digest, base.head_digest),
    ):
        if measured != expected:
            raise ExtensionMismatchError(f"extension's {base.workload} {label} is {measured!r}, not the base run's {expected!r}.")
    return WorkloadRatio(workload=base.workload, family=base.family, pairs=base.pairs + run.pairs, ref_digest=base.ref_digest, head_digest=base.head_digest)
