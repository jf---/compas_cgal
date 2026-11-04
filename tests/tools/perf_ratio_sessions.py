"""Synthetic ``perf-ratio/1`` sessions for unit tests: exact CPUs, no timing.

A session is built from per-workload ``(ref_cpu, head_cpu)`` pairs through the
record's own checked builders (`PairRun.build`), with pair orders alternating as
the harness runs them (`tools.perf_ratio.pair_order`), so every derived value a
test reads is the record's own arithmetic on the given CPUs.

`promotion_link` writes one synthetic anchor-chain link under a test root: the
promoted gate session and its extensions in the promotion layout
(`tools.perf_ratio_gate_registry.promotion_dir`), and a K = 5 A/A calibration at
the link's ``from_ref`` whose every margin is 0.03.

The A/A fixtures the registry, citation and chain tests share live here too:
the five-session calibration over ``w`` and ``v`` (`W_PAIRS`, `V_PAIRS`,
`_aa_sessions`), a registry over it (`_registry`, `_written`), and the hand-edit
helpers that bypass `write_registry` (`_rewrite`, `_dump_unchecked`, `_dirty`).
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
from collections.abc import Callable
from collections.abc import Mapping
from collections.abc import Sequence
from pathlib import Path
from typing import Any
from typing import Optional

from tools.perf_history_store import REPOSITORY
from tools.perf_ratio import pair_order
from tools.perf_ratio_gate_registry import CALIBRATION_SESSIONS
from tools.perf_ratio_gate_registry import AnchorLink
from tools.perf_ratio_gate_registry import Calibration
from tools.perf_ratio_gate_registry import CitedRecord
from tools.perf_ratio_gate_registry import GateRegistry
from tools.perf_ratio_gate_registry import build_registry
from tools.perf_ratio_gate_registry import calibration_dir
from tools.perf_ratio_gate_registry import calibration_of
from tools.perf_ratio_gate_registry import promotion_dir
from tools.perf_ratio_gate_registry import write_registry
from tools.perf_ratio_record import ArmIdentity
from tools.perf_ratio_record import ArmName
from tools.perf_ratio_record import CpuSeconds
from tools.perf_ratio_record import Injection
from tools.perf_ratio_record import LoadOneMinute
from tools.perf_ratio_record import PairRun
from tools.perf_ratio_record import PerfRatioSession
from tools.perf_ratio_record import WorkloadRatio
from tools.run_provenance import sha256_file

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
# The one (ref_cpu, head_cpu) pair of each workload in each of a link's K = 5 A/A sessions: ratios
# 1.0, 1.015, 0.99, 0.995, 1.005 (distinct, so no two records coincide), so max |rho - 1| = 0.015 and
# SAFETY 2 gives m_w = 0.03 (to the float).
LINK_AA_PAIRS = ((100.0, 100.0), (100.0, 101.5), (100.0, 99.0), (100.0, 99.5), (100.0, 100.5))
# The record file name inside every session directory.
RECORD_NAME = "perf_ratio_session.json"


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


def _under(root: Path, directory: Path) -> Path:
    # A repository directory's counterpart under a test root.
    return root / directory.relative_to(REPOSITORY)


def cite(root: Path, path: Path) -> CitedRecord:
    """*path* cited relative to *root*, with its file's SHA-256."""
    return CitedRecord(path=path.relative_to(root).as_posix(), record_sha256=sha256_file(path))


def link_calibration(root: Path, ref: str, workloads: Sequence[str], *, harness: str = HARNESS, cpu_model: str = CPU_MODEL) -> Calibration:
    """The A/A calibration at *ref* over *workloads* (every margin 0.03), its records written under *root*."""
    directory = _under(root, calibration_dir(ref, cpu_model))
    paths = [
        write_record(
            synthetic_session({name: (pair,) for name in workloads}, ref_commit=ref, head_rev=ref, harness=harness, cpu_model=cpu_model),
            directory / f"aa-{index}" / RECORD_NAME,
        )
        for index, pair in enumerate(LINK_AA_PAIRS)
    ]
    return calibration_of(paths, ref_commit=ref, harness_sha256=harness, root=root)


def write_promotion(root: Path, from_ref: str, to_ref: str, gate: PerfRatioSession, extensions: Sequence[PerfRatioSession] = ()) -> tuple[CitedRecord, ...]:
    """Write a promotion's gate session and its extensions in the promotion layout; their citations, in merge order."""
    directory = _under(root, promotion_dir(from_ref, to_ref))
    paths = [write_record(gate, directory / "gate" / RECORD_NAME)]
    paths += [write_record(extension, directory / f"ext-{index}" / RECORD_NAME) for index, extension in enumerate(extensions, start=1)]
    return tuple(cite(root, path) for path in paths)


def promotion_link(
    root: Path,
    from_ref: str,
    to_ref: str,
    gate: Mapping[str, Sequence[tuple[float, float]]],
    extensions: Sequence[tuple[str, Sequence[tuple[float, float]]]] = (),
) -> AnchorLink:
    """An honest link *from_ref* -> *to_ref*: its gate session runs *gate*, each extension one workload's pairs.

    Args:
        root: The registry root every record is written under.
        from_ref: The reference the promoted session measured against.
        to_ref: The promoted HEAD.
        gate: The gate session's pairs per workload, each ``(ref_cpu, head_cpu)``.
        extensions: ``(workload, pairs)`` per extension session, in merge order.

    Returns:
        The link, citing every record it wrote and the A/A calibration at *from_ref* over *gate*'s workloads.
    """
    sessions = [synthetic_session(gate, ref_commit=from_ref, head_rev=to_ref)]
    sessions += [synthetic_session({name: pairs}, ref_commit=from_ref, head_rev=to_ref) for name, pairs in extensions]
    cited = write_promotion(root, from_ref, to_ref, sessions[0], sessions[1:])
    return AnchorLink(from_ref=from_ref, to_ref=to_ref, sessions=cited, calibration=link_calibration(root, from_ref, sorted(gate)))


# (ref_cpu, head_cpu) of each A/A session's single "w" pair: ratios 1.00, 0.98, 1.01, 1.00, 1.005.
W_PAIRS = ((100.0, 100.0), (100.0, 98.0), (100.0, 101.0), (150.0, 150.0), (200.0, 201.0))
# Pairs of the second workload "v", the same in every session (two pairs, ratio 0.99).
V_PAIRS = ((50.0, 49.5), (51.0, 52.0))


def _ratios(index: int) -> dict[str, Sequence[tuple[float, float]]]:
    return {"w": (W_PAIRS[index],), "v": V_PAIRS}


def _aa_sessions(count: int = CALIBRATION_SESSIONS) -> list[PerfRatioSession]:
    return [synthetic_session(_ratios(index)) for index in range(count)]


def _write_all(root: Path, sessions: Sequence[PerfRatioSession], cpu: str = "apple-m1-max") -> list[Path]:
    return [write_record(session, root / cpu / f"aa-{index}" / "perf_ratio_session.json") for index, session in enumerate(sessions)]


def _registry(root: Path) -> GateRegistry:
    calibration = calibration_of(_write_all(root, _aa_sessions()), ref_commit=REF, harness_sha256=HARNESS, root=root)
    return build_registry(ref_commit=REF, harness_sha256=HARNESS, calibrations=[calibration])


def _written(root: Path) -> Path:
    path = root / "perf_ratio_gate.json"
    write_registry(_registry(root), path, root=root)
    return path


def _content_digest(document: Mapping[str, Any]) -> str:
    content = {key: value for key, value in document.items() if key != "registry_sha256"}
    return hashlib.sha256(json.dumps(content, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def _rewrite(path: Path, edit: Callable[[dict[str, Any]], None], *, redigest: bool) -> None:
    document = json.loads(path.read_text(encoding="utf-8"))
    edit(document)
    if redigest:
        document["registry_sha256"] = _content_digest(document)
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _dump_unchecked(registry: GateRegistry, path: Path) -> Path:
    """Write *registry* bypassing `write_registry`'s verification, as a hand edit would."""
    path.write_text(registry.dumps(), encoding="utf-8")
    return path


def _dirty(session: PerfRatioSession) -> PerfRatioSession:
    head = dataclasses.replace(session.arms[ArmName.HEAD], dirty=True, patch_sha256="p" * 64)
    return dataclasses.replace(session, arms={ArmName.REF: session.arms[ArmName.REF], ArmName.HEAD: head})
