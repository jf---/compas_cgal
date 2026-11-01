"""The ratio harness's workloads: registry, each arm child's command, and parsing its result line.

Every workload's child reports its own in-child process CPU. Polygon and Ω*
generation use ``tools.held_generation_clock``'s worker form, the production
certificate uses ``tools.coverage_certificate_budget_gate --worker``, and the
GT pipeline cases run the arm's own ``tools.pipeline_timing`` under a
harness-owned CPU clock (`GT_CPU_CLOCK_SNIPPET`). The result-line markers are
imported from those tools, never copied, so a renamed marker fails loud.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping
from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum
from typing import Any
from typing import NamedTuple

from tools.coverage_certificate_budget_gate import BUDGET_WORKER_MARKER
from tools.held_generation_clock import WORKER_MARKER
from tools.held_generation_clock import GenerationClock
from tools.perf_ratio_record import CpuSeconds
from tools.pipeline_timing import RESULT_LINE_PREFIX


class UnknownWorkloadError(Exception):
    """A workload that no workload set registers."""


class InjectionUnsupportedError(Exception):
    """Injected CPU was requested for a family whose clock the harness does not own (only GT)."""


class MissingChildResultError(Exception):
    """A child's stdout holds zero result lines, or more than one."""


class ChildWorkloadMismatchError(Exception):
    """A child reports a result for a different workload than the one it was launched for."""


class MalformedChildResultError(Exception):
    """A child's result line is not a JSON object, lacks a required key, or holds a wrong-typed value."""


class Family(str, Enum):
    """Which tool's clock measures a workload."""

    POLYGON = "polygon"
    OMEGA = "omega"
    CERTIFICATE = "certificate"
    GT = "gt"


@dataclass(frozen=True)
class Workload:
    """One measured workload.

    Attributes:
        family: The family whose child measures it.
        argument: The child's workload argument (a held workload or a gate case key).
    """

    family: Family
    argument: str

    @property
    def name(self) -> str:
        """The workload's name, ``<family>/<argument>``."""
        return f"{self.family.value}/{self.argument}"


HELD_WORKLOADS = ("figure5_refined", "figure8_crossed_skis_boundary", "figure8_monstera_boundary_qualified_corners")
CERTIFICATE_WORKLOAD = "figure8_monstera_boundary_qualified_corners"
GT_CAP120_KEYS = ("cap-120-default/pipeline/rect_12x8", "cap-120-default/pipeline/rect_20x12", "cap-120-default/pipeline/L_shape")
GT_CAP40_KEYS = ("cap-40-attribution/pipeline/rect_12x8", "cap-40-attribution/pipeline/rect_20x12", "cap-40-attribution/pipeline/L_shape")

DEFAULT_SET_NAME = "default"
GT_CAP40_SET_NAME = "gt-cap40"
WORKLOAD_SETS: Mapping[str, tuple[Workload, ...]] = {
    DEFAULT_SET_NAME: (
        *(Workload(Family.POLYGON, w) for w in HELD_WORKLOADS),
        *(Workload(Family.OMEGA, w) for w in HELD_WORKLOADS),
        Workload(Family.CERTIFICATE, CERTIFICATE_WORKLOAD),
        *(Workload(Family.GT, k) for k in GT_CAP120_KEYS),
    ),
    GT_CAP40_SET_NAME: tuple(Workload(Family.GT, k) for k in GT_CAP40_KEYS),
}

# The certificate row's deterministic fields; its timing and load fields are excluded from the digest.
CERTIFICATE_DIGEST_FIELDS = (
    "workload",
    "circles",
    "connectors",
    "sub_pieces",
    "omega_pieces",
    "vertex_disks",
    "repair_inserted",
    "verdict",
    "witnesses",
    "witnesses_by_kind",
    "undecided_edges",
    "regions",
    "pieces",
    "edges",
)

# Clock reads `tools.pipeline_timing.run_case_here` makes (two `perf_counter` calls around the
# generator); the literal 2 in `GT_CPU_CLOCK_SNIPPET` is this value, kept verbatim in the snippet.
GT_CLOCK_READS = 2

# Runs in an arm's interpreter: argv = [-c, case key, injected CPU seconds]. Imports THAT arm's
# tools.pipeline_timing, gives the module a CPU stand-in for `time`, and spins the injected CPU
# inside the measured span (before the end read returns).
GT_CPU_CLOCK_SNIPPET = """
import sys, time
import tools.pipeline_timing as timed
key, inject = sys.argv[1], float(sys.argv[2])
class _CpuClock:
    reads = 0
    def perf_counter(self):
        _CpuClock.reads += 1
        if _CpuClock.reads == 2 and inject > 0.0:
            until = time.process_time() + inject
            while time.process_time() < until:
                pass
        return time.process_time()
vars(timed)["time"] = _CpuClock()
status = timed.main(["--case", key])
if _CpuClock.reads != 2:
    raise SystemExit(f"tools.pipeline_timing read its clock {_CpuClock.reads} times, not 2: move the reference")
sys.exit(status)
"""

_MARKERS: Mapping[Family, str] = {
    Family.POLYGON: WORKER_MARKER,
    Family.OMEGA: WORKER_MARKER,
    Family.CERTIFICATE: BUDGET_WORKER_MARKER,
    Family.GT: RESULT_LINE_PREFIX,
}


class ChildRun(NamedTuple):
    """One child's measurement.

    Attributes:
        cpu_seconds: The child's in-child process CPU for the workload.
        digest: The child's output digest (determinism and same-output witness).
    """

    cpu_seconds: CpuSeconds
    digest: str


def workload_named(name: str) -> Workload:
    """The registered workload called *name*, searching every set.

    Raises:
        UnknownWorkloadError: No set registers *name*.
    """
    for workloads in WORKLOAD_SETS.values():
        for workload in workloads:
            if workload.name == name:
                return workload
    raise UnknownWorkloadError(f"{name!r} is not a registered workload.")


def child_command(workload: Workload, python: Sequence[str], inject: CpuSeconds) -> tuple[str, ...]:
    """The command that measures *workload* once in a fresh child of the arm interpreter *python*.

    Args:
        workload: A registered workload.
        python: The arm interpreter's command prefix.
        inject: Process CPU seconds to spin inside the measured span; positive only for GT.

    Returns:
        The child's argv.

    Raises:
        UnknownWorkloadError: *workload* is not registered (the certificate worker ignores any
            argument, so an unregistered one would be measured under the wrong name).
        InjectionUnsupportedError: *inject* is positive for a family other than GT.
    """
    workload_named(workload.name)
    if workload.family is not Family.GT and inject > 0.0:
        raise InjectionUnsupportedError(f"{workload.name}: CPU injection needs a harness-owned clock, which only the GT family has.")
    if workload.family in (Family.POLYGON, Family.OMEGA):
        return (*python, "-m", "tools.held_generation_clock", "--domain", workload.family.value, workload.argument)
    if workload.family is Family.CERTIFICATE:
        return (*python, "-m", "tools.coverage_certificate_budget_gate", "--worker")
    return (*python, "-c", GT_CPU_CLOCK_SNIPPET, workload.argument, repr(float(inject)))


def parse_child(workload: Workload, stdout: str) -> ChildRun:
    """The measurement in a child's stdout: exactly one line starting with the family's marker.

    Args:
        workload: The workload the child measured.
        stdout: The child's whole stdout; other lines (log noise) are ignored.

    Returns:
        The child's CPU and output digest.

    Raises:
        MissingChildResultError: Zero or more than one result line.
        MalformedChildResultError: The result line is not a JSON object, lacks a required key, or holds a
            wrong-typed value (CPU not a finite number, digest not the expected object or string).
        ChildWorkloadMismatchError: A certificate row names a workload other than *workload*'s argument.
    """
    marker = _MARKERS[workload.family]
    lines = [line for line in stdout.splitlines() if line.startswith(marker)]
    if len(lines) != 1:
        raise MissingChildResultError(f"{workload.name}: expected one {marker.strip()} line, found {len(lines)}.")
    where = workload.name
    try:
        decoded: object = json.loads(lines[0][len(marker) :])
    except json.JSONDecodeError as error:
        raise MalformedChildResultError(f"{where}: the result line is not JSON: {error}") from error
    raw = _object(decoded, f"{where} result")
    if workload.family in (Family.POLYGON, Family.OMEGA):
        clock_raw = _object(_key(raw, "clock", where), f"{where} clock")
        try:
            clock = GenerationClock.from_json(json.dumps(clock_raw))
        except KeyError as error:
            raise MalformedChildResultError(f"{where}: the clock lacks {error}.") from error
        except (TypeError, ValueError, AttributeError) as error:
            raise MalformedChildResultError(f"{where}: the clock holds a wrong-typed value: {error}") from error
        digest = _object(_key(raw, "digest", where), f"{where} digest")
        return ChildRun(_cpu(clock.cpu_seconds, f"{where} clock CPU"), _sha256_json(digest))
    if workload.family is Family.CERTIFICATE:
        # The worker ignores its argument and certifies its own tree's production workload.
        certified = _key(raw, "workload", where)
        if certified != workload.argument:
            raise ChildWorkloadMismatchError(f"{where}: the child certified {certified!r}, not {workload.argument!r}.")
        cpu = _cpu(_key(raw, "certificate_cpu_seconds", where), f"{where} certificate_cpu_seconds")
        return ChildRun(cpu, _sha256_json({k: _key(raw, k, where) for k in CERTIFICATE_DIGEST_FIELDS}))
    cpu = _cpu(_key(raw, "generation_s", where), f"{where} generation_s")
    digest_text = _key(raw, "digest", where)
    if not isinstance(digest_text, str):
        raise MalformedChildResultError(f"{where} digest is a string, not {digest_text!r}.")
    operations = _key(raw, "operations", where)
    if isinstance(operations, bool) or not isinstance(operations, int):
        raise MalformedChildResultError(f"{where} operations is an integer, not {operations!r}.")
    return ChildRun(cpu, f"{digest_text}:{operations}")


def _object(value: object, where: str) -> Mapping[str, Any]:
    if not isinstance(value, dict):
        raise MalformedChildResultError(f"{where} is a JSON object, not {value!r}.")
    return value


def _key(record: Mapping[str, Any], key: str, where: str) -> Any:
    try:
        return record[key]
    except KeyError as error:
        raise MalformedChildResultError(f"{where}: the result lacks {key!r}.") from error


def _cpu(value: object, where: str) -> CpuSeconds:
    # bool is an int subclass; a JSON true is not a CPU reading.
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise MalformedChildResultError(f"{where} is a finite number, not {value!r}.")
    return CpuSeconds(float(value))


def _sha256_json(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode("utf-8")).hexdigest()
