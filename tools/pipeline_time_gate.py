"""The Task 1.2 generation-time gate (ruling D6), generated from three records.

Per ``pipeline`` gate case the gate requires the 1.2 serial median
(`benchmarks/results/pipeline_timing.json`) to be at most the slowest of the
1.1 tree's serial runs plus the case's median `C_pencil`
(`benchmarks/results/pipeline_pencil_cost.json`, `tools.pipeline_pencil_cost`).

The 1.1 base is the `pipeline_timing.json` committed with the Task 1.1 tree,
`COMMITTED_1_1_REV`, read through ``git show``. Another revision is used only
when named explicitly (``--base-rev``), and a same-session re-timing of the 1.1
tree only when passed explicitly (``--base-session``). The page names the
source, its SHA-256 and its load.

A verdict counts only from quiet-machine evidence measured back to back
(user ruling D2, 2026-10-22; orchestrator rulings on RV25 O1 and RV27 I2), and
the gate checks that from the records' load and `tools.run_provenance`
provenance, never from a flag:

- the 1.1 base is a same-session re-timing (`tools.pipeline_timing_1_1_base`):
  a loaded base inflates its slowest run and loosens the bound, the unsafe
  direction (`CommittedBaseError`);
- every timing record carries load and provenance (`UnrecordedLoadError`);
- the base's commit is the expected 1.1 commit (`BaseRevisionMismatchError`)
  and the base is unmodified (`DirtyBaseError`);
- the base was written by the wrapper: its tool digest is the one recomputed
  from the 1.1 commit's timing sources and this tree's wrapper sources
  (`BaseToolMismatchError`);
- the base timed the 1.1 library: when the measured native sources differ
  (the records' ``native_sources_sha256``, read from the trees actually
  measured, never from a commit id), the base's native identity
  (``native_sha256``, every extension) must differ from the 1.2 records'
  (`ForeignBaseLibraryError`);
- the 1.2 timing and the pencil record describe the same tree
  (`PencilTreeMismatchError`) and ran the same native build
  (`NativeLibraryMismatchError`). The base runs its own build: each record's
  native build is held unchanged within its own run
  (``tools.run_native_identity``), and the base's is compared with the 1.2
  records' only by the foreign-library check above;
- all three ran on one machine (`HostMismatchError`), in one pixi environment
  and Python version (`EnvironmentMismatchError`), with one ``pixi.lock``
  (`LockMismatchError`);
- they ran in order, base, 1.2 timing, pencil, each starting within
  `BACK_TO_BACK_WINDOW_S` of the previous one's end (`NotBackToBackError`);
- every load is within `QUIET_LOAD_AVERAGE_MAX` and the pencil record is not
  indicative (`MachineNotQuietError`).

When several fail, the first in that order is raised, with every reason in
its message. ``--indicative`` renders the comparison anyway, labelled as no
verdict and naming the reasons; with the committed base that is the only
render besides pending.
Without a pencil record (``--pending``) the page renders the inputs it has and
states why the verdict is pending.

``python -m tools.pipeline_time_gate`` writes
`benchmarks/results/pipeline_time_gate.md`, which `docs/refactor/stage1.md`
includes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional
from typing import Sequence

from tools.gate_case_keys import CASE_KEY_SEPARATOR
from tools.held_generation_clock import QUIET_LOAD_AVERAGE_MAX
from tools.held_generation_clock import MachineNotQuietError
from tools.pipeline_pencil_cost_results import PENCIL_GENERATOR
from tools.pipeline_pencil_cost_results import PIPELINE_PENCIL_COST_PATH
from tools.pipeline_pencil_cost_results import PipelinePencilCosts
from tools.pipeline_pencil_cost_results import is_pencil_key
from tools.pipeline_pencil_cost_results import load_pipeline_pencil_costs
from tools.pipeline_pencil_cost_results import rooted
from tools.pipeline_timing_1_1_base import expected_base_tool_sha256
from tools.pipeline_timing_results import LOAD_DECIMALS
from tools.pipeline_timing_results import PIPELINE_TIMING_PATH
from tools.pipeline_timing_results import SCHEMA_WITHOUT_LOAD
from tools.pipeline_timing_results import PipelineTimings
from tools.pipeline_timing_results import format_load
from tools.run_provenance import RunProvenance

ROOT = Path(__file__).resolve().parents[1]
PIPELINE_TIME_GATE_TABLE_PATH = Path("benchmarks") / "results" / "pipeline_time_gate.md"
# The commit whose `pipeline_timing.json` is the Task 1.1 tree's serial timing:
# the last commit touching that file before any Task 1.2 pipeline source commit.
COMMITTED_1_1_REV = "7b001d47"
# Policy, not a derivation: the longest end-to-start gap, in seconds, still called one session
# (base -> 1.2 timing -> pencil). Each run's own quietness is checked separately
# (`QUIET_LOAD_AVERAGE_MAX`); this bounds drift of machine state (thermal and frequency state, background
# jobs, OS updates) between the compared runs. 15 min = one 15-minute load-average time constant.
BACK_TO_BACK_WINDOW_S = 15 * 60
# Rendered precision of seconds: milliseconds (as `tools.pipeline_timing_table`).
RENDERED_SECONDS_DECIMALS = 3
# What to do instead of a refused verdict, named in the refusal.
QUIET_RERUN_HINT = (
    "On a quiet machine, re-time the 1.1 tree (pass it with --base-session), then run `pixi run pipeline-timing` and "
    "`pixi run pipeline-pencil-cost`, back to back; or render without a verdict: python -m tools.pipeline_time_gate --indicative."
)


class StalePencilCostError(Exception):
    """The pencil record was not measured beside the timing document the gate reads (SHA-256 or digest differ)."""


class IncompleteTimeGateInputError(Exception):
    """A ``pipeline`` case the gate judges is missing from the 1.1 base or the 1.2 timing document."""


class UnrecordedLoadError(Exception):
    """A timing record (the 1.2 run or a same-session 1.1 base) holds no load average, so its quietness is unknown: no verdict."""


class CommittedBaseError(Exception):
    """The 1.1 base is a committed record, not a same-session re-timing on the quiet machine: no verdict."""


class BaseRevisionMismatchError(Exception):
    """The same-session base was measured on another commit than the expected 1.1 tree: no verdict."""


class DirtyBaseError(Exception):
    """The same-session base's tree differed from its commit: the 1.1 tree must be timed unmodified."""


class BaseToolMismatchError(Exception):
    """The same-session base's tool digest is not the wrapper's over the 1.1 commit's timing sources: another harness wrote it."""


class ForeignBaseLibraryError(Exception):
    """The base's native identity equals the 1.2 records' although the measured native sources differ: the base timed the 1.2 library."""


class NativeLibraryMismatchError(Exception):
    """The pencil record and the 1.2 timing record ran different native builds (``native_sha256``): no verdict."""


class EnvironmentMismatchError(Exception):
    """The three records were not measured in one pixi environment and Python version: no verdict."""


class PencilTreeMismatchError(Exception):
    """The pencil record and the 1.2 timing record describe different trees (commit or difference digest): no verdict."""


class HostMismatchError(Exception):
    """The three records were not measured on one machine (host name or CPU model differ): no verdict."""


class LockMismatchError(Exception):
    """The three records were not measured with one ``pixi.lock``: no verdict."""


class NotBackToBackError(Exception):
    """The records did not run in order base, 1.2 timing, pencil, each within `BACK_TO_BACK_WINDOW_S` of the previous end: no verdict."""


class SessionBaseOutsideRepositoryError(Exception):
    """A same-session base file lies outside the repository, so the page could not name it reproducibly."""


class TimeGateMode(Enum):
    """What a rendered time gate is.

    VERDICT: from quiet-machine timing and pencil records. INDICATIVE: rendered
    on request from evidence that is loaded, indicative or of unrecorded load;
    no verdict. PENDING: no pencil record yet.
    """

    VERDICT = "verdict"
    INDICATIVE = "indicative"
    PENDING = "pending"


@dataclass(frozen=True)
class TimeGateBase:
    """The 1.1 tree's serial timing the gate compares against, with its provenance.

    Attributes:
        source: Where it came from, as rendered: a committed record at a named
            revision, or a same-session measurement at a named path.
        sha256: SHA-256 of the record's bytes.
        timings: The validated timing document.
        same_session: True for a re-timing measured in this session, which a verdict requires.
    """

    source: str
    sha256: str
    timings: PipelineTimings
    same_session: bool


def committed_time_gate_base(rev: str = COMMITTED_1_1_REV) -> TimeGateBase:
    """The 1.1 base from the `pipeline_timing.json` committed at *rev* (``git show``), `COMMITTED_1_1_REV` by default.

    Raises:
        subprocess.CalledProcessError: *rev* holds no such file.
        MalformedPipelineTimingsError: The committed record breaks its schema.
    """
    spec = f"{rev}:{PIPELINE_TIMING_PATH.as_posix()}"
    data = subprocess.run(["git", "show", spec], cwd=ROOT, capture_output=True, check=True).stdout
    return TimeGateBase(source=f"committed record `{spec}`", sha256=hashlib.sha256(data).hexdigest(), timings=PipelineTimings.build(json.loads(data)), same_session=False)


def session_time_gate_base(path: Path, root: Path = ROOT) -> TimeGateBase:
    """The 1.1 base from a `pipeline-timing` document of the 1.1 tree measured in this session, at *path* inside *root*.

    The page names it by its path relative to *root*.

    Raises:
        SessionBaseOutsideRepositoryError: *path* is not inside *root*.
        MalformedPipelineTimingsError: The document breaks its schema.
    """
    try:
        relative = path.resolve().relative_to(root.resolve())
    except ValueError as error:
        raise SessionBaseOutsideRepositoryError(f"{path} is outside {root}: commit the base record inside the repository.") from error
    data = path.read_bytes()
    return TimeGateBase(
        source=f"same-session measurement `{relative.as_posix()}`",
        sha256=hashlib.sha256(data).hexdigest(),
        timings=PipelineTimings.build(json.loads(data)),
        same_session=True,
    )


@dataclass(frozen=True)
class CaseTimeGate:
    """One ``pipeline`` case's generation-time comparison.

    Attributes:
        key: The case key.
        median_s: The 1.2 serial median (uninstrumented).
        base_slowest_s: The slowest of the 1.1 tree's serial runs.
        c_pencil_s: The median `C_pencil`, or None while pending.
    """

    key: str
    median_s: float
    base_slowest_s: float
    c_pencil_s: Optional[float]

    @property
    def bound_s(self) -> Optional[float]:
        """The slowest 1.1 run plus `C_pencil`, or None while pending."""
        return None if self.c_pencil_s is None else self.base_slowest_s + self.c_pencil_s

    @property
    def holds(self) -> Optional[bool]:
        """Whether the 1.2 median is within the bound (equality holds), or None while pending."""
        bound = self.bound_s
        return None if bound is None else self.median_s <= bound


def _timing_load(timings: PipelineTimings) -> str:
    if timings.load_start is None or timings.load_end is None:
        return f"unrecorded (schema `{SCHEMA_WITHOUT_LOAD}`)"
    return f"{format_load(timings.load_start)} at the start, {format_load(timings.load_end)} at the end"


def _timing_breach(label: str, timings: PipelineTimings) -> Optional[str]:
    """Why *timings* is not quiet-machine evidence (`PipelineTimings.quiet`), or None when it is."""
    peak, end_fifteen = timings.max_one_minute_load, timings.end_fifteen_minute_load
    if peak is None or end_fifteen is None:
        return f"{label} holds no load average (schema `{SCHEMA_WITHOUT_LOAD}`)"
    if peak > QUIET_LOAD_AVERAGE_MAX:
        return f"{label}'s 1-minute load reached {peak:.{LOAD_DECIMALS}f}"
    if end_fifteen > QUIET_LOAD_AVERAGE_MAX:
        return f"{label}'s 15-minute load at the end reached {end_fifteen:.{LOAD_DECIMALS}f}"
    return None


Breach = tuple[type[Exception], str]


def _tree(record: RunProvenance) -> tuple[str, bool, Optional[str]]:
    return record.rev, record.dirty, record.patch_sha256


def _gap_breach(earlier: str, earlier_record: RunProvenance, later: str, later_record: RunProvenance) -> Optional[str]:
    gap_s = (later_record.started - earlier_record.ended).total_seconds()
    if 0 <= gap_s <= BACK_TO_BACK_WINDOW_S:
        return None
    return f"{later} started {gap_s:+.0f} s after {earlier} ended ({earlier_record.ended_utc} -> {later_record.started_utc}; allowed 0 to {BACK_TO_BACK_WINDOW_S} s)"


def evidence_breaches(base: TimeGateBase, timings: PipelineTimings, costs: PipelinePencilCosts, expected_base_rev: str, expected_base_tool_sha256: str) -> tuple[Breach, ...]:
    """Why *base*, *timings* and *costs* are not back-to-back, quiet, same-machine evidence of the expected trees, in raise order; empty when they are.

    Args:
        base: The 1.1 base.
        timings: The 1.2 timing record.
        costs: The pencil record.
        expected_base_rev: The full commit id the base must have been measured on.
        expected_base_tool_sha256: The base's tool digest as the wrapper writes it for that commit.

    Returns:
        (error class, reason) per failed check, in the order the module docstring lists.
    """
    breaches: list[Breach] = []
    if not base.same_session:
        breaches.append((CommittedBaseError, f"the 1.1 base is a {base.source}, not a same-session re-timing (--base-session)"))
    for label, record in (("the 1.1 base", base.timings), ("the 1.2 timing record", timings)):
        if record.provenance is None or record.max_one_minute_load is None:
            breaches.append((UnrecordedLoadError, f"{label} holds no load or provenance (schema `{SCHEMA_WITHOUT_LOAD}`)"))
    base_origin, timing_origin, pencil_origin = base.timings.provenance, timings.provenance, costs.provenance
    if base_origin is not None and base_origin.rev != expected_base_rev:
        breaches.append((BaseRevisionMismatchError, f"the 1.1 base was measured on {base_origin.rev}, not {expected_base_rev}"))
    if base_origin is not None and base_origin.dirty:
        breaches.append((DirtyBaseError, f"the 1.1 base's tree differed from {base_origin.rev} (difference SHA-256 {base_origin.patch_sha256})"))
    if base_origin is not None and base_origin.timing_tool_sha256 != expected_base_tool_sha256:
        breaches.append((BaseToolMismatchError, f"the 1.1 base's tool digest {base_origin.timing_tool_sha256} is not the wrapper's {expected_base_tool_sha256}"))
    head_origins = (pencil_origin,) if timing_origin is None else (timing_origin, pencil_origin)
    if base_origin is not None:
        sources_differ = any(base_origin.native_sources_sha256 != origin.native_sources_sha256 for origin in head_origins)
        if sources_differ and base_origin.native_sha256 in {origin.native_sha256 for origin in head_origins}:
            breaches.append(
                (ForeignBaseLibraryError, f"the 1.1 base's native extensions {base_origin.native_sha256} are the 1.2 records', although the measured native sources differ")
            )
    if timing_origin is not None and _tree(timing_origin) != _tree(pencil_origin):
        breaches.append((PencilTreeMismatchError, f"the 1.2 timing tree {_tree(timing_origin)} is not the pencil tree {_tree(pencil_origin)}"))
    if timing_origin is not None and timing_origin.native_sha256 != pencil_origin.native_sha256:
        breaches.append((NativeLibraryMismatchError, f"the 1.2 timing ran native build {timing_origin.native_sha256}, the pencil run {pencil_origin.native_sha256}"))
    origins = {name: origin for name, origin in (("1.1 base", base_origin), ("1.2 timing", timing_origin), ("pencil", pencil_origin)) if origin is not None}
    machines = {f"{origin.host} ({origin.cpu_model})" for origin in origins.values()}
    if len(machines) > 1:
        breaches.append((HostMismatchError, "the records ran on " + ", ".join(f"{name}: {origin.host} ({origin.cpu_model})" for name, origin in origins.items())))
    if len({(origin.environment, origin.python_version) for origin in origins.values()}) > 1:
        breaches.append(
            (EnvironmentMismatchError, "the records ran in " + ", ".join(f"{name}: {origin.environment} / Python {origin.python_version}" for name, origin in origins.items()))
        )
    if len({origin.pixi_lock_sha256 for origin in origins.values()}) > 1:
        breaches.append((LockMismatchError, "the records ran with different pixi.lock: " + ", ".join(f"{name}: {origin.pixi_lock_sha256}" for name, origin in origins.items())))
    if base_origin is not None and timing_origin is not None:
        gaps = (_gap_breach("the 1.1 base", base_origin, "the 1.2 timing", timing_origin), _gap_breach("the 1.2 timing", timing_origin, "the pencil run", pencil_origin))
        breaches.extend((NotBackToBackError, gap) for gap in gaps if gap is not None)
    for label, record in (("the 1.1 base", base.timings), ("the 1.2 timing record", timings)):
        breach = _timing_breach(label, record)
        if breach is not None and record.max_one_minute_load is not None:
            breaches.append((MachineNotQuietError, breach))
    if costs.indicative:
        breaches.append((MachineNotQuietError, "the pencil record is indicative (measured without the quiet-machine check)"))
    if costs.max_one_minute_load > QUIET_LOAD_AVERAGE_MAX:
        breaches.append((MachineNotQuietError, f"the pencil record's 1-minute load reached {costs.max_one_minute_load:.{LOAD_DECIMALS}f}"))
    if costs.max_fifteen_minute_load > QUIET_LOAD_AVERAGE_MAX:
        breaches.append((MachineNotQuietError, f"the pencil record's 15-minute load reached {costs.max_fifteen_minute_load:.{LOAD_DECIMALS}f}"))
    return tuple(breaches)


@dataclass(frozen=True)
class PipelineTimeGate:
    """The generation-time gate over every measured ``pipeline`` case.

    Raw construction is unchecked; `build` validates.

    Attributes:
        base: The 1.1 base.
        timings: The 1.2 timing document.
        timing_sha256: SHA-256 of the 1.2 timing document.
        costs: The pencil record, or None while pending.
        mode: Verdict, indicative or pending.
        breaches: Why the evidence admits no verdict (`evidence_breaches` reasons); empty for a verdict.
        cases: The per-case comparisons, in gate-case order.
    """

    base: TimeGateBase
    timings: PipelineTimings
    timing_sha256: str
    costs: Optional[PipelinePencilCosts]
    mode: TimeGateMode
    breaches: tuple[str, ...]
    cases: tuple[CaseTimeGate, ...]

    @classmethod
    def build(
        cls,
        base: TimeGateBase,
        timings: PipelineTimings,
        timing_sha256: str,
        costs: Optional[PipelinePencilCosts],
        *,
        expected_base_rev: str,
        expected_base_tool_sha256: str,
        indicative: bool,
    ) -> PipelineTimeGate:
        """Compare every measured ``pipeline`` case of *timings* against *base* plus *costs*.

        Args:
            base: The 1.1 base.
            timings: The 1.2 timing document.
            timing_sha256: SHA-256 of *timings*' bytes.
            costs: The pencil record, or None to render the pending state.
            expected_base_rev: The full commit id a same-session base must have been measured on.
            expected_base_tool_sha256: The tool digest a same-session base must carry
                (`tools.pipeline_timing_1_1_base.expected_base_tool_sha256`).
            indicative: Render evidence that admits no verdict, labelled as such.

        Returns:
            The gate.

        Raises:
            StalePencilCostError: *costs* was measured beside another timing
                document, or a case's digest differs from *timings*'.
            CommittedBaseError, UnrecordedLoadError, BaseRevisionMismatchError,
            DirtyBaseError, BaseToolMismatchError, ForeignBaseLibraryError,
            PencilTreeMismatchError, NativeLibraryMismatchError, HostMismatchError, EnvironmentMismatchError,
            LockMismatchError, NotBackToBackError, MachineNotQuietError: the
                first failed check of `evidence_breaches`, when *indicative* is False.
            IncompleteTimeGateInputError: A case is missing from *base* or *timings*.
        """
        breaches: tuple[str, ...] = ()
        if costs is None:
            keys = tuple(key for key in timings.cases if is_pencil_key(key))
            mode = TimeGateMode.PENDING
        else:
            if costs.timing_sha256 != timing_sha256:
                raise StalePencilCostError(f"the pencil record was measured beside timing document {costs.timing_sha256}, not {timing_sha256}.")
            keys = tuple(costs.cases)
            stale = [key for key in keys if key in timings.cases and timings.cases[key].digest != costs.cases[key].digest]
            if stale:
                raise StalePencilCostError(f"the pencil record's digests differ from the timing document's for {stale}.")
            found = evidence_breaches(base, timings, costs, expected_base_rev, expected_base_tool_sha256)
            breaches = tuple(reason for _error, reason in found)
            if found and not indicative:
                first_error = found[0][0]
                raise first_error(f"no verdict: {'; '.join(breaches)} (quiet limit {QUIET_LOAD_AVERAGE_MAX}). {QUIET_RERUN_HINT}")
            mode = TimeGateMode.INDICATIVE if indicative else TimeGateMode.VERDICT
        missing = [key for key in keys if key not in base.timings.cases or key not in timings.cases]
        if missing:
            raise IncompleteTimeGateInputError(f"no 1.1 base or 1.2 timing for {missing}.")
        cases = tuple(
            CaseTimeGate(
                key=key,
                median_s=timings.cases[key].median_s,
                base_slowest_s=max(base.timings.cases[key].runs_s),
                c_pencil_s=None if costs is None else costs.cases[key].median_decide_s,
            )
            for key in keys
        )
        return cls(base=base, timings=timings, timing_sha256=timing_sha256, costs=costs, mode=mode, breaches=breaches, cases=cases)


def _seconds(value: float) -> str:
    return f"{value:.{RENDERED_SECONDS_DECIMALS}f}"


def _optional_seconds(value: Optional[float]) -> str:
    return "—" if value is None else _seconds(value)


def _verdict_cell(case: CaseTimeGate, mode: TimeGateMode) -> str:
    if case.holds is None:
        return "pending: no `C_pencil` record"
    outcome = "holds" if case.holds else "fails"
    return f"{outcome} (indicative)" if mode is TimeGateMode.INDICATIVE else outcome


def _origin(name: str, origin: Optional[RunProvenance]) -> str:
    if origin is None:
        return f"{name}: unrecorded"
    tree = f"`{origin.rev}`" + (f" dirty (difference SHA-256 `{origin.patch_sha256}`)" if origin.dirty else " clean")
    return f"{name}: {tree}, {origin.host} ({origin.cpu_model}), {origin.started_utc} to {origin.ended_utc}, pixi.lock `{origin.pixi_lock_sha256}`"


def _provenance_line(gate: PipelineTimeGate) -> str:
    records = [_origin("1.1 base", gate.base.timings.provenance), _origin("1.2 timing", gate.timings.provenance)]
    if gate.costs is not None:
        records.append(_origin("pencil", gate.costs.provenance))
    return "Provenance: " + "; ".join(records) + "."


def _gate_status(gate: PipelineTimeGate) -> str:
    timing_load = f"1.1 base load: {_timing_load(gate.base.timings)}; 1.2 timing load: {_timing_load(gate.timings)} (load averages over 1 / 5 / 15 minutes)"
    if gate.costs is None:
        return (
            f"**Verdict pending**: no `C_pencil` record (`{PIPELINE_PENCIL_COST_PATH.as_posix()}`) yet. {timing_load}; a verdict also "
            f"needs a same-session 1.1 base and all three records within a load of {QUIET_LOAD_AVERAGE_MAX}."
        )
    pencil_load = f"pencil record's highest 1-minute load {gate.costs.max_one_minute_load:.{LOAD_DECIMALS}f}"
    if gate.mode is TimeGateMode.INDICATIVE:
        why = "; ".join(gate.breaches) if gate.breaches else "rendered with --indicative"
        return f"**INDICATIVE, not a verdict**: {why} (quiet limit {QUIET_LOAD_AVERAGE_MAX}). {timing_load}; {pencil_load}."
    loads = sorted({format_load(load) for case in gate.costs.cases.values() for load in (case.load_before, case.load_after)})
    return f"**Verdict** from quiet-machine records (quiet limit {QUIET_LOAD_AVERAGE_MAX}). {timing_load}; {pencil_load} (recorded loads: {'; '.join(loads)})."


def render_time_gate_table(gate: PipelineTimeGate) -> str:
    """The Markdown method statement, inputs, measured load and per-case table of *gate*.

    Args:
        gate: The gate.

    Returns:
        The Markdown text, ending in a newline.
    """
    pencil = "none" if gate.costs is None else f"`{PIPELINE_PENCIL_COST_PATH.as_posix()}`, measured beside timing SHA-256 `{gate.costs.timing_sha256}`"
    lines = [
        "<!-- Generated by tools.pipeline_time_gate from the 1.1 base, benchmarks/results/pipeline_timing.json and "
        "benchmarks/results/pipeline_pencil_cost.json. Do not edit; run `python -m tools.pipeline_time_gate`. -->",
        "",
        (
            f"Generation-time gate (Task 1.2, ruling D6): per `{PENCIL_GENERATOR}` case, the 1.2 serial median must be at most the "
            "slowest of the 1.1 tree's serial runs plus the case's median `C_pencil`. "
            f"1.1 base: {gate.base.source} (SHA-256 `{gate.base.sha256}`). "
            f"1.2 medians: `{PIPELINE_TIMING_PATH.as_posix()}` (SHA-256 `{gate.timing_sha256}`). "
            f"`C_pencil`: {pencil}. "
            "`C_pencil` is measured on the code it excuses: any slowdown inside the pencil's decisions, the cap machinery they "
            "share included, is credited automatically, so the gate bounds only the time outside them."
        ),
        "",
        _gate_status(gate),
        "",
        _provenance_line(gate),
        "",
        "| cap | pocket | 1.2 median (s) | 1.1 slowest run (s) | `C_pencil` (s) | bound (s) | verdict |",
        "| :--- | :--- | ---: | ---: | ---: | ---: | :--- |",
    ]
    for case in gate.cases:
        cap_id, _generator, pocket = case.key.split(CASE_KEY_SEPARATOR)
        lines.append(
            f"| {cap_id} | {pocket} | {_seconds(case.median_s)} | {_seconds(case.base_slowest_s)} | {_optional_seconds(case.c_pencil_s)} "
            f"| {_optional_seconds(case.bound_s)} | {_verdict_cell(case, gate.mode)} |"
        )
    return "\n".join(lines) + "\n"


def resolve_commit(rev: str) -> str:
    """The full commit id of *rev*.

    Raises:
        subprocess.CalledProcessError: *rev* names no commit.
    """
    return subprocess.run(["git", "rev-parse", "--verify", f"{rev}^{{commit}}"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Command-line entry: render the generation-time gate from the 1.1 base, the 1.2 timing and the pencil record.

    Relative paths resolve against the repository root.

    Args:
        argv: Arguments, defaulting to ``sys.argv[1:]``.

    Returns:
        The process exit status.
    """
    parser = argparse.ArgumentParser(prog="python -m tools.pipeline_time_gate", description="Render the Task 1.2 generation-time gate (D6).")
    parser.add_argument(
        "--base-rev", metavar="REV", default=COMMITTED_1_1_REV, help=f"the 1.1 commit (default {COMMITTED_1_1_REV}): the committed base's, and the one a session base must carry"
    )
    parser.add_argument("--base-session", type=Path, metavar="PATH", help="a same-session 1.1 base from tools.pipeline_timing_1_1_base, inside the repository")
    parser.add_argument("--timing", type=Path, default=PIPELINE_TIMING_PATH, help="the 1.2 timing document")
    parser.add_argument("--pencil-cost", type=Path, default=PIPELINE_PENCIL_COST_PATH, help="the pencil record")
    parser.add_argument("--pending", action="store_true", help="render without a pencil record: inputs only, verdict pending")
    parser.add_argument("--indicative", action="store_true", help="render evidence that is not quiet, labelled as no verdict")
    parser.add_argument("--write", type=Path, default=PIPELINE_TIME_GATE_TABLE_PATH)
    args = parser.parse_args(argv)
    expected_base_rev = resolve_commit(args.base_rev)
    base = session_time_gate_base(rooted(args.base_session)) if args.base_session is not None else committed_time_gate_base(args.base_rev)
    timing_bytes = rooted(args.timing).read_bytes()
    timings = PipelineTimings.build(json.loads(timing_bytes))
    costs = None if args.pending else load_pipeline_pencil_costs(rooted(args.pencil_cost))
    gate = PipelineTimeGate.build(
        base,
        timings,
        hashlib.sha256(timing_bytes).hexdigest(),
        costs,
        expected_base_rev=expected_base_rev,
        expected_base_tool_sha256=expected_base_tool_sha256(expected_base_rev),
        indicative=args.indicative,
    )
    write = rooted(args.write)
    write.parent.mkdir(parents=True, exist_ok=True)
    write.write_text(render_time_gate_table(gate))
    print(f"wrote {write}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
