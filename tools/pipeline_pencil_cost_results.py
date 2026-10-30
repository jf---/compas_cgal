"""The pencil decision-cost record and the generation-time gate it feeds (Task 1.2, ruling D6).

`benchmarks/results/pipeline_pencil_cost.json` holds, per ``pipeline`` gate
case, each of `RUNS_PER_CASE` fresh-process generations split into three wall
clocks: the whole instrumented generator call (``generation_s``), the pencil's
own decisions (``decide_s`` = `C_pencil`: `TangentPencil.build`,
`RampOut.build`, `RampOut.next_member`) and the first cut (``first_cut_s``:
`generate.enter_first_cut`, helix + ramp-out emission + root loop, which
contains ``decide_s``). It also records the call count of each wrapped
callable on its own, the number of ramp-out members, the operation count, the
ops digest (equal to the uninstrumented digest of the same case in
`benchmarks/results/pipeline_timing.json`, or the tool fails), the 1/5/15
minute load average before and after each case, whether the run was
indicative (measured without the quiet-machine check), the SHA-256 of that
timing document, binding the cost to the timing run it was measured beside,
and the run's `tools.run_provenance.RunProvenance` (tree, tool, environment,
machine, UTC interval), which the time gate checks against the timing records.
`tools.pipeline_pencil_cost` measures it.

The call counts are checked against the first cut's structure, so the record
says what `C_pencil` contains: one `TangentPencil.build`, one `RampOut.build`
and one `enter_first_cut` per generation, and one `RampOut.next_member` per
emitted member plus the one that returns the root.

The generation-time gate that consumes this record is `tools.pipeline_time_gate`.
"""

from __future__ import annotations

import json
import re
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from tools.gate_case_keys import CASE_KEY_SEPARATOR
from tools.gate_case_keys import is_case_key
from tools.gate_measure import case_key
from tools.gate_measure import gate_cases
from tools.held_generation_clock import QUIET_LOAD_AVERAGE_MAX
from tools.held_generation_clock import LoadAverage
from tools.pipeline_timing_results import DIGEST_HEX_CHARS
from tools.pipeline_timing_results import LOAD_AVERAGE_WINDOWS
from tools.pipeline_timing_results import LOAD_DECIMALS
from tools.pipeline_timing_results import RUNS_PER_CASE
from tools.pipeline_timing_results import format_load
from tools.pipeline_timing_results import parse_load_average
from tools.run_provenance import MalformedProvenanceError
from tools.run_provenance import RunProvenance

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "pipeline-pencil-cost/2"
PIPELINE_PENCIL_COST_PATH = Path("benchmarks") / "results" / "pipeline_pencil_cost.json"
PIPELINE_PENCIL_COST_TABLE_PATH = Path("benchmarks") / "results" / "pipeline_pencil_cost.md"
# The gate generator whose first cut runs the tangent pencil.
PENCIL_GENERATOR = "pipeline"
# The wrapped callables, as the record names them. The first three are the
# pencil's decisions (`C_pencil`); the last is the first cut that contains them.
TANGENT_PENCIL_BUILD = "TangentPencil.build"
RAMP_OUT_BUILD = "RampOut.build"
RAMP_OUT_NEXT_MEMBER = "RampOut.next_member"
ENTER_FIRST_CUT = "generate.enter_first_cut"
DECIDE_CALLABLES = (TANGENT_PENCIL_BUILD, RAMP_OUT_BUILD, RAMP_OUT_NEXT_MEMBER)
WRAPPED_CALLABLES = (*DECIDE_CALLABLES, ENTER_FIRST_CUT)
# First cuts per generation: `generate` enters the first chain once.
FIRST_CUTS_PER_GENERATION = 1
# `RampOut.next_member` calls beyond one per emitted member: the call that returns the root.
ROOT_NEXT_MEMBER_CALLS = 1
# Rendered precision of seconds: milliseconds (as `tools.pipeline_timing_table`).
RENDERED_SECONDS_DECIMALS = 3
# Hex characters of a full SHA-256 (a document's content address).
SHA256_HEX_CHARS = 64
_DIGEST_PATTERN = re.compile(rf"[0-9a-f]{{{DIGEST_HEX_CHARS}}}")
_SHA256_PATTERN = re.compile(rf"[0-9a-f]{{{SHA256_HEX_CHARS}}}")
_RUN_KEYS = frozenset({"generation_s", "decide_s", "first_cut_s"})
_CASE_KEYS = frozenset({"runs", "calls", "ramp_out_members", "operations", "digest", "load_before", "load_after"})
_DOCUMENT_KEYS = frozenset({"schema", "runs_per_case", "timing_sha256", "indicative", "cases", "excluded", "provenance"})


class MalformedPencilCostError(Exception):
    """The pencil-cost document breaks its schema."""


class CallCompositionError(Exception):
    """The wrapped call counts do not match the first cut's structure: `C_pencil` would not be what it says it is."""


@dataclass(frozen=True)
class PencilRun:
    """One instrumented generation's three wall clocks, in seconds.

    Attributes:
        generation_s: The whole instrumented generator call.
        decide_s: The pencil's decisions (`C_pencil`), summed over every call.
        first_cut_s: `generate.enter_first_cut`, summed over every call; contains ``decide_s``.
    """

    generation_s: float
    decide_s: float
    first_cut_s: float


def require_call_composition(key: str, calls: Mapping[str, int], ramp_out_members: int) -> None:
    """Refuse call counts that do not match one first cut with *ramp_out_members* emitted members.

    Raises:
        CallCompositionError: A wrapped callable is missing or extra, a build or
            first-cut count is not one, or the ``next_member`` count is not one
            per member plus the root's.
    """
    if set(calls) != set(WRAPPED_CALLABLES):
        raise CallCompositionError(f"{key}: calls name exactly {list(WRAPPED_CALLABLES)}, not {sorted(calls)}.")
    for name in (TANGENT_PENCIL_BUILD, RAMP_OUT_BUILD, ENTER_FIRST_CUT):
        if calls[name] != FIRST_CUTS_PER_GENERATION:
            raise CallCompositionError(f"{key}: {name} ran {calls[name]} time(s), not once per generation.")
    expected = ramp_out_members + ROOT_NEXT_MEMBER_CALLS
    if calls[RAMP_OUT_NEXT_MEMBER] != expected:
        raise CallCompositionError(f"{key}: {RAMP_OUT_NEXT_MEMBER} ran {calls[RAMP_OUT_NEXT_MEMBER]} time(s), not {ramp_out_members} member(s) + the root = {expected}.")


@dataclass(frozen=True)
class CasePencilCost:
    """One ``pipeline`` gate case's pencil cost over `RUNS_PER_CASE` fresh-process runs.

    Attributes:
        runs: The runs, in run order.
        calls: Wrapped callable (`WRAPPED_CALLABLES`) to its call count per run (equal in every run).
        ramp_out_members: Pencil members the first cut emitted (`RampOutRecord.members`).
        operations: The operation count every run produced.
        digest: The ops digest every run produced, equal to the timing document's.
        load_before: Load average (1, 5, 15 minutes) before the case's first run.
        load_after: Load average after its last run.
    """

    runs: tuple[PencilRun, ...]
    calls: Mapping[str, int]
    ramp_out_members: int
    operations: int
    digest: str
    load_before: LoadAverage
    load_after: LoadAverage

    @property
    def median_generation_s(self) -> float:
        """The low median of the runs' ``generation_s``, one of the measured values."""
        return float(statistics.median_low(run.generation_s for run in self.runs))

    @property
    def median_decide_s(self) -> float:
        """The low median of the runs' ``decide_s`` (`C_pencil`), one of the measured values."""
        return float(statistics.median_low(run.decide_s for run in self.runs))

    @property
    def median_first_cut_s(self) -> float:
        """The low median of the runs' ``first_cut_s``, one of the measured values."""
        return float(statistics.median_low(run.first_cut_s for run in self.runs))


def is_pencil_key(key: str) -> bool:
    """Whether *key* is a case key of the ``pipeline`` generator."""
    return is_case_key(key) and key.split(CASE_KEY_SEPARATOR)[1] == PENCIL_GENERATOR


def pencil_cases() -> dict[str, tuple[float, str]]:
    """Every ``pipeline`` gate case key to its (cap degrees, pocket name), in gate-case order."""
    return {case_key(cap_id, generator, pocket): (cap, pocket) for cap_id, cap, generator, pocket in gate_cases() if generator == PENCIL_GENERATOR}


def pencil_case_keys() -> tuple[str, ...]:
    """Every ``pipeline`` gate case key, in gate-case order."""
    return tuple(pencil_cases())


@dataclass(frozen=True)
class PipelinePencilCosts:
    """Every ``pipeline`` gate case's pencil cost, and the ones the admissions exclude.

    Raw construction is unchecked; `build` validates.

    Attributes:
        timing_sha256: SHA-256 of the timing document whose digests the runs matched.
        indicative: True when the run was measured without the quiet-machine check.
        cases: Case key to its cost, in gate-case order.
        excluded: ``pipeline`` case keys the admissions exclude, not generated.
        provenance: Where the run came from.
    """

    timing_sha256: str
    indicative: bool
    cases: Mapping[str, CasePencilCost]
    excluded: tuple[str, ...]
    provenance: RunProvenance

    @classmethod
    def build(cls, raw: object) -> PipelinePencilCosts:
        """Validate a parsed JSON document into pencil costs.

        Args:
            raw: The output of `json.loads`.

        Returns:
            The validated costs.

        Raises:
            MalformedPencilCostError: The document breaks the schema: wrong keys
                or schema tag, a run count other than `RUNS_PER_CASE`, a
                negative or non-numeric time or load, a decision cost above its
                first cut, a missing call count, a non-positive count, a
                malformed digest or SHA-256, or a ``pipeline`` gate case neither
                measured nor excluded.
            CallCompositionError: A case's call counts do not match its first cut.
        """
        if not isinstance(raw, dict) or set(raw) != _DOCUMENT_KEYS:
            raise MalformedPencilCostError(f"a pencil-cost document has exactly the keys {sorted(_DOCUMENT_KEYS)}.")
        if raw["schema"] != SCHEMA:
            raise MalformedPencilCostError(f"schema {raw['schema']!r} is not {SCHEMA!r}.")
        if raw["runs_per_case"] != RUNS_PER_CASE:
            raise MalformedPencilCostError(f"runs_per_case {raw['runs_per_case']!r} is not {RUNS_PER_CASE}.")
        timing_sha256 = raw["timing_sha256"]
        if not isinstance(timing_sha256, str) or _SHA256_PATTERN.fullmatch(timing_sha256) is None:
            raise MalformedPencilCostError(f"timing_sha256 is {SHA256_HEX_CHARS} lowercase hex characters.")
        indicative = raw["indicative"]
        if not isinstance(indicative, bool):
            raise MalformedPencilCostError("indicative is a boolean.")
        excluded = raw["excluded"]
        if not isinstance(excluded, list) or not all(is_case_key(key) for key in excluded):
            raise MalformedPencilCostError("excluded is a list of case keys.")
        cases_raw = raw["cases"]
        if not isinstance(cases_raw, dict):
            raise MalformedPencilCostError("cases is an object keyed by case key.")
        cases = {key: _case_cost(key, value) for key, value in cases_raw.items()}
        expected = pencil_case_keys()
        measured_in_order = [key for key in expected if key in cases]
        excluded_in_order = [key for key in expected if key in excluded]
        if list(cases) != measured_in_order or list(excluded) != excluded_in_order or sorted(measured_in_order + excluded_in_order) != sorted(expected):
            raise MalformedPencilCostError(f"every {PENCIL_GENERATOR} gate case is either measured or excluded, exactly once, in gate-case order.")
        try:
            provenance = RunProvenance.build(raw["provenance"])
        except MalformedProvenanceError as error:
            raise MalformedPencilCostError(f"provenance: {error}") from error
        return cls(timing_sha256=timing_sha256, indicative=indicative, cases=cases, excluded=tuple(excluded), provenance=provenance)

    @property
    def max_one_minute_load(self) -> float:
        """The highest 1-minute load average recorded before or after any case."""
        return max(load[0] for case in self.cases.values() for load in (case.load_before, case.load_after))

    @property
    def max_fifteen_minute_load(self) -> float:
        """The highest 15-minute load average recorded before or after any case."""
        return max(load[2] for case in self.cases.values() for load in (case.load_before, case.load_after))

    def to_json(self) -> dict[str, object]:
        """The JSON-ready document."""
        return {
            "schema": SCHEMA,
            "runs_per_case": RUNS_PER_CASE,
            "timing_sha256": self.timing_sha256,
            "indicative": self.indicative,
            "cases": {
                key: {
                    "runs": [{"generation_s": run.generation_s, "decide_s": run.decide_s, "first_cut_s": run.first_cut_s} for run in case.runs],
                    "calls": {name: case.calls[name] for name in WRAPPED_CALLABLES},
                    "ramp_out_members": case.ramp_out_members,
                    "operations": case.operations,
                    "digest": case.digest,
                    "load_before": list(case.load_before),
                    "load_after": list(case.load_after),
                }
                for key, case in self.cases.items()
            },
            "excluded": list(self.excluded),
            "provenance": self.provenance.to_json(),
        }

    def dumps(self) -> str:
        """The document as the committed text, ending in a newline."""
        return json.dumps(self.to_json(), indent=2) + "\n"


def _seconds_value(key: str, name: str, value: object) -> float:
    if not isinstance(value, (int, float)) or isinstance(value, bool) or value < 0:
        raise MalformedPencilCostError(f"{key}: {name} is a non-negative number.")
    return float(value)


def _count_value(key: str, name: str, value: object) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        raise MalformedPencilCostError(f"{key}: {name} is a positive integer.")
    return value


def _load_value(key: str, name: str, value: object) -> LoadAverage:
    load = parse_load_average(value)
    if load is None:
        raise MalformedPencilCostError(f"{key}: {name} is a list of {LOAD_AVERAGE_WINDOWS} non-negative load averages (1, 5, 15 minutes).")
    return load


def _run(key: str, value: object) -> PencilRun:
    if not isinstance(value, dict) or set(value) != _RUN_KEYS:
        raise MalformedPencilCostError(f"{key}: a run has exactly the keys {sorted(_RUN_KEYS)}.")
    run = PencilRun(**{name: _seconds_value(key, name, value[name]) for name in sorted(_RUN_KEYS)})
    if not run.decide_s <= run.first_cut_s <= run.generation_s:
        raise MalformedPencilCostError(f"{key}: a run nests decide_s <= first_cut_s <= generation_s, not {run!r}.")
    return run


def _calls(key: str, value: object) -> dict[str, int]:
    if not isinstance(value, dict) or set(value) != set(WRAPPED_CALLABLES):
        raise MalformedPencilCostError(f"{key}: calls has exactly the keys {list(WRAPPED_CALLABLES)}.")
    return {name: _count_value(key, f"calls[{name}]", value[name]) for name in WRAPPED_CALLABLES}


def _case_cost(key: str, value: object) -> CasePencilCost:
    if not is_pencil_key(key):
        raise MalformedPencilCostError(f"{key!r} is not a {PENCIL_GENERATOR} case key.")
    if not isinstance(value, dict) or set(value) != _CASE_KEYS:
        raise MalformedPencilCostError(f"{key}: a case has exactly the keys {sorted(_CASE_KEYS)}.")
    runs = value["runs"]
    if not isinstance(runs, list) or len(runs) != RUNS_PER_CASE:
        raise MalformedPencilCostError(f"{key}: runs holds exactly {RUNS_PER_CASE} runs.")
    digest = value["digest"]
    if not isinstance(digest, str) or _DIGEST_PATTERN.fullmatch(digest) is None:
        raise MalformedPencilCostError(f"{key}: digest is {DIGEST_HEX_CHARS} lowercase hex characters.")
    members = value["ramp_out_members"]
    if not isinstance(members, int) or isinstance(members, bool) or members < 0:
        raise MalformedPencilCostError(f"{key}: ramp_out_members is a non-negative integer.")
    calls = _calls(key, value["calls"])
    require_call_composition(key, calls, members)
    return CasePencilCost(
        runs=tuple(_run(key, run) for run in runs),
        calls=calls,
        ramp_out_members=members,
        operations=_count_value(key, "operations", value["operations"]),
        digest=digest,
        load_before=_load_value(key, "load_before", value["load_before"]),
        load_after=_load_value(key, "load_after", value["load_after"]),
    )


def load_pipeline_pencil_costs(path: Path) -> PipelinePencilCosts:
    """Read and validate the pencil-cost document at *path*.

    Raises:
        MalformedPencilCostError: The document breaks the schema.
        CallCompositionError: A case's call counts do not match its first cut.
    """
    return PipelinePencilCosts.build(json.loads(path.read_text()))


def _seconds(value: float) -> str:
    return f"{value:.{RENDERED_SECONDS_DECIMALS}f}"


def _quietness(costs: PipelinePencilCosts) -> str:
    peak = f"highest recorded 1-minute load {costs.max_one_minute_load:.{LOAD_DECIMALS}f}"
    if costs.indicative:
        return f"**INDICATIVE**: measured without the quiet-machine check ({peak}); these seconds yield no time-gate verdict."
    return f"Quiet-machine run: every case started and ended at a 1-minute load of at most {QUIET_LOAD_AVERAGE_MAX} ({peak})."


def render_pipeline_pencil_cost_table(costs: PipelinePencilCosts) -> str:
    """The Markdown method statement, the measured load and the per-case table of *costs*.

    Args:
        costs: The validated costs.

    Returns:
        The Markdown text, ending in a newline.
    """
    lines = [
        "<!-- Generated by tools.pipeline_pencil_cost from benchmarks/results/pipeline_pencil_cost.json. Do not edit; run `pixi run pipeline-pencil-cost`. -->",
        "",
        (
            f"Pencil decision cost `C_pencil` (Task 1.2, ruling D6), measured by `tools.pipeline_pencil_cost`: every `{PENCIL_GENERATOR}` "
            f"gate case, one case at a time, each of {RUNS_PER_CASE} generations in a fresh process; the table reports the low median "
            f"of each clock. `decide_s` (`C_pencil`) sums the wall clock of every `{TANGENT_PENCIL_BUILD}`, `{RAMP_OUT_BUILD}` and "
            f"`{RAMP_OUT_NEXT_MEMBER}` call, each counted on its own: one build of each and one `{ENTER_FIRST_CUT}` per generation, "
            f"and one `{RAMP_OUT_NEXT_MEMBER}` per emitted member plus the root's, or the task fails. `first_cut_s` sums "
            f"`{ENTER_FIRST_CUT}` (helix, ramp-out emission and root loop) and contains `decide_s`; `generation_s` is the whole "
            "instrumented generator call, not the uninstrumented timing median. The wrappers only time: the ops digest equals the "
            "uninstrumented digest of the same case in `benchmarks/results/pipeline_timing.json` (SHA-256 "
            f"`{costs.timing_sha256}`), or the task fails. Seconds are wall clock on the machine that ran the task and are valid for "
            "the Task 1.2 time gate only beside a `pipeline-timing` run of the same session. Load averages (1 / 5 / 15 minutes) "
            "are read before each case's first run and after its last."
        ),
        "",
        _quietness(costs),
        "",
        "| cap | pocket | `decide_s` = `C_pencil` (median) | `first_cut_s` (median) | instrumented `generation_s` (median) | `decide_s` runs (s) "
        f"| calls ({TANGENT_PENCIL_BUILD} / {RAMP_OUT_BUILD} / {RAMP_OUT_NEXT_MEMBER} / {ENTER_FIRST_CUT}) | members | operations | ops digest "
        "| load before | load after |",
        "| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :--- | ---: | ---: |",
    ]
    for key, case in costs.cases.items():
        cap_id, _generator, pocket = key.split(CASE_KEY_SEPARATOR)
        runs = ", ".join(_seconds(run.decide_s) for run in case.runs)
        calls = " / ".join(str(case.calls[name]) for name in WRAPPED_CALLABLES)
        lines.append(
            f"| {cap_id} | {pocket} | {_seconds(case.median_decide_s)} | {_seconds(case.median_first_cut_s)} | {_seconds(case.median_generation_s)} "
            f"| {runs} | {calls} | {case.ramp_out_members} | {case.operations} | `{case.digest}` | {format_load(case.load_before)} | {format_load(case.load_after)} |"
        )
    if costs.excluded:
        lines += ["", "Excluded by the gate admissions, not generated: " + ", ".join(f"`{key}`" for key in costs.excluded) + "."]
    return "\n".join(lines) + "\n"


def rooted(path: Path) -> Path:
    """*path*, resolved against the repository root when relative."""
    return path if path.is_absolute() else ROOT / path
