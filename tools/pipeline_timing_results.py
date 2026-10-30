"""The serial generation-time results file: its schema, validation and serialisation.

`benchmarks/results/pipeline_timing.json` holds, per gate case, the wall-clock
seconds of each of `RUNS_PER_CASE` generations (each in a fresh process), the
operation count and the ops digest the runs agreed on. Schema
``pipeline-timing/2`` also holds the 1/5/15-minute load average at the start
and at the end of the run and its `tools.run_provenance.RunProvenance` (tree,
tool, environment, machine, UTC interval); a ``pipeline-timing/1`` record (written before load
was recorded, e.g. the Task 1.1 record at 7b001d47) still loads, with its load
and provenance unrecorded (None), and serialises back to ``/1`` unchanged.
`tools.pipeline_timing` measures it and `tools.pipeline_timing_table` renders
it for the docs.
"""

from __future__ import annotations

import json
import re
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping
from typing import Optional

from tools.gate_case_keys import is_case_key
from tools.gate_measure import case_key
from tools.gate_measure import gate_cases
from tools.held_generation_clock import QUIET_LOAD_AVERAGE_MAX
from tools.held_generation_clock import LoadAverage
from tools.run_provenance import MalformedProvenanceError
from tools.run_provenance import RunProvenance

SCHEMA = "pipeline-timing/2"
# The schema before load was recorded; read and round-tripped, never written by the tool.
SCHEMA_WITHOUT_LOAD = "pipeline-timing/1"
PIPELINE_TIMING_PATH = Path("benchmarks") / "results" / "pipeline_timing.json"
PIPELINE_TIMING_TABLE_PATH = Path("benchmarks") / "results" / "pipeline_timing.md"
# Generations per case; odd, so the median is one of the measured runs.
RUNS_PER_CASE = 3
# Stored precision of a run's seconds: 0.1 ms, below the run-to-run noise of a
# fresh-process generation.
STORED_SECONDS_DECIMALS = 4
# Hex characters of the truncated SHA-256 ops digest (64 bits).
DIGEST_HEX_CHARS = 16
_DIGEST_PATTERN = re.compile(rf"[0-9a-f]{{{DIGEST_HEX_CHARS}}}")
# Stored and rendered precision of a load average (runnable processes), as `uptime` prints it.
LOAD_DECIMALS = 2
# Windows of a load average: 1, 5 and 15 minutes (`os.getloadavg`).
LOAD_AVERAGE_WINDOWS = 3
_KEYS_WITHOUT_LOAD = frozenset({"schema", "runs_per_case", "cases", "excluded"})
_KEYS_WITH_LOAD = _KEYS_WITHOUT_LOAD | {"load_start", "load_end", "provenance"}


class MalformedPipelineTimingsError(Exception):
    """The timing document breaks its schema."""


def parse_load_average(value: object) -> Optional[LoadAverage]:
    """*value* as a load average if it is a list of `LOAD_AVERAGE_WINDOWS` non-negative numbers, else None (the caller raises)."""
    if not isinstance(value, list) or len(value) != LOAD_AVERAGE_WINDOWS:
        return None
    if not all(isinstance(window, (int, float)) and not isinstance(window, bool) and window >= 0 for window in value):
        return None
    one, five, fifteen = (float(window) for window in value)
    return (one, five, fifteen)


def stored_load(load: LoadAverage) -> LoadAverage:
    """*load* rounded to `LOAD_DECIMALS`, as records store it."""
    one, five, fifteen = (round(window, LOAD_DECIMALS) for window in load)
    return (one, five, fifteen)


def format_load(load: LoadAverage) -> str:
    """*load* as ``1 / 5 / 15``-minute averages at `LOAD_DECIMALS`."""
    return " / ".join(f"{window:.{LOAD_DECIMALS}f}" for window in load)


@dataclass(frozen=True)
class CaseTiming:
    """One gate case's serial timing.

    Attributes:
        runs_s: Wall-clock seconds of each generation, in run order.
        operations: The number of operations every run produced.
        digest: The ops digest every run produced.
    """

    runs_s: tuple[float, ...]
    operations: int
    digest: str

    @property
    def median_s(self) -> float:
        """The median of `runs_s`, one of the measured runs."""
        return float(statistics.median_low(self.runs_s))


@dataclass(frozen=True)
class PipelineTimings:
    """Every gate case's serial timing, and the cases the gate admissions exclude.

    Raw construction is unchecked; `build` validates.

    Attributes:
        cases: Case key to its timing, in gate-case order.
        excluded: Case keys the admissions exclude, which are not generated.
        load_start: Load average (1, 5, 15 minutes) before the first case, or
            None for a ``pipeline-timing/1`` record.
        load_end: Load average after the last case, or None likewise.
        provenance: Where the run came from, or None likewise.
    """

    cases: Mapping[str, CaseTiming]
    excluded: tuple[str, ...]
    load_start: Optional[LoadAverage] = None
    load_end: Optional[LoadAverage] = None
    provenance: Optional[RunProvenance] = None

    @classmethod
    def build(cls, raw: object) -> PipelineTimings:
        """Validate a parsed JSON document into timings.

        Args:
            raw: The output of `json.loads`.

        Returns:
            The validated timings.

        Raises:
            MalformedPipelineTimingsError: The document breaks the schema: wrong
                schema tag or keys, a non-gate case key, a run count other than
                `RUNS_PER_CASE`, a negative or non-numeric time, a malformed
                digest or load, or a gate case neither timed nor excluded.
        """
        if not isinstance(raw, dict) or raw.get("schema") not in (SCHEMA, SCHEMA_WITHOUT_LOAD):
            raise MalformedPipelineTimingsError(f"a timing document's schema is {SCHEMA!r} or {SCHEMA_WITHOUT_LOAD!r}.")
        with_load = raw["schema"] == SCHEMA
        expected_keys = _KEYS_WITH_LOAD if with_load else _KEYS_WITHOUT_LOAD
        if set(raw) != expected_keys:
            raise MalformedPipelineTimingsError(f"a {raw['schema']} document has exactly the keys {sorted(expected_keys)}.")
        load_start = parse_load_average(raw["load_start"]) if with_load else None
        load_end = parse_load_average(raw["load_end"]) if with_load else None
        if with_load and (load_start is None or load_end is None):
            raise MalformedPipelineTimingsError(f"load_start and load_end are lists of {LOAD_AVERAGE_WINDOWS} non-negative load averages (1, 5, 15 minutes).")
        provenance = None
        if with_load:
            try:
                provenance = RunProvenance.build(raw["provenance"])
            except MalformedProvenanceError as error:
                raise MalformedPipelineTimingsError(f"provenance: {error}") from error
        if raw["runs_per_case"] != RUNS_PER_CASE:
            raise MalformedPipelineTimingsError(f"runs_per_case {raw['runs_per_case']!r} is not {RUNS_PER_CASE}.")
        excluded = raw["excluded"]
        if not isinstance(excluded, list) or not all(is_case_key(key) for key in excluded):
            raise MalformedPipelineTimingsError("excluded is a list of case keys.")
        cases_raw = raw["cases"]
        if not isinstance(cases_raw, dict):
            raise MalformedPipelineTimingsError("cases is an object keyed by case key.")
        cases = {key: _case_timing(key, value) for key, value in cases_raw.items()}
        expected = [case_key(cap_id, generator, pocket) for cap_id, _cap, generator, pocket in gate_cases()]
        timed_in_order = [key for key in expected if key in cases]
        excluded_in_order = [key for key in expected if key in excluded]
        if list(cases) != timed_in_order or list(excluded) != excluded_in_order or sorted(timed_in_order + excluded_in_order) != sorted(expected):
            raise MalformedPipelineTimingsError("every gate case is either timed or excluded, exactly once, in gate-case order.")
        return cls(cases=cases, excluded=tuple(excluded), load_start=load_start, load_end=load_end, provenance=provenance)

    @property
    def max_one_minute_load(self) -> Optional[float]:
        """The higher 1-minute load average of the start and the end, or None when unrecorded."""
        if self.load_start is None or self.load_end is None:
            return None
        return max(self.load_start[0], self.load_end[0])

    @property
    def end_fifteen_minute_load(self) -> Optional[float]:
        """The 15-minute load average at the end, spanning a run of up to 15 minutes, or None when unrecorded."""
        return None if self.load_end is None else self.load_end[2]

    @property
    def quiet(self) -> bool:
        """Whether the load was recorded and at most `QUIET_LOAD_AVERAGE_MAX`: 1-minute at the start and the end, and 15-minute at the end.

        The load is read only before the first case and after the last; the
        15-minute window at the end covers the middle of the run as well.
        """
        peak, end_fifteen = self.max_one_minute_load, self.end_fifteen_minute_load
        return peak is not None and end_fifteen is not None and peak <= QUIET_LOAD_AVERAGE_MAX and end_fifteen <= QUIET_LOAD_AVERAGE_MAX

    def to_json(self) -> dict[str, object]:
        """The JSON-ready document: ``pipeline-timing/2`` when load and provenance are recorded, ``/1`` when neither is.

        Raises:
            MalformedPipelineTimingsError: Only some of the ``/2`` fields are set.
        """
        recorded = (self.load_start is not None, self.load_end is not None, self.provenance is not None)
        if any(recorded) and not all(recorded):
            raise MalformedPipelineTimingsError("a pipeline-timing/2 record carries load_start, load_end and provenance together.")
        document: dict[str, object] = {
            "schema": SCHEMA if all(recorded) else SCHEMA_WITHOUT_LOAD,
            "runs_per_case": RUNS_PER_CASE,
            "cases": {key: {"runs_s": list(case.runs_s), "operations": case.operations, "digest": case.digest} for key, case in self.cases.items()},
            "excluded": list(self.excluded),
        }
        if self.load_start is not None and self.load_end is not None and self.provenance is not None:
            document["load_start"] = list(self.load_start)
            document["load_end"] = list(self.load_end)
            document["provenance"] = self.provenance.to_json()
        return document

    def dumps(self) -> str:
        """The document as the committed text, ending in a newline."""
        return json.dumps(self.to_json(), indent=2) + "\n"


def _case_timing(key: str, value: object) -> CaseTiming:
    if not is_case_key(key):
        raise MalformedPipelineTimingsError(f"{key!r} is not a case key.")
    if not isinstance(value, dict) or set(value) != {"runs_s", "operations", "digest"}:
        raise MalformedPipelineTimingsError(f"{key}: a case has exactly the keys runs_s, operations and digest.")
    runs = value["runs_s"]
    if not isinstance(runs, list) or len(runs) != RUNS_PER_CASE:
        raise MalformedPipelineTimingsError(f"{key}: runs_s holds exactly {RUNS_PER_CASE} times.")
    if not all(isinstance(run, (int, float)) and not isinstance(run, bool) and run >= 0 for run in runs):
        raise MalformedPipelineTimingsError(f"{key}: every run is a non-negative number of seconds.")
    operations = value["operations"]
    if not isinstance(operations, int) or isinstance(operations, bool) or operations < 0:
        raise MalformedPipelineTimingsError(f"{key}: operations is a non-negative integer.")
    digest = value["digest"]
    if not isinstance(digest, str) or _DIGEST_PATTERN.fullmatch(digest) is None:
        raise MalformedPipelineTimingsError(f"{key}: digest is {DIGEST_HEX_CHARS} lowercase hex characters.")
    return CaseTiming(runs_s=tuple(float(run) for run in runs), operations=operations, digest=digest)


def load_pipeline_timings(path: Path) -> PipelineTimings:
    """Read and validate the timing document at *path*.

    Raises:
        MalformedPipelineTimingsError: The document breaks the schema.
    """
    return PipelineTimings.build(json.loads(path.read_text()))
