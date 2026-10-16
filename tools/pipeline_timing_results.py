"""The serial generation-time results file: its schema, validation and serialisation.

`benchmarks/results/pipeline_timing.json` holds, per gate case, the wall-clock
seconds of each of `RUNS_PER_CASE` generations (each in a fresh process), the
operation count and the ops digest the runs agreed on. `tools.pipeline_timing`
measures it and `tools.pipeline_timing_table` renders it for the docs.
"""

from __future__ import annotations

import json
import re
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from tools.gate_case_keys import is_case_key
from tools.gate_measure import case_key
from tools.gate_measure import gate_cases

SCHEMA = "pipeline-timing/1"
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


class MalformedPipelineTimingsError(Exception):
    """The timing document breaks its schema."""


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
    """

    cases: Mapping[str, CaseTiming]
    excluded: tuple[str, ...]

    @classmethod
    def build(cls, raw: object) -> PipelineTimings:
        """Validate a parsed JSON document into timings.

        Args:
            raw: The output of `json.loads`.

        Returns:
            The validated timings.

        Raises:
            MalformedPipelineTimingsError: The document breaks the schema: wrong
                schema tag, a non-gate case key, a run count other than
                `RUNS_PER_CASE`, a negative or non-numeric time, a malformed
                digest, or a gate case neither timed nor excluded.
        """
        if not isinstance(raw, dict) or set(raw) != {"schema", "runs_per_case", "cases", "excluded"}:
            raise MalformedPipelineTimingsError("a timing document has exactly the keys schema, runs_per_case, cases and excluded.")
        if raw["schema"] != SCHEMA:
            raise MalformedPipelineTimingsError(f"schema {raw['schema']!r} is not {SCHEMA!r}.")
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
        return cls(cases=cases, excluded=tuple(excluded))

    def to_json(self) -> dict[str, object]:
        """The JSON-ready document."""
        return {
            "schema": SCHEMA,
            "runs_per_case": RUNS_PER_CASE,
            "cases": {key: {"runs_s": list(case.runs_s), "operations": case.operations, "digest": case.digest} for key, case in self.cases.items()},
            "excluded": list(self.excluded),
        }

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
