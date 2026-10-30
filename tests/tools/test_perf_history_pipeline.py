"""The pipeline-record adapters (``tools.perf_history_adapters_pipeline``), registered in ``tools.perf_history_adapters``.

They read ``pipeline_timing`` /1 and /2 and ``pipeline_pencil_cost`` through
their tools' validating loaders. The core's own /1 reader
(``tools.perf_history_adapters.pipeline_timing``) wrote the history before
registration; on /1 records the two must agree, so no history line is
duplicated across the switch.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

from tests.tools.test_pipeline_pencil_cost import _synthetic as synthetic_pencil_costs
from tests.tools.test_pipeline_timing import _with_load as synthetic_timings_with_load
from tools.perf_history_adapters import FAMILIES
from tools.perf_history_adapters import PENDING_FAMILIES
from tools.perf_history_adapters import pipeline_timing as core_pipeline_timing
from tools.perf_history_adapter_types import AdaptedRecord
from tools.perf_history_adapters_pipeline import PIPELINE_FAMILIES
from tools.perf_history_entry import PerfEntry
from tools.perf_history_entry import RecordSource
from tools.perf_history_ingest import record_entries

ROOT = Path(__file__).resolve().parents[2]
TIMING = "benchmarks/results/pipeline_timing.json"
PENCIL = "benchmarks/results/pipeline_pencil_cost.json"
# A commit whose pipeline_timing.json is schema /1 (the last one before the PENCIL unit wrote /2).
SCHEMA_1_REVISION = "1fcda9c4"


def _no_sibling(name: str) -> bytes:
    raise AssertionError(f"no sibling expected, asked for {name}")


def _entries(record_path: str, adapted: AdaptedRecord) -> tuple[PerfEntry, ...]:
    source = RecordSource(family=Path(record_path).stem, record_path=record_path, record_sha256="a" * 64)
    return tuple(PerfEntry.build(source, adapted.provenance, measurement) for measurement in adapted.measurements)


def test_every_pipeline_family_is_registered_with_its_module_adapter() -> None:
    assert set(PIPELINE_FAMILIES) == {"pipeline_timing.json", "pipeline_pencil_cost.json"}
    assert all(FAMILIES[name] is adapter and name not in PENDING_FAMILIES for name, adapter in PIPELINE_FAMILIES.items())


def test_a_quiet_pipeline_timing_2_record_yields_quiet_entries_with_its_provenance() -> None:
    timings = synthetic_timings_with_load()
    entries = _entries(TIMING, PIPELINE_FAMILIES["pipeline_timing.json"](json.loads(timings.dumps()), _no_sibling))
    assert len(entries) == len(timings.cases)
    assert timings.provenance is not None and all(entry.quiet and entry.provenance.commit == timings.provenance.rev for entry in entries)
    assert all(entry.provenance.native_sha256 == timings.provenance.native_sha256 for entry in entries)
    assert all(PerfEntry.from_json_line(entry.to_json_line()) == entry for entry in entries)


def test_a_pencil_cost_record_round_trips_and_its_indicative_flag_is_kept() -> None:
    adapter = PIPELINE_FAMILIES["pipeline_pencil_cost.json"]
    quiet = _entries(PENCIL, adapter(json.loads(synthetic_pencil_costs().dumps()), _no_sibling))
    indicative = _entries(PENCIL, adapter(json.loads(synthetic_pencil_costs(indicative=True).dumps()), _no_sibling))
    assert quiet and all(entry.quiet for entry in quiet)
    assert indicative and not any(entry.quiet for entry in indicative)
    assert all(PerfEntry.from_json_line(entry.to_json_line()) == entry for entry in (*quiet, *indicative))


def test_on_a_committed_schema_1_record_the_module_agrees_with_the_core() -> None:
    """Registering the module must not change a /1 entry's digest: the history stays deduplicated across the switch."""
    committed = subprocess.run(["git", "show", f"{SCHEMA_1_REVISION}:{TIMING}"], cwd=ROOT, capture_output=True, check=True).stdout
    assert json.loads(committed)["schema"] == "pipeline-timing/1"
    registered = record_entries(TIMING, committed, _no_sibling)
    source = registered[0].source
    core = core_pipeline_timing(json.loads(committed), _no_sibling)
    assert tuple(PerfEntry.build(source, core.provenance, measurement) for measurement in core.measurements) == registered
