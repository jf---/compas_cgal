"""Adapters for the product-pipeline timing records: ``pipeline_timing`` /1 and /2, and ``pipeline_pencil_cost``.

These depend on ``tools.pipeline_pencil_cost_results``, ``tools.run_provenance``
and the ``pipeline_timing`` /2 fields of ``tools.pipeline_timing_results``.
``tools.perf_history_adapters`` registers them: ``PIPELINE_FAMILIES`` merges
into its ``FAMILIES``. This module imports only the shared adapter types
(``tools.perf_history_adapter_types``), never the registry.

Both records are read through their tools' validating loaders, so a record the
tool would refuse is refused here too.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any
from typing import Optional

from tools.perf_history_adapter_types import AdaptedRecord
from tools.perf_history_adapter_types import Adapter
from tools.perf_history_adapter_types import SiblingReader
from tools.perf_history_entry import ClockSummary
from tools.perf_history_entry import LoadAverage1m
from tools.perf_history_entry import Measurement
from tools.perf_history_entry import RecordProvenance
from tools.perf_history_json import MalformedTimingRecordError
from tools.pipeline_pencil_cost_results import CallCompositionError
from tools.pipeline_pencil_cost_results import MalformedPencilCostError
from tools.pipeline_pencil_cost_results import PipelinePencilCosts
from tools.pipeline_timing_results import MalformedPipelineTimingsError
from tools.pipeline_timing_results import PipelineTimings
from tools.run_provenance import RunProvenance


def _run_provenance(provenance: Optional[RunProvenance]) -> RecordProvenance:
    if provenance is None:
        return RecordProvenance()
    return RecordProvenance(
        commit=provenance.rev, dirty=provenance.dirty, native_sha256=provenance.native_sha256, machine=provenance.cpu_model, recorded_at=provenance.started_utc
    )


def pipeline_timing(raw: Any, _sibling: SiblingReader) -> AdaptedRecord:
    """Per gate case: the fresh-process generation runs, read through ``PipelineTimings.build`` (/1 or /2).

    Load is the record's 1-minute load before the first case and after the last;
    a record its own ``quiet`` rule refuses (15-minute load at the end) is indicative.
    """
    try:
        timings = PipelineTimings.build(raw)
    except MalformedPipelineTimingsError as error:
        raise MalformedTimingRecordError(f"pipeline_timing: {error}") from error
    load = (
        None if timings.load_start is None else LoadAverage1m(timings.load_start[0]),
        None if timings.load_end is None else LoadAverage1m(timings.load_end[0]),
    )
    downgrade = None if timings.quiet or timings.load_start is None else "the record's own quiet rule (1-minute at both ends and 15-minute at the end) refuses it"
    measurements = tuple(
        Measurement(key, "pipeline_timing/generation", None, ClockSummary.over_runs(case.runs_s), load[0], load[1], downgrade) for key, case in timings.cases.items()
    )
    return AdaptedRecord(_run_provenance(timings.provenance), measurements)


def pipeline_pencil_cost(raw: Any, _sibling: SiblingReader) -> AdaptedRecord:
    """Per ``pipeline`` gate case: generation, pencil decisions and first cut runs, read through ``PipelinePencilCosts.build``."""
    try:
        costs = PipelinePencilCosts.build(raw)
    except (MalformedPencilCostError, CallCompositionError) as error:
        raise MalformedTimingRecordError(f"pipeline_pencil_cost: {error}") from error
    downgrade = "the record marks the run indicative (measured without the quiet-machine check)" if costs.indicative else None
    measurements: list[Measurement] = []
    for key, case in costs.cases.items():
        before, after = LoadAverage1m(case.load_before[0]), LoadAverage1m(case.load_after[0])
        clocks = {
            "generation": [run.generation_s for run in case.runs],
            "decide": [run.decide_s for run in case.runs],
            "first_cut": [run.first_cut_s for run in case.runs],
        }
        for stage, seconds in clocks.items():
            measurements.append(Measurement(key, f"pipeline_pencil_cost/{stage}", None, ClockSummary.over_runs(seconds), before, after, downgrade))
    return AdaptedRecord(_run_provenance(costs.provenance), tuple(measurements))


# Record file name -> adapter; merged into ``tools.perf_history_adapters.FAMILIES``.
PIPELINE_FAMILIES: Mapping[str, Adapter] = {
    "pipeline_pencil_cost.json": pipeline_pencil_cost,
    "pipeline_timing.json": pipeline_timing,
}
