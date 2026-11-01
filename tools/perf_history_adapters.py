"""One adapter per timing-record family under ``benchmarks/results``: a record's bytes to its measured stages.

A family is named by the record's file name (``FAMILIES``). Every JSON file
under ``benchmarks/results`` is either a timing family here or listed in
``NON_TIMING_RECORDS`` with the reason it holds no performance measurement;
anything else raises ``UnknownRecordFamilyError`` (never a silent skip), and
``tests/tools/test_perf_history.py`` pins that every committed record is
classified.

Each adapter reads only what its record states. Stage names are qualified by
the family (``<family>/<stage>``): two tools that both time "generation" count
different work (case load in or out, fresh process or not), so they never
share a baseline. A record that keeps no load average yields indicative
entries (``tools.perf_history_entry.indicative_reason`` says "no load
recorded"). Where a record keeps one load reading per run, the first run's is
the load before and the last run's the load after.

``pipeline_timing`` (/1 and /2) and ``pipeline_pencil_cost`` are read by the
adapters of ``tools.perf_history_adapters_pipeline``, registered here
(``PIPELINE_FAMILIES`` merged into ``FAMILIES``) through their tools'
validating loaders, so a record the tool would refuse is refused here too.
The core's own /1 reader (``pipeline_timing`` below) wrote the history before
that registration and stays as the reference the registered adapter must
match on /1 records. A family whose adapter exists but is not registered is
listed in ``PENDING_FAMILIES`` and raises ``PendingRecordFamilyError``: never
a silent skip. A record schema an adapter does not know raises
``UnknownRecordSchemaError``: a new schema needs its adapter extended, it is
never guessed at.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from collections.abc import Sequence
from pathlib import PurePosixPath
from typing import Any
from typing import Optional

from benchmarks.errors import BenchmarkError
from tools.coverage_certificate_budget_gate import ProductionBudget
from tools.held_generation_clock import MachineNotQuietError
from tools.perf_history_adapter_types import AdaptedRecord
from tools.perf_history_adapter_types import Adapter
from tools.perf_history_adapter_types import SiblingReader as SiblingReader
from tools.perf_history_adapters_pipeline import PIPELINE_FAMILIES
from tools.perf_history_entry import UNRECORDED_MACHINE
from tools.perf_history_entry import ClockStatistic
from tools.perf_history_entry import ClockSummary
from tools.perf_history_entry import LoadAverage1m
from tools.perf_history_entry import Measurement
from tools.perf_history_entry import RecordProvenance
from tools.perf_history_json import MalformedTimingRecordError
from tools.perf_history_json import field
from tools.perf_history_json import integer
from tools.perf_history_json import items
from tools.perf_history_json import number
from tools.perf_history_json import numbers
from tools.perf_history_json import obj
from tools.perf_history_json import optional_bool
from tools.perf_history_json import optional_text
from tools.perf_history_json import text
from tools.pipeline_timing_results import MalformedPipelineTimingsError
from tools.pipeline_timing_results import PipelineTimings

# ``AdaptedRecord``, ``Adapter`` and ``SiblingReader`` live in ``tools.perf_history_adapter_types`` so adapter
# modules registered here never import this registry; ``SiblingReader`` is re-exported (the ``as`` import)
# because ``tools.perf_history_ingest`` reads it from here.

# The provenance envelope written beside a ``benchmarks.cli corpus`` report.
CORPUS_STAMP = "stamp.json"


class UnknownRecordFamilyError(Exception):
    """A record's file name is neither a timing family nor a listed non-timing record."""


class PendingRecordFamilyError(Exception):
    """The record's family has an adapter that is not registered yet (``PENDING_FAMILIES``)."""


class NotATimingRecordError(Exception):
    """The record is a listed non-timing record; it has no performance measurement to ingest."""


class UnknownRecordSchemaError(Exception):
    """A record of a known family carries a schema its adapter does not read."""


class MissingSiblingRecordError(Exception):
    """A record's adapter needs a file beside it that is not there."""


def _machine(identity: object, where: str) -> str:
    if identity is None:
        return UNRECORDED_MACHINE
    brand = obj(identity, where).get("cpu_brand")
    return UNRECORDED_MACHINE if brand is None else text(brand, f"{where}.cpu_brand")


def _build_provenance(record: Mapping[str, Any], identity_key: str) -> RecordProvenance:
    """Provenance of a record carrying ``build`` (commit, dirty, recorded_at) and a machine identity object."""
    build = obj(field(record, "build", "record"), "build")
    return RecordProvenance(
        commit=optional_text(build, "commit", "build"),
        dirty=optional_bool(build, "dirty", "build"),
        source_sha256=optional_text(record, "source_sha256", "record"),
        native_sha256=optional_text(record, "native_sha256", "record"),
        machine=_machine(record.get(identity_key), identity_key),
        recorded_at=optional_text(build, "recorded_at", "build"),
    )


def _load(value: object, where: str) -> Optional[LoadAverage1m]:
    """A 1-minute load: a number, or the first of a (1, 5, 15)-minute triple; None when null."""
    if value is None:
        return None
    if isinstance(value, list):
        return LoadAverage1m(number(field_index(value, 0, where), f"{where}[0]"))
    return LoadAverage1m(number(value, where))


def field_index(value: Sequence[Any], index: int, where: str) -> Any:
    """Raises:
    MalformedTimingRecordError: ``value`` has no element ``index``.
    """
    if len(value) <= index:
        raise MalformedTimingRecordError(f"{where} has no element {index}.")
    return value[index]


def _load_pair(value: object, where: str) -> tuple[Optional[LoadAverage1m], Optional[LoadAverage1m]]:
    """A ``[before, after]`` pair of 1-minute loads."""
    pair = items(value, where)
    if len(pair) != 2:
        raise MalformedTimingRecordError(f"{where} is [before, after], not {pair!r}.")
    return _load(pair[0], f"{where}[0]"), _load(pair[1], f"{where}[1]")


def _first_last(loads: Sequence[Optional[LoadAverage1m]]) -> tuple[Optional[LoadAverage1m], Optional[LoadAverage1m]]:
    return (loads[0], loads[-1]) if loads else (None, None)


def _wall(workload: str, stage: str, wall: ClockSummary, load: tuple[Optional[LoadAverage1m], Optional[LoadAverage1m]] = (None, None)) -> Measurement:
    return Measurement(workload=workload, stage=stage, cpu=None, wall=wall, load_before=load[0], load_after=load[1])


# A stage's (process-CPU key or None, wall key) in a Held workload row.
StageFields = Mapping[str, tuple[Optional[str], str]]


def _cpu(workload: str, stage: str, cpu: ClockSummary, load: tuple[Optional[LoadAverage1m], Optional[LoadAverage1m]] = (None, None)) -> Measurement:
    return Measurement(workload=workload, stage=stage, cpu=cpu, wall=None, load_before=load[0], load_after=load[1])


def _stage_fields(row: Mapping[str, Any], workload: str, family: str, spec: StageFields, load: tuple[Optional[LoadAverage1m], Optional[LoadAverage1m]]) -> list[Measurement]:
    """The named seconds fields of one Held workload row, one single-run stage each; every named field must be present."""
    measurements: list[Measurement] = []
    for stage, (cpu_key, wall_key) in spec.items():
        cpu = None if cpu_key is None else ClockSummary.single(number(field(row, cpu_key, workload), f"{workload}.{cpu_key}"))
        wall = ClockSummary.single(number(field(row, wall_key, workload), f"{workload}.{wall_key}"))
        measurements.append(Measurement(workload, f"{family}/{stage}", cpu, wall, load[0], load[1]))
    return measurements


def _require_schema(record: Mapping[str, Any], accepted: Sequence[Optional[str]], family: str) -> Optional[str]:
    schema = record.get("schema")
    if schema not in accepted:
        raise UnknownRecordSchemaError(f"{family}: schema {schema!r}; the adapter reads {list(accepted)}.")
    return None if schema is None else str(schema)


# --- adapters ------------------------------------------------------------------------------------------


def actual_stock_profile(raw: Any, _sibling: SiblingReader) -> AdaptedRecord:
    """The profiled generation, and the generation before each budget run.

    Both ``generation_s`` are process CPU (``ReplayWorkload.build`` clocks them with
    ``time.process_time``); the record keeps no wall clock for them. Adapter
    version 1 filed them as wall (RV40 I1); version 2 supersedes those entries.

    The budget runs' sweep/query CPU is not a performance measurement: each run
    stops at its budget, so its seconds are the budget, not the work.
    """
    record = obj(raw, "record")
    _require_schema(record, ("actual-stock-profile/1",), "actual_stock_profile")
    family = "actual_stock_profile"
    workload = text(field(record, "workload", "record"), "workload")
    load = (_load(field(record, "load_start", "record"), "load_start"), _load(field(record, "load_end", "record"), "load_end"))
    measurements = [_cpu(workload, f"{family}/generation", ClockSummary.single(number(field(record, "generation_s", "record"), "generation_s")), load)]
    for index, value in enumerate(items(field(record, "budget_runs", "record"), "budget_runs")):
        run = obj(value, f"budget_runs[{index}]")
        where = f"budget_runs[{index}]"
        run_load = (_load(field(run, "load_before", where), f"{where}.load_before"), _load(field(run, "load_after", where), f"{where}.load_after"))
        measurements.append(
            _cpu(
                text(field(run, "workload", where), f"{where}.workload"),
                f"{family}/budget_run_generation",
                ClockSummary.single(number(field(run, "generation_s", where), f"{where}.generation_s")),
                run_load,
            )
        )
    return AdaptedRecord(_build_provenance(record, "machine"), tuple(measurements))


def _runs_with_load(runs: Sequence[Any], where: str) -> tuple[ClockSummary, tuple[Optional[LoadAverage1m], Optional[LoadAverage1m]]]:
    rows = [obj(run, f"{where}[{index}]") for index, run in enumerate(runs)]
    seconds = [number(field(row, "seconds", where), f"{where}.seconds") for row in rows]
    loads = [_load(row.get("load_1min"), f"{where}.load_1min") for row in rows]
    return ClockSummary.over_runs(seconds), _first_last(loads)


def coverage_certificate_acceptance(raw: Any, _sibling: SiblingReader) -> AdaptedRecord:
    """Per workload: the path generation (no load recorded) and the certificate runs."""
    record = obj(raw, "record")
    _require_schema(record, ("coverage-certificate-acceptance/1",), "coverage_certificate_acceptance")
    family = "coverage_certificate_acceptance"
    measurements: list[Measurement] = []
    for index, value in enumerate(items(field(record, "workloads", "record"), "workloads")):
        row = obj(value, f"workloads[{index}]")
        workload = text(field(row, "workload", "workloads[]"), "workload")
        measurements.append(_wall(workload, f"{family}/generation", ClockSummary.single(number(field(row, "generation_s", workload), f"{workload}.generation_s"))))
        clock, load = _runs_with_load(items(field(row, "runs", workload), f"{workload}.runs"), f"{workload}.runs")
        measurements.append(_wall(workload, f"{family}/certificate", clock, load))
    return AdaptedRecord(_build_provenance(record, "identity"), tuple(measurements))


# The circle-core tool's input: the frozen Monstera spike fixture, not a generated path.
CIRCLE_CORE_WORKLOAD = "monstera_circles_fixture"


def coverage_certificate_circle_core(raw: Any, _sibling: SiblingReader) -> AdaptedRecord:
    """The circle-circle core runs, and the audit run, on the frozen Monstera fixture."""
    record = obj(raw, "record")
    _require_schema(record, ("coverage-certificate-circle-core/1",), "coverage_certificate_circle_core")
    family = "coverage_certificate_circle_core"
    clock, load = _runs_with_load(items(field(record, "runs", "record"), "runs"), "runs")
    audit = obj(field(record, "audit", "record"), "audit")
    audit_load = _load(audit.get("load_1min"), "audit.load_1min")
    return AdaptedRecord(
        _build_provenance(record, "identity"),
        (
            _wall(CIRCLE_CORE_WORKLOAD, f"{family}/core", clock, load),
            _wall(CIRCLE_CORE_WORKLOAD, f"{family}/audit", ClockSummary.single(number(field(audit, "seconds", "audit"), "audit.seconds")), (audit_load, audit_load)),
        ),
    )


def coverage_certificate_cost_model(raw: Any, _sibling: SiblingReader) -> AdaptedRecord:
    """The one timing in the cost-model inputs: the path generation (no load recorded)."""
    record = obj(raw, "record")
    _require_schema(record, ("coverage-certificate-cost-model/1",), "coverage_certificate_cost_model")
    inputs = obj(field(record, "inputs", "record"), "inputs")
    workload = text(field(inputs, "workload", "inputs"), "inputs.workload")
    clock = ClockSummary.single(number(field(inputs, "generation_s", "inputs"), "inputs.generation_s"))
    return AdaptedRecord(_build_provenance(record, "identity"), (_wall(workload, "coverage_certificate_cost_model/generation", clock),))


# ``tools.coverage_certificate_omega_star.CASE``: the polygonal Monstera path the diagnostic certifies.
OMEGA_STAR_WORKLOAD = "figure8_monstera"
# The schema whose ``budget`` is the production workload's gate (``tools.coverage_certificate_budget_gate``).
OMEGA_STAR_BUDGET_SCHEMA = "coverage-certificate-omega-star/3"


def _omega_star_budget(record: Mapping[str, Any], family: str) -> Measurement:
    """The /3 production budget, read through ``ProductionBudget.from_raw`` (the gate's own validating loader).

    CPU and wall over the fresh-process runs; load before is the parent's reading
    before the first run, load after the worker's reading after the last certificate.
    """
    try:
        budget = ProductionBudget.from_raw(dict(obj(field(record, "budget", "record"), "budget")))
    except (BenchmarkError, MachineNotQuietError) as error:
        raise MalformedTimingRecordError(f"{family}: budget: {error}") from error
    load = (LoadAverage1m(budget.runs[0].load_1min), LoadAverage1m(budget.runs[-1].certificate_load_1min[1]))
    downgrade = "the record marks the budget runs indicative (measured without the quiet-machine check)" if budget.indicative else None
    return Measurement(
        budget.workload,
        f"{family}/budget",
        ClockSummary.over_runs([run.cpu_s for run in budget.runs]),
        ClockSummary.over_runs([run.wall_s for run in budget.runs]),
        load[0],
        load[1],
        downgrade,
    )


def coverage_certificate_omega_star(raw: Any, _sibling: SiblingReader) -> AdaptedRecord:
    """The diagnostic certificate runs (CPU and wall, load per run) and its one-off phases; from /3 also the production budget.

    /3 keeps the /2 diagnostic fields and adds ``budget``: the production
    workload's fresh-process runs (``_omega_star_budget``).
    """
    record = obj(raw, "record")
    schema = _require_schema(record, ("coverage-certificate-omega-star/1", "coverage-certificate-omega-star/2", OMEGA_STAR_BUDGET_SCHEMA), "coverage_certificate_omega_star")
    family = "coverage_certificate_omega_star"
    cpu = numbers(field(record, "runs_cpu_s", "record"), "runs_cpu_s")
    wall = numbers(field(record, "runs_wall_s", "record"), "runs_wall_s")
    loads = [LoadAverage1m(value) for value in numbers(field(record, "runs_load", "record"), "runs_load")]
    if not (len(cpu) == len(wall) == len(loads)):
        raise MalformedTimingRecordError(f"runs_cpu_s, runs_wall_s and runs_load differ in length: {len(cpu)}, {len(wall)}, {len(loads)}.")
    load = _first_last(loads)
    measurements = [
        Measurement(OMEGA_STAR_WORKLOAD, f"{family}/certificate", ClockSummary.over_runs(cpu), ClockSummary.over_runs(wall), load[0], load[1]),
        _wall(OMEGA_STAR_WORKLOAD, f"{family}/target_build", ClockSummary.single(number(field(record, "target_build_s", "record"), "target_build_s"))),
    ]
    if schema in ("coverage-certificate-omega-star/2", OMEGA_STAR_BUDGET_SCHEMA):
        measurements.append(_wall(OMEGA_STAR_WORKLOAD, f"{family}/attribution", ClockSummary.single(number(field(record, "attribution_s", "record"), "attribution_s"))))
    if schema == OMEGA_STAR_BUDGET_SCHEMA:
        measurements.append(_omega_star_budget(record, family))
    return AdaptedRecord(_build_provenance(record, "identity"), tuple(measurements))


def gate_vectors(raw: Any, _sibling: SiblingReader) -> AdaptedRecord:
    """Per gate case: the rendered generation and gate wall seconds (no load, no machine recorded)."""
    record = obj(raw, "record")
    _require_schema(record, ("gate-vectors/1", "gate-vectors/2"), "gate_vectors")
    measurements: list[Measurement] = []
    for key, value in obj(field(record, "cases", "record"), "cases").items():
        case = obj(value, key)
        for name in ("generation_s", "gate_s"):
            measurements.append(_wall(key, f"gate_vectors/{name[:-2]}", ClockSummary.single(number(field(case, name, key), f"{key}.{name}"))))
    return AdaptedRecord(RecordProvenance(), tuple(measurements))


# ``tools.held_centre_domain_timing.DETERMINISTIC_KEYS`` siblings: the timed phases and their reported clocks.
CENTRE_DOMAIN_PHASES = ("construction", "projection")


def held_centre_domain_timing(raw: Any, _sibling: SiblingReader) -> AdaptedRecord:
    """Per case: construction and projection, the median of ``samples`` fresh processes (CPU and wall); load read at the end only."""
    record = obj(raw, "record")
    samples = integer(field(record, "samples", "record"), "samples")
    after = _load(field(record, "load_average_1m_at_end", "record"), "load_average_1m_at_end")
    measurements: list[Measurement] = []
    for index, value in enumerate(items(field(record, "cases", "record"), "cases")):
        case = obj(value, f"cases[{index}]")
        workload = text(field(case, "case", "cases[]"), "case")
        for phase in CENTRE_DOMAIN_PHASES:
            cpu = ClockSummary.reported_median(number(field(case, f"{phase}_cpu_s", workload), f"{workload}.{phase}_cpu_s"), samples)
            wall = ClockSummary.reported_median(number(field(case, f"{phase}_wall_s", workload), f"{workload}.{phase}_wall_s"), samples)
            measurements.append(Measurement(workload, f"held_centre_domain_timing/{phase}", cpu, wall, None, after))
    return AdaptedRecord(RecordProvenance(), tuple(measurements))


CONNECTOR_ENGAGEMENT_STAGES: StageFields = {"generation": (None, "generation_seconds"), "tiers": (None, "tiers_seconds")}
REPAIR_BUDGET_STAGES: StageFields = {"generation_published": (None, "generation_seconds_published"), "generation_gate_off": (None, "generation_seconds_gate_off")}
OMEGA_CONNECTOR_GATE_STAGES: StageFields = {
    "generation": (None, "generation_seconds"),
    "path_build": (None, "path_build_seconds"),
    "assessment": (None, "assessment_seconds"),
    "gate": (None, "gate_seconds"),
    "repair_stage_gate": (None, "repair_stage_gate_seconds"),
}
# ``generation`` is clocked around the generator call (CPU and wall); ``generation_stage_sum``
# is the path's own ``generation_seconds``, the sum of its stage clocks.
OMEGA_CONNECTOR_REPAIR_STAGES: StageFields = {
    "generation": ("generation_cpu_seconds", "generation_wall_seconds"),
    "generation_stage_sum": (None, "generation_seconds"),
    "connector_repair": (None, "connector_repair_seconds"),
    "repair_gate": (None, "repair_gate_seconds"),
    "repair_pass": (None, "repair_pass_seconds"),
}
OMEGA_COVERAGE_STAGES: StageFields = {"generation": (None, "generation_seconds"), "certificate": ("certificate_cpu_seconds", "certificate_wall_seconds")}


def held_connector_engagement(raw: Any, _sibling: SiblingReader) -> AdaptedRecord:
    """Per workload: generation and the connector tiers (no load, no provenance recorded)."""
    measurements: list[Measurement] = []
    for index, value in enumerate(items(raw, "record")):
        row = obj(value, f"record[{index}]")
        workload = text(field(row, "workload", f"record[{index}]"), "workload")
        measurements.extend(_stage_fields(row, workload, "held_connector_engagement", CONNECTOR_ENGAGEMENT_STAGES, (None, None)))
    return AdaptedRecord(RecordProvenance(), tuple(measurements))


def held_connector_repair_budget(raw: Any, _sibling: SiblingReader) -> AdaptedRecord:
    """Per workload: generation published and with the repair gate off; the tier-2 rate run's phase seconds."""
    record = obj(raw, "record")
    family = "held_connector_repair_budget"
    measurements: list[Measurement] = []
    for index, value in enumerate(items(field(record, "workloads", "record"), "workloads")):
        row = obj(value, f"workloads[{index}]")
        workload = text(field(row, "workload", "workloads[]"), "workload")
        measurements.extend(
            _stage_fields(row, workload, family, REPAIR_BUDGET_STAGES, _load_pair(field(row, "load_1min_before_after", workload), f"{workload}.load_1min_before_after"))
        )
    rates = obj(field(record, "tier2_rates", "record"), "tier2_rates")
    workload = text(field(rates, "workload", "tier2_rates"), "tier2_rates.workload")
    load = _load_pair(field(rates, "load_1min_before_after", "tier2_rates"), "tier2_rates.load_1min_before_after")
    for name, seconds in sorted(obj(field(rates, "seconds", "tier2_rates"), "tier2_rates.seconds").items()):
        measurements.append(_wall(workload, f"{family}/tier2_{name}", ClockSummary.single(number(seconds, f"tier2_rates.seconds.{name}")), load))
    return AdaptedRecord(_build_provenance(record, "machine"), tuple(measurements))


def held_curved_circle_baseline(raw: Any, _sibling: SiblingReader) -> AdaptedRecord:
    """Per case: the sum of the per-piece query medians (each over ``repeats``), and the import (no load recorded)."""
    record = obj(raw, "record")
    if record.get("schema_version") != 1:
        raise UnknownRecordSchemaError(f"held_curved_circle_baseline: schema_version {record.get('schema_version')!r}; the adapter reads 1.")
    family = "held_curved_circle_baseline"
    measurements: list[Measurement] = []
    for workload, value in obj(field(record, "cases", "record"), "cases").items():
        case = obj(value, workload)
        pieces = [obj(piece, f"{workload}.pieces[]") for piece in items(field(case, "pieces", workload), f"{workload}.pieces")]
        if pieces:
            total = sum(number(field(piece, "seconds", workload), f"{workload}.pieces[].seconds") for piece in pieces)
            measurements.append(_wall(workload, f"{family}/queries", ClockSummary.total(total, len(pieces), ClockStatistic.SUM_OF_MEDIANS)))
        measurements.append(_wall(workload, f"{family}/import", ClockSummary.single(number(field(case, "import_seconds", workload), f"{workload}.import_seconds"))))
    return AdaptedRecord(_build_provenance(record, "identity"), tuple(measurements))


def held_native_clearance(raw: Any, _sibling: SiblingReader) -> AdaptedRecord:
    """Per workload and clearance backend: the generation runs (load per run)."""
    record = obj(raw, "record")
    measurements: list[Measurement] = []
    for index, value in enumerate(items(field(record, "generation", "record"), "generation")):
        row = obj(value, f"generation[{index}]")
        workload = text(field(row, "workload", "generation[]"), "workload")
        by_backend: dict[str, list[Mapping[str, Any]]] = {}
        for run_value in items(field(row, "runs", workload), f"{workload}.runs"):
            run = obj(run_value, f"{workload}.runs[]")
            by_backend.setdefault(text(field(run, "backend", workload), f"{workload}.runs[].backend"), []).append(run)
        for backend, runs in sorted(by_backend.items()):
            seconds = [number(field(run, "seconds", workload), f"{workload}.runs[].seconds") for run in runs]
            first = _load_pair(field(runs[0], "load_1min_before_after", workload), f"{workload}.runs[0].load_1min_before_after")
            last = _load_pair(field(runs[-1], "load_1min_before_after", workload), f"{workload}.runs[-1].load_1min_before_after")
            measurements.append(_wall(workload, f"held_native_clearance/generation[{backend}]", ClockSummary.over_runs(seconds), (first[0], last[1])))
    return AdaptedRecord(_build_provenance(record, "machine"), tuple(measurements))


def _held_workload_rows(family: str, spec: StageFields) -> Adapter:
    """An adapter for the Ω* Held records: per workload, the ``spec`` seconds fields, load before/after per workload."""

    def adapt(raw: Any, _sibling: SiblingReader) -> AdaptedRecord:
        record = obj(raw, "record")
        measurements: list[Measurement] = []
        for index, value in enumerate(items(field(record, "workloads", "record"), "workloads")):
            row = obj(value, f"workloads[{index}]")
            workload = text(field(row, "workload", "workloads[]"), "workload")
            measurements.extend(_stage_fields(row, workload, family, spec, _load_pair(field(row, "load_1min_before_after", workload), f"{workload}.load_1min_before_after")))
        return AdaptedRecord(_build_provenance(record, "machine"), tuple(measurements))

    adapt.__doc__ = f"Per workload of ``{family}``: stages {list(spec)}, load before/after."
    return adapt


def held_splice_ab(raw: Any, _sibling: SiblingReader) -> AdaptedRecord:
    """The full refinement, the repair round, and the summed A and B candidate seconds (one load window for the record)."""
    record = obj(raw, "record")
    family = "held_splice_ab"
    workload = text(field(record, "workload", "record"), "workload")
    load = _load_pair(field(record, "load_1min_before_after", "record"), "load_1min_before_after")
    candidates = [obj(value, "candidates[]") for value in items(field(record, "candidates", "record"), "candidates")]
    measurements = [
        _wall(workload, f"{family}/full_refinement", ClockSummary.single(number(field(record, "full_refinement_seconds", "record"), "full_refinement_seconds")), load),
        _wall(workload, f"{family}/repair_round", ClockSummary.single(number(field(record, "repair_round_seconds", "record"), "repair_round_seconds")), load),
    ]
    if candidates:
        for side in ("a", "b"):
            total = sum(number(field(candidate, f"{side}_seconds", "candidates[]"), f"candidates[].{side}_seconds") for candidate in candidates)
            measurements.append(_wall(workload, f"{family}/candidates_{side}", ClockSummary.total(total, len(candidates), ClockStatistic.SUM_OF_SINGLES), load))
    return AdaptedRecord(_build_provenance(record, "machine"), tuple(measurements))


def held_toolpath_baseline(raw: Any, _sibling: SiblingReader) -> AdaptedRecord:
    """Per workload: the generation and each generator stage, one run (no load recorded)."""
    record = obj(raw, "record")
    if record.get("schema_version") != 1:
        raise UnknownRecordSchemaError(f"held_toolpath_baseline: schema_version {record.get('schema_version')!r}; the adapter reads 1.")
    family = "held_toolpath_baseline"
    measurements: list[Measurement] = []
    for workload, value in obj(field(record, "runs", "record"), "runs").items():
        run = obj(value, workload)
        measurements.append(_wall(workload, f"{family}/generation", ClockSummary.single(number(field(run, "generation_seconds", workload), f"{workload}.generation_seconds"))))
        for stage, seconds in obj(field(run, "stage_seconds", workload), f"{workload}.stage_seconds").items():
            measurements.append(_wall(workload, f"{family}/stage_{stage}", ClockSummary.single(number(seconds, f"{workload}.stage_seconds.{stage}"))))
        validation = run.get("validation_seconds")
        if validation is not None:
            measurements.append(_wall(workload, f"{family}/validation", ClockSummary.single(number(validation, f"{workload}.validation_seconds"))))
    return AdaptedRecord(_build_provenance(record, "identity"), tuple(measurements))


# The one ``pipeline_timing`` schema the core's reference reader reads; the registered
# adapter (``tools.perf_history_adapters_pipeline``) reads /1 and /2 (load + provenance).
PIPELINE_TIMING_CORE_SCHEMA = "pipeline-timing/1"
# Where the registered pipeline adapters live.
PIPELINE_ADAPTERS_MODULE = "tools.perf_history_adapters_pipeline"


def pipeline_timing(raw: Any, _sibling: SiblingReader) -> AdaptedRecord:
    """Schema /1 only: per gate case, the fresh-process generation runs (no load, no provenance), read through ``PipelineTimings.build``.

    Not registered: ``FAMILIES`` reads ``pipeline_timing.json`` with ``PIPELINE_ADAPTERS_MODULE``'s
    adapter. This reader wrote the history's /1 entries before that registration and is the
    reference the registered adapter must reproduce on /1 records
    (``tests/tools/test_perf_history_pipeline.py``).

    Raises:
        UnknownRecordSchemaError: any other schema; /2 is read by ``PIPELINE_ADAPTERS_MODULE``.
    """
    record = obj(raw, "record")
    if record.get("schema") != PIPELINE_TIMING_CORE_SCHEMA:
        raise UnknownRecordSchemaError(
            f"pipeline_timing: schema {record.get('schema')!r}; the core reads only {PIPELINE_TIMING_CORE_SCHEMA!r}. "
            f"/2 is read by {PIPELINE_ADAPTERS_MODULE}, registered in tools.perf_history_adapters.FAMILIES."
        )
    try:
        timings = PipelineTimings.build(raw)
    except MalformedPipelineTimingsError as error:
        raise MalformedTimingRecordError(f"pipeline_timing: {error}") from error
    measurements = tuple(Measurement(key, "pipeline_timing/generation", None, ClockSummary.over_runs(case.runs_s), None, None) for key, case in timings.cases.items())
    return AdaptedRecord(RecordProvenance(), measurements)


def benchmark_report(raw: Any, sibling: SiblingReader) -> AdaptedRecord:
    """A ``benchmarks.cli corpus`` report: per case, generation and certification (one run, no load), provenance from the sibling ``stamp.json``."""
    stamp = obj(json.loads(sibling(CORPUS_STAMP)), CORPUS_STAMP)
    provenance = RecordProvenance(
        commit=optional_text(stamp, "commit", CORPUS_STAMP),
        dirty=optional_bool(stamp, "dirty", CORPUS_STAMP),
        recorded_at=optional_text(stamp, "started", CORPUS_STAMP),
    )
    measurements: list[Measurement] = []
    for index, value in enumerate(items(raw, "record")):
        case = obj(value, f"record[{index}]")
        workload = text(field(case, "name", f"record[{index}]"), "name")
        if case.get("error") is not None:
            raise MalformedTimingRecordError(f"benchmark_report: {workload} recorded an error {case['error']!r}; its seconds are not a measurement.")
        for name in ("generate", "certify"):
            measurements.append(_wall(workload, f"benchmark_report/{name}", ClockSummary.single(number(field(case, f"{name}_seconds", workload), f"{workload}.{name}_seconds"))))
    return AdaptedRecord(provenance, tuple(measurements))


# Record file name -> its adapter; the pipeline families come from ``PIPELINE_FAMILIES``.
FAMILIES: Mapping[str, Adapter] = {
    "actual_stock_profile.json": actual_stock_profile,
    "benchmark_report.json": benchmark_report,
    "coverage_certificate_acceptance.json": coverage_certificate_acceptance,
    "coverage_certificate_circle_core.json": coverage_certificate_circle_core,
    "coverage_certificate_cost_model.json": coverage_certificate_cost_model,
    "coverage_certificate_omega_star.json": coverage_certificate_omega_star,
    "gate_vectors.json": gate_vectors,
    "held_centre_domain_timing.json": held_centre_domain_timing,
    "held_connector_engagement.json": held_connector_engagement,
    "held_connector_repair_budget.json": held_connector_repair_budget,
    "held_curved_circle_baseline.json": held_curved_circle_baseline,
    "held_native_clearance.json": held_native_clearance,
    "held_omega_connector_gate.json": _held_workload_rows("held_omega_connector_gate", OMEGA_CONNECTOR_GATE_STAGES),
    "held_omega_connector_repair.json": _held_workload_rows("held_omega_connector_repair", OMEGA_CONNECTOR_REPAIR_STAGES),
    "held_omega_coverage.json": _held_workload_rows("held_omega_coverage", OMEGA_COVERAGE_STAGES),
    "held_splice_ab.json": held_splice_ab,
    "held_toolpath_baseline.json": held_toolpath_baseline,
} | dict(PIPELINE_FAMILIES)

# Families whose adapter exists but is not registered yet: record file name -> where it is.
# Registration is explicit, in this module and nowhere else: import the adapter module's
# families, merge them into ``FAMILIES`` and remove them here. A pending record raises
# ``PendingRecordFamilyError``. Empty since the pipeline families were registered.
PENDING_FAMILIES: Mapping[str, str] = {}

# Family -> the version of its adapter. Bump a family's version whenever its adapter
# changes what it reads from an unchanged record (a clock relabelled, a stage
# renamed or added): re-running ``backfill`` then appends the corrected entries
# beside the old ones (the history is never rewritten), and readers keep, per
# record version, only the entries of the newest adapter version
# (``tools.perf_history_baseline.current_entries``).
ADAPTER_VERSIONS: Mapping[str, int] = {name[: -len(".json")]: 1 for name in FAMILIES} | {
    # 2: generation_s is process CPU, not wall (RV40 I1, 2026-10-23).
    "actual_stock_profile": 2,
}

# Record file name -> why it holds no performance measurement.
NON_TIMING_RECORDS: Mapping[str, str] = {
    "coverage_certificate_raster.json": "raster coverage areas; no timing",
    "held_boundary_census.json": "boundary primitive census; no timing",
    "held_connector_engagement_tier2.json": "per-connector engagement bounds; no timing",
    "held_declared_boundary.json": "declared-boundary census; no timing",
    "held_path_digests.json": "path digests; no timing",
    "performance_lever_shares.json": "profiler time shares (cProfile-inflated, percentages); not a clocked run",
    "performance_levers.json": "lever-to-function mapping; an input, no timing",
    "performance_time_buckets.json": "profiler time buckets (cProfile-inflated totals); not a clocked run",
    "perf_history.jsonl": "this history itself",
    "perf_ratio_session.json": "same-moment HEAD/reference CPU ratios; owned by tools.perf_ratio_history, never the absolute history",
    CORPUS_STAMP: "provenance envelope of a corpus report; read by the benchmark_report adapter",
}


def family_of(record_path: str) -> str:
    """The family of the record at ``record_path`` (its file name).

    Raises:
        NotATimingRecordError: the file is a listed non-timing record.
        PendingRecordFamilyError: the family's adapter is not registered yet.
        UnknownRecordFamilyError: the file name is none of these.
    """
    name = PurePosixPath(record_path).name
    if name in FAMILIES:
        return name[: -len(".json")]
    if name in NON_TIMING_RECORDS:
        raise NotATimingRecordError(f"{record_path}: {NON_TIMING_RECORDS[name]}.")
    if name in PENDING_FAMILIES:
        raise PendingRecordFamilyError(f"{record_path}: {PENDING_FAMILIES[name]}.")
    raise UnknownRecordFamilyError(f"{record_path}: no timing-record adapter for {name!r} and not a listed non-timing record; add one to tools.perf_history_adapters.")


def adapt(record_path: str, data: bytes, sibling: SiblingReader) -> AdaptedRecord:
    """Every measured stage of the record ``data`` read from ``record_path``.

    Raises:
        NotATimingRecordError: the file is a listed non-timing record.
        UnknownRecordFamilyError: its file name has no adapter.
        UnknownRecordSchemaError: its schema is not one its adapter reads.
        MalformedTimingRecordError: it is not the shape its adapter reads.
        MissingSiblingRecordError: its adapter needs a file beside it that is absent.
    """
    family_of(record_path)
    try:
        raw = json.loads(data)
    except json.JSONDecodeError as error:
        raise MalformedTimingRecordError(f"{record_path}: not JSON: {error}") from error
    return FAMILIES[PurePosixPath(record_path).name](raw, sibling)
