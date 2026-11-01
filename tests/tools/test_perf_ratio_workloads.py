"""The ratio workload registry, child commands, child-output parsing and the GT CPU clock."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any
from typing import Callable

import pytest

from tools.coverage_certificate_budget_gate import BUDGET_WORKER_MARKER
from tools.held_generation_clock import WORKER_MARKER
from tools.perf_ratio_record import CpuSeconds
from tools.perf_ratio_workloads import CERTIFICATE_WORKLOAD
from tools.perf_ratio_workloads import ChildWorkloadMismatchError
from tools.perf_ratio_workloads import DEFAULT_SET_NAME
from tools.perf_ratio_workloads import GT_CAP40_KEYS
from tools.perf_ratio_workloads import GT_CAP40_SET_NAME
from tools.perf_ratio_workloads import GT_CAP120_KEYS
from tools.perf_ratio_workloads import HELD_WORKLOADS
from tools.perf_ratio_workloads import WORKLOAD_SETS
from tools.perf_ratio_workloads import Family
from tools.perf_ratio_workloads import InjectionUnsupportedError
from tools.perf_ratio_workloads import MalformedChildResultError
from tools.perf_ratio_workloads import MissingChildResultError
from tools.perf_ratio_workloads import UnknownWorkloadError
from tools.perf_ratio_workloads import Workload
from tools.perf_ratio_workloads import child_command
from tools.perf_ratio_workloads import parse_child
from tools.perf_ratio_workloads import workload_named
from tools.pipeline_timing import RESULT_LINE_PREFIX
from tools.pipeline_timing import gate_timing_keys

# Process CPU the injection probe spins, s: well above process_time's resolution, short enough for a unit test.
INJECT_PROBE_S = 0.2
# Iterations of the stub generator's busy loop: enough arithmetic for a non-zero process_time span.
STUB_BUSY_ITERATIONS = 200_000

_GT_KEY = GT_CAP120_KEYS[0]
_MARKERS = {
    Family.POLYGON: WORKER_MARKER,
    Family.OMEGA: WORKER_MARKER,
    Family.CERTIFICATE: BUDGET_WORKER_MARKER,
    Family.GT: RESULT_LINE_PREFIX,
}
_ONE_PER_FAMILY = (
    Workload(Family.POLYGON, HELD_WORKLOADS[0]),
    Workload(Family.OMEGA, HELD_WORKLOADS[0]),
    Workload(Family.CERTIFICATE, CERTIFICATE_WORKLOAD),
    Workload(Family.GT, _GT_KEY),
)


def _stub_pipeline_timing(tree: Path, clock_reads: int) -> None:
    """A stub ``tools/pipeline_timing.py`` whose ``main`` busy-loops between its clock reads."""
    tools = tree / "tools"
    tools.mkdir()
    (tools / "__init__.py").write_text("")
    extra_reads = "\n".join("    time.perf_counter()" for _ in range(clock_reads - 2))
    (tools / "pipeline_timing.py").write_text(
        "import json\n"
        "import time\n"
        "def main(argv):\n"
        "    started = time.perf_counter()\n"
        f"    total = sum(i * i for i in range({STUB_BUSY_ITERATIONS}))\n"
        f"{extra_reads}\n"
        "    elapsed = time.perf_counter() - started\n"
        f"    print('noise before the result line')\n"
        f"    print({RESULT_LINE_PREFIX!r} + json.dumps({{'generation_s': elapsed, 'operations': 7, 'digest': 'd', 'total': total}}))\n"
        "    return 0\n"
    )


def _run_child(tree: Path, workload: Workload, inject: float) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        child_command(workload, (sys.executable,), CpuSeconds(inject)),
        cwd=tree,
        env={**os.environ, "PYTHONPATH": str(tree)},
        capture_output=True,
        text=True,
        check=False,
    )


def _certificate_row(**overrides: object) -> dict[str, object]:
    row: dict[str, object] = {
        "workload": CERTIFICATE_WORKLOAD,
        "circles": 10,
        "connectors": 4,
        "sub_pieces": 6,
        "omega_pieces": 3,
        "vertex_disks": 2,
        "repair_inserted": 1,
        "generation_seconds": 1.5,
        "verdict": "COVERED",
        "witnesses": 0,
        "witnesses_by_kind": {},
        "undecided_edges": 0,
        "regions": 13,
        "pieces": 20,
        "edges": 40,
        "certificate_wall_seconds": 2.5,
        "certificate_cpu_seconds": 2.0,
        "load_1min_before_after": [3.0, 4.0],
    }
    row.update(overrides)
    return row


def test_the_default_set_is_the_specs_ten_workloads() -> None:
    expected = (
        [f"polygon/{w}" for w in HELD_WORKLOADS]
        + [f"omega/{w}" for w in HELD_WORKLOADS]
        + ["certificate/figure8_monstera_boundary_qualified_corners"]
        + [f"gt/{k}" for k in GT_CAP120_KEYS]
    )
    assert [w.name for w in WORKLOAD_SETS[DEFAULT_SET_NAME]] == expected
    assert len(expected) == 10
    assert WORKLOAD_SETS[GT_CAP40_SET_NAME] == tuple(Workload(Family.GT, k) for k in GT_CAP40_KEYS)


def test_the_gt_keys_are_real_gate_cases() -> None:
    assert set(GT_CAP120_KEYS) | set(GT_CAP40_KEYS) <= set(gate_timing_keys()[0])


def test_held_and_certificate_commands_are_the_tools_worker_forms() -> None:
    assert child_command(Workload(Family.OMEGA, "figure5_refined"), ("py",), CpuSeconds(0.0)) == (
        "py",
        "-m",
        "tools.held_generation_clock",
        "--domain",
        "omega",
        "figure5_refined",
    )
    assert child_command(Workload(Family.POLYGON, "figure5_refined"), ("py",), CpuSeconds(0.0)) == (
        "py",
        "-m",
        "tools.held_generation_clock",
        "--domain",
        "polygon",
        "figure5_refined",
    )
    assert child_command(Workload(Family.CERTIFICATE, CERTIFICATE_WORKLOAD), ("py",), CpuSeconds(0.0)) == (
        "py",
        "-m",
        "tools.coverage_certificate_budget_gate",
        "--worker",
    )


def test_injection_is_refused_outside_gt() -> None:
    for workload in _ONE_PER_FAMILY:
        if workload.family is Family.GT:
            continue
        with pytest.raises(InjectionUnsupportedError):
            child_command(workload, ("py",), CpuSeconds(1.0))


def test_an_unregistered_workload_has_no_child_command() -> None:
    # The certificate worker always certifies its own production workload, so an argument it
    # ignores would be measured under the wrong name.
    for workload in (Workload(Family.CERTIFICATE, "figure5_refined"), Workload(Family.POLYGON, "nope"), Workload(Family.GT, "gt-nope")):
        with pytest.raises(UnknownWorkloadError):
            child_command(workload, ("py",), CpuSeconds(0.0))


def test_the_gt_snippet_reports_cpu_from_the_arms_own_module(tmp_path: Path) -> None:
    _stub_pipeline_timing(tmp_path, clock_reads=2)
    workload = Workload(Family.GT, _GT_KEY)
    completed = _run_child(tmp_path, workload, 0.0)
    assert completed.returncode == 0, completed.stderr
    run = parse_child(workload, completed.stdout)
    assert run.digest == "d:7"
    assert run.cpu_seconds > 0.0


def test_the_gt_snippet_spins_the_injected_cpu_inside_the_span(tmp_path: Path) -> None:
    _stub_pipeline_timing(tmp_path, clock_reads=2)
    workload = Workload(Family.GT, _GT_KEY)
    completed = _run_child(tmp_path, workload, INJECT_PROBE_S)
    assert completed.returncode == 0, completed.stderr
    assert parse_child(workload, completed.stdout).cpu_seconds >= INJECT_PROBE_S


def test_a_moved_clock_fails_loud(tmp_path: Path) -> None:
    _stub_pipeline_timing(tmp_path, clock_reads=3)
    completed = _run_child(tmp_path, Workload(Family.GT, _GT_KEY), 0.0)
    assert completed.returncode != 0
    assert "not 2" in completed.stderr


def test_zero_or_two_result_lines_are_refused() -> None:
    for workload in _ONE_PER_FAMILY:
        line = _MARKERS[workload.family] + "{}"
        for stdout in ("log noise only\n", "", f"{line}\nnoise\n{line}\n"):
            with pytest.raises(MissingChildResultError):
                parse_child(workload, stdout)


def test_the_certificate_digest_ignores_timing_and_load() -> None:
    workload = Workload(Family.CERTIFICATE, CERTIFICATE_WORKLOAD)

    def digest_of(row: dict[str, object]) -> str:
        return parse_child(workload, f"noise\n{BUDGET_WORKER_MARKER}{json.dumps(row)}\n").digest

    base = digest_of(_certificate_row())
    retimed = digest_of(_certificate_row(certificate_cpu_seconds=9.0, certificate_wall_seconds=9.5, generation_seconds=7.0, load_1min_before_after=[30.0, 40.0]))
    assert base == retimed
    assert digest_of(_certificate_row(edges=41)) != base
    assert parse_child(workload, f"{BUDGET_WORKER_MARKER}{json.dumps(_certificate_row())}").cpu_seconds == 2.0


def test_an_unknown_workload_name_is_refused() -> None:
    with pytest.raises(UnknownWorkloadError):
        workload_named("gt/nope")
    assert workload_named(f"gt/{GT_CAP40_KEYS[1]}") == Workload(Family.GT, GT_CAP40_KEYS[1])


def _held_payload() -> dict[str, Any]:
    clock = {"stage_cpu_seconds": {"stage": 1.0}, "stage_wall_seconds": {"stage": 1.25}, "exit_cpu_seconds": 0.5, "exit_wall_seconds": 0.75}
    return {"digest": {"operations": 3}, "clock": clock}


def _gt_payload() -> dict[str, Any]:
    return {"generation_s": 1.5, "operations": 7, "digest": "d"}


def _without(key: str) -> Callable[[dict[str, Any]], None]:
    return lambda payload: payload.pop(key)


def _with(key: str, value: object) -> Callable[[dict[str, Any]], None]:
    return lambda payload: payload.__setitem__(key, value)


def _without_clock_key(payload: dict[str, Any]) -> None:
    payload["clock"].pop("exit_cpu_seconds")


def _with_clock_value(payload: dict[str, Any]) -> None:
    payload["clock"]["exit_cpu_seconds"] = "slow"


_FAMILY_PAYLOADS: dict[Family, Callable[[], dict[str, Any]]] = {
    Family.POLYGON: _held_payload,
    Family.OMEGA: _held_payload,
    Family.CERTIFICATE: _certificate_row,
    Family.GT: _gt_payload,
}
# (family, case, mutation of a valid payload, whether the refusal chains an original exception)
_MALFORMED_PAYLOADS = [
    *((family, "missing clock", _without("clock"), True) for family in (Family.POLYGON, Family.OMEGA)),
    *((family, "missing clock field", _without_clock_key, True) for family in (Family.POLYGON, Family.OMEGA)),
    *((family, "missing digest", _without("digest"), True) for family in (Family.POLYGON, Family.OMEGA)),
    *((family, "clock field not a number", _with_clock_value, True) for family in (Family.POLYGON, Family.OMEGA)),
    *((family, "clock not an object", _with("clock", [1.0]), False) for family in (Family.POLYGON, Family.OMEGA)),
    *((family, "digest not an object", _with("digest", "abc"), False) for family in (Family.POLYGON, Family.OMEGA)),
    (Family.CERTIFICATE, "missing cpu", _without("certificate_cpu_seconds"), True),
    (Family.CERTIFICATE, "missing digest field", _without("edges"), True),
    (Family.CERTIFICATE, "cpu a string", _with("certificate_cpu_seconds", "2.0"), False),
    (Family.CERTIFICATE, "cpu a boolean", _with("certificate_cpu_seconds", True), False),
    (Family.CERTIFICATE, "cpu not finite", _with("certificate_cpu_seconds", float("inf")), False),
    (Family.GT, "missing cpu", _without("generation_s"), True),
    (Family.GT, "missing digest", _without("digest"), True),
    (Family.GT, "missing operations", _without("operations"), True),
    (Family.GT, "cpu a string", _with("generation_s", "1.5"), False),
    (Family.GT, "cpu not finite", _with("generation_s", float("nan")), False),
    (Family.GT, "digest not a string", _with("digest", 5), False),
    (Family.GT, "operations not an integer", _with("operations", 7.0), False),
]


def _workload_of(family: Family) -> Workload:
    return next(workload for workload in _ONE_PER_FAMILY if workload.family is family)


@pytest.mark.parametrize("family", list(Family), ids=lambda family: family.value)
def test_a_valid_payload_parses_for_every_family(family: Family) -> None:
    run = parse_child(_workload_of(family), _MARKERS[family] + json.dumps(_FAMILY_PAYLOADS[family]()))
    assert run.cpu_seconds > 0.0
    assert run.digest


@pytest.mark.parametrize("family", list(Family), ids=lambda family: family.value)
@pytest.mark.parametrize("payload", ["{not json", "[1, 2]", "null"], ids=["bad-json", "array", "null"])
def test_a_non_object_payload_is_malformed(family: Family, payload: str) -> None:
    with pytest.raises(MalformedChildResultError) as refused:
        parse_child(_workload_of(family), _MARKERS[family] + payload)
    if payload == "{not json":
        assert isinstance(refused.value.__cause__, json.JSONDecodeError)


@pytest.mark.parametrize(
    ("family", "mutate", "chained"),
    [(family, mutate, chained) for family, _case, mutate, chained in _MALFORMED_PAYLOADS],
    ids=[f"{family.value}-{case}" for family, case, _mutate, _chained in _MALFORMED_PAYLOADS],
)
def test_a_malformed_payload_is_refused(family: Family, mutate: Callable[[dict[str, Any]], None], chained: bool) -> None:
    payload = _FAMILY_PAYLOADS[family]()
    mutate(payload)
    with pytest.raises(MalformedChildResultError) as refused:
        parse_child(_workload_of(family), _MARKERS[family] + json.dumps(payload))
    assert (refused.value.__cause__ is not None) == chained


def test_a_certificate_row_for_another_workload_is_refused() -> None:
    # The worker certifies the arm's own production workload whatever the harness asked for.
    workload = Workload(Family.CERTIFICATE, CERTIFICATE_WORKLOAD)
    with pytest.raises(ChildWorkloadMismatchError):
        parse_child(workload, BUDGET_WORKER_MARKER + json.dumps(_certificate_row(workload="figure5_refined")))
