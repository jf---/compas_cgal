"""The ``perf-ratio/1`` session record: round-trip, derived values, and its named refusals."""

from __future__ import annotations

import json
import math
import statistics
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest
from tools.perf_history_json import MalformedTimingRecordError
from tools.perf_ratio_record import ArmIdentity
from tools.perf_ratio_record import ArmName
from tools.perf_ratio_record import CpuSeconds
from tools.perf_ratio_record import InconsistentRatioRecordError
from tools.perf_ratio_record import LoadOneMinute
from tools.perf_ratio_record import MalformedRatioRecordError
from tools.perf_ratio_record import NonPositiveCpuError
from tools.perf_ratio_record import PairRun
from tools.perf_ratio_record import PerfRatioSession
from tools.perf_ratio_record import WorkloadRatio
from tools.perf_ratio_record import arms_to_json
from tools.perf_ratio_record import build_arms
from tools.perf_ratio_record import load_session
from tools.perf_ratio_record import stage_of

REF_COMMIT = "15fbed70ebb016bdd8ef49453e04606d03e0b1fe"
HEAD_COMMIT = "0123456789abcdef0123456789abcdef01234567"
# (ref, head) CPU seconds per pair, chosen so min, median and pair-ratio median all differ.
PAIR_CPU = ((2.0, 3.0), (1.0, 4.0), (4.0, 2.0))


def _arm(rev: str, name: str) -> ArmIdentity:
    return ArmIdentity(
        rev=rev,
        dirty=False,
        patch_sha256=None,
        native_sha256=f"native-{name}",
        native_sources_sha256=f"sources-{name}",
        environment="default",
        python_version="3.12.7",
        stock_2_sha256=f"stock-{name}",
        pixi_lock_sha256=f"lock-{name}",
        compas_cgal_file=f"/trees/{name}/src/compas_cgal/__init__.py",
        stock_2_origin=f"/trees/{name}/build/_stock_2.so",
        editable_build_dir=f"/trees/{name}/build",
        interpreter_prefix=f"/trees/{name}/.pixi/envs/default",
    )


def _pairs() -> tuple[PairRun, ...]:
    orders = ((ArmName.REF, ArmName.HEAD), (ArmName.HEAD, ArmName.REF))
    return tuple(
        PairRun.build(
            order=orders[index % 2],
            ref_cpu=CpuSeconds(ref),
            head_cpu=CpuSeconds(head),
            load_before=LoadOneMinute(100.0 + index),
            load_after=LoadOneMinute(110.0 + index),
        )
        for index, (ref, head) in enumerate(PAIR_CPU)
    )


def _session(head_digest: str = "digest-a") -> PerfRatioSession:
    workload = WorkloadRatio(workload="gt-cap120-case-1", family="gt", pairs=_pairs(), ref_digest="digest-a", head_digest=head_digest)
    return PerfRatioSession(
        ref_commit=REF_COMMIT,
        arms={ArmName.REF: _arm(REF_COMMIT, "ref"), ArmName.HEAD: _arm(HEAD_COMMIT, "head")},
        workloads=(workload,),
        host="bench-host",
        cpu_model="Apple M1 Max",
        started_utc="2026-11-01T10:00:00+00:00",
        ended_utc="2026-11-01T10:30:00+00:00",
        harness_sha256="harness-digest",
        injected=None,
    )


def _raw(session: PerfRatioSession) -> Any:
    return json.loads(session.dumps())


def test_a_session_round_trips_through_json() -> None:
    session = _session()
    assert PerfRatioSession.build(json.loads(session.dumps())) == session


def test_the_arms_and_a_workload_round_trip_through_their_public_builders() -> None:
    session = _session()
    assert build_arms(json.loads(json.dumps(arms_to_json(session.arms)))) == session.arms
    workload = session.workloads[0]
    assert WorkloadRatio.build(json.loads(json.dumps(workload.to_json()))) == workload
    assert ArmIdentity.build(session.arms[ArmName.HEAD].to_json()) == session.arms[ArmName.HEAD]
    assert stage_of(REF_COMMIT) == session.stage == "ratio-vs-15fbed70"


def test_the_public_builders_refuse_with_the_record_errors() -> None:
    session = _session()
    with pytest.raises(MalformedRatioRecordError):
        build_arms({"ref": session.arms[ArmName.REF].to_json()})
    with pytest.raises(MalformedRatioRecordError):
        ArmIdentity.build("an arm")
    tampered = dict(session.workloads[0].to_json(), ratio=0.25)
    with pytest.raises(InconsistentRatioRecordError):
        WorkloadRatio.build(tampered)
    with pytest.raises(MalformedRatioRecordError):
        WorkloadRatio.build([])


def test_ratio_is_min_head_over_min_ref_and_pair_median_beside_it() -> None:
    workload = _session().workloads[0]
    assert workload.ratio == 2.0 / 1.0
    assert workload.pair_ratio_median == statistics.median([1.5, 4.0, 0.5])
    assert workload.ref_min == 1.0
    assert workload.head_median == 3.0


def test_a_tampered_derived_value_is_refused() -> None:
    raw = _raw(_session())
    raw["workloads"][0]["ratio"] = raw["workloads"][0]["ratio"] * 0.5
    with pytest.raises(InconsistentRatioRecordError):
        PerfRatioSession.build(raw)


def test_a_non_positive_cpu_run_is_refused() -> None:
    order = (ArmName.REF, ArmName.HEAD)
    load = LoadOneMinute(1.0)
    with pytest.raises(NonPositiveCpuError):
        PairRun.build(order=order, ref_cpu=CpuSeconds(0.0), head_cpu=CpuSeconds(1.0), load_before=load, load_after=load)
    with pytest.raises(NonPositiveCpuError):
        PairRun.build(order=order, ref_cpu=CpuSeconds(1.0), head_cpu=CpuSeconds(-1.0), load_before=load, load_after=load)
    raw = _raw(_session())
    raw["workloads"][0]["pairs"][0]["ref_cpu"] = 0.0
    with pytest.raises(NonPositiveCpuError):
        PerfRatioSession.build(raw)


def test_another_schema_is_refused() -> None:
    raw = _raw(_session())
    raw["schema"] = "perf-ratio/0"
    with pytest.raises(MalformedRatioRecordError):
        PerfRatioSession.build(raw)


def test_the_stage_names_the_reference_short_rev() -> None:
    assert _session().stage == "ratio-vs-15fbed70"


def test_both_arms_are_required() -> None:
    raw = _raw(_session())
    del raw["arms"]["head"]
    with pytest.raises(MalformedRatioRecordError):
        PerfRatioSession.build(raw)


def test_same_output_is_digest_equality() -> None:
    same = _session(head_digest="digest-a")
    different = _session(head_digest="digest-b")
    assert same.workloads[0].same_output is True
    assert different.workloads[0].same_output is False
    assert PerfRatioSession.build(_raw(different)) == different


def test_the_reference_commit_is_the_ref_arm_rev() -> None:
    raw = _raw(_session())
    raw["ref_commit"] = HEAD_COMMIT
    with pytest.raises(InconsistentRatioRecordError):
        PerfRatioSession.build(raw)


def test_a_raw_workload_without_pairs_is_refused() -> None:
    with pytest.raises(MalformedRatioRecordError):
        WorkloadRatio(workload="w", family="gt", pairs=(), ref_digest="a", head_digest="a")


@pytest.mark.parametrize(
    ("key", "value"),
    [("host", 5), ("workloads", {}), ("arms", []), ("ref_commit", None)],
)
def test_an_accessor_type_error_surfaces_as_a_ratio_record_error(key: str, value: object) -> None:
    raw = _raw(_session())
    raw[key] = value
    with pytest.raises(MalformedRatioRecordError) as caught:
        PerfRatioSession.build(raw)
    assert not isinstance(caught.value, MalformedTimingRecordError)
    assert isinstance(caught.value.__cause__, MalformedTimingRecordError)


def test_load_session_round_trips_a_written_session(tmp_path: Path) -> None:
    session = _session()
    path = tmp_path / "perf_ratio_session.json"
    path.write_text(session.dumps(), encoding="utf-8")
    assert load_session(path) == session


@pytest.mark.parametrize("content", ["{not json", "[]", '{"schema": "perf-ratio/1"}'])
def test_load_session_refuses_a_malformed_file(tmp_path: Path, content: str) -> None:
    path = tmp_path / "perf_ratio_session.json"
    path.write_text(content, encoding="utf-8")
    with pytest.raises(MalformedRatioRecordError):
        load_session(path)


Mutation = Callable[[Any], None]


def _set(*path: Any, value: object) -> Mutation:
    def mutate(raw: Any) -> None:
        target = raw
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value

    return mutate


def _delete(*path: Any) -> Mutation:
    def mutate(raw: Any) -> None:
        target = raw
        for key in path[:-1]:
            target = target[key]
        del target[path[-1]]

    return mutate


def _duplicate_first_workload(raw: Any) -> None:
    raw["workloads"].append(json.loads(json.dumps(raw["workloads"][0])))


def _with_patch_digest(*, dirty: object) -> Mutation:
    # A patch digest is present, so the dirty-xor-patch check accepts a truthy dirty and only the boolean check refuses it.
    def mutate(raw: Any) -> None:
        raw["arms"]["head"]["dirty"] = dirty
        raw["arms"]["head"]["patch_sha256"] = "patch"

    return mutate


PAIR0 = ("workloads", 0, "pairs", 0)

MALFORMED_MUTATIONS = {
    "unexpected session key": _set("extra", value=1),
    "unexpected arm key": _set("arms", "ref", "extra", value="x"),
    "unexpected arm name": _set("arms", "base", value={}),
    "unexpected workload key": _set("workloads", 0, "pair_ratio_mean", value=1.0),
    "missing workload key": _delete("workloads", 0, "family"),
    "unexpected pair key": _set(*PAIR0, "wall", value=1.0),
    "unexpected injected key": _set("injected", value={"arm": "head", "seconds": 1.0, "extra": 1}),
    "pair order runs one arm twice": _set(*PAIR0, "order", value=["ref", "ref"]),
    "pair order names three arms": _set(*PAIR0, "order", value=["ref", "head", "ref"]),
    "pair order names an unknown arm": _set(*PAIR0, "order", value=["ref", "base"]),
    "infinite cpu": _set(*PAIR0, "ref_cpu", value=math.inf),
    "nan load": _set(*PAIR0, "load_after", value=math.nan),
    "decimal-string cpu": _set(*PAIR0, "head_cpu", value="3.0"),
    "decimal-string derived value": _set("workloads", 0, "ratio", value="2.0"),
    "boolean cpu": _set(*PAIR0, "ref_cpu", value=True),
    "empty pairs": _set("workloads", 0, "pairs", value=[]),
    "empty workloads": _set("workloads", value=[]),
    "duplicate workload names": _duplicate_first_workload,
    "dirty without patch digest": _set("arms", "head", "dirty", value=True),
    "clean with patch digest": _set("arms", "head", "patch_sha256", value="patch"),
    "non-boolean dirty": _with_patch_digest(dirty=1),
    "non-boolean same_output": _set("workloads", 0, "same_output", value=1),
    "start without utc offset": _set("started_utc", value="2026-11-01T10:00:00"),
    "end not iso 8601": _set("ended_utc", value="yesterday"),
    "zero injected seconds": _set("injected", value={"arm": "head", "seconds": 0.0}),
    "negative injected seconds": _set("injected", value={"arm": "ref", "seconds": -1.0}),
    "injected unknown arm": _set("injected", value={"arm": "base", "seconds": 1.0}),
}


@pytest.mark.parametrize("mutate", list(MALFORMED_MUTATIONS.values()), ids=list(MALFORMED_MUTATIONS))
def test_a_malformed_record_is_refused(mutate: Mutation) -> None:
    raw = _raw(_session())
    mutate(raw)
    with pytest.raises(MalformedRatioRecordError):
        PerfRatioSession.build(raw)


@pytest.mark.parametrize(
    ("order", "ref_cpu", "load"),
    [
        ((ArmName.HEAD, ArmName.HEAD), 1.0, 1.0),
        ((ArmName.REF, ArmName.HEAD), math.inf, 1.0),
        ((ArmName.REF, ArmName.HEAD), 1.0, math.inf),
    ],
    ids=["order runs one arm twice", "infinite cpu", "infinite load"],
)
def test_a_malformed_pair_is_refused_at_build(order: tuple[ArmName, ArmName], ref_cpu: float, load: float) -> None:
    with pytest.raises(MalformedRatioRecordError):
        PairRun.build(order=order, ref_cpu=CpuSeconds(ref_cpu), head_cpu=CpuSeconds(1.0), load_before=LoadOneMinute(load), load_after=LoadOneMinute(1.0))
