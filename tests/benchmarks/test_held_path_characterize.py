from __future__ import annotations

import importlib
from dataclasses import replace
from typing import cast

import numpy as np
import pytest
from compas.geometry import Line
from compas.geometry import Point
from compas.geometry import Polygon

from benchmarks.errors import MutatedHeldToolpathError
from benchmarks.held_path_evidence import HeldFigure5Characterization
from benchmarks.held_path_snapshot import HeldOperationSnapshot
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.quality_observations import QualityEvidence
from benchmarks.spec import PocketSpec
from benchmarks.survey import PathSurvey
from benchmarks.units import Seconds
from compas_cgal.engagement import EngagementReport
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult


def _result() -> ToolpathResult:
    operation = ToolpathOperation(
        geometry=Line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
        operation=OperationType.CUT,
        path_index=0,
        clockwise=False,
        start_tangent=np.array([1.0, 0.0, 0.0]),
        end_tangent=np.array([1.0, 0.0, 0.0]),
    )
    return ToolpathResult(
        operations=[operation],
        polyline=np.empty((0, 3), dtype=np.float64),
    )


def _sentinel_report() -> EngagementReport:
    return EngagementReport(2.0, 1.0, [], 0.0, 0, 0)


def test_engagement_adapter_translates_exact_values_and_identity(monkeypatch: pytest.MonkeyPatch) -> None:
    adapters = importlib.import_module("benchmarks.held_consumer_adapters")
    base_spec = load_held_reference_case("figure5").pocket_spec()
    hole = Polygon([Point(1.0, 1.0), Point(2.0, 1.0), Point(1.0, 2.0)])
    spec = replace(base_spec, holes=(hole,))
    result = _result()
    report = _sentinel_report()
    received: list[tuple[object, ...]] = []

    def audit(*args: object) -> EngagementReport:
        received.append(args)
        return report

    monkeypatch.setattr(adapters, "audit_toolpath_engagement", audit)

    actual = adapters.audit_figure5_engagement(spec, result)

    assert actual is report
    assert len(received) == 1
    polygon, received_result, diameter, cap, holes = received[0]
    assert polygon is spec.polygon
    assert received_result is result
    assert diameter == spec.tool_diameter
    assert cap == spec.tea_cap_rad
    assert type(holes) is list
    assert holes is not spec.holes
    assert len(cast(list[object], holes)) == len(spec.holes)
    assert all(actual_hole is expected_hole for actual_hole, expected_hole in zip(cast(list[object], holes), spec.holes))


def test_engagement_adapter_propagates_exception_object(monkeypatch: pytest.MonkeyPatch) -> None:
    adapters = importlib.import_module("benchmarks.held_consumer_adapters")
    failure = RuntimeError("sentinel audit failure")

    def audit(*args: object) -> EngagementReport:
        raise failure

    monkeypatch.setattr(adapters, "audit_toolpath_engagement", audit)

    with pytest.raises(RuntimeError) as caught:
        adapters.audit_figure5_engagement(load_held_reference_case("figure5").pocket_spec(), _result())

    assert caught.value is failure


class _CaseSpy:
    def __init__(self, spec: PocketSpec, events: list[str]) -> None:
        self._spec = spec
        self._events = events

    def pocket_spec(self) -> PocketSpec:
        self._events.append("pocket_spec")
        return self._spec


def test_characterization_has_exact_order_identities_and_durations(monkeypatch: pytest.MonkeyPatch) -> None:
    characterize = importlib.import_module("benchmarks.held_path_characterize")
    events: list[str] = []
    spec = load_held_reference_case("figure5").pocket_spec()
    case = cast(HeldReferenceCase, _CaseSpy(spec, events))
    result = _result()
    snapshot = cast(tuple[HeldOperationSnapshot, ...], (object(),))
    report = _sentinel_report()
    survey = cast(PathSurvey, object())
    quality = cast(QualityEvidence, object())
    characterization = cast(HeldFigure5Characterization, object())
    times = iter((Seconds(1.0), Seconds(3.0), Seconds(10.0), Seconds(14.0), Seconds(20.0), Seconds(25.0), Seconds(30.0), Seconds(37.0)))
    generator_calls = 0

    def load(name: str) -> HeldReferenceCase:
        assert name == "figure5"
        events.append("load")
        return case

    def phase_observer(phase: object) -> None:
        events.append(f"phase:{phase}")

    def clock() -> Seconds:
        events.append("clock")
        return next(times)

    def generator(received_spec: PocketSpec) -> ToolpathResult:
        nonlocal generator_calls
        generator_calls += 1
        assert received_spec is spec
        events.append("generate")
        return result

    def take_snapshot(received_result: ToolpathResult) -> tuple[HeldOperationSnapshot, ...]:
        assert received_result is result
        events.append("snapshot")
        return snapshot

    def assert_unchanged(received_result: ToolpathResult, received_snapshot: tuple[HeldOperationSnapshot, ...]) -> None:
        assert received_result is result
        assert received_snapshot is snapshot
        events.append("assert")

    def audit(received_spec: PocketSpec, received_result: ToolpathResult) -> EngagementReport:
        assert received_spec is spec
        assert received_result is result
        events.append("audit")
        return report

    def survey_path(received_spec: PocketSpec, received_result: ToolpathResult) -> PathSurvey:
        assert received_spec is spec
        assert received_result is result
        events.append("survey")
        return survey

    def reduce(
        received_spec: PocketSpec,
        received_snapshot: tuple[HeldOperationSnapshot, ...],
        received_survey: PathSurvey,
    ) -> QualityEvidence:
        assert received_spec is spec
        assert received_snapshot is snapshot
        assert received_survey is survey
        events.append("reduce")
        return quality

    class BuildBoundary:
        @staticmethod
        def build(**values: object) -> HeldFigure5Characterization:
            events.append("build")
            assert values == {
                "case": case,
                "snapshot": snapshot,
                "audit": report,
                "survey": survey,
                "quality": quality,
                "generation_seconds": Seconds(2.0),
                "audit_seconds": Seconds(4.0),
                "survey_seconds": Seconds(5.0),
                "reduction_seconds": Seconds(7.0),
            }
            assert values["case"] is case
            assert values["snapshot"] is snapshot
            assert values["audit"] is report
            assert values["survey"] is survey
            assert values["quality"] is quality
            return characterization

    monkeypatch.setattr(characterize, "load_held_reference_case", load)
    monkeypatch.setattr(characterize, "snapshot_toolpath", take_snapshot)
    monkeypatch.setattr(characterize, "assert_toolpath_matches_snapshot", assert_unchanged)
    monkeypatch.setattr(characterize, "HeldFigure5Characterization", BuildBoundary)

    actual = characterize.characterize_figure5(
        generator,
        audit,
        survey_path,
        reduce,
        phase_observer=phase_observer,
        clock=clock,
    )

    assert actual is characterization
    assert generator_calls == 1
    assert events == [
        "load",
        "pocket_spec",
        "phase:generation",
        "clock",
        "generate",
        "clock",
        "snapshot",
        "assert",
        "phase:guarded_audit",
        "clock",
        "assert",
        "audit",
        "clock",
        "assert",
        "phase:survey",
        "clock",
        "assert",
        "survey",
        "clock",
        "assert",
        "phase:quality_reduction",
        "clock",
        "assert",
        "reduce",
        "clock",
        "assert",
        "build",
    ]


def _constant_clock() -> Seconds:
    return Seconds(0.0)


@pytest.mark.parametrize("mutating_stage", ["audit", "survey", "reduce"])
def test_characterization_rejects_consumer_mutation(
    monkeypatch: pytest.MonkeyPatch,
    mutating_stage: str,
) -> None:
    characterize = importlib.import_module("benchmarks.held_path_characterize")
    result = _result()
    report = _sentinel_report()
    survey = cast(PathSurvey, object())
    quality = cast(QualityEvidence, object())
    phases: list[object] = []
    builds: list[object] = []

    def mutate() -> None:
        result.operations[0].operation = OperationType.RETRACT

    def audit(spec: PocketSpec, received_result: ToolpathResult) -> EngagementReport:
        if mutating_stage == "audit":
            mutate()
        return report

    def survey_path(spec: PocketSpec, received_result: ToolpathResult) -> PathSurvey:
        if mutating_stage == "survey":
            mutate()
        return survey

    def reduce(spec: PocketSpec, snapshot: tuple[HeldOperationSnapshot, ...], received_survey: PathSurvey) -> QualityEvidence:
        if mutating_stage == "reduce":
            mutate()
        return quality

    class BuildBoundary:
        @staticmethod
        def build(**values: object) -> HeldFigure5Characterization:
            builds.append(values)
            return cast(HeldFigure5Characterization, object())

    monkeypatch.setattr(characterize, "HeldFigure5Characterization", BuildBoundary)

    with pytest.raises(MutatedHeldToolpathError):
        characterize.characterize_figure5(
            lambda spec: result,
            audit,
            survey_path,
            reduce,
            phase_observer=phases.append,
            clock=_constant_clock,
        )

    expected_phases = {
        "audit": ["generation", "guarded_audit"],
        "survey": ["generation", "guarded_audit", "survey"],
        "reduce": ["generation", "guarded_audit", "survey", "quality_reduction"],
    }
    assert phases == expected_phases[mutating_stage]
    assert builds == []


@pytest.mark.parametrize("mutator", ["observer", "clock"])
@pytest.mark.parametrize("affected_phase", ["guarded_audit", "survey", "quality_reduction"])
def test_characterization_rejects_callback_mutation_before_affected_consumer(
    monkeypatch: pytest.MonkeyPatch,
    mutator: str,
    affected_phase: str,
) -> None:
    characterize = importlib.import_module("benchmarks.held_path_characterize")
    result = _result()
    report = _sentinel_report()
    survey = cast(PathSurvey, object())
    quality = cast(QualityEvidence, object())
    current_phase = ""
    start_clock_pending = False
    consumer_calls: list[str] = []
    builds: list[object] = []

    def mutate() -> None:
        result.operations[0].operation = OperationType.RETRACT

    def phase_observer(phase: str) -> None:
        nonlocal current_phase, start_clock_pending
        current_phase = phase
        start_clock_pending = True
        if mutator == "observer" and phase == affected_phase:
            mutate()

    def clock() -> Seconds:
        nonlocal start_clock_pending
        if start_clock_pending:
            start_clock_pending = False
            if mutator == "clock" and current_phase == affected_phase:
                mutate()
        return Seconds(0.0)

    def audit(spec: PocketSpec, received_result: ToolpathResult) -> EngagementReport:
        consumer_calls.append("guarded_audit")
        return report

    def survey_path(spec: PocketSpec, received_result: ToolpathResult) -> PathSurvey:
        consumer_calls.append("survey")
        return survey

    def reduce(spec: PocketSpec, snapshot: tuple[HeldOperationSnapshot, ...], received_survey: PathSurvey) -> QualityEvidence:
        consumer_calls.append("quality_reduction")
        return quality

    class BuildBoundary:
        @staticmethod
        def build(**values: object) -> HeldFigure5Characterization:
            builds.append(values)
            return cast(HeldFigure5Characterization, object())

    monkeypatch.setattr(characterize, "HeldFigure5Characterization", BuildBoundary)

    with pytest.raises(MutatedHeldToolpathError):
        characterize.characterize_figure5(
            lambda spec: result,
            audit,
            survey_path,
            reduce,
            phase_observer=phase_observer,
            clock=clock,
        )

    ordered_consumers = ["guarded_audit", "survey", "quality_reduction"]
    assert consumer_calls == ordered_consumers[: ordered_consumers.index(affected_phase)]
    assert builds == []


@pytest.mark.parametrize("failing_stage", ["generation", "audit", "survey", "reduce"])
def test_characterization_propagates_stage_exception_without_partial_build(
    monkeypatch: pytest.MonkeyPatch,
    failing_stage: str,
) -> None:
    characterize = importlib.import_module("benchmarks.held_path_characterize")
    failure = RuntimeError(f"sentinel {failing_stage} failure")
    result = _result()
    report = _sentinel_report()
    survey = cast(PathSurvey, object())
    quality = cast(QualityEvidence, object())
    calls: list[str] = []
    builds: list[object] = []

    def stage(name: str, returned: object) -> object:
        calls.append(name)
        if failing_stage == name:
            raise failure
        return returned

    def generator(spec: PocketSpec) -> ToolpathResult:
        return cast(ToolpathResult, stage("generation", result))

    def audit(spec: PocketSpec, received_result: ToolpathResult) -> EngagementReport:
        return cast(EngagementReport, stage("audit", report))

    def survey_path(spec: PocketSpec, received_result: ToolpathResult) -> PathSurvey:
        return cast(PathSurvey, stage("survey", survey))

    def reduce(spec: PocketSpec, snapshot: tuple[HeldOperationSnapshot, ...], received_survey: PathSurvey) -> QualityEvidence:
        return cast(QualityEvidence, stage("reduce", quality))

    class BuildBoundary:
        @staticmethod
        def build(**values: object) -> HeldFigure5Characterization:
            builds.append(values)
            return cast(HeldFigure5Characterization, object())

    monkeypatch.setattr(characterize, "HeldFigure5Characterization", BuildBoundary)

    with pytest.raises(RuntimeError) as caught:
        characterize.characterize_figure5(
            generator,
            audit,
            survey_path,
            reduce,
            phase_observer=lambda phase: None,
            clock=_constant_clock,
        )

    stage_order = ["generation", "audit", "survey", "reduce"]
    assert caught.value is failure
    assert calls == stage_order[: stage_order.index(failing_stage) + 1]
    assert builds == []
