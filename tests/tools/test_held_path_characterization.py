from __future__ import annotations

from datetime import datetime
from datetime import timezone
from pathlib import Path
from typing import cast

import pytest

import tools.held_path_characterization as tool_module
from benchmarks.held_path_characterize import CharacterizationPhase
from tests.benchmarks.test_held_post_qualification import _characterization


UTC_INSTANT = datetime(2026, 9, 5, 8, 9, 10, tzinfo=timezone.utc)


class PendingSpyPath:
    def __init__(self, destination: "SpyPath") -> None:
        self.destination = destination

    def write_text(self, markdown: str, *, encoding: str) -> int:
        self.destination.calls.append(("pending-write", markdown, encoding))
        self.destination.pending_content = markdown
        if self.destination.failure_stage == "write":
            raise self.destination.failure
        return len(markdown)

    def replace(self, destination: object) -> object:
        assert destination is self.destination
        self.destination.calls.append(("replace",))
        if self.destination.failure_stage == "replace":
            raise self.destination.failure
        self.destination.content = self.destination.pending_content
        self.destination.pending_content = None
        return destination


class SpyPath:
    def __init__(self) -> None:
        self.content = "prior report\n"
        self.pending_content: str | None = None
        self.calls: list[tuple[object, ...]] = []
        self.failure_stage: str | None = None
        self.failure = RuntimeError("publication failure")

    def with_suffix(self, suffix: str) -> PendingSpyPath:
        assert suffix == ".pending.md"
        self.calls.append(("pending-path", suffix))
        return PendingSpyPath(self)


def test_writer_requires_an_explicit_phase_observer() -> None:
    with pytest.raises(TypeError):
        tool_module.write_held_figure5_report(pixi_command="run")  # type: ignore[call-arg]


def test_writer_wires_production_dependencies_observer_context_and_one_write(monkeypatch: pytest.MonkeyPatch) -> None:
    characterization = _characterization(all_open=True)
    destination = SpyPath()
    characterize_calls: list[tuple[object, ...]] = []
    rendered: list[tuple[object, object]] = []
    phases: list[CharacterizationPhase] = []
    monkeypatch.setattr(tool_module, "_print_phase", phases.append)

    def characterize(*args: object, phase_observer: object) -> object:
        characterize_calls.append((*args, phase_observer))
        for phase in cast(tuple[CharacterizationPhase, ...], ("generation", "guarded_audit", "survey", "quality_reduction")):
            cast(object, phase_observer)(phase)  # type: ignore[operator]
        return characterization

    def render(value: object, context: object) -> str:
        rendered.append((value, context))
        return "report\n"

    monkeypatch.setattr(tool_module, "characterize_figure5", characterize)
    monkeypatch.setattr(tool_module, "render_held_figure5_2d_path", render)
    monkeypatch.setattr(tool_module, "_utc_now", lambda: UTC_INSTANT)

    observer = phases.append
    returned = tool_module.write_held_figure5_report(
        pixi_command="pixi run held-figure5-characterize",
        phase_observer=observer,
        path=cast(Path, destination),
    )

    assert returned is characterization
    assert characterize_calls == [
        (
            tool_module.generate_toolpath,
            tool_module.audit_figure5_engagement,
            tool_module.survey_path,
            tool_module.reduce_figure5_quality,
            observer,
        )
    ]
    assert phases == ["generation", "guarded_audit", "survey", "quality_reduction"]
    assert rendered[0][0] is characterization
    context = rendered[0][1]
    assert context.generated_at_utc is UTC_INSTANT
    assert context.pixi_command == "pixi run held-figure5-characterize"
    assert context.generator_policy_name == "benchmarks.runner.generate_toolpath"
    assert destination.calls == [
        ("pending-path", ".pending.md"),
        ("pending-write", "report\n", "utf-8"),
        ("replace",),
    ]
    assert destination.content == "report\n"
    assert destination.pending_content is None


@pytest.mark.parametrize("failure_stage", ["characterize", "clock", "render", "write", "replace"])
def test_boundary_failures_propagate_without_later_side_effects(monkeypatch: pytest.MonkeyPatch, failure_stage: str) -> None:
    characterization = _characterization()
    destination = SpyPath()
    calls: list[str] = []
    failure = RuntimeError(failure_stage)

    def characterize(*args: object, **kwargs: object) -> object:
        calls.append("characterize")
        if failure_stage == "characterize":
            raise failure
        return characterization

    def clock() -> datetime:
        calls.append("clock")
        if failure_stage == "clock":
            raise failure
        return UTC_INSTANT

    def render(value: object, context: object) -> str:
        calls.append("render")
        if failure_stage == "render":
            raise failure
        return "report\n"

    destination.failure_stage = failure_stage
    destination.failure = failure

    monkeypatch.setattr(tool_module, "characterize_figure5", characterize)
    monkeypatch.setattr(tool_module, "_utc_now", clock)
    monkeypatch.setattr(tool_module, "render_held_figure5_2d_path", render)

    with pytest.raises(RuntimeError) as exc_info:
        tool_module.write_held_figure5_report(
            pixi_command="run",
            phase_observer=lambda phase: None,
            path=cast(Path, destination),
        )
    assert exc_info.value is failure
    if failure_stage == "characterize":
        assert calls == ["characterize"]
    elif failure_stage == "clock":
        assert calls == ["characterize", "clock"]
    elif failure_stage == "render":
        assert calls == ["characterize", "clock", "render"]
    else:
        assert calls == ["characterize", "clock", "render"]
        assert destination.content == "prior report\n"
        assert destination.pending_content == "report\n"


def test_print_phase_has_closed_stable_spelling(capsys: pytest.CaptureFixture[str]) -> None:
    for phase in cast(tuple[CharacterizationPhase, ...], ("generation", "guarded_audit", "survey", "quality_reduction")):
        tool_module._print_phase(phase)
    assert capsys.readouterr().out.splitlines() == [
        "PHASE: generation",
        "PHASE: guarded_audit",
        "PHASE: survey",
        "PHASE: quality_reduction",
    ]


def test_main_reports_completion_path_and_canonical_ineligible_verdict(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    characterization = _characterization(all_open=True)
    writes: list[tuple[str, object, Path]] = []
    failure_calls: list[object] = []

    def writer(*, pixi_command: str, phase_observer: object, path: Path) -> object:
        writes.append((pixi_command, phase_observer, path))
        return characterization

    def failures(value: object) -> tuple[str, ...]:
        failure_calls.append(value)
        return ("sentinel refusal",)

    monkeypatch.setattr(tool_module, "write_held_figure5_report", writer)
    monkeypatch.setattr(tool_module, "post_qualification_failures", failures)
    tool_module.main()

    assert writes == [("pixi run held-figure5-characterize", tool_module._print_phase, tool_module.DEFAULT_REPORT_PATH)]
    assert failure_calls == [characterization]
    output = capsys.readouterr().out
    assert str(tool_module.DEFAULT_REPORT_PATH) in output
    assert "CHARACTERIZATION COMPLETED" in output
    assert "not eligible for postprocessor qualification: sentinel refusal" in output
    assert "Exit zero means evidence completion, not path eligibility." in output


def test_main_reports_canonical_eligible_verdict(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    characterization = _characterization()
    monkeypatch.setattr(tool_module, "write_held_figure5_report", lambda **kwargs: characterization)
    monkeypatch.setattr(tool_module, "post_qualification_failures", lambda value: ())
    tool_module.main()
    assert "eligible for postprocessor qualification evidence entry" in capsys.readouterr().out
