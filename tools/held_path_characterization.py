"""Write the deterministic Held Figure 5 evidence report."""

from __future__ import annotations

from collections.abc import Callable
from datetime import datetime
from datetime import timezone
from pathlib import Path

from typing_extensions import TypeAlias

from benchmarks.held_consumer_adapters import audit_figure5_engagement
from benchmarks.held_consumer_adapters import reduce_figure5_quality
from benchmarks.held_path_characterize import CharacterizationPhase
from benchmarks.held_path_characterize import characterize_figure5
from benchmarks.held_path_evidence import HeldFigure5Characterization
from benchmarks.held_path_report import HeldPathReportContext
from benchmarks.held_path_report import render_held_figure5_2d_path
from benchmarks.held_post_qualification import post_qualification_failures
from benchmarks.runner import generate_toolpath
from benchmarks.survey import survey_path

DEFAULT_REPORT_PATH = Path("docs/benchmarks/held_figure5_2d_path.md")
GENERATOR_POLICY_NAME = "benchmarks.runner.generate_toolpath"
UtcClock: TypeAlias = Callable[[], datetime]


def _utc_now() -> datetime:
    return datetime.now(timezone.utc)


def _print_phase(phase: CharacterizationPhase) -> None:
    print(f"PHASE: {phase}")


def write_held_figure5_report(
    *,
    pixi_command: str,
    phase_observer: Callable[[CharacterizationPhase], None],
    path: Path = DEFAULT_REPORT_PATH,
) -> HeldFigure5Characterization:
    """Characterize once and write exactly one report."""
    characterization = characterize_figure5(
        generate_toolpath,
        audit_figure5_engagement,
        survey_path,
        reduce_figure5_quality,
        phase_observer=phase_observer,
    )
    context = HeldPathReportContext.build(
        generated_at_utc=_utc_now(),
        pixi_command=pixi_command,
        generator_policy_name=GENERATOR_POLICY_NAME,
    )
    markdown = render_held_figure5_2d_path(characterization, context)
    pending_path = path.with_suffix(".pending.md")
    pending_path.write_text(markdown, encoding="utf-8")
    pending_path.replace(path)
    return characterization


def main() -> None:
    characterization = write_held_figure5_report(
        pixi_command="pixi run held-figure5-characterize",
        phase_observer=_print_phase,
        path=DEFAULT_REPORT_PATH,
    )
    failures = post_qualification_failures(characterization)
    print(f"REPORT: {DEFAULT_REPORT_PATH}")
    print("CHARACTERIZATION COMPLETED")
    if failures:
        print("not eligible for postprocessor qualification: " + "; ".join(failures))
    else:
        print("eligible for postprocessor qualification evidence entry")
    print("Exit zero means evidence completion, not path eligibility.")


if __name__ == "__main__":
    main()
