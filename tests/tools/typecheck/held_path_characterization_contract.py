from collections.abc import Callable
from datetime import datetime
from datetime import timezone
from pathlib import Path

from typing_extensions import assert_type

from benchmarks.held_path_evidence import HeldFigure5Characterization
from benchmarks.held_path_report import HeldPathReportContext
from benchmarks.held_path_report import render_held_figure5_2d_path
from benchmarks.held_post_qualification import HeldPostQualificationCandidate
from tools.held_path_characterization import DEFAULT_REPORT_PATH
from tools.held_path_characterization import _utc_now
from tools.held_path_characterization import write_held_figure5_report


assert_type(DEFAULT_REPORT_PATH, Path)
assert_type(_utc_now(), datetime)
writer: Callable[..., HeldFigure5Characterization] = write_held_figure5_report
assert_type(writer, Callable[..., HeldFigure5Characterization])
characterization = assert_type(
    write_held_figure5_report(pixi_command="pixi run held-figure5-characterize", path=Path("report.md")),
    HeldFigure5Characterization,
)
candidate: HeldPostQualificationCandidate = characterization  # type: ignore[assignment]
missing: None = characterization  # type: ignore[assignment]
context = assert_type(
    HeldPathReportContext.build(
        generated_at_utc=datetime(2026, 9, 5, tzinfo=timezone.utc),
        pixi_command="pixi run held-figure5-characterize",
        generator_policy_name="benchmarks.runner.generate_toolpath",
    ),
    HeldPathReportContext,
)
assert_type(render_held_figure5_2d_path(characterization, context), str)
assert_type(candidate, HeldPostQualificationCandidate)
assert_type(missing, None)
