from __future__ import annotations

from datetime import datetime
from datetime import timedelta
from datetime import timezone

import pytest

import benchmarks.held_path_report as report_module
from benchmarks.errors import InvalidHeldPathReportContextError
from benchmarks.held_path_report import HeldPathReportContext
from benchmarks.held_path_report import render_held_figure5_2d_path
from benchmarks.quality_observations import CRITERION_NAMES
from tests.benchmarks.test_held_path_evidence import _build
from tests.benchmarks.test_held_path_evidence import _inputs
from tests.benchmarks.test_held_path_evidence import _sample
from tests.benchmarks.test_held_post_qualification import _characterization


UTC_INSTANT = datetime(2026, 9, 5, 8, 9, 10, tzinfo=timezone.utc)


def _context(*, command: str = "pixi run held-figure5-characterize", policy: str = "benchmarks.runner.generate_toolpath") -> HeldPathReportContext:
    return HeldPathReportContext.build(generated_at_utc=UTC_INSTANT, pixi_command=command, generator_policy_name=policy)


def test_context_has_one_strict_utc_factory() -> None:
    with pytest.raises(TypeError):
        HeldPathReportContext()

    context = _context()
    assert context.generated_at_utc is UTC_INSTANT
    assert context.pixi_command == "pixi run held-figure5-characterize"
    assert context.generator_policy_name == "benchmarks.runner.generate_toolpath"


@pytest.mark.parametrize(
    ("generated_at", "command", "policy"),
    [
        (datetime(2026, 9, 5), "run", "policy"),
        (datetime(2026, 9, 5, tzinfo=timezone(timedelta(hours=1))), "run", "policy"),
        ("2026-09-05", "run", "policy"),
        (UTC_INSTANT, "", "policy"),
        (UTC_INSTANT, "  ", "policy"),
        (UTC_INSTANT, 2, "policy"),
        (UTC_INSTANT, "run", ""),
        (UTC_INSTANT, "run", "\t"),
        (UTC_INSTANT, "run", None),
    ],
)
def test_context_rejects_invalid_metadata(generated_at: object, command: object, policy: object) -> None:
    with pytest.raises(InvalidHeldPathReportContextError):
        HeldPathReportContext.build(  # type: ignore[arg-type]
            generated_at_utc=generated_at,
            pixi_command=command,
            generator_policy_name=policy,
        )


def test_report_is_deterministic_complete_and_bounded() -> None:
    characterization = _characterization(all_open=True)
    first = render_held_figure5_2d_path(characterization, _context())
    second = render_held_figure5_2d_path(characterization, _context())

    assert first == second
    assert first.endswith("\n")
    assert first.startswith("# Held Figure 5 2D path characterization\n\n**HISTORICAL 2D BENCHMARK CHARACTERIZATION - NOT A MANUFACTURING RELEASE**")
    for required in (
        "2026-09-05T08:09:10+00:00",
        "Benchmark coordinates and lengths are normalized-scale millimetres; they are not a machine setup.",
        "Sampled negative observations are bounded evidence, not global proofs.",
        "The exercised guarded replay is not evidence of native audit compatibility.",
        "No controller, machine, setup, or manufacturing-release claim is made.",
        "figure5",
        "| Tool diameter (benchmark-normalized mm) | 2.000000 |",
        "| TEA cap (degrees) | 80.000000 |",
        "| Reference primitives | 31 |",
        "| Projection vertices | 65 |",
        "| Total source operations | 3 |",
        "| TEA-audited lateral operations | 2 |",
        "| TEA-audit-excluded operations | 1 |",
        "| Sampled material-contact operations |",
        "| Generation | 1.000000 |",
        "| Guarded audit | 2.000000 |",
        "| Survey | 3.000000 |",
        "| Coverage plus reduction | 4.000000 |",
        "| certified |",
        "| demonstrated_exceeded |",
        "| unresolved |",
        "## Report-only PathQuality",
        "not eligible for postprocessor qualification",
        "This verdict is not a manufacturing release.",
    ):
        assert required in first
    assert "G-code" not in first
    assert "cycle-time estimate" not in first
    assert "outperforms Held" not in first
    assert "native audit compatible" not in first


def test_witness_and_criterion_rows_preserve_canonical_order() -> None:
    samples = (
        _sample(cap=True, engagement=11.0, x=3.0, y=4.0),
        _sample(cap=True, engagement=12.0, x=1.0, y=2.0),
    )
    characterization = _build(_inputs(samples=samples, second_certified=False, cap_violations=1))
    markdown = render_held_figure5_2d_path(characterization, _context())

    first_witness = "| 2 | 3.000000 | 4.000000 |"
    second_witness = "| 2 | 1.000000 | 2.000000 |"
    assert markdown.index(first_witness) < markdown.index(second_witness)
    criterion_positions = [markdown.index(f"| {name} |") for name in CRITERION_NAMES]
    assert criterion_positions == sorted(criterion_positions)
    assert markdown.count("| no_failure_observed |") >= 1
    assert markdown.count("| criterion_satisfied |") >= 1
    assert markdown.count("| within_declared_tolerance |") >= 1


def test_report_only_rows_are_literal_and_gate_fields_are_not_duplicated() -> None:
    markdown = render_held_figure5_2d_path(_characterization(), _context())
    report_only_names = (
        "gouge_free",
        "rapid_safety",
        "marginal_loops",
        "recut_fraction",
        "max_engagement_deg",
        "engagement_p95_deg",
        "engagement_variance_deg2",
        "max_chip_thickness_ratio",
        "low_chip_thickness_ratio",
        "max_engagement_gradient_deg_per_length",
        "immersion_steady_fraction",
        "immersion_at_design_fraction",
        "immersion_excursions",
        "mean_radial_depth",
        "radial_depth_variance",
        "wall_scallop_height",
        "loop_radius_cv",
        "cutting_length",
        "air_length",
        "air_fraction",
        "max_curvature",
        "curvature_breaks",
        "direction_reversals",
        "retract_count",
        "reentry_count",
        "material_entries",
        "cut_air_alternations",
        "alternations_per_length",
        "block_count",
        "cut_blocks",
        "min_block_length",
        "median_block_length",
        "short_block_length",
        "block_length_cv",
        "blocks_per_unit_length",
        "arc_length_fraction",
        "cut_operations",
        "path_length",
    )
    for name in report_only_names:
        assert markdown.count(f"| {name} |") == 1
    assert "engagement_length_histogram[0]" in markdown
    for name in (
        "uncut_fraction",
        "gouging_motions",
        "unsafe_rapids",
        "continuity_breaks",
        "zero_length_motions",
        "degenerate_loops",
        "redundant_operations",
        "cap_exceedances",
        "slotting_motions",
        "max_engagement_step_deg",
        "max_loop_radius_step",
        "tangent_breaks",
    ):
        assert f"| {name} |" not in markdown


def test_renderer_uses_canonical_failure_collector_once(monkeypatch: pytest.MonkeyPatch) -> None:
    characterization = _characterization()
    calls: list[object] = []

    def failures(value: object) -> tuple[str, ...]:
        calls.append(value)
        return ("first failure", "second failure")

    monkeypatch.setattr(report_module, "post_qualification_failures", failures)
    markdown = render_held_figure5_2d_path(characterization, _context())

    assert calls == [characterization]
    assert markdown.index("first failure") < markdown.index("second failure")
    assert "not eligible for postprocessor qualification" in markdown

    monkeypatch.setattr(report_module, "post_qualification_failures", lambda value: ())
    eligible = render_held_figure5_2d_path(characterization, _context())
    assert "eligible for postprocessor qualification evidence entry" in eligible
    assert "This verdict is not a manufacturing release." in eligible


def test_dynamic_cells_escape_backslash_pipe_and_newlines(monkeypatch: pytest.MonkeyPatch) -> None:
    characterization = _characterization()
    monkeypatch.setattr(report_module, "post_qualification_failures", lambda value: ("bad\\pipe|one\r\ntwo\rthree\nfour",))
    context = _context(command="run\\this|now\r\nnext", policy="policy|name")

    markdown = render_held_figure5_2d_path(characterization, context)

    assert "run\\\\this\\|now<br>next" in markdown
    assert "policy\\|name" in markdown
    assert "bad\\\\pipe\\|one<br>two<br>three<br>four" in markdown
