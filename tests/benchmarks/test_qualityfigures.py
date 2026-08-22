"""Contract tests for the machining-quality figures.

No golden images: a figure's bytes change when matplotlib changes, and a test
that fails on a library upgrade teaches nobody anything. What is pinned instead
is (1) that every writer produces a real file with LIVE TEXT, and (2) that the
DATA behind each figure is what its caption claims -- because a caption is a
claim, and a figure whose claim has quietly stopped being true is worse than no
figure at all.

Generating and surveying a tool path costs seconds, so the two paths every
figure reads are built once for the whole module.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import matplotlib
import pytest

# The figure modules draw, so they must draw without a display here too.
matplotlib.use("Agg")

from benchmarks.coverage import measure_coverage  # noqa: E402
from benchmarks.figures import regenerate_every_figure  # noqa: E402
from benchmarks.palette import Theme  # noqa: E402
from benchmarks.pathmetrics import entry_cut_indices  # noqa: E402
from benchmarks.quality import chip_thickness_ratio  # noqa: E402
from benchmarks.quality import measure_quality  # noqa: E402
from benchmarks.qualityfigures import FIGURE_SMALL_POCKET  # noqa: E402
from benchmarks.qualityfigures import QUALITY_FIGURES  # noqa: E402
from benchmarks.qualityfigures import _cumulative_samples  # noqa: E402
from benchmarks.qualityfigures import corner_motions  # noqa: E402
from benchmarks.qualityfigures import draw_chip_thinning  # noqa: E402
from benchmarks.qualityfigures import draw_corner_defect  # noqa: E402
from benchmarks.qualityfigures import draw_coverage_residual  # noqa: E402
from benchmarks.qualityfigures import draw_curvature_feed  # noqa: E402
from benchmarks.qualityfigures import draw_engagement_along_path  # noqa: E402
from benchmarks.qualityfigures import draw_engagement_histogram  # noqa: E402
from benchmarks.qualityfigures import draw_engagement_map  # noqa: E402
from benchmarks.qualityfigures import draw_length_vs_time  # noqa: E402
from benchmarks.qualityfigures import draw_material_entries  # noqa: E402
from benchmarks.qualityfigures import measured_path  # noqa: E402
from benchmarks.qualityfigures import operating_band  # noqa: E402
from benchmarks.qualityfigures import per_operation_engagement  # noqa: E402
from benchmarks.qualityfigures import residual_points  # noqa: E402
from benchmarks.qualityfigures import write_chip_thinning  # noqa: E402
from benchmarks.qualityfigures import write_engagement_histogram  # noqa: E402

# One degree either side of the plateau corner, for the monotonicity sweep.
PLATEAU_PROBE_DEG = (90.0, 91.0, 120.0, 150.0, 179.0, 180.0)


@pytest.fixture(scope="module")
def measured() -> Any:
    """The 20x12 path every measured figure is drawn on, generated once."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return measured_path()


@pytest.fixture(scope="module")
def small() -> Any:
    """The 12x8 path, for the figures that compare two instances."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return measured_path(FIGURE_SMALL_POCKET)


# ---------------------------------------------------------------------------
# The claims each caption makes.
# ---------------------------------------------------------------------------


def test_the_corner_carries_a_degenerate_loop_a_link_and_a_proper_trochoid(measured: Any) -> None:
    """The figure annotates three motions; if the generator stops emitting them it must fail loudly."""
    loop, link, trochoid = corner_motions(measured.survey, measured.spec)
    assert loop.loop_radius is not None and loop.loop_radius <= measured.spec.tool_radius
    assert link.loop_radius is None
    assert trochoid.loop_radius is not None and trochoid.loop_radius > measured.spec.tool_radius
    assert link.index == loop.index + 1, "the slot the caption describes is the motion immediately after the plunge"


def test_the_residual_shading_counts_exactly_what_the_metric_counts(small: Any) -> None:
    """The shaded squares and the percentage beside them must be one measurement."""
    estimate = measure_coverage(small.spec, small.survey.final_stock)
    xs, ys = residual_points(small.spec, small.survey.final_stock)
    assert len(xs) == len(ys) == estimate.uncut_reachable_samples


def test_the_engagement_trace_runs_the_whole_cut(measured: Any) -> None:
    distances, engagements = _cumulative_samples(measured.survey)
    assert len(distances) == len(engagements)
    assert distances == sorted(distances)
    assert distances[-1] <= measured.survey.cut_length


def test_the_chip_curve_never_falls_and_is_flat_past_the_plateau() -> None:
    """The physics trap: `min(sin θ, 1)` would fall past 90° and invert the argument."""
    values = [chip_thickness_ratio(float(degrees)) for degrees in range(0, 181)]
    assert all(later >= earlier - 1e-12 for earlier, later in zip(values, values[1:]))
    for degrees in PLATEAU_PROBE_DEG:
        assert chip_thickness_ratio(degrees) == pytest.approx(1.0)


def test_the_operating_band_is_inside_the_engagement_the_path_actually_reaches(small: Any) -> None:
    low, high = operating_band(small.survey)
    peaks = [motion.peak_engagement_deg for motion in small.survey.motions if motion.is_engaged]
    assert 0.0 < low <= high <= max(peaks)


def test_every_operation_gets_an_engagement_or_an_explicit_none(measured: Any) -> None:
    """`ColourBy.ENGAGEMENT` needs exactly one entry per operation, or it refuses to draw."""
    peaks = per_operation_engagement(measured)
    assert len(peaks) == len(measured.result.operations)
    measured_indices = {motion.index for motion in measured.survey.motions}
    assert {index for index, value in enumerate(peaks) if value is not None} == measured_indices


def test_the_worst_curvature_is_the_tightest_loop(measured: Any) -> None:
    radii = [motion.loop_radius for motion in measured.survey.motions if motion.loop_radius]
    assert max(motion.curvature for motion in measured.survey.motions) == pytest.approx(1.0 / min(radii))


def test_the_feed_limited_time_is_never_shorter_than_the_ideal(measured: Any) -> None:
    """The curvature bound can only ever slow the machine down."""
    from benchmarks.models import MachineModel
    from benchmarks.quality import machine_outcome

    model = MachineModel.build()
    outcome = machine_outcome(measured.survey, model)
    assert outcome.cutting_seconds >= measured.survey.cut_length / model.feed_mm_per_s - 1e-9


def test_the_histogram_bands_account_for_the_whole_cut(measured: Any) -> None:
    quality = measure_quality(measured.spec, measured.result)
    banded = sum(length for _low, _high, length in quality.longevity.engagement_length_histogram)
    assert banded == pytest.approx(measured.survey.cut_length)


def test_every_marked_entry_is_an_entry_the_metric_counts(measured: Any) -> None:
    entries = entry_cut_indices(measured.result)
    marked = [motion for motion in measured.survey.motions if motion.index in entries]
    assert len(marked) == measure_quality(measured.spec, measured.result).longevity.material_entries


# ---------------------------------------------------------------------------
# Every figure draws, and draws differently in each theme.
# ---------------------------------------------------------------------------


def test_every_figure_draws_on_the_dark_surface(measured: Any, small: Any) -> None:
    """A crash in any of the nine is caught here.

    DARK only: the light surface is rendered end to end by the composite
    regeneration test below, and rendering both here doubled the module's run
    time to no extra effect.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for theme in (Theme.DARK,):
            assert draw_corner_defect(measured, theme=theme).figure is not None
            assert draw_engagement_along_path(measured, theme=theme).figure is not None
            assert draw_chip_thinning(small, theme=theme).figure is not None
            assert draw_engagement_map(measured, theme=theme).figure is not None
            assert draw_curvature_feed(measured, theme=theme).figure is not None
            assert draw_length_vs_time([small], theme=theme).figure is not None
            assert draw_engagement_histogram(measured, theme=theme).figure is not None
            assert draw_material_entries(measured, theme=theme).figure is not None
            assert draw_coverage_residual(small, theme=theme).figure is not None


def test_a_written_svg_carries_live_text_and_no_raster(tmp_path: Path) -> None:
    """`svg.fonttype: none` keeps real <text>; the default converts glyphs to outlines."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        written = write_chip_thinning(tmp_path, formats=("svg",))
    body = written[0].read_bytes()
    assert written[0].exists() and len(body) > 0
    assert b"<text" in body
    assert b"<image" not in body


def test_the_dark_variant_lands_on_its_own_file(tmp_path: Path) -> None:
    """The docs link the light figure by path; a dark redraw must not overwrite it."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        light = write_engagement_histogram(tmp_path, formats=("svg",), theme=Theme.LIGHT)
        dark = write_engagement_histogram(tmp_path, formats=("svg",), theme=Theme.DARK)
    assert light[0] != dark[0]
    assert light[0].read_bytes() != dark[0].read_bytes()


def test_the_quality_registry_names_every_figure_exactly_once() -> None:
    """A figure absent from the registry cannot be regenerated at all."""
    names = [name for name, _writer in QUALITY_FIGURES]
    assert len(names) == len(set(names)) == 9


def test_the_composite_entry_point_writes_both_families(tmp_path: Path) -> None:
    """`regenerate_every_figure` is the one command, so it must reach BOTH registries.

    Run on one format and one theme: what is under test is that both registries
    are walked, not how many surfaces they are walked for.
    """
    from benchmarks.figures import TOOLPATHS_FIGURE_NAME

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        written = regenerate_every_figure(tmp_path, formats=("svg",), themes=(Theme.LIGHT,))
    stems = {path.stem for path in written}
    assert TOOLPATHS_FIGURE_NAME in stems, "the tool-path comparison must still be drawn"
    assert {name for name, _writer in QUALITY_FIGURES} <= stems, "every quality figure must be drawn"
    for path in written:
        assert path.exists() and path.stat().st_size > 0
