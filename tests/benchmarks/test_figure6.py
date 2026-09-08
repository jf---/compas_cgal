from __future__ import annotations

import json
from pathlib import Path

import pytest

import benchmarks.figure6 as figure6_module
from benchmarks.families.analytic import rectangle
from benchmarks.figure6 import NO_BASELINE_CELL
from benchmarks.figure6 import Figure6Point
from benchmarks.figure6 import Figure6Run
from benchmarks.figure6 import figure6_payload
from benchmarks.figure6 import figure6_points
from benchmarks.figure6 import render_figure6_markdown
from benchmarks.figure6 import reference_pocket
from benchmarks.figure6 import run_figure6
from benchmarks.figure6 import write_figure6
from benchmarks.mathsm import MathsmPoint
from benchmarks.pathmetrics import PathMetrics
from benchmarks.held_reference_cases import load_held_reference_case

# Small enough to run both generators several times per test, large enough that
# the generator emits more than one chain and therefore more than one entry cut.
TEST_WIDTH = 8.0
TEST_HEIGHT = 6.0
TEST_TOOL = 2.0

# Two caps far enough apart that the advance regulation has to react, so a flat
# result is a real finding and not two samples of the same operating point.
TEST_CAPS = (60.0, 160.0)

# Two trial spacings a factor of two apart: enough for a baseline to exist and be
# selected, without paying for the twelve-point default sweep in a unit test.
TEST_SPACINGS = (0.2, 0.4)

# Banned as literal text: the practice is to ORDER the document conclusion-first,
# not to announce that you have (CLAUDE.md, documentation format).
BANNED_LABELS = ("BLUF", "TL;DR", "Bottom Line", "## Summary")


def test_figure6_uses_the_reconstructed_figure5_pocket() -> None:
    spec = reference_pocket()
    expected = load_held_reference_case("figure5").pocket_spec()

    assert tuple(tuple(point) for point in spec.polygon.points) == tuple(tuple(point) for point in expected.polygon.points)
    assert spec.tool_diameter == pytest.approx(2.0)


def test_figure6_reference_pocket_rebuilds_the_requested_cap() -> None:
    spec = reference_pocket(tea_cap_deg=140.0)

    assert spec.tea_cap_deg == pytest.approx(140.0)
    assert spec.name == "figure5"


def _spec(cap_deg: float = 120.0):
    """The pocket every test in this module measures."""
    return rectangle(width=TEST_WIDTH, height=TEST_HEIGHT, tool_diameter=TEST_TOOL, tea_cap_deg=cap_deg)


def _metrics(length: float, max_tea_after_entry_deg: float) -> PathMetrics:
    """Path metrics with hand-chosen numbers, so rendering is tested without generating."""
    return PathMetrics(length=length, cut_motions=20, entry_cuts=2, max_tea_deg=360.0, max_tea_after_entry_deg=max_tea_after_entry_deg)


def _trial(spacing: float, length: float, max_tea_after_entry_deg: float) -> MathsmPoint:
    return MathsmPoint(spacing_tool_diameters=spacing, metrics=_metrics(length, max_tea_after_entry_deg))


def _point(cap_deg: float, ours: float, theirs: float | None, ours_tea: float | None = None, exceedances: int = 0) -> Figure6Point:
    baseline = None if theirs is None else _trial(0.2, theirs, cap_deg)
    return Figure6Point(
        cap_deg=cap_deg,
        controlled=_metrics(ours, cap_deg if ours_tea is None else ours_tea),
        controlled_exceedances_after_entry=exceedances,
        mathsm=baseline,
    )


def _run(points: list[Figure6Point], trials: list[MathsmPoint] | None = None) -> Figure6Run:
    return Figure6Run(spec=_spec(), points=tuple(points), trials=tuple(trials or ()))


def test_points_cover_every_requested_cap_in_order() -> None:
    spec = _spec()
    trials = [_trial(0.2, 400.0, 90.0), _trial(0.4, 250.0, 150.0)]
    points = figure6_points(spec, caps=TEST_CAPS, trials=trials)
    assert [p.cap_deg for p in points] == list(TEST_CAPS)


def test_the_cap_axis_actually_moves_our_path_length() -> None:
    """Without this the figure would be two samples of one operating point.

    A benchmark whose x-axis does not move its y-axis cannot show the effect it
    exists to show, and the failure is invisible in a table of plausible numbers.
    """
    spec = _spec()
    tight, loose = figure6_points(spec, caps=TEST_CAPS, trials=[])
    assert tight.controlled.length > loose.controlled.length


def test_the_controlled_path_honours_its_cap_away_from_the_entry_on_this_pocket() -> None:
    """The regulation works here, and the measurement has to be able to see that.

    On this pocket the advance search finds an admissible step everywhere, so the
    measured maximum tracks the cap from below (49.5 of 60, 158.8 of 160) and the
    exact predicate never fires away from an entry. A measurement pipeline that
    reported otherwise would be reading the wrong motions.
    """
    spec = _spec()
    for point in figure6_points(spec, caps=TEST_CAPS, trials=[]):
        assert point.controlled.entry_cuts > 0
        assert point.controlled.max_tea_deg == pytest.approx(360.0)
        assert point.controlled_meets_cap, point.cap_deg
        assert point.controlled_exceedances_after_entry == 0, point.cap_deg


def test_a_looser_cap_lets_the_tool_engage_more() -> None:
    """The cap has to bind the engagement, not just the length."""
    tight, loose = figure6_points(_spec(), caps=TEST_CAPS, trials=[])
    assert tight.controlled.max_tea_after_entry_deg < loose.controlled.max_tea_after_entry_deg


def test_the_baseline_is_selected_per_cap_from_one_shared_sweep() -> None:
    spec = _spec()
    trials = [_trial(0.2, 400.0, 90.0), _trial(0.4, 250.0, 150.0)]
    tight, loose = figure6_points(spec, caps=(100.0, 160.0), trials=trials)
    assert tight.mathsm is not None and tight.mathsm.spacing_tool_diameters == pytest.approx(0.2)
    assert loose.mathsm is not None and loose.mathsm.spacing_tool_diameters == pytest.approx(0.4)


def test_length_ratio_is_none_without_a_baseline() -> None:
    assert _point(40.0, ours=300.0, theirs=None).length_ratio is None


def test_markdown_states_the_comparison_before_the_table() -> None:
    text = render_figure6_markdown(_run([_point(40.0, 300.0, 400.0), _point(80.0, 200.0, 260.0)]))
    assert text.index("shorter") < text.index("| cap (deg) |")


def test_no_banned_summary_labels_appear() -> None:
    text = render_figure6_markdown(_run([_point(40.0, 300.0, 400.0)]))
    for banned in BANNED_LABELS:
        assert banned not in text


def test_markdown_reports_a_missing_baseline_honestly() -> None:
    text = render_figure6_markdown(_run([_point(40.0, 300.0, None)]))
    assert NO_BASELINE_CELL in text
    assert "no compliant spacing" in text


def test_length_ratio_is_reported_per_cap() -> None:
    point = _point(40.0, ours=300.0, theirs=400.0)
    text = render_figure6_markdown(_run([point]))
    assert point.length_ratio == pytest.approx(0.75)
    assert "0.75" in text


def test_markdown_names_both_generators_so_neither_reads_as_the_paper_s() -> None:
    """The paper's curves are not reproducible here; claiming otherwise is the risk."""
    text = render_figure6_markdown(_run([_point(40.0, 300.0, 400.0)]))
    assert "engagement_controlled_toolpath" in text
    assert "trochoidal_mat_toolpath_circular" in text
    assert "Neither is one of the paper's curves." in text


def test_markdown_discloses_that_our_path_may_miss_the_cap_the_baseline_met() -> None:
    """The comparison is asymmetric and the report must not let that pass silently.

    The baseline is SELECTED for compliance; the controlled generator is merely
    ASKED. A ratio between them at a cap the controlled path misses compares a
    compliant path against one that is not, and the row has to say so.
    """
    missed = _point(40.0, ours=300.0, theirs=400.0, ours_tea=127.1)
    met = _point(160.0, ours=200.0, theirs=400.0, ours_tea=120.0)
    text = render_figure6_markdown(_run([missed, met]))
    assert missed.controlled_meets_cap is False
    assert met.controlled_meets_cap is True
    rows = [line for line in text.splitlines() if line.startswith("| 40 |") or line.startswith("| 160 |")]
    assert rows[0].split("|")[4].strip() == "no"
    assert rows[1].split("|")[4].strip() == "yes"
    assert "127.1" in rows[0]
    assert "not compliant" in text


def test_the_conclusion_counts_the_caps_our_own_path_actually_met() -> None:
    met = _point(160.0, ours=300.0, theirs=400.0, ours_tea=120.0)
    missed = _point(40.0, ours=900.0, theirs=1000.0, ours_tea=127.1)
    conclusion = render_figure6_markdown(_run([missed, met])).split("## Per-cap comparison")[0]
    assert "**1 of 2** caps" in conclusion


def test_the_conclusion_singles_out_the_rows_where_both_paths_met_the_cap() -> None:
    """The rest of the table compares a selected path against a merely asked-for one.

    A mean ratio taken over every comparable cap silently mixes those two kinds of
    row, so the subset where both paths were measured under their cap is stated on
    its own.
    """
    met = _point(160.0, ours=300.0, theirs=400.0, ours_tea=120.0)
    missed = _point(40.0, ours=900.0, theirs=1000.0, ours_tea=127.1)
    conclusion = render_figure6_markdown(_run([missed, met])).split("## Per-cap comparison")[0]
    assert "**1** cap(s) where BOTH paths met the cap it is **0.75x**" in conclusion


def test_the_conclusion_says_when_no_row_is_like_for_like() -> None:
    missed = _point(40.0, ours=900.0, theirs=1000.0, ours_tea=127.1)
    conclusion = render_figure6_markdown(_run([missed])).split("## Per-cap comparison")[0]
    assert "No cap had both paths meet it" in conclusion


def test_a_run_with_no_baseline_anywhere_says_so_instead_of_claiming_a_ratio() -> None:
    conclusion = render_figure6_markdown(_run([_point(40.0, 300.0, None)])).split("## Per-cap comparison")[0]
    assert "no length comparison" in conclusion


def test_the_trial_table_is_emitted_as_evidence_for_the_brute_force_search() -> None:
    trials = [_trial(0.2, 400.0, 90.0), _trial(0.4, 250.0, 150.0)]
    text = render_figure6_markdown(_run([_point(120.0, 300.0, 400.0)], trials))
    assert "Constant-spacing trials" in text
    assert text.index("| cap (deg) |") < text.index("| spacing (tool diam.) |")


def test_trial_table_does_not_infer_a_pattern_from_monotone_measurements() -> None:
    trials = [
        _trial(0.2, 400.0, 90.0),
        _trial(0.4, 300.0, 120.0),
        _trial(0.6, 200.0, 150.0),
    ]
    text = render_figure6_markdown(_run([_point(160.0, 300.0, 400.0)], trials))

    assert "falls and rises" not in text
    assert "evaluates every compliant trial" in text
    assert "without assuming spacing orders engagement" in text
    assert "| 0.200 | 400.0 | 20 | 90.0 | 360.0 |" in text
    assert "| 0.400 | 300.0 | 20 | 120.0 | 360.0 |" in text
    assert "| 0.600 | 200.0 | 20 | 150.0 | 360.0 |" in text


def test_payload_round_trips_through_json() -> None:
    run = _run([_point(120.0, 300.0, 400.0), _point(40.0, 900.0, None)], [_trial(0.2, 400.0, 90.0)])
    payload = json.loads(json.dumps(figure6_payload(run)))
    assert [p["cap_deg"] for p in payload["points"]] == [120.0, 40.0]
    assert payload["points"][0]["length_ratio"] == pytest.approx(0.75)
    assert payload["points"][1]["mathsm_length"] is None
    assert payload["pocket"]["name"] == run.spec.name
    assert payload["pocket"]["holes"] == []
    assert len(payload["spacing_trials"]) == 1


def test_write_figure6_emits_both_artifacts(tmp_path) -> None:
    md_path, json_path = write_figure6(_run([_point(120.0, 300.0, 400.0)]), tmp_path)
    assert md_path.exists() and json_path.exists()
    assert "| cap (deg) |" in md_path.read_text()
    assert json.loads(json_path.read_text())["points"][0]["cap_deg"] == 120.0
    assert md_path.read_text(encoding="utf-8").startswith("# Figure 6 reproduction")
    assert json.loads(json_path.read_text(encoding="utf-8"))["pocket"]["holes"] == []


def test_run_figure6_carries_its_own_trials_for_the_report() -> None:
    run = run_figure6(_spec(), caps=(160.0,), spacings_tool_diameters=TEST_SPACINGS)
    assert [t.spacing_tool_diameters for t in run.trials] == list(TEST_SPACINGS)
    assert len(run.points) == 1
    assert run.spec.name == _spec().name
