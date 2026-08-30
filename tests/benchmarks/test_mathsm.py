from __future__ import annotations

import pytest

from benchmarks.families.analytic import rectangle
from benchmarks.mathsm import SPACING_SWEEP_TOOL_DIAMETERS
from benchmarks.mathsm import MathsmPoint
from benchmarks.mathsm import shortest_within_cap
from benchmarks.mathsm import sweep_spacing
from benchmarks.pathmetrics import PathMetrics

# The pocket every timed test in this module uses: two chains, a handful of
# seconds to sweep, and wide enough that the tool has room at every spacing.
TEST_WIDTH = 8.0
TEST_HEIGHT = 6.0
TEST_TOOL = 2.0

# Two spacings a factor of two apart, so the length ordering they must produce is
# unambiguous rather than a near-tie the audit could reorder.
COARSE_PAIR = (0.2, 0.4)

# Four spacings spanning the useful range of the test pocket, used only for the
# does-the-axis-move check where the SPREAD of the measurements is the assertion.
SPREAD_SWEEP = (0.1, 0.2, 0.4, 0.6)

# Smallest spread, in degrees, that counts as the spacing axis genuinely moving
# the engagement rather than jittering. Ten degrees is far outside anything
# station placement can produce and far inside the ~85 degree spread measured.
MEANINGFUL_TEA_SPREAD_DEG = 10.0


def _point(spacing: float, length: float, max_tea_after_entry_deg: float, max_tea_deg: float = 360.0) -> MathsmPoint:
    """A trial with hand-chosen numbers, so selection is tested without generating."""
    metrics = PathMetrics(length=length, cut_motions=10, entry_cuts=1, max_tea_deg=max_tea_deg, max_tea_after_entry_deg=max_tea_after_entry_deg)
    return MathsmPoint(spacing_tool_diameters=spacing, metrics=metrics)


def test_the_default_sweep_is_strictly_increasing_and_scale_free() -> None:
    """Spacings are tool diameters, so the sweep means the same at every scale."""
    assert list(SPACING_SWEEP_TOOL_DIAMETERS) == sorted(SPACING_SWEEP_TOOL_DIAMETERS)
    assert len(set(SPACING_SWEEP_TOOL_DIAMETERS)) == len(SPACING_SWEEP_TOOL_DIAMETERS)
    assert all(s > 0.0 for s in SPACING_SWEEP_TOOL_DIAMETERS)


def test_shortest_within_cap_picks_the_shortest_compliant_trial() -> None:
    points = [
        _point(0.1, length=500.0, max_tea_after_entry_deg=60.0),
        _point(0.3, length=300.0, max_tea_after_entry_deg=110.0),
        _point(0.6, length=200.0, max_tea_after_entry_deg=150.0),
    ]
    chosen = shortest_within_cap(points, cap_deg=120.0)
    assert chosen is not None
    assert chosen.spacing_tool_diameters == pytest.approx(0.3)


def test_shortest_within_cap_returns_none_when_nothing_complies() -> None:
    assert shortest_within_cap([_point(0.6, length=200.0, max_tea_after_entry_deg=150.0)], cap_deg=120.0) is None


def test_selection_reads_the_after_entry_maximum_not_the_raw_one() -> None:
    """Filtering on the raw maximum would empty the comparison for every cap.

    Every trial's raw maximum is the full slot its entry cut takes, so a cap below
    360 degrees would reject all of them and the baseline curve would silently not
    exist. The regression this guards is a comparison that cannot show a result.
    """
    points = [_point(0.4, length=200.0, max_tea_after_entry_deg=90.0, max_tea_deg=360.0)]
    chosen = shortest_within_cap(points, cap_deg=120.0)
    assert chosen is not None
    assert chosen.spacing_tool_diameters == pytest.approx(0.4)


def test_a_cap_at_the_measured_maximum_is_compliant() -> None:
    """The cap is an upper bound the path may reach, so the comparison is closed."""
    points = [_point(0.4, length=200.0, max_tea_after_entry_deg=120.0)]
    assert shortest_within_cap(points, cap_deg=120.0) is not None


def test_sweep_measures_every_requested_spacing_in_order() -> None:
    spec = rectangle(width=TEST_WIDTH, height=TEST_HEIGHT, tool_diameter=TEST_TOOL, tea_cap_deg=120.0)
    points = sweep_spacing(spec, COARSE_PAIR)
    assert [p.spacing_tool_diameters for p in points] == list(COARSE_PAIR)
    assert all(p.metrics.length > 0.0 for p in points)
    assert all(p.metrics.cut_motions > 0 for p in points)


def test_a_coarser_spacing_gives_a_shorter_path() -> None:
    """Fewer machining circles cover the same pocket, which is the whole trade."""
    spec = rectangle(width=TEST_WIDTH, height=TEST_HEIGHT, tool_diameter=TEST_TOOL, tea_cap_deg=120.0)
    fine, coarse = sweep_spacing(spec, COARSE_PAIR)
    assert coarse.metrics.length < fine.metrics.length


def test_the_spacing_axis_actually_moves_the_measured_engagement() -> None:
    """Without this the baseline curve would be a constant dressed as a sweep.

    This test asserts meaningful spread only; it does not establish whether the
    sampled relation is monotone. `shortest_within_cap` minimises over all
    compliant trials because the generic protocol assumes no spacing order.
    """
    spec = rectangle(width=TEST_WIDTH, height=TEST_HEIGHT, tool_diameter=TEST_TOOL, tea_cap_deg=120.0)
    measured = [p.metrics.max_tea_after_entry_deg for p in sweep_spacing(spec, SPREAD_SWEEP)]
    assert max(measured) - min(measured) > MEANINGFUL_TEA_SPREAD_DEG


def test_every_trial_is_measured_at_the_same_reference_cap() -> None:
    """A trial's engagement must be a property of its spacing, not of a cap.

    The audit refines line motions against the cap it is given, so measuring each
    spacing at the cap it will later be compared against would let the x-axis of
    the baseline curve shift with the y-axis it is plotted against.
    """
    tight = rectangle(width=TEST_WIDTH, height=TEST_HEIGHT, tool_diameter=TEST_TOOL, tea_cap_deg=30.0)
    loose = rectangle(width=TEST_WIDTH, height=TEST_HEIGHT, tool_diameter=TEST_TOOL, tea_cap_deg=170.0)
    assert sweep_spacing(tight, COARSE_PAIR)[0].metrics == sweep_spacing(loose, COARSE_PAIR)[0].metrics
