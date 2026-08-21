from __future__ import annotations

import math

import numpy as np
import pytest
from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Frame
from compas.geometry import Line

from benchmarks.families.analytic import rectangle
from benchmarks.pathmetrics import REFERENCE_CAP_DEG
from benchmarks.pathmetrics import UnmeasurableOperationLengthError
from benchmarks.pathmetrics import entry_cut_indices
from benchmarks.pathmetrics import max_tea_after_entry
from benchmarks.pathmetrics import measure_path
from benchmarks.pathmetrics import path_length
from compas_cgal.engagement import OperationEngagement
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult
from compas_cgal.toolpath import trochoidal_mat_toolpath_circular

# A pocket small enough to audit in about a second, large enough that the
# generator emits several chains and therefore several entry cuts.
TEST_WIDTH = 8.0
TEST_HEIGHT = 6.0
TEST_TOOL = 2.0

# Constant spacing, in tool diameters, whose entry cut is a full slot on the test
# pocket while its steady-state cutting is far below one. Measured, not guessed.
FULL_SLOT_SPACING = 0.2

# The tessellated polyline chords a circle from inside, so the analytic length is
# the larger of the two. 0.5% is two orders above the measured 0.04% gap at the
# generator's default of ten samples per radian.
TESSELLATION_SLACK = 0.005


def _op(geometry, operation: OperationType) -> ToolpathOperation:
    """A toolpath operation with the metadata the classifier never reads."""
    return ToolpathOperation(geometry=geometry, operation=operation, path_index=0)


def _plunge(x: float, y: float, z_top: float = 4.0) -> ToolpathOperation:
    """A downward bore from the clearance height to the cutting plane."""
    return _op(Line([x, y, z_top], [x, y, 0.0]), OperationType.PLUNGE)


def _cut_circle(x: float, y: float, radius: float = 1.0) -> ToolpathOperation:
    """A machining circle at the cutting plane."""
    return _op(Circle(radius, frame=Frame([x, y, 0.0])), OperationType.CUT)


def _retract(x: float, y: float, z_top: float = 4.0) -> ToolpathOperation:
    """An upward rapid to the clearance height."""
    return _op(Line([x, y, 0.0], [x, y, z_top]), OperationType.RETRACT)


def _result(operations: list[ToolpathOperation]) -> ToolpathResult:
    """A synthetic toolpath; the polyline is unused by every metric under test."""
    return ToolpathResult(operations=operations, polyline=np.zeros((0, 3), dtype=float))


def _engagement(op_index: int, max_tea_deg: float, certified: bool = True) -> OperationEngagement:
    """One audited operation, built directly so no toolpath has to be generated."""
    return OperationEngagement(op_index=op_index, operation=OperationType.CUT, max_tea=math.radians(max_tea_deg), cap_certified=certified, stations=20)


def test_the_first_cut_after_a_plunge_is_the_entry() -> None:
    result = _result([_plunge(0.0, 0.0), _cut_circle(0.0, 0.0), _cut_circle(1.0, 0.0), _cut_circle(2.0, 0.0)])
    assert entry_cut_indices(result) == frozenset({1})


def test_each_plunge_arms_exactly_one_further_entry() -> None:
    result = _result([_plunge(0.0, 0.0), _cut_circle(0.0, 0.0), _retract(0.0, 0.0), _plunge(5.0, 0.0), _cut_circle(5.0, 0.0), _cut_circle(6.0, 0.0)])
    assert entry_cut_indices(result) == frozenset({1, 4})


def test_rapid_moves_between_the_plunge_and_the_cut_do_not_disarm_it() -> None:
    """A rapid removes no material, so the cut after it still meets virgin stock."""
    result = _result([_plunge(0.0, 0.0), _retract(0.0, 0.0), _cut_circle(0.0, 0.0)])
    assert entry_cut_indices(result) == frozenset({2})


def test_a_path_that_never_plunges_has_no_entry_cut() -> None:
    result = _result([_cut_circle(0.0, 0.0), _cut_circle(1.0, 0.0)])
    assert entry_cut_indices(result) == frozenset()


def test_path_length_sums_the_exact_primitive_lengths() -> None:
    """Analytic, never the tessellation: a circle contributes its circumference."""
    result = _result([_op(Line([0.0, 0.0, 0.0], [3.0, 4.0, 0.0]), OperationType.CUT), _cut_circle(0.0, 0.0, radius=2.0)])
    assert path_length(result) == pytest.approx(5.0 + 4.0 * math.pi)


def test_path_length_counts_an_arc_by_radius_times_swept_angle() -> None:
    quarter = Arc(radius=2.0, start_angle=0.0, end_angle=0.5 * math.pi, frame=Frame([0.0, 0.0, 0.0]))
    assert path_length(_result([_op(quarter, OperationType.CUT)])) == pytest.approx(math.pi)


def test_an_unknown_primitive_fails_loudly_rather_than_measuring_zero() -> None:
    with pytest.raises(UnmeasurableOperationLengthError):
        path_length(_result([_op(Frame.worldXY(), OperationType.CUT)]))


def test_max_tea_after_entry_ignores_the_entry_operations() -> None:
    operations = [_engagement(1, 360.0), _engagement(2, 95.0), _engagement(3, 40.0)]
    assert math.degrees(max_tea_after_entry(operations, frozenset({1}))) == pytest.approx(95.0)


def test_max_tea_after_entry_is_zero_when_every_cut_is_an_entry() -> None:
    assert max_tea_after_entry([_engagement(1, 360.0)], frozenset({1})) == pytest.approx(0.0)


def test_measure_path_reports_the_analytic_length_the_tessellation_approximates() -> None:
    """A contract test at the seam: the polyline is a chorded copy of the path.

    The two are computed by different machinery, so agreement to a tessellation
    bound is real evidence that neither has silently changed meaning.
    """
    spec = rectangle(width=TEST_WIDTH, height=TEST_HEIGHT, tool_diameter=TEST_TOOL, tea_cap_deg=REFERENCE_CAP_DEG)
    result = _constant_spacing_path(spec, FULL_SLOT_SPACING)
    metrics = measure_path(spec, result)
    chords = float(np.linalg.norm(np.diff(np.asarray(result.polyline, dtype=float), axis=0), axis=1).sum())
    assert metrics.length >= chords
    assert metrics.length == pytest.approx(chords, rel=TESSELLATION_SLACK)


def test_the_entry_cut_is_what_pins_the_raw_maximum_at_a_full_slot() -> None:
    """The load-bearing claim behind every after-entry number in the corpus.

    A tool entering solid stock cuts a full slot on its first circle whatever the
    stepover is, so the raw maximum carries no information about the stepover
    control and a comparison filtered on it is empty for every generator. The
    after-entry maximum is the one that moves.
    """
    spec = rectangle(width=TEST_WIDTH, height=TEST_HEIGHT, tool_diameter=TEST_TOOL, tea_cap_deg=REFERENCE_CAP_DEG)
    metrics = measure_path(spec, _constant_spacing_path(spec, FULL_SLOT_SPACING))
    assert metrics.entry_cuts > 0
    assert metrics.max_tea_deg == pytest.approx(360.0)
    assert metrics.max_tea_after_entry_deg < metrics.max_tea_deg


def _constant_spacing_path(spec, spacing_tool_diameters: float) -> ToolpathResult:
    """The unregulated generator at one constant spacing.

    Called directly rather than through `benchmarks.mathsm` so these tests hold
    the path-metric contract on its own, independent of the baseline that uses it.
    """
    return trochoidal_mat_toolpath_circular(
        spec.polygon,
        tool_diameter=spec.tool_diameter,
        stepover=spacing_tool_diameters * spec.tool_diameter,
        holes=list(spec.holes),
        clearance_z=2.0 * spec.tool_diameter,
    )
