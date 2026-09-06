"""Consumer contracts for exact world-XY millimetre boundary contacts."""

from __future__ import annotations

from collections.abc import Sequence
from math import sqrt

import numpy as np
import pytest

from compas_cgal import _coverage_2


def _cycle(shift_mm: float = 0.0) -> _coverage_2.ReachableBoundaryCycle2:
    # Integer rotation of an L: offsetting its diagonal sides by 1 mm creates
    # irrational contacts (coordinate displacements of sqrt(1/2) mm).
    boundary = np.asarray(
        [(0, 0, 0), (12, 12, 0), (6, 18, 0), (0, 12, 0), (-6, 18, 0), (-12, 12, 0)],
        dtype=np.float64,
    )
    boundary[:, 0] += shift_mm
    return _coverage_2.build_center_boundary_cycle(boundary, [], 1.0)


def _assert_retained(
    actual: Sequence[_coverage_2.ReachableBoundaryPrimitive2],
    expected: Sequence[_coverage_2.ReachableBoundaryPrimitive2],
) -> None:
    assert len(actual) == len(expected)
    for result, source in zip(actual, expected):
        assert result.kind == source.kind
        assert result.start == source.start
        assert result.end == source.end
        assert result.source_piece_records == source.source_piece_records
        assert result.source_piece_records
    assert all(first.end == second.start for first, second in zip(actual, actual[1:]))


def test_native_contacts_preserve_algebraic_line_arc_transition() -> None:
    cycle = _cycle()
    primitives = cycle.primitives
    assert cycle.counterclockwise
    assert {item.kind for item in primitives} == {"line", "arc"}
    arc = next(item for item in primitives if item.kind == "arc")
    # A geometric witness, not a tolerance used by the transition: this endpoint
    # has an irrational coordinate and cannot survive rational double injection.
    assert abs(arc.start_mm[1] - arc.arc_center_mm[1]) == pytest.approx(sqrt(0.5))
    assert arc.start.reporting_xy_mm == arc.start_mm
    assert arc.end.reporting_xy_mm == arc.end_mm
    assert isinstance(arc.start, _coverage_2.WorldXYBoundaryPointMm)
    assert all(first.end == second.start for first, second in zip(primitives, (*primitives[1:], primitives[0])))
    assert primitives[0].start != primitives[0].end

    # Every start is handed back to exact native membership, including both
    # algebraic arc contacts; converting any contact to doubles would fail.
    for primitive in primitives:
        transition = cycle.ccw_transition(primitive.start, primitive.end)
        _assert_retained(transition, [primitive])

    start = primitives[-1].start
    end = primitives[1].end
    transition = cycle.ccw_transition(start, end)
    _assert_retained(transition, [primitives[-1], *primitives[:2]])
    assert transition[0].start == start
    assert transition[-1].end == end


def test_native_points_are_values_across_independent_cycles() -> None:
    first = _cycle()
    second = _cycle()
    assert first.primitives[0].start == second.primitives[0].start
    transition = first.ccw_transition(second.primitives[0].start, second.primitives[0].end)
    _assert_retained(transition, [first.primitives[0]])


def test_ccw_transition_rejects_zero_progress_and_foreign_boundary() -> None:
    cycle = _cycle()
    start = cycle.primitives[0].start
    with pytest.raises(_coverage_2.ReachableArrangementTopologyError, match="positive boundary progress"):
        cycle.ccw_transition(start, start)
    foreign = _cycle(100.0).primitives[0].start
    with pytest.raises(_coverage_2.ReachableArrangementTopologyError, match="does not lie"):
        cycle.ccw_transition(start, foreign)
    with pytest.raises(_coverage_2.ReachableArrangementTopologyError, match="does not lie"):
        cycle.ccw_transition(foreign, start)


def test_ccw_transition_rejects_raw_coordinates_and_point_construction() -> None:
    cycle = _cycle()
    primitive = cycle.primitives[0]
    with pytest.raises(TypeError):
        cycle.ccw_transition(primitive.start_mm, primitive.end)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        cycle.ccw_transition(primitive.start, primitive.end_mm)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        _coverage_2.WorldXYBoundaryPointMm()  # type: ignore[call-arg]
