from __future__ import annotations

import numpy as np
import pytest

from benchmarks.depletion import replay_depletion
from benchmarks.errors import UnreplayableOperationError
from benchmarks.families.analytic import rectangle
from benchmarks.instrument import probe_digits
from benchmarks.instrument import probe_size
from benchmarks.runner import generate_toolpath
from benchmarks.spec import PocketSpec
from compas.geometry import Line
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult

# The synthetic operations below never reach the polyline field, which exists only
# for visualization; an empty (0, 3) array is its well-formed zero value.
EMPTY_POLYLINE = np.zeros((0, 3), dtype=np.float64)

# Height of the synthetic clearance-plane link, in the pocket's own units. Any
# value strictly above the cut plane exercises the same classification branch.
LINK_Z = 2.0


def _pocket() -> PocketSpec:
    """The smallest instance that still produces a many-operation toolpath."""
    return rectangle(width=8.0, height=6.0, tool_diameter=2.0, tea_cap_deg=120.0)


def _cut_across_the_middle() -> ToolpathOperation:
    """A cut-plane linear motion that removes a capsule of material."""
    return ToolpathOperation(geometry=Line([-2.0, 0.0, 0.0], [2.0, 0.0, 0.0]), operation=OperationType.CUT, path_index=0)


def test_depletion_grows_the_arrangement_and_the_exact_coordinates() -> None:
    # The regression this file exists for: reading the kernel probes on a fresh,
    # undepleted stock reports the virgin values forever, so the two diagnostics
    # can never move. Both must be strictly larger after the toolpath is replayed.
    spec = _pocket()
    virgin = Stock(spec.polygon, list(spec.holes))
    depleted = replay_depletion(spec, generate_toolpath(spec))
    assert probe_size(depleted).vertices > probe_size(virgin).vertices
    assert probe_digits(depleted).max_digits > probe_digits(virgin).max_digits


def test_clearance_height_links_remove_no_material() -> None:
    spec = _pocket()
    cut = _cut_across_the_middle()
    link = ToolpathOperation(geometry=Line([2.0, 0.0, LINK_Z], [-2.0, 2.0, LINK_Z]), operation=OperationType.LINK, path_index=0)
    with_link = replay_depletion(spec, ToolpathResult(operations=[cut, link], polyline=EMPTY_POLYLINE))
    without_link = replay_depletion(spec, ToolpathResult(operations=[cut], polyline=EMPTY_POLYLINE))
    assert with_link.exactly_equals(without_link)


def test_retracts_remove_no_material() -> None:
    spec = _pocket()
    cut = _cut_across_the_middle()
    retract = ToolpathOperation(geometry=Line([2.0, 0.0, 0.0], [2.0, 0.0, LINK_Z]), operation=OperationType.RETRACT, path_index=0)
    with_retract = replay_depletion(spec, ToolpathResult(operations=[cut, retract], polyline=EMPTY_POLYLINE))
    without_retract = replay_depletion(spec, ToolpathResult(operations=[cut], polyline=EMPTY_POLYLINE))
    assert with_retract.exactly_equals(without_retract)


def test_a_plunge_removes_exactly_the_tool_disk() -> None:
    spec = _pocket()
    plunge = ToolpathOperation(geometry=Line([0.0, 0.0, LINK_Z], [0.0, 0.0, 0.0]), operation=OperationType.PLUNGE, path_index=0)
    depleted = replay_depletion(spec, ToolpathResult(operations=[plunge], polyline=EMPTY_POLYLINE))
    expected = Stock(spec.polygon, list(spec.holes))
    expected.subtract_disk(0.0, 0.0, spec.tool_radius)
    assert depleted.exactly_equals(expected)


def test_a_ramped_cut_is_refused_rather_than_mis_depleted() -> None:
    # A line that descends while travelling in XY is outside the cut-plane replay
    # model; projecting it to a flat capsule would remove material the machine
    # never cut, so the replay must refuse instead of guessing.
    spec = _pocket()
    ramp = ToolpathOperation(geometry=Line([-2.0, 0.0, LINK_Z], [2.0, 0.0, 0.0]), operation=OperationType.CUT, path_index=0)
    with pytest.raises(UnreplayableOperationError):
        replay_depletion(spec, ToolpathResult(operations=[ramp], polyline=EMPTY_POLYLINE))
