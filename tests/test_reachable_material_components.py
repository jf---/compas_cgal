"""The reachable-material predicate reports the cutter-centre domain's component count."""

from __future__ import annotations

import numpy as np
import pytest

from compas_cgal import _coverage_2

# A 12 x 8 rectangle and a tool of radius 1: the centre domain is the inner
# 10 x 6 rectangle, one component.
RECTANGLE = np.array([[0.0, 0.0], [12.0, 0.0], [12.0, 8.0], [0.0, 8.0]], dtype=float)
TOOL_RADIUS = 1.0
# One unit past the rectangle's half-height of 4: no centre is a tool radius
# clear of both long walls, so the centre domain is empty.
SLOT_FILLING_RADIUS = 4.0 + 1.0


def test_a_rectangle_has_one_centre_domain_component() -> None:
    predicate = _coverage_2.ReachableMaterialPredicate2.build(RECTANGLE, [], TOOL_RADIUS)
    assert predicate.center_component_count() == 1


def test_a_tool_too_large_for_the_pocket_is_refused_by_name_before_any_count() -> None:
    with pytest.raises(_coverage_2.PocketNotMachinableError):
        _coverage_2.ReachableMaterialPredicate2.build(RECTANGLE, [], SLOT_FILLING_RADIUS)
