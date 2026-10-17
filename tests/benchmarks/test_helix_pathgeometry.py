"""A HELIX operation in the figure path: drawn from its own exact parametrisation, with its own mark (Task 1.2).

`benchmarks.pathgeometry` refused a helix (`"helix"` was not a known class) and
`benchmarks.marks` fell through to a ``"<name> (ramped)"`` secondary-ink stroke
for anything it did not name. A helix is a cutting entry with its own mark; an
unknown class raises by name.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any
from typing import List

import pytest
from compas.geometry import Line
from compas.tolerance import TOL

from benchmarks.errors import UnknownOperationClassError
from benchmarks.marks import ColourBy
from benchmarks.marks import encode
from benchmarks.palette import LIGHT
from benchmarks.pathgeometry import ARC_DEGREES_PER_SAMPLE
from benchmarks.pathgeometry import KNOWN_OPERATIONS
from benchmarks.pathgeometry import Motion
from benchmarks.pathgeometry import drawable_motions
from compas_cgal.helix import HELIX_CLOSING_TURNS
from compas_cgal.helix import Helix
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from tests.benchmarks.helix_fixtures import emitted
from tests.benchmarks.helix_fixtures import make_helix

CENTRE = (1.0, 2.0)
FOOTPRINT_RADIUS = 1.0
# An off-axis phase, so the exact entry is not a trivially representable point.
PHASE = 0.7


@dataclass(frozen=True)
class _Result:
    operations: List[Any]


def _helix_operation(clockwise: bool) -> tuple[ToolpathOperation, Helix]:
    helix = make_helix(clockwise, centre=CENTRE, radius=FOOTPRINT_RADIUS, phase=PHASE)
    return emitted(helix), helix


def _signed_turn(points: tuple[tuple[float, float], ...]) -> float:
    """Total signed angle swept about CENTRE, radians (positive = counter-clockwise)."""
    angles = [math.atan2(y - CENTRE[1], x - CENTRE[0]) for x, y in points]
    return sum(math.remainder(b - a, math.tau) for a, b in zip(angles, angles[1:]))


def test_every_operation_type_is_a_known_drawing_class() -> None:
    assert KNOWN_OPERATIONS == frozenset(operation.value for operation in OperationType)


@pytest.mark.parametrize("clockwise", [True, False])
def test_drawable_motions_samples_a_helix(clockwise: bool) -> None:
    operation, helix = _helix_operation(clockwise)
    (motion,) = drawable_motions(_Result(operations=[operation]))
    entry = (helix.entry_xy.x, helix.entry_xy.y)
    assert motion.name == "helix"
    assert motion.is_cutting and not motion.is_point
    assert motion.points[0] == entry and motion.points[-1] == entry
    assert len(motion.points) == math.ceil(360.0 * (helix.turns + HELIX_CLOSING_TURNS) / ARC_DEGREES_PER_SAMPLE) + 1
    assert all(TOL.is_close(math.hypot(x - CENTRE[0], y - CENTRE[1]), FOOTPRINT_RADIUS) for x, y in motion.points)
    assert TOL.is_close(_signed_turn(motion.points), (-1.0 if clockwise else 1.0) * math.tau * (helix.turns + HELIX_CLOSING_TURNS))


def test_a_helix_has_its_own_mark_and_unknown_ops_raise() -> None:
    (motion,) = drawable_motions(_Result(operations=[_helix_operation(True)[0]]))
    mark = encode([motion], ColourBy.OPERATION, LIGHT, None)[motion.index]
    assert mark.label == "helix"
    assert "(ramped)" not in mark.label
    unknown = Motion(index=0, name="dwell", path_index=0, points=((0.0, 0.0), (1.0, 0.0)), length=1.0, is_point=False, is_cutting=False)  # type: ignore[arg-type]
    with pytest.raises(UnknownOperationClassError):
        encode([unknown], ColourBy.OPERATION, LIGHT, None)


def test_a_plunge_with_extent_in_the_plane_is_still_named_a_ramp() -> None:
    ramp = ToolpathOperation(geometry=Line([0.0, 0.0, 4.0], [1.0, 0.0, 0.0]), operation=OperationType.PLUNGE, path_index=0)
    (motion,) = drawable_motions(_Result(operations=[ramp]))
    assert encode([motion], ColourBy.OPERATION, LIGHT, None)[0].label == "plunge (ramped)"
