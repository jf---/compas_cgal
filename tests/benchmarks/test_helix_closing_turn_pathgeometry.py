"""The HELIX op's flat closing turn in the figure path (Task 1.2, user ruling Q6).

The HELIX operation is ``turns`` ramp turns plus one flat closing turn at
``bottom_z``, in the same op. `benchmarks.pathgeometry` drew the ramp only: the
plan view swept ``turns`` turns and nothing showed the flat turn. The plan view
now sweeps the closing turn too, in the helix's own sense, from and to the
exact entry; the side elevation (Z against plan arc length) shows the ramp and
the flat turn, its last ring at ``bottom_z`` exactly.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any
from typing import List

import pytest
from compas.tolerance import TOL

from benchmarks.pathgeometry import ARC_DEGREES_PER_SAMPLE
from benchmarks.pathgeometry import drawable_motions
from benchmarks.pathgeometry import helix_elevation
from compas_cgal.helix import HELIX_CLOSING_TURNS
from tests.benchmarks.helix_fixtures import emitted
from tests.benchmarks.helix_fixtures import make_helix

CENTRE = (1.0, 2.0)
FOOTPRINT_RADIUS = 1.0
# An off-axis phase, so the exact entry is not a trivially representable point.
PHASE = 0.7


@dataclass(frozen=True)
class _Result:
    operations: List[Any]


def _signed_turn(points: tuple[tuple[float, float], ...]) -> float:
    """Total signed angle swept about CENTRE, radians (positive = counter-clockwise)."""
    angles = [math.atan2(y - CENTRE[1], x - CENTRE[0]) for x, y in points]
    return sum(math.remainder(b - a, math.tau) for a, b in zip(angles, angles[1:]))


def test_the_closing_turn_is_one_flat_turn() -> None:
    assert HELIX_CLOSING_TURNS == 1


@pytest.mark.parametrize("clockwise", [True, False])
def test_the_plan_view_sweeps_the_ramp_and_the_closing_turn(clockwise: bool) -> None:
    helix = make_helix(clockwise, centre=CENTRE, radius=FOOTPRINT_RADIUS, phase=PHASE)
    (motion,) = drawable_motions(_Result(operations=[emitted(helix)]))
    swept = helix.turns + HELIX_CLOSING_TURNS
    entry = (helix.entry_xy.x, helix.entry_xy.y)
    assert motion.points[0] == entry and motion.points[-1] == entry
    assert len(motion.points) == math.ceil(360.0 * swept / ARC_DEGREES_PER_SAMPLE) + 1
    assert TOL.is_close(_signed_turn(motion.points), (-1.0 if clockwise else 1.0) * math.tau * swept)


@pytest.mark.parametrize("clockwise", [True, False])
def test_the_elevation_ramps_then_closes_flat_at_bottom_z(clockwise: bool) -> None:
    helix = make_helix(clockwise, centre=CENTRE, radius=FOOTPRINT_RADIUS, phase=PHASE)
    profile = helix_elevation(helix)
    top, bottom = float(helix.top_z.value), float(helix.bottom_z.value)
    circumference = math.tau * float(helix.footprint.radius.value)
    ramp_end = helix.turns * circumference
    assert profile[0] == (0.0, top)
    assert profile[-1] == ((helix.turns + HELIX_CLOSING_TURNS) * circumference, bottom)
    assert all(b[0] > a[0] for a, b in zip(profile, profile[1:]))
    # The last ring -- everything past the ramp -- is at bottom_z exactly.
    assert [z for s, z in profile if s >= ramp_end] == [bottom, bottom]
    assert all(top >= z >= bottom for _, z in profile)
