"""Exact NOTHING-circle removal: a circle whose swept annulus lies in the union of earlier circles' swept annuli."""

import math

import pytest

from benchmarks.held_nothing_circles import remove_nothing_circles
from benchmarks.held_nothing_circles import removes_nothing
from benchmarks.held_nothing_circles import union_covers_annulus
from benchmarks.held_standard_placement import PaperCircleCandidate
from benchmarks.held_standard_placement import exact_predecessor_engagement
from compas_cgal.motion import EngagementCap
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

TOOL = ToolRadius.build(1.0)
TIGHT = EngagementCap.build(math.radians(80.0))
LOOSE = EngagementCap.build(math.radians(160.0))


def _circle(x: float, y: float, radius: float) -> PaperCircleCandidate:
    return PaperCircleCandidate.build(
        center=Point2[WorldXY].build(x, y),
        guide_radius=GuideRadius.build(radius),
        contact_point=Point2[WorldXY].build(x, y - radius),
    )


def _identity(circle: PaperCircleCandidate) -> PaperCircleCandidate:
    return circle


# Annulus [2, 4] about the origin.
BIG = _circle(0.0, 0.0, 3.0)
# Annulus [0.5, 2.5] about the origin: with BIG the union is the annulus [0.5, 4].
SMALL = _circle(0.0, 0.0, 1.5)
# Annulus [1.25, 3.25] about (0.5, 0): in [0.5, 4] about the origin, in neither BIG nor SMALL alone.
INSIDE_UNION = _circle(0.5, 0.0, 2.25)
# Its successor: within 80 deg after INSIDE_UNION, beyond 80 deg (within 160 deg) after SMALL.
SUCCESSOR = _circle(1.0, 0.0, 2.25)


def test_exact_duplicate_removes_nothing_and_is_removed() -> None:
    second = _circle(0.5, 0.0, 3.0)
    duplicate = _circle(0.5, 0.0, 3.0)
    assert removes_nothing(duplicate, (BIG, second), TOOL)
    result = remove_nothing_circles((BIG, second, duplicate), _identity, TOOL, TIGHT)
    assert result.circles == (BIG, second)
    assert result.removed == (2,)


def test_one_ulp_shift_leaves_a_sliver_and_is_kept() -> None:
    # Shifted by the smallest positive double: the outer rim pokes out of BIG's annulus.
    shifted = _circle(math.nextafter(0.0, 1.0), 0.0, 3.0)
    assert not removes_nothing(shifted, (BIG,), TOOL)
    result = remove_nothing_circles((BIG, shifted), _identity, TOOL, TIGHT)
    assert result.circles == (BIG, shifted)
    assert result.removed == ()


def test_diagonal_sliver_is_refused_by_the_exact_union_test() -> None:
    # A shift along (1, 1): the exact arrangement test, not only the point witnesses, must refuse it.
    delta = 2.0**-40
    shifted = _circle(delta, delta, 3.0)
    assert not union_covers_annulus(shifted, (BIG,), TOOL)


def test_circle_covered_only_by_the_union_of_two_predecessors_is_removed() -> None:
    assert not removes_nothing(INSIDE_UNION, (BIG,), TOOL)
    assert not removes_nothing(INSIDE_UNION, (SMALL,), TOOL)
    assert removes_nothing(INSIDE_UNION, (BIG, SMALL), TOOL)


def test_successor_cap_is_re_decided_against_the_new_predecessor() -> None:
    # Fixture premises, decided by the same exact predicate the pass uses.
    assert exact_predecessor_engagement(INSIDE_UNION, SUCCESSOR, TOOL).within(TIGHT)
    assert not exact_predecessor_engagement(SMALL, SUCCESSOR, TOOL).within(TIGHT)
    assert exact_predecessor_engagement(SMALL, SUCCESSOR, TOOL).within(LOOSE)
    circles = (BIG, SMALL, INSIDE_UNION, SUCCESSOR)

    tight = remove_nothing_circles(circles, _identity, TOOL, TIGHT)
    assert tight.circles == circles
    assert tight.removed == ()
    assert tight.kept_for_successor_cap == (2,)

    loose = remove_nothing_circles(circles, _identity, TOOL, LOOSE)
    assert loose.circles == (BIG, SMALL, SUCCESSOR)
    assert loose.removed == (2,)
    assert loose.kept_for_successor_cap == ()


def test_first_circle_is_never_removed() -> None:
    result = remove_nothing_circles((BIG, _circle(0.0, 0.0, 3.0)), _identity, TOOL, TIGHT)
    assert result.removed == (1,)
    assert result.circles[0] is BIG


@pytest.mark.parametrize("count", [0, 1])
def test_paths_without_a_predecessor_pass_through(count: int) -> None:
    circles = (BIG,)[:count]
    result = remove_nothing_circles(circles, _identity, TOOL, TIGHT)
    assert result.circles == circles
    assert result.removed == ()
