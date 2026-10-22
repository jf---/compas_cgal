"""Exact NOTHING-circle removal: a circle whose swept annulus lies in the union of earlier circles' swept annuli."""

import math
from fractions import Fraction

import numpy as np
import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks import held_nothing_circles
from benchmarks.held_nothing_circles import material_witness
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
    result = remove_nothing_circles((BIG, second, duplicate), _identity, TOOL, TIGHT, closed=False)
    assert result.circles == (BIG, second)
    assert result.removed == (2,)


def test_one_ulp_shift_leaves_a_sliver_and_is_kept() -> None:
    # Shifted by the smallest positive double: the outer rim pokes out of BIG's annulus.
    shifted = _circle(math.nextafter(0.0, 1.0), 0.0, 3.0)
    assert not removes_nothing(shifted, (BIG,), TOOL)
    result = remove_nothing_circles((BIG, shifted), _identity, TOOL, TIGHT, closed=False)
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

    tight = remove_nothing_circles(circles, _identity, TOOL, TIGHT, closed=False)
    assert tight.circles == circles
    assert tight.removed == ()
    assert tight.kept_for_successor_cap == (2,)

    loose = remove_nothing_circles(circles, _identity, TOOL, LOOSE, closed=False)
    assert loose.circles == (BIG, SMALL, SUCCESSOR)
    assert loose.removed == (2,)
    assert loose.kept_for_successor_cap == ()


def test_first_circle_is_never_removed() -> None:
    result = remove_nothing_circles((BIG, _circle(0.0, 0.0, 3.0)), _identity, TOOL, TIGHT, closed=False)
    assert result.removed == (1,)
    assert result.circles[0] is BIG


@pytest.mark.parametrize("count", [0, 1])
def test_paths_without_a_predecessor_pass_through(count: int) -> None:
    circles = (BIG,)[:count]
    result = remove_nothing_circles(circles, _identity, TOOL, TIGHT, closed=False)
    assert result.circles == circles
    assert result.removed == ()


# Closed paths: the seam runs from the last kept circle back to circle 0.
# SMALL is not NOTHING after (SUCCESSOR, INSIDE_UNION): (0, 0.5) lies in neither annulus.
SEAM_DUPLICATE = _circle(0.5, 0.0, 2.25)


def test_closed_path_keeps_a_last_duplicate_whose_removal_puts_the_seam_over_the_cap() -> None:
    # Fixture premises, decided by the same exact predicate the pass uses.
    assert not removes_nothing(INSIDE_UNION, (SUCCESSOR,), TOOL)
    assert not removes_nothing(SMALL, (SUCCESSOR, INSIDE_UNION), TOOL)
    assert removes_nothing(SEAM_DUPLICATE, (SUCCESSOR, INSIDE_UNION, SMALL), TOOL)
    assert exact_predecessor_engagement(SEAM_DUPLICATE, SUCCESSOR, TOOL).within(TIGHT)
    assert not exact_predecessor_engagement(SMALL, SUCCESSOR, TOOL).within(TIGHT)
    circles = (SUCCESSOR, INSIDE_UNION, SMALL, SEAM_DUPLICATE)

    closed = remove_nothing_circles(circles, _identity, TOOL, TIGHT, closed=True)
    assert closed.circles == circles
    assert closed.removed == ()
    assert closed.kept_for_successor_cap == (3,)

    # Open, the last circle has no successor: the same duplicate is removed.
    opened = remove_nothing_circles(circles, _identity, TOOL, TIGHT, closed=False)
    assert opened.circles == (SUCCESSOR, INSIDE_UNION, SMALL)
    assert opened.removed == (3,)


def test_closed_path_removes_a_last_duplicate_when_the_new_seam_stays_within_the_cap() -> None:
    assert exact_predecessor_engagement(INSIDE_UNION, SUCCESSOR, TOOL).within(TIGHT)
    circles = (SUCCESSOR, INSIDE_UNION, SEAM_DUPLICATE)
    result = remove_nothing_circles(circles, _identity, TOOL, TIGHT, closed=True)
    assert result.circles == (SUCCESSOR, INSIDE_UNION)
    assert result.removed == (2,)
    assert result.kept_for_successor_cap == ()


def test_closed_path_seam_is_re_decided_at_the_looser_cap() -> None:
    circles = (SUCCESSOR, INSIDE_UNION, SMALL, SEAM_DUPLICATE)
    result = remove_nothing_circles(circles, _identity, TOOL, LOOSE, closed=True)
    assert result.circles == (SUCCESSOR, INSIDE_UNION, SMALL)
    assert result.removed == (3,)


def test_closed_path_never_removes_circle_zero_and_decides_a_self_seam_exactly() -> None:
    # Circle 0 is cut first: no earlier union exists, so it is never removed.
    # Removing its exact duplicate leaves the seam circle 0 -> circle 0, which the
    # predicate decides exactly (theta = 0: the successor disk lies in the predecessor's).
    duplicate = _circle(0.0, 0.0, 3.0)
    assert exact_predecessor_engagement(BIG, BIG, TOOL).within(TIGHT)
    result = remove_nothing_circles((BIG, duplicate), _identity, TOOL, TIGHT, closed=True)
    assert result.circles == (BIG,)
    assert result.removed == (1,)


def test_closed_is_a_required_keyword() -> None:
    with pytest.raises(TypeError, match="closed"):
        remove_nothing_circles((BIG, SMALL), _identity, TOOL, TIGHT)  # type: ignore[call-arg]


# The keep screen: a circle is kept without the exact union test only on an exact
# witness, a rational point of its annulus outside every earlier annulus.


def _in_annulus(circle: PaperCircleCandidate, point: tuple[Fraction, Fraction]) -> bool:
    squared = (point[0] - Fraction(circle.center.x)) ** 2 + (point[1] - Fraction(circle.center.y)) ** 2
    guide, tool = Fraction(circle.guide_radius.value), Fraction(TOOL.value)
    return (guide - tool) ** 2 * (guide > tool) <= squared <= (guide + tool) ** 2


def test_a_witness_is_an_exact_point_of_the_annulus_outside_every_earlier_annulus() -> None:
    target = _circle(1.0, 0.0, 3.0)
    witness = material_witness(target, (BIG, SMALL), TOOL)
    assert witness is not None
    assert all(isinstance(value, Fraction) for value in witness)
    assert _in_annulus(target, witness)
    assert not _in_annulus(BIG, witness)
    assert not _in_annulus(SMALL, witness)


def test_a_covered_annulus_has_no_witness() -> None:
    assert material_witness(INSIDE_UNION, (BIG, SMALL), TOOL) is None


def test_a_sliver_below_every_sample_grid_is_kept_on_a_rim_gap_witness() -> None:
    # Shifted by 2^-40 along (1, 1): no sample grid sees the sliver, the rim-gap search does.
    delta = 2.0**-40
    shifted = _circle(delta, delta, 3.0)
    witness = material_witness(shifted, (BIG,), TOOL)
    assert witness is not None
    assert _in_annulus(shifted, witness)
    assert not _in_annulus(BIG, witness)
    result = remove_nothing_circles((BIG, shifted), _identity, TOOL, TIGHT, closed=False)
    assert result.removed == ()
    assert result.union_tests == 0


# Dyadic coordinates so exact coincidences, tangencies and duplicates occur.
_COORDINATE = st.integers(-12, 12).map(lambda n: n / 4.0)
_RADIUS = st.integers(1, 16).map(lambda n: n / 4.0)
_CIRCLES = st.lists(st.builds(_circle, _COORDINATE, _COORDINATE, _RADIUS), min_size=2, max_size=6)


@settings(max_examples=60, deadline=None)
@given(circles=_CIRCLES, data=st.data())
def test_the_screen_never_changes_the_verdict(circles: list[PaperCircleCandidate], data: st.DataObject) -> None:
    index = data.draw(st.integers(1, len(circles) - 1))
    earlier, target = tuple(circles[:index]), circles[index]
    witness = material_witness(target, earlier, TOOL)
    covered = union_covers_annulus(target, earlier, TOOL)
    if witness is not None:
        assert _in_annulus(target, witness)
        assert not any(_in_annulus(circle, witness) for circle in earlier)
        assert not covered
    assert removes_nothing(target, earlier, TOOL) == covered


@settings(max_examples=40, deadline=None)
@given(circles=st.lists(st.builds(_circle, _COORDINATE, _COORDINATE, _RADIUS), min_size=2, max_size=2 * held_nothing_circles.FIRST_UNION_QUERY))
def test_the_batched_screen_never_changes_the_pass(circles: list[PaperCircleCandidate]) -> None:
    screened = remove_nothing_circles(circles, _identity, TOOL, LOOSE, closed=True)
    # Reference: no screen at all, every circle goes to the exact union test.
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(held_nothing_circles, "_batch_witnessed", lambda geometry, *_: np.zeros(len(geometry), dtype=bool))
        patch.setattr(held_nothing_circles, "_witness", lambda *_, **__: None)
        exact = remove_nothing_circles(circles, _identity, TOOL, LOOSE, closed=True)
    assert screened.removed == exact.removed
    assert screened.kept_for_successor_cap == exact.kept_for_successor_cap
    assert screened.union_tests <= exact.union_tests


# Exact boundaries of the witness checks. Earlier E: annulus [2, 4] about the origin;
# target T: disc of radius 2 about (5, 0). (4, 0) lies on E's outer rim and inside T.
RIM_EARLIER = _circle(0.0, 0.0, 3.0)
RIM_TARGET = _circle(5.0, 0.0, 1.0)


def _rim_fixture() -> tuple[np.ndarray, float, held_nothing_circles._DyadicImage, held_nothing_circles._Selection]:
    geometry = held_nothing_circles._arrays((RIM_EARLIER, RIM_TARGET))
    tool = float(TOOL.value)
    selection = held_nothing_circles._Selection.build(geometry, tool)
    assert list(selection.of(1)) == [0]
    return geometry, tool, held_nothing_circles._DyadicImage.build(geometry, tool), selection


@pytest.mark.parametrize(
    ("candidate", "witnessed"),
    [
        ((4.0, 0.0), False),
        ((8.0, 0.0), False),
        ((6.0, 0.0), True),
        ((4.5, 0.0), True),
        ((6.5, 0.0), True),
        ((3.5, 0.0), False),
        ((4.5, 0.25), True),
        ((3.75, 1.5), True),
        ((7.0, 0.0), True),
    ],
    # The fixture's path is all-integer (dyadic scale 0): (4.5, 0) must be checked at its own value, not snapped onto E's rim.
    ids=[
        "on-earlier-rim",
        "outside-target",
        "inside-target-outside-earlier",
        "off-the-path-scale",
        "off-scale-inside-target",
        "off-scale-inside-earlier",
        "off-scale-both-coordinates",
        "off-scale-near-both-rims",
        "on-target-rim",
    ],
)
def test_the_batched_check_treats_both_annuli_as_closed(monkeypatch: pytest.MonkeyPatch, candidate: tuple[float, float], witnessed: bool) -> None:
    geometry, tool, image, selection = _rim_fixture()
    slot = np.full((2, 2), np.nan)
    slot[1] = candidate
    monkeypatch.setattr(held_nothing_circles, "_own_rim_candidates", lambda *_: [slot])
    assert bool(held_nothing_circles._batch_witnessed(geometry, image, selection, tool)[1]) is witnessed


def test_the_rational_check_treats_the_earlier_rim_as_covered() -> None:
    _, _, image, selection = _rim_fixture()
    tiny = Fraction(1, 1 << 80)
    # Sub-ulp dyadic points: exactly on E's outer rim (covered), and just outside it (a witness).
    assert not image.is_rational_witness(1, selection.of(1), Fraction(4) + 0 * tiny, Fraction(0))
    assert image.is_rational_witness(1, selection.of(1), Fraction(4) + tiny, Fraction(0))
    assert not image.is_rational_witness(1, selection.of(1), Fraction(4) - tiny, Fraction(0))


def test_a_witness_off_a_coarse_dyadic_path_is_not_snapped_onto_a_rim() -> None:
    # An all-integer path has dyadic scale 0; a binary64 candidate must still be checked at its own value.
    assert material_witness(RIM_TARGET, (RIM_EARLIER,), TOOL) is not None
    result = remove_nothing_circles((RIM_EARLIER, RIM_TARGET), _identity, TOOL, TIGHT, closed=False)
    assert result.union_tests == 0


# The union test's early exit runs only past FIRST_UNION_QUERY sweeps.
FAR = tuple(_circle(40.0 + 10.0 * k, 0.0, 3.0) for k in range(held_nothing_circles.FIRST_UNION_QUERY))


def test_the_early_exit_returns_a_nothing_prefix_verdict() -> None:
    # The exact duplicate comes first: the query after FIRST_UNION_QUERY sweeps proves NOTHING.
    assert union_covers_annulus(BIG, (_circle(0.0, 0.0, 3.0), *FAR, *FAR), TOOL)


def test_the_early_exit_does_not_return_a_partial_verdict() -> None:
    # Covered only by the sweep after the first query: the prefix verdict is not NOTHING, the full one is.
    assert len(FAR) == held_nothing_circles.FIRST_UNION_QUERY
    assert not union_covers_annulus(BIG, FAR, TOOL)
    assert union_covers_annulus(BIG, (*FAR, _circle(0.0, 0.0, 3.0)), TOOL)


@settings(max_examples=40, deadline=None)
@given(circles=st.lists(st.builds(_circle, _COORDINATE, _COORDINATE, _RADIUS), min_size=2, max_size=30))
def test_the_batched_selection_equals_the_per_circle_selection(circles: list[PaperCircleCandidate]) -> None:
    geometry = held_nothing_circles._arrays(circles)
    tool = float(TOOL.value)
    selection = held_nothing_circles._Selection.build(geometry, tool)
    for index in range(len(circles)):
        assert np.array_equal(selection.of(index), held_nothing_circles._overlapping_earlier(geometry, index, tool))
