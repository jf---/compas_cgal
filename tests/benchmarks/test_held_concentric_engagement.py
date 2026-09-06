"""Concentric growth must not bypass the predecessor engagement model."""

import math

import pytest

from benchmarks.held_standard_placement import PaperCircleCandidate
from benchmarks.held_standard_placement import maximum_predecessor_engagement
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


def _circle(radius: float) -> PaperCircleCandidate:
    return PaperCircleCandidate.build(
        center=Point2[WorldXY].build(0.0, 0.0),
        guide_radius=GuideRadius.build(radius),
        contact_point=Point2[WorldXY].build(radius, 0.0),
    )


@pytest.mark.parametrize(
    ("previous_radius", "next_radius", "expected"),
    (
        (1.0, 1.0, 0.0),
        (2.0, 1.0, 0.0),
        (1.0, 1.5, math.acos(0.25)),
        (1.0, 3.0, math.pi),
        (1.0, 4.0, math.pi),
    ),
)
def test_concentric_engagement_respects_successor_radius(previous_radius: float, next_radius: float, expected: float) -> None:
    # Section 2.3 uses the previous swept disk. For Rprev=2, rho=1.5,
    # tool=1, the cosine-law intersection has cosine (4-2.25-1)/3=1/4.
    actual = maximum_predecessor_engagement(_circle(previous_radius), _circle(next_radius), ToolRadius.build(1.0))
    assert float(actual) == pytest.approx(expected)


@pytest.mark.parametrize("offset", (0.5, 1.0))
def test_displaced_successor_inside_predecessor_swept_disk_has_zero_engagement(offset: float) -> None:
    successor = PaperCircleCandidate.build(
        center=Point2[WorldXY].build(offset, 0.0),
        guide_radius=GuideRadius.build(1.0),
        contact_point=Point2[WorldXY].build(offset + 1.0, 0.0),
    )
    actual = maximum_predecessor_engagement(_circle(2.0), successor, ToolRadius.build(1.0))
    assert float(actual) == 0.0
