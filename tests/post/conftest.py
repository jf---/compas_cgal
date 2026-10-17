"""Shared fixtures for the test-cut post: the real rect_12x8 120-degree pipeline output and a default setup."""

from fractions import Fraction

import pytest

from benchmarks.gate import gate_pocket
from benchmarks.gate import pipeline
from benchmarks.pathmetrics import CLEARANCE_Z_TOOL_DIAMETERS
from benchmarks.spec import PocketSpec
from compas_cgal.pipeline.result import PipelineResult
from compas_cgal.post.machine import PostSetup
from compas_cgal.post.testcut_defaults import build_testcut_defaults
from compas_cgal.post.mapping import SourcePath
from compas_cgal.units import ClearanceZ
from compas_cgal.units import CutZ
from compas_cgal.units import ToolRadius
from tests.post.support import ACCEPTED_PLUNGES
from tests.post.support import TEST_ENVELOPE


@pytest.fixture(scope="session")
def rect_spec() -> PocketSpec:
    return gate_pocket("rect_12x8")


@pytest.fixture(scope="session")
def rect_result(rect_spec: PocketSpec) -> PipelineResult:
    result = pipeline(rect_spec)
    assert isinstance(result, PipelineResult)
    return result


@pytest.fixture(scope="session")
def rect_source(rect_spec: PocketSpec) -> SourcePath:
    return SourcePath(
        cut=CutZ.build(0.0),
        clearance=ClearanceZ.build(CLEARANCE_Z_TOOL_DIAMETERS * rect_spec.tool_diameter),
        tool=ToolRadius.build(rect_spec.tool_radius),
    )


@pytest.fixture(scope="session")
def rect_min_corner(rect_spec: PocketSpec) -> tuple[Fraction, Fraction]:
    points = rect_spec.polygon.points
    return min(Fraction(float(p[0])) for p in points), min(Fraction(float(p[1])) for p in points)


@pytest.fixture(scope="session")
def setup(rect_min_corner: tuple[Fraction, Fraction]) -> PostSetup:
    return build_testcut_defaults(envelope=TEST_ENVELOPE, pocket_min_corner=rect_min_corner, accepted_plunges=ACCEPTED_PLUNGES)
