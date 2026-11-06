"""The exact NOTHING pass inside select-certify-repair: every certificate run, and the repair, see the final circle set.

Fixture: the square / fan of ``tests/benchmarks/test_held_coverage_repair.py``
(centre domain (0,0)-(6,6), r = 1, 180-degree cap), built here from the same
dyadic values.
"""

from __future__ import annotations

import math
from fractions import Fraction

import pytest

from benchmarks.held_coverage_repair import certify_boundary_ordered_selection
from benchmarks.held_coverage_repair import repair_coverage
from benchmarks.held_coverage_repair_cover import contains_exactly
from benchmarks.held_coverage_repair_invariant_errors import CertifiedNothingMismatchError
from benchmarks.held_coverage_repair_types import BoundaryOrderedEmission
from benchmarks.held_coverage_repair_types import CertifiedBoundaryOrderedPath
from benchmarks.held_coverage_repair_types import EngagementRefinementSettings
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_engagement_refinement import TracedEngagementRefinement
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_nothing_circles import NothingCirclePass
from compas_cgal.coverage_certificate.certificate import CoverageCertificate
from compas_cgal.coverage_certificate.sweeps import BoundaryConnectorSweep
from compas_cgal.coverage_certificate.sweeps import FullCircleSweep
from compas_cgal.coverage_certificate.target import CentreDomain
from compas_cgal.motion import EngagementCap
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY
from tests.benchmarks.held_part_fixture import part_around

SQUARE = tuple(Point2[WorldXY].build(x, y) for x, y in ((0, 0), (6, 0), (6, 6), (0, 6)))
TOOL = ToolRadius.build(1)
FAN_CAP = EngagementCap.build(math.pi)  # the fan's engagements are exactly 180 degrees
SETTINGS = EngagementRefinementSettings.build(TOOL, FAN_CAP, part=part_around(SQUARE, TOOL))
HALF = Fraction(1, 2)
# Index of the emitted circle duplicated by the defect injection: an interior
# circle, so the duplicate's successor pair is an emitted (certified) pair.
DUPLICATED = 1


def _circle(ordinal: int, segment: int, contact: tuple[float, float], centre: tuple[float, float], radius: float, run: int = 0) -> Figure5CounterclockwiseCircle:
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(segment), HALF, 4)
    return Figure5CounterclockwiseCircle(
        GuideRunId(run), GuideRunStationOrdinal(ordinal), site, site, Point2[WorldXY].build(*centre), Point2[WorldXY].build(*contact), GuideRadius.build(radius)
    )


FIXTURE = (
    _circle(0, 0, (3, 0), (3, 1), 1),
    _circle(1, 0, (3, 0), (3, 2), 2),
    _circle(2, 0, (3, 0), (3, 3), 3, run=1),
    _circle(3, 2, (3, 6), (3, 5), 1),
)
RUNS = ((GuideRunId(0),), (GuideRunId(0),), (GuideRunId(1),), (GuideRunId(0),))


def _certified() -> CertifiedBoundaryOrderedPath:
    return certify_boundary_ordered_selection(FIXTURE, RUNS, SQUARE, SETTINGS, CoverageCertificate.build)


def _with_duplicate(emission: BoundaryOrderedEmission) -> BoundaryOrderedEmission:
    """The emission with an exact duplicate of circle ``DUPLICATED`` inserted right after it."""
    refinement = emission.refinement
    circles = (*refinement.circles[: DUPLICATED + 1], refinement.circles[DUPLICATED], *refinement.circles[DUPLICATED + 1 :])
    indices = tuple(k if k <= DUPLICATED else k + 1 for k in refinement.original_indices)
    return BoundaryOrderedEmission.build(emission.selected, TracedEngagementRefinement(circles, indices, refinement.trace))


def test_the_certified_path_emits_the_pass_result_and_exposes_it_as_circles() -> None:
    certified = _certified()
    assert isinstance(certified.nothing, NothingCirclePass)
    assert certified.circles == certified.nothing.circles
    assert certified.ledger.emitted_after == len(certified.circles)


def test_every_certificate_run_sees_the_pruned_circles_and_the_duplicate_is_not_emitted() -> None:
    reference = _certified()
    seen: list[int] = []

    def recording(domain: CentreDomain, tool: ToolRadius, cuts: tuple[FullCircleSweep, ...], connectors: tuple[BoundaryConnectorSweep, ...]) -> CoverageCertificate:
        seen.append(len(cuts))
        return CoverageCertificate.build(domain, tool, cuts, connectors)

    injected = _with_duplicate(reference.emission)
    certified = repair_coverage(reference.sources, SQUARE, SETTINGS, injected, recording)
    shifted = tuple(k + 1 if k > DUPLICATED else k for k in reference.nothing.removed)
    assert certified.nothing.removed == tuple(sorted((DUPLICATED + 1, *shifted)))
    assert certified.circles == reference.circles
    assert seen == [len(reference.circles)]
    assert (certified.ledger.emitted_before, certified.ledger.emitted_after) == (len(reference.circles), len(reference.circles))
    assert certified.emission == injected, "the repair's working emission keeps every refined circle"


@pytest.mark.parametrize("forge", ("kept_mismatch", "circle_zero_removed", "kept_for_cap_removed"))
def test_build_rejects_a_pass_result_that_is_not_the_emission_minus_its_removed_circles(forge: str) -> None:
    certified = _certified()
    refined, nothing = certified.emission.refinement.circles, certified.nothing
    removed = nothing.removed
    forged = {
        "kept_mismatch": NothingCirclePass(certified.circles[1:], removed, (), 0, 0),
        "circle_zero_removed": NothingCirclePass(tuple(c for k, c in enumerate(refined) if k not in {0, *removed}), (0, *removed), (), 0, 0),
        "kept_for_cap_removed": NothingCirclePass(certified.circles, removed, (removed[0],), 0, 0),
    }[forge]
    with pytest.raises(CertifiedNothingMismatchError):
        CertifiedBoundaryOrderedPath.build(certified.sources, certified.emission, forged, certified.certificate, certified.ledger)


def test_every_repair_witness_lies_in_a_kept_circle_of_the_emitted_path() -> None:
    """No witness lies in a removed annulus: each removed annulus is in the union of the kept cuts before it."""
    certified = _certified()
    assert certified.ledger.insertions
    for insertion in certified.ledger.insertions:
        point = insertion.witness.uncovered_point
        assert any(contains_exactly(circle, TOOL, point) for circle in certified.circles)
