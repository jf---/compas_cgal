"""Ω* connector repair: refused connectors repaired by ∂Ω* midpoint circles, convex corners over the cap refused by name.

Fixtures are declared pockets with exact Ω* (r = 1 mm) and binary64 circles
contained in Ω* exactly. The hexagon has exterior angles 53.13° and 73.74°
(all under the 80° cap), so its seam lap is repairable; the L-shape has
right-angle convex corners, where the convex-corner lemma forbids any
certificate under a cap below 90°.
"""

from __future__ import annotations

import math
from fractions import Fraction
from functools import cache

import pytest

from benchmarks.held_connector_engagement import ConnectorTier
from benchmarks.held_connector_engagement_exact import ExactConnectorAssessment
from benchmarks.held_connector_engagement_range import RangeEngagementLemma
from benchmarks.held_connector_engagement_range import SweepKind
from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_declared_boundary import DeclaredLine
from benchmarks.held_declared_boundary import PythagoreanDirection
from benchmarks.held_declared_clearance import DeclaredClearance
from benchmarks.held_exact_centre_domain import CentreDomain
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_omega_circles import OmegaCircle
from benchmarks.held_omega_circles import contain
from benchmarks.held_omega_circles import place_local
from benchmarks.held_omega_connector_engagement import lemma_inputs
from benchmarks.held_omega_connector_engagement import omega_connector_path
from benchmarks.held_omega_connector_gate import OmegaConnectorGate
from benchmarks.held_omega_connector_engagement import exact_omega_connector
from benchmarks.held_omega_connector_engagement import prior_circle
from benchmarks.held_omega_connector_repair import OmegaConnectorRepair
from benchmarks.held_omega_connector_repair import OmegaConnectorRepairBudgetError
from benchmarks.held_omega_connector_repair import OmegaConnectorRepairContractError
from benchmarks.held_omega_connector_repair import OmegaConnectorRepairDepthError
from benchmarks.held_omega_connector_repair import OmegaConvexCornerEngagementError
from benchmarks.held_omega_connector_repair import _require_paired
from benchmarks.held_omega_connector_repair import repair_omega_connectors
from benchmarks.held_omega_connectors import OmegaConnectors
from benchmarks.held_omega_connectors import connect
from benchmarks.held_omega_convex_corner import exterior_angle_exceeds
from benchmarks.held_omega_convex_corner import sharp_vertex
from benchmarks.held_omega_gouge_audit import audit_omega_path
from benchmarks.held_omega_refinement import transition_engagement
from benchmarks.held_omega_sites import OmegaLoop
from benchmarks.held_omega_sites import OmegaSite
from benchmarks.held_omega_toolpath import OmegaGeneratedToolpath
from benchmarks.held_omega_toolpath import _containment_oracle
from benchmarks.held_omega_toolpath import generate_omega
from compas_cgal.motion import EngagementCap
from compas_cgal.units import ExactMillimetre
from compas_cgal.units import ToolRadius
from tests.benchmarks.test_held_exact_centre_domain import _point
from tests.benchmarks.test_held_exact_centre_domain import _Pocket
from tests.benchmarks.test_held_exact_centre_domain import l_shape_pocket

F = Fraction
TOOL = ToolRadius.build(1.0)
CAP = EngagementCap.build(math.radians(80))
# Admissible gap (mm) between a fixture contact and ∂Ω*: covers binary64 rounding of the contacts.
GAP = ExactMillimetre(F(1, 4))
# Hexagon with a convex-arc bottom: arc (0,0) -> (10,0) about (5, 3.75), ρ = 6.25; its Ω* arc has ρ' = 5.25.
ARC_CENTRE = (5.0, 3.75)
OFFSET_ARC_RADIUS = 5.25
# Angular position (rad from straight down) of the first arc contact; the second lies 1.5 mm further along Ω*.
ARC_START = -0.3
ARC_STEP_MM = 1.5


def _hexagon(convex_arc_bottom: bool) -> DeclaredBoundary:
    """CCW hexagon, edges at 0°, 53.13°, 126.87°, 180°, 233.13°, 306.87° (lengths 10, 5, 5, 10, 5, 5); optionally a tangent convex-arc bottom."""
    pocket = _Pocket(_point(0, 0))
    pocket = pocket.arc(_point(10, 0), F(-1, 2), 0, 1) if convex_arc_bottom else pocket.line(0, 10)
    for t, quarter, length in ((F(1, 2), 0, 5), (F(-1, 2), 2, 5), (F(0), 2, 10), (F(1, 2), 2, 5), (F(-1, 2), 0, 5)):
        line = DeclaredLine.build(pocket.cursor, PythagoreanDirection.build(t, quarter), F(length), len(pocket.primitives))
        pocket.primitives.append(line)
        pocket.cursor = line.end
    return pocket.close("hexagon_arc" if convex_arc_bottom else "hexagon")


class _Fixture:
    def __init__(self, boundary: DeclaredBoundary, contacts: list[tuple[float, float]], radii: tuple[float, ...]) -> None:
        self.clearance = DeclaredClearance.build(boundary, F(1))
        self.domain = CentreDomain.build(boundary, F(1))
        self.loop = OmegaLoop.build(self.domain)
        sites = [OmegaSite.decided(self.loop, decision) for decision in self.domain.decide_contacts(contacts, GAP).contacts]
        self.circles = tuple(contain(place_local(site, radius, _provenance(k)), self.clearance, self.domain) for k, (site, radius) in enumerate(zip(sites, radii)))

    def assessment(self, circles: tuple[OmegaCircle, ...], cap: EngagementCap = CAP) -> ExactConnectorAssessment:
        return ExactConnectorAssessment.build(omega_connector_path(circles, OmegaConnectors.build(self.loop.pieces, circles)), F(1), cap.chord_ratio)

    def repair(self, cap: EngagementCap = CAP, **budget: int) -> OmegaConnectorRepair:
        return repair_omega_connectors(self.circles, self.loop, self.clearance, TOOL, cap, **budget)


def _provenance(index: int) -> tuple[GuideRunId, GuideRunStationOrdinal, ProjectionBoundarySite]:
    return GuideRunId(0), GuideRunStationOrdinal(index), ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), F(0), 1)


@cache
def line_fixture() -> _Fixture:
    """Two circles on the hexagon's bottom Ω* line (y = 1), ρ 2 -> 1, 1 mm apart: pair within 80°, connector refused on its line."""
    return _Fixture(_hexagon(False), [(4.0, 1.0), (5.0, 1.0)], (2.0, 1.0))


@cache
def arc_fixture() -> _Fixture:
    """Two circles on the hexagon's offset convex arc, ρ 2 -> 1, 1.5 mm apart: refused on its arc and entry junction."""
    angles = (ARC_START, ARC_START + ARC_STEP_MM / OFFSET_ARC_RADIUS)
    contacts = [(ARC_CENTRE[0] + OFFSET_ARC_RADIUS * math.sin(a), ARC_CENTRE[1] - OFFSET_ARC_RADIUS * math.cos(a)) for a in angles]
    return _Fixture(_hexagon(True), contacts, (2.0, 1.0))


@cache
def corner_fixture() -> _Fixture:
    """L-shape: one circle on the bottom Ω* line, one on the right one; connector 0 crosses the right-angle vertex (19, 1)."""
    return _Fixture(l_shape_pocket(), [(17.0, 1.0), (19.0, 3.0)], (1.0, 1.0))


def _assert_renumbered(result: OmegaConnectorRepair) -> None:
    """Carried-over verdicts name their connector and only prior circles; tier 1 names the circle just cut."""
    for j, verdict in enumerate(result.verdicts):
        assert verdict.connector == j
        assert all(sub.circle <= j for sub in verdict.sub_sweeps)
        if verdict.tier is ConnectorTier.TIER_1:
            assert {sub.circle for sub in verdict.sub_sweeps} == {j}


def _assert_repaired(fixture: _Fixture, result: OmegaConnectorRepair) -> None:
    fresh = fixture.assessment(result.circles)
    assert fresh.refused == ()
    assert [verdict.tier for verdict in result.verdicts] == [verdict.tier for verdict in fresh.verdicts]
    _assert_renumbered(result)
    audit = audit_omega_path(fixture.clearance, result.circles, result.connectors.connectors)
    assert not audit.gouging_circles and not audit.gouging_connectors
    for first, second in zip(result.circles, result.circles[1:]):
        assert transition_engagement(first, second, TOOL, 0).within(CAP)


def test_the_fixtures_refuse_a_connector_whose_pair_is_within_the_cap() -> None:
    for fixture, kinds in ((line_fixture(), {SweepKind.OFFSET_LINE}), (arc_fixture(), {SweepKind.OFFSET_ARC, SweepKind.CIRCLE_JUNCTION})):
        first, second = fixture.circles
        assert transition_engagement(first, second, TOOL, 0).within(CAP)
        verdict = fixture.assessment(fixture.circles).verdicts[0]
        assert verdict.tier is ConnectorTier.REFUSED
        assert {refusal.kind for refusal in verdict.refusals} == kinds


def test_the_gate_decides_every_connector_as_the_post_gate() -> None:
    fixture = line_fixture()
    circles = fixture.repair().circles
    post = fixture.assessment(circles)
    gate = OmegaConnectorGate(TOOL, CAP)
    priors = [prior_circle(circle) for circle in circles]
    count = len(circles)
    for k, expected in enumerate(post.verdicts):
        leaving, entering = circles[k], circles[(k + 1) % count]
        exact = exact_omega_connector(leaving, entering, connect(fixture.loop.pieces, leaving, entering))
        decided = gate.decide(priors[: k + 1], exact, k, seam=expected.seam)
        assert decided == expected
    refused_input = fixture.assessment(fixture.circles).verdicts[0]
    first, second = fixture.circles
    exact = exact_omega_connector(first, second, connect(fixture.loop.pieces, first, second))
    assert gate.decide([prior_circle(first)], exact, 0, seam=False) == refused_input


def test_a_refused_line_connector_is_repaired_by_one_midpoint() -> None:
    fixture = line_fixture()
    result = fixture.repair()
    _assert_repaired(fixture, result)
    line_action = result.actions[0]
    assert (line_action.connector, line_action.seam, line_action.inserted) == (0, False, 1)
    first, second = fixture.circles
    assert result.circles[0] is first and result.circles[2] is second
    midpoint = result.circles[1].site
    assert first.site.key < midpoint.key < second.site.key


def test_a_refused_arc_connector_is_repaired_and_its_junction_with_it() -> None:
    fixture = arc_fixture()
    result = fixture.repair()
    _assert_repaired(fixture, result)
    arc_action = result.actions[0]
    assert (arc_action.connector, arc_action.seam, arc_action.inserted) == (0, False, 1)
    assert {refusal.kind for refusal in arc_action.refusals} == {SweepKind.OFFSET_ARC, SweepKind.CIRCLE_JUNCTION}


def test_the_seam_is_repaired_by_appending_and_circle_zero_is_never_recut() -> None:
    fixture = line_fixture()
    result = fixture.repair()
    seam_action = result.actions[-1]
    assert seam_action.seam and seam_action.connector == 1 and seam_action.inserted > 0
    assert result.positions == (0, 2)
    assert result.circles[0] is fixture.circles[0]
    assert len(result.circles) == 3 + seam_action.inserted
    assert result.verdicts[-1].seam and result.verdicts[-1].tier is not ConnectorTier.REFUSED
    assert result.inserted == sum(action.inserted for action in result.actions)


def test_a_certified_path_is_returned_unchanged() -> None:
    fixture = line_fixture()
    repaired = fixture.repair().circles
    again = repair_omega_connectors(repaired, fixture.loop, fixture.clearance, TOOL, CAP)
    assert again.circles is repaired and again.actions == () and again.positions == tuple(range(len(repaired)))


def test_depth_exhaustion_names_the_connector_and_its_refused_sweeps() -> None:
    with pytest.raises(OmegaConnectorRepairDepthError) as raised:
        line_fixture().repair(max_depth=0)
    assert (raised.value.connector, raised.value.seam) == (0, False)
    assert {refusal.kind for refusal in raised.value.refusals} == {SweepKind.OFFSET_LINE}


def test_the_circle_budget_fails_loud() -> None:
    with pytest.raises(OmegaConnectorRepairBudgetError) as raised:
        line_fixture().repair(max_circles=0)
    assert (raised.value.connector, raised.value.seam, raised.value.budget) == (0, False, 0)


@cache
def budget_fixture() -> _Fixture:
    """Four circles on the hexagon's bottom Ω* line, ρ 2, 1, 2, 1 at x = 2..5: connectors 0..3 refused; repair inserts 1, 6, 1, then the seam lap."""
    return _Fixture(_hexagon(False), [(2.0, 1.0), (3.0, 1.0), (4.0, 1.0), (5.0, 1.0)], (2.0, 1.0, 2.0, 1.0))


# Insertions per input connector when the budget is ample, in connector order (the last is the seam lap).
BUDGET_FIXTURE_INSERTED = (1, 6, 1)


def test_the_circle_budget_is_global_across_connectors() -> None:
    fixture = budget_fixture()
    full = fixture.repair()
    assert tuple(action.inserted for action in full.actions[:3]) == BUDGET_FIXTURE_INSERTED
    # 7 circles cover connectors 0 and 1 (1 + 6) but not connector 2 as well: a per-connector budget would
    # only run out on the seam lap.
    budget = sum(BUDGET_FIXTURE_INSERTED[:2])
    with pytest.raises(OmegaConnectorRepairBudgetError) as raised:
        fixture.repair(max_circles=budget)
    assert (raised.value.connector, raised.value.seam, raised.value.budget) == (2, False, budget)


def test_verdicts_that_do_not_pair_with_the_repaired_connectors_fail_loud() -> None:
    fixture = line_fixture()
    result = fixture.repair()
    verdicts = list(result.verdicts)
    _require_paired(verdicts, result.connectors, True)
    with pytest.raises(OmegaConnectorRepairContractError):
        _require_paired(verdicts[:-1], result.connectors, True)
    with pytest.raises(OmegaConnectorRepairContractError):
        _require_paired(verdicts, result.connectors, False)
    with pytest.raises(OmegaConnectorRepairContractError):
        _require_paired([verdicts[1], verdicts[0], *verdicts[2:]], result.connectors, True)


def test_a_right_angle_corner_under_an_80_degree_cap_is_refused_by_name_before_any_midpoint(monkeypatch: pytest.MonkeyPatch) -> None:
    def no_midpoint(*_: object) -> OmegaCircle:
        raise AssertionError("a midpoint was tried at a convex corner over the cap")

    monkeypatch.setattr("benchmarks.held_omega_connector_repair.omega_midpoint_circle", no_midpoint)
    fixture = corner_fixture()
    with pytest.raises(OmegaConvexCornerEngagementError) as raised:
        fixture.repair()
    first, second = fixture.circles
    assert (raised.value.connector, raised.value.seam, raised.value.piece) == (0, False, second.site.index)


def test_the_corner_decision_compares_the_exterior_angle_with_the_cap_exactly() -> None:
    fixture = corner_fixture()
    first, second = fixture.circles
    for degrees, sharp in ((80, True), (89, True), (91, False), (100, False)):
        lemma = RangeEngagementLemma.build(F(1), EngagementCap.build(math.radians(degrees)).chord_ratio)
        found = sharp_vertex(fixture.loop, [first], second, lemma.cos_sq_half_cap, F(1))
        assert (found == second.site.index) is sharp, degrees
        incoming, outgoing = fixture.loop.pieces[first.site.index].support, fixture.loop.pieces[second.site.index].support
        assert exterior_angle_exceeds(incoming, outgoing, lemma.cos_sq_half_cap) is sharp
        assert not exterior_angle_exceeds(outgoing, incoming, lemma.cos_sq_half_cap)


def test_hexagon_corners_are_never_sharp_under_80_degrees() -> None:
    loop = line_fixture().loop
    lemma = RangeEngagementLemma.build(F(1), CAP.chord_ratio)
    pieces = loop.pieces
    turns = [(pieces[k - 1].support, pieces[k].support) for k in range(len(pieces))]
    assert all(exterior_angle_exceeds(a, b, lemma.cos_sq_half_cap) is False for a, b in turns)


@cache
def _figure5() -> tuple[OmegaGeneratedToolpath, OmegaConnectorRepair]:
    path = generate_omega("figure5_refined")
    cap = EngagementCap.build(math.radians(float(path.case.tea_cap)))
    oracle = _containment_oracle(path.declared, F(float(path.case.tool_radius.value)), path.clearance, path.generation_backend)
    return path, repair_omega_connectors(path.nothing.circles, path.loop, oracle, path.case.tool_radius, cap)


def test_figure5_on_omega_is_fully_certified_after_repair_and_gouge_free() -> None:
    path, result = _figure5()
    assert [action.inserted for action in result.actions] == [1]
    assert len(result.circles) == len(path.nothing.circles) + 1
    assert result.circles == path.circles  # the generator's connector_repair stage emits exactly this
    tool, ratio = lemma_inputs(path)
    fresh = ExactConnectorAssessment.build(omega_connector_path(result.circles, result.connectors), tool, ratio)
    assert fresh.refused == ()
    assert [verdict.tier for verdict in fresh.verdicts] == [verdict.tier for verdict in result.verdicts]
    _assert_renumbered(result)
    audit = audit_omega_path(path.clearance, result.circles, result.connectors.connectors)
    assert not audit.gouging_circles and not audit.gouging_connectors
