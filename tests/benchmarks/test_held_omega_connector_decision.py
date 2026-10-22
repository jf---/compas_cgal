"""One decision, one junction model: the public connector decision and the Ω* adapter's single connector builder.

``decide_connector`` is the certificate module's per-connector decision
(tier 1 with the circle just cut, else tier 1b over the priors); the
assessment and the repair gate both call it, so they cannot drift apart.
The Ω* adapter builds every exact connector (junction lines included) with
one function, which the post-gate path and the repair gate share.
"""

from __future__ import annotations

from fractions import Fraction

import benchmarks.held_omega_connector_gate as gate_module
from benchmarks.held_connector_engagement_exact import CandidateIndex
from benchmarks.held_connector_engagement_exact import ExactConnectorAssessment
from benchmarks.held_connector_engagement_exact import decide_connector
from benchmarks.held_connector_engagement_range import RangeEngagementLemma
from benchmarks.held_omega_connector_engagement import exact_omega_connector
from benchmarks.held_omega_connector_engagement import junction_point
from benchmarks.held_omega_connector_engagement import omega_connector_path
from benchmarks.held_omega_connector_engagement import prior_circle
from benchmarks.held_omega_connector_gate import OmegaConnectorGate
from benchmarks.held_omega_connectors import OmegaConnectors
from tests.benchmarks.test_held_omega_connector_repair import CAP
from tests.benchmarks.test_held_omega_connector_repair import TOOL
from tests.benchmarks.test_held_omega_connector_repair import arc_fixture
from tests.benchmarks.test_held_omega_connector_repair import line_fixture


def _paths():
    """Refused inputs (line and arc) and a repaired multi-circle path with a seam."""
    for fixture in (line_fixture(), arc_fixture()):
        yield fixture, fixture.circles
    fixture = line_fixture()
    yield fixture, fixture.repair().circles


def test_decide_connector_is_the_assessment_verdict_on_all_circles_and_on_the_priors_alone() -> None:
    lemma = RangeEngagementLemma.build(Fraction(1), CAP.chord_ratio)
    for fixture, circles in _paths():
        path = omega_connector_path(circles, OmegaConnectors.build(fixture.loop.pieces, circles))
        assessment = ExactConnectorAssessment.build(path, Fraction(1), CAP.chord_ratio)
        index = CandidateIndex(path.circles, lemma.tool_radius)
        for k, (connector, expected) in enumerate(zip(path.connectors, assessment.verdicts)):
            assert decide_connector(lemma, path.circles, connector, k, last=k, seam=expected.seam, index=index) == expected
            assert decide_connector(lemma, path.circles[: k + 1], connector, k, last=k, seam=expected.seam) == expected


def test_the_repair_gate_is_the_public_decision() -> None:
    fixture = line_fixture()
    circles = fixture.repair().circles
    path = omega_connector_path(circles, OmegaConnectors.build(fixture.loop.pieces, circles))
    gate = OmegaConnectorGate(TOOL, CAP)
    for k, connector in enumerate(path.connectors):
        seam = k == len(circles) - 1
        expected = decide_connector(gate.lemma, path.circles[: k + 1], connector, k, last=k, seam=seam)
        assert gate.decide(list(path.circles[: k + 1]), connector, k, seam=seam) == expected


def test_the_post_gate_path_is_built_by_the_single_connector_builder() -> None:
    for fixture, circles in _paths():
        connectors = OmegaConnectors.build(fixture.loop.pieces, circles)
        path = omega_connector_path(circles, connectors)
        count = len(circles)
        assert path.junctions == tuple(junction_point(circle) for circle in circles)
        assert path.circles == tuple(prior_circle(circle) for circle in circles)
        for i, connector in enumerate(connectors.connectors):
            assert path.connectors[i] == exact_omega_connector(circles[i], circles[(i + 1) % count], connector)


def test_the_repair_gate_has_no_junction_model_of_its_own() -> None:
    for name in ("_junction", "junction_point", "exact_omega_connector", "prior_circle", "circle_junction_point"):
        assert name not in vars(gate_module), name
