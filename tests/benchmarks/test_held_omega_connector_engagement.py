"""Contract tests of the Ω* connector-engagement adapter on a generated Ω* path (Figure 5).

The adapter (:mod:`benchmarks.held_omega_connector_engagement`) turns the
emitted Ω* path (binary64 circles, exact ∂Ω* connectors between exact sites)
into the certificate's exact connector path under OC's junction model:
``J_k`` is rational and exactly on circle ``k``, and connector ``i`` is
``J_i -> q_i``, the emitted ∂Ω* pieces, ``q_{i+1} -> J_{i+1}``. These tests
pin that model on real generated geometry, a certified and a refused
fragment, a fan and the seam, and the gate-results tool's junction check.
"""

from __future__ import annotations

import dataclasses
from fractions import Fraction

import pytest

from benchmarks.held_connector_engagement import ConnectorTier
from benchmarks.held_connector_engagement_arcs import ArcSweep
from benchmarks.held_connector_engagement_exact import ExactConnectorAssessment
from benchmarks.held_connector_engagement_exact import ExactConnectorCertificate
from benchmarks.held_connector_engagement_exact import ExactConnectorPath
from benchmarks.held_connector_engagement_exact import UncertifiedExactConnectorError
from benchmarks.held_connector_engagement_lines import LineSweep
from benchmarks.held_connector_engagement_range import RangeEngagementLemma
from benchmarks.held_connector_engagement_range import SweepKind
from benchmarks.held_centre_domain_supports import OffsetLineSupport
from benchmarks.held_centre_domain_supports import PieceKind
from benchmarks.held_figure5_boundary_path import Figure5SeamClosure
from benchmarks.held_omega_circles import OmegaCircle
from benchmarks.held_omega_connector_engagement import assess_omega_connectors
from benchmarks.held_omega_connector_engagement import circle_junction_point
from benchmarks.held_omega_connector_engagement import lemma_inputs
from benchmarks.held_omega_connector_engagement import omega_connector_path
from benchmarks.held_omega_connector_engagement import require_exact_omega_connector_certificate
from benchmarks.held_omega_connectors import OmegaConnector
from benchmarks.held_omega_connectors import OmegaConnectors
from benchmarks.held_omega_connectors import connect
from benchmarks.held_omega_toolpath import OmegaGeneratedToolpath
from benchmarks.held_omega_toolpath import generate_omega
from benchmarks.held_one_root import OneRootPoint2
from compas_cgal.units import Millimetre
from tools.held_omega_connector_gate_results import junction_model_check
from tools.held_omega_connector_gate_results import render

_ARC_KIND = {PieceKind.OFFSET_ARC: SweepKind.OFFSET_ARC, PieceKind.JUNCTION_ARC: SweepKind.JUNCTION_ARC}


@pytest.fixture(scope="module")
def figure5() -> OmegaGeneratedToolpath:
    return generate_omega("figure5_refined")


@pytest.fixture(scope="module")
def figure5_exact(figure5: OmegaGeneratedToolpath) -> ExactConnectorPath:
    return omega_connector_path(figure5.circles, figure5.connectors)


@pytest.fixture(scope="module")
def figure5_assessment(figure5: OmegaGeneratedToolpath) -> ExactConnectorAssessment:
    return assess_omega_connectors(figure5)


@pytest.fixture(scope="module")
def figure5_unrepaired(figure5: OmegaGeneratedToolpath) -> OmegaGeneratedToolpath:
    """Figure 5 as emitted before connector repair (the NOTHING pass's circles), which still refuses a connector to pin."""
    circles = figure5.nothing.circles
    connectors = OmegaConnectors.build(figure5.loop.pieces, circles)
    tool, ratio = lemma_inputs(figure5)
    verdicts = ExactConnectorAssessment.build(omega_connector_path(circles, connectors), tool, ratio).verdicts
    repair = dataclasses.replace(figure5.repair, circles=circles, positions=tuple(range(len(circles))), connectors=connectors, verdicts=verdicts, actions=())
    return dataclasses.replace(figure5, repair=repair)


@pytest.fixture(scope="module")
def figure5_unrepaired_assessment(figure5_unrepaired: OmegaGeneratedToolpath) -> ExactConnectorAssessment:
    return assess_omega_connectors(figure5_unrepaired)


def _junction(circle: OmegaCircle) -> OneRootPoint2:
    return circle_junction_point((float(circle.center.x), float(circle.center.y)), float(circle.radius.value), circle.site.point())


def _exactly_on(circle: OmegaCircle, point: OneRootPoint2) -> bool:
    cx, cy, radius = Fraction(float(circle.center.x)), Fraction(float(circle.center.y)), Fraction(float(circle.radius.value))
    return point.x.is_rational and point.y.is_rational and (point.x.rational - cx) ** 2 + (point.y.rational - cy) ** 2 == radius * radius


def _expected_sweeps(leaving: OmegaCircle, entering: OmegaCircle, connector: OmegaConnector) -> list[tuple[SweepKind, OneRootPoint2, OneRootPoint2]]:
    """OC's junction model, spelled out independently of the adapter: (kind, start, end) per sweep."""
    first, last = _junction(leaving), _junction(entering)
    expected = [(SweepKind.CIRCLE_JUNCTION, first, connector.start)] if first != connector.start else []
    for piece in connector.pieces:
        kind = SweepKind.OFFSET_LINE if isinstance(piece.piece.support, OffsetLineSupport) else _ARC_KIND[piece.piece.kind]
        expected.append((kind, piece.start, piece.end))
    if connector.end != last:
        expected.append((SweepKind.CIRCLE_JUNCTION, connector.end, last))
    return expected


def _fragment(path: OmegaGeneratedToolpath, index: int) -> tuple[tuple[OmegaCircle, ...], OmegaConnectors]:
    """Circles ``index, index + 1`` and the one emitted connector between them, as an open two-circle path."""
    circles = (path.circles[index], path.circles[index + 1])
    connector = path.connectors.connectors[index]
    return circles, OmegaConnectors((connector,), Figure5SeamClosure.CONNECTOR, Millimetre(0.0))


def test_every_junction_is_rational_exactly_on_its_binary64_circle(figure5: OmegaGeneratedToolpath, figure5_exact: ExactConnectorPath) -> None:
    assert len(figure5_exact.junctions) == len(figure5.circles)
    for circle, junction in zip(figure5.circles, figure5_exact.junctions):
        assert junction == _junction(circle)
        assert _exactly_on(circle, junction)


def test_every_connector_is_its_junction_lines_around_the_emitted_pieces(figure5: OmegaGeneratedToolpath, figure5_exact: ExactConnectorPath) -> None:
    circles, emitted = figure5.circles, figure5.connectors.connectors
    assert figure5_exact.closed and len(figure5_exact.connectors) == len(emitted) == len(circles)
    for i, (connector, modelled) in enumerate(zip(emitted, figure5_exact.connectors)):
        expected = _expected_sweeps(circles[i], circles[(i + 1) % len(circles)], connector)
        assert [(sweep.kind, sweep.start, sweep.end) for sweep in modelled.sweeps] == expected, f"connector {i}"
        pieces = iter(connector.pieces)
        for sweep in modelled.sweeps:
            if sweep.kind is SweepKind.CIRCLE_JUNCTION:
                assert isinstance(sweep, LineSweep)
                continue
            support = next(pieces).piece.support
            if isinstance(support, OffsetLineSupport):
                assert isinstance(sweep, LineSweep)
            else:
                assert isinstance(sweep, ArcSweep)
                assert (sweep.centre, sweep.radius, sweep.sense) == (support.centre, support.radius, support.sense)


def test_the_seam_runs_back_to_the_first_junction(figure5: OmegaGeneratedToolpath, figure5_exact: ExactConnectorPath) -> None:
    assert figure5.connectors.seam is Figure5SeamClosure.CONNECTOR
    seam = figure5_exact.connectors[-1]
    assert seam.sweeps[-1].end == figure5_exact.junctions[0]
    assert seam.sweeps[0].start == figure5_exact.junctions[-1]


def test_a_certified_fragment_is_tier_1_with_its_junction_lines(figure5: OmegaGeneratedToolpath, figure5_assessment: ExactConnectorAssessment) -> None:
    def chosen(i: int) -> bool:
        connector = figure5.connectors.connectors[i]
        has_arc = any(not isinstance(piece.piece.support, OffsetLineSupport) for piece in connector.pieces)
        return figure5_assessment.verdicts[i].tier is ConnectorTier.TIER_1 and has_arc and len(_expected_sweeps(figure5.circles[i], figure5.circles[i + 1], connector)) >= 3

    index = next(i for i in range(len(figure5.circles) - 1) if chosen(i))
    circles, connectors = _fragment(figure5, index)
    exact = omega_connector_path(circles, connectors)
    assert not exact.closed
    kinds = [sweep.kind for sweep in exact.connectors[0].sweeps]
    assert kinds[0] is SweepKind.CIRCLE_JUNCTION and kinds[-1] is SweepKind.CIRCLE_JUNCTION
    tool, ratio = lemma_inputs(figure5)
    certificate = ExactConnectorCertificate.build(exact, tool, ratio)
    (verdict,) = certificate.verdicts
    assert verdict.tier is ConnectorTier.TIER_1
    assert [sub.circle for sub in verdict.sub_sweeps] == [0] * len(kinds)


def test_a_refused_fragment_is_named_and_its_tier_1_bound_exceeds_the_cap(
    figure5_unrepaired: OmegaGeneratedToolpath, figure5_unrepaired_assessment: ExactConnectorAssessment
) -> None:
    figure5, figure5_assessment = figure5_unrepaired, figure5_unrepaired_assessment
    assert figure5_assessment.refused, "Figure 5 on Ω* before repair has a refused connector to pin"
    index = figure5_assessment.refused[0]
    circles, connectors = _fragment(figure5, index)
    exact = omega_connector_path(circles, connectors)
    tool, ratio = lemma_inputs(figure5)
    with pytest.raises(UncertifiedExactConnectorError) as raised:
        ExactConnectorCertificate.build(exact, tool, ratio)
    fragment = {(refusal.sweep, refusal.kind) for refusal in raised.value.refusals}
    full = {(refusal.sweep, refusal.kind) for refusal in figure5_assessment.verdicts[index].refusals}
    assert raised.value.connectors == (0,)
    assert full <= fragment  # fewer priors can only refuse more
    lemma = RangeEngagementLemma.build(tool, ratio)
    cap = float(figure5.case.tea_cap)
    for sweep, _ in full:
        assert lemma.bound_deg(exact.connectors[0].sweeps[sweep], exact.circles[0]) > cap  # reporting cross-check of "not proved"


def test_require_raises_with_exactly_the_assessed_refusals(figure5_unrepaired: OmegaGeneratedToolpath, figure5_unrepaired_assessment: ExactConnectorAssessment) -> None:
    with pytest.raises(UncertifiedExactConnectorError) as raised:
        require_exact_omega_connector_certificate(figure5_unrepaired)
    assert raised.value.refusals == figure5_unrepaired_assessment.refusals
    assert raised.value.connectors == figure5_unrepaired_assessment.refused


def test_a_fan_is_two_junction_lines_out_and_back(figure5: OmegaGeneratedToolpath) -> None:
    circle = next(c for c in figure5.circles if _junction(c) != c.site.point())
    twin = dataclasses.replace(circle)
    fan = connect(figure5.loop.pieces, circle, twin)
    assert fan.within_fan and fan.pieces == ()
    exact = omega_connector_path((circle, twin), OmegaConnectors((fan,), Figure5SeamClosure.CONNECTOR, Millimetre(0.0)))
    out, back = exact.connectors[0].sweeps
    site, junction = circle.site.point(), _junction(circle)
    assert (out.kind, out.start, out.end) == (SweepKind.CIRCLE_JUNCTION, junction, site)
    assert (back.kind, back.start, back.end) == (SweepKind.CIRCLE_JUNCTION, site, junction)
    tool, ratio = lemma_inputs(figure5)
    assert ExactConnectorAssessment.build(exact, tool, ratio).verdicts[0].tier is ConnectorTier.TIER_1


def test_junction_check_passes_the_adapter_and_flags_a_swapped_connector(figure5: OmegaGeneratedToolpath, figure5_exact: ExactConnectorPath) -> None:
    report = junction_model_check(figure5, figure5_exact)
    assert report["structure_mismatches"] == []
    assert report["junctions_exactly_on_circle"] == len(figure5.circles)
    swapped = list(figure5_exact.connectors)
    swapped[3], swapped[4] = swapped[4], swapped[3]
    broken = ExactConnectorPath(figure5_exact.circles, tuple(swapped), figure5_exact.junctions)
    assert junction_model_check(figure5, broken)["structure_mismatches"] == [3, 4]


def test_render_writes_one_row_per_workload_in_each_table() -> None:
    gaps = {"circles": 1, "max_junction_to_site_mm": 1e-14, "max_contact_to_site_mm": 1e-14, "max_junction_to_contact_mm": 2e-14}
    row = {
        "workload": "w",
        "cap_deg": 80.0,
        "circles": 3,
        "connectors": 3,
        "certified": 2,
        "tiers": {"stationary": 0, "tier 1": 2, "tier 1b": 0, "refused": 1},
        "seam_tier": "tier 1",
        "refusals_by_kind": {"offset_arc": 1},
        "refusals_by_reason": {"not proved at the bisection depth": 1},
        "tier_1_bound": {"refused_sweeps_deg": [81.0, 81.0], "tier_1_connectors_worst_sweep_deg": None},
        "gate_seconds": 1.0,
        "generation_seconds": 2.0,
        "gate_over_generation": 0.5,
        "load_1min_before_after": [10.0, 11.0],
        "junction_model": {
            "structure_mismatches": [],
            "junctions_exactly_on_circle": 3,
            "junction_equals_site": 0,
            "contact_point_exactly_on_circle": 0,
            "junction_lines": 6,
            "fan_connectors": 0,
            "gaps": {"free": gaps, "limited": {**gaps, "circles": 0}},
        },
    }
    record = {"build": {"recorded_at": "t", "commit": "0123456789ab", "dirty": False}, "machine": {"cpu_brand": "cpu"}, "workloads": [row]}
    text = render(record)
    assert text.count("| `w` |") == 3
    assert "| `w` | 2 / 3 | 2 / 0 / 1 | tier 1 | offset_arc 1 |" in text
    assert "81.00–81.00 | — |" in text
