"""The connector post-gate of the Held toolpath benchmark: every connector, seam included, certified or a named error."""

from __future__ import annotations

import math
from typing import Any

import pytest

from benchmarks import held_toolpath_benchmark as tb
from benchmarks.held_centre_domain_supports import PieceKind
from benchmarks.held_connector_engagement import ConnectorEngagementAssessment
from benchmarks.held_connector_engagement import UncertifiedConnectorEngagementError
from benchmarks.held_figure5_boundary_path import Figure5BoundaryConnectors
from benchmarks.held_figure5_boundary_path import Figure5SeamClosure
from benchmarks.held_figure5_engagement_refinement import Figure5EngagementRefinement
from benchmarks.held_figure5_engagement_refinement import refine_figure5_engagement
from benchmarks.held_omega_toolpath import PathDomain
from benchmarks.held_omega_toolpath import generate_omega
from compas_cgal.motion import EngagementCap

pytestmark = pytest.mark.xdist_group("held_toolpath_connector_gate")

# Keeping every second refined circle opens gaps tiers 1 and 1b cannot close on
# Figure 5 (measured: 3 refused connectors), a deliberate certificate defect.
THINNING_STRIDE = 2


def test_a_real_polygon_workload_passes_the_gate_as_its_own_stage() -> None:
    path = tb.generate("figure5_refined")
    assert "connector_certificate" in path.stage_seconds
    assert path.connectors.seam is Figure5SeamClosure.CONNECTOR
    cap = EngagementCap.build(math.radians(float(path.case.tea_cap)))
    assessment = ConnectorEngagementAssessment.build(path.circles, path.connectors, path.case.tool_radius, cap)
    assert not assessment.refused
    assert assessment.verdicts[-1].seam


def test_a_deliberately_uncertified_connector_makes_generation_raise(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[dict[str, Any]] = []

    def thinned(*args: Any, **kwargs: Any) -> Figure5EngagementRefinement:
        calls.append(kwargs)
        result = refine_figure5_engagement(*args, **kwargs)
        return Figure5EngagementRefinement(result.circles[::THINNING_STRIDE], result.original_indices)

    monkeypatch.setattr(tb, "refine_figure5_engagement", thinned)
    with pytest.raises(UncertifiedConnectorEngagementError) as raised:
        tb.generate("figure5_refined")
    assert calls and calls[0].get("close_seam") is True
    assert raised.value.connectors
    # The refused set is exactly the independent assessment's on the same thinned path.
    monkeypatch.undo()
    reference = tb.generate("figure5_refined")
    circles = reference.circles[::THINNING_STRIDE]
    connectors = Figure5BoundaryConnectors.build(reference.boundary, circles)
    cap = EngagementCap.build(math.radians(float(reference.case.tea_cap)))
    assessment = ConnectorEngagementAssessment.build(circles, connectors, reference.case.tool_radius, cap)
    assert raised.value.connectors == assessment.refused


def test_boundary_spacing_repairs_the_seam(monkeypatch: pytest.MonkeyPatch) -> None:
    class Sentinel(Exception):
        """Stops the generator once the refinement call is observed."""

    seen: list[dict[str, Any]] = []

    def spy(*args: Any, **kwargs: Any) -> Any:
        seen.append(kwargs)
        raise Sentinel

    monkeypatch.setattr(tb, "refine_boundary_ordered_path", spy)
    with pytest.raises(Sentinel):
        tb.generate("figure8_crossed_skis_boundary")
    assert seen[0].get("close_seam") is True


def test_omega_connectors_without_a_sweep_type_raise_the_named_error() -> None:
    path = generate_omega("figure5_refined")
    moving = tuple(k for k, connector in enumerate(path.connectors.connectors) if connector.pieces)
    arcs = tuple(k for k in moving if any(piece.piece.kind is not PieceKind.OFFSET_LINE for piece in path.connectors.connectors[k].pieces))
    assert arcs
    with pytest.raises(tb.UncertifiableOmegaConnectorError) as raised:
        tb.require_omega_connector_certificate(path)
    assert isinstance(raised.value, UncertifiedConnectorEngagementError)
    assert raised.value.connectors == moving
    assert raised.value.arc_connectors == arcs
    assert set(raised.value.arc_connectors) | set(raised.value.one_root_line_connectors) | set(raised.value.rational_line_connectors) == set(moving)


def test_generate_on_omega_runs_the_gate() -> None:
    with pytest.raises(tb.UncertifiableOmegaConnectorError):
        tb.generate_on("figure5_refined", PathDomain.OMEGA)
