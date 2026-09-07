"""Small emitted-diagnostic contracts, including native corner events."""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmarks import held_native_boundary_draft as draft
from compas_cgal import _circle_geometry_2
from compas_cgal import _coverage_2
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


def _rectangle(monkeypatch: pytest.MonkeyPatch) -> None:
    case = SimpleNamespace(
        name="rectangle",
        projection=SimpleNamespace(points=tuple(Point2[WorldXY].build(x, y) for x, y in ((0, 0), (12, 0), (12, 8), (0, 8)))),
        tool_radius=ToolRadius.build(1.0),
    )
    monkeypatch.setattr(draft, "load_held_reference_case", lambda _: case)


def test_rectangle_preserves_stationary_corners_and_closes_native_transitions(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    _rectangle(monkeypatch)
    output = tmp_path / "rectangle.png"
    draft.render_native_boundary_draft("rectangle", output=output)
    report = json.loads(output.with_suffix(".json").read_text())
    assert report["complete_construction"] is True
    assert report["offset_primitives"] == 4
    assert report["sampled_events"] == 8
    assert report["stationary_events"] == 4
    assert report["transition_pieces"] == 8
    assert output.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")


def test_native_failure_publishes_partial_diagnostic_then_raises(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    _rectangle(monkeypatch)
    original = _coverage_2.boundary_circle_at_contact
    calls = 0

    def fail_third(
        owner: _circle_geometry_2.BoundaryNormalCircle2, contact: _coverage_2.WorldXYBoundaryPointMm, radius: float
    ) -> _circle_geometry_2.BoundaryNormalCircleProposal2:
        nonlocal calls
        calls += 1
        if calls == 3:
            raise _coverage_2.BoundaryContactConstructionError("deliberate third-contact failure")
        return original(owner, contact, radius)

    monkeypatch.setattr(_coverage_2, "boundary_circle_at_contact", fail_third)
    output = tmp_path / "partial.png"
    with pytest.raises(_coverage_2.BoundaryContactConstructionError, match="third-contact"):
        draft.render_native_boundary_draft("rectangle", output=output)
    report = json.loads(output.with_suffix(".json").read_text())
    assert report["complete_construction"] is False
    assert report["sampled_events"] == 2
    assert report["failure"]["type"] == "BoundaryContactConstructionError"
    assert output.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
