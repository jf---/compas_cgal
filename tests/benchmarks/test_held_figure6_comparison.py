from __future__ import annotations

import math

import pytest

import benchmarks.held_figure6_comparison as comparison_module
from benchmarks.held_figure6_comparison import COMPLIANCE_UNAUDITED
from benchmarks.held_figure6_comparison import REPOSITORY_FIGURE6_CAPS
from benchmarks.held_figure6_comparison import InvalidFigure6RepositoryMeasurementError
from benchmarks.held_figure6_comparison import RepositoryFigure6Point
from benchmarks.held_figure6_comparison import measure_repository_figure6
from benchmarks.held_figure6_publisher import FIGURE6_PUBLISHER_LABELS
from benchmarks.held_figure6_publisher import InvalidFigure6PublisherEvidenceError
from benchmarks.held_figure6_publisher import PublisherFigure6Axes
from benchmarks.held_figure6_publisher import PublisherFigure6Point
from benchmarks.held_figure6_publisher import PublisherFigure6Series
from benchmarks.held_figure6_publisher import PublisherGraphicalPathLength
from benchmarks.held_figure6_publisher import load_figure6_publisher_evidence
from benchmarks.tools.held_figure6_same_axes import build_figure
from benchmarks.tools.held_figure6_same_axes import build_literal_overlay
from benchmarks.units import Degrees
from compas_cgal.adaptive.units import Millimetre


def _axes() -> PublisherFigure6Axes:
    return PublisherFigure6Axes.build(
        x_min_deg=0.0,
        x_max_deg=200.0,
        x_tick_deg=20.0,
        y_min=10.0,
        y_max=100_000.0,
        pixel_left=329,
        pixel_right=1140,
        pixel_top=212,
        pixel_bottom=710,
    )


def test_publisher_axis_maps_pixel_frame_to_labeled_log_axes() -> None:
    axes = _axes()
    assert axes.point_from_pixel(329, 212) == PublisherFigure6Point(Degrees(0.0), PublisherGraphicalPathLength(100_000.0))
    assert axes.point_from_pixel(1140, 710) == PublisherFigure6Point(Degrees(200.0), PublisherGraphicalPathLength(10.0))


@pytest.mark.parametrize(
    ("field", "value"),
    (("x_max_deg", math.inf), ("x_tick_deg", 0.0), ("y_min", -1.0), ("pixel_right", 329)),
)
def test_publisher_axes_reject_invalid_bounds(field: str, value: float) -> None:
    values: dict[str, float | int] = {
        "x_min_deg": 0.0,
        "x_max_deg": 200.0,
        "x_tick_deg": 20.0,
        "y_min": 10.0,
        "y_max": 100_000.0,
        "pixel_left": 329,
        "pixel_right": 1140,
        "pixel_top": 212,
        "pixel_bottom": 710,
    }
    values[field] = value
    with pytest.raises(InvalidFigure6PublisherEvidenceError):
        PublisherFigure6Axes.build(**values)


def test_publisher_series_requires_expected_label_and_strict_x_order() -> None:
    points = (
        PublisherFigure6Point.build(20.0, 1000.0, axes=_axes()),
        PublisherFigure6Point.build(20.0, 900.0, axes=_axes()),
    )
    with pytest.raises(InvalidFigure6PublisherEvidenceError, match="strictly increasing"):
        PublisherFigure6Series.build("standard", points)
    with pytest.raises(InvalidFigure6PublisherEvidenceError, match="label"):
        PublisherFigure6Series.build("unknown", points[:1])


def test_tracked_publisher_evidence_has_three_distinct_ordered_series() -> None:
    evidence = load_figure6_publisher_evidence()
    assert tuple(series.label for series in evidence.series) == FIGURE6_PUBLISHER_LABELS
    assert all(len(series.points) >= 40 for series in evidence.series)
    assert evidence.provenance == "publisher Figure 6 pixel-digitized graphical observations"
    assert evidence.length_unit == "publisher graphical path-length unit"


def test_repository_adapter_uses_real_requested_caps_without_audit(monkeypatch: pytest.MonkeyPatch) -> None:
    spec = object()
    paths = {float(cap): object() for cap in REPOSITORY_FIGURE6_CAPS}
    events: list[tuple[str, object]] = []
    monkeypatch.setattr(comparison_module, "reference_pocket", lambda: spec)

    def controlled(received_spec: object, cap: float) -> object:
        assert received_spec is spec
        events.append(("generate", cap))
        return paths[cap]

    def length(path: object) -> float:
        events.append(("length", path))
        return 1000.0 + 10.0 * list(paths.values()).index(path)

    monkeypatch.setattr(comparison_module, "controlled_path", controlled)
    monkeypatch.setattr(comparison_module, "path_length", length)
    result = measure_repository_figure6()

    assert tuple(point.requested_cap_deg for point in result.points) == REPOSITORY_FIGURE6_CAPS
    assert all(point.compliance == COMPLIANCE_UNAUDITED for point in result.points)
    assert result.constant_spacing_available is False
    assert events == [
        ("generate", 80.0),
        ("length", paths[80.0]),
        ("generate", 120.0),
        ("length", paths[120.0]),
        ("generate", 160.0),
        ("length", paths[160.0]),
    ]


def test_repository_point_rejects_invalid_length() -> None:
    with pytest.raises(InvalidFigure6RepositoryMeasurementError):
        RepositoryFigure6Point.build(Degrees(80.0), Millimetre(math.nan))


def test_render_contract_has_exact_axes_and_truthful_labels() -> None:
    publisher = load_figure6_publisher_evidence()
    repository = comparison_module.RepositoryFigure6Comparison.build(
        (
            RepositoryFigure6Point.build(Degrees(80.0), Millimetre(1000.0)),
            RepositoryFigure6Point.build(Degrees(120.0), Millimetre(800.0)),
            RepositoryFigure6Point.build(Degrees(160.0), Millimetre(700.0)),
        )
    )
    figure = build_figure(publisher, repository)
    axis = figure.axes[0]
    assert axis.get_xlim() == (0.0, 200.0)
    assert axis.get_ylim() == (10.0, 100_000.0)
    assert axis.get_yscale() == "log"
    assert tuple(axis.get_xticks()) == tuple(float(value) for value in range(0, 201, 20))
    labels = tuple(line.get_label() for line in axis.lines)
    assert labels == (
        "publisher standard (digitized)",
        "publisher MATHSM (digitized)",
        "publisher contour (digitized)",
        "repository measured requested-cap; compliance unaudited",
    )
    caption = " ".join(text.get_text() for text in figure.texts).lower()
    assert "constant-spacing repository curve unavailable" in caption
    assert "no unpublished numeric parity" in caption
    assert "implementation identity" not in caption.lower()


def test_literal_overlay_places_repository_points_on_publisher_pixels() -> None:
    publisher = load_figure6_publisher_evidence()
    repository = comparison_module.RepositoryFigure6Comparison.build(
        (
            RepositoryFigure6Point.build(Degrees(80.0), Millimetre(10_000.0)),
            RepositoryFigure6Point.build(Degrees(120.0), Millimetre(1_000.0)),
            RepositoryFigure6Point.build(Degrees(160.0), Millimetre(100.0)),
        )
    )
    figure = build_literal_overlay(publisher, repository)
    axis = figure.axes[0]
    line = axis.lines[0]

    assert tuple(line.get_xdata()) == pytest.approx((498.4, 660.6, 822.8))
    assert tuple(line.get_ydata()) == pytest.approx((161.5, 286.0, 410.5))
    assert line.get_label() == "repository requested cap / path length in mm; compliance unaudited"
    assert len(axis.images) == 1
    caption = " ".join(text.get_text() for text in figure.texts).lower()
    assert "literal same-numeric-axis inspection" in caption
    assert "not common-unit or numeric-parity evidence" in caption
