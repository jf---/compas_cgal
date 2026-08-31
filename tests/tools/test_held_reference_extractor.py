from __future__ import annotations

from pathlib import Path

import pytest

from benchmarks.errors import AmbiguousPublishedBoundaryError
from benchmarks.errors import MissingPublishedBoundaryMarkerError
from benchmarks.errors import MissingPublishedToolCircleError
from benchmarks.errors import UnsupportedPdfBoundaryOperatorError
from benchmarks.held_reference_geometry import PdfPoint2
from benchmarks.held_reference_geometry import SourceCubic
from benchmarks.held_reference_geometry import SourceLine
from tools.held_reference_extractor import AffineTransform
from tools.held_reference_extractor import FIGURE_CROPS
from tools.held_reference_extractor import FigureCrop
from tools.held_reference_extractor import PdfCrop
from tools.held_reference_extractor import PDF_COORDINATE_QUANTUM_PT
from tools.held_reference_extractor import canonicalize_degree_one_endpoints
from tools.held_reference_extractor import extract_case_from_svg
from tools.held_reference_extractor import parse_pdf_svg_path
from tools.held_reference_extractor import select_boundary_paths
from tools.held_reference_extractor import select_tool_circle

GREEN = "rgb(17.999268%, 54.499817%, 34.098816%)"
RED = "rgb(100%, 0%, 0%)"
BLACK = "rgb(0%, 0%, 0%)"


def _path(
    path_data: str,
    *,
    stroke: str = GREEN,
    width: float = 2.0,
    transform: str = "matrix(1, 0, 0, 1, 0, 0)",
) -> str:
    return f'<path fill="none" stroke-width="{width}" stroke="{stroke}" stroke-linejoin="round" d="{path_data}" transform="{transform}"/>'


def _write_svg(path: Path, *paths: str) -> None:
    path.write_text('<svg xmlns="http://www.w3.org/2000/svg">' + "".join(paths) + "</svg>")


def _rectangle(x0: float, y0: float, x1: float, y1: float) -> tuple[str, ...]:
    return (
        _path(f"M {x0} {y0} L {x1} {y0}"),
        _path(f"M {x1} {y0} L {x1} {y1}"),
        _path(f"M {x1} {y1} L {x0} {y1}"),
        _path(f"M {x0} {y1} L {x0} {y0}"),
    )


def _black_marker(x: float, y: float) -> str:
    radius = 0.125
    handle = 0.069036
    return (
        f'<path fill-rule="evenodd" fill="{BLACK}" fill-opacity="1" '
        f'd="M {x + radius} {y} '
        f"C {x + radius} {y + handle} {x + handle} {y + radius} {x} {y + radius} "
        f"C {x - handle} {y + radius} {x - radius} {y + handle} {x - radius} {y} "
        f"C {x - radius} {y - handle} {x - handle} {y - radius} {x} {y - radius} "
        f'C {x + handle} {y - radius} {x + radius} {y - handle} {x + radius} {y}"/>'
    )


def _rectangle_markers(x0: float, y0: float, x1: float, y1: float) -> tuple[str, ...]:
    return (
        _black_marker(x0, y0),
        _black_marker(x1, y0),
        _black_marker(x1, y1),
        _black_marker(x0, y1),
    )


def _crop(*, expected_tool_circles: int = 1) -> FigureCrop:
    return FigureCrop.build(
        name="fixture",
        page=12,
        crop=PdfCrop.build(0.0, 0.0, 10.0, 10.0),
        transform=AffineTransform.identity(),
        expected_line_count=4,
        expected_cubic_count=0,
        expected_tool_circle_count=expected_tool_circles,
        normalization_radius_local=1.0,
    )


def test_parser_applies_affine_transform_to_line() -> None:
    primitives = parse_pdf_svg_path(
        "M 1 2 L 3 4",
        AffineTransform.build(2.0, 0.0, 0.0, -2.0, 10.0, 20.0),
    )

    assert primitives == (
        SourceLine.build(
            PdfPoint2.build(12.0, 16.0),
            PdfPoint2.build(16.0, 12.0),
        ),
    )


def test_parser_splits_compound_cubics_at_authored_junctions() -> None:
    primitives = parse_pdf_svg_path(
        "M 0 0 C 1 0 2 1 3 1 C 4 1 5 0 6 0",
        AffineTransform.identity(),
    )

    assert primitives == (
        SourceCubic.build(
            PdfPoint2.build(0.0, 0.0),
            PdfPoint2.build(1.0, 0.0),
            PdfPoint2.build(2.0, 1.0),
            PdfPoint2.build(3.0, 1.0),
        ),
        SourceCubic.build(
            PdfPoint2.build(3.0, 1.0),
            PdfPoint2.build(4.0, 1.0),
            PdfPoint2.build(5.0, 0.0),
            PdfPoint2.build(6.0, 0.0),
        ),
    )


@pytest.mark.parametrize(
    "path_data",
    [
        "M 0 0 L 1 0 Z",
        "m 0 0 l 1 0",
        "M 0 0 L 1 0 C 2 0 3 1 4 1",
        "M 0 0 C 1 0 2 1",
        "M 0 0 L nan 1",
        "M 0 0 C 1 0 2 1 3 1 C 4 1 5 0 6 0 C 7 0 8 1 9 1 C 10 1 11 0 12 0 C 13 0 14 1 15 1",
    ],
)
def test_parser_rejects_unsupported_or_malformed_path(path_data: str) -> None:
    with pytest.raises(UnsupportedPdfBoundaryOperatorError):
        parse_pdf_svg_path(path_data, AffineTransform.identity())


def test_crop_rejects_colour_only_distractors(tmp_path: Path) -> None:
    svg = tmp_path / "fixture.svg"
    _write_svg(
        svg,
        *_rectangle(1.0, 1.0, 4.0, 4.0),
        *_rectangle_markers(1.0, 1.0, 4.0, 4.0),
        _path("M 20 20 L 21 20"),
        _path("M 6 6 L 7 6", width=1.0),
        _path("M 7 7 L 8 7"),
        _path("M 1 1 L 4 1", transform="matrix(1, 0, 0, 1, 0.5, 0)"),
    )

    selected = select_boundary_paths(svg, _crop())

    assert len(selected) == 4
    assert all(path.stroke_width == pytest.approx(2.0) for path in selected)
    assert sum(len(path.primitives) for path in selected) == 4


def test_crop_rejects_two_valid_boundary_cycles(tmp_path: Path) -> None:
    svg = tmp_path / "ambiguous.svg"
    _write_svg(
        svg,
        *_rectangle(1.0, 1.0, 3.0, 3.0),
        *_rectangle(6.0, 6.0, 8.0, 8.0),
        *_rectangle_markers(1.0, 1.0, 3.0, 3.0),
        *_rectangle_markers(6.0, 6.0, 8.0, 8.0),
    )

    with pytest.raises(AmbiguousPublishedBoundaryError):
        select_boundary_paths(svg, _crop())


def test_tool_circle_uses_only_red_width_point_eight_closed_four_cubic_path(
    tmp_path: Path,
) -> None:
    svg = tmp_path / "circle.svg"
    circle = "M 5 4 C 5.552285 4 6 4.447715 6 5 C 6 5.552285 5.552285 6 5 6 C 4.447715 6 4 5.552285 4 5 C 4 4.447715 4.447715 4 5 4"
    _write_svg(
        svg,
        _path(circle, stroke=RED, width=0.8),
        _path(circle, stroke=RED, width=2.0),
    )

    selected = select_tool_circle(svg, _crop())

    assert selected.centre == PdfPoint2.build(5.0, 5.0)
    assert float(selected.radius) == pytest.approx(1.0)


def test_tool_circle_fails_when_scale_circle_is_ambiguous(tmp_path: Path) -> None:
    svg = tmp_path / "circles.svg"
    first = "M 2 1 C 2.5 1 3 1.5 3 2 C 3 2.5 2.5 3 2 3 C 1.5 3 1 2.5 1 2 C 1 1.5 1.5 1 2 1"
    second = "M 7 6 C 7.5 6 8 6.5 8 7 C 8 7.5 7.5 8 7 8 C 6.5 8 6 7.5 6 7 C 6 6.5 6.5 6 7 6"
    _write_svg(svg, _path(first, stroke=RED, width=0.8), _path(second, stroke=RED, width=0.8))

    with pytest.raises(MissingPublishedToolCircleError):
        select_tool_circle(svg, _crop(expected_tool_circles=1))


def test_tool_circle_accepts_one_quantized_extremum_discrepancy(tmp_path: Path) -> None:
    svg = tmp_path / "quantized-circle.svg"
    circle = "M 5 4 C 5.552285 4 6 4.447715 6 5.001953 C 6 5.554238 5.552285 6.003906 5 6.003906 C 4.447715 6.003906 4 5.554238 4 5.001953 C 4 4.447715 4.447715 4 5 4"
    _write_svg(svg, _path(circle, stroke=RED, width=0.8))

    selected = select_tool_circle(svg, _crop())

    assert float(selected.radius) == pytest.approx(1.0)


def test_tool_circle_radius_is_rotation_invariant(tmp_path: Path) -> None:
    svg = tmp_path / "rotated-circle.svg"
    circle = (
        "M 5.866025404 5.5 C 5.589883029 5.978292623 4.978292623 6.142167779 4.5 5.866025404 "
        "C 4.021707377 5.589883029 3.857832221 4.978292623 4.133974596 4.5 "
        "C 4.410116971 4.021707377 5.021707377 3.857832221 5.5 4.133974596 "
        "C 5.978292623 4.410116971 6.142167779 5.021707377 5.866025404 5.5"
    )
    _write_svg(svg, _path(circle, stroke=RED, width=0.8))

    selected = select_tool_circle(svg, _crop())

    assert float(selected.radius) == pytest.approx(1.0, abs=0.000001)


def test_tool_circle_rejects_radially_displaced_handle(tmp_path: Path) -> None:
    svg = tmp_path / "malformed-circle.svg"
    malformed = "M 5 4 C 5.2 4.2 6 4.447715 6 5 C 6 5.552285 5.552285 6 5 6 C 4.447715 6 4 5.552285 4 5 C 4 4.447715 4.447715 4 5 4"
    _write_svg(svg, _path(malformed, stroke=RED, width=0.8))

    with pytest.raises(MissingPublishedToolCircleError):
        select_tool_circle(svg, _crop())


def test_unique_coordinate_quantum_seams_are_canonicalized() -> None:
    primitives = (
        SourceLine.build(PdfPoint2.build(0.0, 0.0), PdfPoint2.build(1.0, 0.0)),
        SourceLine.build(PdfPoint2.build(1.0, 0.003906466), PdfPoint2.build(5.0, 0.0)),
        SourceLine.build(PdfPoint2.build(10.0, 0.0), PdfPoint2.build(11.0, 0.0)),
        SourceLine.build(PdfPoint2.build(11.003906, 0.00006), PdfPoint2.build(15.0, 0.0)),
        SourceLine.build(PdfPoint2.build(20.0, 0.0), PdfPoint2.build(21.0, 0.0)),
        SourceLine.build(PdfPoint2.build(21.003906, 0.0), PdfPoint2.build(25.0, 0.0)),
    )

    canonical = canonicalize_degree_one_endpoints(primitives)

    assert canonical[0].end == canonical[1].start
    assert canonical[2].end == canonical[3].start
    assert canonical[4].end == canonical[5].start
    assert canonical[1].end != canonical[2].start


def test_diagonal_gap_outside_euclidean_quantum_is_not_canonicalized() -> None:
    component = 0.75 * float(PDF_COORDINATE_QUANTUM_PT)
    primitives = (
        SourceLine.build(PdfPoint2.build(0.0, 0.0), PdfPoint2.build(1.0, 1.0)),
        SourceLine.build(PdfPoint2.build(1.0 + component, 1.0 + component), PdfPoint2.build(2.0, 2.0)),
    )

    canonical = canonicalize_degree_one_endpoints(primitives)

    assert canonical[0].end != canonical[1].start


def test_boundary_rejects_unrelated_black_markers(tmp_path: Path) -> None:
    svg = tmp_path / "unrelated-markers.svg"
    _write_svg(
        svg,
        *_rectangle(1.0, 1.0, 4.0, 4.0),
        *_rectangle_markers(5.0, 5.0, 8.0, 8.0),
    )

    with pytest.raises(MissingPublishedBoundaryMarkerError):
        select_boundary_paths(svg, _crop())


def test_boundary_marker_crop_membership_uses_validated_centre(tmp_path: Path) -> None:
    svg = tmp_path / "marker-centres.svg"
    _write_svg(
        svg,
        *_rectangle(1.0, 0.0, 4.0, 4.0),
        *_rectangle_markers(1.0, 0.0, 4.0, 4.0),
    )

    selected = select_boundary_paths(svg, _crop())

    assert len(selected) == 4


def test_approved_crops_hold_measured_publisher_families() -> None:
    observed = {
        crop.name: (
            crop.page,
            tuple(float(value) for value in crop.crop.bounds),
            crop.expected_line_count,
            crop.expected_cubic_count,
            float(crop.normalization_radius_local),
        )
        for crop in FIGURE_CROPS
    }

    assert observed == {
        "figure-5": (12, (40.0, 70.0, 265.0, 223.0), 5, 10, 6.769349),
        "figure-8-upper": (16, (132.0, 70.0, 410.0, 259.0), 4, 22, 5.977675),
        "figure-8-skis": (16, (132.0, 253.0, 407.0, 322.0), 2, 28, 5.118423),
        "figure-8-monstera": (16, (133.0, 329.0, 407.0, 617.0), 103, 107, 4.863683),
    }


def test_case_extraction_orders_the_selected_cycle(tmp_path: Path) -> None:
    svg = tmp_path / "case.svg"
    circle = "M 5 4 C 5.552285 4 6 4.447715 6 5 C 6 5.552285 5.552285 6 5 6 C 4.447715 6 4 5.552285 4 5 C 4 4.447715 4.447715 4 5 4"
    boundary = _rectangle(1.0, 1.0, 4.0, 4.0)
    _write_svg(
        svg,
        boundary[2],
        boundary[0],
        boundary[3],
        boundary[1],
        *_rectangle_markers(1.0, 1.0, 4.0, 4.0),
        _path(circle, stroke=RED, width=0.8),
    )

    extracted = extract_case_from_svg(svg, _crop())

    assert len(extracted.sources) == 4
    assert len(extracted.boundary_markers) == 4
    assert all(current.end == following.start for current, following in zip(extracted.sources, (*extracted.sources[1:], extracted.sources[0])))
    assert float(extracted.tool_circle.radius) == pytest.approx(1.0)
