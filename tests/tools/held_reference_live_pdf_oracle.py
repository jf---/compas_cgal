from __future__ import annotations

from pathlib import Path

from benchmarks.held_reference_geometry import PdfPointUnit
from benchmarks.held_reference_geometry import SourceCubic
from benchmarks.held_reference_geometry import SourceLine
from tools.held_reference_extractor import extract_reference_sources


def test_live_publisher_pdf_matches_approved_vector_census(
    held_publisher_pdf: Path,
) -> None:
    cases = extract_reference_sources(held_publisher_pdf)

    assert tuple((case.name, case.page) for case in cases) == (
        ("figure-5", 12),
        ("figure-8-upper", 16),
        ("figure-8-skis", 16),
        ("figure-8-monstera", 16),
    )
    expected = (
        (5, 10, 12, 1, 1),
        (4, 22, 16, 1, 1),
        (2, 28, 19, 1, 1),
        (103, 107, 192, 1, 1),
    )
    for case, census in zip(cases, expected):
        assert all(current.end == following.start for current, following in zip(case.sources, (*case.sources[1:], case.sources[0])))
        assert isinstance(case.boundary_stroke_width, float)
        assert isinstance(case.tool_circle.radius, float)
        assert PdfPointUnit(float(case.boundary_stroke_width)) > 0.0
        assert PdfPointUnit(float(case.tool_circle.radius)) > 0.0
        assert (
            sum(isinstance(source, SourceLine) for source in case.sources),
            sum(isinstance(source, SourceCubic) for source in case.sources),
            len(case.boundary_markers),
            1,
            len(case.start_markers),
        ) == census
