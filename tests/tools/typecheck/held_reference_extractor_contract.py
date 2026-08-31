from __future__ import annotations

from typing import assert_type

from benchmarks.held_reference_geometry import PdfPoint2
from benchmarks.held_reference_geometry import PdfPointUnit
from tools.held_reference_extractor import AffineTransform
from tools.held_reference_extractor import BOUNDARY_STROKE_WIDTH
from tools.held_reference_extractor import ExtractedCase
from tools.held_reference_extractor import PdfCrop
from tools.held_reference_extractor import PublishedToolCircle
from tools.held_reference_extractor import SvgPath
from tools.held_reference_extractor import TOOL_CIRCLE_STROKE_WIDTH

transform = AffineTransform.build(
    1.0,
    0.0,
    0.0,
    1.0,
    PdfPointUnit(2.0),
    PdfPointUnit(3.0),
)
assert_type(transform.e, PdfPointUnit)
assert_type(transform.f, PdfPointUnit)
assert_type(transform.point(PdfPointUnit(1.0), PdfPointUnit(2.0)), PdfPoint2)

crop = PdfCrop.build(
    PdfPointUnit(0.0),
    PdfPointUnit(1.0),
    PdfPointUnit(2.0),
    PdfPointUnit(3.0),
)
assert_type(crop.bounds, tuple[PdfPointUnit, PdfPointUnit, PdfPointUnit, PdfPointUnit])
assert_type(BOUNDARY_STROKE_WIDTH, PdfPointUnit)
assert_type(TOOL_CIRCLE_STROKE_WIDTH, PdfPointUnit)

circle = PublishedToolCircle.build(PdfPoint2.build(0.0, 0.0), PdfPointUnit(1.0))
assert_type(circle.radius, PdfPointUnit)


def consume_path(path: SvgPath) -> None:
    assert_type(path.stroke_width, PdfPointUnit)


def consume_case(case: ExtractedCase) -> None:
    assert_type(case.boundary_stroke_width, PdfPointUnit)
