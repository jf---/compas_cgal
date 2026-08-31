import math
from typing import assert_type

from benchmarks.held_reference_geometry import MillimetresPerPdfPoint
from benchmarks.held_reference_geometry import PdfPoint2
from benchmarks.held_reference_geometry import PolygonProjection
from benchmarks.held_reference_geometry import ReferenceArc
from benchmarks.held_reference_geometry import ReferenceBoundary
from benchmarks.held_reference_geometry import ReferenceLine
from benchmarks.held_reference_geometry import SourceCubic
from benchmarks.held_reference_geometry import SourceLine
from benchmarks.held_reference_geometry import SourceToWorld
from benchmarks.held_reference_geometry import project_boundary
from benchmarks.held_reference_geometry import reconstruct_cubic
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Radian
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY

pdf_start = assert_type(PdfPoint2.build(1.0, 0.0), PdfPoint2)
pdf_end = assert_type(PdfPoint2.build([0.0, 1.0]), PdfPoint2)
pdf_control1 = PdfPoint2.build(1.0, 0.5)
pdf_control2 = PdfPoint2.build(0.5, 1.0)
source_line = assert_type(SourceLine.build(pdf_start, pdf_end), SourceLine)
source_cubic = assert_type(
    SourceCubic.build(pdf_start, pdf_control1, pdf_control2, pdf_end),
    SourceCubic,
)
transform = assert_type(
    SourceToWorld.build(
        source_origin=PdfPoint2.build(0.0, 0.0),
        world_origin=Point2[WorldXY].build(0.0, 0.0),
        scale=MillimetresPerPdfPoint(1.0),
    ),
    SourceToWorld,
)
world_start = assert_type(transform.point(pdf_start), Point2[WorldXY])
world_end = transform.point(pdf_end)
reference_line = assert_type(ReferenceLine.build(world_start, world_end), ReferenceLine)
reference_arc = assert_type(
    ReferenceArc.build(
        Point2[WorldXY].build(1.0, 0.0),
        Point2[WorldXY].build(0.0, 1.0),
        Point2[WorldXY].build(0.0, 0.0),
        Radian(math.pi / 2.0),
    ),
    ReferenceArc,
)
corners = (
    Point2[WorldXY].build(0.0, 0.0),
    Point2[WorldXY].build(1.0, 0.0),
    Point2[WorldXY].build(1.0, 1.0),
    Point2[WorldXY].build(0.0, 1.0),
)
closed_lines = tuple(ReferenceLine.build(start, end) for start, end in zip(corners, (*corners[1:], corners[0])))
boundary = assert_type(
    ReferenceBoundary.build(
        closed_lines,
        ToolRadius.build(1.0),
        Millimetre(0.1),
    ),
    ReferenceBoundary,
)
polygon_projection = assert_type(
    PolygonProjection.build(
        corners,
        Millimetre(0.01),
        Millimetre(0.0),
    ),
    PolygonProjection,
)
projection = assert_type(project_boundary(boundary, Millimetre(0.01)), PolygonProjection)
arcs = assert_type(
    reconstruct_cubic(source_cubic, transform, Millimetre(0.01)),
    tuple[ReferenceArc, ...],
)
