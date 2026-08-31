from typing import assert_type
from typing import cast

from benchmarks.held_reference_geometry import MillimetresPerPdfPoint
from benchmarks.held_reference_geometry import PdfPoint2
from benchmarks.held_reference_geometry import PolygonProjection
from benchmarks.held_reference_geometry import ReferenceArc
from benchmarks.held_reference_geometry import ReferenceBoundary
from benchmarks.held_reference_geometry import SourceCubic
from benchmarks.held_reference_geometry import SourceToWorld
from benchmarks.held_reference_geometry import project_boundary
from benchmarks.held_reference_geometry import reconstruct_cubic
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY

pdf_start = assert_type(PdfPoint2.build(1.0, 0.0), PdfPoint2)
pdf_control1 = PdfPoint2.build(1.0, 0.5)
pdf_control2 = PdfPoint2.build(0.5, 1.0)
pdf_end = PdfPoint2.build(0.0, 1.0)
source = assert_type(
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
arcs = assert_type(
    reconstruct_cubic(source, transform, Millimetre(0.01)),
    tuple[ReferenceArc, ...],
)

# The projection's consumer type remains explicit even though full boundary
# fixtures belong in runtime tests rather than this static contract.
boundary = cast(ReferenceBoundary, object())
projection = assert_type(project_boundary(boundary, Millimetre(0.01)), PolygonProjection)
tool_radius = assert_type(ToolRadius.build(1.0), ToolRadius)
