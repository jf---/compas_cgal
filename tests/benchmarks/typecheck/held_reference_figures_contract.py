from pathlib import Path
from typing import Literal
from typing import assert_type

from benchmarks.held_reference_cases import Figure7Observation
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import AxisAlignedDisplayAffine
from benchmarks.held_reference_figures import Figure7PanelEvidence
from benchmarks.held_reference_figures import InwardOffsetEvidence
from benchmarks.held_reference_figures import OverlayMark
from benchmarks.held_reference_figures import SourceCrop
from benchmarks.held_reference_figures import figure7_inward_offset
from benchmarks.held_reference_figures import measure_figure7_observation
from benchmarks.held_reference_figures import measure_figure7_panels
from benchmarks.held_reference_figures import reference_overlay_caption
from benchmarks.held_reference_figures import reference_overlay_caption_lines
from benchmarks.held_reference_figures import reference_overlay_marks
from benchmarks.held_reference_figures import render_reference_overlay
from benchmarks.held_reference_geometry import MillimetresPerPdfPoint
from benchmarks.held_reference_geometry import PdfPoint2
from benchmarks.held_reference_geometry import SourceLine
from benchmarks.held_reference_geometry import SourceToWorld
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY

start = PdfPoint2.build(0.0, 0.0)
end = PdfPoint2.build(1.0, 0.0)
source = assert_type(
    SourceCrop.build(
        name="figure5",
        primitives=(SourceLine.build(start, end),),
        transform=SourceToWorld.build(
            source_origin=start,
            world_origin=Point2[WorldXY].build(0.0, 0.0),
            scale=MillimetresPerPdfPoint(1.0),
        ),
    ),
    SourceCrop,
)
case = load_held_reference_case("figure5")
assert_type(reference_overlay_marks(case, source), tuple[OverlayMark, ...])
assert_type(reference_overlay_caption(case), str)
assert_type(reference_overlay_caption_lines(case), tuple[str, str])
assert_type(render_reference_overlay(case, source, Path("figure5.png")), None)
inward = assert_type(figure7_inward_offset(case), InwardOffsetEvidence)
assert_type(inward.component_provenance, Literal["certified_polygon_projection"])
assert_type(inward.support_provenance, Literal["analytic_boundary_extrema"])
panels = assert_type(measure_figure7_panels(Path("paper.pdf"), case), tuple[Figure7PanelEvidence, ...])
assert_type(measure_figure7_observation(Path("paper.pdf"), case), Figure7Observation)
assert_type(panels[0].registration, AxisAlignedDisplayAffine)
