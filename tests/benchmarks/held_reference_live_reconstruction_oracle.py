from __future__ import annotations

from pathlib import Path

from benchmarks.held_reference_geometry import MillimetresPerPdfPoint
from benchmarks.held_reference_geometry import ReferenceArc
from benchmarks.held_reference_geometry import ReferenceBoundary
from benchmarks.held_reference_geometry import ReferenceLine
from benchmarks.held_reference_geometry import SourceToWorld
from benchmarks.held_reference_geometry import reconstruct_source_path
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY
from tools.held_reference_extractor import extract_reference_sources


def test_live_publisher_cycles_close_the_reconstruction_contract(
    held_publisher_pdf: Path,
) -> None:
    cases = extract_reference_sources(held_publisher_pdf)
    census: list[tuple[str, int, int, int]] = []
    for case in cases:
        scale = 1.0 / float(case.tool_circle.radius)
        transform = SourceToWorld.build(
            source_origin=case.sources[0].start,
            world_origin=Point2[WorldXY].build(0.0, 0.0),
            scale=MillimetresPerPdfPoint(scale),
            reflect_source_y=True,
        )
        limit = Millimetre(float(case.boundary_stroke_width) * scale / 4.0)

        reconstruction = reconstruct_source_path(case.sources, transform, limit)
        transformed_starts = tuple(transform.point(source.start) for source in case.sources)
        twice_signed_area = sum(
            float(start.x) * float(end.y) - float(end.x) * float(start.y) for start, end in zip(transformed_starts, (*transformed_starts[1:], transformed_starts[0]))
        )
        boundary = ReferenceBoundary.build(
            reconstruction.primitives,
            ToolRadius.build(1.0),
            Millimetre(float(case.boundary_stroke_width) * scale),
        )

        assert twice_signed_area > 0.0
        assert boundary.primitives == reconstruction.primitives
        assert float(reconstruction.deviation_upper_bound) <= float(limit)
        assert all(
            current.end == following.start
            for current, following in zip(
                reconstruction.primitives,
                (*reconstruction.primitives[1:], reconstruction.primitives[0]),
            )
        )
        census.append(
            (
                case.name,
                sum(isinstance(primitive, ReferenceLine) for primitive in reconstruction.primitives),
                sum(isinstance(primitive, ReferenceArc) for primitive in reconstruction.primitives),
                len(reconstruction.primitives),
            )
        )

    assert tuple(census) == (
        ("figure-5", 5, 26, 31),
        ("figure-8-upper", 4, 66, 70),
        ("figure-8-skis", 2, 56, 58),
        ("figure-8-monstera", 103, 214, 317),
    )
