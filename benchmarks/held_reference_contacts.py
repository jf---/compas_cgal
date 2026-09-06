"""Qualify draft offset contacts without erasing rejected input evidence."""

from __future__ import annotations

from fractions import Fraction

from benchmarks.errors import BenchmarkError
from benchmarks.held_figure5_raw_guide import ProjectionAdmissibleBoundaryHypothesis
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from compas_cgal import _circle_geometry_2
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY


class MissingReferenceContactRunError(BenchmarkError):
    """The qualified offset contacts cannot represent every input guide run."""


def qualify_reference_contacts(
    hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...],
    boundary: tuple[Point2[WorldXY], ...],
    evidence_bound: Millimetre,
) -> tuple[dict[ProjectionAdmissibleBoundaryHypothesis, ProjectionBoundarySite], tuple[ProjectionAdmissibleBoundaryHypothesis, ...]]:
    """Separate bounded contacts from native distance rejections.

    CGAL owns projection, tie handling and distance decisions. Returned site
    parameters are approximate construction views. A rejected hypothesis is
    preserved verbatim; ambiguous projections and malformed geometry raise.
    Empty input or loss of an entire source run also raises.
    """
    if not hypotheses:
        raise MissingReferenceContactRunError("Contact qualification requires input guide hypotheses.")
    polygon = [(float(point.x), float(point.y)) for point in boundary]
    sites = {}
    rejected = []
    for hypothesis in hypotheses:
        try:
            side, parameter = _circle_geometry_2.project_boundary_contact(polygon, (float(hypothesis.contact_point.x), float(hypothesis.contact_point.y)), float(evidence_bound))
        except _circle_geometry_2.BoundaryProjectionDistanceError:
            rejected.append(hypothesis)
        else:
            # Representation bridge only; projection arithmetic stays native.
            sites[hypothesis] = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(side), Fraction(parameter), len(boundary))
    missing = {hypothesis.run_id for hypothesis in hypotheses} - {hypothesis.run_id for hypothesis in sites}
    if missing:
        raise MissingReferenceContactRunError(f"No bounded offset contact for guide runs {sorted(map(int, missing))}.")
    return sites, tuple(rejected)
