"""Ordered paper-derived machining-circle proposals along one MAT edge."""

from compas_cgal.adaptive.candidates import MathsmCircleProposal
from compas_cgal.adaptive.errors import InvalidMathsmProposalError
from compas_cgal.adaptive.medial_axis import MatEdge
from compas_cgal.adaptive.medial_axis import MatSite
from compas_cgal.adaptive.medial_axis import MedialAxis
from compas_cgal.adaptive.policy import CircleOrientation


def maximum_radius_edge_proposals(
    *,
    axis: MedialAxis,
    edge: MatEdge,
    generator_site: MatSite,
    circle_orientation: CircleOrientation,
) -> tuple[MathsmCircleProposal, ...]:
    """Emit maximum-radius, boundary-phase circles in MAT sample order.

    Args:
        axis: Exact typed medial-axis owner.
        edge: One edge owned by `axis`.
        generator_site: The boundary side used for the one-sided middle curve.
        circle_orientation: Travel direction for every full circle.

    Returns:
        One proposal per native sample, ordered by its edge ordinal.

    Raises:
        InvalidMathsmProposalError: Inputs are foreign, the edge has no samples,
            or a sample has no positive paper-derived machining-circle radius.
    """
    if type(axis) is not MedialAxis:
        raise InvalidMathsmProposalError(
            "MATHSM edge path requires one exact typed MAT owner.",
        )
    if type(edge) is not MatEdge or axis.edge_by_id.get(edge.identity) != edge:
        raise InvalidMathsmProposalError(
            "MATHSM edge path requires an edge owned by its MAT.",
        )
    if type(generator_site) is not MatSite or axis.site_by_id.get(generator_site.identity) != generator_site or generator_site.identity not in edge.generator_site_ids:
        raise InvalidMathsmProposalError(
            "MATHSM edge path generator site must be owned by and bound to its edge.",
        )
    samples = tuple(
        sorted(
            (sample for sample in axis.samples if sample.edge_id == edge.identity),
            key=lambda sample: sample.ordinal_on_edge,
        )
    )
    if not samples:
        raise InvalidMathsmProposalError(
            "MATHSM edge path requires at least one owned native sample.",
        )
    return tuple(
        MathsmCircleProposal.build_maximum(
            generator_site=generator_site,
            middle_point=sample.point,
            tool_radius=axis.tool_radius,
            circle_orientation=circle_orientation,
        )
        for sample in samples
    )
