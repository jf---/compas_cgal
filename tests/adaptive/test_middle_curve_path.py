from fractions import Fraction

import numpy as np
import pytest

from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.errors import InvalidMathsmProposalError
from compas_cgal.adaptive.medial_axis import MatEdge
from compas_cgal.adaptive.medial_axis import MatSite
from compas_cgal.adaptive.medial_axis import MedialAxis
from compas_cgal.adaptive.middle_curve_path import maximum_radius_edge_proposals
from compas_cgal.adaptive.policy import CircleOrientation
from compas_cgal.adaptive.units import ChordBound
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Spacing
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


L_SHAPE = np.asarray(
    (
        (0.0, 0.0),
        (6.0, 0.0),
        (6.0, 2.0),
        (2.0, 2.0),
        (2.0, 6.0),
        (0.0, 6.0),
    ),
    dtype=np.float64,
)


def _axis() -> MedialAxis:
    boundary = CanonicalRingV1.build_outer(
        tuple(Point2[WorldXY].build(x, y) for x, y in L_SHAPE),
    )
    return MedialAxis.build(
        design_boundary=boundary,
        holes=(),
        tool_radius=ToolRadius.build(0.5),
        station_spacing=Spacing.build(0.75),
        max_sagitta=ChordBound.build(0.02),
        max_refinement_depth=32,
    )


def _constant_clearance_edge(axis: MedialAxis) -> tuple[MatEdge, MatSite]:
    samples_by_edge = {edge.identity: tuple(sample for sample in axis.samples if sample.edge_id == edge.identity) for edge in axis.edges}
    edge = next(
        edge
        for edge in axis.edges
        if edge.curve_kind == "line" and len(samples_by_edge[edge.identity]) == 5 and samples_by_edge[edge.identity][0].point == Point2[WorldXY].build(2.0, 1.0)
    )
    site = next(axis.site_by_id[site_id] for site_id in edge.generator_site_ids if axis.site_by_id[site_id].kind == "open-segment" and axis.site_by_id[site_id].source.y == 0.0)
    return edge, site


def test_maximum_radius_edge_proposals_follow_mat_progress_not_policy_order() -> None:
    axis = _axis()
    edge, site = _constant_clearance_edge(axis)
    edge_samples = tuple(
        sorted(
            (sample for sample in axis.samples if sample.edge_id == edge.identity),
            key=lambda sample: sample.ordinal_on_edge,
        )
    )

    proposals = maximum_radius_edge_proposals(
        axis=axis,
        edge=edge,
        generator_site=site,
        circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
    )

    assert tuple(proposal.middle_point for proposal in proposals) == tuple(sample.point for sample in edge_samples)
    assert all(proposal.middle_point.y == 1.0 for proposal in proposals)
    assert all(proposal.motion.center.y == 0.75 for proposal in proposals)
    assert all(proposal.guide_radius == Fraction(1, 4) for proposal in proposals)
    assert all(
        Point2[WorldXY].build(
            proposal.motion.center.x + proposal.motion.phase_vector.x,
            proposal.motion.center.y + proposal.motion.phase_vector.y,
        )
        == proposal.mathsm_contact_point
        for proposal in proposals
    )
    assert proposals == maximum_radius_edge_proposals(
        axis=axis,
        edge=edge,
        generator_site=site,
        circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
    )


def test_maximum_radius_edge_proposals_reject_foreign_generator_site() -> None:
    axis = _axis()
    edge, _ = _constant_clearance_edge(axis)
    foreign_site = next(site for site in axis.sites if site.identity not in edge.generator_site_ids)

    with pytest.raises(InvalidMathsmProposalError, match="generator site"):
        maximum_radius_edge_proposals(
            axis=axis,
            edge=edge,
            generator_site=foreign_site,
            circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
        )
