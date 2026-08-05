from dataclasses import replace

import numpy as np
import pytest

from compas_cgal import _containment_2
from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.containment import GougeContainment
from compas_cgal.adaptive.errors import ContainmentConstructionError
from compas_cgal.adaptive.errors import GougeContainmentError
from compas_cgal.adaptive.errors import InvalidContainmentCertificateError
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.reachable_domain import ReachableDomain
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import Vector2
from compas_cgal.adaptive.units import WorldXY

OUTER = (
    (0.0, 0.0),
    (20.0, 0.0),
    (20.0, 20.0),
    (0.0, 20.0),
)
SQUARE_ISLAND = (
    (8.0, 8.0),
    (12.0, 8.0),
    (12.0, 12.0),
    (8.0, 12.0),
)
DIAMOND_ISLAND = (
    (10.0, 8.0),
    (12.0, 10.0),
    (10.0, 12.0),
    (8.0, 10.0),
)
L_SHAPE = (
    (0.0, 0.0),
    (6.0, 0.0),
    (6.0, 2.0),
    (2.0, 2.0),
    (2.0, 6.0),
    (0.0, 6.0),
)


def _ring(
    points: tuple[tuple[float, float], ...],
    *,
    outer: bool,
) -> CanonicalRingV1:
    framed = tuple(Point2[WorldXY].build(x, y) for x, y in points)
    return CanonicalRingV1.build_outer(framed) if outer else CanonicalRingV1.build_hole(framed)


def _domain(
    boundary: tuple[tuple[float, float], ...] = OUTER,
    holes: tuple[tuple[tuple[float, float], ...], ...] = (),
    *,
    tool_radius: float,
) -> ReachableDomain:
    return ReachableDomain.build(
        design_boundary=_ring(boundary, outer=True),
        holes=tuple(_ring(hole, outer=False) for hole in holes),
        tool_radius=ToolRadius.build(tool_radius),
    )


def test_annular_sweep_around_island_does_not_require_filled_outer_disk() -> None:
    domain = _domain(OUTER, (SQUARE_ISLAND,), tool_radius=0.5)
    containment = GougeContainment.build(domain)
    motion = ExactCircleMotion.build(
        Point2[WorldXY].build(10.0, 10.0),
        Vector2[WorldXY].build(4.0, 0.0),
        False,
    )

    certificate = containment.certify_full_circle(
        motion,
        ToolRadius.build(0.5),
    )

    assert certificate.motion is motion
    assert certificate.domain_digest == domain.certificate.digest
    assert certificate.guide_anchor_in_center_domain
    assert not certificate.disk_sweep
    assert not certificate.outer_disk_contained
    assert certificate.native_strategy_version
    assert certificate.native_structural_record
    assert len(certificate.digest) == 32


def test_circle_disk_case_and_exact_outer_tangency_are_accepted() -> None:
    disk_domain = _domain(tool_radius=1.0)
    disk_motion = ExactCircleMotion.build(
        Point2[WorldXY].build(5.0, 5.0),
        Vector2[WorldXY].build(0.5, 0.0),
        True,
    )
    disk_certificate = GougeContainment.build(disk_domain).certify_full_circle(
        disk_motion,
        ToolRadius.build(1.0),
    )

    tangent_domain = _domain(
        (
            (0.0, 0.0),
            (10.0, 0.0),
            (10.0, 10.0),
            (0.0, 10.0),
        ),
        tool_radius=2.0,
    )
    tangent_motion = ExactCircleMotion.build(
        Point2[WorldXY].build(5.0, 5.0),
        Vector2[WorldXY].build(3.0, 0.0),
        False,
    )
    tangent_certificate = GougeContainment.build(tangent_domain).certify_full_circle(
        tangent_motion,
        ToolRadius.build(2.0),
    )

    assert disk_certificate.disk_sweep
    assert disk_certificate.outer_disk_contained
    assert not tangent_certificate.disk_sweep
    assert tangent_certificate.outer_disk_contained

    across = ExactCircleMotion.build(
        tangent_motion.center,
        Vector2[WorldXY].build(np.nextafter(3.0, 4.0), 0.0),
        False,
    )
    with pytest.raises(GougeContainmentError, match="full-circle"):
        GougeContainment.build(tangent_domain).certify_full_circle(
            across,
            ToolRadius.build(2.0),
        )


def test_segment_capsule_accepts_equality_and_rejects_one_binary_quantum_gouge() -> None:
    domain = _domain(
        (
            (0.0, 0.0),
            (10.0, 0.0),
            (10.0, 10.0),
            (0.0, 10.0),
        ),
        tool_radius=1.0,
    )
    containment = GougeContainment.build(domain)
    tangent = ExactSegmentMotion.build(
        Point2[WorldXY].build(1.0, 1.0),
        Point2[WorldXY].build(9.0, 1.0),
    )

    certificate = containment.certify_segment(
        tangent,
        ToolRadius.build(1.0),
    )

    assert certificate.endpoints_in_center_domain
    assert certificate.domain_digest == domain.certificate.digest

    across_y = np.nextafter(1.0, 0.0)
    across = ExactSegmentMotion.build(
        Point2[WorldXY].build(1.0, across_y),
        Point2[WorldXY].build(9.0, across_y),
    )
    with pytest.raises(GougeContainmentError, match="segment"):
        containment.certify_segment(across, ToolRadius.build(1.0))


def test_segment_capsule_accepts_a_lower_dimensional_center_locus() -> None:
    """Certify an exact tool-diameter arm even when 2D erosion drops its axis.

    A radius-1 cutter fits the width-2 horizontal arm exactly. Its center locus
    is one-dimensional and therefore absent from the regularized area-valued
    center domain, but the complete closed capsule remains an exact subset of
    the L-shaped design and must retain independent gouge evidence.
    """
    domain = _domain(L_SHAPE, tool_radius=1.0)
    motion = ExactSegmentMotion.build(
        Point2[WorldXY].build(2.0, 1.0),
        Point2[WorldXY].build(5.0, 1.0),
    )

    certificate = GougeContainment.build(domain).certify_segment(
        motion,
        ToolRadius.build(1.0),
    )

    assert certificate.native_record.contained
    assert not certificate.endpoints_in_center_domain


def test_reflex_crossing_is_rejected_even_when_motion_is_nondegenerate() -> None:
    domain = _domain(L_SHAPE, tool_radius=0.5)
    containment = GougeContainment.build(domain)
    safe = ExactSegmentMotion.build(
        Point2[WorldXY].build(0.5, 0.5),
        Point2[WorldXY].build(5.5, 0.5),
    )
    crossing = ExactSegmentMotion.build(
        Point2[WorldXY].build(1.0, 3.0),
        Point2[WorldXY].build(3.0, 3.0),
    )

    assert containment.certify_segment(
        safe,
        ToolRadius.build(0.5),
    ).endpoints_in_center_domain
    with pytest.raises(GougeContainmentError, match="segment"):
        containment.certify_segment(crossing, ToolRadius.build(0.5))


def test_island_boundary_contact_is_closed_but_one_quantum_incursion_is_not() -> None:
    domain = _domain(OUTER, (DIAMOND_ISLAND,), tool_radius=0.5)
    containment = GougeContainment.build(domain)
    tangent = ExactCircleMotion.build(
        Point2[WorldXY].build(10.0, 10.0),
        Vector2[WorldXY].build(2.5, 0.0),
        False,
    )
    certificate = containment.certify_full_circle(
        tangent,
        ToolRadius.build(0.5),
    )

    assert certificate.guide_anchor_in_center_domain
    assert not certificate.outer_disk_contained

    incursion = ExactCircleMotion.build(
        tangent.center,
        Vector2[WorldXY].build(np.nextafter(2.5, 2.0), 0.0),
        False,
    )
    with pytest.raises(GougeContainmentError, match="full-circle"):
        containment.certify_full_circle(
            incursion,
            ToolRadius.build(0.5),
        )


def test_containment_certificate_cannot_be_relabelled_to_another_motion() -> None:
    domain = _domain(tool_radius=0.5)
    motion = ExactCircleMotion.build(
        Point2[WorldXY].build(10.0, 10.0),
        Vector2[WorldXY].build(2.0, 0.0),
        False,
    )
    certificate = GougeContainment.build(domain).certify_full_circle(
        motion,
        ToolRadius.build(0.5),
    )
    changed_motion = ExactCircleMotion.build(
        motion.center,
        Vector2[WorldXY].build(2.25, 0.0),
        False,
    )

    with pytest.raises(InvalidContainmentCertificateError, match="native record"):
        replace(certificate, motion=changed_motion)


def test_containment_certificate_cannot_be_relabelled_to_another_domain() -> None:
    domain = _domain(tool_radius=0.5)
    other_domain = _domain(
        (
            (0.0, 0.0),
            (24.0, 0.0),
            (24.0, 24.0),
            (0.0, 24.0),
        ),
        tool_radius=0.5,
    )
    containment = GougeContainment.build(domain)
    segment = ExactSegmentMotion.build(
        Point2[WorldXY].build(1.0, 1.0),
        Point2[WorldXY].build(2.0, 1.0),
    )
    circle = ExactCircleMotion.build(
        Point2[WorldXY].build(10.0, 10.0),
        Vector2[WorldXY].build(2.0, 0.0),
        False,
    )

    certificates = (
        containment.certify_segment(segment, ToolRadius.build(0.5)),
        containment.certify_full_circle(circle, ToolRadius.build(0.5)),
    )

    for certificate in certificates:
        with pytest.raises(
            InvalidContainmentCertificateError,
            match="native record",
        ):
            replace(
                certificate,
                domain_digest=other_domain.certificate.digest,
            )


def test_native_containment_requires_one_sha256_authority_digest() -> None:
    domain = _domain(tool_radius=0.5)

    with pytest.raises(ContainmentConstructionError, match="SHA-256"):
        _containment_2.evaluate_exact_segment_containment(
            domain.design_region,
            domain.center_domain,
            b"not-a-digest",
            1.0,
            1.0,
            2.0,
            1.0,
            0.5,
        )
