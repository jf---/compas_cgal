"""Ordered physical trial for one exact zero-guide segment advance."""

from compas_cgal.adaptive.containment import GougeContainment
from compas_cgal.adaptive.containment import SegmentContainmentCertificate
from compas_cgal.adaptive.coverage import CoverageLedger
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion_certificate import MotionCertifier
from compas_cgal.adaptive.operation import AdvanceSegmentOperation
from compas_cgal.adaptive.policy import DepletionPolicy
from compas_cgal.adaptive.replay_trace import ReplayLateralWitness
from compas_cgal.adaptive.stock_area import Stock2Area
from compas_cgal.adaptive.units import ToolRadius


def evaluate_advancing_segment_trial(
    *,
    containment_authority: GougeContainment,
    stock: Stock2Area,
    coverage: CoverageLedger,
    operation_index: int,
    operation: AdvanceSegmentOperation,
    tool_radius: ToolRadius,
    user_cap: EngagementCap,
    effective_cap: EngagementCap,
    depletion_policy: DepletionPolicy,
) -> ReplayLateralWitness:
    """Evaluate one advancing cut against immutable pre-mutation evidence.

    The order is part of the safety contract: exact capsule containment and
    swept-prefix TEA are proved on the pre-state before trial-local stock and
    coverage owners are mutated.

    Args:
        containment_authority: Exact reachable-domain containment authority.
        stock: Trial-local exact stock owner.
        coverage: Trial-local exact coverage ledger.
        operation_index: Canonical operation ordinal.
        operation: Typed advancing segment and traversal decision.
        tool_radius: Typed cutter radius.
        user_cap: User engagement limit.
        effective_cap: Policy-derived engagement limit.
        depletion_policy: Exact under-covering depletion policy.

    Returns:
        Cross-validated replay witness after both trial-local mutations.

    Raises:
        GougeContainmentError: If the exact cutter capsule leaves the design.
        InvalidMotionCertificateError: If native proof structure contradicts
            the submitted motion or stock identity.
        UnresolvedMotionEventError: If the clear-start exact-pi theorem cannot
            certify this advancing cut.
        InvalidStockAreaError: If exact depletion breaks owned chronology.
        InvalidCoverageSweepError: If exact coverage evidence diverges from
            the trial ledger.
    """
    containment: SegmentContainmentCertificate = containment_authority.certify_segment(
        operation.motion,
        tool_radius,
    )
    certifier = MotionCertifier.build(
        stock=stock,
        tool_radius=tool_radius,
    )
    motion_witness = certifier.certify_swept_prefix_segment(
        operation_index=operation_index,
        motion=operation.motion,
        user_cap=user_cap,
        effective_cap=effective_cap,
    )
    depletion = stock.deplete(
        operation.motion,
        tool_radius,
        depletion_policy,
    )
    sweep = coverage.add_sweep(
        operation.motion,
        tool_radius,
    )
    return ReplayLateralWitness(
        operation_index=operation_index,
        operation=operation,
        effective_cap_decision=operation.effective_cap_decision,
        stock_boundary_digest=certifier.canonical_boundary_digest,
        containment_certificate=containment,
        motion_witness=motion_witness,
        depletion_witness=depletion,
        sweep_witness=sweep,
    )
