"""Ordered physical trial shared by exact swept-prefix segment operations."""

from typing import TypeAlias

from compas_cgal.adaptive.containment import GougeContainment
from compas_cgal.adaptive.coverage import CoverageLedger
from compas_cgal.adaptive.errors import InvalidRetraceSegmentOperationError
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion_certificate import MotionCertifier
from compas_cgal.adaptive.operation import AdvanceSegmentOperation
from compas_cgal.adaptive.operation import RetraceSegmentOperation
from compas_cgal.adaptive.policy import DepletionPolicy
from compas_cgal.adaptive.replay_trace import ReplayLateralWitness
from compas_cgal.adaptive.stock_area import Stock2Area
from compas_cgal.adaptive.units import ToolRadius

SweptPrefixSegmentOperation: TypeAlias = AdvanceSegmentOperation | RetraceSegmentOperation


def _evaluate_swept_prefix_segment_trial(
    *,
    containment_authority: GougeContainment,
    stock: Stock2Area,
    coverage: CoverageLedger,
    operation_index: int,
    operation: SweptPrefixSegmentOperation,
    tool_radius: ToolRadius,
    user_cap: EngagementCap,
    effective_cap: EngagementCap,
    depletion_policy: DepletionPolicy,
) -> ReplayLateralWitness:
    """Evaluate the single ordered exact swept-prefix physical proof path."""
    containment = containment_authority.certify_segment(
        operation.motion,
        tool_radius,
    )
    certifier = MotionCertifier.build(stock=stock, tool_radius=tool_radius)
    motion_witness = certifier.certify_swept_prefix_segment(
        operation_index=operation_index,
        motion=operation.motion,
        user_cap=user_cap,
        effective_cap=effective_cap,
    )
    depletion = stock.deplete(operation.motion, tool_radius, depletion_policy)
    sweep = coverage.add_sweep(operation.motion, tool_radius)
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


def evaluate_retrace_segment_trial(
    *,
    containment_authority: GougeContainment,
    stock: Stock2Area,
    coverage: CoverageLedger,
    operation_index: int,
    operation: RetraceSegmentOperation,
    tool_radius: ToolRadius,
    user_cap: EngagementCap,
    effective_cap: EngagementCap,
    depletion_policy: DepletionPolicy,
) -> ReplayLateralWitness:
    """Evaluate one retrace cut against immutable pre-mutation evidence.

    Exact capsule containment and the two-stratum swept-prefix TEA theorem are
    proved on the pre-state before trial-local stock and coverage owners are
    mutated. The generic event certifier is not part of this deciding path.

    Args:
        containment_authority: Exact reachable-domain containment authority.
        stock: Trial-local exact stock owner.
        coverage: Trial-local exact coverage ledger.
        operation_index: Canonical operation ordinal.
        operation: Source-derived exact retrace segment.
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
            certify this retrace cut.
        InvalidStockAreaError: If exact depletion breaks owned chronology.
        InvalidCoverageSweepError: If exact coverage evidence diverges from
            the trial ledger.
    """
    if type(operation) is not RetraceSegmentOperation:
        raise InvalidRetraceSegmentOperationError(
            "retrace trial requires one exact RetraceSegmentOperation.",
        )
    return _evaluate_swept_prefix_segment_trial(
        containment_authority=containment_authority,
        stock=stock,
        coverage=coverage,
        operation_index=operation_index,
        operation=operation,
        tool_radius=tool_radius,
        user_cap=user_cap,
        effective_cap=effective_cap,
        depletion_policy=depletion_policy,
    )
