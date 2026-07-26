import math

import numpy as np
import pytest

from compas_cgal import _stock_2
from compas_cgal.adaptive.errors import InvalidMotionCertificateError
from compas_cgal.adaptive.errors import UnresolvedMotionEventError
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.motion_certificate import MotionCertifier
from compas_cgal.adaptive.stock_area import Stock2Area
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import Vector2
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.toolpath import OperationType

SQUARE = np.array(
    [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]],
    dtype=np.float64,
)


def _certifier() -> MotionCertifier:
    return MotionCertifier.build(
        stock=Stock2Area(_stock_2.Stock2(SQUARE, []), ()),
        tool_radius=ToolRadius.build(0.5),
    )


def _circle() -> ExactCircleMotion:
    return ExactCircleMotion.build(
        Point2[WorldXY].build(5.0, 5.0),
        Vector2[WorldXY].build(1.0, 0.0),
        False,
    )


def _segment() -> ExactSegmentMotion:
    return ExactSegmentMotion.build(
        Point2[WorldXY].build(4.0, 5.0),
        Point2[WorldXY].build(6.0, 5.0),
    )


def test_certifier_rejects_invalid_local_contract_before_native_dispatch() -> None:
    certifier = _certifier()

    with pytest.raises(InvalidMotionCertificateError, match="nonnegative"):
        certifier.certify(
            operation_index=-1,
            operation_kind=OperationType.CUT,
            motion=_circle(),
            user_cap=EngagementCap.build(math.pi),
            effective_cap=EngagementCap.build(math.pi),
        )
    with pytest.raises(InvalidMotionCertificateError, match="exceeds"):
        certifier.certify(
            operation_index=0,
            operation_kind=OperationType.CUT,
            motion=_circle(),
            user_cap=EngagementCap.build(math.pi / 2.0),
            effective_cap=EngagementCap.build(math.pi),
        )
    with pytest.raises(InvalidMotionCertificateError, match="incompatible"):
        certifier.certify(
            operation_index=0,
            operation_kind=OperationType.LINK,
            motion=_circle(),
            user_cap=EngagementCap.build(math.pi),
            effective_cap=EngagementCap.build(math.pi),
        )
    with pytest.raises(InvalidMotionCertificateError, match="lateral"):
        certifier.certify(
            operation_index=0,
            operation_kind=OperationType.PLUNGE,
            motion=_segment(),
            user_cap=EngagementCap.build(math.pi),
            effective_cap=EngagementCap.build(math.pi),
        )


def test_real_native_unresolved_is_named_and_const() -> None:
    certifier = _certifier()
    boundary_digest = certifier.canonical_boundary_digest
    lineage_digest = certifier.stock_lineage_digest
    contains = certifier.contains(5.0, 5.0)

    with pytest.raises(UnresolvedMotionEventError, match="unresolved"):
        certifier.certify(
            operation_index=0,
            operation_kind=OperationType.CUT,
            motion=_circle(),
            user_cap=EngagementCap.build(math.pi),
            effective_cap=EngagementCap.build(math.pi),
        )

    assert certifier.canonical_boundary_digest == boundary_digest
    assert certifier.stock_lineage_digest == lineage_digest
    assert certifier.contains(5.0, 5.0) is contains


def test_segment_dependency_gap_fails_loudly() -> None:
    with pytest.raises(UnresolvedMotionEventError, match="segment oracle"):
        _certifier().certify(
            operation_index=0,
            operation_kind=OperationType.LINK,
            motion=_segment(),
            user_cap=EngagementCap.build(math.pi),
            effective_cap=EngagementCap.build(math.pi),
        )


def test_witness_type_is_frozen() -> None:
    from compas_cgal.adaptive.motion_certificate import MotionWitness

    assert MotionWitness.__dataclass_params__.frozen
