from __future__ import annotations

import hashlib

import pytest

from compas_cgal.adaptive.units import Radian
from compas_cgal.engagement_audit.errors import InvalidMeasuredOperationAuditError
from compas_cgal.engagement_audit.errors import InvalidMotionVerdictError
from compas_cgal.engagement_audit.errors import InvalidNonEngagingOperationAuditError
from compas_cgal.engagement_audit.records import MeasuredOperationAudit
from compas_cgal.engagement_audit.records import NonEngagingOperationAudit
from compas_cgal.engagement_audit.records import OperationDigest


def _digest(seed: bytes) -> bytes:
    return hashlib.sha256(seed).digest()


def _measured(**changes: object) -> MeasuredOperationAudit:
    arguments: dict[str, object] = {
        "operation_index": 0,
        "operation_digest": _digest(b"cut"),
        "verdict": "certified",
        "max_tea": Radian(0.5),
        "station_count": 1,
        "pre_motion_stock_lineage": _digest(b"stock"),
        "motion_certificate_digest": _digest(b"certificate"),
    }
    arguments.update(changes)
    return MeasuredOperationAudit.build(**arguments)  # type: ignore[arg-type]


def test_non_engaging_record_has_no_certificate_verdict() -> None:
    record = NonEngagingOperationAudit.build(
        operation_index=0,
        operation_digest=OperationDigest(_digest(b"retract")),
        reason="vertical_retract",
    )

    assert not hasattr(record, "verdict")
    assert not hasattr(record, "max_tea")
    assert not hasattr(record, "station_count")


def test_measured_record_rejects_foreign_verdict() -> None:
    with pytest.raises(InvalidMotionVerdictError, match="foreign"):
        _measured(verdict="foreign")


@pytest.mark.parametrize("station_count", [0, -1, True])
def test_measured_record_requires_positive_exact_station_count(
    station_count: object,
) -> None:
    with pytest.raises(InvalidMeasuredOperationAuditError, match="station count"):
        _measured(station_count=station_count)


def test_unresolved_measurement_remains_a_measured_record() -> None:
    record = _measured(verdict="unresolved", max_tea=Radian(1.25))

    assert record.verdict == "unresolved"
    assert record.max_tea == Radian(1.25)
    assert bytes(record.digest) == hashlib.sha256(record.canonical_bytes).digest()


def test_measurement_digest_binds_pre_motion_stock_and_certificate() -> None:
    original = _measured()
    changed_stock = _measured(pre_motion_stock_lineage=_digest(b"other-stock"))
    changed_certificate = _measured(motion_certificate_digest=_digest(b"other-certificate"))

    assert len({original.digest, changed_stock.digest, changed_certificate.digest}) == 3


def test_non_engaging_record_rejects_foreign_reason() -> None:
    with pytest.raises(InvalidNonEngagingOperationAuditError, match="foreign"):
        NonEngagingOperationAudit.build(
            operation_index=0,
            operation_digest=OperationDigest(_digest(b"transport")),
            reason="foreign",  # type: ignore[arg-type]
        )
