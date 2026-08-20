from __future__ import annotations

import json

import pytest

from benchmarks.errors import MalformedRecordError
from benchmarks.measurement import MeasurementRecord


def _record() -> MeasurementRecord:
    """A fully populated record standing in for a real measurement."""
    return MeasurementRecord(
        name="rect_20x10",
        family="analytic",
        params={"width": 20.0, "height": 10.0},
        tool_diameter=1.0,
        tea_cap_deg=120.0,
        generate_seconds=0.01,
        certify_seconds=1.5,
        operations=100,
        cut_operations=90,
        stations=1200,
        max_tea_deg=187.5,
        cap_violations=4,
        unresolved=0,
        arrangement_vertices_final=900,
        max_coordinate_digits=64,
        error=None,
    )


def test_record_round_trips_through_json() -> None:
    record = _record()
    restored = MeasurementRecord.from_dict(json.loads(json.dumps(record.to_dict())))
    assert restored == record


def test_failed_record_carries_the_error_and_zero_timings() -> None:
    record = MeasurementRecord.failed(name="bad", family="necks", params={"pinch": 1.0}, tool_diameter=1.0, tea_cap_deg=120.0, error="DegeneratePocketError: too narrow")
    assert record.error is not None
    assert record.certify_seconds == 0.0
    assert record.cap_violations == 0


def test_failed_record_round_trips_through_json() -> None:
    record = MeasurementRecord.failed(name="bad", family="necks", params={"pinch": 1.0}, tool_diameter=1.0, tea_cap_deg=120.0, error="DegeneratePocketError: too narrow")
    assert MeasurementRecord.from_dict(json.loads(json.dumps(record.to_dict()))) == record


def test_from_dict_rejects_an_unknown_column() -> None:
    payload = _record().to_dict()
    payload["certify_millis"] = 1500.0
    with pytest.raises(MalformedRecordError):
        MeasurementRecord.from_dict(payload)


def test_from_dict_rejects_a_missing_column() -> None:
    payload = _record().to_dict()
    del payload["stations"]
    with pytest.raises(MalformedRecordError):
        MeasurementRecord.from_dict(payload)


def test_from_dict_accepts_an_omitted_error_column() -> None:
    payload = _record().to_dict()
    del payload["error"]
    assert MeasurementRecord.from_dict(payload).error is None
