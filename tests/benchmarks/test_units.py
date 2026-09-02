from __future__ import annotations

import math

import pytest

from benchmarks.errors import InvalidHeldPathEvidenceError
from benchmarks.units import closed_unit_fraction
from benchmarks.units import degrees_value
from benchmarks.units import motion_count
from benchmarks.units import operation_index
from benchmarks.units import seconds_value
from benchmarks.units import tool_radius_multiple


@pytest.mark.parametrize("value", [0, 1.25])
def test_seconds_value_accepts_non_negative_finite_values(value: float) -> None:
    assert seconds_value(value, name="audit") == value


@pytest.mark.parametrize("value", [-20, 0.0, 360])
def test_degrees_value_accepts_finite_values(value: float) -> None:
    assert degrees_value(value, name="reported angle") == value


@pytest.mark.parametrize("value", [0.0, 0.25, 1.0])
def test_closed_unit_fraction_accepts_both_boundaries(value: float) -> None:
    assert closed_unit_fraction(value, name="uncut fraction") == value


@pytest.mark.parametrize("value", [0, 17])
def test_motion_count_accepts_non_negative_integers(value: int) -> None:
    assert motion_count(value, name="gouging motions") == value


@pytest.mark.parametrize("value", [0.0, 2, 3.5])
def test_tool_radius_multiple_accepts_non_negative_finite_values(value: float) -> None:
    assert tool_radius_multiple(value, name="step length") == value


def test_operation_index_accepts_index_in_operation_stream() -> None:
    assert operation_index(3, operation_count=4) == 3


@pytest.mark.parametrize(
    ("validator", "kwargs"),
    [
        (seconds_value, {"name": "audit"}),
        (degrees_value, {"name": "reported angle"}),
        (closed_unit_fraction, {"name": "uncut fraction"}),
        (motion_count, {"name": "gouging motions"}),
        (tool_radius_multiple, {"name": "step length"}),
    ],
)
def test_observation_unit_validators_reject_bool(validator: object, kwargs: dict[str, str]) -> None:
    with pytest.raises(InvalidHeldPathEvidenceError):
        validator(True, **kwargs)  # type: ignore[operator]


@pytest.mark.parametrize("value", [math.nan, math.inf, -math.inf])
@pytest.mark.parametrize(
    ("validator", "name"),
    [
        (seconds_value, "audit"),
        (degrees_value, "reported angle"),
        (closed_unit_fraction, "uncut fraction"),
        (tool_radius_multiple, "step length"),
    ],
)
def test_float_observation_unit_validators_reject_non_finite_values(
    validator: object,
    name: str,
    value: float,
) -> None:
    with pytest.raises(InvalidHeldPathEvidenceError):
        validator(value, name=name)  # type: ignore[operator]


def test_seconds_value_rejects_negative_value() -> None:
    with pytest.raises(InvalidHeldPathEvidenceError):
        seconds_value(-0.01, name="audit")


@pytest.mark.parametrize("value", [-0.01, 1.01])
def test_closed_unit_fraction_rejects_value_outside_closed_interval(value: float) -> None:
    with pytest.raises(InvalidHeldPathEvidenceError):
        closed_unit_fraction(value, name="uncut fraction")


@pytest.mark.parametrize("value", [True, -1, 1.0, math.nan, math.inf, -math.inf])
def test_motion_count_rejects_non_count_value(value: object) -> None:
    with pytest.raises(InvalidHeldPathEvidenceError):
        motion_count(value, name="gouging motions")  # type: ignore[arg-type]


def test_tool_radius_multiple_rejects_negative_value() -> None:
    with pytest.raises(InvalidHeldPathEvidenceError):
        tool_radius_multiple(-0.01, name="step length")


@pytest.mark.parametrize(
    ("value", "operation_count"),
    [
        (True, 4),
        (1.0, 4),
        (math.nan, 4),
        (math.inf, 4),
        (-math.inf, 4),
        (-1, 4),
        (4, 4),
        (0, True),
        (0, 0),
    ],
)
def test_operation_index_rejects_value_outside_operation_stream(value: object, operation_count: object) -> None:
    with pytest.raises(InvalidHeldPathEvidenceError):
        operation_index(value, operation_count=operation_count)  # type: ignore[arg-type]
