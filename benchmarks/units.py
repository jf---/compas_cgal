from __future__ import annotations

import math
from typing import NewType

from benchmarks.errors import InvalidHeldPathEvidenceError

Degrees = NewType("Degrees", float)
Seconds = NewType("Seconds", float)
UnitFraction = NewType("UnitFraction", float)
MotionCount = NewType("MotionCount", int)
ToolRadiusMultiple = NewType("ToolRadiusMultiple", float)
OperationIndex = NewType("OperationIndex", int)


def _finite(value: float, *, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise InvalidHeldPathEvidenceError(f"{name} must be a finite real number.")
    try:
        numeric = float(value)
    except OverflowError:
        raise InvalidHeldPathEvidenceError(f"{name} exceeds the binary64 range.") from None
    if not math.isfinite(numeric):
        raise InvalidHeldPathEvidenceError(f"{name} must be finite.")
    return numeric


def _non_negative(value: float, *, name: str) -> float:
    numeric = _finite(value, name=name)
    if numeric < 0.0:
        raise InvalidHeldPathEvidenceError(f"{name} must be non-negative.")
    return numeric


def seconds_value(value: float, *, name: str) -> Seconds:
    return Seconds(_non_negative(value, name=name))


def degrees_value(value: float, *, name: str) -> Degrees:
    return Degrees(_finite(value, name=name))


def closed_unit_fraction(value: float, *, name: str) -> UnitFraction:
    numeric = _finite(value, name=name)
    if not 0.0 <= numeric <= 1.0:
        raise InvalidHeldPathEvidenceError(f"{name} must lie in [0, 1].")
    return UnitFraction(numeric)


def motion_count(value: int, *, name: str) -> MotionCount:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise InvalidHeldPathEvidenceError(f"{name} must be a non-negative integer.")
    return MotionCount(value)


def tool_radius_multiple(value: float, *, name: str) -> ToolRadiusMultiple:
    return ToolRadiusMultiple(_non_negative(value, name=name))


def operation_index(value: int, *, operation_count: int) -> OperationIndex:
    if isinstance(operation_count, bool) or not isinstance(operation_count, int) or operation_count <= 0:
        raise InvalidHeldPathEvidenceError("operation count must be a positive integer.")
    if isinstance(value, bool) or not isinstance(value, int) or not 0 <= value < operation_count:
        raise InvalidHeldPathEvidenceError(f"operation index must lie in [0, {operation_count}), got {value!r}.")
    return OperationIndex(value)
