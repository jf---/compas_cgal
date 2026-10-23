"""Typed reads of a parsed timing record, each failing with ``MalformedTimingRecordError``.

The adapters in ``tools.perf_history_adapters`` read many hand-shaped JSON
records; these accessors make every read checked, so a record whose shape
drifted fails loudly naming the field instead of yielding a wrong entry.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from collections.abc import Sequence
from typing import Any
from typing import Optional


class MalformedTimingRecordError(Exception):
    """A timing record of a known family does not have the shape its adapter reads."""


def obj(value: object, where: str) -> Mapping[str, Any]:
    """Raises:
    MalformedTimingRecordError: ``value`` is not a JSON object.
    """
    if not isinstance(value, dict):
        raise MalformedTimingRecordError(f"{where} is an object, not {type(value).__name__}.")
    return value


def items(value: object, where: str) -> Sequence[Any]:
    """Raises:
    MalformedTimingRecordError: ``value`` is not a JSON list.
    """
    if not isinstance(value, list):
        raise MalformedTimingRecordError(f"{where} is a list, not {type(value).__name__}.")
    return value


def field(record: Mapping[str, Any], key: str, where: str) -> Any:
    """Raises:
    MalformedTimingRecordError: ``key`` is missing.
    """
    if key not in record:
        raise MalformedTimingRecordError(f"{where} has no {key!r}.")
    return record[key]


def number(value: object, where: str) -> float:
    """A finite JSON number, or a decimal string (``gate_vectors`` renders its seconds as strings).

    Raises:
        MalformedTimingRecordError: ``value`` is neither.
    """
    if isinstance(value, str):
        try:
            parsed = float(value)
        except ValueError as error:
            raise MalformedTimingRecordError(f"{where} is a number, not {value!r}.") from error
    elif isinstance(value, bool) or not isinstance(value, (int, float)):
        raise MalformedTimingRecordError(f"{where} is a number, not {value!r}.")
    else:
        parsed = float(value)
    if not math.isfinite(parsed):
        raise MalformedTimingRecordError(f"{where} is finite, not {value!r}.")
    return parsed


def numbers(value: object, where: str) -> list[float]:
    """Raises:
    MalformedTimingRecordError: ``value`` is not a list of numbers.
    """
    return [number(item, f"{where}[{index}]") for index, item in enumerate(items(value, where))]


def text(value: object, where: str) -> str:
    """Raises:
    MalformedTimingRecordError: ``value`` is not a string.
    """
    if not isinstance(value, str):
        raise MalformedTimingRecordError(f"{where} is a string, not {value!r}.")
    return value


def optional_text(record: Mapping[str, Any], key: str, where: str) -> Optional[str]:
    """The string at ``key``, or None when absent or null.

    Raises:
        MalformedTimingRecordError: present and not a string.
    """
    value = record.get(key)
    return None if value is None else text(value, f"{where}.{key}")


def optional_bool(record: Mapping[str, Any], key: str, where: str) -> Optional[bool]:
    """Raises:
    MalformedTimingRecordError: present and not a boolean.
    """
    value = record.get(key)
    if value is not None and not isinstance(value, bool):
        raise MalformedTimingRecordError(f"{where}.{key} is a boolean, not {value!r}.")
    return value


def integer(value: object, where: str) -> int:
    """Raises:
    MalformedTimingRecordError: ``value`` is not a JSON integer.
    """
    if isinstance(value, bool) or not isinstance(value, int):
        raise MalformedTimingRecordError(f"{where} is an integer, not {value!r}.")
    return value
