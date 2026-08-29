"""Strict JSON decoding for measurement-claim artifacts."""

from __future__ import annotations

import json
import math
from typing import Dict
from typing import List
from typing import TypeVar
from typing import cast

ErrorT = TypeVar("ErrorT", bound=Exception)


def validate_finite(value: object, field: str, error: type[ErrorT]) -> None:
    """Reject non-finite floats recursively without changing the decoded tree."""
    if type(value) is float and not math.isfinite(value):
        raise error(f"{field}: must contain only finite numbers")
    if type(value) is list:
        for index, item in enumerate(cast(List[object], value)):
            validate_finite(item, f"{field}[{index}]", error)
    elif type(value) is dict:
        for key, item in cast(Dict[str, object], value).items():
            validate_finite(item, f"{field}.{key}", error)


def decode_strict(data: bytes, field: str, error: type[ErrorT]) -> object:
    """Decode duplicate-free, standards-conforming UTF-8 JSON."""

    def pairs(items: List[tuple[str, object]]) -> Dict[str, object]:
        result: Dict[str, object] = {}
        for key, value in items:
            if key in result:
                raise error(f"duplicate JSON key: {key}")
            result[key] = value
        return result

    def constant(value: str) -> object:
        raise error(f"non-standard JSON constant: {value}")

    try:
        value = json.loads(data.decode("utf-8"), object_pairs_hook=pairs, parse_constant=constant)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise error(f"{field}: not strict UTF-8 JSON") from exc
    validate_finite(value, field, error)
    return value
