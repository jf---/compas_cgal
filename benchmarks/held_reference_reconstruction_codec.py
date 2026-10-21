"""Lossless JSON codec for a certified :class:`ReferenceReconstruction`.

The reconstruction carries binary64 values only (``Point2`` coordinates,
``Radian`` sweeps, the ``Millimetre`` deviation bound). Each is written as
``float.hex()``: the exact binary64 value, including the sign of zero, which
a numerator/denominator pair would lose (``(-0.0).as_integer_ratio() == (0, 1)``)
and which ``atan2`` downstream distinguishes.
"""

from __future__ import annotations

import math
from collections.abc import Mapping

from benchmarks.errors import BenchmarkError
from benchmarks.held_reference_geometry import ReferenceArc
from benchmarks.held_reference_geometry import ReferenceLine
from benchmarks.held_reference_geometry import ReferencePrimitive
from benchmarks.held_reference_geometry import ReferenceReconstruction
from benchmarks.held_reference_reconstruction_cache_errors import MalformedReconstructionEncodingError
from compas_cgal.errors import InvalidUnitValueError
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import WorldXY

EncodedFloat = str
EncodedPoint = list[EncodedFloat]
EncodedPrimitive = dict[str, "EncodedPoint | EncodedFloat"]
EncodedReconstruction = dict[str, "list[EncodedPrimitive] | EncodedFloat"]

_LINE_KEYS = frozenset({"kind", "start", "end"})
_ARC_KEYS = frozenset({"kind", "start", "end", "centre", "sweep"})
_RECONSTRUCTION_KEYS = frozenset({"primitives", "deviation_upper_bound"})


def _encode_float(value: float) -> EncodedFloat:
    return float(value).hex()


def _encode_point(point: Point2[WorldXY]) -> EncodedPoint:
    return [_encode_float(point.x), _encode_float(point.y)]


def _encode_primitive(primitive: ReferencePrimitive) -> EncodedPrimitive:
    if isinstance(primitive, ReferenceLine):
        return {"kind": "line", "start": _encode_point(primitive.start), "end": _encode_point(primitive.end)}
    return {
        "kind": "arc",
        "start": _encode_point(primitive.start),
        "end": _encode_point(primitive.end),
        "centre": _encode_point(primitive.centre),
        "sweep": _encode_float(primitive.sweep),
    }


def encode_reconstruction(reconstruction: ReferenceReconstruction) -> EncodedReconstruction:
    """Encode a reconstruction as JSON-ready data, bit for bit.

    Args:
        reconstruction: The certified reconstruction to encode.

    Returns:
        A JSON-serializable mapping whose decoding is bitwise equal to the input.
    """
    return {
        "primitives": [_encode_primitive(primitive) for primitive in reconstruction.primitives],
        "deviation_upper_bound": _encode_float(reconstruction.deviation_upper_bound),
    }


def _decode_float(value: object, name: str) -> float:
    if not isinstance(value, str):
        raise MalformedReconstructionEncodingError(f"{name} must be a float.hex() string, got {type(value).__name__}.")
    try:
        decoded = float.fromhex(value)
    except ValueError as error:
        raise MalformedReconstructionEncodingError(f"{name} is not a float.hex() string: {value!r}.") from error
    if not math.isfinite(decoded) or decoded.hex() != value:
        raise MalformedReconstructionEncodingError(f"{name} is not the canonical hex of a finite binary64: {value!r}.")
    return decoded


def _decode_point(value: object, name: str) -> Point2[WorldXY]:
    if not isinstance(value, list) or len(value) != 2:
        raise MalformedReconstructionEncodingError(f"{name} must be a two-element list.")
    return Point2[WorldXY].build(_decode_float(value[0], f"{name}.x"), _decode_float(value[1], f"{name}.y"))


def _mapping(value: object, keys: frozenset[str], name: str) -> Mapping[str, object]:
    if not isinstance(value, dict) or set(value) != keys:
        raise MalformedReconstructionEncodingError(f"{name} must be an object with exactly the keys {sorted(keys)}.")
    return value


def _decode_primitive(value: object, index: int) -> ReferencePrimitive:
    name = f"primitives[{index}]"
    kind = value.get("kind") if isinstance(value, dict) else None
    if kind == "line":
        line = _mapping(value, _LINE_KEYS, name)
        return ReferenceLine.build(_decode_point(line["start"], f"{name}.start"), _decode_point(line["end"], f"{name}.end"))
    if kind == "arc":
        arc = _mapping(value, _ARC_KEYS, name)
        return ReferenceArc.build(
            _decode_point(arc["start"], f"{name}.start"),
            _decode_point(arc["end"], f"{name}.end"),
            _decode_point(arc["centre"], f"{name}.centre"),
            Radian(_decode_float(arc["sweep"], f"{name}.sweep")),
        )
    raise MalformedReconstructionEncodingError(f"{name} has unknown kind {kind!r}.")


def decode_reconstruction(value: object) -> ReferenceReconstruction:
    """Decode data written by :func:`encode_reconstruction`.

    Every primitive is rebuilt through its validating factory, so a decoded
    value satisfies the same invariants as a freshly computed one.

    Args:
        value: Parsed JSON data.

    Returns:
        The reconstruction, bitwise equal to the encoded one.

    Raises:
        MalformedReconstructionEncodingError: The data is not a valid encoding.
    """
    document = _mapping(value, _RECONSTRUCTION_KEYS, "reconstruction")
    primitives = document["primitives"]
    if not isinstance(primitives, list):
        raise MalformedReconstructionEncodingError("reconstruction.primitives must be a list.")
    try:
        return ReferenceReconstruction.build(
            tuple(_decode_primitive(primitive, index) for index, primitive in enumerate(primitives)),
            Millimetre(_decode_float(document["deviation_upper_bound"], "deviation_upper_bound")),
        )
    except (BenchmarkError, InvalidUnitValueError) as error:
        if isinstance(error, MalformedReconstructionEncodingError):
            raise
        raise MalformedReconstructionEncodingError(f"encoded reconstruction violates an invariant: {error}") from error
