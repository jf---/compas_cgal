"""Pin the two-cap operation-stream premise used by the quality product gate.

Three pockets by two generators is the smallest product in which a defect can
be attributed. At the default cap both generators intentionally emit the same
stream; at the attribution cap each pocket must expose a generator difference.
"""

from __future__ import annotations

import struct
import warnings
from typing import Iterable
from typing import Literal
from typing import NewType
from typing import Optional
from typing import Tuple
from typing import Union
from typing import cast

import pytest
from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Line

from benchmarks.gate import GATE_ATTRIBUTION_CAP_DEG
from benchmarks.gate import GATE_CAP_DEG
from benchmarks.gate import GATE_GENERATOR_NAMES
from benchmarks.gate import GATE_GENERATORS
from benchmarks.gate import GATE_POCKET_NAMES
from benchmarks.gate import GateCapDegrees
from benchmarks.gate import gate_pocket
from compas_cgal.engagement_toolpath import UnavoidableEngagementWarning
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult


class UnsupportedGateOperationGeometryError(TypeError):
    """The C3 witness has no defining-geometry contract for this primitive."""


class MalformedGateOperationWitnessError(ValueError):
    """A C3 geometry coordinate or tangent does not have three components."""


Binary64Bits = NewType("Binary64Bits", bytes)
Binary64Vector3 = Tuple[Binary64Bits, Binary64Bits, Binary64Bits]
LineGeometryWitness = Tuple[Literal["Line"], Binary64Vector3, Binary64Vector3]
CircleGeometryWitness = Tuple[
    Literal["Circle"],
    Binary64Vector3,
    Binary64Vector3,
    Binary64Vector3,
    Binary64Bits,
]
GateGeometryWitness = Union[LineGeometryWitness, CircleGeometryWitness]
GateOperationWitness = Tuple[
    GateGeometryWitness,
    str,
    str,
    int,
    bool,
    Optional[Binary64Vector3],
    Optional[Binary64Vector3],
]
GateOperationStreamWitness = Tuple[GateOperationWitness, ...]


def _binary64_bits(value: float) -> Binary64Bits:
    """Return the network-order binary64 representation of one scalar."""
    return Binary64Bits(struct.pack("!d", float(value)))


def _binary64_vector3(values: Iterable[float], *, field: str) -> Binary64Vector3:
    """Pack one three-component coordinate or tangent without normalization."""
    components = tuple(values)
    if len(components) != 3:
        raise MalformedGateOperationWitnessError(f"C3 witness requires exactly three {field} components; got {len(components)}.")
    return (
        _binary64_bits(components[0]),
        _binary64_bits(components[1]),
        _binary64_bits(components[2]),
    )


def _qualified_type_name(value: object) -> str:
    """Return a stable diagnostic name for an unsupported runtime type."""
    value_type = type(value)
    return f"{value_type.__module__}.{value_type.__qualname__}"


def _unsupported_geometry(geometry: object) -> UnsupportedGateOperationGeometryError:
    """Build the named fail-loud boundary error for C3's closed witness."""
    return UnsupportedGateOperationGeometryError(f"C3's deliberately bounded Line/Circle contract does not support {_qualified_type_name(geometry)}.")


def _geometry_witness(geometry: Union[Line, Arc, Circle]) -> GateGeometryWitness:
    """Witness defining Line/Circle geometry; reject every other exact type."""
    if type(geometry) is Line:
        line_tag: Literal["Line"] = "Line"
        return (
            line_tag,
            _binary64_vector3(geometry.start, field="line-start coordinate"),
            _binary64_vector3(geometry.end, field="line-end coordinate"),
        )
    if type(geometry) is Circle:
        circle_tag: Literal["Circle"] = "Circle"
        return (
            circle_tag,
            _binary64_vector3(geometry.frame.point, field="circle-frame point"),
            _binary64_vector3(geometry.frame.xaxis, field="circle-frame x-axis"),
            _binary64_vector3(geometry.frame.yaxis, field="circle-frame y-axis"),
            _binary64_bits(geometry.radius),
        )
    if type(geometry) is Arc:
        raise _unsupported_geometry(geometry)
    raise _unsupported_geometry(geometry)


def _optional_tangent_witness(
    tangent: Optional[Iterable[float]],
    *,
    field: str,
) -> Optional[Binary64Vector3]:
    """Pack one optional tangent without truncating or normalizing it."""
    if tangent is None:
        return None
    return _binary64_vector3(tangent, field=field)


def _operation_witness(operation: ToolpathOperation) -> GateOperationWitness:
    """Witness one emitted operation's defining geometry and metadata."""
    return (
        _geometry_witness(operation.geometry),
        type(operation.operation).__qualname__,
        operation.operation.value,
        operation.path_index,
        operation.clockwise,
        _optional_tangent_witness(operation.start_tangent, field="start-tangent"),
        _optional_tangent_witness(operation.end_tangent, field="end-tangent"),
    )


def _operation_stream_witness(result: ToolpathResult) -> GateOperationStreamWitness:
    """Witness every operation in emission order, retaining duplicates."""
    return tuple(_operation_witness(operation) for operation in result.operations)


def _generated_stream(
    generator_name: str,
    pocket_name: str,
    tea_cap_deg: GateCapDegrees,
) -> GateOperationStreamWitness:
    """Generate and witness one fresh gate instance under exact warning policy."""
    spec = gate_pocket(pocket_name, tea_cap_deg=tea_cap_deg)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        warnings.filterwarnings("ignore", category=UnavoidableEngagementWarning)
        result = GATE_GENERATORS[generator_name](spec)
    return _operation_stream_witness(result)


_FIRST_GENERATOR, _SECOND_GENERATOR = GATE_GENERATOR_NAMES


@pytest.mark.parametrize("pocket_name", GATE_POCKET_NAMES, ids=GATE_POCKET_NAMES)
def test_default_cap_operation_streams_are_identical(pocket_name: str) -> None:
    """Catch a default-cap generator change that invalidates saturated attribution."""
    first = _generated_stream(_FIRST_GENERATOR, pocket_name, GATE_CAP_DEG)
    second = _generated_stream(_SECOND_GENERATOR, pocket_name, GATE_CAP_DEG)
    assert first == second


@pytest.mark.parametrize("pocket_name", GATE_POCKET_NAMES, ids=GATE_POCKET_NAMES)
def test_attribution_cap_operation_streams_diverge(pocket_name: str) -> None:
    """Catch a cap-40 change that leaves the radius ladder inert on one pocket."""
    first = _generated_stream(_FIRST_GENERATOR, pocket_name, GATE_ATTRIBUTION_CAP_DEG)
    second = _generated_stream(_SECOND_GENERATOR, pocket_name, GATE_ATTRIBUTION_CAP_DEG)
    assert first != second


def test_arc_geometry_is_rejected_by_named_error() -> None:
    """Catch accidental expansion of C3's deliberately closed primitive contract."""
    operation = ToolpathOperation(
        geometry=Arc(radius=1.0, start_angle=0.0, end_angle=1.0),
        operation=OperationType.CUT,
        path_index=0,
    )
    with pytest.raises(UnsupportedGateOperationGeometryError) as raised:
        _operation_witness(operation)
    assert "Arc" in str(raised.value)
    assert "Line/Circle" in str(raised.value)


def test_unknown_geometry_is_rejected_by_named_error() -> None:
    """Catch fallback serialization of a primitive outside the closed witness."""
    operation = ToolpathOperation(
        geometry=cast(Union[Line, Arc, Circle], object()),
        operation=OperationType.CUT,
        path_index=0,
    )
    with pytest.raises(UnsupportedGateOperationGeometryError) as raised:
        _operation_witness(operation)
    assert "builtins.object" in str(raised.value)
