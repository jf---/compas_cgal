import math
import struct
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import Generic
from typing import NewType
from typing import Self
from typing import TypeVar
from typing import overload

from compas_cgal.adaptive.errors import InvalidUnitValueError

Millimetre = NewType("Millimetre", float)
DirectionComponent = NewType("DirectionComponent", float)
ExactMillimetre = NewType("ExactMillimetre", Fraction)
SquaredMillimetre = NewType("SquaredMillimetre", Fraction)
Radian = NewType("Radian", float)


class WorldXY:
    """Phantom frame tag for world-XY geometry."""


class WorldXYZ:
    """Phantom frame tag for world-XYZ ingress geometry."""


FrameT = TypeVar("FrameT")


def _finite(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise InvalidUnitValueError(f"{name} must be a finite real number.")
    try:
        numeric = float(value)
    except OverflowError:
        raise InvalidUnitValueError(f"{name} exceeds the binary64 range.") from None
    if not math.isfinite(numeric):
        raise InvalidUnitValueError(f"{name} must be finite.")
    return numeric


def _positive(value: object, name: str) -> float:
    numeric = _finite(value, name)
    if numeric <= 0.0:
        raise InvalidUnitValueError(f"{name} must be positive.")
    return numeric


def _non_negative(value: object, name: str) -> float:
    numeric = _finite(value, name)
    if numeric < 0.0:
        raise InvalidUnitValueError(f"{name} must be non-negative.")
    return numeric


def _components(args: tuple[object, ...], name: str) -> tuple[float, float]:
    if len(args) == 1 and isinstance(args[0], Sequence) and not isinstance(args[0], (str, bytes)):
        values = args[0]
        if len(values) != 2:
            raise InvalidUnitValueError(f"{name} requires exactly two components.")
        return _finite(values[0], f"{name}.x"), _finite(values[1], f"{name}.y")
    if len(args) == 2:
        return _finite(args[0], f"{name}.x"), _finite(args[1], f"{name}.y")
    raise InvalidUnitValueError(f"{name} requires two scalars or one two-component sequence.")


def _components3(args: tuple[object, ...], name: str) -> tuple[float, float, float]:
    if len(args) == 1 and isinstance(args[0], Sequence) and not isinstance(args[0], (str, bytes)):
        values = args[0]
        if len(values) != 3:
            raise InvalidUnitValueError(f"{name} requires exactly three components.")
        return (
            _finite(values[0], f"{name}.x"),
            _finite(values[1], f"{name}.y"),
            _finite(values[2], f"{name}.z"),
        )
    if len(args) == 3:
        return (
            _finite(args[0], f"{name}.x"),
            _finite(args[1], f"{name}.y"),
            _finite(args[2], f"{name}.z"),
        )
    raise InvalidUnitValueError(f"{name} requires three scalars or one three-component sequence.")


@dataclass(frozen=True)
class Point2(Generic[FrameT]):
    x: Millimetre
    y: Millimetre

    def __post_init__(self) -> None:
        object.__setattr__(self, "x", Millimetre(_finite(self.x, "Point2.x")))
        object.__setattr__(self, "y", Millimetre(_finite(self.y, "Point2.y")))

    @classmethod
    @overload
    def build(cls, x: float, y: float, /) -> Self: ...

    @classmethod
    @overload
    def build(cls, components: Sequence[float], /) -> Self: ...

    @classmethod
    def build(cls, *args: object) -> Self:  # pyright: ignore[reportInconsistentOverload]
        x, y = _components(args, cls.__name__)
        return cls(Millimetre(x), Millimetre(y))


@dataclass(frozen=True)
class Vector2(Generic[FrameT]):
    x: Millimetre
    y: Millimetre

    def __post_init__(self) -> None:
        object.__setattr__(self, "x", Millimetre(_finite(self.x, "Vector2.x")))
        object.__setattr__(self, "y", Millimetre(_finite(self.y, "Vector2.y")))

    @classmethod
    @overload
    def build(cls, x: float, y: float, /) -> Self: ...

    @classmethod
    @overload
    def build(cls, components: Sequence[float], /) -> Self: ...

    @classmethod
    def build(cls, *args: object) -> Self:  # pyright: ignore[reportInconsistentOverload]
        x, y = _components(args, cls.__name__)
        return cls(Millimetre(x), Millimetre(y))


@dataclass(frozen=True)
class Point3(Generic[FrameT]):
    x: Millimetre
    y: Millimetre
    z: Millimetre

    def __post_init__(self) -> None:
        object.__setattr__(self, "x", Millimetre(_finite(self.x, "Point3.x")))
        object.__setattr__(self, "y", Millimetre(_finite(self.y, "Point3.y")))
        object.__setattr__(self, "z", Millimetre(_finite(self.z, "Point3.z")))

    @classmethod
    @overload
    def build(cls, x: float, y: float, z: float, /) -> Self: ...

    @classmethod
    @overload
    def build(cls, components: Sequence[float], /) -> Self: ...

    @classmethod
    def build(cls, *args: object) -> Self:  # pyright: ignore[reportInconsistentOverload]
        x, y, z = _components3(args, cls.__name__)
        return cls(Millimetre(x), Millimetre(y), Millimetre(z))


@dataclass(frozen=True)
class Direction3(Generic[FrameT]):
    x: DirectionComponent
    y: DirectionComponent
    z: DirectionComponent

    def __post_init__(self) -> None:
        object.__setattr__(self, "x", DirectionComponent(_finite(self.x, "Direction3.x")))
        object.__setattr__(self, "y", DirectionComponent(_finite(self.y, "Direction3.y")))
        object.__setattr__(self, "z", DirectionComponent(_finite(self.z, "Direction3.z")))

    @classmethod
    @overload
    def build(cls, x: float, y: float, z: float, /) -> Self: ...

    @classmethod
    @overload
    def build(cls, components: Sequence[float], /) -> Self: ...

    @classmethod
    def build(cls, *args: object) -> Self:  # pyright: ignore[reportInconsistentOverload]
        x, y, z = _components3(args, cls.__name__)
        return cls(DirectionComponent(x), DirectionComponent(y), DirectionComponent(z))


@dataclass(frozen=True)
class ToolRadius:
    value: Millimetre

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", Millimetre(_positive(self.value, "tool radius")))

    @classmethod
    def build(cls, value: float) -> Self:
        return cls(Millimetre(_positive(value, "tool radius")))


@dataclass(frozen=True)
class EntryRadius:
    value: Millimetre

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", Millimetre(_positive(self.value, "entry radius")))

    @classmethod
    def build(cls, value: float) -> Self:
        return cls(Millimetre(_positive(value, "entry radius")))


@dataclass(frozen=True)
class Clearance:
    value: Millimetre

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", Millimetre(_non_negative(self.value, "clearance")))

    @classmethod
    def build(cls, value: float) -> Self:
        return cls(Millimetre(_non_negative(value, "clearance")))


@dataclass(frozen=True)
class GuideRadius:
    value: Millimetre

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", Millimetre(_positive(self.value, "guide radius")))

    @classmethod
    def build(cls, value: float) -> Self:
        return cls(Millimetre(_positive(value, "guide radius")))


@dataclass(frozen=True)
class ObservedRadius:
    """Finite millimetre observation awaiting native geometric validation."""

    value: Millimetre

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", Millimetre(_finite(self.value, "observed radius")))

    @classmethod
    def build(cls, value: float) -> Self:
        return cls(Millimetre(_finite(value, "observed radius")))


@dataclass(frozen=True)
class Spacing:
    value: Millimetre

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", Millimetre(_positive(self.value, "spacing")))

    @classmethod
    def build(cls, value: float) -> Self:
        return cls(Millimetre(_positive(value, "spacing")))


@dataclass(frozen=True)
class ChordBound:
    value: Millimetre

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", Millimetre(_positive(self.value, "chord bound")))

    @classmethod
    def build(cls, value: float) -> Self:
        return cls(Millimetre(_positive(value, "chord bound")))

    @property
    def exact_bytes(self) -> bytes:
        return struct.pack(">d", self.value)


@dataclass(frozen=True)
class CutZ:
    value: Millimetre

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", Millimetre(_finite(self.value, "cut Z")))

    @classmethod
    def build(cls, value: float) -> Self:
        return cls(Millimetre(_finite(value, "cut Z")))


@dataclass(frozen=True)
class ClearanceZ:
    value: Millimetre

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", Millimetre(_finite(self.value, "clearance Z")))

    @classmethod
    def build(cls, value: float) -> Self:
        return cls(Millimetre(_finite(value, "clearance Z")))


@dataclass(frozen=True)
class CutPlane:
    cut_z: CutZ
    clearance_z: ClearanceZ

    def __post_init__(self) -> None:
        if not isinstance(self.cut_z, CutZ) or not isinstance(self.clearance_z, ClearanceZ):
            raise InvalidUnitValueError("cut plane requires typed cut and clearance Z values.")
        if self.clearance_z.value <= self.cut_z.value:
            raise InvalidUnitValueError("clearance Z must be strictly above cut Z.")

    @classmethod
    def build(cls, cut_z: CutZ, clearance_z: ClearanceZ) -> Self:
        return cls(cut_z, clearance_z)
