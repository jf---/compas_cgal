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

from compas_cgal.errors import InvalidUnitValueError

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
    # Fast path: a finite exact ``float`` is already its own validated value
    # (``float(value)`` would return it unchanged). Everything else, including
    # float subclasses, ints and the error cases, takes the full check below.
    if type(value) is float and math.isfinite(value):
        return value
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


# The four vector units' ``__post_init__`` writes a component back only when
# validation changed it (an int widened to float): ``_finite`` returns a finite
# ``float`` itself, and rewriting the same object is a no-op that costs a quarter
# of the construction.


@dataclass(frozen=True)
class Point2(Generic[FrameT]):
    x: Millimetre
    y: Millimetre

    def __post_init__(self) -> None:
        x = _finite(self.x, "Point2.x")
        if x is not self.x:
            object.__setattr__(self, "x", Millimetre(x))
        y = _finite(self.y, "Point2.y")
        if y is not self.y:
            object.__setattr__(self, "y", Millimetre(y))

    @classmethod
    @overload
    def build(cls, x: float, y: float, /) -> Self: ...

    @classmethod
    @overload
    def build(cls, components: Sequence[float], /) -> Self: ...

    @classmethod
    def build(cls, *args: object) -> Self:  # pyright: ignore[reportInconsistentOverload]
        if len(args) == 2:
            x, y = args
            if type(x) is float and type(y) is float and math.isfinite(x) and math.isfinite(y):
                if cls is Point2:
                    # Two exact finite doubles are validated here, and are exactly
                    # what ``__post_init__`` would keep: fill the frozen slots
                    # directly (the hot path, ~10^5 points per Figure 5 run).
                    point = object.__new__(cls)
                    object.__setattr__(point, "x", x)
                    object.__setattr__(point, "y", y)
                    return point
                # A subclass may extend ``__post_init__``: it is the one validation.
                return cls(Millimetre(x), Millimetre(y))
        x, y = _components(args, cls.__name__)
        return cls(Millimetre(x), Millimetre(y))


@dataclass(frozen=True)
class Vector2(Generic[FrameT]):
    x: Millimetre
    y: Millimetre

    def __post_init__(self) -> None:
        x = _finite(self.x, "Vector2.x")
        if x is not self.x:
            object.__setattr__(self, "x", Millimetre(x))
        y = _finite(self.y, "Vector2.y")
        if y is not self.y:
            object.__setattr__(self, "y", Millimetre(y))

    @classmethod
    @overload
    def build(cls, x: float, y: float, /) -> Self: ...

    @classmethod
    @overload
    def build(cls, components: Sequence[float], /) -> Self: ...

    @classmethod
    def build(cls, *args: object) -> Self:  # pyright: ignore[reportInconsistentOverload]
        if len(args) == 2:
            x, y = args
            # Two finite doubles: ``__post_init__`` is the one validation.
            if type(x) is float and type(y) is float and math.isfinite(x) and math.isfinite(y):
                return cls(Millimetre(x), Millimetre(y))
        x, y = _components(args, cls.__name__)
        return cls(Millimetre(x), Millimetre(y))


@dataclass(frozen=True)
class Point3(Generic[FrameT]):
    x: Millimetre
    y: Millimetre
    z: Millimetre

    def __post_init__(self) -> None:
        x = _finite(self.x, "Point3.x")
        if x is not self.x:
            object.__setattr__(self, "x", Millimetre(x))
        y = _finite(self.y, "Point3.y")
        if y is not self.y:
            object.__setattr__(self, "y", Millimetre(y))
        z = _finite(self.z, "Point3.z")
        if z is not self.z:
            object.__setattr__(self, "z", Millimetre(z))

    @classmethod
    @overload
    def build(cls, x: float, y: float, z: float, /) -> Self: ...

    @classmethod
    @overload
    def build(cls, components: Sequence[float], /) -> Self: ...

    @classmethod
    def build(cls, *args: object) -> Self:  # pyright: ignore[reportInconsistentOverload]
        if len(args) == 3:
            x, y, z = args
            # Three finite doubles: ``__post_init__`` is the one validation.
            if type(x) is float and type(y) is float and type(z) is float and math.isfinite(x) and math.isfinite(y) and math.isfinite(z):
                return cls(Millimetre(x), Millimetre(y), Millimetre(z))
        x, y, z = _components3(args, cls.__name__)
        return cls(Millimetre(x), Millimetre(y), Millimetre(z))


@dataclass(frozen=True)
class Direction3(Generic[FrameT]):
    x: DirectionComponent
    y: DirectionComponent
    z: DirectionComponent

    def __post_init__(self) -> None:
        x = _finite(self.x, "Direction3.x")
        if x is not self.x:
            object.__setattr__(self, "x", DirectionComponent(x))
        y = _finite(self.y, "Direction3.y")
        if y is not self.y:
            object.__setattr__(self, "y", DirectionComponent(y))
        z = _finite(self.z, "Direction3.z")
        if z is not self.z:
            object.__setattr__(self, "z", DirectionComponent(z))

    @classmethod
    @overload
    def build(cls, x: float, y: float, z: float, /) -> Self: ...

    @classmethod
    @overload
    def build(cls, components: Sequence[float], /) -> Self: ...

    @classmethod
    def build(cls, *args: object) -> Self:  # pyright: ignore[reportInconsistentOverload]
        if len(args) == 3:
            x, y, z = args
            # Three finite doubles: ``__post_init__`` is the one validation.
            if type(x) is float and type(y) is float and type(z) is float and math.isfinite(x) and math.isfinite(y) and math.isfinite(z):
                return cls(DirectionComponent(x), DirectionComponent(y), DirectionComponent(z))
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
        # ``__post_init__`` is the one validation (same check, same message).
        return cls(Millimetre(value))


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


@dataclass(frozen=True)
class CutDepth:
    """Depth of one cut level below the stock top, along world -Z, mm (positive).

    Frame: world Z with the stock top as datum (`CutDepth.cut_plane`).
    """

    value: Millimetre

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", Millimetre(_positive(self.value, "cut depth")))

    @classmethod
    def build(cls, value: float) -> Self:
        return cls(Millimetre(_positive(value, "cut depth")))

    def cut_plane(self, stock_top: ClearanceZ) -> CutPlane:
        """The cut plane of this depth under *stock_top*: cut at ``stock_top - depth``, a ramp starting at *stock_top*.

        The subtraction is one double operation; every downstream decision
        reads the resulting doubles exactly.

        Args:
            stock_top: The stock top, where a ramp into the stock starts.

        Returns:
            The cut plane.
        """
        return CutPlane(CutZ(Millimetre(stock_top.value - self.value)), stock_top)
