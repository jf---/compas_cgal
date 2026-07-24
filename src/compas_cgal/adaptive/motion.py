import struct
from dataclasses import dataclass
from dataclasses import field
from fractions import Fraction
from typing import Self

from compas_cgal import _stock_2
from compas_cgal.adaptive.errors import DegenerateCircleMotionError
from compas_cgal.adaptive.errors import DegenerateSegmentMotionError
from compas_cgal.adaptive.errors import InvalidEngagementCapError
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Radian
from compas_cgal.adaptive.units import SquaredMillimetre
from compas_cgal.adaptive.units import Vector2
from compas_cgal.adaptive.units import WorldXY


@dataclass(frozen=True)
class ExactSegmentMotion:
    start: Point2[WorldXY]
    end: Point2[WorldXY]

    def __post_init__(self) -> None:
        if not isinstance(self.start, Point2) or not isinstance(self.end, Point2):
            raise DegenerateSegmentMotionError("segment endpoints must be world-XY points.")
        if self.start == self.end:
            raise DegenerateSegmentMotionError("segment motion requires nonzero exact progress.")

    @classmethod
    def build(cls, start: Point2[WorldXY], end: Point2[WorldXY]) -> Self:
        return cls(start, end)


@dataclass(frozen=True)
class ExactCircleMotion:
    center: Point2[WorldXY]
    phase_vector: Vector2[WorldXY]
    clockwise: bool

    def __post_init__(self) -> None:
        if not isinstance(self.center, Point2) or not isinstance(self.phase_vector, Vector2):
            raise DegenerateCircleMotionError("circle motion requires a world-XY center and phase vector.")
        if self.phase_vector.x == 0.0 and self.phase_vector.y == 0.0:
            raise DegenerateCircleMotionError("circle motion requires a nonzero exact phase vector.")
        if type(self.clockwise) is not bool:
            raise DegenerateCircleMotionError("circle orientation must be explicitly clockwise or counterclockwise.")

    @classmethod
    def build(
        cls,
        center: Point2[WorldXY],
        phase_vector: Vector2[WorldXY],
        clockwise: bool,
    ) -> Self:
        return cls(center, phase_vector, clockwise)

    @property
    def squared_radius(self) -> SquaredMillimetre:
        dx = Fraction.from_float(self.phase_vector.x)
        dy = Fraction.from_float(self.phase_vector.y)
        return SquaredMillimetre(dx * dx + dy * dy)


@dataclass(frozen=True)
class EngagementCap:
    theta: Radian
    chord_ratio: float = field(init=False)
    chord_ratio_bytes: bytes = field(init=False)

    def __post_init__(self) -> None:
        if isinstance(self.theta, bool) or not isinstance(self.theta, (int, float)):
            raise InvalidEngagementCapError("engagement cap must be a real angle.")
        try:
            theta = float(self.theta)
        except OverflowError:
            raise InvalidEngagementCapError("engagement cap exceeds the binary64 range.") from None
        try:
            native_ratio = _stock_2.cap_chord_ratio(theta)
        except ValueError as error:
            raise InvalidEngagementCapError(str(error)) from None
        object.__setattr__(self, "theta", Radian(theta))
        object.__setattr__(self, "chord_ratio", native_ratio)
        object.__setattr__(self, "chord_ratio_bytes", struct.pack(">d", native_ratio))

    @classmethod
    def build(cls, theta: float) -> Self:
        return cls(Radian(theta))
