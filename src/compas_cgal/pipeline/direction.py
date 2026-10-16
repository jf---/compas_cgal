"""A dimensionless planar direction, frame-tagged: the unit lives in the type, not in a comment."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Generic
from typing import Self
from typing import overload

from compas_cgal.units import DirectionComponent
from compas_cgal.units import FrameT
from compas_cgal.units import _components


@dataclass(frozen=True)
class Direction2(Generic[FrameT]):
    """A planar direction in frame *FrameT*: dimensionless components, not millimetres.

    The components are a double construction (a normalised difference or a
    rotation of one); they place geometry and decide nothing.

    Attributes:
        x: X component.
        y: Y component.
    """

    x: DirectionComponent
    y: DirectionComponent

    @classmethod
    @overload
    def build(cls, x: float, y: float, /) -> Self: ...

    @classmethod
    @overload
    def build(cls, components: Sequence[float], /) -> Self: ...

    @classmethod
    def build(cls, *args: object) -> Self:  # pyright: ignore[reportInconsistentOverload]
        """A direction from two finite components or one two-component sequence.

        Raises:
            InvalidUnitValueError: The arguments are not two finite scalars.
        """
        x, y = _components(args, cls.__name__)
        return cls(DirectionComponent(x), DirectionComponent(y))

    def reversed(self) -> Direction2[FrameT]:
        """The opposite direction (exact negation)."""
        return Direction2[FrameT](DirectionComponent(-self.x), DirectionComponent(-self.y))
