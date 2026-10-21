"""Native declared-boundary clearance oracle: batch calls into ``_held_declared_clearance_2``.

The C++ oracle (src/held_declared_clearance_2/) decides exactly what
:class:`benchmarks.held_declared_clearance.DeclaredClearance` decides, with an
AABB tree over the declared primitives and Epeck / Sqrt_extension predicates,
and selects the reported primitive the same way. Every call is one batch: no
per-element Python crossing.

Owner-by-construction TOUCH: an exact element may name the declared LINE or
ARC primitives it offsets (a ∂Ω* piece's owners). The native oracle accepts an
owner only where a stated lemma proves ``dist == r`` and checks every hypothesis
exactly (src/held_declared_clearance_2/owner_lemma_2.h); every other pair, and
every non-owner primitive, is decided in full.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import Self

from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_declared_boundary import DeclaredLine
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_clearance import CircleElement
from benchmarks.held_declared_clearance import ClearanceVerdict
from benchmarks.held_declared_clearance import Element
from benchmarks.held_declared_clearance_wire import binary64_circle
from benchmarks.held_declared_clearance_wire import binary64_point
from benchmarks.held_declared_clearance_wire import element_wire
from benchmarks.held_declared_clearance_wire import primitive_wires
from benchmarks.held_declared_clearance_wire import rational_wire
from benchmarks.held_gouge_audit import Clearance
from benchmarks.held_gouge_boundary import PrimitiveKind
from compas_cgal import _held_declared_clearance_2
from compas_cgal.units import WorldXY


@dataclass(frozen=True)
class NativeCounters:
    """Reported work of the last batch: elements, pruned-in pairs, owner-lemma touches, exact crossings."""

    elements: int
    pairs: int
    owner_touches: int
    crossings: int


def _verdict(wire: tuple[int, int, float]) -> ClearanceVerdict:
    clearance, primitive, depth = wire
    return ClearanceVerdict(Clearance(clearance), None if primitive < 0 else primitive, depth)


@dataclass(frozen=True)
class NativeDeclaredClearance:
    """The native oracle over one declared boundary and one tool radius."""

    boundary: DeclaredBoundary
    tool_radius: Fraction
    handle: _held_declared_clearance_2.DeclaredClearance2

    @classmethod
    def build(cls, boundary: DeclaredBoundary, tool_radius: Fraction) -> Self:
        """Index the declared primitives natively.

        Raises:
            NonPositiveClearanceRadiusError: If ``tool_radius`` is not a positive binary64 (native).
        """
        handle = _held_declared_clearance_2.DeclaredClearance2.build(primitive_wires(boundary), rational_wire(tool_radius))
        return cls(boundary, Fraction(tool_radius), handle)

    @property
    def primitive_count(self) -> int:
        return self.handle.primitive_count

    def circles(self, elements: Sequence[CircleElement]) -> list[ClearanceVerdict]:
        """Verdicts of binary64 machining circles, one batch.

        Raises:
            NonBinary64ClearanceInputError: If a centre or radius is not a binary64.
        """
        values = [binary64_circle(element) for element in elements]
        wires = self.handle.circles_binary64([x for x, _, _ in values], [y for _, y, _ in values], [radius for _, _, radius in values])
        return [_verdict(wire) for wire in wires]

    def segments_binary64(self, segments: Sequence[tuple[ExactPoint2[WorldXY], ExactPoint2[WorldXY]]]) -> list[ClearanceVerdict]:
        """Verdicts of polyline pieces between distinct binary64 points (``SegmentElement.between``), one batch.

        Raises:
            NonBinary64ClearanceInputError: If a coordinate is not a binary64.
        """
        starts = [binary64_point(start) for start, _ in segments]
        ends = [binary64_point(end) for _, end in segments]
        wires = self.handle.segments_binary64([x for x, _ in starts], [y for _, y in starts], [x for x, _ in ends], [y for _, y in ends])
        return [_verdict(wire) for wire in wires]

    def elements(self, elements: Sequence[Element], owners: Sequence[Sequence[int]] | None = None) -> list[ClearanceVerdict]:
        """Verdicts of exact segment and arc elements, one batch; ``owners[i]`` are element ``i``'s declared owners."""
        listed: Sequence[Sequence[int]] = owners if owners is not None else [()] * len(elements)
        wires = self.handle.elements([element_wire(element, owned) for element, owned in zip(elements, listed, strict=True)])
        return [_verdict(wire) for wire in wires]

    def primitive_kind(self, index: int) -> PrimitiveKind:
        """``"line"``, ``"convex-arc"`` or ``"concave-arc"`` (the gouge audit's vocabulary)."""
        primitive = self.boundary.primitives[index]
        if isinstance(primitive, DeclaredLine):
            return "line"
        return "convex-arc" if primitive.sense == 1 else "concave-arc"


def last_counters() -> NativeCounters:
    """The reported work of the last native batch."""
    return NativeCounters(*_held_declared_clearance_2.last_counters())
