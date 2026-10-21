"""Declared Held boundary: the Pythagorean (offset-rational) boundary contract.

The declared boundary replaces the stored float primitives as the authority for
exact Held geometry (docs/held_exact_centre_domain.md, "Boundary contract").
It sits beside :class:`benchmarks.held_reference_geometry.ReferenceBoundary`;
the stored representation is unchanged and consumers switch separately.

Every number is a :class:`~fractions.Fraction` reached from binary64 inputs by
field operations only:

- a direction is ``i^quarter · ((1 - t²), 2t) / (1 + t²)`` with ``t`` on the
  dyadic grid ``2^-56``, so it is an exactly unit rational vector;
- a line is ``start + λ·u``, ``λ ∈ [0, L]``, ``L`` a binary64;
- an arc is ``(S, E, start tangent T, sense σ)``; with ``N = σ·T^⊥`` the unit
  normal towards the centre, ``s = |E - S|² / (2 N·(E - S))`` gives the
  rational centre ``S + s·N`` and rational radius ``s``.

G1 is not imposed: a junction may be an exact kink of any size. Closure is
exact and is the one invariant :meth:`DeclaredBoundary.build` enforces.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable
from dataclasses import dataclass
from dataclasses import field
from fractions import Fraction
from typing import Generic
from typing import Literal
from typing import Self
from typing import TypeAlias
from typing import TypeVar

from benchmarks.errors import BenchmarkError
from compas_cgal.units import ExactMillimetre
from compas_cgal.units import SquaredMillimetre
from compas_cgal.units import WorldXY

# Grid of the half-angle tangent t: 2^-56 gives 2^-55 rad (2.8e-17) angular
# resolution, below binary64 direction error, and caps every direction
# denominator at 113 bits (measured; unrounded t reaches 211 bits).
T_GRID_BITS = 56
T_GRID_DENOMINATOR = 2**T_GRID_BITS
# |t| = |tan(phi / 2)| with |phi| <= pi/4 after the quarter-turn reduction, so
# |t| <= tan(pi/8) = 0.41421...; 1/2 is the declared rational cap (grid rounding
# adds at most 2^-57).
MAX_ABS_HALF_ANGLE_TANGENT = Fraction(1, 2)
QUARTER_TURNS = 4
# Version of the declared-boundary contract embedded in every identity hash.
DECLARED_BOUNDARY_CONTRACT_VERSION = 1

Sense: TypeAlias = Literal[1, -1]
FrameT = TypeVar("FrameT")
ValueT = TypeVar("ValueT")


@dataclass
class _Memo(Generic[ValueT]):
    """One derived value of a frozen owner, computed on first read.

    The owner holds exactly one per value, created in its ``__init__``
    (``init=False``, excluded from eq/hash/repr), so the frozen owner's own
    attributes are never written after construction and
    ``dataclasses.replace`` never shares a memo between instances.
    """

    value: ValueT | None = None

    def get(self, compute: Callable[[], ValueT]) -> ValueT:
        if self.value is None:
            self.value = compute()
        return self.value


class DeclaredBoundaryError(BenchmarkError):
    """Base error of the declared (Pythagorean) Held boundary."""


class InvalidPythagoreanDirectionError(DeclaredBoundaryError):
    """A direction's t is off the 2^-56 grid, exceeds the half-angle cap, or its quarter is not 0..3."""


class InvalidDeclaredLineError(DeclaredBoundaryError):
    """A declared line has a non-positive or non-binary64 length."""


class DegenerateArcChordError(DeclaredBoundaryError):
    """A declared arc's chord does not advance towards its centre (N·(E - S) <= 0)."""


class OpenDeclaredCycleError(DeclaredBoundaryError):
    """A declared primitive does not end exactly where its successor starts."""


class InvalidDeclaredProvenanceError(DeclaredBoundaryError):
    """Declared primitives do not map one-to-one, in order, onto the source primitive indices."""


@dataclass(frozen=True)
class ExactPoint2(Generic[FrameT]):
    """A point with exact rational millimetre coordinates in frame ``FrameT``."""

    x: ExactMillimetre
    y: ExactMillimetre

    @classmethod
    def build(cls, x: Fraction, y: Fraction) -> Self:
        """Build from two rationals.

        Args:
            x: Exact x coordinate in millimetres.
            y: Exact y coordinate in millimetres.

        Returns:
            The exact point.
        """
        return cls(ExactMillimetre(Fraction(x)), ExactMillimetre(Fraction(y)))

    @classmethod
    def from_binary64(cls, x: float, y: float) -> Self:
        """Inject two binary64 coordinates exactly (every double is a rational).

        Args:
            x: Binary64 x coordinate in millimetres.
            y: Binary64 y coordinate in millimetres.

        Returns:
            The exact point equal to the doubles.
        """
        return cls.build(Fraction(x), Fraction(y))

    def offset(self, direction: tuple[Fraction, Fraction], distance: Fraction) -> Self:
        """Return ``self + distance·direction`` exactly."""
        return type(self).build(self.x + distance * direction[0], self.y + distance * direction[1])

    def minus(self, other: ExactPoint2[FrameT]) -> tuple[Fraction, Fraction]:
        """Return the exact vector ``self - other``."""
        return Fraction(self.x - other.x), Fraction(self.y - other.y)

    def as_floats(self) -> tuple[float, float]:
        """Nearest binary64 coordinates (reporting only)."""
        return float(self.x), float(self.y)


@dataclass(frozen=True)
class PythagoreanDirection:
    """Exactly unit rational direction ``i^quarter · ((1 - t²), 2t) / (1 + t²)``."""

    t: Fraction
    quarter: int
    _vector: _Memo[tuple[Fraction, Fraction]] = field(default_factory=_Memo, init=False, repr=False, compare=False)

    @classmethod
    def build(cls, t: Fraction, quarter: int) -> Self:
        """Validate the grid, the half-angle cap and the quarter index.

        Args:
            t: Half-angle tangent on the ``2^-56`` grid, ``|t| <= 1/2``.
            quarter: Number of counter-clockwise quarter turns, ``0..3``.

        Returns:
            The direction.

        Raises:
            InvalidPythagoreanDirectionError: On any contract breach.
        """
        t = Fraction(t)
        if T_GRID_DENOMINATOR % t.denominator != 0:
            raise InvalidPythagoreanDirectionError(f"t = {t} is not on the 2^-{T_GRID_BITS} grid.")
        if abs(t) > MAX_ABS_HALF_ANGLE_TANGENT:
            raise InvalidPythagoreanDirectionError(f"|t| = {abs(t)} exceeds the half-angle cap {MAX_ABS_HALF_ANGLE_TANGENT}.")
        if isinstance(quarter, bool) or quarter not in range(QUARTER_TURNS):
            raise InvalidPythagoreanDirectionError(f"quarter = {quarter!r} is not one of 0..{QUARTER_TURNS - 1}.")
        return cls(t, quarter)

    @property
    def vector(self) -> tuple[Fraction, Fraction]:
        """The exact unit vector (``x² + y² == 1``)."""
        return self._vector.get(self._unit_vector)

    def _unit_vector(self) -> tuple[Fraction, Fraction]:
        denominator = 1 + self.t * self.t
        x, y = (1 - self.t * self.t) / denominator, 2 * self.t / denominator
        for _ in range(self.quarter):
            x, y = -y, x
        return x, y

    @property
    def left_normal(self) -> tuple[Fraction, Fraction]:
        """The vector turned a quarter counter-clockwise (exactly unit)."""
        x, y = self.vector
        return -y, x


@dataclass(frozen=True)
class DeclaredLine:
    """Segment ``start + λ·u``, ``λ ∈ [0, length]``, derived from source primitive ``source_index``."""

    start: ExactPoint2[WorldXY]
    direction: PythagoreanDirection
    length: ExactMillimetre
    source_index: int
    _end: _Memo[ExactPoint2[WorldXY]] = field(default_factory=_Memo, init=False, repr=False, compare=False)

    @classmethod
    def build(
        cls,
        start: ExactPoint2[WorldXY],
        direction: PythagoreanDirection,
        length: Fraction,
        source_index: int,
    ) -> Self:
        """Validate that the length is a positive binary64.

        Args:
            start: Exact start (the predecessor's exact end).
            direction: Exact unit direction.
            length: Segment length, a binary64 value.
            source_index: Index of the stored primitive this line re-emits.

        Returns:
            The declared line.

        Raises:
            InvalidDeclaredLineError: If the length is not a positive binary64.
        """
        length = Fraction(length)
        if length <= 0 or Fraction(float(length)) != length:
            raise InvalidDeclaredLineError(f"Line {source_index}: length {length} is not a positive binary64.")
        return cls(start, direction, ExactMillimetre(length), source_index)

    @property
    def end(self) -> ExactPoint2[WorldXY]:
        """Exact end ``start + length·u``."""
        return self._end.get(lambda: self.start.offset(self.direction.vector, self.length))

    @property
    def start_tangent(self) -> tuple[Fraction, Fraction]:
        """Exact unit tangent at the start."""
        return self.direction.vector

    @property
    def end_tangent(self) -> tuple[Fraction, Fraction]:
        """Exact unit tangent at the end."""
        return self.direction.vector


@dataclass(frozen=True)
class DeclaredArc:
    """Arc from ``start`` to ``end`` leaving ``start`` along ``start_tangent`` in sense ``sense``."""

    start: ExactPoint2[WorldXY]
    end: ExactPoint2[WorldXY]
    start_tangent_direction: PythagoreanDirection
    sense: Sense
    source_index: int
    _radius: _Memo[ExactMillimetre] = field(default_factory=_Memo, init=False, repr=False, compare=False)
    _centre: _Memo[ExactPoint2[WorldXY]] = field(default_factory=_Memo, init=False, repr=False, compare=False)

    @classmethod
    def build(
        cls,
        start: ExactPoint2[WorldXY],
        end: ExactPoint2[WorldXY],
        start_tangent: PythagoreanDirection,
        sense: Sense,
        source_index: int,
    ) -> Self:
        """Validate that the chord advances towards the centre.

        Args:
            start: Exact start junction.
            end: Exact end junction.
            start_tangent: Exact unit tangent at ``start``.
            sense: ``1`` counter-clockwise, ``-1`` clockwise.
            source_index: Index of the stored primitive this arc re-emits.

        Returns:
            The declared arc.

        Raises:
            DegenerateArcChordError: If ``N·(E - S) <= 0``.
        """
        if sense not in (1, -1):
            raise DegenerateArcChordError(f"Arc {source_index}: sense {sense!r} is not 1 or -1.")
        arc = cls(start, end, start_tangent, sense, source_index)
        if _dot(arc.inward_normal, end.minus(start)) <= 0:
            raise DegenerateArcChordError(f"Arc {source_index}: the chord does not advance towards the centre.")
        return arc

    @property
    def inward_normal(self) -> tuple[Fraction, Fraction]:
        """Exact unit normal at ``start`` pointing to the centre (``σ·T^⊥``)."""
        x, y = self.start_tangent_direction.left_normal
        return self.sense * x, self.sense * y

    @property
    def radius(self) -> ExactMillimetre:
        """Exact rational radius ``|E - S|² / (2 N·(E - S))``."""
        return self._radius.get(self._chord_radius)

    def _chord_radius(self) -> ExactMillimetre:
        chord = self.end.minus(self.start)
        return ExactMillimetre(_dot(chord, chord) / (2 * _dot(self.inward_normal, chord)))

    @property
    def centre(self) -> ExactPoint2[WorldXY]:
        """Exact rational centre ``S + radius·N``."""
        return self._centre.get(lambda: self.start.offset(self.inward_normal, self.radius))

    @property
    def radius_squared(self) -> SquaredMillimetre:
        """Exact squared radius."""
        return SquaredMillimetre(self.radius * self.radius)

    @property
    def start_tangent(self) -> tuple[Fraction, Fraction]:
        """Exact unit tangent at the start."""
        return self.start_tangent_direction.vector

    @property
    def end_tangent(self) -> tuple[Fraction, Fraction]:
        """Exact unit tangent at the end, ``σ·(E - c)^⊥ / ρ`` (rational)."""
        x, y = self.end.minus(self.centre)
        return -self.sense * y / self.radius, self.sense * x / self.radius


DeclaredPrimitive: TypeAlias = "DeclaredLine | DeclaredArc"


@dataclass(frozen=True)
class JunctionCensus:
    """Exact classification of every junction by the sign of ``cross(T_in, T_out)``."""

    exact_g1: int
    left_turns: int
    right_turns: int
    reversals: int


@dataclass(frozen=True)
class DeclaredBoundaryIdentity:
    """Content address of a declared boundary: contract versions plus input and output hashes."""

    contract_version: int
    t_grid_bits: int
    case_name: str
    primitive_count: int
    source_sha256: str
    declared_sha256: str


@dataclass(frozen=True)
class DeclaredBoundary:
    """Exactly closed cycle of declared primitives, in source order."""

    case_name: str
    primitives: tuple[DeclaredPrimitive, ...]
    source_sha256: str
    _identity: _Memo[DeclaredBoundaryIdentity] = field(default_factory=_Memo, init=False, repr=False, compare=False)

    @classmethod
    def build(cls, case_name: str, primitives: tuple[DeclaredPrimitive, ...], source_sha256: str) -> Self:
        """Validate provenance and exact closure.

        Args:
            case_name: Name of the reference case.
            primitives: Declared primitives; primitive ``i`` re-emits source primitive ``i``.
            source_sha256: Hash of the stored primitives the boundary was derived from.

        Returns:
            The declared boundary.

        Raises:
            InvalidDeclaredProvenanceError: If the source indices are not ``0..n-1`` in order.
            OpenDeclaredCycleError: If any primitive's exact end differs from its successor's start.
        """
        if not primitives or tuple(primitive.source_index for primitive in primitives) != tuple(range(len(primitives))):
            raise InvalidDeclaredProvenanceError(f"{case_name}: declared primitives must re-emit source indices 0..n-1 in order.")
        for current, following in zip(primitives, (*primitives[1:], primitives[0])):
            if current.end != following.start:
                raise OpenDeclaredCycleError(
                    f"{case_name}: primitive {current.source_index} ends at {current.end.as_floats()}, primitive {following.source_index} starts at {following.start.as_floats()}."
                )
        return cls(case_name, primitives, source_sha256)

    def junction_census(self) -> JunctionCensus:
        """Classify each junction exactly (G1 is not imposed, so kinks are counted, not refused)."""
        exact_g1 = left = right = reversals = 0
        for incoming, outgoing in zip((self.primitives[-1], *self.primitives[:-1]), self.primitives):
            turn = _cross(incoming.end_tangent, outgoing.start_tangent)
            if turn > 0:
                left += 1
            elif turn < 0:
                right += 1
            elif _dot(incoming.end_tangent, outgoing.start_tangent) > 0:
                exact_g1 += 1
            else:
                reversals += 1
        return JunctionCensus(exact_g1, left, right, reversals)

    def canonical_bytes(self) -> bytes:
        """Deterministic serialisation of every exact value (the hashed content)."""
        records: list[dict[str, object]] = []
        for primitive in self.primitives:
            record: dict[str, object] = {
                "source_index": primitive.source_index,
                "start": [_fraction_text(primitive.start.x), _fraction_text(primitive.start.y)],
            }
            if isinstance(primitive, DeclaredLine):
                record |= {
                    "kind": "line",
                    "t": _fraction_text(primitive.direction.t),
                    "quarter": primitive.direction.quarter,
                    "length": _fraction_text(primitive.length),
                }
            else:
                record |= {
                    "kind": "arc",
                    "end": [_fraction_text(primitive.end.x), _fraction_text(primitive.end.y)],
                    "t": _fraction_text(primitive.start_tangent_direction.t),
                    "quarter": primitive.start_tangent_direction.quarter,
                    "sense": primitive.sense,
                }
            records.append(record)
        document = {
            "contract_version": DECLARED_BOUNDARY_CONTRACT_VERSION,
            "t_grid_bits": T_GRID_BITS,
            "case_name": self.case_name,
            "primitives": records,
        }
        return json.dumps(document, sort_keys=True, separators=(",", ":")).encode()

    @property
    def identity(self) -> DeclaredBoundaryIdentity:
        """Content address of this boundary."""
        return self._identity.get(self._content_identity)

    def _content_identity(self) -> DeclaredBoundaryIdentity:
        return DeclaredBoundaryIdentity(
            contract_version=DECLARED_BOUNDARY_CONTRACT_VERSION,
            t_grid_bits=T_GRID_BITS,
            case_name=self.case_name,
            primitive_count=len(self.primitives),
            source_sha256=self.source_sha256,
            declared_sha256=hashlib.sha256(self.canonical_bytes()).hexdigest(),
        )


def _fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def _dot(left: tuple[Fraction, Fraction], right: tuple[Fraction, Fraction]) -> Fraction:
    return left[0] * right[0] + left[1] * right[1]


def _cross(left: tuple[Fraction, Fraction], right: tuple[Fraction, Fraction]) -> Fraction:
    return left[0] * right[1] - left[1] * right[0]
