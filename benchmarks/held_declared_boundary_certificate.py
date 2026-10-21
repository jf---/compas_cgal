"""Compose a bound (B + δ) on a declared Held boundary's distance to the publisher PDF evidence.

This is a composed bound, not a stand-alone certificate: it is conditional on
one open item, named in the type as ``OpenLemma.BIARC_POINT_SET_BOUND`` (see
below). The stored primitives carry the existing reconstruction certificate: the
loader (:func:`benchmarks.held_reference_cases.validate_case_payload`) re-runs
:func:`benchmarks.held_reference_geometry.reconstruct_source_path` on the
publisher source curves and accepts the case only if it reproduces the
recorded primitives and their deviation bound ``B``. That bound is a two-sided
point-set bound between each source curve and its stored primitive (radial
annulus plus monotone polar angle for circle spans, parameter correspondence
for biarcs, control hull for lines).

Open item (``OpenLemma.BIARC_POINT_SET_BOUND``): for biarc spans, that B is a
two-sided point-set bound is read from the construction (a parameter
correspondence evaluated on float inputs), not proved. Circle spans are
proved (exact hull/annulus enclosures rounded outward). Until the lemma is
proved, or the per-primitive certificate is ported to rationals, B + δ is a
bound composed under that lemma.

B enters exact-land once, as the next binary64 above the recorded double
(declared API semantic): the loader guarantees the record is no lower than
its recomputation, and reading it one step up also covers a record that
rounded a true bound to nearest from below.

The declared boundary's bound is then composed by the triangle inequality:
``dist(source, declared) <= B + δ`` with ``δ`` a rigorous rational upper bound
on the two-sided Hausdorff distance between every stored primitive and its
declared re-emission. The decision ``B + δ <= limit`` is an exact rational
comparison; square roots enter only as outward rational enclosures.

Stored arc convention (the same as the gouge audit's): the circle about
``c0`` of radius ``ρ0 = |S0 - c0|`` over the directions from ``dir(S0 - c0)``
to ``dir(E0 - c0)`` in the sense of the stored sweep.

Arc move bound. With ``c, ρ`` the declared centre and radius,
``Δ = c - c0`` and ``pow0(P) = |P - c0|² - ρ0²``:

- every declared point ``P = c + ρz`` has ``pow0(P) = K + 2ρ z·Δ`` with
  ``K = |Δ|² + ρ² - ρ0²``, linear in the unit vector ``z``, so its extreme
  over the declared arc is at an end (``|S - c0|² - ρ0²``,
  ``|E - c0|² - ρ0²``, rational) or at ``z = ±Δ/|Δ|`` when that lies in the
  arc's wedge (decided by exact orientation);
- ``G = max|pow0| / ρ0`` bounds the radial gap ``||P - c0| - ρ0||``;
- with ``c0`` strictly inside the declared circle, rays from ``c0`` pair the
  two arcs one-to-one over their common directions, at distance ``<= G``;
- the arcs end on the same exact point, so their direction ranges differ only
  by a start sliver; with ``s = |S - S0|`` its points are within
  ``max(2G + s, G + 2s)`` of the other arc.

Hence ``δ_arc <= G`` when ``S = S0`` and ``δ_arc <= 2G + 2s`` otherwise. A line
moves by at most ``max(|S - S0|, |E - E0|)`` (convex correspondence).
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum
from fractions import Fraction
from typing import Literal
from typing import Self

from benchmarks.held_declared_boundary import DeclaredArc
from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_declared_boundary import DeclaredBoundaryError
from benchmarks.held_declared_boundary import DeclaredBoundaryIdentity
from benchmarks.held_declared_boundary import DeclaredLine
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_boundary_emission import stored_primitives_sha256
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_reference_geometry import ReferenceArc
from benchmarks.held_reference_geometry import ReferenceLine
from benchmarks.held_reference_geometry import ReferencePrimitive
from compas_cgal.units import ExactMillimetre
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY

# Square-root enclosures are dyadic with 2^-192 absolute width: the quantities
# bounded are >= 1e-30 mm wide in practice, so the enclosure slack is < 1e-27
# of any reported bound.
SQRT_ENCLOSURE_BITS = 192
# The loader requires the normalised stroke width to be four reconstruction
# limits (benchmarks/held_reference_cases.py, "one quarter of the stroke").
STROKE_WIDTHS_PER_RECONSTRUCTION_LIMIT = 4

_Vector = tuple[Fraction, Fraction]
_IntegerVector = tuple[int, int]


class UncertifiableArcMoveError(DeclaredBoundaryError):
    """A stored/declared arc pair breaks a precondition of the exact move bound."""


class DeclaredProvenanceMismatchError(DeclaredBoundaryError):
    """A declared boundary was not re-emitted from the case's stored primitives."""


class DeclaredBoundaryRecertificationError(DeclaredBoundaryError):
    """The composed bound ``B + δ`` exceeds the case's reconstruction limit."""


class OpenLemma(Enum):
    """An unproved statement a composed bound is conditional on."""

    BIARC_POINT_SET_BOUND = "biarc-point-set-bound"
    """The recorded B is a two-sided point-set bound for biarc spans (read from the construction, not proved)."""


@dataclass(frozen=True)
class PrimitiveMove:
    """Rigorous upper bound on the Hausdorff move of one primitive."""

    source_index: int
    kind: Literal["line", "arc"]
    bound: ExactMillimetre
    radial_gap: ExactMillimetre
    start_gap: ExactMillimetre


@dataclass(frozen=True)
class DeclaredBoundaryCertificate:
    """The declared boundary's composed (B + δ) deviation bound from the publisher source curves.

    δ is proved here (exact rationals, outward square-root enclosures); B is the
    loader's reconstruction bound, injected one binary64 step up. The whole
    bound holds conditional on ``conditional_on``; it is not a certificate for
    biarc-derived spans until ``OpenLemma.BIARC_POINT_SET_BOUND`` is proved.
    """

    boundary: DeclaredBoundaryIdentity
    stored_bound: ExactMillimetre
    maximum_move: ExactMillimetre
    declared_bound: ExactMillimetre
    reconstruction_limit: ExactMillimetre
    moves: tuple[PrimitiveMove, ...]
    conditional_on: tuple[OpenLemma, ...]

    @classmethod
    def build(cls, case: HeldReferenceCase, declared: DeclaredBoundary) -> Self:
        """Compose the loader-verified reconstruction bound with the exact move bound.

        The recorded B is read as ``nextafter(B, +inf)`` (API semantic, see the
        module docstring) and the result is conditional on
        ``OpenLemma.BIARC_POINT_SET_BOUND``.

        Args:
            case: A loaded case (its reconstruction certificate was re-run on load).
            declared: The declared boundary re-emitted from ``case.boundary``.

        Returns:
            The certificate.

        Raises:
            DeclaredProvenanceMismatchError: If ``declared`` was not re-emitted from this case.
            DeclaredBoundaryRecertificationError: If ``B + δ`` exceeds the reconstruction limit.
        """
        stored = case.boundary.primitives
        if declared.case_name != case.name or declared.source_sha256 != stored_primitives_sha256(stored):
            raise DeclaredProvenanceMismatchError(f"Declared boundary {declared.case_name!r} was not re-emitted from case {case.name!r}.")
        moves = measure_moves(stored, declared)
        maximum_move = max(move.bound for move in moves)
        stored_bound = Fraction(math.nextafter(float(case.reconstruction.deviation_upper_bound), math.inf))
        declared_bound = stored_bound + maximum_move
        limit = Fraction(float(case.boundary.boundary_stroke_width)) / STROKE_WIDTHS_PER_RECONSTRUCTION_LIMIT
        if declared_bound > limit:
            raise DeclaredBoundaryRecertificationError(f"{case.name}: B + δ = {float(declared_bound)!r} mm exceeds the reconstruction limit {float(limit)!r} mm.")
        return cls(
            declared.identity,
            ExactMillimetre(stored_bound),
            ExactMillimetre(maximum_move),
            ExactMillimetre(declared_bound),
            ExactMillimetre(limit),
            moves,
            (OpenLemma.BIARC_POINT_SET_BOUND,),
        )

    @property
    def identity_sha256(self) -> str:
        """Content address of the certified boundary: declared hash plus exact bounds."""
        document = {
            "declared_sha256": self.boundary.declared_sha256,
            "source_sha256": self.boundary.source_sha256,
            "contract_version": self.boundary.contract_version,
            "stored_bound": _text(self.stored_bound),
            "maximum_move": _text(self.maximum_move),
            "declared_bound": _text(self.declared_bound),
            "reconstruction_limit": _text(self.reconstruction_limit),
            "conditional_on": [lemma.value for lemma in self.conditional_on],
        }
        return hashlib.sha256(json.dumps(document, sort_keys=True).encode()).hexdigest()


def measure_moves(stored: Sequence[ReferencePrimitive], declared: DeclaredBoundary) -> tuple[PrimitiveMove, ...]:
    """Bound the move of every declared primitive from its stored source.

    Args:
        stored: Stored primitives, indexed as the declared provenance.
        declared: Declared boundary.

    Returns:
        One move bound per primitive.

    Raises:
        DeclaredProvenanceMismatchError: If kinds or counts disagree.
    """
    if len(stored) != len(declared.primitives):
        raise DeclaredProvenanceMismatchError("Stored and declared primitive counts differ.")
    moves: list[PrimitiveMove] = []
    for primitive in declared.primitives:
        source = stored[primitive.source_index]
        if isinstance(primitive, DeclaredLine) and isinstance(source, ReferenceLine):
            moves.append(line_move_bound(source, primitive))
        elif isinstance(primitive, DeclaredArc) and isinstance(source, ReferenceArc):
            moves.append(arc_move_bound(source, primitive))
        else:
            raise DeclaredProvenanceMismatchError(f"Primitive {primitive.source_index}: stored and declared kinds differ.")
    return tuple(moves)


def line_move_bound(stored: ReferenceLine, declared: DeclaredLine) -> PrimitiveMove:
    """Hausdorff bound between a stored segment and its declared re-emission.

    Args:
        stored: Stored segment.
        declared: Declared segment.

    Returns:
        ``max(|S - S0|, |E - E0|)`` enclosed from above.
    """
    start_gap = sqrt_upper(_squared(declared.start.minus(_exact(stored.start))))
    end_gap = sqrt_upper(_squared(declared.end.minus(_exact(stored.end))))
    return PrimitiveMove(declared.source_index, "line", ExactMillimetre(max(start_gap, end_gap)), ExactMillimetre(Fraction(0)), ExactMillimetre(start_gap))


def arc_move_bound(stored: ReferenceArc, declared: DeclaredArc) -> PrimitiveMove:
    """Hausdorff bound between a stored arc and its declared re-emission (module docstring).

    Args:
        stored: Stored arc.
        declared: Declared arc.

    Returns:
        The move bound with its radial-gap and start-sliver components.

    Raises:
        UncertifiableArcMoveError: If a precondition of the bound fails.
    """
    index = declared.source_index
    stored_sense = 1 if float(stored.sweep) > 0.0 else -1
    if stored_sense != declared.sense:
        raise UncertifiableArcMoveError(f"Arc {index}: stored and declared senses differ.")
    if declared.end != _exact(stored.end):
        raise UncertifiableArcMoveError(f"Arc {index}: the declared arc does not end on the stored end.")
    # One exact integer frame per arc: every coordinate (and the radius) times
    # one positive integer ``scale``, refined by 2^SQRT_ENCLOSURE_BITS so that
    # integer square-root enclosures are at least that fine in millimetres.
    exact_values = (
        declared.centre.x,
        declared.centre.y,
        declared.start.x,
        declared.start.y,
        declared.end.x,
        declared.end.y,
        Fraction(float(stored.centre.x)),
        Fraction(float(stored.centre.y)),
        Fraction(float(stored.start.x)),
        Fraction(float(stored.start.y)),
        declared.radius,
    )
    scale = math.lcm(*(value.denominator for value in exact_values)) << SQRT_ENCLOSURE_BITS
    cx, cy, sx, sy, ex, ey, c0x, c0y, s0x, s0y, radius = (int(value * scale) for value in exact_values)
    stored_radius_squared = (s0x - c0x) ** 2 + (s0y - c0y) ** 2
    shift = (cx - c0x, cy - c0y)
    shift_squared = shift[0] ** 2 + shift[1] ** 2
    if shift_squared >= radius * radius:
        raise UncertifiableArcMoveError(f"Arc {index}: the stored centre is not strictly inside the declared circle.")
    start_radial = (sx - cx, sy - cy)
    end_radial = (ex - cx, ey - cy)
    declared_major = _declared_is_major(start_radial, end_radial, declared.sense)
    # Float decision on the stored binary64 sweep against fl(pi): a precondition
    # guard on recorded data. Only a sweep within an ulp of pi could be misjudged;
    # the four cases' stored sweeps all sit >= 2.35 rad from pi (measured 2026-10-21).
    if declared_major is not None and declared_major != (abs(float(stored.sweep)) > math.pi):
        raise UncertifiableArcMoveError(f"Arc {index}: stored and declared arcs disagree on major/minor.")

    power_extremes = [
        abs((sx - c0x) ** 2 + (sy - c0y) ** 2 - stored_radius_squared),
        abs((ex - c0x) ** 2 + (ey - c0y) ** 2 - stored_radius_squared),
    ]
    if shift_squared > 0:
        constant = shift_squared + radius * radius - stored_radius_squared
        shift_lower = math.isqrt(shift_squared)
        shift_upper = shift_lower if shift_lower * shift_lower == shift_squared else shift_lower + 1
        for orientation in (1, -1):
            if _in_wedge(start_radial, end_radial, (orientation * shift[0], orientation * shift[1]), declared.sense):
                twice_radius = 2 * orientation * radius
                power_extremes.append(max(abs(constant + twice_radius * shift_lower), abs(constant + twice_radius * shift_upper)))
    # ||P - c0| - ρ0| <= |pow0(P)| / ρ0; in the frame pow0 scales by scale² and ρ0 by scale.
    radial_gap = Fraction(max(power_extremes), scale * math.isqrt(stored_radius_squared))

    start_offset = (sx - s0x, sy - s0y)
    if start_offset == (0, 0):
        return PrimitiveMove(index, "arc", ExactMillimetre(radial_gap), ExactMillimetre(radial_gap), ExactMillimetre(Fraction(0)))
    if (sx - c0x) * (s0x - c0x) + (sy - c0y) * (s0y - c0y) <= 0:
        raise UncertifiableArcMoveError(f"Arc {index}: the start sliver spans a quarter turn or more.")
    offset_squared = start_offset[0] ** 2 + start_offset[1] ** 2
    offset_root = math.isqrt(offset_squared)
    start_gap = Fraction(offset_root if offset_root * offset_root == offset_squared else offset_root + 1, scale)
    return PrimitiveMove(index, "arc", ExactMillimetre(2 * radial_gap + 2 * start_gap), ExactMillimetre(radial_gap), ExactMillimetre(start_gap))


def sqrt_upper(value: Fraction) -> Fraction:
    """Dyadic upper enclosure of ``√value`` (width ``2^-SQRT_ENCLOSURE_BITS``).

    Args:
        value: Non-negative rational.

    Returns:
        A rational ``>= √value``.
    """
    scaled = value * 4**SQRT_ENCLOSURE_BITS
    ceiling = -((-scaled.numerator) // scaled.denominator)
    root = math.isqrt(ceiling)
    if root * root < ceiling:
        root += 1
    return Fraction(root, 2**SQRT_ENCLOSURE_BITS)


def sqrt_lower(value: Fraction) -> Fraction:
    """Dyadic lower enclosure of ``√value`` (width ``2^-SQRT_ENCLOSURE_BITS``).

    Args:
        value: Non-negative rational.

    Returns:
        A rational ``<= √value``.
    """
    scaled = value * 4**SQRT_ENCLOSURE_BITS
    return Fraction(math.isqrt(scaled.numerator // scaled.denominator), 2**SQRT_ENCLOSURE_BITS)


def _in_wedge(start: _IntegerVector, end: _IntegerVector, probe: _IntegerVector, sense: int) -> bool:
    """Whether ``probe`` lies in the closed ``sense`` arc of directions from ``start`` to ``end``.

    Exactly-antipodal ends are classified conservatively (a superset), which
    can only enlarge the bound.
    """
    if sense < 0:
        start, end = end, start
    turn = start[0] * end[1] - start[1] * end[0]
    probe_after_start = start[0] * probe[1] - start[1] * probe[0]
    probe_before_end = probe[0] * end[1] - probe[1] * end[0]
    if turn > 0:
        return probe_after_start >= 0 and probe_before_end >= 0
    if turn < 0:
        return probe_after_start >= 0 or probe_before_end >= 0
    return True


def _declared_is_major(start: _IntegerVector, end: _IntegerVector, sense: int) -> bool | None:
    """Exact major/minor class of the declared arc; ``None`` for an exact half turn."""
    turn = sense * (start[0] * end[1] - start[1] * end[0])
    if turn == 0 and start[0] * end[0] + start[1] * end[1] < 0:
        return None
    return turn < 0


def _exact(point: Point2[WorldXY]) -> ExactPoint2[WorldXY]:
    return ExactPoint2[WorldXY].from_binary64(float(point.x), float(point.y))


def _squared(vector: _Vector) -> Fraction:
    return vector[0] * vector[0] + vector[1] * vector[1]


def _text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"
