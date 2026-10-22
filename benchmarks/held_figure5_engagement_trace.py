"""What each refinement decision read beyond its own circle pair: the tier-1b filter's piece bounds and the trace.

``EngagementTrace`` records, per decision, the inputs that make it depend on
emitted circles other than its predecessor and target; the coverage repair's
splice stops a continuation early from it (``continue_with_early_stop``).
``tier_1b_kept`` is the gate's superset filter and the trace's influence test
alike: one implementation, so replay identity never rests on two copies.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction

import numpy as np
import numpy.typing as npt

from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from compas_cgal.motion import EngagementCap
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

# binary64 round-to-nearest unit roundoff u: |fl(x op y) - x op y| <= u |x op y|.
_UNIT_ROUNDOFF = 2.0**-53

FloatArray = npt.NDArray[np.float64]
BoolArray = npt.NDArray[np.bool_]


def rational_point(point: Point2[WorldXY]) -> tuple[Fraction, Fraction]:
    return Fraction(float(point.x)), Fraction(float(point.y))


def _down(value: FloatArray) -> FloatArray:
    """The next binary64 below: a strict lower bound on the exact result of the op that produced ``value``."""
    return np.nextafter(value, -np.inf)


def _up(value: FloatArray) -> FloatArray:
    """The next binary64 above: a strict upper bound on the exact result of the op that produced ``value``."""
    return np.nextafter(value, np.inf)


def tier_1b_reach(circle: Figure5CounterclockwiseCircle, tool_radius: ToolRadius) -> float:
    """The binary64 reach ``fl(rho + 2 r)`` the tier-1b filter tests a circle with.

    The single source for the gate (``_ConnectorGate._nearby``) and the
    trace's influence test (``EngagementTrace.last_influenced``): the early
    stop is exact because both replay the same operations on the same doubles.
    """
    return float(circle.radius.value) + 2 * float(tool_radius.value)


def floor_reach(tool_radius: ToolRadius) -> Fraction:
    """The exact ``2 r`` of the convex-corner floor hypothesis, for the gate and the trace alike."""
    return 2 * Fraction(float(tool_radius.value))


@dataclass(frozen=True, eq=False)
class PieceBounds:
    """Outward binary64 bounds of one connector's sample pieces, the per-piece half of ``_tier_1b_neighbours``.

    Per piece ``[A, B]``: ``mid_lo <= M <= mid_hi`` per axis (``M`` the exact
    midpoint), ``half >= H`` (the exact half-length) and ``extent = 8 P``
    (``P = max |A, B coordinates|``, exact: a power-of-two scaling of a double).
    The arrays are read-only; equality and hash are bit identity of all four
    (the replay identity the early stop relies on), so footprints and traces
    compare and hash by value.
    """

    mid_lo: FloatArray
    mid_hi: FloatArray
    half: FloatArray
    extent: FloatArray

    def __post_init__(self) -> None:
        for array in (self.mid_lo, self.mid_hi, self.half, self.extent):
            array.setflags(write=False)

    def _bits(self) -> tuple[bytes, ...]:
        return tuple(array.tobytes() for array in (self.mid_lo, self.mid_hi, self.half, self.extent))

    def __eq__(self, other: object) -> bool:
        return isinstance(other, PieceBounds) and self._bits() == other._bits()

    def __hash__(self) -> int:
        return hash(self._bits())

    @classmethod
    def of(cls, samples: Sequence[tuple[float, float]]) -> PieceBounds:
        points = np.array(samples, dtype=np.float64).reshape(-1, 2)
        a, b = points[:-1], points[1:]
        total = a + b
        width = _up(np.abs(b - a))
        half = _up(_up(np.sqrt(_up(_up(width[:, 0] * width[:, 0]) + _up(width[:, 1] * width[:, 1])))) / 2.0)
        return cls(_down(_down(total) / 2.0), _up(_up(total) / 2.0), half, 8.0 * np.maximum(np.abs(a), np.abs(b)).max(axis=1, initial=0.0))

    @classmethod
    def concatenate(cls, parts: Sequence[PieceBounds]) -> PieceBounds:
        if not parts:
            return cls.of(())
        return cls(*(np.concatenate([getattr(part, name) for part in parts]) for name in ("mid_lo", "mid_hi", "half", "extent")))

    def __len__(self) -> int:
        return len(self.half)


def tier_1b_kept(centres: FloatArray, reach: FloatArray, pieces: PieceBounds) -> BoolArray:
    """``(circles, pieces)``: circle ``i`` is not provably beyond ``_CandidateIndex.near`` on piece ``j``.

    The test of ``_tier_1b_neighbours`` (which states and proves it), every
    binary64 operation elementwise and one ulp outward; broadcasting does not
    change a single rounding.
    """
    cx, cy = centres[:, 0:1], centres[:, 1:2]
    spread = _up(2.0 * np.maximum(np.abs(cx), np.abs(cy)))
    lx = np.maximum(0.0, np.maximum(_down(cx - pieces.mid_hi[:, 0]), _down(pieces.mid_lo[:, 0] - cx)))
    ly = np.maximum(0.0, np.maximum(_down(cy - pieces.mid_hi[:, 1]), _down(pieces.mid_lo[:, 1] - cy)))
    lower = _down(_down(lx * lx) + _down(ly * ly))
    span = _up(reach[:, None] + pieces.half)
    error = _up(_UNIT_ROUNDOFF * _up(_up(_up(6.0 * span) + pieces.extent) + spread))
    bound = _up(span + error)
    return ~(lower > _up(bound * bound))


def _tier_1b_neighbours(centres: FloatArray, reach: FloatArray, samples: Sequence[tuple[float, float]]) -> BoolArray:
    """Circles not provably beyond ``_CandidateIndex.near`` for every dyadic sub-piece of every piece.

    ``reach`` is the candidate filter's own binary64 ``rho + 2 r``. For a piece
    ``[A, B]`` (exact doubles) with midpoint ``M``, half-length ``H``,
    ``P = max |A, B coordinates|`` and, per circle, ``C = max |c coordinates|``,
    the circle is kept unless ``|c - M| > R + H + E`` provably, with
    ``E = u (6 (R + H) + 8 P + 2 C)``.

    Why that is a superset. ``_CandidateIndex.near`` asks about a sub-piece
    ``[p, q]`` of the dyadic grid, which lies exactly on ``[A, B]`` (the grid
    keeps ``MAX_BISECTION_DEPTH`` spare bits, so ``_at`` never rounds); its
    midpoint ``m`` and half-length ``h`` satisfy ``|m - M| <= H - h``. It
    keeps a circle only if ``within`` (an outward bracket proving
    ``|c - m| <= R_low``, ``R_low`` one ulp below ``R``) or ``_within_exactly``
    (``|c - m| <= R* + h`` decided by squaring on the grid integers), where
    ``R* = rho + 2 r`` is the exact sum of the two doubles. No libm function
    is involved, so either way a kept circle has ``|c - m| <= R* + h`` exactly,
    and the triangle inequality gives ``|c - M| <= R* + H``. ``R`` here is the
    binary64 ``fl(rho + 2 r)`` (the same ``float(rho) + 2 float(r)``), one
    correctly rounded addition (the doubling is exact), so
    ``R* <= R / (1 - u) < R + u R (1 + 2 u)``. Hence a kept circle has
    ``|c - M| < R + H + u R (1 + 2 u) <= R + H + E``, since
    ``E >= 6 u R``. ``E``'s other terms (``8 P``, ``2 C``, the ``6 u H``) were
    sized for an earlier binary64 ``hypot`` test of ``near`` and are no longer
    needed; they are kept (``tests/benchmarks/test_held_tier_1b_error_term_pins.py``
    pins them) because extra width only admits extra circles, which the
    exact gate then decides on, and never drops one.

    Why the test is outward. Round-to-nearest errs by at most half an ulp, so
    nudging every op's result one ulp outward (``_down`` / ``_up``) brackets the
    exact value: ``lower <= |c - M|**2`` (``lx, ly`` lower-bound the distance
    from ``c`` to the box ``[mid_lo, mid_hi]`` that contains ``M``),
    ``half >= H`` and ``upper >= (R + H + E)**2``, and a circle is dropped only
    when ``lower > upper``. Dropping a circle could only lose a certificate,
    never create one; keeping a superset makes the gate decide exactly as over
    every emitted circle.
    """
    return tier_1b_kept(centres, reach, PieceBounds.of(samples)).any(axis=1)


@dataclass(frozen=True)
class EngagementFootprint:
    """One refinement decision that read emitted circles beyond its own pair, at emission ``position``.

    ``position`` is ``len(emitted)`` when it was decided. ``pieces`` is set
    for a tier-1b connector decision (the gate's ``_tier_1b_neighbours``
    bounds): a circle ``tier_1b_kept`` by none of them is not on its grid, so it
    cannot change the decision. ``vertex`` is set for a convex-corner floor
    hypothesis (``_ConnectorGate._floor_hypothesis``), which skips every
    circle with ``|c - V| > rho + 2 r`` exactly. Every other decision reads
    only its predecessor, its target and the boundary.
    """

    position: int
    pieces: PieceBounds | None = None
    vertex: int | None = None

    def shifted(self, delta: int) -> EngagementFootprint:
        return EngagementFootprint(self.position + delta, self.pieces, self.vertex)


@dataclass(frozen=True)
class EngagementRefinementKey:
    """Everything besides the sources and the circle budget that a refinement's decisions read."""

    boundary: tuple[Point2[WorldXY], ...]
    tool_radius: ToolRadius
    cap: EngagementCap
    max_depth: int
    corner_approaches: bool
    close_seam: bool


@dataclass(frozen=True)
class EngagementTrace:
    """The footprints of every decision of one refinement, in emission order, and the settings they were made under."""

    key: EngagementRefinementKey
    footprints: tuple[EngagementFootprint, ...]

    def last_influenced(self, circles: Sequence[Figure5CounterclockwiseCircle], since: int) -> int:
        """The last footprint position ``>= since`` whose decision ``circles`` could change, or -1.

        Exact: a circle changes a tier-1b decision only if the gate's own
        superset filter keeps it for one of the decision's pieces (``tier_1b_kept``,
        the same binary64 operations on the same doubles), and a floor
        hypothesis only if ``|c - V|^2 <= (rho + 2 r)^2`` in ``Fraction``.
        """
        later = [footprint for footprint in self.footprints if footprint.position >= since]
        if not circles or not later:
            return -1
        centres = np.array([(float(c.center.x), float(c.center.y)) for c in circles], dtype=np.float64).reshape(-1, 2)
        reach = np.array([tier_1b_reach(c, self.key.tool_radius) for c in circles], dtype=np.float64)
        tier_1b = [(footprint.position, footprint.pieces) for footprint in later if footprint.pieces is not None]
        owners = np.repeat(np.array([position for position, _ in tier_1b], dtype=np.int64), [len(pieces) for _, pieces in tier_1b])
        hit = tier_1b_kept(centres, reach, PieceBounds.concatenate([pieces for _, pieces in tier_1b])).any(axis=0)
        last = int(owners[hit].max()) if hit.any() else -1
        reach_2r = floor_reach(self.key.tool_radius)
        for footprint in later:
            if footprint.vertex is None or footprint.position <= last:
                continue
            vx, vy = rational_point(self.key.boundary[footprint.vertex])
            for circle in circles:
                (x, y), rho = rational_point(circle.center), Fraction(float(circle.radius.value))
                if (x - vx) * (x - vx) + (y - vy) * (y - vy) <= (rho + reach_2r) * (rho + reach_2r):
                    last = footprint.position
                    break
        return last
