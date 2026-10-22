"""Exact removal of NOTHING circles: circles whose swept region removes no material.

A full machining circle of guide radius rho about c, swept by a tool of radius
r, sweeps the closed annulus A = {p : max(rho - r, 0) <= |p - c| <= rho + r}.
A circle removes NOTHING when A lies in the union of the regions swept before
it. This module proves a sufficient form of that statement exactly: A lies in
the union of the swept annuli of EARLIER CIRCLES (connectors ignored). Leaving
connectors out only shrinks the union, so every circle proven here removes
nothing on the true stock as well; circles covered only with a connector's help
are kept.

Decision (exact over the injected doubles): a local exact stock
(``_stock_2.Stock2``, Epeck circle-segment traits) is seeded with an axis square
containing A, every selected earlier annulus is subtracted exactly
(``subtract_circle_sweep``), and ``Stock2.circle_sweep_removal`` decides A
against the rest (a regularized interior intersection). ``NOTHING`` is the
proof.

Selection is not decision. Which earlier circles enter the union is chosen in
binary64 (``_overlapping_earlier``); leaving a sweep out can only keep a
circle, never remove one.

The keep screen is exact. A circle skips the exact union test only on an exact
witness: a rational point of its closed annulus outside every selected closed
annulus (``_DyadicImage``, integer arithmetic over the injected doubles).
Such a point proves material is left: the complement of the closed union is
open, so a disc about the witness misses the union, and that disc meets the
annulus's interior in positive area (a closed annulus with r > 0 is the closure
of its interior); that is exactly a non-NOTHING verdict of the regularized
removal test. Binary64 only PROPOSES candidate points (the
uncovered angle gaps of the annulus's own rims for the whole path at once,
``_batch_witnessed``; then of every boundary circle met by the annulus,
``_rim_gaps``, in binary64 and at sub-ulp dyadic offsets);
a candidate that fails the exact check is discarded, so a float error can only
send a circle on to the exact union test, never keep or remove one. With no
witness the exact union test decides, subtracting the selected annuli nearest
first and stopping as soon as NOTHING is proven (the stock only shrinks, so a
NOTHING verdict on a prefix of the sweeps holds for all of them); the verdict
is the one the full subtraction gives.

Removing circle k makes circle k+1's predecessor the last kept circle; that
pair's engagement is re-decided by the exact standard-model predicate
(``exact_predecessor_engagement(...).within(cap)``) before k is removed, and a
NOTHING circle whose removal would put its successor over the cap is kept.
On a closed path the last circle's successor is circle 0 across the seam, so
removing the current last circle re-decides the seam pair (last kept, circle 0)
the same way (see :func:`remove_nothing_circles` for the exact rule).
Connectors are built after this pass from the surviving circles, so the
connector from the last kept circle to k+1 is formed by the ordinary connector
builder.
"""

from __future__ import annotations

from collections.abc import Callable
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import Generic
from typing import TypeVar

import numpy as np

from benchmarks.errors import BenchmarkError
from benchmarks.held_standard_placement import PaperCircleCandidate
from benchmarks.held_standard_placement import StandardPlacementFragmentationError
from benchmarks.held_standard_placement import exact_predecessor_engagement
from compas_cgal import _stock_2
from compas_cgal.motion import EngagementCap
from compas_cgal.units import ToolRadius

C = TypeVar("C")

# Rim-gap candidates: the widest uncovered gaps tried per circle, and the radial
# offsets tried off each gap's boundary circle, as fractions of the binary64
# clearance at the gap's midpoint (half of it, then 1/64 of it). Proposal only:
# every candidate is checked exactly. Measured on the 10 workload paths
# (build/nu/stages.py): the target's own rims yield an exact witness for every
# circle that is not proven NOTHING or a sub-ulp sliver.
RIM_GAPS_TRIED = 4
RIM_GAP_OFFSETS = (0.5, 2.0**-6)

# Sub-ulp rational offsets off a rim gap, as 2**-k of the rim radius (one ulp of
# a coordinate is about 2**-52 of it). Proposal only.
DEEP_OFFSET_EXPONENTS = (48, 54, 60, 70, 85, 100)

# Selection and own-rim gap proposals run for the whole path in row blocks of
# this many circles (memory: a block holds block x n binary64 distances).
SELECTION_BLOCK = 256

# Angle offset separating the rows of the batched gap search: a row's angles lie
# in [-1, 4 pi], so rows offset by 8 pi never interleave in one sort.
ROW_SPAN = 8.0 * np.pi

# Exact union test: subtract this many sweeps (nearest first) before the first
# NOTHING query, doubling the count between queries. A query costs about two
# subtractions (build/nu/prof.py), so doubling keeps the queries below half the
# subtraction work.
FIRST_UNION_QUERY = 8

# Half-width factor of the local stock's seed square about the annulus's outer
# radius. The square's exact containment of the annulus is checked, not assumed
# (``_seed_square``); the factor only keeps that check from ever failing (the
# corners' binary64 rounding is ulps, far below 1/16 of the radius). A tighter
# square keeps less of the earlier sweeps' arrangement: 1.0625 against 1.25 cut
# the upper_lane union tests from 0.83 s to 0.72 s CPU (build/nu/prof2.py).
SEED_SQUARE_FACTOR = 1.0625


class NothingCircleError(BenchmarkError):
    """Base error of the NOTHING-circle pass."""


class NonDyadicWitnessError(NothingCircleError):
    """A proposed witness point is not a dyadic rational (proposals are built from doubles and powers of two)."""


class SeedSquareContainmentError(NothingCircleError):
    """The local stock's binary64 seed square does not contain the annulus exactly."""


@dataclass(frozen=True)
class _Annulus:
    """A swept full circle, exact: centre (x, y), guide radius, tool radius."""

    x: Fraction
    y: Fraction
    guide: Fraction
    tool: Fraction

    @classmethod
    def of(cls, circle: PaperCircleCandidate, tool_radius: ToolRadius) -> _Annulus:
        return cls(Fraction(float(circle.center.x)), Fraction(float(circle.center.y)), Fraction(float(circle.guide_radius.value)), Fraction(float(tool_radius.value)))

    def contains(self, px: Fraction, py: Fraction) -> bool:
        """Closed membership of the exact point (px, py)."""
        squared = (px - self.x) ** 2 + (py - self.y) ** 2
        outer = self.guide + self.tool
        if squared > outer * outer:
            return False
        if self.guide <= self.tool:
            return True
        inner = self.guide - self.tool
        return squared >= inner * inner


def _arrays(circles: Sequence[PaperCircleCandidate]) -> np.ndarray:
    return np.array([(float(c.center.x), float(c.center.y), float(c.guide_radius.value)) for c in circles], dtype=np.float64).reshape(-1, 3)


def _overlapping_earlier(geometry: np.ndarray, index: int, tool: float) -> np.ndarray:
    """Indices below ``index`` whose annulus may meet ``index``'s (binary64 selection, not a decision)."""
    x, y, rho = geometry[index]
    earlier = geometry[:index]
    distance = np.hypot(earlier[:, 0] - x, earlier[:, 1] - y)
    inner = np.maximum(rho - tool, 0.0)
    earlier_inner = np.maximum(earlier[:, 2] - tool, 0.0)
    reach = (distance <= earlier[:, 2] + rho + 2.0 * tool) & (distance >= earlier_inner - rho - tool) & (distance >= inner - earlier[:, 2] - tool)
    return np.flatnonzero(reach)


@dataclass(frozen=True)
class _DyadicImage:
    """The path's doubles as integers at one dyadic scale: value = integer / 2**scale, exactly.

    Every finite double is a dyadic rational, so one scale makes every centre
    coordinate, radius and the tool radius an integer, and squared distances and
    squared radii exact integers. Membership of an integer point is then decided
    exactly with integer comparisons.
    """

    scale: int
    x: tuple[int, ...]
    y: tuple[int, ...]
    outer: tuple[int, ...]
    inner: tuple[int, ...]

    @classmethod
    def build(cls, geometry: np.ndarray, tool: float) -> _DyadicImage:
        values = [float(v) for v in geometry.ravel()] + [tool]
        scale = max(v.as_integer_ratio()[1].bit_length() - 1 for v in values)

        def lift(v: float) -> int:
            numerator, denominator = v.as_integer_ratio()
            return numerator << (scale - denominator.bit_length() + 1)

        t = lift(tool)
        guides = [lift(float(rho)) for rho in geometry[:, 2]]
        return cls(
            scale,
            tuple(lift(float(v)) for v in geometry[:, 0]),
            tuple(lift(float(v)) for v in geometry[:, 1]),
            tuple((g + t) * (g + t) for g in guides),
            tuple((g - t) * (g - t) if g > t else -1 for g in guides),
        )

    def lift(self, value: float) -> tuple[int, int]:
        """(integer, extra) with ``value == integer / 2**(scale + extra)`` exactly; extra > 0 only for a double finer than the scale."""
        numerator, denominator = value.as_integer_ratio()
        extra = max(0, denominator.bit_length() - 1 - self.scale)
        return numerator << (self.scale + extra - denominator.bit_length() + 1), extra

    def is_rational_witness(self, index: int, selected: np.ndarray, px: Fraction, py: Fraction) -> bool:
        """The dyadic point (px, py) lies in closed annulus ``index`` and outside every selected closed annulus, exact.

        A point finer than the scale is checked against the image refined by the extra bits.
        """
        extra = max(0, max(px.denominator, py.denominator).bit_length() - 1 - self.scale)
        unit = 1 << (self.scale + extra)
        point = (px.numerator * unit // px.denominator, py.numerator * unit // py.denominator)
        if Fraction(point[0], unit) != px or Fraction(point[1], unit) != py:
            raise NonDyadicWitnessError(f"witness ({px}, {py}) is not dyadic.")

        def contains(j: int) -> bool:
            dx, dy = point[0] - (self.x[j] << extra), point[1] - (self.y[j] << extra)
            squared = dx * dx + dy * dy
            return self.inner[j] << (2 * extra) <= squared <= self.outer[j] << (2 * extra) if self.inner[j] >= 0 else squared <= self.outer[j] << (2 * extra)

        return contains(index) and not any(contains(int(j)) for j in selected)


def _clearance(px: np.ndarray, py: np.ndarray, other: np.ndarray, tool: float) -> np.ndarray:
    """Binary64 distance of each point outside its nearest selected annulus (negative inside). Proposal only."""
    distance = np.hypot(px[:, None] - other[None, :, 0], py[:, None] - other[None, :, 1])
    outside = np.maximum(distance - (other[None, :, 2] + tool), np.where(other[None, :, 2] > tool, (other[None, :, 2] - tool) - distance, -np.inf))
    return outside.min(axis=1) if other.shape[0] else np.full(px.shape, np.inf)


def _arc_bounds(bx: np.ndarray, by: np.ndarray, radius: np.ndarray, ax: np.ndarray, ay: np.ndarray, arho: np.ndarray, tool: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """For circles about (bx, by) of ``radius`` and closed annuli (ax, ay, arho), broadcast: (phi, alpha, beta).

    The circle lies inside the annulus at angles phi + u with beta <= |u| <= alpha;
    alpha < 0: never inside; beta > pi: never inside. Binary64.
    """
    dx, dy = ax - bx, ay - by
    distance = np.hypot(dx, dy)
    phi = np.arctan2(dy, dx)
    outer = arho + tool
    inner = np.maximum(arho - tool, 0.0)
    r = radius
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        cos_outer = (r * r + distance * distance - outer * outer) / (2.0 * r * distance)
        cos_inner = (r * r + distance * distance - inner * inner) / (2.0 * r * distance)
    concentric = distance == 0.0
    member = (inner <= r) & (r <= outer)
    alpha = np.where(cos_outer <= -1.0, np.pi, np.where(cos_outer > 1.0, -1.0, np.arccos(np.clip(cos_outer, -1.0, 1.0))))
    beta = np.where(cos_inner >= 1.0, 0.0, np.where(cos_inner < -1.0, 4.0, np.arccos(np.clip(cos_inner, -1.0, 1.0))))
    alpha = np.where(concentric, np.where(member, np.pi, -1.0), alpha)
    beta = np.where(concentric, 0.0, beta)
    return phi, alpha, beta


def _membership_intervals(phi: np.ndarray, alpha: np.ndarray, beta: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(starts, ends, valid) of the two angle intervals per (phi, alpha, beta), stacked on a new last axis. Binary64."""
    valid = alpha >= beta
    single = beta == 0.0
    starts = np.stack([np.where(single, phi - alpha, phi + beta), phi - alpha], axis=-1)
    ends = np.stack([phi + alpha, phi - beta], axis=-1)
    return starts, ends, np.stack([valid, valid & ~single], axis=-1)


def _widest_gaps(starts: np.ndarray, ends: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per row: the widest angle gap on [0, 2 pi) outside every interval (start, end), end - start < 2 pi; invalid rows hold start = 2 pi. Binary64."""
    turn = 2.0 * np.pi
    wrap = np.maximum((ends - turn).max(axis=1, initial=0.0), 0.0)
    order = np.argsort(starts, axis=1)
    starts = np.take_along_axis(starts, order, axis=1)
    ends = np.take_along_axis(ends, order, axis=1)
    reach = np.maximum.accumulate(np.concatenate([wrap[:, None], ends], axis=1), axis=1)
    lower = reach
    upper = np.concatenate([starts, np.full((starts.shape[0], 1), turn)], axis=1)
    width = upper - lower
    best = np.argmax(width, axis=1)
    rows = np.arange(starts.shape[0])
    return lower[rows, best], upper[rows, best]


def _rim_gaps(geometry: np.ndarray, index: int, selected: np.ndarray, tool: float, *, own_rims: bool) -> list[_RimGap]:
    """The widest uncovered gaps of the boundary circles met by annulus ``index``. Proposal only.

    The boundary of any uncovered part of the annulus lies on the annulus's
    own rims or on a selected annulus's rims, so each uncovered part of
    positive area shows as an angle gap on one of those circles: inside the
    target annulus and outside every other selected annulus. A point just off
    the gap, on the uncovered side, is a candidate (``_rim_gap_points``). ``own_rims`` searches the
    target's two rims only (the new material usually touches them, and the
    search is then linear in the selected annuli rather than quadratic).
    """
    x, y, rho = geometry[index]
    others = geometry[selected]
    m = others.shape[0]
    centres = np.concatenate([[[x, y, rho]], others])
    bx = np.concatenate([centres[:, 0], centres[:, 0]])
    by = np.concatenate([centres[:, 1], centres[:, 1]])
    radius = np.concatenate([centres[:, 2] + tool, centres[:, 2] - tool])
    owner = np.concatenate([np.arange(-1, m), np.arange(-1, m)])
    # Off the target's rims: into the band; off another annulus's rims: out of it.
    side = np.concatenate([[-1.0], np.ones(m), [1.0], -np.ones(m)])
    live = (radius > 0.0) & ((owner < 0) | (not own_rims))
    bx, by, radius, owner, side = bx[live], by[live], radius[live], owner[live], side[live]
    turn = 2.0 * np.pi
    starts, ends, use = _membership_intervals(*_arc_bounds(bx[:, None], by[:, None], radius[:, None], others[None, :, 0], others[None, :, 1], others[None, :, 2], tool))
    use &= (owner[:, None] != np.arange(m)[None, :])[:, :, None]
    starts, ends, use = starts.reshape(len(bx), -1), ends.reshape(len(bx), -1), use.reshape(len(bx), -1)
    # Outside the target annulus (rows on other annuli's rims only): covered as well.
    tphi, talpha, tbeta = _arc_bounds(bx, by, radius, x, y, rho, tool)
    foreign = owner >= 0
    starts = np.concatenate([starts, (tphi + talpha)[:, None], (tphi - tbeta)[:, None]], axis=1)
    ends = np.concatenate([ends, (tphi + turn - talpha)[:, None], (tphi + tbeta)[:, None]], axis=1)
    use = np.concatenate([use, (foreign & (talpha < np.pi))[:, None], (foreign & (tbeta > 0.0))[:, None]], axis=1)
    full = (use & (ends - starts >= turn)).any(axis=1)
    length = ends - starts
    starts = np.where(use, np.mod(starts, turn), turn)
    ends = np.where(use, starts + length, 0.0)
    lower, upper = _widest_gaps(starts, ends)
    score = np.where(full, 0.0, (upper - lower) * radius)
    gaps: list[_RimGap] = []
    for row in np.argsort(-score)[:RIM_GAPS_TRIED]:
        if score[row] <= 0.0:
            break
        gaps.append(_RimGap(float(bx[row]), float(by[row]), float(radius[row]), float(side[row]), int(owner[row]), 0.5 * float(lower[row] + upper[row])))
    return gaps


@dataclass(frozen=True)
class _RimGap:
    """An uncovered angle gap's midpoint on a boundary circle, and the side of the circle the uncovered part lies on (binary64)."""

    x: float
    y: float
    radius: float
    side: float
    owner: int
    angle: float


def _rim_gap_points(geometry: np.ndarray, index: int, selected: np.ndarray, tool: float, gaps: Sequence[_RimGap]) -> list[tuple[float, float]]:
    """Binary64 candidates just off each rim gap, by fractions of the binary64 clearance there. Proposal only."""
    target = geometry[index : index + 1]
    others = geometry[selected]
    points: list[tuple[float, float]] = []
    for gap in gaps:
        ux, uy = float(np.cos(gap.angle)), float(np.sin(gap.angle))
        qx, qy = np.array([gap.x + gap.radius * ux]), np.array([gap.y + gap.radius * uy])
        covering = others if gap.owner < 0 else np.delete(others, gap.owner, axis=0)
        room = float(_clearance(qx, qy, covering, tool)[0])
        # Depth inside the target annulus: its own rims offset into the band (at most half its width).
        depth = tool if gap.owner < 0 else -float(_clearance(qx, qy, target, tool)[0])
        reach = min(room, depth)
        if reach <= 0.0:
            continue
        for fraction in RIM_GAP_OFFSETS:
            offset = gap.radius + gap.side * fraction * reach
            points.append((gap.x + offset * ux, gap.y + offset * uy))
    return points


def _deep_rim_points(gaps: Sequence[_RimGap]) -> list[tuple[Fraction, Fraction]]:
    """Rational candidates off each rim gap by 2**-k of the radius, below binary64 resolution. Proposal only.

    A sliver thinner than one ulp of its coordinates has no binary64 point; a
    dyadic point between the rim and the sliver's far side can still exist.
    """
    points: list[tuple[Fraction, Fraction]] = []
    for gap in gaps:
        ux, uy = Fraction(float(np.cos(gap.angle))), Fraction(float(np.sin(gap.angle)))
        bx, by, radius = Fraction(gap.x), Fraction(gap.y), Fraction(gap.radius)
        for depth in DEEP_OFFSET_EXPONENTS:
            offset = radius * (1 + Fraction(int(gap.side), 1 << depth))
            points.append((bx + offset * ux, by + offset * uy))
    return points


def _witness(geometry: np.ndarray, image: _DyadicImage, index: int, selected: np.ndarray, tool: float, *, own_rims_tried: bool) -> tuple[Fraction, Fraction] | None:
    """An exact witness that annulus ``index`` leaves material outside the selected annuli, or None (then the exact test decides).

    Proposals, cheapest first: the own-rim gaps (skipped when the batched screen
    already tried them), then every boundary circle's gaps in binary64, then the
    same gaps at sub-ulp rational offsets.
    """
    if selected.size == 0:
        x, y, rho = geometry[index]
        point = Fraction(float(x)), Fraction(float(y + rho))
        return point if image.is_rational_witness(index, selected, *point) else None
    own = [] if own_rims_tried else _rim_gap_points(geometry, index, selected, tool, _rim_gaps(geometry, index, selected, tool, own_rims=True))
    for px, py in own:
        point = Fraction(px), Fraction(py)
        if image.is_rational_witness(index, selected, *point):
            return point
    gaps = _rim_gaps(geometry, index, selected, tool, own_rims=False)
    for px, py in _rim_gap_points(geometry, index, selected, tool, gaps):
        point = Fraction(px), Fraction(py)
        if image.is_rational_witness(index, selected, *point):
            return point
    for wx, wy in _deep_rim_points(gaps):
        if image.is_rational_witness(index, selected, wx, wy):
            return wx, wy
    return None


@dataclass(frozen=True)
class _Selection:
    """Every circle's selected earlier circles (``_overlapping_earlier``, the same binary64 expressions), as CSR rows."""

    offsets: np.ndarray
    columns: np.ndarray

    @classmethod
    def build(cls, geometry: np.ndarray, tool: float) -> _Selection:
        n = geometry.shape[0]
        inner = np.maximum(geometry[:, 2] - tool, 0.0)
        rows: list[np.ndarray] = []
        cols: list[np.ndarray] = []
        for start in range(0, n, SELECTION_BLOCK):
            stop = min(start + SELECTION_BLOCK, n)
            block = geometry[start:stop]
            earlier = geometry[:stop]
            distance = np.hypot(earlier[None, :, 0] - block[:, None, 0], earlier[None, :, 1] - block[:, None, 1])
            reach = (
                (distance <= earlier[None, :, 2] + block[:, None, 2] + 2.0 * tool)
                & (distance >= inner[None, :stop] - block[:, None, 2] - tool)
                & (distance >= inner[start:stop, None] - earlier[None, :, 2] - tool)
            )
            reach &= np.arange(stop)[None, :] < np.arange(start, stop)[:, None]
            row, col = np.nonzero(reach)
            rows.append(row + start)
            cols.append(col)
        row = np.concatenate(rows)
        offsets = np.zeros(n + 1, dtype=np.int64)
        np.add.at(offsets, row + 1, 1)
        return cls(np.cumsum(offsets), np.concatenate(cols))

    def of(self, index: int) -> np.ndarray:
        return self.columns[self.offsets[index] : self.offsets[index + 1]]

    def rows(self) -> np.ndarray:
        return np.repeat(np.arange(len(self.offsets) - 1), np.diff(self.offsets))


def _own_rim_candidates(geometry: np.ndarray, selection: _Selection, tool: float) -> list[np.ndarray]:
    """Per candidate slot, (n, 2) binary64 points (NaN: none) off the widest uncovered gap of each circle's own rims. Proposal only.

    The batched form of ``_rim_gap_points(..., own_rims=True)``: rows are
    (circle, rim), columns the circle's selected annuli; rows are offset by
    ``ROW_SPAN`` so one sort and one running maximum serve the whole path.
    """
    n = geometry.shape[0]
    turn = 2.0 * np.pi
    pair_row = selection.rows()
    pair_col = selection.columns
    slots: list[np.ndarray] = []
    for rim, side in ((1.0, -1.0), (-1.0, 1.0)):
        radius = geometry[:, 2] + rim * tool
        phi, alpha, beta = _arc_bounds(geometry[pair_row, 0], geometry[pair_row, 1], radius[pair_row], geometry[pair_col, 0], geometry[pair_col, 1], geometry[pair_col, 2], tool)
        starts, ends, use = _membership_intervals(phi, alpha, beta)
        owner = np.repeat(pair_row, 2)
        starts, ends, use = starts.ravel(), ends.ravel(), use.ravel()
        full = np.zeros(n, dtype=bool)
        full[owner[use & (ends - starts >= turn)]] = True
        owner, starts, length = owner[use], np.mod(starts[use], turn), (ends - starts)[use]
        ends = starts + length
        wrap = np.zeros(n)
        np.maximum.at(wrap, owner, ends - turn)
        # Per row: a head interval (-1, wrap) and a tail (2 pi, 2 pi) bracket the real ones.
        rows = np.concatenate([np.arange(n), owner, np.arange(n)])
        lo = np.concatenate([np.full(n, -1.0), starts, np.full(n, turn)]) + rows * ROW_SPAN
        hi = np.concatenate([wrap, ends, np.full(n, turn)]) + rows * ROW_SPAN
        order = np.argsort(lo, kind="stable")
        rows, lo, hi = rows[order], lo[order], hi[order]
        reach = np.maximum.accumulate(hi)
        width = np.full(lo.shape, -np.inf)
        width[1:] = lo[1:] - reach[:-1]
        width[np.concatenate([[True], rows[1:] != rows[:-1]])] = -np.inf
        best = np.full(n, -np.inf)
        np.maximum.at(best, rows, width)
        winner = np.nonzero((width == best[rows]) & (width > 0.0))[0]
        first = winner[np.unique(rows[winner], return_index=True)[1]]
        circle = rows[first]
        middle = 0.5 * (reach[first - 1] + lo[first]) - circle * ROW_SPAN
        live = ~full[circle] & (radius[circle] > 0.0)
        circle, middle = circle[live], middle[live]
        ux, uy = np.cos(middle), np.sin(middle)
        qx, qy = geometry[circle, 0] + radius[circle] * ux, geometry[circle, 1] + radius[circle] * uy
        room = _pair_clearance(geometry, selection, circle, qx, qy, tool)
        reach_in = np.minimum(room, tool)
        for fraction in RIM_GAP_OFFSETS:
            slot = np.full((n, 2), np.nan)
            ok = reach_in > 0.0
            offset = radius[circle] + side * fraction * reach_in
            slot[circle[ok], 0] = (geometry[circle, 0] + offset * ux)[ok]
            slot[circle[ok], 1] = (geometry[circle, 1] + offset * uy)[ok]
            slots.append(slot)
    return slots


def _pair_clearance(geometry: np.ndarray, selection: _Selection, circle: np.ndarray, qx: np.ndarray, qy: np.ndarray, tool: float) -> np.ndarray:
    """Binary64 clearance of point q[i] outside every annulus selected for circle[i] (inf when none selected)."""
    counts = np.diff(selection.offsets)[circle]
    owner = np.repeat(np.arange(len(circle)), counts)
    cols = np.concatenate([selection.of(int(c)) for c in circle]) if len(circle) else np.zeros(0, dtype=np.int64)
    distance = np.hypot(qx[owner] - geometry[cols, 0], qy[owner] - geometry[cols, 1])
    rho = geometry[cols, 2]
    outside = np.maximum(distance - (rho + tool), np.where(rho > tool, (rho - tool) - distance, -np.inf))
    room = np.full(len(circle), np.inf)
    np.minimum.at(room, owner, outside)
    return room


def _batch_witnessed(geometry: np.ndarray, image: _DyadicImage, selection: _Selection, tool: float) -> np.ndarray:
    """Per circle: an exact witness was found among its own-rim candidates (the exact check is integer, per point and per selected annulus)."""
    n = geometry.shape[0]
    witnessed = np.zeros(n, dtype=bool)
    witnessed[0] = True
    x, y = np.array(image.x, dtype=object), np.array(image.y, dtype=object)
    outer, inner = np.array(image.outer, dtype=object), np.array(image.inner, dtype=object)
    for slot in _own_rim_candidates(geometry, selection, tool):
        circle = np.nonzero(~witnessed & ~np.isnan(slot[:, 0]))[0]
        if circle.size == 0:
            continue
        # Each candidate double lifted exactly; the slot shares the finest scale any of them needs.
        lifted = [(image.lift(float(slot[c, 0])), image.lift(float(slot[c, 1]))) for c in circle]
        extra = max(max(ex, ey) for (_, ex), (_, ey) in lifted)
        px = np.array([nx << (extra - ex) for (nx, ex), _ in lifted], dtype=object)
        py = np.array([ny << (extra - ey) for _, (ny, ey) in lifted], dtype=object)
        unit, area = 1 << extra, 1 << (2 * extra)
        dx, dy = px - x[circle] * unit, py - y[circle] * unit
        squared = dx * dx + dy * dy
        inside = (squared <= outer[circle] * area) & (squared >= inner[circle] * area)
        counts = np.diff(selection.offsets)[circle]
        owner = np.repeat(np.arange(circle.size), counts)
        cols = np.concatenate([selection.of(int(c)) for c in circle])
        dx, dy = px[owner] - x[cols] * unit, py[owner] - y[cols] * unit
        squared = dx * dx + dy * dy
        covered = np.zeros(circle.size, dtype=np.int64)
        np.add.at(covered, owner, ((squared <= outer[cols] * area) & (squared >= inner[cols] * area)).astype(np.int64))
        witnessed[circle[inside.astype(bool) & (covered == 0)]] = True
    return witnessed


def material_witness(target: PaperCircleCandidate, earlier: Sequence[PaperCircleCandidate], tool_radius: ToolRadius) -> tuple[Fraction, Fraction] | None:
    """An exact rational point of ``target``'s swept annulus outside every overlapping earlier annulus, or None.

    A returned point proves ``target`` removes material from the union of
    ``earlier``'s annuli. None proves nothing: the exact union test decides.
    """
    geometry = _arrays((*earlier, target))
    tool = float(tool_radius.value)
    image = _DyadicImage.build(geometry, tool)
    return _witness(geometry, image, len(earlier), _overlapping_earlier(geometry, len(earlier), tool), tool, own_rims_tried=False)


def _seed_square(target: _Annulus) -> np.ndarray:
    x, y = float(target.x), float(target.y)
    half = SEED_SQUARE_FACTOR * float(target.guide + target.tool)
    corners = ((x - half, y - half), (x + half, y - half), (x + half, y + half), (x - half, y + half))
    outer = target.guide + target.tool
    xs = sorted({Fraction(cx) for cx, _ in corners})
    ys = sorted({Fraction(cy) for _, cy in corners})
    if xs[0] > target.x - outer or xs[-1] < target.x + outer or ys[0] > target.y - outer or ys[-1] < target.y + outer:
        raise SeedSquareContainmentError(f"seed square of half-width {half!r} does not contain the annulus about ({x!r}, {y!r}).")
    return np.array([(cx, cy, 0.0) for cx, cy in corners], dtype=np.float64)


def _union_covers(target: _Annulus, earlier: Sequence[_Annulus]) -> bool:
    """Exact: the target annulus lies in the union of ``earlier`` (in the order given).

    The NOTHING query runs after ``FIRST_UNION_QUERY`` subtractions, then after
    each doubling, and after the last. Subtracting only shrinks the stock, so a
    NOTHING verdict after a prefix is the verdict of the full union.
    """
    stock = _stock_2.Stock2(_seed_square(target), [])
    query = FIRST_UNION_QUERY
    for count, annulus in enumerate(earlier, start=1):
        stock.subtract_circle_sweep(float(annulus.x), float(annulus.y), float(annulus.guide), float(annulus.tool))
        if count == query and count < len(earlier):
            if _removes_nothing_from(stock, target):
                return True
            query *= 2
    return _removes_nothing_from(stock, target)


def _removes_nothing_from(stock: _stock_2.Stock2, target: _Annulus) -> bool:
    verdict = stock.circle_sweep_removal(float(target.x), float(target.y), float(target.guide), float(target.tool))
    return verdict == _stock_2.SweepRemovalVerdict.NOTHING


def _nearest_first(geometry: np.ndarray, index: int, selected: np.ndarray) -> np.ndarray:
    """``selected`` ordered by binary64 centre distance plus guide-radius difference to ``index`` (order only, not a decision)."""
    x, y, rho = geometry[index]
    other = geometry[selected]
    key = np.hypot(other[:, 0] - x, other[:, 1] - y) + np.abs(other[:, 2] - rho)
    return selected[np.argsort(key, kind="stable")]


def union_covers_annulus(target: PaperCircleCandidate, earlier: Sequence[PaperCircleCandidate], tool_radius: ToolRadius) -> bool:
    """Step 2 alone, exact: ``target``'s swept annulus lies in the union of ``earlier``'s swept annuli."""
    return _union_covers(_Annulus.of(target, tool_radius), [_Annulus.of(circle, tool_radius) for circle in earlier])


def removes_nothing(target: PaperCircleCandidate, earlier: Sequence[PaperCircleCandidate], tool_radius: ToolRadius) -> bool:
    """Exact sufficient test: True proves ``target``'s swept annulus lies in the union of ``earlier``'s.

    False means either material is left or the proof was not attempted (the
    selection and the screen only ever reject).
    """
    geometry = _arrays((*earlier, target))
    tool = float(tool_radius.value)
    annuli = [_Annulus.of(circle, tool_radius) for circle in (*earlier, target)]
    return _decide(geometry, _DyadicImage.build(geometry, tool), len(earlier), annuli, tool, _overlapping_earlier(geometry, len(earlier), tool), own_rims_tried=False)[0]


def _decide(geometry: np.ndarray, image: _DyadicImage, index: int, annuli: Sequence[_Annulus], tool: float, selected: np.ndarray, *, own_rims_tried: bool) -> tuple[bool, bool]:
    """(NOTHING proven, exact union test run) for circle ``index`` against its ``selected`` earlier circles."""
    if _witness(geometry, image, index, selected, tool, own_rims_tried=own_rims_tried) is not None:
        return False, False
    return _union_covers(annuli[index], [annuli[int(j)] for j in _nearest_first(geometry, index, selected)]), True


@dataclass(frozen=True)
class NothingCirclePass(Generic[C]):
    """The surviving circles and the indices (into the input) the pass removed or had to keep.

    ``kept_for_successor_cap``: proven NOTHING, kept because the successor
    would exceed the cap after the new predecessor (or the model has no
    engagement for that pair).
    """

    circles: tuple[C, ...]
    removed: tuple[int, ...]
    kept_for_successor_cap: tuple[int, ...]
    proven_nothing: int
    union_tests: int


def _successor_within_cap(predecessor: PaperCircleCandidate, successor: PaperCircleCandidate, tool_radius: ToolRadius, cap: EngagementCap) -> bool:
    try:
        return exact_predecessor_engagement(predecessor, successor, tool_radius).within(cap)
    except StandardPlacementFragmentationError:
        return False


def remove_nothing_circles(
    circles: Sequence[C],
    paper: Callable[[C], PaperCircleCandidate],
    tool_radius: ToolRadius,
    cap: EngagementCap,
    *,
    closed: bool,
) -> NothingCirclePass[C]:
    """Remove, in path order, every circle proven NOTHING whose successor stays within ``cap`` after the new predecessor.

    Earlier circles are all the input circles before k: a removed circle's
    annulus lies in the union of the kept ones before it, so that union equals
    the union of all earlier circles.

    Successor rule (exact, ``exact_predecessor_engagement(...).within(cap)``;
    a fragmentation error counts as not within). A circle k proven NOTHING is
    removed only when the pair (last kept circle before k, successor of k) is
    within ``cap``, where the successor of k is circle k+1 for k < n-1; for
    k = n-1 it is circle 0 when ``closed`` (the seam) and there is none when
    open. Circles are decided in path order, so when k = n-1 is decided the
    last kept circle already reflects every earlier removal; while circle n-1
    is kept the seam pair (n-1, 0) is unchanged by any removal before it, so
    the seam needs re-deciding only at k = n-1. When every circle between is
    removed the seam pair is (0, 0), decided by the same predicate.

    Circle 0 is never removed, open or closed: it is cut first, so it has no
    earlier union to be covered by, and on a closed path it anchors the seam
    (it is not re-cut; the seam ends on it).

    Args:
        circles: The circles in path order.
        paper: Each circle's paper candidate (centre, guide radius).
        tool_radius: The tool radius.
        cap: The engagement cap the successor (and the seam, when closed) must keep.
        closed: The path closes with a seam from the last circle back to circle 0.
            Required: every generated path is closed, and an open default
            would silently leave the seam undecided.
    """
    candidates = [paper(circle) for circle in circles]
    if len(candidates) < 2:
        return NothingCirclePass(tuple(circles), (), (), 0, 0)
    geometry = _arrays(candidates)
    annuli = [_Annulus.of(circle, tool_radius) for circle in candidates]
    tool = float(tool_radius.value)
    image = _DyadicImage.build(geometry, tool)
    selection = _Selection.build(geometry, tool)
    witnessed = _batch_witnessed(geometry, image, selection, tool)
    last = len(candidates) - 1
    kept = [0]
    removed: list[int] = []
    refused: list[int] = []
    union_tests = 0
    for index in range(1, len(candidates)):
        nothing, tested = (False, False) if witnessed[index] else _decide(geometry, image, index, annuli, tool, selection.of(index), own_rims_tried=True)
        union_tests += tested
        if not nothing:
            kept.append(index)
            continue
        successor = index + 1 if index < last else (0 if closed else None)
        if successor is not None and not _successor_within_cap(candidates[kept[-1]], candidates[successor], tool_radius, cap):
            refused.append(index)
            kept.append(index)
            continue
        removed.append(index)
    return NothingCirclePass(tuple(circles[k] for k in kept), tuple(removed), tuple(refused), len(removed) + len(refused), union_tests)
