"""Exact continuous engagement certificate for Held connectors (tiers 1 and 1b).

Doctrine and derivation: ``docs/held_connector_engagement.md``. While connector
``i`` runs, the remaining stock is a subset of the pocket minus the swept
annulus of ANY circle ``k <= i`` (stock only shrinks). The tool rim minus one
annulus is an outer arc (outside the disk ``(c_k, R_k)``, ``R_k = rho_k + r``)
and a hole arc (inside the open disk ``(c_k, h_k)``, ``h_k = rho_k - r``). With
``d = |p - c_k|`` the outer angle is non-decreasing in ``d`` and the hole angle
is largest where ``kappa_h`` is smallest, at an end of ``d``'s range or at
``d* = sqrt(r^2 - h^2)``; ``d`` is convex along a straight piece. So two or
three rational tests bound the engagement over the whole piece, continuously.

Every decision is an exact sign test on rationals: the binary64 coordinates
and radii are injected exactly, and the cap enters only as ``EngagementCap``'s
exact chord surrogate ``s = 4 sin^2(theta / 2)``, so ``C = 1 - s / 4`` is
``cos^2(theta / 2)`` in surrogate semantics. Every test is homogeneous in
length, so the assessment runs them on integers: all inputs on one common
dyadic grid (``_DyadicGrid``), the same exact decisions without ``Fraction``'s
gcd normalisation. Degrees are reported for humans
only (:meth:`ConnectorEngagementLemma.bound_deg`) and never feed a decision.

A refusal means "not proved", never "violated"; violations are tier 2's to
show (``tools/held_connector_engagement_tier2.py``).
"""

from __future__ import annotations

import math
import time
from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum
from fractions import Fraction
from typing import Self

import numpy as np
import numpy.typing as npt

from benchmarks.errors import BenchmarkError
from benchmarks.held_figure5_boundary_path import Figure5BoundaryConnectors
from benchmarks.held_figure5_boundary_path import Figure5BoundaryTransition
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_boundary_path import Figure5SeamClosure
from compas_cgal.motion import EngagementCap
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

# Bisection depth per piece in tier 1b: a work limit, not a tolerance.
# Exhaustion means "not certified", never a verdict. The design's value.
MAX_BISECTION_DEPTH = 12
# Tiers 1 + 1b over every connector of one workload, seconds, on the Apple
# M1 Max: set before code in the design for Monstera (4,366 connectors).
TIERS_BUDGET_S = 0.5
# The chord surrogate's domain: s = 4 sin^2(theta / 2) lies in (0, 4] for a
# cap theta in (0, pi]. Exact binary64 bounds, not tolerances.
CHORD_RATIO_MAX = 4.0

# An exact rational: a Fraction in mm, or an integer on a common dyadic grid
# (mm * 2**shift). Every lemma test is homogeneous in length, so either works.
Exact = Fraction | int
ExactPoint = tuple[Exact, Exact]
RationalPoint = tuple[Fraction, Fraction]
FloatArray = npt.NDArray[np.float64]


class ConnectorEngagementError(BenchmarkError):
    """Base failure of the connector-engagement certificate."""


class InvalidConnectorEngagementInputError(ConnectorEngagementError):
    """The cap surrogate or the path shape violates the certificate's contract."""


class UncertifiedConnectorEngagementError(ConnectorEngagementError):
    """At least one connector is left without a tier-1 or tier-1b certificate."""

    def __init__(self, connectors: tuple[int, ...]) -> None:
        self.connectors = connectors
        super().__init__(f"{len(connectors)} connector(s) without an engagement certificate: {list(connectors)}")


def _rational(point: Point2[WorldXY]) -> RationalPoint:
    return Fraction(float(point.x)), Fraction(float(point.y))


def _squared_distance(point: ExactPoint, centre: ExactPoint) -> Exact:
    dx, dy = point[0] - centre[0], point[1] - centre[1]
    return dx * dx + dy * dy


def foot_squared_distance(a: ExactPoint, b: ExactPoint, centre: ExactPoint) -> Exact:
    """Exact squared distance from ``centre`` to the closed segment ``a b`` (``a != b``).

    Interior foot: ``cross(w, e)^2 / |e|^2``, the one division, kept rational.
    """
    ex, ey = b[0] - a[0], b[1] - a[1]
    wx, wy = centre[0] - a[0], centre[1] - a[1]
    dot = wx * ex + wy * ey
    if dot <= 0:
        return wx * wx + wy * wy
    length_sq = ex * ex + ey * ey
    if dot >= length_sq:
        return _squared_distance(b, centre)
    cross = wx * ey - wy * ex
    return Fraction(cross * cross, length_sq)


@dataclass(frozen=True)
class LemmaCircle:
    """One cut circle as the lemma reads it: exact centre and exact guide radius."""

    centre: ExactPoint
    guide_radius: Exact

    @classmethod
    def from_circle(cls, circle: Figure5CounterclockwiseCircle) -> Self:
        return cls(_rational(circle.center), Fraction(float(circle.radius.value)))


@dataclass(frozen=True)
class ConnectorEngagementLemma:
    """The single-sweep bound with its exact decision; reported degrees are separate.

    ``tool_radius`` is the exact rational tool radius ``r``; ``cos_sq_half_cap``
    is ``C = 1 - s / 4`` for the cap's binary64 chord surrogate ``s``.
    """

    tool_radius: Exact
    cos_sq_half_cap: Fraction

    @classmethod
    def build(cls, tool: ToolRadius, cap: EngagementCap) -> Self:
        return cls.from_chord_ratio(tool, cap.chord_ratio)

    @classmethod
    def from_chord_ratio(cls, tool: ToolRadius, chord_ratio: float) -> Self:
        """Inject a chord surrogate ``s`` directly.

        Raises:
            InvalidConnectorEngagementInputError: ``s`` is not a finite binary64 in ``(0, 4]``.
        """
        if not math.isfinite(chord_ratio) or not 0.0 < chord_ratio <= CHORD_RATIO_MAX:
            raise InvalidConnectorEngagementInputError(f"chord surrogate must lie in (0, 4], got {chord_ratio!r}")
        return cls(Fraction(float(tool.value)), 1 - Fraction(chord_ratio) / 4)

    def outer_within(self, squared_distance: Exact, outer_radius: Exact) -> bool:
        """Outer arc at ``d^2 = squared_distance`` within the cap: ``N_O <= 0`` and ``N_O^2 >= 4 C D r^2``."""
        r_sq = self.tool_radius * self.tool_radius
        n = squared_distance + r_sq - outer_radius * outer_radius
        c = self.cos_sq_half_cap
        return n <= 0 and n * n * c.denominator >= 4 * c.numerator * squared_distance * r_sq

    def hole_within(self, squared_distance: Exact, inner_radius: Exact) -> bool:
        """Hole arc at ``d^2 = squared_distance > 0`` within the cap; ``D = 0`` is refused."""
        if squared_distance == 0:
            return False
        r_sq = self.tool_radius * self.tool_radius
        n = squared_distance + r_sq - inner_radius * inner_radius
        c = self.cos_sq_half_cap
        return n >= 0 and n * n * c.denominator >= 4 * c.numerator * squared_distance * r_sq

    def hole_interval_within(self, d_min_sq: Exact, d_max_sq: Exact, inner_radius: Exact) -> bool:
        """Hole arc within the cap for every ``d^2`` in ``[d_min_sq, d_max_sq]``: both ends and ``d*`` inside."""
        if not (self.hole_within(d_min_sq, inner_radius) and self.hole_within(d_max_sq, inner_radius)):
            return False
        r_sq = self.tool_radius * self.tool_radius
        star = r_sq - inner_radius * inner_radius
        if star > 0 and d_min_sq < star < d_max_sq:
            c = self.cos_sq_half_cap
            return star * c.denominator >= c.numerator * r_sq
        return True

    def piece_certified_by(self, a: ExactPoint, b: ExactPoint, circle: LemmaCircle) -> bool:
        """Exact: every engaged run along the straight piece ``a b`` is within the cap, bounded by ``circle``'s annulus alone."""
        r = self.tool_radius
        rho = circle.guide_radius
        d_max_sq = max(_squared_distance(a, circle.centre), _squared_distance(b, circle.centre))
        if not self.outer_within(d_max_sq, rho + r):
            return False
        if rho <= r:
            return True
        d_min_sq = foot_squared_distance(a, b, circle.centre)
        if d_min_sq >= rho * rho:
            return True
        return self.hole_interval_within(d_min_sq, d_max_sq, rho - r)

    def bound_deg(self, a: ExactPoint, b: ExactPoint, circle: LemmaCircle) -> float:
        """REPORTING ONLY: the lemma's bound over the piece, degrees (outer at the far end, hole at its extremes)."""
        r = float(self.tool_radius)
        rho = float(circle.guide_radius)
        d_max_sq = max(_squared_distance(a, circle.centre), _squared_distance(b, circle.centre))
        bound = outer_arc_deg(math.sqrt(float(d_max_sq)), r, rho + r)
        if rho > r:
            d_min_sq = foot_squared_distance(a, b, circle.centre)
            if d_min_sq < circle.guide_radius * circle.guide_radius:
                h = rho - r
                distances = [math.sqrt(float(d_min_sq)), math.sqrt(float(d_max_sq))]
                if r > h and float(d_min_sq) < r * r - h * h < float(d_max_sq):
                    distances.append(math.sqrt(r * r - h * h))
                bound = max(bound, *(hole_arc_deg(d, r, h) for d in distances))
        return bound


def outer_arc_deg(d: float, r: float, outer: float) -> float:
    """REPORTING ONLY: outer-arc angle, degrees, at distance ``d`` from a centre with outer radius ``outer``."""
    if d == 0.0:
        return 0.0
    kappa = (d * d + r * r - outer * outer) / (2 * d * r)
    return 0.0 if kappa <= -1 else 360.0 if kappa >= 1 else math.degrees(2 * math.pi - 2 * math.acos(kappa))


def hole_arc_deg(d: float, r: float, inner: float) -> float:
    """REPORTING ONLY: hole-arc angle, degrees, at distance ``d`` from a centre with hole radius ``inner``."""
    if d == 0.0:
        return 360.0 if r < inner else 0.0
    kappa = (d * d + r * r - inner * inner) / (2 * d * r)
    return 0.0 if kappa >= 1 else 360.0 if kappa <= -1 else math.degrees(2 * math.acos(kappa))


def prior_circle_indices(connector: int, circle_count: int) -> range:
    """Circles cut before connector ``connector``, most recent first: ``k <= i``.

    ``transitions[i]`` follows ``circles[i]``, so this is ``i, i - 1, ..., 0``.
    The seam is the last connector, ``i = N - 1``, so every circle is a prior.
    A circle after ``i`` would be a false certificate.
    """
    if not 0 <= connector < circle_count:
        raise InvalidConnectorEngagementInputError(f"connector {connector} has no circle among {circle_count}")
    return range(connector, -1, -1)


class ConnectorTier(Enum):
    """How a connector is certified, or that it is not."""

    STATIONARY = "stationary"
    TIER_1 = "tier 1"
    TIER_1B = "tier 1b"
    REFUSED = "refused"


@dataclass(frozen=True)
class SubPieceCertificate:
    """Piece ``piece`` between exact parameters ``start`` and ``end`` is certified by circle ``circle``."""

    piece: int
    start: Fraction
    end: Fraction
    circle: int


@dataclass(frozen=True)
class ConnectorVerdict:
    """One connector's verdict; ``sub_pieces`` cover every moving piece exactly when certified."""

    connector: int
    seam: bool
    tier: ConnectorTier
    pieces: int
    sub_pieces: tuple[SubPieceCertificate, ...]

    @property
    def certified(self) -> bool:
        return self.tier is not ConnectorTier.REFUSED


def moving_pieces(transition: Figure5BoundaryTransition) -> tuple[tuple[RationalPoint, RationalPoint], ...]:
    """The connector's non-degenerate straight pieces, exactly."""
    vertices = [_rational(point) for point in transition.samples]
    return tuple((a, b) for a, b in zip(vertices, vertices[1:]) if a != b)


class _DyadicGrid:
    """One integer grid for a whole path: every binary64 input times ``2**shift`` is an integer.

    Every lemma test is homogeneous in length (degree 2 or 4 on both sides), so
    the common scale cancels: integers decide exactly what the rationals do. The
    shift keeps ``MAX_BISECTION_DEPTH`` spare bits, so every bisection point
    ``a + (b - a) k / 2**depth`` is an integer too.
    """

    def __init__(self, values: Iterable[float]) -> None:
        bits = max((value.as_integer_ratio()[1].bit_length() - 1 for value in values), default=0)
        self.shift = bits + MAX_BISECTION_DEPTH
        self.scale = 1 << self.shift

    def length(self, value: float) -> int:
        numerator, denominator = value.as_integer_ratio()
        return numerator << (self.shift - (denominator.bit_length() - 1))

    def point(self, point: Point2[WorldXY]) -> tuple[int, int]:
        return self.length(float(point.x)), self.length(float(point.y))


BoolArray = npt.NDArray[np.bool_]


def _down(value: FloatArray) -> FloatArray:
    return np.nextafter(value, -np.inf)


def _up(value: FloatArray) -> FloatArray:
    return np.nextafter(value, np.inf)


def _below(value: float) -> float:
    return math.nextafter(value, -math.inf)


def _above(value: float) -> float:
    return math.nextafter(value, math.inf)


def _magnitude_bounds(low: FloatArray, high: FloatArray) -> tuple[FloatArray, FloatArray]:
    """Bounds of ``|v|`` for ``v`` in ``[low, high]``."""
    return np.where(low > 0, low, np.where(high < 0, -high, 0.0)), np.maximum(np.abs(low), np.abs(high))


def reach_bracket(centres: FloatArray, reach_low: FloatArray, reach_high: FloatArray, x: tuple[float, float], y: tuple[float, float]) -> tuple[BoolArray, BoolArray]:
    """``(provably beyond, provably within)``: ``|p - c|`` against an exact reach in ``[reach_low, reach_high]``, ``p`` in the box ``x × y``.

    Only correctly rounded IEEE ``+ - ×`` and ``nextafter`` are used: each
    round-to-nearest result is off by at most half an ulp, so nudging it one
    ulp outward brackets the exact value. No ``hypot``, ``sqrt`` or other libm
    call, so the brackets hold on every IEEE platform. A circle in neither
    mask is undecided and goes to an exact stage.
    """
    dx_low, dx_high = _down(centres[:, 0] - x[1]), _up(centres[:, 0] - x[0])
    dy_low, dy_high = _down(centres[:, 1] - y[1]), _up(centres[:, 1] - y[0])
    small_x, large_x = _magnitude_bounds(dx_low, dx_high)
    small_y, large_y = _magnitude_bounds(dy_low, dy_high)
    distance_low = _down(_down(small_x * small_x) + _down(small_y * small_y))
    distance_high = _up(_up(large_x * large_x) + _up(large_y * large_y))
    return distance_low > _up(reach_high * reach_high), distance_high <= _down(reach_low * reach_low)


class _CandidateIndex:
    """Exact filter of the circles that could certify a sub-piece: ``|c - m| <= rho + 2 r + h``.

    ``m`` is the sub-piece's midpoint and ``h`` its half-length. A circle
    certifies only if the rim does not lie wholly outside its outer disk, so
    both ends, hence the midpoint, are within ``rho + 2 r`` of its centre, and
    the half-length keeps the test a superset. The filter only chooses which
    exact questions to ask: dropping a circle can lose a certificate, never
    create one. The test is decided exactly and without any libm function:
    an outward-rounded binary64 bracket (:func:`reach_bracket`, with
    ``h <= (|dx| + |dy|) / 2`` for "beyond" and ``h >= 0`` for "within")
    decides the clear cases; the rest is decided on the gate's integer grid
    by squaring, ``d <= R + h`` iff ``d² - R² - h² <= 0`` or
    ``(d² - R² - h²)² <= 4 R² h²`` (all quantities doubled to stay integral).
    """

    def __init__(self, circles: tuple[Figure5CounterclockwiseCircle, ...], tool: ToolRadius, grid: _DyadicGrid) -> None:
        self._centres: FloatArray = np.array([(float(c.center.x), float(c.center.y)) for c in circles], dtype=np.float64).reshape(-1, 2)
        reach = np.array([float(c.radius.value) + 2 * float(tool.value) for c in circles], dtype=np.float64)
        self._reach_low, self._reach_high = _down(reach), _up(reach)
        self._exact_centres = [grid.point(c.center) for c in circles]
        tool_length = grid.length(float(tool.value))
        self._exact_reach = [grid.length(float(c.radius.value)) + 2 * tool_length for c in circles]
        self._scale = grid.scale

    def _float(self, value: int) -> tuple[float, float]:
        rounded = value / self._scale  # int / int: correctly rounded
        return _below(rounded), _above(rounded)

    def _within_exactly(self, k: int, a: tuple[int, int], b: tuple[int, int]) -> bool:
        cx, cy = self._exact_centres[k]
        distance_sq = (2 * cx - a[0] - b[0]) ** 2 + (2 * cy - a[1] - b[1]) ** 2
        reach_sq = (2 * self._exact_reach[k]) ** 2
        half_sq = (b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2
        excess = distance_sq - reach_sq - half_sq
        return excess <= 0 or excess * excess <= 4 * reach_sq * half_sq

    def near(self, a: tuple[int, int], b: tuple[int, int], prior: range) -> list[int]:
        stop = max(prior[0], prior[-1]) + 1
        centres = self._centres[:stop]
        ax, ay, bx, by = (self._float(value) for value in (*a, *b))
        mx = (_below(ax[0] + bx[0]) / 2, _above(ax[1] + bx[1]) / 2)
        my = (_below(ay[0] + by[0]) / 2, _above(ay[1] + by[1]) / 2)
        width = _above(max(abs(bx[1] - ax[0]), abs(bx[0] - ax[1])))
        height = _above(max(abs(by[1] - ay[0]), abs(by[0] - ay[1])))
        half_high = _above(width + height) / 2
        beyond, _ = reach_bracket(centres, self._reach_low[:stop], _up(self._reach_high[:stop] + half_high), mx, my)
        _, within = reach_bracket(centres, self._reach_low[:stop], self._reach_high[:stop], mx, my)
        kept = np.nonzero(~beyond)[0][::-1]
        return [int(k) for k in kept if int(k) in prior and (within[k] or self._within_exactly(int(k), a, b))]


# Bisection parameters are integers k / 2**MAX_BISECTION_DEPTH.
_FULL_PIECE = 1 << MAX_BISECTION_DEPTH


def _at(a: tuple[int, int], b: tuple[int, int], k: int) -> tuple[int, int]:
    """``a + (b - a) k / 2**MAX_BISECTION_DEPTH``, exact: the grid keeps those bits spare."""
    return a[0] + (((b[0] - a[0]) * k) >> MAX_BISECTION_DEPTH), a[1] + (((b[1] - a[1]) * k) >> MAX_BISECTION_DEPTH)


def _bisect(
    lemma: ConnectorEngagementLemma,
    circles: tuple[LemmaCircle, ...],
    index: _CandidateIndex,
    prior: range,
    piece: int,
    a: tuple[int, int],
    b: tuple[int, int],
) -> tuple[SubPieceCertificate, ...] | None:
    """Tier 1b on one piece: any prior circle per sub-piece, exact dyadic midpoints, depth-first."""
    certified: list[SubPieceCertificate] = []
    stack: list[tuple[int, int, int]] = [(0, _FULL_PIECE, 0)]
    while stack:
        start, end, depth = stack.pop()
        sub_a, sub_b = _at(a, b, start), _at(a, b, end)
        winner = next((k for k in index.near(sub_a, sub_b, prior) if lemma.piece_certified_by(sub_a, sub_b, circles[k])), None)
        if winner is not None:
            certified.append(SubPieceCertificate(piece, Fraction(start, _FULL_PIECE), Fraction(end, _FULL_PIECE), winner))
            continue
        if depth >= MAX_BISECTION_DEPTH:
            return None
        middle = (start + end) >> 1
        stack.append((middle, end, depth + 1))
        stack.append((start, middle, depth + 1))
    return tuple(certified)


def _inputs(circles: tuple[Figure5CounterclockwiseCircle, ...], connectors: Figure5BoundaryConnectors, tool: ToolRadius) -> Iterable[float]:
    yield float(tool.value)
    for circle in circles:
        yield from (float(circle.center.x), float(circle.center.y), float(circle.radius.value))
    for transition in connectors.transitions:
        for point in transition.samples:
            yield from (float(point.x), float(point.y))


def _closed(circles: tuple[Figure5CounterclockwiseCircle, ...], connectors: Figure5BoundaryConnectors) -> bool:
    if not circles:
        raise InvalidConnectorEngagementInputError("a path requires at least one circle")
    expected = {Figure5SeamClosure.CONNECTOR: len(circles), Figure5SeamClosure.COINCIDENT: len(circles) - 1}[connectors.seam]
    if len(connectors.transitions) != expected:
        raise InvalidConnectorEngagementInputError(f"a {connectors.seam.value} seam over {len(circles)} circles takes {expected} connectors, not {len(connectors.transitions)}")
    for k, transition in enumerate(connectors.transitions):
        successor = circles[(k + 1) % len(circles)]
        if len(transition.samples) < 2 or transition.samples[0] != circles[k].contact_point or transition.samples[-1] != successor.contact_point:
            raise InvalidConnectorEngagementInputError(f"connector {k} does not run between the contacts of circles {k} and {(k + 1) % len(circles)}")
    return connectors.seam is Figure5SeamClosure.CONNECTOR


@dataclass(frozen=True)
class ConnectorEngagementAssessment:
    """Tier 1 and 1b verdicts for every connector, refusals included; never raises on a refusal."""

    verdicts: tuple[ConnectorVerdict, ...]
    seconds: float

    @classmethod
    def build(
        cls,
        circles: tuple[Figure5CounterclockwiseCircle, ...],
        connectors: Figure5BoundaryConnectors,
        tool: ToolRadius,
        cap: EngagementCap,
    ) -> Self:
        """Decide every connector: tier 1 with the circle just cut, else tier 1b over every prior circle.

        Raises:
            InvalidConnectorEngagementInputError: the connector count does not
                match the seam, or a connector does not join its circles' contacts.
        """
        started = time.perf_counter()
        closed = _closed(circles, connectors)
        surrogate = ConnectorEngagementLemma.build(tool, cap)
        grid = _DyadicGrid(_inputs(circles, connectors, tool))
        lemma = ConnectorEngagementLemma(grid.length(float(tool.value)), surrogate.cos_sq_half_cap)
        exact = tuple(LemmaCircle(grid.point(circle.center), grid.length(float(circle.radius.value))) for circle in circles)
        index = _CandidateIndex(circles, tool, grid)
        seam_index = len(circles) - 1 if closed else None
        verdicts = []
        for k, transition in enumerate(connectors.transitions):
            vertices = [grid.point(point) for point in transition.samples]
            pieces = [(a, b) for a, b in zip(vertices, vertices[1:]) if a != b]
            seam = k == seam_index
            if not pieces:
                verdicts.append(ConnectorVerdict(k, seam, ConnectorTier.STATIONARY, 0, ()))
                continue
            first_try = [lemma.piece_certified_by(a, b, exact[k]) for a, b in pieces]
            if all(first_try):
                verdicts.append(
                    ConnectorVerdict(k, seam, ConnectorTier.TIER_1, len(pieces), tuple(SubPieceCertificate(p, Fraction(0), Fraction(1), k) for p in range(len(pieces))))
                )
                continue
            prior = prior_circle_indices(k, len(circles))
            sub_pieces: list[SubPieceCertificate] = []
            refused = False
            for p, ((a, b), certified_by_own) in enumerate(zip(pieces, first_try)):
                if certified_by_own:
                    sub_pieces.append(SubPieceCertificate(p, Fraction(0), Fraction(1), k))
                    continue
                found = _bisect(lemma, exact, index, prior, p, a, b)
                if found is None:
                    refused = True
                    break
                sub_pieces.extend(found)
            tier = ConnectorTier.REFUSED if refused else ConnectorTier.TIER_1B
            verdicts.append(ConnectorVerdict(k, seam, tier, len(pieces), () if refused else tuple(sub_pieces)))
        return cls(tuple(verdicts), time.perf_counter() - started)

    @property
    def refused(self) -> tuple[int, ...]:
        return tuple(verdict.connector for verdict in self.verdicts if not verdict.certified)

    def count(self, tier: ConnectorTier) -> int:
        return sum(1 for verdict in self.verdicts if verdict.tier is tier)


@dataclass(frozen=True)
class ConnectorEngagementCertificate:
    """Every connector of a path certified within the cap, continuously, by tier 1 or 1b."""

    verdicts: tuple[ConnectorVerdict, ...]

    @classmethod
    def build(
        cls,
        circles: tuple[Figure5CounterclockwiseCircle, ...],
        connectors: Figure5BoundaryConnectors,
        tool: ToolRadius,
        cap: EngagementCap,
    ) -> Self:
        """Post-gate: certify all N connectors (N - 1 with a ``COINCIDENT`` seam), the seam included.

        Raises:
            InvalidConnectorEngagementInputError: see :meth:`ConnectorEngagementAssessment.build`.
            UncertifiedConnectorEngagementError: a connector is left without a
                certificate; it names every such connector.
        """
        assessment = ConnectorEngagementAssessment.build(circles, connectors, tool, cap)
        if assessment.refused:
            raise UncertifiedConnectorEngagementError(assessment.refused)
        return cls(assessment.verdicts)
