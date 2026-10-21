"""Pure predecessor-only placement for the Held Figure 5(a) standard path.

The paper's standard model evaluates a candidate machining circle only against
the immediately preceding swept disk.  This module implements that local model
on already transformed ``(c, rho, q)`` circles.  It neither owns nor consults a
global depleted-stock contour.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from fractions import Fraction
from typing import Self
from typing import Sequence

from compas.tolerance import TOL

from benchmarks.errors import BenchmarkError
from compas_cgal import _circle_geometry_2
from compas_cgal.motion import EngagementCap
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

# Held and Pfeiffer, Section 2.4: stop once theta_max - 0.001 <= theta <= theta_max.
# Decided exactly, never on acos: the band floor is the declared binary64
# ``cap.theta - 0.001`` (one IEEE subtraction of the two declared parameters),
# compared through its EngagementCap chord-ratio surrogate 4 sin^2(floor / 2).
# The sub-ulp gap to the transcendental band edge is that API semantic.
PAPER_BISECTION_ANGLE_TOLERANCE = Radian(0.001)


class InvalidPaperCircleCandidateError(ValueError):
    """A transformed ``(c, rho, q)`` tuple does not define one circle."""


class InvalidStandardPlacementInputError(ValueError):
    """The predecessor placement model received invalid typed input."""


class StandardPlacementResolutionError(RuntimeError):
    """The ordered candidates cannot resolve the paper's next placement."""


class StandardPlacementFragmentationError(RuntimeError):
    """The next candidate violates the paper's swept-disk overlap bound."""


class MixedRadicandError(BenchmarkError):
    """Quadratic-surd arithmetic was attempted across two different radicands."""


class ExactEngagementCosineRangeError(BenchmarkError):
    """An exact standard-model engagement cosine left [-1, 1], contradicting its branch."""


@dataclass(frozen=True)
class PaperCircleCandidate:
    """One already transformed paper circle ``(c, rho, q)`` in world XY."""

    center: Point2[WorldXY]
    guide_radius: GuideRadius
    contact_point: Point2[WorldXY]

    def __post_init__(self) -> None:
        if type(self.center) is not Point2 or type(self.guide_radius) is not GuideRadius or type(self.contact_point) is not Point2:
            raise InvalidPaperCircleCandidateError("Paper circle requires world-XY c/q points and one typed guide radius.")
        radius = float(self.guide_radius.value)
        if not TOL.is_between(
            math.dist(
                (float(self.center.x), float(self.center.y)),
                (float(self.contact_point.x), float(self.contact_point.y)),
            ),
            radius,
            radius,
        ):
            raise InvalidPaperCircleCandidateError("Paper circle contact q must lie on the circle (c, rho).")

    @classmethod
    def build(
        cls,
        *,
        center: Point2[WorldXY],
        guide_radius: GuideRadius,
        contact_point: Point2[WorldXY],
    ) -> Self:
        return cls(center, guide_radius, contact_point)


@dataclass(frozen=True)
class StandardPlacement:
    """Selected successor and paper-derived spacing diagnostics."""

    candidate: PaperCircleCandidate
    candidate_advance: int
    center_spacing: Millimetre
    guide_progress: Millimetre
    maximum_engagement: Radian
    cap_shortfall: Radian
    cap_excess: Radian
    overlap_margin: Millimetre
    reached_terminal: bool
    forced: bool


def _radius(candidate: PaperCircleCandidate) -> float:
    return float(candidate.guide_radius.value)


def _center_spacing(first: PaperCircleCandidate, second: PaperCircleCandidate) -> float:
    return math.dist(
        (float(first.center.x), float(first.center.y)),
        (float(second.center.x), float(second.center.y)),
    )


def _successor_disk_is_cleared(predecessor: PaperCircleCandidate, candidate: PaperCircleCandidate) -> bool:
    # The common tool radius cancels; CGAL owns the exact tangency decision.
    return bool(
        _circle_geometry_2.disk_contains_disk(
            (float(predecessor.center.x), float(predecessor.center.y)),
            _radius(predecessor),
            (float(candidate.center.x), float(candidate.center.y)),
            _radius(candidate),
        )
    )


def _resolvability_floor(
    first: PaperCircleCandidate,
    second: PaperCircleCandidate,
    tool_radius: ToolRadius,
) -> float:
    values = (
        float(first.center.x),
        float(first.center.y),
        float(second.center.x),
        float(second.center.y),
        _radius(first),
        _radius(second),
        float(tool_radius.value),
    )
    return max(math.ulp(abs(value)) for value in values)


def _overlap_margin(
    predecessor: PaperCircleCandidate,
    candidate: PaperCircleCandidate,
    tool_radius: ToolRadius,
) -> float:
    return 2.0 * float(tool_radius.value) - _center_spacing(predecessor, candidate) - _radius(candidate) + _radius(predecessor)


def predecessor_overlap_margin(
    predecessor: PaperCircleCandidate,
    candidate: PaperCircleCandidate,
    tool_radius: ToolRadius,
) -> Millimetre:
    """Return the paper Eq. 4 overlap margin for one ordered pair."""
    if type(predecessor) is not PaperCircleCandidate or type(candidate) is not PaperCircleCandidate or type(tool_radius) is not ToolRadius:
        raise InvalidStandardPlacementInputError("Overlap margin requires two paper circles and a tool radius.")
    return Millimetre(_overlap_margin(predecessor, candidate, tool_radius))


def _clamped_unit(value: float) -> float:
    lower_construction_bound = math.nextafter(-1.0, -math.inf)
    upper_construction_bound = math.nextafter(1.0, math.inf)
    if not lower_construction_bound <= value <= upper_construction_bound:
        raise StandardPlacementFragmentationError("Paper circle construction produced an invalid angular ratio.")
    return min(1.0, max(-1.0, value))


def float_predecessor_engagement(
    predecessor: PaperCircleCandidate,
    candidate: PaperCircleCandidate,
    tool_radius: ToolRadius,
) -> Radian:
    """Binary64 evaluation of Section 2.3's standard-model maximum engagement.

    Ablation witness only: its branch decisions are taken in floating point
    and misread tangencies (q_y^2 rounded below zero reports pi). No decision
    or report uses it; ``maximum_predecessor_engagement`` reports through
    ``exact_predecessor_engagement``.
    """
    if type(predecessor) is not PaperCircleCandidate or type(candidate) is not PaperCircleCandidate or type(tool_radius) is not ToolRadius:
        raise InvalidStandardPlacementInputError("Predecessor engagement requires two paper circles and a tool radius.")
    distance = _center_spacing(predecessor, candidate)
    tool = float(tool_radius.value)
    previous_radius = _radius(predecessor)
    radius = _radius(candidate)
    if _successor_disk_is_cleared(predecessor, candidate):
        return Radian(0.0)
    if distance == 0.0:
        previous_outer_radius = previous_radius + tool
        if radius - tool >= previous_outer_radius:
            return Radian(math.pi)
        # Concentric growth still cuts a new band. Equation 7 applies directly
        # with b = previous swept radius; the displaced-circle correction below
        # would divide by the zero center spacing.
        cosine = (previous_outer_radius * previous_outer_radius - tool * tool - radius * radius) / (2.0 * tool * radius)
        return Radian(math.acos(_clamped_unit(cosine)))

    if _overlap_margin(predecessor, candidate, tool_radius) <= 0.0:
        return Radian(math.pi)

    previous_outer_radius = previous_radius + tool
    current_outer_radius = radius + tool
    b_coordinate = previous_outer_radius - distance
    if b_coordinate <= 0.0:
        return Radian(math.pi)

    q_x = (b_coordinate * b_coordinate - tool * tool + radius * radius) / (2.0 * b_coordinate)
    q_y_squared = radius * radius - q_x * q_x
    if q_y_squared < 0.0:
        return Radian(math.pi)

    outer_scale = current_outer_radius / radius
    w_to_previous_squared = current_outer_radius * current_outer_radius + 2.0 * distance * outer_scale * q_x + distance * distance
    if w_to_previous_squared > previous_outer_radius * previous_outer_radius:
        cosine = (b_coordinate * b_coordinate - tool * tool - radius * radius) / (2.0 * tool * radius)
        return Radian(math.acos(_clamped_unit(cosine)))

    # Keep the corrected chord geometry in CGAL: expanding its squared
    # distance cancels the square roots before the final reporting conversion.
    try:
        cosine_squared = _circle_geometry_2.corrected_engagement_cosine_squared(
            (float(predecessor.center.x), float(predecessor.center.y)),
            previous_radius,
            (float(candidate.center.x), float(candidate.center.y)),
            radius,
            tool,
        )
    except _circle_geometry_2.NoCircleIntersectionError as error:
        raise StandardPlacementFragmentationError("Paper overlap correction has no real swept-disk intersection.") from error
    return Radian(2.0 * math.acos(math.sqrt(cosine_squared)))


@dataclass(frozen=True)
class _Surd:
    """Exact ``rational + radical * sqrt(radicand)`` over the integers, one fixed radicand >= 0."""

    rational: int
    radical: int
    radicand: int

    def _same(self, other: _Surd) -> None:
        if self.radicand != other.radicand:
            raise MixedRadicandError("Surd arithmetic requires one shared radicand.")

    def __add__(self, other: _Surd) -> _Surd:
        self._same(other)
        return _Surd(self.rational + other.rational, self.radical + other.radical, self.radicand)

    def __sub__(self, other: _Surd) -> _Surd:
        self._same(other)
        return _Surd(self.rational - other.rational, self.radical - other.radical, self.radicand)

    def __mul__(self, other: _Surd) -> _Surd:
        self._same(other)
        return _Surd(
            self.rational * other.rational + self.radical * other.radical * self.radicand,
            self.rational * other.radical + self.radical * other.rational,
            self.radicand,
        )

    def scaled(self, factor: int) -> _Surd:
        return _Surd(self.rational * factor, self.radical * factor, self.radicand)

    def plus(self, value: int) -> _Surd:
        return _Surd(self.rational + value, self.radical, self.radicand)

    def sign(self) -> int:
        """Exact sign: equal-sign parts decide directly, opposite signs compare squares."""
        rational = (self.rational > 0) - (self.rational < 0)
        radical = 0 if self.radicand == 0 else (self.radical > 0) - (self.radical < 0)
        if radical == 0 or rational == radical:
            return rational if rational != 0 else radical
        if rational == 0:
            return radical
        squares = self.rational * self.rational - self.radical * self.radical * self.radicand
        return rational if squares > 0 else radical if squares < 0 else 0

    def reported_ratio(self, denominator: int) -> float:
        # Reporting conversion only; never feeds a decision. Fractions keep
        # huge scaled integers out of float overflow.
        magnitude = math.sqrt(float(Fraction(self.radical * self.radical * self.radicand, denominator * denominator)))
        radical = -magnitude if self.radical < 0 else magnitude
        return float(Fraction(self.rational, denominator)) + radical


@dataclass(frozen=True)
class PredecessorEngagementCosine:
    """Exact cosine ``numerator / denominator`` of the Section 2.3 maximum engagement.

    The branches mirror ``float_predecessor_engagement``; every branch
    decision and the cap verdict are exact over the injected doubles. Only
    ``reported_angle`` leaves exact arithmetic.
    """

    numerator: _Surd
    denominator: int  # > 0

    def _surrogate_sign(self, cap: EngagementCap) -> int:
        """Exact sign of ``cos theta - (1 - chord_ratio / 2)``: +1 iff theta is below ``cap``.

        ``cap.chord_ratio`` is the exact surrogate ``4 sin^2(cap / 2)`` in
        (0, 4]; cosine is decreasing on [0, pi], so the comparison is exact
        for the surrogate. The sub-ulp gap to the transcendental angle is the
        EngagementCap API semantic, not a correction here.
        """
        ratio, ratio_denominator = cap.chord_ratio.as_integer_ratio()
        return self.numerator.scaled(2 * ratio_denominator).plus(-self.denominator * (2 * ratio_denominator - ratio)).sign()

    def within(self, cap: EngagementCap) -> bool:
        """Decide ``theta <= cap`` as ``cos theta >= 1 - chord_ratio / 2``."""
        return self._surrogate_sign(cap) >= 0

    def reaches(self, floor: EngagementCap) -> bool:
        """Decide ``theta >= floor`` as ``cos theta <= 1 - chord_ratio / 2``."""
        return self._surrogate_sign(floor) <= 0

    def reported_angle(self) -> Radian:
        return Radian(math.acos(min(1.0, max(-1.0, self.numerator.reported_ratio(self.denominator)))))


def _common_integers(values: tuple[float, ...]) -> tuple[int, ...]:
    """Exactly inject doubles as integers over their common power-of-two denominator."""
    ratios = tuple(value.as_integer_ratio() for value in values)
    scale = max(denominator for _, denominator in ratios)
    return tuple(numerator * (scale // denominator) for numerator, denominator in ratios)


def exact_predecessor_engagement(
    predecessor: PaperCircleCandidate,
    candidate: PaperCircleCandidate,
    tool_radius: ToolRadius,
) -> PredecessorEngagementCosine:
    """Exact standard-model engagement cosine for one ordered pair.

    Inputs are dyadic, so one power-of-two scale makes every length an
    integer; the model is homogeneous, so the scale cancels from each test.
    Center spacing d = sqrt(D) is carried symbolically: every Equation 7
    quantity lies in Z[sqrt(D)], and the Fig. 4(d) correction is rational,
    as in ``corrected_engagement_cosine_squared``. Positive denominators are
    multiplied through, never divided.

    Raises:
        StandardPlacementFragmentationError: The corrected swept disks have no
            real intersection (same failure as ``float_predecessor_engagement``).
        ExactEngagementCosineRangeError: An exact cosine left [-1, 1].
    """
    px, py, cx, cy, previous_radius, radius, tool = _common_integers(
        (
            float(predecessor.center.x),
            float(predecessor.center.y),
            float(candidate.center.x),
            float(candidate.center.y),
            float(predecessor.guide_radius.value),
            float(candidate.guide_radius.value),
            float(tool_radius.value),
        )
    )
    squared_spacing = (cx - px) ** 2 + (cy - py) ** 2
    zero = _Surd(0, 0, squared_spacing)
    # Successor disk inside the predecessor disk (tool radius cancels): theta = 0.
    if previous_radius >= radius and squared_spacing <= (previous_radius - radius) ** 2:
        return PredecessorEngagementCosine(zero.plus(1), 1)
    full_slot = PredecessorEngagementCosine(zero.plus(-1), 1)
    previous_outer = previous_radius + tool
    if squared_spacing == 0:
        if radius - tool >= previous_outer:
            return full_slot
        return _checked(zero.plus(previous_outer * previous_outer - tool * tool - radius * radius), 2 * tool * radius)
    spacing = _Surd(0, 1, squared_spacing)
    # Eq. 4 overlap margin 2r - d - rho + rho_prev <= 0.
    if (zero.plus(2 * tool - radius + previous_radius) - spacing).sign() <= 0:
        return full_slot
    b_coordinate = zero.plus(previous_outer) - spacing
    if b_coordinate.sign() <= 0:
        return full_slot
    # n = 2 b q_x; q_y^2 < 0 iff |n| > 2 b rho (b > 0).
    contact_numerator = (b_coordinate * b_coordinate).plus(radius * radius - tool * tool)
    bound = b_coordinate.scaled(2 * radius)
    if (contact_numerator - bound).sign() > 0 or (contact_numerator + bound).sign() < 0:
        return full_slot
    outer = radius + tool
    # |w - c_prev|^2 > (rho_prev + r)^2, multiplied through by rho * b > 0:
    # (W^2 + D - P^2) rho b + W d n > 0.
    outside = b_coordinate.scaled((outer * outer + squared_spacing - previous_outer * previous_outer) * radius) + (spacing * contact_numerator).scaled(outer)
    if outside.sign() > 0:
        return _checked((b_coordinate * b_coordinate).plus(-tool * tool - radius * radius), 2 * tool * radius)
    numerator = previous_outer * previous_outer - outer * outer - squared_spacing
    if 4 * squared_spacing * outer * outer - numerator * numerator < 0:
        raise StandardPlacementFragmentationError("Paper overlap correction has no real swept-disk intersection.")
    # Q = W q^2 = W rho^2 + rho n + W D; chord C = W (q^2 + r^2 - P^2).
    weighted_q_squared = outer * radius * radius + radius * numerator + outer * squared_spacing
    if weighted_q_squared <= 0:
        raise StandardPlacementFragmentationError("Paper overlap correction has no chord direction.")
    chord = weighted_q_squared + outer * (tool * tool - previous_outer * previous_outer)
    # cos(theta / 2)^2 = C^2 / (4 r^2 Q W); cos theta = 2 cos(theta / 2)^2 - 1.
    denominator = 4 * tool * tool * weighted_q_squared * outer
    if chord * chord > denominator:
        raise StandardPlacementFragmentationError("Paper overlap correction has no real contour intersection.")
    return _checked(zero.plus(2 * chord * chord - denominator), denominator)


def _checked(numerator: _Surd, denominator: int) -> PredecessorEngagementCosine:
    if numerator.plus(-denominator).sign() > 0 or numerator.plus(denominator).sign() < 0:
        raise ExactEngagementCosineRangeError("Exact standard-model engagement cosine is outside [-1, 1].")
    return PredecessorEngagementCosine(numerator, denominator)


def maximum_predecessor_engagement(
    predecessor: PaperCircleCandidate,
    candidate: PaperCircleCandidate,
    tool_radius: ToolRadius,
) -> Radian:
    """Report Section 2.3's standard-model maximum engagement.

    Branches are decided exactly by ``exact_predecessor_engagement``; only
    the final ``acos`` is binary64. All circle geometry is taken from the
    transformed candidate; no baseline radius or stock contour is substituted.

    Raises:
        InvalidStandardPlacementInputError: Untyped circles or tool radius.
        StandardPlacementFragmentationError: The corrected swept disks have no
            real intersection.
    """
    if type(predecessor) is not PaperCircleCandidate or type(candidate) is not PaperCircleCandidate or type(tool_radius) is not ToolRadius:
        raise InvalidStandardPlacementInputError("Predecessor engagement requires two paper circles and a tool radius.")
    return exact_predecessor_engagement(predecessor, candidate, tool_radius).reported_angle()


@dataclass(frozen=True)
class _BisectionStopBand:
    """Section 2.4 stopping band ``[cap - PAPER_BISECTION_ANGLE_TOLERANCE, cap]``, decided exactly.

    ``floor`` is the declared binary64 band floor as an EngagementCap
    surrogate, or None when that floor is not a positive angle: every
    engagement (theta >= 0) then reaches the band.
    """

    floor: EngagementCap | None

    @classmethod
    def build(cls, cap: EngagementCap) -> Self:
        floor = float(cap.theta) - float(PAPER_BISECTION_ANGLE_TOLERANCE)
        return cls(EngagementCap.build(floor) if floor > 0.0 else None)

    def stops(self, cosine: PredecessorEngagementCosine) -> bool:
        """Decide ``theta >= floor`` exactly; the caller has already decided ``theta <= cap``."""
        return self.floor is None or cosine.reaches(self.floor)


def _validated_spacings(
    predecessor: PaperCircleCandidate,
    candidates: Sequence[PaperCircleCandidate],
    tool_radius: ToolRadius,
) -> tuple[float, ...]:
    """Centre spacing of every consecutive ordered pair, each one resolvable."""
    if not candidates:
        raise InvalidStandardPlacementInputError("Standard placement requires at least one successor candidate.")
    spacings = []
    previous = predecessor
    for candidate in candidates:
        if type(candidate) is not PaperCircleCandidate:
            raise InvalidStandardPlacementInputError("Standard placement candidates must be typed paper circles.")
        spacing = _center_spacing(previous, candidate)
        if spacing <= _resolvability_floor(previous, candidate, tool_radius):
            raise StandardPlacementResolutionError("Ordered paper candidates contain a zero or binary64-unresolvable restart.")
        spacings.append(spacing)
        previous = candidate
    return tuple(spacings)


def _running_progress(spacings: Sequence[float]) -> float:
    # Left-to-right from 0.0, the order the per-step progress always used.
    total = 0.0
    for spacing in spacings:
        total += spacing
    return total


def _validate_ordered_progress(
    predecessor: PaperCircleCandidate,
    candidates: Sequence[PaperCircleCandidate],
    tool_radius: ToolRadius,
) -> tuple[float, ...]:
    progress = []
    total = 0.0
    for spacing in _validated_spacings(predecessor, candidates, tool_radius):
        total += spacing
        progress.append(total)
    return tuple(progress)


@dataclass(frozen=True)
class OrderedPlacementSpan:
    """One ordered candidate run between two placement stops, validated once.

    ``candidates[0]`` is the span's first predecessor and ``spacings[i]`` the
    centre spacing from ``candidates[i]`` to ``candidates[i + 1]``. Every later
    predecessor in the span sees a suffix of the same ordered pairs, so one
    validation here decides exactly what validating each step's remaining
    slice decided, in O(L) instead of O(L^2) per span.
    """

    candidates: tuple[PaperCircleCandidate, ...]
    spacings: tuple[float, ...]
    tool_radius: ToolRadius
    cap: EngagementCap

    @classmethod
    def build(
        cls,
        candidates: Sequence[PaperCircleCandidate],
        tool_radius: ToolRadius,
        cap: EngagementCap,
    ) -> Self:
        if not candidates or type(candidates[0]) is not PaperCircleCandidate or type(tool_radius) is not ToolRadius or type(cap) is not EngagementCap:
            raise InvalidStandardPlacementInputError("Standard placement requires a paper predecessor, tool radius, and cap.")
        spacings = _validated_spacings(candidates[0], candidates[1:], tool_radius)
        return cls(tuple(candidates), spacings, tool_radius, cap)


def select_next_in_span(span: OrderedPlacementSpan, cursor: int) -> StandardPlacement:
    """Select the next standard placement after ``span.candidates[cursor]``.

    Decides and reports exactly what ``select_next_standard_candidate`` does on
    ``span.candidates[cursor]`` and ``span.candidates[cursor + 1:]``.
    """
    if type(span) is not OrderedPlacementSpan or not 0 <= cursor < len(span.spacings):
        raise InvalidStandardPlacementInputError("Span placement requires a predecessor with at least one successor in its span.")
    return _bisect_next_candidate(
        span.candidates[cursor],
        span.candidates,
        cursor + 1,
        len(span.candidates) - cursor - 1,
        span.tool_radius,
        span.cap,
        lambda selected: _running_progress(span.spacings[cursor : cursor + selected + 1]),
    )


def select_next_standard_candidate(
    predecessor: PaperCircleCandidate,
    candidates: Sequence[PaperCircleCandidate],
    tool_radius: ToolRadius,
    cap: EngagementCap,
) -> StandardPlacement:
    """Bisect ordered candidates for Section 2.4's next standard placement."""
    if type(predecessor) is not PaperCircleCandidate or type(tool_radius) is not ToolRadius or type(cap) is not EngagementCap:
        raise InvalidStandardPlacementInputError("Standard placement requires a paper predecessor, tool radius, and cap.")
    progress = _validate_ordered_progress(predecessor, candidates, tool_radius)
    return _bisect_next_candidate(predecessor, candidates, 0, len(candidates), tool_radius, cap, progress.__getitem__)


def _bisect_next_candidate(
    predecessor: PaperCircleCandidate,
    candidates: Sequence[PaperCircleCandidate],
    first: int,
    count: int,
    tool_radius: ToolRadius,
    cap: EngagementCap,
    progress_at: Callable[[int], float],
) -> StandardPlacement:
    """Bisect ``candidates[first:first + count]``; indices below are relative to ``first``.

    The cap test and the paper's stopping band are exact decisions on the
    engagement cosine (``within`` / ``_BisectionStopBand``); only the
    reported diagnostics go through ``acos``.
    """
    band = _BisectionStopBand.build(cap)
    low = 0
    high = count - 1
    best: tuple[int, PredecessorEngagementCosine] | None = None
    while low <= high:
        middle = (low + high) // 2
        cosine = exact_predecessor_engagement(
            predecessor,
            candidates[first + middle],
            tool_radius,
        )
        if cosine.within(cap):
            best = middle, cosine
            if band.stops(cosine):
                break
            low = middle + 1
        else:
            high = middle - 1

    if best is None:
        forced = True
        first_margin = _overlap_margin(predecessor, candidates[first], tool_radius)
        if first_margin <= 0.0:
            raise StandardPlacementFragmentationError("The first successor exceeds the paper swept-disk overlap bound.")
        selected_index = 0
        cosine = exact_predecessor_engagement(
            predecessor,
            candidates[first + selected_index],
            tool_radius,
        )
    else:
        forced = False
        selected_index, cosine = best

    # Reporting only: no decision reads these angles.
    cap_angle = float(cap.theta)
    engagement = cosine.reported_angle()
    reached_terminal = selected_index == count - 1
    cap_shortfall = max(0.0, cap_angle - float(engagement))
    cap_excess = max(0.0, float(engagement) - cap_angle)
    candidate = candidates[first + selected_index]
    spacing = _center_spacing(predecessor, candidate)
    margin = _overlap_margin(predecessor, candidate, tool_radius)
    if margin <= 0.0:
        raise StandardPlacementFragmentationError("Selected successor violates the paper swept-disk overlap bound.")
    return StandardPlacement(
        candidate=candidate,
        candidate_advance=selected_index + 1,
        center_spacing=Millimetre(spacing),
        guide_progress=Millimetre(progress_at(selected_index)),
        maximum_engagement=engagement,
        cap_shortfall=Radian(cap_shortfall),
        cap_excess=Radian(cap_excess),
        overlap_margin=Millimetre(margin),
        reached_terminal=reached_terminal,
        forced=forced,
    )
