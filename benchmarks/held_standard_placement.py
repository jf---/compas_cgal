"""Pure predecessor-only placement for the Held Figure 5(a) standard path.

The paper's standard model evaluates a candidate machining circle only against
the immediately preceding swept disk.  This module implements that local model
on already transformed ``(c, rho, q)`` circles.  It neither owns nor consults a
global depleted-stock contour.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Self
from typing import Sequence

from compas.tolerance import TOL

from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Radian
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY

# Held and Pfeiffer, Section 2.4: stop once theta_max - 0.001 <= theta <= theta_max.
PAPER_BISECTION_ANGLE_TOLERANCE = Radian(0.001)


class InvalidPaperCircleCandidateError(ValueError):
    """A transformed ``(c, rho, q)`` tuple does not define one circle."""


class InvalidStandardPlacementInputError(ValueError):
    """The predecessor placement model received invalid typed input."""


class StandardPlacementResolutionError(RuntimeError):
    """The ordered candidates cannot resolve the paper's next placement."""


class StandardPlacementFragmentationError(RuntimeError):
    """The next candidate violates the paper's swept-disk overlap bound."""


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


def maximum_predecessor_engagement(
    predecessor: PaperCircleCandidate,
    candidate: PaperCircleCandidate,
    tool_radius: ToolRadius,
) -> Radian:
    """Compute Section 2.3's standard-model maximum engagement.

    The reported angle follows the paper's two geometric branches.  All circle
    geometry is taken from the transformed candidate; no baseline radius or
    stock contour is substituted.
    """
    if type(predecessor) is not PaperCircleCandidate or type(candidate) is not PaperCircleCandidate or type(tool_radius) is not ToolRadius:
        raise InvalidStandardPlacementInputError("Predecessor engagement requires two paper circles and a tool radius.")
    distance = _center_spacing(predecessor, candidate)
    tool = float(tool_radius.value)
    previous_radius = _radius(predecessor)
    radius = _radius(candidate)
    if distance == 0.0:
        if radius <= previous_radius:
            return Radian(0.0)
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

    if distance > previous_outer_radius + current_outer_radius or distance < abs(previous_outer_radius - current_outer_radius):
        raise StandardPlacementFragmentationError("Paper overlap correction requires intersecting predecessor/current swept disks.")
    w_x = (previous_outer_radius * previous_outer_radius - current_outer_radius * current_outer_radius - distance * distance) / (2.0 * distance)
    w_y_squared = current_outer_radius * current_outer_radius - w_x * w_x
    if w_y_squared < 0.0:
        raise StandardPlacementFragmentationError("Paper overlap correction has no real swept-disk intersection.")
    q_scale = radius / current_outer_radius
    q_corrected_x = q_scale * w_x
    q_corrected_y_squared = q_scale * q_scale * w_y_squared
    q_to_previous_squared = (q_corrected_x + distance) * (q_corrected_x + distance) + q_corrected_y_squared
    q_to_previous = math.sqrt(q_to_previous_squared)
    chord_offset_ratio = (q_to_previous_squared + tool * tool - previous_outer_radius * previous_outer_radius) / (2.0 * q_to_previous * tool)
    return Radian(2.0 * math.acos(abs(_clamped_unit(chord_offset_ratio))))


def _validate_ordered_progress(
    predecessor: PaperCircleCandidate,
    candidates: Sequence[PaperCircleCandidate],
    tool_radius: ToolRadius,
) -> tuple[float, ...]:
    if not candidates:
        raise InvalidStandardPlacementInputError("Standard placement requires at least one successor candidate.")
    progress = []
    previous = predecessor
    total = 0.0
    for candidate in candidates:
        if type(candidate) is not PaperCircleCandidate:
            raise InvalidStandardPlacementInputError("Standard placement candidates must be typed paper circles.")
        spacing = _center_spacing(previous, candidate)
        if spacing <= _resolvability_floor(previous, candidate, tool_radius):
            raise StandardPlacementResolutionError("Ordered paper candidates contain a zero or binary64-unresolvable restart.")
        total += spacing
        progress.append(total)
        previous = candidate
    return tuple(progress)


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
    cap_angle = float(cap.theta)
    low = 0
    high = len(candidates) - 1
    best: tuple[int, Radian] | None = None
    while low <= high:
        middle = (low + high) // 2
        engagement = maximum_predecessor_engagement(
            predecessor,
            candidates[middle],
            tool_radius,
        )
        if float(engagement) <= cap_angle:
            best = middle, engagement
            if cap_angle - float(engagement) <= float(PAPER_BISECTION_ANGLE_TOLERANCE):
                break
            low = middle + 1
        else:
            high = middle - 1

    if best is None:
        forced = True
        first_margin = _overlap_margin(predecessor, candidates[0], tool_radius)
        if first_margin <= 0.0:
            raise StandardPlacementFragmentationError("The first successor exceeds the paper swept-disk overlap bound.")
        selected_index = 0
        engagement = maximum_predecessor_engagement(
            predecessor,
            candidates[selected_index],
            tool_radius,
        )
    else:
        forced = False
        selected_index, engagement = best

    reached_terminal = selected_index == len(candidates) - 1
    cap_shortfall = max(0.0, cap_angle - float(engagement))
    cap_excess = max(0.0, float(engagement) - cap_angle)
    candidate = candidates[selected_index]
    spacing = _center_spacing(predecessor, candidate)
    margin = _overlap_margin(predecessor, candidate, tool_radius)
    if margin <= 0.0:
        raise StandardPlacementFragmentationError("Selected successor violates the paper swept-disk overlap bound.")
    return StandardPlacement(
        candidate=candidate,
        candidate_advance=selected_index + 1,
        center_spacing=Millimetre(spacing),
        guide_progress=Millimetre(progress[selected_index]),
        maximum_engagement=engagement,
        cap_shortfall=Radian(cap_shortfall),
        cap_excess=Radian(cap_excess),
        overlap_margin=Millimetre(margin),
        reached_terminal=reached_terminal,
        forced=forced,
    )
