#include "boundary_circle_engagement_2.h"

#include <CGAL/squared_distance_2.h>

#include <algorithm>
#include <cmath>

namespace boundary_normal {
namespace {

std::pair<double, bool> report(const FT& cosine, const FT& cap_cosine)
{
    if (CGAL::compare(cosine, FT(-1)) == CGAL::SMALLER
        || CGAL::compare(cosine, FT(1)) == CGAL::LARGER) {
        throw BoundaryCircleEngagementGeometryError(
            "Exact predecessor engagement cosine is outside [-1, 1].");
    }
    const bool exceeds = CGAL::compare(cosine, cap_cosine) == CGAL::SMALLER;
    // Clamp only the final approximate reporting conversion, never a decision.
    return {std::acos(std::clamp(CGAL::to_double(cosine), -1.0, 1.0)), exceeds};
}

FT corrected_cosine(const FT& previous_outer_squared,
                    const FT& outer, const FT& radius,
                    const FT& tool, const FT& distance_squared)
{
    const FT numerator = previous_outer_squared - outer * outer - distance_squared;
    if (CGAL::sign(4 * distance_squared * outer * outer - numerator * numerator)
        == CGAL::NEGATIVE) {
        throw BoundaryCircleEngagementGeometryError(
            "Corrected predecessor outer circles have no real intersection.");
    }
    // Fig. 4(d): radially project the outer-circle intersection onto M_i.
    // Expansion cancels the intersection height and its square root, exactly
    // as in corrected_engagement_cosine_squared in the older reporting API.
    const FT q_squared = radius * radius + (radius / outer) * numerator
        + distance_squared;
    if (CGAL::sign(q_squared) != CGAL::POSITIVE) {
        throw BoundaryCircleEngagementGeometryError(
            "Corrected predecessor contact has no chord direction.");
    }
    const FT chord_numerator = q_squared + tool * tool - previous_outer_squared;
    const FT half_angle_cosine_squared = chord_numerator * chord_numerator
        / (4 * tool * tool * q_squared);
    return 2 * half_angle_cosine_squared - 1;
}

} // namespace

std::pair<double, bool> boundary_circle_engagement(
    const BoundaryNormalCircleProposal2& previous,
    const BoundaryNormalCircleProposal2& current,
    double cap_chord_ratio)
{
    if (!std::isfinite(cap_chord_ratio)
        || CGAL::compare(FT(cap_chord_ratio), FT(0)) == CGAL::SMALLER
        || CGAL::compare(FT(cap_chord_ratio), FT(4)) == CGAL::LARGER) {
        throw InvalidBoundaryEngagementCapError(
            "Engagement cap requires a finite squared-chord ratio in [0, 4].");
    }
    const FT tool = previous.exact_tool_radius();
    if (CGAL::sign(tool) != CGAL::POSITIVE
        || CGAL::compare(tool, current.exact_tool_radius()) != CGAL::EQUAL) {
        throw BoundaryCircleToolMismatchError(
            "Predecessor and successor require the same positive exact tool radius.");
    }
    const FT cap_cosine = 1 - FT(cap_chord_ratio) / 2;
    const FT previous_radius = previous.exact_guide_radius();
    const FT radius = current.exact_guide_radius();
    const FT previous_outer = previous_radius + tool;
    const FT outer = radius + tool;
    const FT distance_squared = CGAL::squared_distance(
        previous.exact_center(), current.exact_center());
    const FT distance = CGAL::sqrt(distance_squared);

    if (CGAL::compare(distance + outer, previous_outer) != CGAL::LARGER) {
        return report(FT(1), cap_cosine);
    }
    if (CGAL::is_zero(radius)) {
        throw UncoveredStationaryCircleError(
            "An uncovered stationary successor needs explicit entry or transition semantics.");
    }
    // Eq. 4 keeps the prior disk and successor sweep connected without a hole.
    const auto spacing = CGAL::compare(distance + radius - previous_radius, 2 * tool);
    if (spacing == CGAL::LARGER) {
        throw BoundaryCircleSpacingError(
            "Successor violates Held Eq. 4 swept-disk continuity.");
    }
    if (spacing == CGAL::EQUAL) return report(FT(-1), cap_cosine);

    const FT previous_outer_squared = previous_outer * previous_outer;
    if (CGAL::is_zero(distance_squared)) {
        return report((previous_outer_squared - tool * tool - radius * radius)
                      / (2 * tool * radius), cap_cosine);
    }

    const FT b = previous_outer - distance;
    // These branches retain the existing conservative full-slot bound. They
    // are not advertised as an exact maximum outside the Fig. 4 domain.
    if (CGAL::sign(b) != CGAL::POSITIVE) return report(FT(-1), cap_cosine);
    const FT q_x = (b * b - tool * tool + radius * radius) / (2 * b);
    if (CGAL::sign(radius * radius - q_x * q_x) == CGAL::NEGATIVE) {
        return report(FT(-1), cap_cosine);
    }
    const FT w_squared = outer * outer
        + 2 * distance * (outer / radius) * q_x + distance_squared;
    if (CGAL::compare(w_squared, previous_outer_squared) == CGAL::LARGER) {
        return report((b * b - tool * tool - radius * radius)
                      / (2 * tool * radius), cap_cosine);
    }
    return report(corrected_cosine(previous_outer_squared, outer, radius,
                                   tool, distance_squared), cap_cosine);
}

} // namespace boundary_normal
