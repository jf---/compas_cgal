#include "held_disk_contour_2.h"

#include <cmath>
#include <optional>

namespace {
using FT = Epeck::FT;
using Coord = GpsPoint::CoordNT;

Coord point_distance_squared(const GpsPoint& point, const EPoint& center)
{
    // One point's coordinates share a root; subtracting rational center
    // coordinates never mixes extensions.
    const Coord dx = point.x() - Coord(center.x());
    const Coord dy = point.y() - Coord(center.y());
    return dx * dx + dy * dy;
}

Coord arc_distance_squared(const GpsXCurve& arc, const EPoint& center)
{
    const Coord start = point_distance_squared(arc.source(), center);
    const Coord end = point_distance_squared(arc.target(), center);
    Coord minimum = CGAL::compare(start, end) == CGAL::SMALLER ? start : end;
    if (!arc.is_circular()) {
        throw BrokenHeldContourBoundaryError("Filled disk union has a non-circular boundary.");
    }
    const ECircle support = arc.supporting_circle();
    const auto direction = center - support.center();
    const FT distance = direction.squared_length();
    if (CGAL::is_zero(distance)) return Coord(support.squared_radius());

    const FT root = support.squared_radius() / distance;
    const GpsPoint nearest(Coord(support.center().x(), direction.x(), root),
                            Coord(support.center().y(), direction.y(), root));
    const GpsTraits traits;
    if (arc.is_in_x_range(nearest)
        && traits.compare_y_at_x_2_object()(nearest, arc) == CGAL::EQUAL) {
        // (sqrt(distance) - sqrt(radius²))², represented in a single root.
        minimum = Coord(distance + support.squared_radius(), FT(-2),
                        distance * support.squared_radius());
    }
    return minimum;
}
} // namespace

std::pair<GpsPoint::CoordNT, bool> HeldDiskContour2::engagement_bound(
    const XY& candidate_center, double guide_radius, double cap_chord_ratio) const
{
    if (!std::isfinite(candidate_center[0]) || !std::isfinite(candidate_center[1])
        || !std::isfinite(guide_radius) || guide_radius <= 0
        || !std::isfinite(cap_chord_ratio) || cap_chord_ratio < 0 || cap_chord_ratio > 4) {
        throw InvalidHeldContourInputError("Contour bound requires finite XY, positive guide radius and chord ratio in [0,4].");
    }
    const EPoint center(candidate_center[0], candidate_center[1]);
    const FT guide(guide_radius);
    const Coord cap_cosine(FT(1) - FT(cap_chord_ratio) / FT(2));
    auto result = [&](const Coord& cosine) {
        return std::pair{cosine, CGAL::compare(cosine, cap_cosine) == CGAL::SMALLER};
    };
    if (contour_.oriented_side(GpsPoint(center.x(), center.y())) != CGAL::ON_POSITIVE_SIDE) {
        return result(Coord(FT(-1)));
    }
    std::optional<Coord> clearance_squared;
    const auto& arrangement = contour_.arrangement();
    for (auto edge = arrangement.edges_begin(); edge != arrangement.edges_end(); ++edge) {
        if (edge->face()->contained() == edge->twin()->face()->contained()) continue;
        const Coord distance = arc_distance_squared(edge->curve(), center);
        if (!clearance_squared || CGAL::compare(distance, *clearance_squared) == CGAL::SMALLER) {
            clearance_squared = distance;
        }
    }
    if (!clearance_squared) {
        throw BrokenHeldContourBoundaryError("Finite filled disk union has no exposed boundary.");
    }
    if (CGAL::compare(*clearance_squared, Coord(CGAL::square(guide + tool_radius_))) != CGAL::SMALLER) {
        return result(Coord(FT(1)));
    }
    if (CGAL::compare(*clearance_squared, Coord(CGAL::square(guide - tool_radius_))) != CGAL::LARGER) {
        return result(Coord(FT(-1)));
    }
    // Disk(center, clearance) is wholly cleared. On the advancing half-rim,
    // |p-center|² = guide² + tool² + 2*guide*tool*cos(alpha), alpha in [0,pi].
    // Uncut points must precede the disk crossing, so its angle bounds EVERY
    // connected engagement and their total, for every pose around the orbit.
    const Coord cosine = (*clearance_squared - Coord(CGAL::square(guide) + CGAL::square(tool_radius_)))
        / Coord(FT(2) * guide * tool_radius_);
    return result(cosine);
}
