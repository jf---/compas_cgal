#include "held_disk_contour_2.h"

#include <cmath>
#include <optional>

namespace {
using FT = Epeck::FT;
using Coord = GpsPoint::CoordNT;

EPoint checked_point(const HeldDiskContour2::XY& xy)
{
    if (!std::isfinite(xy[0]) || !std::isfinite(xy[1])) {
        throw InvalidHeldContourInputError("Held contour coordinates must be finite.");
    }
    return EPoint(xy[0], xy[1]);
}

FT checked_radius(double radius)
{
    if (!std::isfinite(radius) || radius <= 0) {
        throw InvalidHeldContourInputError("Held contour radii must be finite and positive.");
    }
    return FT(radius);
}

// On one circle, CCW order needs only exact coordinate comparisons. In
// particular, do not form cross products between different root extensions.
bool ccw_less(const GpsPoint& a, const GpsPoint& b, const EPoint& center)
{
    auto upper = [&](const GpsPoint& p) {
        const auto y = CGAL::compare(p.y(), Coord(center.y()));
        return y == CGAL::LARGER
            || (y == CGAL::EQUAL
                && CGAL::compare(p.x(), Coord(center.x())) == CGAL::LARGER);
    };
    const bool a_upper = upper(a), b_upper = upper(b);
    if (a_upper != b_upper) return a_upper;
    const auto x = CGAL::compare(a.x(), b.x());
    return a_upper ? x == CGAL::LARGER : x == CGAL::SMALLER;
}

bool on_ccw_arc(const GpsPoint& point, const GpsPoint& start,
                 const GpsPoint& end, const EPoint& center)
{
    const bool after_start = !ccw_less(point, start, center);
    const bool before_end = !ccw_less(end, point, center);
    return ccw_less(start, end, center)
        ? after_start && before_end : after_start || before_end;
}
} // namespace

HeldDiskContour2::HeldDiskContour2(const XY& center, double guide_radius,
                                  double tool_radius)
    : HeldDiskContour2(checked_point(center), checked_radius(guide_radius),
                       checked_radius(tool_radius))
{
}

HeldDiskContour2::HeldDiskContour2(const EPoint& center, const FT& guide_radius,
                                  const FT& tool_radius)
    : tool_radius_(tool_radius),
      predecessor_(center, CGAL::square(guide_radius + tool_radius))
{
    contour_.join(disk_polygon(center, guide_radius + tool_radius));
}

void HeldDiskContour2::append(const XY& center, double guide_radius)
{
    const EPoint point = checked_point(center);
    const FT outer = checked_radius(guide_radius) + tool_radius_;
    const ECircle next(point, CGAL::square(outer));
    contour_.join(disk_polygon(point, outer));
    predecessor_ = next;
}

std::pair<GpsPoint, bool> HeldDiskContour2::contact_toward(
    const XY& candidate_center) const
{
    const EPoint candidate = checked_point(candidate_center);
    const EPoint center = predecessor_.center();
    const auto direction = candidate - center;
    const FT distance_squared = direction.squared_length();
    if (CGAL::is_zero(distance_squared)) {
        throw UndefinedPredecessorDirectionError("Concentric candidate has no unique critical-point direction.");
    }
    const FT root = predecessor_.squared_radius() / distance_squared;
    const GpsPoint b(Coord(center.x(), direction.x(), root),
                     Coord(center.y(), direction.y(), root));

    std::optional<GpsPoint> before_b, last;
    const auto& arrangement = contour_.arrangement();
    for (auto edge = arrangement.edges_begin(); edge != arrangement.edges_end(); ++edge) {
        if (edge->face()->contained() == edge->twin()->face()->contained()) continue;
        const GpsXCurve& arc = edge->curve();
        if (!arc.is_circular()) continue;
        const ECircle support = arc.supporting_circle();
        if (support.center() != center
            || support.squared_radius() != predecessor_.squared_radius()) continue;

        GpsPoint start = arc.source(), end = arc.target();
        if (arc.orientation() == CGAL::CLOCKWISE) std::swap(start, end);
        if (on_ccw_arc(b, start, end, center)) return {b, false};

        // Enter an exposed arc at its CCW end while travelling clockwise.
        // If no endpoint precedes b, cross the +x seam to the last endpoint.
        if (!last || ccw_less(*last, end, center)) last = end;
        if (!ccw_less(b, end, center)
            && (!before_b || ccw_less(*before_b, end, center))) before_b = end;
    }
    if (before_b) return {*before_b, true};
    if (last) return {*last, true};
    throw NoExposedPredecessorArcError("Latest outer disk has no exposed contour arc.");
}
