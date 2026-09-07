#include "native_boundary_curve_2.h"
#include "reachable_boundary_sampling_2.h"

#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/squared_distance_2.h>
#include <cmath>
#include <optional>

namespace {
using Point = ReachKernelPoint;
using Vector = ReachKernelVector;
using FT = ReachFT;
using Segment = ReachKernel::Segment_2;

bool interior_foot(const Point& point, const Segment& segment)
{
    const FT projection = CGAL::scalar_product(point - segment.source(), segment.to_vector());
    return CGAL::sign(projection) == CGAL::POSITIVE
        && CGAL::compare(projection, segment.squared_length()) == CGAL::SMALLER;
}

bool on_arc(const Point& point, const ReachCurve& arc, bool strict)
{
    const Point start = reachable_kernel_point(arc.source()), end = reachable_kernel_point(arc.target());
    if (point == start || point == end) return !strict;
    const auto circle = arc.supporting_circle();
    if (!circle.has_on_boundary(point)) return false;
    const auto direction = (point - circle.center()).direction();
    const auto first = (start - circle.center()).direction();
    const auto last = (end - circle.center()).direction();
    return arc.orientation() == CGAL::COUNTERCLOCKWISE
        ? direction.counterclockwise_in_between(first, last)
        : direction.counterclockwise_in_between(last, first);
}

FT boundary_distance_squared(const Point& point, const ReachCurve& curve)
{
    const Point start = reachable_kernel_point(curve.source()), end = reachable_kernel_point(curve.target());
    if (curve.is_linear()) return CGAL::squared_distance(point, Segment(start, end));
    const auto circle = curve.supporting_circle();
    const Vector radial = point - circle.center();
    const FT distance_squared = radial.squared_length();
    if (CGAL::is_zero(distance_squared)) return circle.squared_radius();
    const Point foot = circle.center() + radial * CGAL::sqrt(circle.squared_radius() / distance_squared);
    if (on_arc(foot, curve, false)) {
        const FT delta = CGAL::sqrt(distance_squared) - CGAL::sqrt(circle.squared_radius());
        return delta * delta;
    }
    const FT a = CGAL::squared_distance(point, start), b = CGAL::squared_distance(point, end);
    return CGAL::compare(a, b) == CGAL::SMALLER ? a : b;
}
}

boundary_normal::BoundaryNormalCircleProposal2 NativeBoundary2::circle_on_piece(
    std::int64_t piece_index, double parameter, double tool_radius) const
{
    if (piece_index < 0 || static_cast<std::size_t>(piece_index) >= cycle_.curves.size()
        || !std::isfinite(parameter) || !std::isfinite(tool_radius)
        || parameter < 0 || parameter > 1 || tool_radius <= 0) {
        throw InvalidNativeBoundaryMedialInputError("Query requires a valid piece, parameter in [0,1], and positive finite tool radius.");
    }
    const auto& source_piece = cycle_.curves[static_cast<std::size_t>(piece_index)];
    const Point p = reachable_kernel_point(sample_reachable_boundary(source_piece, parameter));
    Vector normal;
    FT normal_length;
    std::optional<FT> first;
    std::optional<Point> focal;
    const auto consider = [&](const FT& t) {
        if (CGAL::sign(t) == CGAL::POSITIVE
            && (!first || CGAL::compare(t, *first) == CGAL::SMALLER)) first = t;
    };
    if (source_piece.curve.is_linear()) {
        normal = (reachable_kernel_point(source_piece.curve.target())
                - reachable_kernel_point(source_piece.curve.source())).perpendicular(cycle_.orientation);
        normal_length = CGAL::sqrt(normal.squared_length());
    } else {
        const auto source_circle = source_piece.curve.supporting_circle();
        const bool convex = source_piece.curve.orientation() == cycle_.orientation;
        normal = convex ? source_circle.center() - p : p - source_circle.center();
        normal_length = CGAL::sqrt(source_circle.squared_radius());
        // At a convex arc's center its entire support arc ties the source.
        // The ordinary signed-circle equations are identities here, so the
        // positive focal event must be supplied explicitly (t=1 for this n).
        if (convex) { focal = source_circle.center(); consider(FT(1)); }
    }
    const FT normal_squared = normal_length * normal_length;
    for (const auto& feature : curves_) {
        const auto& curve = feature.curve();
        const Point vertex = reachable_kernel_point(curve.source());
        const Vector delta = vertex - p;
        const FT point_denominator = FT(2) * CGAL::scalar_product(normal, delta);
        if (CGAL::sign(point_denominator) == CGAL::POSITIVE) {
            consider(delta.squared_length() / point_denominator);
        }
        if (curve.is_linear()) {
            const Segment segment(vertex, reachable_kernel_point(curve.target()));
            const Vector other_normal = segment.to_vector().perpendicular(CGAL::COUNTERCLOCKWISE);
            const FT at_p = CGAL::scalar_product(other_normal, p - segment.source());
            const FT slope = CGAL::scalar_product(other_normal, normal);
            const FT norm_product = CGAL::sqrt(normal_squared * other_normal.squared_length());
            for (int sign : {-1, 1}) {
                const FT denominator = FT(sign) * norm_product - slope;
                if (CGAL::is_zero(denominator)) continue;
                const FT t = at_p / denominator;
                if (CGAL::sign(t) == CGAL::POSITIVE && interior_foot(p + t * normal, segment)) consider(t);
            }
        } else {
            const auto circle = curve.supporting_circle();
            const FT radius = CGAL::sqrt(circle.squared_radius());
            const Vector delta_center = p - circle.center();
            const FT numerator = circle.squared_radius() - delta_center.squared_length();
            for (int sign : {-1, 1}) {
                const FT denominator = FT(2) * (CGAL::scalar_product(delta_center, normal)
                    - FT(sign) * radius * normal_length);
                if (CGAL::is_zero(denominator)) continue;
                const FT t = numerator / denominator;
                if (CGAL::sign(t) != CGAL::POSITIVE) continue;
                const FT distance = radius + FT(sign) * t * normal_length;
                // Keep the unsquared branch: |m-o|=R±t|n| must be nonnegative.
                // The spurious swallowing branch |m-o|=t|n|-R is not a contact.
                if (CGAL::sign(distance) == CGAL::NEGATIVE) continue;
                const Point candidate = p + t * normal;
                if (CGAL::compare(CGAL::squared_distance(candidate, circle.center()), distance * distance) != CGAL::EQUAL) continue;
                if (CGAL::is_zero(distance)) consider(t);
                else {
                    const Point foot = circle.center() + (radius / distance) * (candidate - circle.center());
                    if (on_arc(foot, curve, true)) consider(t);
                }
            }
        }
    }
    if (!first) throw NativeBoundaryMedialConstructionError("Boundary normal has no positive medial contact.");
    const Point m = focal && CGAL::compare(*first, FT(1)) == CGAL::EQUAL
        ? *focal : p + *first * normal;
    const FT clearance = *first * normal_length, squared_clearance = clearance * clearance;
    if (design_.oriented_side(ReachPoint(m.x(), m.y())) != CGAL::ON_POSITIVE_SIDE) {
        throw NativeBoundaryMedialConstructionError("First medial candidate is not inside its exact design.");
    }
    std::vector<std::size_t> vertices, segments, arcs;
    for (std::size_t i = 0; i < curves_.size(); ++i) {
        const auto& curve = curves_[i].curve();
        const Point start = reachable_kernel_point(curve.source()), end = reachable_kernel_point(curve.target());
        const auto relation = CGAL::compare(boundary_distance_squared(m, curve), squared_clearance);
        if (relation == CGAL::SMALLER) {
            throw NativeBoundaryMedialConstructionError("Medial disk crosses a closer boundary feature; one-sided normal may not be admissible.");
        }
        if (start != p && CGAL::compare(CGAL::squared_distance(m, start), squared_clearance) == CGAL::EQUAL) vertices.push_back(i);
        if (relation != CGAL::EQUAL) continue;
        if (curve.is_linear()) {
            const Segment segment(start, end);
            if (interior_foot(m, segment)) {
                const Point foot = segment.supporting_line().projection(m);
                if (foot != p) segments.push_back(i);
            }
        } else {
            const auto circle = curve.supporting_circle();
            const Vector radial = m - circle.center();
            if (radial == Vector(CGAL::NULL_VECTOR)) arcs.push_back(i);
            else {
                const Point foot = circle.center() + radial * CGAL::sqrt(circle.squared_radius() / radial.squared_length());
                if (foot != p && on_arc(foot, curve, true)) arcs.push_back(i);
            }
        }
    }
    if (vertices.empty() && segments.empty() && arcs.empty()) {
        throw NativeBoundaryMedialConstructionError("Medial contact has no distinct competing boundary contact.");
    }
    const FT tool(tool_radius);
    if (CGAL::compare(clearance, tool) == CGAL::SMALLER) {
        throw NoPositiveNativeBoundaryCircleError("Medial clearance is smaller than the tool radius.");
    }
    const Point q = p + (tool / normal_length) * normal;
    const Point center = CGAL::midpoint(q, m);
    const FT guide = (clearance - tool) / FT(2);
    if (CGAL::compare(CGAL::squared_distance(m, p), squared_clearance) != CGAL::EQUAL
        || CGAL::compare(CGAL::squared_distance(q, p), tool * tool) != CGAL::EQUAL
        || CGAL::compare(CGAL::squared_distance(center, q), guide * guide) != CGAL::EQUAL
        || CGAL::compare(CGAL::squared_distance(center, m), guide * guide) != CGAL::EQUAL) {
        throw NativeBoundaryMedialConstructionError("Coupled boundary-circle diameter invariant failed.");
    }
    return boundary_normal::BoundaryNormalCircleProposal2(p, m, q, center, guide, clearance,
        std::move(vertices), std::move(segments), std::move(arcs));
}
