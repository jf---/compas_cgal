#include "native_boundary_curve_2.h"
#include "reachable_boundary_sampling_2.h"

#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/determinant.h>
#include <CGAL/squared_distance_2.h>
#include <cmath>
#include <optional>
#include <utility>

// Decision discipline for this translation unit. The kernel is CORE::Expr:
// its floating-point filter decides every generic sign in nanoseconds, but a
// value that is identically zero can only be certified by refining to the
// root-separation bound, which grows exponentially with the number of
// independent square roots in the expression. Measured 2026-09-07: one such
// identity cost 1.5 ms to 2 s on a shallow scene and the Figure 5 query
// exceeded 337 s, while the same query without identities ran in 23 us.
// Therefore: (1) a quantity that holds by construction is recorded, never
// re-decided; (2) coincidences that would make a decision identically zero
// (shared support, parallel or collinear segments, a competitor tangent to
// the source tangent line at p) are detected first by rational or structural
// tests; (3) only genuine geometric ties, measure-zero for generic input,
// reach the exact path.

namespace {
using Point = ReachKernelPoint;
using Vector = ReachKernelVector;
using FT = ReachFT;
using Segment = ReachKernel::Segment_2;
using Circle = ReachKernel::Circle_2;
using Direction = ReachKernel::Direction_2;

enum class ContactKind { Focal, Vertex, Segment, Arc };
struct Winner {
    FT t;
    ContactKind kind;
    std::size_t index;  // meaningless for Focal
};

bool interior_foot(const Point& point, const Segment& segment)
{
    const FT projection = CGAL::scalar_product(point - segment.source(), segment.to_vector());
    return CGAL::sign(projection) == CGAL::POSITIVE
        && CGAL::compare(projection, segment.squared_length()) == CGAL::SMALLER;
}

// Angular membership of a radial direction on a trimmed arc. Every caller
// derives `radial` from a point that lies on the supporting circle by
// construction, so incidence is never re-decided here.
bool on_arc_direction(const Vector& radial, const ReachCurve& arc, bool strict)
{
    const auto circle = arc.supporting_circle();
    const Direction direction = radial.direction();
    const Direction first = (reachable_kernel_point(arc.source()) - circle.center()).direction();
    const Direction last = (reachable_kernel_point(arc.target()) - circle.center()).direction();
    if (direction == first || direction == last) return !strict;
    return arc.orientation() == CGAL::COUNTERCLOCKWISE
        ? direction.counterclockwise_in_between(first, last)
        : direction.counterclockwise_in_between(last, first);
}

bool same_supporting_line(const Point& a, const Point& b, const ReachCurve& other)
{
    return other.is_linear()
        && CGAL::collinear(a, b, reachable_kernel_point(other.source()))
        && CGAL::collinear(a, b, reachable_kernel_point(other.target()));
}

bool same_supporting_circle(const Circle& circle, const ReachCurve& other)
{
    if (other.is_linear()) return false;
    const auto candidate = other.supporting_circle();
    return candidate.center() == circle.center() && candidate.squared_radius() == circle.squared_radius();
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
    FT normal_squared;
    std::optional<Point> focal;
    std::optional<Circle> source_circle;
    std::optional<std::pair<Point, Point>> source_line;
    std::optional<Winner> first;
    const auto consider = [&](const FT& t, ContactKind kind, std::size_t index) {
        if (CGAL::sign(t) == CGAL::POSITIVE
            && (!first || CGAL::compare(t, first->t) == CGAL::SMALLER)) first = Winner{t, kind, index};
    };
    if (source_piece.curve.is_linear()) {
        const Point a = reachable_kernel_point(source_piece.curve.source());
        const Point b = reachable_kernel_point(source_piece.curve.target());
        normal = (b - a).perpendicular(cycle_.orientation);
        normal_squared = normal.squared_length();
        source_line = std::make_pair(a, b);
    } else {
        const Circle circle = source_piece.curve.supporting_circle();
        const bool convex = source_piece.curve.orientation() == cycle_.orientation;
        normal = convex ? circle.center() - p : p - circle.center();
        // p lies on the supporting circle by construction, so |normal|^2 is
        // the rational squared radius; the radical is never re-derived.
        normal_squared = circle.squared_radius();
        source_circle = circle;
        // At a convex arc's center its entire support arc ties the source.
        // The ordinary signed-circle equations are identities here, so the
        // positive focal event must be supplied explicitly (t=1 for this n).
        if (convex) { focal = circle.center(); consider(FT(1), ContactKind::Focal, 0); }
    }
    for (std::size_t i = 0; i < curves_.size(); ++i) {
        const auto& curve = curves_[i].curve();
        const Point vertex = reachable_kernel_point(curve.source());
        // A vertex on the source circle (every join vertex of an arc) ties the
        // focal event at exactly t = 1 for a convex source and is never a
        // positive candidate for a concave one: decided on rational data here,
        // not through the sample point's radicals.
        const bool vertex_on_source_circle = source_circle && source_circle->has_on_boundary(vertex);
        if (!vertex_on_source_circle) {
            const Vector delta = vertex - p;
            const FT point_denominator = FT(2) * CGAL::scalar_product(normal, delta);
            if (CGAL::sign(point_denominator) == CGAL::POSITIVE) {
                consider(delta.squared_length() / point_denominator, ContactKind::Vertex, i);
            }
        }
        if (curve.is_linear()) {
            // The source's own supporting line meets the normal only at p.
            if (source_line && same_supporting_line(source_line->first, source_line->second, curve)) continue;
            const Segment segment(vertex, reachable_kernel_point(curve.target()));
            const Vector other_normal = segment.to_vector().perpendicular(CGAL::COUNTERCLOCKWISE);
            const FT at_p = CGAL::scalar_product(other_normal, p - segment.source());
            const FT slope = CGAL::scalar_product(other_normal, normal);
            if (CGAL::is_zero(CGAL::determinant(normal.x(), normal.y(), other_normal.x(), other_normal.y()))) {
                // Parallel supports: |n||n'| equals |slope| exactly, so one
                // branch is identically empty and the other is rational.
                const FT t = -at_p / (FT(2) * slope);
                if (CGAL::sign(t) == CGAL::POSITIVE && interior_foot(p + t * normal, segment)) {
                    consider(t, ContactKind::Segment, i);
                }
                continue;
            }
            const FT norm_product = CGAL::sqrt(normal_squared * other_normal.squared_length());
            for (int sign : {-1, 1}) {
                const FT denominator = FT(sign) * norm_product - slope;  // nonzero: supports are not parallel
                const FT t = at_p / denominator;
                if (CGAL::sign(t) == CGAL::POSITIVE && interior_foot(p + t * normal, segment)) {
                    consider(t, ContactKind::Segment, i);
                }
            }
        } else {
            // On the source's own supporting circle the numerator and one
            // denominator vanish identically; it is handled as the focal event.
            if (source_circle && same_supporting_circle(*source_circle, curve)) continue;
            const Circle circle = curve.supporting_circle();
            const Vector delta_center = p - circle.center();
            const FT numerator = circle.squared_radius() - delta_center.squared_length();
            const FT projection = CGAL::scalar_product(delta_center, normal);
            const FT product_squared = circle.squared_radius() * normal_squared;  // (R_i |n|)^2, no radical
            // A competitor tangent to the source tangent line at p makes one
            // branch's denominator vanish identically; that branch is empty.
            const bool tangent_at_p = CGAL::compare(projection * projection, product_squared) == CGAL::EQUAL;
            const int empty_sign = tangent_at_p ? static_cast<int>(CGAL::sign(projection)) : 0;
            const FT root = CGAL::sqrt(product_squared);  // the only radical this competitor contributes
            for (int sign : {-1, 1}) {
                if (sign == empty_sign) continue;
                const FT t = numerator / (FT(2) * (projection - FT(sign) * root));
                if (CGAL::sign(t) != CGAL::POSITIVE) continue;
                // |m-o| = R_i + sign t |n| must be nonnegative. For sign = -1
                // that is t^2 |n|^2 <= R_i^2, decided without a new radical.
                // The spurious swallowing branch |m-o| = t|n| - R_i is not a contact.
                if (sign < 0) {
                    const auto reach = CGAL::compare(t * t * normal_squared, circle.squared_radius());
                    if (reach == CGAL::LARGER) continue;
                    if (reach == CGAL::EQUAL) { consider(t, ContactKind::Arc, i); continue; }  // m at the competitor center
                }
                const Vector radial = (p + t * normal) - circle.center();
                if (on_arc_direction(radial, curve, true)) consider(t, ContactKind::Arc, i);
            }
        }
    }
    if (!first) throw NativeBoundaryMedialConstructionError("Boundary normal has no positive medial contact.");
    const Point m = first->kind == ContactKind::Focal ? *focal : p + first->t * normal;
    const FT squared_clearance = first->t * first->t * normal_squared;
    if (design_.oriented_side(ReachPoint(m.x(), m.y())) != CGAL::ON_POSITIVE_SIDE) {
        throw NativeBoundaryMedialConstructionError("First medial candidate is not inside its exact design.");
    }
    std::vector<std::size_t> vertices, segments, arcs;
    for (std::size_t i = 0; i < curves_.size(); ++i) {
        const auto& curve = curves_[i].curve();
        const bool is_winner = first->kind != ContactKind::Focal && i == first->index;
        const Point start = reachable_kernel_point(curve.source());
        // Every boundary vertex is the start of exactly one curve, so this
        // test sees each vertex once; curve interiors are judged separately
        // below and never fall back to an endpoint distance, which for the
        // winning vertex would re-decide the tie the disk was built on.
        if (source_circle && source_circle->has_on_boundary(start)) {
            // On the source circle the tie with the focal disk holds by
            // construction, and no other winner can reach such a vertex.
            if (first->kind == ContactKind::Focal && start != p) vertices.push_back(i);
        } else if (is_winner && first->kind == ContactKind::Vertex) {
            vertices.push_back(i);  // holds by construction
        } else if (start != p) {
            const auto vertex_relation = CGAL::compare(CGAL::squared_distance(m, start), squared_clearance);
            if (vertex_relation == CGAL::SMALLER) {
                throw NativeBoundaryMedialConstructionError("Medial disk contains a boundary vertex; one-sided normal may not be admissible.");
            }
            if (vertex_relation == CGAL::EQUAL) vertices.push_back(i);
        }
        if (is_winner) {
            // The winning contact holds by construction: recorded, not re-decided.
            if (first->kind == ContactKind::Segment) segments.push_back(i);
            if (first->kind == ContactKind::Arc) arcs.push_back(i);
            continue;
        }
        const bool source_support = (source_line && same_supporting_line(source_line->first, source_line->second, curve))
            || (source_circle && same_supporting_circle(*source_circle, curve));
        if (source_support) {
            // The disk is tangent to the source support at p, so that support
            // never enters the disk; it is a distinct contact only at the
            // focal event, where the whole support ties.
            if (first->kind == ContactKind::Focal) arcs.push_back(i);
            continue;
        }
        // Interior of the curve only. If the nearest point is an endpoint,
        // the vertex test above has already judged it.
        FT interior_distance_squared;
        if (curve.is_linear()) {
            const Segment segment(start, reachable_kernel_point(curve.target()));
            if (!interior_foot(m, segment)) continue;
            interior_distance_squared = CGAL::squared_distance(m, segment.supporting_line());
        } else {
            const Circle circle = curve.supporting_circle();
            const Vector radial = m - circle.center();
            const FT distance_squared = radial.squared_length();
            if (CGAL::is_zero(distance_squared)) {
                interior_distance_squared = circle.squared_radius();
            } else {
                if (!on_arc_direction(radial, curve, true)) continue;
                // (|radial| - R)^2 with a single radical node.
                interior_distance_squared = distance_squared + circle.squared_radius()
                    - FT(2) * CGAL::sqrt(distance_squared * circle.squared_radius());
            }
        }
        const auto relation = CGAL::compare(interior_distance_squared, squared_clearance);
        if (relation == CGAL::SMALLER) {
            throw NativeBoundaryMedialConstructionError("Medial disk crosses a closer boundary feature; one-sided normal may not be admissible.");
        }
        if (relation != CGAL::EQUAL) continue;  // a genuine tie is the only way past this point
        (curve.is_linear() ? segments : arcs).push_back(i);
    }
    if (vertices.empty() && segments.empty() && arcs.empty()) {
        throw NativeBoundaryMedialConstructionError("Medial contact has no distinct competing boundary contact.");
    }
    const FT tool(tool_radius);
    if (CGAL::compare(squared_clearance, tool * tool) == CGAL::SMALLER) {
        throw NoPositiveNativeBoundaryCircleError("Medial clearance is smaller than the tool radius.");
    }
    // Outputs only from here on: one radical for the reported lengths.
    const FT normal_length = CGAL::sqrt(normal_squared);
    const FT clearance = first->t * normal_length;
    const Point q = p + (tool / normal_length) * normal;
    const Point center = CGAL::midpoint(q, m);
    const FT guide = (clearance - tool) / FT(2);
    // |m-p| = clearance, |q-p| = tool and |c-q| = |c-m| = guide hold by
    // construction; the medial tests witness them on reporting values.
    return boundary_normal::BoundaryNormalCircleProposal2(p, m, q, center, guide, clearance,
        std::move(vertices), std::move(segments), std::move(arcs));
}
