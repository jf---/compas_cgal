#include "boundary_normal_circle_2.h"

#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_2.h>

#include <cmath>
#include <optional>
#include <utility>

namespace boundary_normal {
namespace {

XY report(const Point& point)
{
    return {CGAL::to_double(point.x()), CGAL::to_double(point.y())};
}

bool interior_foot(const Point& point, const Kernel::Segment_2& segment)
{
    const auto direction = segment.to_vector();
    const FT projection = CGAL::scalar_product(point - segment.source(), direction);
    return CGAL::sign(projection) == CGAL::POSITIVE
        && CGAL::compare(projection, direction.squared_length()) == CGAL::SMALLER;
}

} // namespace

BoundaryNormalCircle2::BoundaryNormalCircle2(const std::vector<XY>& boundary)
{
    if (boundary.size() < 3) {
        throw InvalidBoundaryPolygonError("A boundary polygon requires at least three vertices.");
    }
    for (const XY& xy : boundary) {
        if (!std::isfinite(xy[0]) || !std::isfinite(xy[1])) {
            throw InvalidBoundaryPolygonError("Boundary coordinates must be finite world-XY millimetres.");
        }
        polygon_.push_back(Point(xy[0], xy[1]));
    }
    for (std::size_t index = 0; index < polygon_.size(); ++index) {
        if (polygon_[index] == polygon_[(index + 1) % polygon_.size()]) {
            throw InvalidBoundaryPolygonError("Boundary contains a zero-length segment.");
        }
    }
    if (!polygon_.is_simple() || polygon_.orientation() == CGAL::COLLINEAR) {
        throw InvalidBoundaryPolygonError("Boundary must be a simple polygon with nonzero area.");
    }
}

BoundaryNormalCircleProposal2 BoundaryNormalCircle2::query(
    std::int64_t source_segment, double parameter, double tool_radius) const
{
    if (source_segment < 0 || static_cast<std::size_t>(source_segment) >= polygon_.size()
        || !std::isfinite(parameter) || !std::isfinite(tool_radius)) {
        throw InvalidBoundaryNormalInputError("Query requires a valid segment index and finite parameter/radius.");
    }
    const FT u(parameter), tool(tool_radius);
    if (CGAL::sign(u) == CGAL::NEGATIVE || CGAL::compare(u, FT(1)) == CGAL::LARGER
        || CGAL::sign(tool) != CGAL::POSITIVE) {
        throw InvalidBoundaryNormalInputError("Parameter must lie in [0,1] and cutter radius must be positive.");
    }
    if (CGAL::is_zero(u) || CGAL::compare(u, FT(1)) == CGAL::EQUAL) {
        throw BoundaryVertexQueryUnsupportedError("A vertex requires a separate normal-cone query; use a strict segment interior.");
    }
    const auto source = static_cast<std::size_t>(source_segment);
    const Point a = polygon_[source], b = polygon_[(source + 1) % polygon_.size()];
    const Point p = CGAL::barycenter(a, FT(1) - u, b, u);
    const auto normal = (b - a).perpendicular(polygon_.orientation());
    return construct(p, normal, tool, CGAL::sqrt(normal.squared_length()));
}

BoundaryNormalCircleProposal2 BoundaryNormalCircle2::query_vertex(
    std::int64_t vertex, const XY& inward_direction, double tool_radius) const
{
    if (vertex < 0 || static_cast<std::size_t>(vertex) >= polygon_.size()
        || !std::isfinite(tool_radius) || CGAL::sign(FT(tool_radius)) != CGAL::POSITIVE) {
        throw InvalidBoundaryNormalInputError("Vertex query requires a valid index and positive finite cutter radius.");
    }
    if (!std::isfinite(inward_direction[0]) || !std::isfinite(inward_direction[1])) {
        throw InvalidBoundaryVertexDirectionError("Vertex direction must be finite.");
    }
    const auto index = static_cast<std::size_t>(vertex);
    const Point p = polygon_[index];
    const Point previous = polygon_[(index + polygon_.size() - 1) % polygon_.size()];
    const Point next = polygon_[(index + 1) % polygon_.size()];
    const auto turn = CGAL::orientation(previous, p, next);
    const Kernel::Vector_2 direction(inward_direction[0], inward_direction[1]);
    // A reflex point owns precisely the directions whose projections onto both
    // incident outgoing edges are nonpositive. Equality includes sector ends.
    if (turn == CGAL::COLLINEAR || turn == polygon_.orientation()
        || direction == Kernel::Vector_2(CGAL::NULL_VECTOR)
        || CGAL::sign(CGAL::scalar_product(direction, previous - p)) == CGAL::POSITIVE
        || CGAL::sign(CGAL::scalar_product(direction, next - p)) == CGAL::POSITIVE) {
        throw InvalidBoundaryVertexDirectionError("Direction must lie in a reflex vertex's point-site normal cone.");
    }
    return construct(p, direction, FT(tool_radius), CGAL::sqrt(direction.squared_length()));
}

BoundaryNormalCircleProposal2 BoundaryNormalCircle2::at_contact(
    const Point& q, double tool_radius) const
{
    if (!std::isfinite(tool_radius) || CGAL::sign(FT(tool_radius)) != CGAL::POSITIVE) {
        throw InvalidBoundaryCircleContactError("Contact requires a positive finite cutter radius.");
    }
    const FT tool(tool_radius), squared = tool * tool;
    if (polygon_.bounded_side(q) != CGAL::ON_BOUNDED_SIDE) {
        throw InvalidBoundaryCircleContactError("Contact is not inside its boundary polygon.");
    }
    std::optional<Point> source;
    std::vector<std::size_t> vertices, segments;
    for (std::size_t index = 0; index < polygon_.size(); ++index) {
        const Kernel::Segment_2 segment(polygon_[index], polygon_[(index + 1) % polygon_.size()]);
        const auto relation = CGAL::compare(CGAL::squared_distance(q, segment), squared);
        if (relation == CGAL::SMALLER) {
            throw InvalidBoundaryCircleContactError("Contact cutter crosses the polygon boundary.");
        }
        if (relation != CGAL::EQUAL) continue;
        const auto direction = segment.to_vector();
        const FT projection = CGAL::scalar_product(q - segment.source(), direction);
        const bool at_source = CGAL::sign(projection) != CGAL::POSITIVE;
        const bool at_target = CGAL::compare(projection, direction.squared_length()) != CGAL::SMALLER;
        const Point foot = at_source ? segment.source() : at_target ? segment.target()
            : segment.source() + (projection / direction.squared_length()) * direction;
        if (!source) source = foot;
        else if (foot != *source) {
            if (at_source) vertices.push_back(index);
            else if (at_target) vertices.push_back((index + 1) % polygon_.size());
            else segments.push_back(index);
        }
    }
    if (!source) {
        throw InvalidBoundaryCircleContactError("Contact does not lie on this radius's cutter-centre boundary.");
    }
    if (!vertices.empty() || !segments.empty()) {
        // Full scan proved B(q,r) contained with two distinct boundary feet.
        // Along p->q, expanding tangent disks nest inside B(q,r): no earlier
        // competing contact is possible. Thus q itself is the first medial
        // point and the q-m diameter vanishes. Same-foot incident-edge ties
        // at reflex sector joins do not authorize this stationary event.
        return BoundaryNormalCircleProposal2(*source, q, q, q, FT(0), tool,
                                              std::move(vertices), std::move(segments));
    }
    auto result = construct(*source, q - *source, tool, tool, true);
    if (result.exact_contact() != q) {
        throw BoundaryNormalConstructionError("Inverse contact construction changed the exact contact.");
    }
    return result;
}

BoundaryNormalCircleProposal2 BoundaryNormalCircle2::construct(
    const Point& p, const Kernel::Vector_2& normal, const FT& tool, const FT& normal_length, bool allow_stationary) const
{
    // The inverse query has already proved |q-p| equals the cutter radius.
    // Retain that exact length instead of rebuilding nested square roots.
    const FT normal_squared = normal_length * normal_length;
    std::optional<FT> first;
    const auto consider = [&](const FT& t) {
        if (CGAL::sign(t) == CGAL::POSITIVE
            && (!first || CGAL::compare(t, *first) == CGAL::SMALLER)) {
            first = t;
        }
    };

    // An expanding disk centered at p+t*n remains tangent to the source.
    // A point v ties it when |v-p|² = 2*t*n.(v-p).
    for (const Point& vertex : polygon_) {
        const auto delta = vertex - p;
        const FT denominator = FT(2) * CGAL::scalar_product(normal, delta);
        if (CGAL::sign(denominator) == CGAL::POSITIVE) {
            const FT t = delta.squared_length() / denominator;
            consider(t);
        }
    }
    for (std::size_t index = 0; index < polygon_.size(); ++index) {
        const Kernel::Segment_2 segment(polygon_[index], polygon_[(index + 1) % polygon_.size()]);
        const auto other_normal = segment.to_vector().perpendicular(CGAL::COUNTERCLOCKWISE);
        const FT at_p = CGAL::scalar_product(other_normal, p - segment.source());
        const FT slope = CGAL::scalar_product(other_normal, normal);
        const FT norm_product = CGAL::sqrt(normal_squared * other_normal.squared_length());
        // Signed supporting-line distance equals either +t*|n| or -t*|n|.
        // Endpoints are handled above; an infinite supporting line alone cannot
        // authorize an event whose perpendicular foot lies outside its segment.
        for (int sign : {-1, 1}) {
            const FT denominator = FT(sign) * norm_product - slope;
            if (CGAL::is_zero(denominator)) continue;
            const FT t = at_p / denominator;
            if (CGAL::sign(t) != CGAL::POSITIVE) continue;
            const Point candidate = p + t * normal;
            if (interior_foot(candidate, segment)) consider(t);
        }
    }
    if (!first) {
        throw BoundaryNormalConstructionError("A bounded polygon normal has no positive medial contact.");
    }
    const Point m = p + *first * normal;
    const FT clearance = *first * normal_length;
    const FT clearance_squared = clearance * clearance;
    std::vector<std::size_t> vertices, segments;
    if (polygon_.bounded_side(m) != CGAL::ON_BOUNDED_SIDE) {
        throw BoundaryNormalConstructionError("First medial contact is not inside its polygon.");
    }
    for (std::size_t index = 0; index < polygon_.size(); ++index) {
        const Kernel::Segment_2 segment(polygon_[index], polygon_[(index + 1) % polygon_.size()]);
        const auto distance = CGAL::compare(CGAL::squared_distance(m, segment), clearance_squared);
        if (distance == CGAL::SMALLER) {
            throw BoundaryNormalConstructionError("Proposed medial disk crosses a closer boundary feature.");
        }
        if (distance == CGAL::EQUAL && interior_foot(m, segment)) {
            const auto direction = segment.to_vector();
            const FT projection = CGAL::scalar_product(m - segment.source(), direction) / direction.squared_length();
            if (segment.source() + projection * direction != p) segments.push_back(index);
        }
        if (polygon_[index] != p
            && CGAL::compare(CGAL::squared_distance(m, polygon_[index]), clearance_squared) == CGAL::EQUAL) vertices.push_back(index);
    }
    if (vertices.empty() && segments.empty()) {
        throw BoundaryNormalConstructionError("Medial contact lacks a competing boundary feature.");
    }
    if (CGAL::compare(clearance, tool) == CGAL::SMALLER
        || (!allow_stationary && CGAL::compare(clearance, tool) == CGAL::EQUAL)) {
        throw NoPositiveBoundaryCircleError("Medial clearance does not exceed the cutter radius.");
    }
    const Point q = p + (tool / normal_length) * normal;
    const Point center = CGAL::midpoint(q, m);
    const FT guide = (clearance - tool) / FT(2);
    // The complete outer cutter disk is nested in the boundary-tangent medial
    // disk: |m-center| + guide + tool = clearance. All checks stay exact.
    if (CGAL::compare(CGAL::squared_distance(m, p), clearance_squared) != CGAL::EQUAL
        || CGAL::compare(CGAL::squared_distance(q, p), tool * tool) != CGAL::EQUAL
        || CGAL::compare(CGAL::squared_distance(center, q), guide * guide) != CGAL::EQUAL
        || CGAL::compare(CGAL::squared_distance(center, m), guide * guide) != CGAL::EQUAL
        || CGAL::compare(FT(2) * guide + tool, clearance) != CGAL::EQUAL) {
        throw BoundaryNormalConstructionError("Boundary-circle diameter or clearance invariant failed.");
    }
    return BoundaryNormalCircleProposal2(p, m, q, center, guide, clearance, std::move(vertices), std::move(segments));
}

BoundaryNormalCircleProposal2::BoundaryNormalCircleProposal2(
    Point p, Point m, Point q, Point center, FT guide_radius, FT clearance,
    std::vector<std::size_t> vertices, std::vector<std::size_t> segments, std::vector<std::size_t> arcs)
    : p_(std::move(p)), m_(std::move(m)), q_(std::move(q)), center_(std::move(center)),
      guide_radius_(std::move(guide_radius)), clearance_(std::move(clearance)),
      vertices_(std::move(vertices)), segments_(std::move(segments)), arcs_(std::move(arcs)) {}

XY BoundaryNormalCircleProposal2::p_mm() const { return report(p_); }
XY BoundaryNormalCircleProposal2::m_mm() const { return report(m_); }
XY BoundaryNormalCircleProposal2::q_mm() const { return report(q_); }
XY BoundaryNormalCircleProposal2::center_mm() const { return report(center_); }
double BoundaryNormalCircleProposal2::guide_radius_mm() const { return CGAL::to_double(guide_radius_); }
double BoundaryNormalCircleProposal2::clearance_mm() const { return CGAL::to_double(clearance_); }
const std::vector<std::size_t>& BoundaryNormalCircleProposal2::competing_vertex_indices() const { return vertices_; }
const std::vector<std::size_t>& BoundaryNormalCircleProposal2::competing_segment_indices() const { return segments_; }

} // namespace boundary_normal
