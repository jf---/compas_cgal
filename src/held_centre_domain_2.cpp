#include "held_centre_domain_2.h"

#include "gps_ccb_repair_2.h"

#include <CGAL/Fraction_traits.h>

#include <algorithm>
#include <chrono>
#include <limits>
#include <map>
#include <iterator>
#include <sstream>
#include <utility>
#include <variant>

namespace compas_cgal::held_centre_domain {

namespace {

using ET = FT::ET;
using FractionTraits = CGAL::Fraction_traits<ET>;
using Integer = FractionTraits::Numerator_type;
using Circle = Kernel::Circle_2;
using Vector = Kernel::Vector_2;
using CoordNT = OneRootPoint::CoordNT;
using MakeXMonotoneResult = std::variant<OneRootPoint, XCurve>;

Integer integer_from_text(const std::string& text)
{
    if (text.empty()) throw InvalidCentreDomainInputError("empty integer text");
    const std::size_t first = text[0] == '-' ? 1 : 0;
    if (first == text.size()) throw InvalidCentreDomainInputError("integer text has no digits");
    for (std::size_t i = first; i < text.size(); ++i) {
        if (text[i] < '0' || text[i] > '9') throw InvalidCentreDomainInputError("integer text must be signed decimal");
    }
    return Integer(text);
}

std::string integer_to_text(const Integer& value)
{
    std::ostringstream stream;
    stream << value;
    return stream.str();
}

OneRootText one_root_to_text(const CoordNT& value)
{
    if (!value.is_extended()) {
        return {rational_to_text(value.a0()), rational_to_text(FT(0)), rational_to_text(FT(0))};
    }
    return {rational_to_text(value.a0()), rational_to_text(value.a1()), rational_to_text(value.root())};
}

PointText point_to_text(const OneRootPoint& point) { return {one_root_to_text(point.x()), one_root_to_text(point.y())}; }

OneRootPoint lift(const RationalPoint& point) { return OneRootPoint(point.x(), point.y()); }

/// centre + factor * (point - centre): the radial image of a point of a
/// rational circle on the concentric circle of radius factor * radius.
RationalPoint radial(const RationalPoint& centre, const RationalPoint& point, const FT& factor)
{
    return centre + (point - centre) * factor;
}

void append_segment(Polygon& polygon, const RationalPoint& source, const RationalPoint& target)
{
    if (source == target) throw InvalidCentreDomainInputError("zero-length segment in a centre-domain member");
    polygon.push_back(XCurve(source, target));
}

void append_arc(Polygon& polygon, const Traits& traits, const RationalPoint& centre, const FT& radius,
                CGAL::Orientation orientation, const RationalPoint& source, const RationalPoint& target)
{
    const Curve curve(Circle(centre, radius * radius, orientation), lift(source), lift(target));
    std::vector<MakeXMonotoneResult> parts;
    traits.make_x_monotone_2_object()(curve, std::back_inserter(parts));
    for (const MakeXMonotoneResult& part : parts) {
        const XCurve* x_curve = std::get_if<XCurve>(&part);
        if (x_curve == nullptr) throw InvalidCentreDomainInputError("arc degenerates to an isolated point");
        polygon.push_back(*x_curve);
    }
}

Polygon disk(const RationalPoint& centre, const FT& radius)
{
    // Split at the two x-extreme points, which are rational because the radius is.
    const Circle circle(centre, radius * radius, CGAL::COUNTERCLOCKWISE);
    const OneRootPoint left(centre.x() - radius, centre.y());
    const OneRootPoint right(centre.x() + radius, centre.y());
    Polygon polygon;
    polygon.push_back(XCurve(circle, left, right, CGAL::COUNTERCLOCKWISE));
    polygon.push_back(XCurve(circle, right, left, CGAL::COUNTERCLOCKWISE));
    return polygon;
}

/// The rectangle of centres whose nearest point on the line primitive is
/// interior to it: between the two r-offsets, closed by the end normals (each
/// a diameter of the junction disk there).
Polygon line_band(const Primitive& line, const Vector& normal, const FT& r)
{
    const Vector offset = normal * r;
    Polygon polygon;
    append_segment(polygon, line.start - offset, line.end - offset);
    append_segment(polygon, line.end - offset, line.end + offset);
    append_segment(polygon, line.end + offset, line.start + offset);
    append_segment(polygon, line.start + offset, line.start - offset);
    return polygon;
}

/// The annular sector max(rho - r, 0) <= |x - c| <= rho + r over the arc's
/// sweep, closed by the radial segments at its ends (diameters of the junction
/// disks). rho <= r gives a pie sector with its apex at the centre: an exact
/// comparison, structural, not a tolerance.
Polygon arc_band(const Primitive& arc, const FT& r, const Traits& traits)
{
    const FT outer = arc.radius + r;
    const bool pie = CGAL::compare(arc.radius, r) != CGAL::LARGER;
    const FT inner = pie ? FT(0) : FT(arc.radius - r);
    const RationalPoint start_outer = radial(arc.centre, arc.start, outer / arc.radius);
    const RationalPoint end_outer = radial(arc.centre, arc.end, outer / arc.radius);
    const RationalPoint start_inner = pie ? arc.centre : radial(arc.centre, arc.start, inner / arc.radius);
    const RationalPoint end_inner = pie ? arc.centre : radial(arc.centre, arc.end, inner / arc.radius);
    // CCW traversal of the sector: outer flank CCW, inner flank CW.
    const bool ccw = arc.sense > 0;
    const RationalPoint& first_outer = ccw ? start_outer : end_outer;
    const RationalPoint& last_outer = ccw ? end_outer : start_outer;
    const RationalPoint& first_inner = ccw ? end_inner : start_inner;
    const RationalPoint& last_inner = ccw ? start_inner : end_inner;
    Polygon polygon;
    append_arc(polygon, traits, arc.centre, outer, CGAL::COUNTERCLOCKWISE, first_outer, last_outer);
    if (pie) {
        append_segment(polygon, last_outer, arc.centre);
        append_segment(polygon, arc.centre, first_outer);
        return polygon;
    }
    append_segment(polygon, last_outer, first_inner);
    append_arc(polygon, traits, arc.centre, inner, CGAL::CLOCKWISE, first_inner, last_inner);
    append_segment(polygon, last_inner, first_outer);
    return polygon;
}

const Vector& line_normal(const Primitive& line)
{
    // The interior normal must be exactly unit and exactly the left normal of
    // the travel direction.
    const Vector& normal = line.normal;
    const Vector travel = line.end - line.start;
    if (CGAL::compare(normal.squared_length(), FT(1)) != CGAL::EQUAL) {
        throw InvalidCentreDomainInputError("line normal is not exactly unit");
    }
    if (CGAL::sign(travel * normal) != CGAL::ZERO ||
        CGAL::orientation(travel, normal) != CGAL::LEFT_TURN) {
        throw InvalidCentreDomainInputError("line normal is not the exact left normal of its travel direction");
    }
    return normal;
}

struct DirectedEdge {
    const XCurve* curve;
    OneRootPoint source;
    OneRootPoint target;
    int sense;
};

template <typename Circulator>
std::vector<DirectedEdge> directed_edges(Circulator first)
{
    std::vector<DirectedEdge> edges;
    Circulator current = first;
    do {
        const XCurve& curve = current->curve();
        const bool along = (current->direction() == CGAL::ARR_LEFT_TO_RIGHT) == curve.is_directed_right();
        int sense = 0;
        if (curve.is_circular()) {
            const bool ccw = curve.orientation() == CGAL::COUNTERCLOCKWISE;
            sense = (ccw == along) ? 1 : -1;
        }
        edges.push_back({&curve, current->source()->point(), current->target()->point(), sense});
    } while (++current != first);
    for (std::size_t i = 0; i < edges.size(); ++i) {
        if (edges[i].target != edges[(i + 1) % edges.size()].source) {
            throw BrokenCentreDomainChainError("centre-domain boundary cycle does not chain exactly");
        }
    }
    return edges;
}

bool same_support(const DirectedEdge& first, const DirectedEdge& second)
{
    return first.sense == second.sense && first.curve->has_same_supporting_curve(*second.curve);
}

std::array<RationalText, 3> support_text(const XCurve& curve)
{
    if (curve.is_circular()) {
        const Circle circle = curve.supporting_circle();
        return {rational_to_text(circle.center().x()), rational_to_text(circle.center().y()),
                rational_to_text(circle.squared_radius())};
    }
    const auto line = curve.supporting_line();
    return {rational_to_text(line.a()), rational_to_text(line.b()), rational_to_text(line.c())};
}

/// Exact (centre, R^2) -> rational R for every circle the construction used,
/// so a circle piece of the result knows its rational radius without a root.
class RadiusIndex {
public:
    void add(const RationalPoint& centre, const FT& radius) { radii_.emplace(Key{centre, radius * radius}, radius); }

    [[nodiscard]] const FT& radius_of(const Circle& circle) const
    {
        const auto found = radii_.find(Key{circle.center(), circle.squared_radius()});
        if (found == radii_.end()) throw BrokenCentreDomainChainError("a circle of the centre domain was never constructed");
        return found->second;
    }

private:
    struct Key {
        RationalPoint centre;
        FT radius_squared;
    };
    struct Less {
        bool operator()(const Key& first, const Key& second) const
        {
            const CGAL::Comparison_result by_point = CGAL::compare_xy(first.centre, second.centre);
            if (by_point != CGAL::EQUAL) return by_point == CGAL::SMALLER;
            return CGAL::compare(first.radius_squared, second.radius_squared) == CGAL::SMALLER;
        }
    };
    std::map<Key, FT, Less> radii_;
};

OneRoot one_root(const CoordNT& value)
{
    if (!value.is_extended()) return {value.a0(), FT(0), FT(0)};
    return {value.a0(), value.a1(), value.root()};
}

void include(QueryPiece& piece, const std::pair<double, double>& x, const std::pair<double, double>& y)
{
    piece.xmin = std::min(piece.xmin, x.first);
    piece.xmax = std::max(piece.xmax, x.second);
    piece.ymin = std::min(piece.ymin, y.first);
    piece.ymax = std::max(piece.ymax, y.second);
}

/// Grow the certified box of `piece` by one x-monotone edge: its endpoints,
/// and for an arc the top or bottom of its circle (an x-monotone arc spans
/// exactly one of the two half circles, and only its extreme y can leave the
/// endpoint box).
void include_edge(QueryPiece& piece, const DirectedEdge& edge, const RadiusIndex& radii)
{
    for (const OneRootPoint* point : {&edge.source, &edge.target}) {
        include(piece, enclosure(one_root(point->x())), enclosure(one_root(point->y())));
    }
    if (!edge.curve->is_circular()) return;
    const XCurve& curve = *edge.curve;
    const Circle circle = curve.supporting_circle();
    const bool upper = (curve.orientation() == CGAL::COUNTERCLOCKWISE) != curve.is_directed_right();
    const FT extreme = upper ? FT(circle.center().y() + radii.radius_of(circle)) : FT(circle.center().y() - radii.radius_of(circle));
    const std::pair<double, double> y = CGAL::to_interval(extreme);
    piece.ymin = std::min(piece.ymin, y.first);
    piece.ymax = std::max(piece.ymax, y.second);
}

QueryPiece query_piece(const std::vector<DirectedEdge>& edges, std::size_t first, std::size_t count,
                       const RadiusIndex& radii)
{
    const std::size_t size = edges.size();
    const DirectedEdge& head = edges[first % size];
    const DirectedEdge& tail = edges[(first + count - 1) % size];
    QueryPiece piece;
    PieceNumbers<FT>& n = piece.exact;
    piece.circular = head.curve->is_circular();
    piece.full = count == size && head.source == tail.target && piece.circular;
    piece.sense = head.sense;
    n.sx = one_root(head.source.x());
    n.sy = one_root(head.source.y());
    n.tx = one_root(tail.target.x());
    n.ty = one_root(tail.target.y());
    n.a = n.b = n.c = n.norm_squared = n.dx = n.dy = n.cx = n.cy = n.radius = FT(0);
    const double infinity = std::numeric_limits<double>::infinity();
    piece.xmin = piece.ymin = infinity;
    piece.xmax = piece.ymax = -infinity;
    for (std::size_t k = 0; k < count; ++k) include_edge(piece, edges[(first + k) % size], radii);
    if (piece.circular) {
        const Circle circle = head.curve->supporting_circle();
        n.cx = circle.center().x();
        n.cy = circle.center().y();
        n.radius = radii.radius_of(circle);
        finish(piece);
        return piece;
    }
    const auto line = head.curve->supporting_line();
    n.a = line.a();
    n.b = line.b();
    n.c = line.c();
    n.norm_squared = n.a * n.a + n.b * n.b;
    // Travel direction: (-b, a) or its opposite, whichever advances from source to target.
    const FT ux = -n.b;
    const FT uy = n.a;
    const OneRoot along_target{ux * n.tx.a + uy * n.ty.a, ux * n.tx.b + uy * n.ty.b,
                               CGAL::is_zero(n.tx.b) ? n.ty.r : n.tx.r};
    const OneRoot along_source{ux * n.sx.a + uy * n.sy.a, ux * n.sx.b + uy * n.sy.b,
                               CGAL::is_zero(n.sx.b) ? n.sy.r : n.sx.r};
    const int advance = compare(along_target, along_source);
    if (advance == 0) throw BrokenCentreDomainChainError("a line piece of the centre domain has zero length");
    n.dx = advance > 0 ? ux : FT(-ux);
    n.dy = advance > 0 ? uy : FT(-uy);
    finish(piece);
    return piece;
}

/// Merge consecutive halfedges on one support into maximal pieces, cyclically.
/// Appends the query form of every piece to `queries`; piece k of the cycle has
/// global index base + k, and its target is the next piece's source vertex.
template <typename Circulator>
Cycle extract_cycle(Circulator first, const RadiusIndex& radii, std::vector<QueryPiece>& queries)
{
    const std::vector<DirectedEdge> edges = directed_edges(first);
    const std::size_t count = edges.size();
    const std::size_t base = queries.size();
    std::size_t start = 0;
    while (start < count && same_support(edges[(start + count - 1) % count], edges[start])) ++start;
    Cycle cycle;
    if (start == count) {
        // The whole cycle lies on one support (a full circle).
        const DirectedEdge& edge = edges.front();
        cycle.push_back({edge.curve->is_circular(), support_text(*edge.curve), edge.sense, point_to_text(edge.source),
                         point_to_text(edge.source)});
        QueryPiece piece = query_piece(edges, 0, count, radii);
        piece.target_vertex = base;
        queries.push_back(std::move(piece));
        return cycle;
    }
    std::size_t i = 0;
    while (i < count) {
        const DirectedEdge& head = edges[(start + i) % count];
        std::size_t j = i + 1;
        while (j < count && same_support(edges[(start + j - 1) % count], edges[(start + j) % count])) ++j;
        const DirectedEdge& tail = edges[(start + j - 1) % count];
        cycle.push_back({head.curve->is_circular(), support_text(*head.curve), head.sense, point_to_text(head.source),
                         point_to_text(tail.target)});
        queries.push_back(query_piece(edges, start + i, j - i, radii));
        i = j;
    }
    for (std::size_t k = 0; k < cycle.size(); ++k) queries[base + k].target_vertex = base + (k + 1) % cycle.size();
    return cycle;
}

void validate(const std::vector<Primitive>& primitives, const FT& tool_radius)
{
    if (primitives.size() < 2) throw InvalidCentreDomainInputError("the declared cycle needs at least two primitives");
    if (CGAL::sign(tool_radius) != CGAL::POSITIVE) throw InvalidCentreDomainInputError("tool radius must be positive");
    for (std::size_t i = 0; i < primitives.size(); ++i) {
        const Primitive& primitive = primitives[i];
        if (primitive.end != primitives[(i + 1) % primitives.size()].start) {
            throw InvalidCentreDomainInputError("the declared cycle does not close exactly");
        }
        if (primitive.start == primitive.end) throw InvalidCentreDomainInputError("zero-length primitive");
        if (primitive.kind != PrimitiveKind::Arc) continue;
        if (CGAL::sign(primitive.radius) != CGAL::POSITIVE || (primitive.sense != 1 && primitive.sense != -1)) {
            throw InvalidCentreDomainInputError("arc needs a positive radius and a sense of +1 or -1");
        }
        const FT radius_squared = primitive.radius * primitive.radius;
        if (CGAL::squared_distance(primitive.start, primitive.centre) != radius_squared ||
            CGAL::squared_distance(primitive.end, primitive.centre) != radius_squared) {
            throw InvalidCentreDomainInputError("arc ends are not exactly on the arc's circle");
        }
    }
}

}  // namespace

FT rational_from_text(const std::string& numerator, const std::string& denominator)
{
    const Integer num = integer_from_text(numerator);
    const Integer den = integer_from_text(denominator);
    if (CGAL::sign(den) != CGAL::POSITIVE) throw InvalidCentreDomainInputError("rational denominator must be positive");
    return FT(FractionTraits::Compose()(num, den));
}

RationalText rational_to_text(const FT& value)
{
    Integer num;
    Integer den;
    FractionTraits::Decompose()(CGAL::exact(value), num, den);
    return {integer_to_text(num), integer_to_text(den)};
}

std::unique_ptr<CentreDomain2> CentreDomain2::build(const std::vector<Primitive>& primitives, const FT& tool_radius)
{
    const auto started = std::chrono::steady_clock::now();
    validate(primitives, tool_radius);
    const Traits traits;
    RadiusIndex radii;
    Polygon pocket;
    std::vector<Polygon> members;
    members.reserve(2 * primitives.size());
    for (const Primitive& primitive : primitives) {
        // Junction disk at the start of every primitive: every junction, so
        // every band end cap (the end normals are its diameters) is closed.
        members.push_back(disk(primitive.start, tool_radius));
        radii.add(primitive.start, tool_radius);
        if (primitive.kind == PrimitiveKind::Line) {
            append_segment(pocket, primitive.start, primitive.end);
            members.push_back(line_band(primitive, line_normal(primitive), tool_radius));
            continue;
        }
        append_arc(pocket, traits, primitive.centre, primitive.radius,
                   primitive.sense > 0 ? CGAL::COUNTERCLOCKWISE : CGAL::CLOCKWISE, primitive.start, primitive.end);
        members.push_back(arc_band(primitive, tool_radius, traits));
        radii.add(primitive.centre, primitive.radius);
        radii.add(primitive.centre, primitive.radius + tool_radius);
        if (CGAL::compare(primitive.radius, tool_radius) == CGAL::LARGER) radii.add(primitive.centre, primitive.radius - tool_radius);
    }
    std::unique_ptr<CentreDomain2> result(new CentreDomain2());
    Gps removed;
    removed.join(members.begin(), members.end());
    result->domain_ = std::make_unique<Gps>(pocket);
    result->domain_->difference(removed);
    // The Gps boolean can file a hole as a second outer CCB (gps_ccb_repair_2.h);
    // the cycle extraction below and every later reader need the planar
    // invariant, so it is restored exactly where it breaks.
    result->refiled_outer_ccbs_ = refile_misfiled_outer_ccbs<CentreDomainCcbRepairDomain>(*result->domain_);
    result->member_count_ = members.size();

    const auto& arrangement = result->domain_->arrangement();
    for (auto face = arrangement.faces_begin(); face != arrangement.faces_end(); ++face) {
        if (!face->contained()) continue;
        if (!face->has_outer_ccb()) throw BrokenCentreDomainChainError("a contained face of the centre domain is unbounded");
        if (face->isolated_vertices_begin() != face->isolated_vertices_end()) {
            throw BrokenCentreDomainChainError("a contained face of the centre domain has an isolated vertex");
        }
        // Post-condition of the repair above: `outer_ccb()` requires exactly one
        // outer CCB and that precondition is compiled out of release builds, so
        // a planar-invariant breach would otherwise silently drop a cycle.
        if (face->number_of_outer_ccbs() != 1) {
            throw BrokenCentreDomainChainError("a contained face of the centre domain does not have exactly one outer boundary cycle");
        }
        Component component;
        component.cycles.push_back(extract_cycle(face->outer_ccb(), radii, result->query_pieces_));
        for (auto hole = face->inner_ccbs_begin(); hole != face->inner_ccbs_end(); ++hole) {
            component.cycles.push_back(extract_cycle(*hole, radii, result->query_pieces_));
        }
        result->components_.push_back(std::move(component));
    }
    result->construction_seconds_ =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    return result;
}

NearestDecision CentreDomain2::decide_nearest(const RationalPoint& point, const FT* gap_squared,
                                              QueryWorkspace& workspace) const
{
    return held_centre_domain::decide_nearest(query_pieces_, point, gap_squared, workspace);
}

NearestValues CentreDomain2::nearest_values(const RationalPoint& point, const NearestDecision& decision) const
{
    return held_centre_domain::nearest_values(query_pieces_, point, decision);
}

Location CentreDomain2::locate(const RationalPoint& point) const
{
    switch (domain_->oriented_side(lift(point))) {
    case CGAL::ON_POSITIVE_SIDE:
        return Location::Interior;
    case CGAL::ON_ORIENTED_BOUNDARY:
        return Location::Boundary;
    default:
        return Location::Exterior;
    }
}

}  // namespace compas_cgal::held_centre_domain
