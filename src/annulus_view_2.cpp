#include "annulus_view_2.h"

#include "engagement_stations_2.h"
#include "exact_rim_order_2.h"
#include "stock_2.h"
#include "swept_prefix_engagement_2.h"

#include <CGAL/AABB_segment_primitive_2.h>
#include <CGAL/AABB_traits_2.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

namespace {

using FT = Epeck::FT;
using CoordNT = GpsPoint::CoordNT;
using RimArc = std::pair<GpsPoint, GpsPoint>;

// A rim-cell witness refines by dyadic bisection of a rational half-angle, which
// reaches any non-empty cell of algebraic endpoints. Hitting this bound is a
// construction breach and raises. It decides nothing about the cell.
constexpr int WITNESS_REFINEMENT_LIMIT = 4096;

// ---- candidate search: outward-rounded boxes in CGAL AABB trees ----------------
//
// Every box coordinate is an injected double or a to_interval bound, so a box miss
// proves an exact miss. A hit only makes a candidate: selection is never a decision.

using SearchKernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using SearchBox = SearchKernel::Iso_rectangle_2;
using SearchSegments = std::vector<SearchKernel::Segment_2>;
using EdgePrimitive =
    CGAL::AABB_segment_primitive_2<SearchKernel, SearchSegments::const_iterator>;
using EdgeTree = CGAL::AABB_tree<CGAL::AABB_traits_2<SearchKernel, EdgePrimitive>>;

struct AnnulusBox2 {
    std::size_t id;
    SearchBox box;
};

struct AnnulusBoxPrimitive {
    using Id = std::size_t;
    using Point = SearchKernel::Point_2;
    using Datum = SearchBox;

    Id box_id;
    SearchBox box;

    explicit AnnulusBoxPrimitive(std::vector<AnnulusBox2>::const_iterator item)
        : box_id(item->id), box(item->box)
    {
    }

    Id id() const { return box_id; }
    const Datum& datum() const { return box; }
    Point reference_point() const { return box.min(); }
};
using AnnulusTree =
    CGAL::AABB_tree<CGAL::AABB_traits_2<SearchKernel, AnnulusBoxPrimitive>>;

// The closed swept annulus of one emitted circle, on full_turn_annulus_bounds' exact radii.
struct ExactAnnulus2 {
    EPoint centre;
    FT inner_squared;
    FT outer_squared;
    bool has_inner;
    // Outward-rounded enclosures of the values above, for the certified shortcut.
    double centre_x;
    double centre_y;
    CGAL::Interval_nt<false> inner;
    CGAL::Interval_nt<false> outer;
};

// How a rim lies relative to a candidate annulus, certified by outward-rounded
// intervals. Each of the three definite answers is a strict inequality on enclosures
// of the exact quantities, so it holds exactly (the Filtered_predicate pattern); the
// distance from the annulus centre to the rim runs over [|d - r|, d + r].
//   CLEAR     d + r < inner, or d - r > outer: the rim misses the closed annulus.
//   COVERED   d - r > inner and d + r < outer: the rim is inside it, no material.
//   UNCERTAIN everything else, decided by the exact path.
enum class RimAnnulusRelation2 { CLEAR, COVERED, UNCERTAIN };

RimAnnulusRelation2 certified_rim_relation(
    const CGAL::Interval_nt<false>& station_x,
    const CGAL::Interval_nt<false>& station_y,
    const CGAL::Interval_nt<false>& radius,
    const ExactAnnulus2& annulus)
{
    using Interval = CGAL::Interval_nt<false>;
    const Interval dx = station_x - Interval(annulus.centre_x);
    const Interval dy = station_y - Interval(annulus.centre_y);
    const Interval squared = CGAL::square(dx) + CGAL::square(dy);
    const auto squared_below = [&squared](const Interval& bound) {
        return bound.inf() > 0.0 && squared.sup() < (bound * bound).inf();
    };
    const auto squared_above = [&squared](const Interval& bound) {
        return bound.inf() > 0.0 && squared.inf() > (bound * bound).sup();
    };
    const Interval far = annulus.outer + radius;
    if (squared_above(far)) return RimAnnulusRelation2::CLEAR;
    if (annulus.has_inner && squared_below(annulus.inner - radius)) {
        return RimAnnulusRelation2::CLEAR;
    }
    if (squared_below(annulus.outer - radius)
        && (!annulus.has_inner || squared_above(annulus.inner + radius))) {
        return RimAnnulusRelation2::COVERED;
    }
    return RimAnnulusRelation2::UNCERTAIN;
}

struct PartEdge2 {
    EPoint start;
    EPoint end;
    EVector direction;
    FT squared_length;
};

// One rim event. `upper` is the half the exact rim order puts it in. The two
// x-extreme events get it by construction: their y equals the station's identically,
// so comparing it would be an identically-zero decision, which no interval filter
// resolves and which forces the station's whole exact value.
struct RimEvent2 {
    GpsPoint point;
    bool wall;
    bool upper;
};

bool rim_event_less(const RimEvent2& lhs, const RimEvent2& rhs)
{
    if (lhs.upper != rhs.upper) return lhs.upper;
    const CGAL::Comparison_result x = CGAL::compare(lhs.point.x(), rhs.point.x());
    return lhs.upper ? x == CGAL::LARGER : x == CGAL::SMALLER;
}

// Whether the rational rim point lies strictly inside the CCW cell (start, end).
// Both extremes are events, so a cell lies in the one half its start belongs to and
// is x-monotone there.
bool strictly_inside_rim_cell(
    const EPoint& point,
    const RimEvent2& start,
    const RimEvent2& end,
    const EPoint& station)
{
    const CGAL::Comparison_result side = CGAL::compare(point.y(), station.y());
    if (side != (start.upper ? CGAL::LARGER : CGAL::SMALLER)) return false;
    const CoordNT x(point.x());
    const GpsPoint& low = start.upper ? end.point : start.point;
    const GpsPoint& high = start.upper ? start.point : end.point;
    return CGAL::compare(x, low.x()) == CGAL::LARGER
        && CGAL::compare(x, high.x()) == CGAL::SMALLER;
}

bool in_closed_annulus(const EPoint& point, const ExactAnnulus2& annulus)
{
    if (CGAL::compare_squared_distance(point, annulus.centre, annulus.outer_squared)
        == CGAL::LARGER) {
        return false;
    }
    return !annulus.has_inner
        || CGAL::compare_squared_distance(point, annulus.centre, annulus.inner_squared)
        != CGAL::SMALLER;
}

// One rim against one annulus centre: the quantities both of the annulus's circles share.
struct RimCentrePair2 {
    const EPoint& station;
    EVector v;
    FT distance_squared;
    FT radius_ratio;  // r^2 / d^2
    FT half_inverse;  // 1 / (2 d^2)
    FT shifted;       // d^2 + r^2
};

// The rim's intersections with the circle of squared radius `squared_radius` about
// the pair's centre. With v = centre - s and d^2 = |v|^2 the points are
// s + a v +- sqrt(q) v_perp, where a = (d^2 + r^2 - R^2) / (2 d^2) and
// q = r^2 / d^2 - a^2: one root for both coordinates.
std::size_t add_circle_events(
    std::vector<RimEvent2>& events,
    const RimCentrePair2& pair,
    const FT& squared_radius)
{
    const FT a = (pair.shifted - squared_radius) * pair.half_inverse;
    const FT q = pair.radius_ratio - a * a;
    const CGAL::Sign sign = CGAL::sign(q);
    if (sign == CGAL::NEGATIVE) return 0;
    const FT mid_x = pair.station.x() + a * pair.v.x();
    const FT mid_y = pair.station.y() + a * pair.v.y();
    if (sign == CGAL::ZERO) {
        GpsPoint point(mid_x, mid_y);
        const bool upper = rim_upper(point, pair.station);
        events.push_back({std::move(point), false, upper});
        return 1;
    }
    for (const int branch : {1, -1}) {
        const FT unit(branch);
        GpsPoint point(
            CoordNT(mid_x, -unit * pair.v.y(), q),
            CoordNT(mid_y, unit * pair.v.x(), q));
        const bool upper = rim_upper(point, pair.station);
        events.push_back({std::move(point), false, upper});
    }
    return 2;
}

// The rim's intersections with one edge a + t (b - a), t in [0, 1]. With L = |b - a|^2,
// B = (b - a) . (a - s) and C = |a - s|^2 - r^2, t = (-B +- sqrt(B^2 - L C)) / L.
void add_edge_events(
    std::vector<RimEvent2>& events,
    const EPoint& station,
    const FT& rim_squared,
    const PartEdge2& edge)
{
    const EVector offset = edge.start - station;
    const FT b = edge.direction * offset;
    const FT c = offset.squared_length() - rim_squared;
    const FT disc = b * b - edge.squared_length * c;
    const CGAL::Sign sign = CGAL::sign(disc);
    if (sign == CGAL::NEGATIVE) return;
    const FT x_slope = edge.direction.x() / edge.squared_length;
    const FT y_slope = edge.direction.y() / edge.squared_length;
    const FT mid_x = edge.start.x() - b * x_slope;
    const FT mid_y = edge.start.y() - b * y_slope;
    if (sign == CGAL::ZERO) {
        const CGAL::Sign below = CGAL::sign(-b);
        if (below != CGAL::NEGATIVE
            && CGAL::compare(-b, edge.squared_length) != CGAL::LARGER) {
            GpsPoint point(mid_x, mid_y);
            const bool upper = rim_upper(point, station);
            events.push_back({std::move(point), true, upper});
        }
        return;
    }
    const CoordNT zero{FT(0)};
    const CoordNT length{edge.squared_length};
    for (const int branch : {-1, 1}) {
        const FT unit(branch);
        const CoordNT parameter(-b, unit, disc);
        if (CGAL::compare(parameter, zero) == CGAL::SMALLER
            || CGAL::compare(parameter, length) == CGAL::LARGER) {
            continue;
        }
        GpsPoint point(
            CoordNT(mid_x, unit * x_slope, disc),
            CoordNT(mid_y, unit * y_slope, disc));
        const bool upper = rim_upper(point, station);
        events.push_back({std::move(point), true, upper});
    }
}

SearchBox rim_search_box(const EPoint& station, const FT& radius)
{
    const FT min_x = station.x() - radius;
    const FT min_y = station.y() - radius;
    const FT max_x = station.x() + radius;
    const FT max_y = station.y() + radius;
    return SearchBox(
        CGAL::to_interval(min_x).first,
        CGAL::to_interval(min_y).first,
        CGAL::to_interval(max_x).second,
        CGAL::to_interval(max_y).second);
}

void validate_range(const AnnulusViewRange2& range, std::size_t size)
{
    if (range.first > range.end || range.end > size) {
        throw InvalidAnnulusViewRangeError(
            "annulus view range must satisfy first <= end <= size");
    }
}

bool same_point_run(const std::vector<RimArc>& lhs, const std::vector<RimArc>& rhs)
{
    if (lhs.empty() || rhs.empty()) return lhs.empty() && rhs.empty();
    const bool lhs_closed = lhs.front().first == lhs.back().second;
    const bool rhs_closed = rhs.front().first == rhs.back().second;
    if (lhs_closed || rhs_closed) return lhs_closed && rhs_closed;
    return lhs.front().first == rhs.front().first
        && lhs.back().second == rhs.back().second;
}

StationAgreement2 classifications_agree(
    const AuditExactStationClassification2& lhs,
    const AuditExactStationClassification2& rhs)
{
    const auto& lhs_runs = lhs.true_run_arcs();
    const auto& rhs_runs = rhs.true_run_arcs();
    bool runs = lhs_runs.size() == rhs_runs.size();
    for (const auto& run : lhs_runs) {
        if (!runs) break;
        runs = std::any_of(
            rhs_runs.begin(), rhs_runs.end(),
            [&run](const auto& other) { return same_point_run(run, other); });
    }
    return {lhs.disposition() == rhs.disposition(), runs};
}

} // namespace

struct AnnulusView2::Index {
    FT tool_radius;
    double tool_radius_binary64;
    FT rim_squared;
    double part_max_x;
    std::vector<PartEdge2> edges;
    SearchSegments edge_segments;
    std::vector<ExactAnnulus2> annuli;
    std::vector<AnnulusBox2> annulus_boxes;
    EdgeTree edge_tree;
    AnnulusTree annulus_tree;

    bool strictly_inside_part(const EPoint& point) const;
    EPoint cell_witness(
        const RimEvent2& start, const RimEvent2& end, const EPoint& station) const;
    std::vector<RimArc> harvest(const EPoint& station, AnnulusViewRange2 range) const;
    AuditExactStationClassification2 classify_static(
        const EPoint& station, AnnulusViewRange2 range,
        const FT& cap_chord_ratio, const FT& gap_close_ratio) const;
};

// Crossing parity along the +x ray, half-open in y. The caller guarantees that the
// point is on no edge (an edge through a witness would be an event inside its
// cell); a collinear crossing edge is that breach and raises.
bool AnnulusView2::Index::strictly_inside_part(const EPoint& point) const
{
    const std::pair<double, double> x_range = CGAL::to_interval(point.x());
    if (x_range.first > part_max_x) return false;
    const std::pair<double, double> y_range = CGAL::to_interval(point.y());
    std::vector<EdgeTree::Primitive_id> hits;
    edge_tree.all_intersected_primitives(
        SearchBox(x_range.first, y_range.first, part_max_x, y_range.second),
        std::back_inserter(hits));
    bool inside = false;
    for (const EdgeTree::Primitive_id& hit : hits) {
        const PartEdge2& edge =
            edges[static_cast<std::size_t>(std::distance(edge_segments.cbegin(), hit))];
        const bool start_above = CGAL::compare_y(edge.start, point) == CGAL::LARGER;
        const bool end_above = CGAL::compare_y(edge.end, point) == CGAL::LARGER;
        if (start_above == end_above) continue;
        const CGAL::Orientation side = CGAL::orientation(edge.start, edge.end, point);
        if (side == CGAL::COLLINEAR) {
            throw AnnulusViewWitnessError("a rim-cell witness lies on the part boundary");
        }
        if (side == (end_above ? CGAL::LEFT_TURN : CGAL::RIGHT_TURN)) {
            inside = !inside;
        }
    }
    return inside;
}

// A rational point strictly inside the CCW rim cell (start, end): the rim point
// along an exact rational unit steered by the double middle angle, checked exactly;
// on a cell too thin for that hint, a rational half-angle rotation of that unit,
// refined by exact rim comparisons. Rationals on the rim stay on the rim exactly.
EPoint AnnulusView2::Index::cell_witness(
    const RimEvent2& start, const RimEvent2& end, const EPoint& station) const
{
    // The hint reads interval bounds, not to_double: Lazy_exact_nt's to_double forces
    // the exact value of anything near zero, and y - station.y is near zero at every
    // rim extreme. The angles are taken inside the cell's half, so they increase along
    // the cell with no wrap, and a cell too thin for doubles reads as one angle beside
    // the cell. A hint needs no accuracy: the check below decides.
    const double centre_x = CGAL::to_interval(station.x()).first;
    const double centre_y = CGAL::to_interval(station.y()).first;
    const double half = start.upper ? 1.0 : -1.0;
    const auto angle_in_half = [&](const GpsPoint& point) {
        return std::atan2(
            half * std::abs(CGAL::to_interval(point.y()).first - centre_y),
            CGAL::to_interval(point.x()).first - centre_x);
    };
    const double start_angle = angle_in_half(start.point);
    const double end_angle = angle_in_half(end.point);
    const double middle = 0.5 * (start_angle + end_angle);
    const EVector unit = exact_unit_along(std::cos(middle), std::sin(middle));
    const EPoint hinted = station + tool_radius * unit;
    if (strictly_inside_rim_cell(hinted, start, end, station)) return hinted;
    const EVector normal = unit.perpendicular(CGAL::COUNTERCLOCKWISE);
    const auto rotated = [&](const FT& t) {
        const FT squared = t * t;
        const FT denominator = FT(1) + squared;
        return station
            + (tool_radius / denominator)
            * ((FT(1) - squared) * unit + (FT(2) * t) * normal);
    };
    FT lower(-1);
    FT upper(1);
    EPoint lower_point = rotated(lower);
    for (int step = 0; step < WITNESS_REFINEMENT_LIMIT; ++step) {
        const FT middle_t = (lower + upper) / FT(2);
        const EPoint candidate = rotated(middle_t);
        if (strictly_inside_rim_cell(candidate, start, end, station)) return candidate;
        const GpsPoint on_rim(candidate.x(), candidate.y());
        const bool before_cell = on_rim == start.point
            || (on_rim != end.point
                && strictly_inside_ccw_arc(
                    on_rim, GpsPoint(lower_point.x(), lower_point.y()), start.point,
                    station));
        if (before_cell) {
            lower = middle_t;
            lower_point = candidate;
        } else {
            upper = middle_t;
        }
    }
    throw AnnulusViewWitnessError(
        "no rational witness found inside a rim cell: start angle "
        + std::to_string(start_angle) + ", end angle " + std::to_string(end_angle));
}

std::vector<RimArc> AnnulusView2::Index::harvest(
    const EPoint& station, AnnulusViewRange2 range) const
{
    const SearchBox rim_box = rim_search_box(station, tool_radius);
    const FT extreme_x_low = station.x() - tool_radius;
    const FT extreme_x_high = station.x() + tool_radius;
    std::vector<RimEvent2> events;
    events.push_back({GpsPoint(extreme_x_high, station.y()), false, true});
    events.push_back({GpsPoint(extreme_x_low, station.y()), false, false});
    const EPoint rim_probe(extreme_x_high, station.y());

    // Annuli: an annulus with no event has the whole rim inside it (no material at
    // all, including a rim overlapping a boundary) or none of it, decided at one
    // rim point. Only the annuli the rim crosses stay in the cell tests.
    std::vector<AnnulusTree::Primitive_id> hits;
    annulus_tree.all_intersected_primitives(rim_box, std::back_inserter(hits));
    std::sort(hits.begin(), hits.end(), std::greater<>());
    std::vector<const ExactAnnulus2*> uncertain;
    {
        CGAL::Protect_FPU_rounding<true> rounding;
        const CGAL::Interval_nt<false> station_x(CGAL::to_interval(station.x()));
        const CGAL::Interval_nt<false> station_y(CGAL::to_interval(station.y()));
        const CGAL::Interval_nt<false> radius(CGAL::to_interval(tool_radius));
        for (const AnnulusTree::Primitive_id hit : hits) {
            if (hit < range.first || hit >= range.end) continue;
            const ExactAnnulus2& annulus = annuli[hit];
            const RimAnnulusRelation2 relation =
                certified_rim_relation(station_x, station_y, radius, annulus);
            if (relation == RimAnnulusRelation2::COVERED) return {};
            if (relation == RimAnnulusRelation2::UNCERTAIN) uncertain.push_back(&annulus);
        }
    }
    std::vector<const ExactAnnulus2*> crossing;
    for (const ExactAnnulus2* candidate : uncertain) {
        const ExactAnnulus2& annulus = *candidate;
        std::size_t count = 0;
        const EVector v = annulus.centre - station;
        const FT distance_squared = v.squared_length();
        // A concentric annulus has no event (no common point, or the rim itself, which
        // the membership below decides). CGAL::sign, never CGAL::is_zero: on
        // Lazy_exact_nt is_zero evaluates the exact value unconditionally, and sign is
        // filtered.
        if (CGAL::sign(distance_squared) != CGAL::ZERO) {
            const RimCentrePair2 pair{
                station, v, distance_squared, rim_squared / distance_squared,
                FT(1) / (FT(2) * distance_squared), distance_squared + rim_squared};
            count = add_circle_events(events, pair, annulus.outer_squared);
            if (annulus.has_inner) {
                count += add_circle_events(events, pair, annulus.inner_squared);
            }
        }
        if (count > 0) {
            crossing.push_back(&annulus);
        } else if (in_closed_annulus(rim_probe, annulus)) {
            return {};
        }
    }

    // Walls: with no wall event the whole rim is on one side of P.
    std::vector<EdgeTree::Primitive_id> edge_hits;
    edge_tree.all_intersected_primitives(rim_box, std::back_inserter(edge_hits));
    const std::size_t events_before_walls = events.size();
    for (const EdgeTree::Primitive_id& hit : edge_hits) {
        add_edge_events(
            events, station, rim_squared,
            edges[static_cast<std::size_t>(std::distance(edge_segments.cbegin(), hit))]);
    }
    const bool any_wall = events.size() != events_before_walls;
    if (!any_wall && !strictly_inside_part(rim_probe)) return {};

    std::sort(events.begin(), events.end(), rim_event_less);
    std::vector<RimEvent2> ordered;
    ordered.reserve(events.size());
    for (RimEvent2& event : events) {
        if (!ordered.empty() && ordered.back().point == event.point) {
            ordered.back().wall = ordered.back().wall || event.wall;
        } else {
            ordered.push_back(std::move(event));
        }
    }

    std::vector<RimArc> arcs;
    const std::size_t count = ordered.size();
    // A wall event may change P membership; between two of them it is constant.
    bool part_known = !any_wall;
    bool in_part = !any_wall;
    for (std::size_t index = 0; index < count; ++index) {
        const RimEvent2& start = ordered[index];
        const RimEvent2& end = ordered[(index + 1) % count];
        if (start.wall) part_known = false;
        if (crossing.empty() && part_known) {
            if (in_part) arcs.emplace_back(start.point, end.point);
            continue;
        }
        const EPoint witness = cell_witness(start, end, station);
        const bool removed = std::any_of(
            crossing.begin(), crossing.end(),
            [&witness](const ExactAnnulus2* annulus) {
                return in_closed_annulus(witness, *annulus);
            });
        if (removed) continue;
        if (!part_known) {
            in_part = strictly_inside_part(witness);
            part_known = true;
        }
        if (in_part) arcs.emplace_back(start.point, end.point);
    }
    return arcs;
}

AuditExactStationClassification2 AnnulusView2::Index::classify_static(
    const EPoint& station, AnnulusViewRange2 range,
    const FT& cap_chord_ratio, const FT& gap_close_ratio) const
{
    return classify_audit_filtered_station_exact(
        station, tool_radius, cap_chord_ratio, gap_close_ratio,
        harvest(station, range));
}

AnnulusView2 AnnulusView2::build(
    Eigen::Ref<const compas::RowMatrixXd> polygon,
    Eigen::Ref<const compas::RowMatrixXd> annuli,
    double tool_radius)
{
    if (!std::isfinite(tool_radius)) {
        throw NonFiniteAnnulusViewInputError("annulus view tool radius must be finite");
    }
    if (!(tool_radius > 0.0)) {
        throw InvalidAnnulusViewRadiiError("annulus view tool radius must be positive");
    }
    if (polygon.rows() < 3 || polygon.cols() < 2) {
        throw InvalidAnnulusViewPolygonError(
            "annulus view polygon needs at least three vertices");
    }
    if (annuli.rows() > 0 && annuli.cols() != 3) {
        throw InvalidAnnulusViewRadiiError(
            "annulus view annuli must be rows of (cx, cy, guide_radius)");
    }
    CGAL::Polygon_2<Epeck> part;
    std::vector<SearchKernel::Point_2> search_vertices;
    double part_max_x = polygon(0, 0);
    for (Eigen::Index row = 0; row < polygon.rows(); ++row) {
        const double x = polygon(row, 0);
        const double y = polygon(row, 1);
        if (!std::isfinite(x) || !std::isfinite(y)) {
            throw NonFiniteAnnulusViewInputError(
                "annulus view polygon vertices must be finite");
        }
        part.push_back(EPoint(x, y));
        search_vertices.emplace_back(x, y);
        part_max_x = std::max(part_max_x, x);
    }
    if (!part.is_simple()) {
        throw InvalidAnnulusViewPolygonError("annulus view polygon must be simple");
    }

    auto index = std::make_unique<Index>();
    index->tool_radius = FT(tool_radius);
    index->tool_radius_binary64 = tool_radius;
    index->rim_squared = index->tool_radius * index->tool_radius;
    index->part_max_x = part_max_x;
    const std::size_t vertex_count = part.size();
    index->edges.reserve(vertex_count);
    index->edge_segments.reserve(vertex_count);
    for (std::size_t i = 0; i < vertex_count; ++i) {
        const std::size_t next = (i + 1) % vertex_count;
        const EPoint& start = part[i];
        const EPoint& end = part[next];
        const EVector direction = end - start;
        index->edges.push_back({start, end, direction, direction.squared_length()});
        index->edge_segments.emplace_back(search_vertices[i], search_vertices[next]);
    }
    index->edge_tree.rebuild(index->edge_segments.cbegin(), index->edge_segments.cend());

    index->annuli.reserve(static_cast<std::size_t>(annuli.rows()));
    index->annulus_boxes.reserve(static_cast<std::size_t>(annuli.rows()));
    for (Eigen::Index row = 0; row < annuli.rows(); ++row) {
        const double cx = annuli(row, 0);
        const double cy = annuli(row, 1);
        const double rho = annuli(row, 2);
        if (!std::isfinite(cx) || !std::isfinite(cy) || !std::isfinite(rho)) {
            throw NonFiniteAnnulusViewInputError(
                "annulus view centres and guide radii must be finite");
        }
        if (rho < 0.0) {
            throw InvalidAnnulusViewRadiiError(
                "annulus view guide radii must be nonnegative");
        }
        const FT guide(rho);
        const FT outer = guide + index->tool_radius;
        const bool has_inner =
            CGAL::compare(guide, index->tool_radius) == CGAL::LARGER;
        const FT inner = has_inner ? FT(guide - index->tool_radius) : FT(0);
        const EPoint centre(cx, cy);
        index->annuli.push_back(
            {centre, inner * inner, outer * outer, has_inner, cx, cy,
             CGAL::Interval_nt<false>(CGAL::to_interval(inner)),
             CGAL::Interval_nt<false>(CGAL::to_interval(outer))});
        const std::size_t id = static_cast<std::size_t>(row);
        index->annulus_boxes.push_back({id, rim_search_box(centre, outer)});
    }
    index->annulus_tree.rebuild(index->annulus_boxes.cbegin(), index->annulus_boxes.cend());
    return AnnulusView2(std::unique_ptr<const Index>(std::move(index)));
}

AnnulusView2::AnnulusView2(std::unique_ptr<const Index> index)
    : index_(std::move(index))
{
}

AnnulusView2::AnnulusView2(AnnulusView2&&) noexcept = default;
AnnulusView2& AnnulusView2::operator=(AnnulusView2&&) noexcept = default;
AnnulusView2::~AnnulusView2() = default;

std::size_t AnnulusView2::size() const
{
    return index_->annuli.size();
}

EngagementSample AnnulusView2::station(
    AnnulusViewRange2 range, double cx, double cy,
    double cap_chord_ratio, double gap_close_ratio) const
{
    validate_range(range, size());
    validate_engagement_input_binary64(
        cx, cy, index_->tool_radius_binary64, cap_chord_ratio, gap_close_ratio);
    const EPoint centre{FT(cx), FT(cy)};
    return engagement_sample_at_station(
        index_->classify_static(centre, range, FT(cap_chord_ratio), FT(gap_close_ratio)),
        centre);
}

std::vector<EngagementSample> AnnulusView2::circle_prefixes(
    AnnulusViewRange2 range, double guide_cx, double guide_cy,
    double guide_radius, double entry_x, double entry_y, bool clockwise,
    const std::vector<double>& parameters, double cap_chord_ratio,
    double gap_close_ratio) const
{
    validate_range(range, size());
    const CirclePrefixStations2 placed = circle_prefix_stations(
        guide_cx, guide_cy, guide_radius, entry_x, entry_y, clockwise, parameters,
        index_->tool_radius_binary64, cap_chord_ratio, gap_close_ratio);
    const FT ratio(cap_chord_ratio);
    const FT gap(gap_close_ratio);
    std::vector<EngagementSample> samples;
    samples.reserve(placed.stations.size());
    for (const ArcPrefixStation2& at : placed.stations) {
        samples.push_back(engagement_sample_at_station(
            clip_arc_swept_prefix_station_exact(
                index_->classify_static(at.centre, range, ratio, gap),
                placed.guide, placed.start, at.centre, placed.clockwise, at.extent,
                index_->tool_radius, ratio, gap),
            at.centre));
    }
    return samples;
}

bool AnnulusView2::circle_exceeded(
    AnnulusViewRange2 range, double guide_cx, double guide_cy,
    double guide_radius, double entry_x, double entry_y, bool clockwise,
    const std::vector<double>& parameters, double cap_chord_ratio,
    double gap_close_ratio) const
{
    validate_range(range, size());
    const CirclePrefixStations2 placed = circle_prefix_stations(
        guide_cx, guide_cy, guide_radius, entry_x, entry_y, clockwise, parameters,
        index_->tool_radius_binary64, cap_chord_ratio, gap_close_ratio);
    const FT ratio(cap_chord_ratio);
    const FT gap(gap_close_ratio);
    for (const std::size_t probe : station_probe_order(placed.stations.size())) {
        const ArcPrefixStation2& at = placed.stations[probe];
        const AuditExactStationClassification2 base =
            index_->classify_static(at.centre, range, ratio, gap);
        // CLIP MONOTONICITY (swept_prefix_engagement_2.h): a static WITHIN_CAP is the
        // in-process verdict.
        if (base.disposition() != AuditExactStationDisposition2::CAP_EXCEEDED) continue;
        if (clip_arc_swept_prefix_station_exact(
                base, placed.guide, placed.start, at.centre, placed.clockwise,
                at.extent, index_->tool_radius, ratio, gap)
                .disposition()
            == AuditExactStationDisposition2::CAP_EXCEEDED) {
            return true;
        }
    }
    return false;
}

StationAgreement2 AnnulusView2::station_agrees_with_stock(
    const Stock2& stock, AnnulusViewRange2 range, double cx, double cy,
    double cap_chord_ratio, double gap_close_ratio) const
{
    validate_range(range, size());
    const EPoint centre{FT(cx), FT(cy)};
    const FT ratio(cap_chord_ratio);
    const FT gap(gap_close_ratio);
    return classifications_agree(
        index_->classify_static(centre, range, ratio, gap),
        classify_audit_unguarded_station_exact(
            stock, centre, index_->tool_radius, ratio, gap));
}

std::vector<StationAgreement2> AnnulusView2::circle_agrees_with_stock(
    const Stock2& stock, AnnulusViewRange2 range, double guide_cx,
    double guide_cy, double guide_radius, double entry_x, double entry_y,
    bool clockwise, const std::vector<double>& parameters,
    double cap_chord_ratio, double gap_close_ratio) const
{
    validate_range(range, size());
    const CirclePrefixStations2 placed = circle_prefix_stations(
        guide_cx, guide_cy, guide_radius, entry_x, entry_y, clockwise, parameters,
        index_->tool_radius_binary64, cap_chord_ratio, gap_close_ratio);
    const FT ratio(cap_chord_ratio);
    const FT gap(gap_close_ratio);
    std::vector<StationAgreement2> agreements;
    agreements.reserve(placed.stations.size());
    for (const ArcPrefixStation2& at : placed.stations) {
        agreements.push_back(classifications_agree(
            clip_arc_swept_prefix_station_exact(
                index_->classify_static(at.centre, range, ratio, gap),
                placed.guide, placed.start, at.centre, placed.clockwise, at.extent,
                index_->tool_radius, ratio, gap),
            classify_arc_swept_prefix_station_exact(
                stock, placed.guide, placed.start, at.centre, placed.clockwise,
                at.extent, index_->tool_radius, ratio, gap)));
    }
    return agreements;
}
