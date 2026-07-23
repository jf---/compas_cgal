#include "toolpath.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <tuple>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes;
typedef CGAL::Straight_skeleton_2<K> Ss;

typedef std::shared_ptr<Ss> SsPtr;

namespace
{

typedef K::Vector_2 Vector_2;
typedef K::Segment_2 Segment_2;
typedef K::Circle_2 Circle_2;
typedef K::Direction_2 Direction_2;

// ---------------------------------------------------------------------------
// Named numerical floors and tolerances.
//
// Geometric *decisions* use exact CGAL predicates.  The constants below guard
// inherently inexact double *constructions* (normalization, tessellation) or
// define output granularity; each states its unit and derivation.
// ---------------------------------------------------------------------------

// Fraction of the tool radius below which two trochoid stations are merged.
// Below this spacing a separate circle removes no measurable material and the
// external-tangent unit normal (denominator = spacing) loses precision.
constexpr double STATION_MERGE_FRACTION = 1e-6;

// Floor for 1/length normalizations of tangent/direction vectors (absolute,
// model units).  Pure division guard, not a geometric tolerance.
constexpr double DIRECTION_NORM_FLOOR = 1e-12;

// Consecutive tessellated output points closer than this merge (absolute,
// model units).  Output granularity only; never used in decisions upstream.
constexpr double OUTPUT_DEDUP_TOL = 1e-9;

// Arrival/departure points on a circle closer than this squared distance skip
// the repositioning arc (= OUTPUT_DEDUP_TOL^2; same granularity rationale).
constexpr double REPOSITION_GATE_SQ = OUTPUT_DEDUP_TOL * OUTPUT_DEDUP_TOL;

// Bridge/lead certification refinement: rounds of midpoint insertion before a
// run is split at the offending gap.  Spacing halves per round, so 8 rounds
// resolve clearance dips down to spacing/256.
constexpr int MAX_BRIDGE_REFINE_ROUNDS = 8;

// Lead lines that gouge are halved at most this many times before being
// dropped (lead length resolution = requested/16).
constexpr int MAX_LEAD_HALVINGS = 4;

inline double approx_length(const Vector_2& v) {
    return std::sqrt(std::max(0.0, CGAL::to_double(v.squared_length())));
}
inline double approx_distance(const Point_2& a, const Point_2& b) {
    return std::sqrt(std::max(0.0, CGAL::to_double(CGAL::squared_distance(a, b))));
}
inline double approx_radius(const Circle_2& c) {
    return std::sqrt(std::max(0.0, CGAL::to_double(c.squared_radius())));
}

struct TrochoidArc {
    Circle_2 circle;   // center + squared_radius + orientation; degenerate = line
    Point_2  start;
    Point_2  end;

    bool is_line() const { return circle.is_degenerate(); }
    bool is_clockwise() const { return circle.orientation() == CGAL::CLOCKWISE; }

    Segment_2 as_segment() const { return Segment_2(start, end); }

    Vector_2 start_tangent() const
    {
        if (is_line()) return end - start;
        return (start - circle.center()).perpendicular(circle.orientation());
    }

    Vector_2 end_tangent() const
    {
        if (is_line()) return end - start;
        return (end - circle.center()).perpendicular(circle.orientation());
    }

    TrochoidArc reversed() const
    {
        return TrochoidArc{is_line() ? circle : circle.opposite(), end, start};
    }

    double sweep() const
    {
        if (is_line()) return 0.0;
        if (start == end) return 2.0 * std::numbers::pi;  // full circle

        const Point_2& c = circle.center();
        const double ux = CGAL::to_double(start.x() - c.x());
        const double uy = CGAL::to_double(start.y() - c.y());
        const double vx = CGAL::to_double(end.x() - c.x());
        const double vy = CGAL::to_double(end.y() - c.y());

        const double cross = ux * vy - uy * vx;
        const double dot_val = ux * vx + uy * vy;
        double ccw = std::atan2(cross, dot_val);
        if (ccw < 0.0) ccw += 2.0 * std::numbers::pi;
        return is_clockwise() ? 2.0 * std::numbers::pi - ccw : ccw;
    }

    double radius() const { return approx_radius(circle); }

    static TrochoidArc make_line(const Point_2& s, const Point_2& e)
    {
        return TrochoidArc{Circle_2(s), s, e};
    }

    static TrochoidArc make_arc(const Point_2& center, double r, const Point_2& s, const Point_2& e, bool cw)
    {
        return TrochoidArc{Circle_2(center, r * r, cw ? CGAL::CLOCKWISE : CGAL::COUNTERCLOCKWISE), s, e};
    }

    static TrochoidArc make_circle(const Point_2& center, double r, const Point_2& on_circle, CGAL::Orientation ori)
    {
        return TrochoidArc{Circle_2(center, r * r, ori), on_circle, on_circle};
    }

};

// ---------------------------------------------------------------------------
// Boundary: outer polygon plus optional holes, flattened to edge segments for
// exact clearance queries.
// ---------------------------------------------------------------------------

struct Boundary {
    std::vector<Segment_2> edges;

    void add_polygon(const Polygon_2& polygon)
    {
        for (auto edge_iter = polygon.edges_begin(); edge_iter != polygon.edges_end(); ++edge_iter) {
            edges.push_back(*edge_iter);
        }
    }

    K::FT squared_distance_to(const Point_2& point) const
    {
        K::FT best = CGAL::squared_distance(point, edges.front());
        for (std::size_t i = 1; i < edges.size(); ++i) {
            const K::FT sq = CGAL::squared_distance(point, edges[i]);
            if (sq < best) best = sq;
        }
        return best;
    }

    K::FT squared_distance_to(const Segment_2& segment) const
    {
        K::FT best = CGAL::squared_distance(segment, edges.front());
        for (std::size_t i = 1; i < edges.size(); ++i) {
            const K::FT sq = CGAL::squared_distance(segment, edges[i]);
            if (sq < best) best = sq;
        }
        return best;
    }

    double distance_to(const Point_2& point) const
    {
        return std::sqrt(std::max(0.0, CGAL::to_double(squared_distance_to(point))));
    }

    // Certified gouge-free test: every point of *segment* keeps at least
    // *clearance* to the boundary.
    bool segment_clear(const Segment_2& segment, double clearance) const
    {
        return !(squared_distance_to(segment) < K::FT(clearance) * K::FT(clearance));
    }
};

// ---------------------------------------------------------------------------
// Stations: cutter-center samples along a skeleton edge, each carrying its
// EXACT boundary clearance.  No interpolation of clearance anywhere — this is
// the structural fix for the reflex-vertex gouge.
// ---------------------------------------------------------------------------

struct Station {
    Point_2 center;
    double clearance;  // exact distance to boundary at *center*
    double radius;     // trochoid radius derived from *clearance*
};

struct RadiusModel {
    double tool_radius;
    double radial_clearance;
    double mat_scale;
    double min_trochoid_radius;
    double max_trochoid_radius;

    double min_centerline() const { return tool_radius + radial_clearance; }

    double radius_from_clearance(double clearance) const
    {
        const double available = mat_scale * std::max(0.0, clearance - tool_radius - radial_clearance);
        if (available <= 0.0) return 0.0;
        double r = std::min(max_trochoid_radius, available);
        if (min_trochoid_radius > 0.0) {
            r = std::min(available, std::max(r, min_trochoid_radius));
        }
        return r;
    }
};

Station make_station(const Point_2& center, const Boundary& boundary, const RadiusModel& model)
{
    const double clearance = boundary.distance_to(center);
    return Station{center, clearance, model.radius_from_clearance(clearance)};
}

// Split edge samples into maximal runs of stations where the cutter fits.
std::vector<std::vector<Station>>
valid_runs(const std::vector<Station>& stations, const RadiusModel& model)
{
    std::vector<std::vector<Station>> runs;
    std::vector<Station> current;
    for (const auto& s : stations) {
        if (s.clearance >= model.min_centerline()) {
            current.push_back(s);
        } else if (!current.empty()) {
            runs.push_back(std::move(current));
            current.clear();
        }
    }
    if (!current.empty()) runs.push_back(std::move(current));
    return runs;
}

// Enforce the engagement bound: per cycle, center advance plus radius growth
// must not exceed *stepover* (the crescent of new material cut by the next
// circle has radial width <= spacing + max(0, r_next - r_prev)).  Clearance is
// 1-Lipschitz so refinement by midpoint insertion converges geometrically.
void refine_engagement(std::vector<Station>& run, const Boundary& boundary,
                       const RadiusModel& model, double stepover)
{
    for (int round = 0; round < MAX_BRIDGE_REFINE_ROUNDS; ++round) {
        bool inserted = false;
        std::vector<Station> refined;
        refined.reserve(run.size() * 2);
        for (std::size_t i = 0; i < run.size(); ++i) {
            refined.push_back(run[i]);
            if (i + 1 >= run.size()) continue;
            const double spacing = approx_distance(run[i].center, run[i + 1].center);
            const double growth = std::max(0.0, run[i + 1].radius - run[i].radius);
            const double merge_floor = STATION_MERGE_FRACTION * model.tool_radius;
            if (spacing + growth > stepover && spacing > 2.0 * merge_floor) {
                const Point_2 mid = CGAL::barycenter(run[i].center, 0.5, run[i + 1].center, 0.5);
                Station s = make_station(mid, boundary, model);
                if (s.clearance >= model.min_centerline()) {
                    refined.push_back(s);
                    inserted = true;
                }
            }
        }
        run.swap(refined);
        if (!inserted) return;
    }
}

// Drop stations that add nothing: near-coincident neighbours (keep the larger
// radius) and disks fully contained in a neighbour's disk.  Guarantees the
// external-tangent precondition |r1 - r0| < distance for consecutive stations,
// which removes the former clamp-induced coordinate blow-up structurally.
void dedup_stations(std::vector<Station>& run, const RadiusModel& model)
{
    const double merge_floor = STATION_MERGE_FRACTION * model.tool_radius;
    const K::FT merge_floor_sq(merge_floor * merge_floor);

    bool changed = true;
    while (changed && run.size() > 1) {
        changed = false;
        std::vector<Station> kept;
        kept.reserve(run.size());
        kept.push_back(run.front());
        for (std::size_t i = 1; i < run.size(); ++i) {
            Station& prev = kept.back();
            const Station& cur = run[i];
            // Near-coincident: keep the larger disk.
            if (CGAL::compare_squared_distance(prev.center, cur.center, merge_floor_sq) != CGAL::LARGER) {
                if (cur.radius > prev.radius) prev = cur;
                changed = true;
                continue;
            }
            // Containment with a conditioning margin: require
            // distance > |dr| + merge_floor so the external-tangent normal
            // (denominator = distance - |dr|) stays well conditioned in double
            // arithmetic downstream.
            const double dr_margin = std::abs(cur.radius - prev.radius) + merge_floor;
            const K::FT dr_margin_sq(dr_margin * dr_margin);
            if (CGAL::compare_squared_distance(prev.center, cur.center, dr_margin_sq) != CGAL::LARGER) {
                // One disk (nearly) contains the other: keep the larger.
                if (cur.radius > prev.radius) prev = cur;
                changed = true;
                continue;
            }
            kept.push_back(cur);
        }
        run.swap(kept);
    }
}

// Choose the tangent side once per chain, from the run direction and milling
// convention.  Single source of truth shared by chain assembly and bridge
// certification.
bool
external_tangent(const Circle_2& c0, const Circle_2& c1, const Vector_2& edge_direction,
                 bool climb_milling, Segment_2& tangent)
{
    if (c0.center() == c1.center()) return false;

    const Vector_2 d = c1.center() - c0.center();
    const double length = approx_length(d);
    const double r0 = approx_radius(c0);
    const double r1 = approx_radius(c1);
    const double delta = r1 - r0;
    if (!(std::abs(delta) < length)) {
        // Structurally unreachable after dedup_stations(); fail loud, never
        // fabricate geometry (the old clamp here caused 6e6-unit excursions).
        throw std::logic_error("external_tangent: contained circles reached tangent construction");
    }

    const double ux = CGAL::to_double(d.x()) / length;
    const double uy = CGAL::to_double(d.y()) / length;
    const double vx = -uy;
    const double vy = ux;
    const double m = -delta / length;
    const double h = std::sqrt(std::max(0.0, 1.0 - m * m));

    const double n1x = m * ux + h * vx;
    const double n1y = m * uy + h * vy;
    const double n2x = m * ux - h * vx;
    const double n2y = m * uy - h * vy;

    const double c0x = CGAL::to_double(c0.center().x());
    const double c0y = CGAL::to_double(c0.center().y());
    const double c1x = CGAL::to_double(c1.center().x());
    const double c1y = CGAL::to_double(c1.center().y());

    const Segment_2 ta(Point_2(c0x + r0 * n1x, c0y + r0 * n1y), Point_2(c1x + r1 * n1x, c1y + r1 * n1y));
    const Segment_2 tb(Point_2(c0x + r0 * n2x, c0y + r0 * n2y), Point_2(c1x + r1 * n2x, c1y + r1 * n2y));

    const Point_2& ci = c0.center();
    const auto orient_a = CGAL::orientation(ci, ci + edge_direction, ta.source());
    tangent = climb_milling
        ? (orient_a == CGAL::LEFT_TURN ? ta : tb)
        : (orient_a == CGAL::RIGHT_TURN ? ta : tb);
    return true;
}

std::vector<TrochoidArc>
trochoid_chain_from_stations(
    const std::vector<Station>& run,
    const Vector_2& edge_direction,
    bool climb_milling)
{
    if (run.size() < 2) return {};

    // Climb milling: tool rotation matches feed direction.
    const CGAL::Orientation arc_ori = climb_milling ? CGAL::CLOCKWISE : CGAL::COUNTERCLOCKWISE;
    const bool cw = (arc_ori == CGAL::CLOCKWISE);

    std::vector<TrochoidArc> chain;
    chain.reserve(run.size() * 3);

    Point_2 prev_arrival;
    bool has_prev = false;

    for (std::size_t i = 0; i + 1 < run.size(); ++i) {
        const Circle_2 ci(run[i].center, run[i].radius * run[i].radius);
        const Circle_2 cj(run[i + 1].center, run[i + 1].radius * run[i + 1].radius);

        Segment_2 tangent;
        if (!external_tangent(ci, cj, edge_direction, climb_milling, tangent)) {
            continue;
        }

        if (!ci.is_degenerate()) {
            const double ri = run[i].radius;
            if (has_prev &&
                CGAL::compare_squared_distance(prev_arrival, tangent.source(), K::FT(REPOSITION_GATE_SQ)) == CGAL::LARGER) {
                // Varying radius: arrival != departure.  Full circle at nominal
                // winding for complete material removal, then short repositioning
                // arc (also nominal winding) to reach the tangent departure point.
                chain.push_back(TrochoidArc::make_circle(run[i].center, ri, prev_arrival, arc_ori));
                chain.push_back(TrochoidArc::make_arc(run[i].center, ri, prev_arrival, tangent.source(), cw));
            } else {
                // First circle or constant radius: full 360 degrees.
                const Point_2& cp = has_prev ? prev_arrival : tangent.source();
                chain.push_back(TrochoidArc::make_circle(run[i].center, ri, cp, arc_ori));
            }
        }

        // Tangent line to next circle (skip if zero-length / overlapping).
        if (tangent.source() != tangent.target()) {
            chain.push_back(TrochoidArc::make_line(tangent.source(), tangent.target()));
        }

        prev_arrival = tangent.target();
        has_prev = true;
    }

    // Full circle on the last circle, at the last tangent arrival point.
    const Station& last = run.back();
    if (has_prev && last.radius > 0.0) {
        chain.push_back(TrochoidArc::make_circle(last.center, last.radius, prev_arrival, arc_ori));
    }

    return chain;
}

// Certify every bridge (line primitive) of a station run against the boundary
// at tool-radius clearance.  Gouging bridges trigger midpoint-station
// insertion; a persistent violation splits the run at that gap (the corridor
// between the stations is genuinely unreachable).
std::vector<std::vector<Station>>
certify_bridges(std::vector<Station> run, const Boundary& boundary,
                const RadiusModel& model, const Vector_2& edge_direction,
                bool climb_milling)
{
    const double merge_floor = STATION_MERGE_FRACTION * model.tool_radius;

    for (int round = 0; round < MAX_BRIDGE_REFINE_ROUNDS; ++round) {
        dedup_stations(run, model);
        if (run.size() < 2) return {std::move(run)};

        bool inserted = false;
        std::vector<Station> refined;
        refined.reserve(run.size() * 2);

        for (std::size_t i = 0; i + 1 < run.size(); ++i) {
            refined.push_back(run[i]);
            const Circle_2 ci(run[i].center, run[i].radius * run[i].radius);
            const Circle_2 cj(run[i + 1].center, run[i + 1].radius * run[i + 1].radius);
            Segment_2 tangent;
            if (!external_tangent(ci, cj, edge_direction, climb_milling, tangent)) continue;
            if (boundary.segment_clear(tangent, model.tool_radius)) continue;

            const double spacing = approx_distance(run[i].center, run[i + 1].center);
            if (spacing > 2.0 * merge_floor) {
                const Point_2 mid = CGAL::barycenter(run[i].center, 0.5, run[i + 1].center, 0.5);
                Station s = make_station(mid, boundary, model);
                if (s.clearance >= model.min_centerline()) {
                    refined.push_back(s);
                    inserted = true;
                    continue;
                }
            }
            // Unreachable corridor between stations: split the run here.
            refined.push_back(Station{run[i].center, -1.0, -1.0});  // split marker
        }
        refined.push_back(run.back());
        run.swap(refined);
        if (!inserted) break;
    }

    // Partition at split markers (clearance < 0).
    std::vector<std::vector<Station>> parts;
    std::vector<Station> current;
    for (const auto& s : run) {
        if (s.clearance < 0.0) {
            if (current.size() >= 2) parts.push_back(std::move(current));
            current.clear();
        } else {
            current.push_back(s);
        }
    }
    if (current.size() >= 2) parts.push_back(std::move(current));
    return parts;
}

std::vector<Point_2>
tessellate_chain(const std::vector<TrochoidArc>& chain, int samples_per_arc)
{
    std::vector<Point_2> points;
    points.reserve(chain.size() * samples_per_arc);

    for (const auto& arc : chain) {
        if (arc.is_line()) {
            points.push_back(arc.start);
        } else {
            const double r = arc.radius();
            const Point_2& c = arc.circle.center();
            const double cx = CGAL::to_double(c.x());
            const double cy = CGAL::to_double(c.y());

            const double sx = CGAL::to_double(arc.start.x()) - cx;
            const double sy = CGAL::to_double(arc.start.y()) - cy;
            const double start_angle = std::atan2(sy, sx);

            const double sw = arc.sweep();
            const double signed_sweep = arc.is_clockwise() ? -sw : sw;

            const int n = std::max(2, samples_per_arc);
            for (int i = 0; i < n; ++i) {
                const double t = static_cast<double>(i) / static_cast<double>(n - 1);
                const double theta = start_angle + signed_sweep * t;
                points.emplace_back(cx + r * std::cos(theta), cy + r * std::sin(theta));
            }
        }
    }

    // Ensure we end at the last arc's endpoint
    if (!chain.empty()) {
        points.push_back(chain.back().end);
    }

    return points;
}

Polygon_2
data_to_polygon(Eigen::Ref<const compas::RowMatrixXd> vertices)
{
    if (vertices.rows() < 3) {
        throw std::invalid_argument("Expected at least three polygon vertices.");
    }

    Polygon_2 polygon;
    polygon.reserve(vertices.rows());

    for (int i = 0; i < vertices.rows(); ++i) {
        polygon.push_back(Point_2(vertices(i, 0), vertices(i, 1)));
    }

    if (!polygon.is_simple()) {
        throw std::invalid_argument(
            "Polygon boundary must be simple (no self-intersections, no repeated vertices).");
    }

    if (polygon.is_clockwise_oriented()) {
        polygon.reverse_orientation();
    }

    return polygon;
}

// Assemble outer boundary + holes: validates simplicity, hole containment and
// pairwise hole disjointness (CGAL straight-skeleton preconditions), and
// returns both the CGAL input and the flattened clearance boundary.
std::tuple<Polygon_with_holes, Boundary>
assemble_domain(
    Eigen::Ref<const compas::RowMatrixXd> vertices,
    const std::vector<compas::RowMatrixXd>& holes)
{
    Polygon_2 outer = data_to_polygon(vertices);

    Boundary boundary;
    boundary.add_polygon(outer);

    Polygon_with_holes domain(outer);
    std::vector<Polygon_2> hole_polygons;
    hole_polygons.reserve(holes.size());

    for (const auto& hole_data : holes) {
        Polygon_2 hole = data_to_polygon(hole_data);  // validates simplicity, makes CCW
        for (auto vertex_iter = hole.vertices_begin(); vertex_iter != hole.vertices_end(); ++vertex_iter) {
            if (outer.bounded_side(*vertex_iter) != CGAL::ON_BOUNDED_SIDE) {
                throw std::invalid_argument("Hole polygon must lie strictly inside the outer boundary.");
            }
        }
        for (const auto& other : hole_polygons) {
            for (auto ea = hole.edges_begin(); ea != hole.edges_end(); ++ea) {
                for (auto eb = other.edges_begin(); eb != other.edges_end(); ++eb) {
                    if (CGAL::do_intersect(*ea, *eb)) {
                        throw std::invalid_argument("Hole polygons must be pairwise disjoint.");
                    }
                }
            }
        }
        boundary.add_polygon(hole);
        hole_polygons.push_back(hole);
        hole.reverse_orientation();  // CGAL interior skeleton expects CW holes
        domain.add_hole(hole);
    }

    return std::make_tuple(domain, boundary);
}

SsPtr
build_skeleton(const Polygon_with_holes& domain)
{
    SsPtr skeleton = CGAL::create_interior_straight_skeleton_2(domain);
    if (!skeleton) {
        throw std::invalid_argument(
            "Straight skeleton construction failed: polygon is degenerate or violates preconditions.");
    }
    return skeleton;
}

void
validate_toolpath_params(
    double tool_diameter, double stepover, double pitch,
    double min_trochoid_radius, double max_trochoid_radius,
    int samples_per_cycle, int max_passes)
{
    if (tool_diameter <= 0.0) throw std::invalid_argument("tool_diameter should be positive.");
    if (stepover <= 0.0) throw std::invalid_argument("stepover should be positive.");
    if (pitch <= 0.0) throw std::invalid_argument("pitch should be positive.");
    if (min_trochoid_radius < 0.0) throw std::invalid_argument("min_trochoid_radius should be >= 0.");
    if (max_trochoid_radius < 0.0) throw std::invalid_argument("max_trochoid_radius should be >= 0.");
    if (samples_per_cycle < 4) throw std::invalid_argument("samples_per_cycle should be at least 4.");
    if (max_passes <= 0) throw std::invalid_argument("max_passes should be positive.");
}

// Walk the straight skeleton and return certified trochoid chains per edge,
// plus the number of qualifying edges skipped once max_passes was reached.
std::tuple<std::vector<std::vector<TrochoidArc>>, int>
skeleton_edge_chains(
    const Boundary& boundary,
    const SsPtr& skeleton,
    const RadiusModel& model,
    double pitch,
    double stepover,
    bool climb_milling,
    int max_passes)
{
    const double station_spacing = std::min(pitch, stepover);

    std::vector<std::vector<TrochoidArc>> chains;
    chains.reserve(static_cast<std::size_t>(max_passes));
    int skipped_edges = 0;

    for (auto edge_iter = skeleton->halfedges_begin();
         edge_iter != skeleton->halfedges_end(); ++edge_iter) {
        if (!(edge_iter->is_bisector() || edge_iter->is_inner_bisector())) continue;

        const auto v0 = edge_iter->vertex();
        const auto v1 = edge_iter->opposite()->vertex();
        // Deterministic halfedge-pair canonicalization by vertex id (allocation
        // addresses are not reproducible across runs/platforms).
        if (v0->id() >= v1->id()) continue;

        Point_2 p0(v0->point().x(), v0->point().y());
        Point_2 p1(v1->point().x(), v1->point().y());
        if (p0 == p1) continue;

        // Exact early rejection: both endpoints too close to the boundary.
        const double min_cd = model.min_centerline();
        const K::FT min_cd_sq(min_cd * min_cd);
        if (boundary.squared_distance_to(p0) < min_cd_sq &&
            boundary.squared_distance_to(p1) < min_cd_sq) {
            continue;
        }

        // Canonical direction: narrower end first (exact compare, id tiebreak
        // handled above by the halfedge selection).
        if (boundary.squared_distance_to(p0) > boundary.squared_distance_to(p1)) {
            std::swap(p0, p1);
        }
        const Vector_2 edge_dir = p1 - p0;

        // Sample stations at engagement-bounded spacing with EXACT clearance.
        const double length = approx_distance(p0, p1);
        const int n_spans = std::max(1, static_cast<int>(std::ceil(length / station_spacing)));
        std::vector<Station> stations;
        stations.reserve(n_spans + 1);
        for (int i = 0; i <= n_spans; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(n_spans);
            stations.push_back(make_station(CGAL::barycenter(p0, 1.0 - t, p1, t), boundary, model));
        }

        for (auto& run : valid_runs(stations, model)) {
            if (run.size() < 2) continue;
            refine_engagement(run, boundary, model, stepover);
            for (auto& part : certify_bridges(std::move(run), boundary, model, edge_dir, climb_milling)) {
                dedup_stations(part, model);
                if (part.size() < 2) continue;
                auto chain = trochoid_chain_from_stations(part, edge_dir, climb_milling);
                if (chain.empty()) continue;
                if (static_cast<int>(chains.size()) >= max_passes) {
                    skipped_edges += 1;
                } else {
                    chains.push_back(std::move(chain));
                }
            }
        }
    }

    return std::make_tuple(std::move(chains), skipped_edges);
}

void
deduplicate_consecutive_points(std::vector<Point_2>& points, double tol)
{
    if (points.size() < 2) {
        return;
    }
    const K::FT tol_sq(tol * tol);
    auto last = std::unique(points.begin(), points.end(),
        [&tol_sq](const Point_2& a, const Point_2& b) {
            return CGAL::compare_squared_distance(a, b, tol_sq) != CGAL::LARGER;
        });
    points.erase(last, points.end());
}

compas::RowMatrixXd
points_to_matrix(const std::vector<Point_2>& points)
{
    compas::RowMatrixXd matrix(points.size(), 3);
    for (std::size_t i = 0; i < points.size(); ++i) {
        matrix(i, 0) = points[i].x();
        matrix(i, 1) = points[i].y();
        matrix(i, 2) = 0.0;
    }
    return matrix;
}

} // namespace

std::tuple<compas::RowMatrixXd, compas::RowMatrixXd>
pmp_polygon_skeleton_clearance(
    Eigen::Ref<const compas::RowMatrixXd> vertices,
    const std::vector<compas::RowMatrixXd>& holes)
{
    auto [domain, boundary] = assemble_domain(vertices, holes);
    SsPtr skeleton = build_skeleton(domain);

    std::vector<Point_2> points;
    std::vector<double> radii;
    points.reserve(skeleton->size_of_vertices());
    radii.reserve(skeleton->size_of_vertices());

    for (auto vertex_iter = skeleton->vertices_begin(); vertex_iter != skeleton->vertices_end(); ++vertex_iter) {
        const auto& point = vertex_iter->point();
        Point_2 point_xy(point.x(), point.y());
        const K::FT sq_dist = boundary.squared_distance_to(point_xy);

        // Keep only interior skeleton samples (exact comparison against zero).
        if (sq_dist > K::FT(0)) {
            points.push_back(point_xy);
            radii.push_back(std::sqrt(std::max(0.0, CGAL::to_double(sq_dist))));
        }
    }

    compas::RowMatrixXd points_matrix(points.size(), 3);
    compas::RowMatrixXd radii_matrix(radii.size(), 1);
    for (std::size_t i = 0; i < points.size(); ++i) {
        points_matrix(i, 0) = points[i].x();
        points_matrix(i, 1) = points[i].y();
        points_matrix(i, 2) = 0.0;
        radii_matrix(i, 0) = radii[i];
    }
    return std::make_tuple(points_matrix, radii_matrix);
}

std::tuple<std::vector<compas::RowMatrixXd>, int>
pmp_trochoidal_mat_toolpath(
    Eigen::Ref<const compas::RowMatrixXd> vertices,
    const std::vector<compas::RowMatrixXd>& holes,
    double tool_diameter,
    double stepover,
    double pitch,
    double min_trochoid_radius,
    double max_trochoid_radius,
    double mat_scale,
    double radial_clearance,
    int samples_per_cycle,
    int max_passes,
    bool climb)
{
    validate_toolpath_params(tool_diameter, stepover, pitch,
        min_trochoid_radius, max_trochoid_radius, samples_per_cycle, max_passes);

    auto [domain, boundary] = assemble_domain(vertices, holes);
    SsPtr skeleton = build_skeleton(domain);

    const RadiusModel model{0.5 * tool_diameter, radial_clearance, mat_scale,
                            min_trochoid_radius, max_trochoid_radius};
    auto [chains, skipped] = skeleton_edge_chains(boundary, skeleton, model,
                                                  pitch, stepover, climb, max_passes);

    std::vector<compas::RowMatrixXd> toolpaths;
    toolpaths.reserve(chains.size());
    for (auto& chain : chains) {
        auto pts = tessellate_chain(chain, samples_per_cycle);
        deduplicate_consecutive_points(pts, OUTPUT_DEDUP_TOL);
        if (pts.size() >= 2) {
            toolpaths.push_back(points_to_matrix(pts));
        }
    }
    return std::make_tuple(std::move(toolpaths), skipped);
}

// ============================================================================
// Linked circular toolpath with path ordering, leads, links, retract/plunge
// ============================================================================

namespace
{

struct ToolpathPrimitive {
    TrochoidArc arc;
    int path_index;
    int operation;  // 0=cut 1=lead_in 2=lead_out 3=link 4=retract 5=plunge
    double z_start;
    double z_end;
};

double
path_length(const std::vector<TrochoidArc>& path)
{
    double length = 0.0;
    for (const auto& a : path) {
        if (a.is_line()) {
            length += approx_distance(a.start, a.end);
        } else {
            length += a.radius() * a.sweep();
        }
    }
    return length;
}

compas::RowMatrixXd
tessellate_operations(const std::vector<ToolpathPrimitive>& ops, double samples_per_radian)
{
    if (ops.empty()) return compas::RowMatrixXd(0, 3);

    std::vector<std::array<double, 3>> pts;
    pts.reserve(ops.size() * 32);

    for (std::size_t oi = 0; oi < ops.size(); ++oi) {
        const auto& op = ops[oi];
        const bool first = pts.empty();

        if (op.arc.is_line()) {
            if (first) {
                pts.push_back({CGAL::to_double(op.arc.start.x()),
                               CGAL::to_double(op.arc.start.y()),
                               op.z_start});
            }
            pts.push_back({CGAL::to_double(op.arc.end.x()),
                           CGAL::to_double(op.arc.end.y()),
                           op.z_end});
        } else {
            const double r = op.arc.radius();
            const double cx = CGAL::to_double(op.arc.circle.center().x());
            const double cy = CGAL::to_double(op.arc.circle.center().y());
            const double sx = CGAL::to_double(op.arc.start.x()) - cx;
            const double sy = CGAL::to_double(op.arc.start.y()) - cy;
            const double start_angle = std::atan2(sy, sx);

            const double sw = op.arc.sweep();
            const double signed_sweep = op.arc.is_clockwise() ? -sw : sw;
            const int n = std::max(2, static_cast<int>(std::ceil(sw * samples_per_radian)));

            // Skip first sample (junction duplicate) unless first operation
            const int i0 = first ? 0 : 1;
            for (int i = i0; i < n - 1; ++i) {
                const double t = static_cast<double>(i) / static_cast<double>(n - 1);
                const double theta = start_angle + signed_sweep * t;
                const double z = op.z_start + t * (op.z_end - op.z_start);
                pts.push_back({cx + r * std::cos(theta),
                               cy + r * std::sin(theta), z});
            }
            // Exact endpoint — avoids cos/sin drift
            pts.push_back({CGAL::to_double(op.arc.end.x()),
                           CGAL::to_double(op.arc.end.y()),
                           op.z_end});
        }
    }

    const int m = static_cast<int>(pts.size());
    compas::RowMatrixXd result(m, 3);
    for (int i = 0; i < m; ++i) {
        result(i, 0) = pts[i][0];
        result(i, 1) = pts[i][1];
        result(i, 2) = pts[i][2];
    }
    return result;
}

// Shrink an engaged lead segment until it is certified gouge-free; empty
// optional means the lead cannot fit and is dropped.
bool
certified_lead(const Point_2& anchor, const Vector_2& direction_unit, double length,
               const Boundary& boundary, double tool_radius, bool outward,
               Point_2& lead_point)
{
    double len = length;
    for (int i = 0; i <= MAX_LEAD_HALVINGS; ++i) {
        const double sign = outward ? 1.0 : -1.0;
        const Point_2 candidate(
            CGAL::to_double(anchor.x()) + sign * len * CGAL::to_double(direction_unit.x()),
            CGAL::to_double(anchor.y()) + sign * len * CGAL::to_double(direction_unit.y()));
        if (boundary.segment_clear(Segment_2(anchor, candidate), tool_radius)) {
            lead_point = candidate;
            return true;
        }
        len *= 0.5;
    }
    return false;
}

} // namespace

std::tuple<compas::RowMatrixXd, compas::RowMatrixXd, compas::RowMatrixXd, compas::RowMatrixXd, compas::RowMatrixXd, compas::RowMatrixXd, compas::RowMatrixXd, compas::RowMatrixXd, int>
pmp_trochoidal_mat_toolpath_circular(
    Eigen::Ref<const compas::RowMatrixXd> vertices,
    const std::vector<compas::RowMatrixXd>& holes,
    double tool_diameter,
    double stepover,
    double pitch,
    double min_trochoid_radius,
    double max_trochoid_radius,
    double mat_scale,
    double radial_clearance,
    int samples_per_cycle,
    int max_passes,
    bool climb,
    double lead_in,
    double lead_out,
    bool link_paths,
    bool optimize_order,
    double cut_z,
    double clearance_z,
    bool has_clearance_z,
    bool retract_at_end,
    double samples_per_radian)
{
    validate_toolpath_params(tool_diameter, stepover, pitch,
        min_trochoid_radius, max_trochoid_radius, samples_per_cycle, max_passes);

    auto [domain, boundary] = assemble_domain(vertices, holes);
    SsPtr skeleton = build_skeleton(domain);

    const double tool_radius = 0.5 * tool_diameter;
    const RadiusModel model{tool_radius, radial_clearance, mat_scale,
                            min_trochoid_radius, max_trochoid_radius};
    auto [chains, skipped] = skeleton_edge_chains(boundary, skeleton, model,
                                                  pitch, stepover, climb, max_passes);

    // Filter chains without arcs, build paths vector
    std::vector<std::vector<TrochoidArc>> paths;
    paths.reserve(chains.size());
    for (auto& chain : chains) {
        bool has_arc = std::any_of(chain.begin(), chain.end(),
            [](const TrochoidArc& a) { return !a.is_line(); });
        if (has_arc && !chain.empty()) {
            paths.push_back(std::move(chain));
        }
    }

    if (paths.empty()) {
        compas::RowMatrixXd empty_meta(0, 4);
        compas::RowMatrixXd empty3(0, 3);
        compas::RowMatrixXd empty1(0, 1);
        return std::make_tuple(empty_meta, empty3, empty3, empty3, empty1, empty3, empty3, empty3, skipped);
    }

    // Greedy nearest-neighbor path ordering (documented heuristic; optimal
    // ordering is a TSP and out of scope).
    if (optimize_order && paths.size() > 1) {
        double max_len = -1.0;
        int start_idx = 0;
        for (int i = 0; i < static_cast<int>(paths.size()); ++i) {
            double len = path_length(paths[i]);
            if (len > max_len) { max_len = len; start_idx = i; }
        }

        std::vector<bool> used(paths.size(), false);
        std::vector<std::vector<TrochoidArc>> ordered;
        ordered.reserve(paths.size());
        ordered.push_back(paths[start_idx]);
        used[start_idx] = true;

        for (std::size_t step = 1; step < paths.size(); ++step) {
            const Point_2& cur_end = ordered.back().back().end;
            int best_idx = -1;

            for (int i = 0; i < static_cast<int>(paths.size()); ++i) {
                if (used[i]) continue;
                if (best_idx < 0 ||
                    CGAL::has_smaller_distance_to_point(cur_end, paths[i].front().start, paths[best_idx].front().start)) {
                    best_idx = i;
                }
            }

            ordered.push_back(paths[best_idx]);
            used[best_idx] = true;
        }
        paths.swap(ordered);
    }

    // Build linked operation stream using ToolpathPrimitive
    const bool use_clearance = has_clearance_z && (clearance_z > cut_z);
    const double safe_z = has_clearance_z ? clearance_z : cut_z;

    std::vector<ToolpathPrimitive> operations;
    std::size_t total_arcs = 0;
    for (const auto& p : paths) total_arcs += p.size();
    operations.reserve(total_arcs + paths.size() * 6);

    auto make_tp_line = [](const Point_2& s, double sz, const Point_2& e, double ez, int op, int pidx) {
        ToolpathPrimitive tp;
        tp.arc = TrochoidArc::make_line(s, e);
        tp.path_index = pidx;
        tp.operation = op;
        tp.z_start = sz;
        tp.z_end = ez;
        return tp;
    };

    Point_2 cur_xy;
    double cur_z = 0.0;
    bool has_current = false;

    for (int pidx = 0; pidx < static_cast<int>(paths.size()); ++pidx) {
        auto& path = paths[pidx];
        if (path.empty()) continue;

        const Point_2& path_start = path.front().start;
        const Point_2& path_end = path.back().end;

        // Certified lead-in start point via start tangent (shrinks to fit,
        // dropped when no gouge-free length exists).
        const Vector_2 st = path.front().start_tangent();
        const double st_len = approx_length(st);
        Point_2 lead_in_pt = path_start;  // default: no lead-in
        bool has_lead_in = false;
        if (lead_in > 0.0 && st_len > DIRECTION_NORM_FLOOR) {
            const double inv = 1.0 / st_len;
            const Vector_2 st_unit(CGAL::to_double(st.x()) * inv, CGAL::to_double(st.y()) * inv);
            has_lead_in = certified_lead(path_start, st_unit, lead_in, boundary,
                                         tool_radius, /*outward=*/false, lead_in_pt);
            if (!has_lead_in) lead_in_pt = path_start;
        }

        // Connect to this path
        if (!has_current) {
            if (use_clearance) {
                operations.push_back(make_tp_line(lead_in_pt, safe_z, lead_in_pt, cut_z, 5 /*plunge*/, pidx));
            }
            cur_xy = lead_in_pt; cur_z = cut_z;
            has_current = true;
        } else if (link_paths) {
            if (use_clearance) {
                // retract
                operations.push_back(make_tp_line(cur_xy, cur_z, cur_xy, safe_z, 4 /*retract*/, pidx));
                cur_z = safe_z;
                // link at clearance
                if (cur_xy != lead_in_pt) {
                    operations.push_back(make_tp_line(cur_xy, safe_z, lead_in_pt, safe_z, 3 /*link*/, pidx));
                }
                // plunge
                operations.push_back(make_tp_line(lead_in_pt, safe_z, lead_in_pt, cut_z, 5 /*plunge*/, pidx));
                cur_xy = lead_in_pt; cur_z = cut_z;
            } else {
                if (cur_xy != lead_in_pt) {
                    // A flat link is an engaged move: it must be certified
                    // against the walls.  Without a clearance plane there is no
                    // safe fallback — fail loud.
                    if (!boundary.segment_clear(Segment_2(cur_xy, lead_in_pt), tool_radius)) {
                        throw std::invalid_argument(
                            "Flat link between paths would gouge the boundary; "
                            "provide clearance_z for safe Z-linking.");
                    }
                    operations.push_back(make_tp_line(cur_xy, cut_z, lead_in_pt, cut_z, 3 /*link*/, pidx));
                    cur_xy = lead_in_pt; cur_z = cut_z;
                }
            }
        } else {
            cur_xy = lead_in_pt; cur_z = cut_z;
        }

        // Lead-in
        if (has_lead_in && lead_in_pt != path_start) {
            operations.push_back(make_tp_line(lead_in_pt, cut_z, path_start, cut_z, 1 /*lead_in*/, pidx));
            cur_xy = path_start; cur_z = cut_z;
        }

        // Cut primitives
        for (const auto& arc : path) {
            ToolpathPrimitive tp;
            tp.arc = arc;
            tp.path_index = pidx;
            tp.operation = 0;  // cut
            tp.z_start = cut_z;
            tp.z_end = cut_z;
            operations.push_back(tp);
            cur_xy = arc.end;
            cur_z = cut_z;
        }

        // Certified lead-out via end tangent
        {
            const Vector_2 et = path.back().end_tangent();
            const double et_len = approx_length(et);
            if (lead_out > 0.0 && et_len > DIRECTION_NORM_FLOOR) {
                const double inv = 1.0 / et_len;
                const Vector_2 et_unit(CGAL::to_double(et.x()) * inv, CGAL::to_double(et.y()) * inv);
                Point_2 lo_end = path_end;
                if (certified_lead(path_end, et_unit, lead_out, boundary,
                                   tool_radius, /*outward=*/true, lo_end) &&
                    lo_end != path_end) {
                    operations.push_back(make_tp_line(path_end, cut_z, lo_end, cut_z, 2 /*lead_out*/, pidx));
                    cur_xy = lo_end; cur_z = cut_z;
                }
            }
        }
    }

    // Final retract
    if (use_clearance && retract_at_end && has_current) {
        if (cur_z != safe_z) {
            operations.push_back(make_tp_line(cur_xy, cur_z, cur_xy, safe_z, 4 /*retract*/, static_cast<int>(paths.size())));
        }
    }

    // Serialize to matrices (meta Nx4)
    const int n = static_cast<int>(operations.size());
    compas::RowMatrixXd meta(n, 4);
    compas::RowMatrixXd starts(n, 3);
    compas::RowMatrixXd ends(n, 3);
    compas::RowMatrixXd centers_out(n, 3);
    compas::RowMatrixXd radii_out(n, 1);
    compas::RowMatrixXd start_tangents(n, 3);
    compas::RowMatrixXd end_tangents(n, 3);
    for (int i = 0; i < n; ++i) {
        const auto& op = operations[i];
        bool is_circle = !op.arc.is_line() && op.arc.start == op.arc.end;

        meta(i, 0) = static_cast<double>(op.path_index);
        meta(i, 1) = op.arc.is_line() ? 0.0 : (is_circle ? 2.0 : 1.0);
        meta(i, 2) = op.arc.is_clockwise() ? 1.0 : 0.0;
        meta(i, 3) = static_cast<double>(op.operation);

        starts(i, 0) = CGAL::to_double(op.arc.start.x());
        starts(i, 1) = CGAL::to_double(op.arc.start.y());
        starts(i, 2) = op.z_start;
        ends(i, 0) = CGAL::to_double(op.arc.end.x());
        ends(i, 1) = CGAL::to_double(op.arc.end.y());
        ends(i, 2) = op.z_end;
        centers_out(i, 0) = CGAL::to_double(op.arc.circle.center().x());
        centers_out(i, 1) = CGAL::to_double(op.arc.circle.center().y());
        centers_out(i, 2) = op.z_start;  // center z matches start z
        radii_out(i, 0) = op.arc.radius();

        // Unit tangent vectors (XY, z=0)
        const Vector_2 st = op.arc.start_tangent();
        const double st_len = approx_length(st);
        if (st_len > DIRECTION_NORM_FLOOR) {
            const double inv = 1.0 / st_len;
            start_tangents(i, 0) = CGAL::to_double(st.x()) * inv;
            start_tangents(i, 1) = CGAL::to_double(st.y()) * inv;
        } else {
            start_tangents(i, 0) = 0.0;
            start_tangents(i, 1) = 0.0;
        }
        start_tangents(i, 2) = 0.0;

        const Vector_2 et = op.arc.end_tangent();
        const double et_len = approx_length(et);
        if (et_len > DIRECTION_NORM_FLOOR) {
            const double inv = 1.0 / et_len;
            end_tangents(i, 0) = CGAL::to_double(et.x()) * inv;
            end_tangents(i, 1) = CGAL::to_double(et.y()) * inv;
        } else {
            end_tangents(i, 0) = 0.0;
            end_tangents(i, 1) = 0.0;
        }
        end_tangents(i, 2) = 0.0;
    }

    auto polyline = tessellate_operations(operations, samples_per_radian);
    return std::make_tuple(meta, starts, ends, centers_out, radii_out, polyline, start_tangents, end_tangents, skipped);
}

NB_MODULE(_toolpath, m)
{
    m.def(
        "polygon_skeleton_clearance",
        &pmp_polygon_skeleton_clearance,
        "Interior straight-skeleton vertices with exact boundary clearance.\n\n"
        "The straight skeleton coincides with the medial axis for convex\n"
        "polygons only; radii are exact clearances at the returned loci.\n\n"
        "Parameters\n"
        "----------\n"
        "vertices : array-like\n"
        "    Matrix of outer boundary vertices (Nx3, float64)\n"
        "holes : list of array-like\n"
        "    Optional hole polygons (each Mx3, float64)\n"
        "\n"
        "Returns\n"
        "-------\n"
        "tuple\n"
        "    - skeleton points (Mx3, float64)\n"
        "    - exact clearance radii (Mx1, float64)",
        "vertices"_a,
        "holes"_a);

    m.def(
        "trochoidal_mat_toolpath",
        &pmp_trochoidal_mat_toolpath,
        "Certified gouge-free trochoidal pocket toolpaths along the straight\n"
        "skeleton, with exact per-station clearance radii.\n\n"
        "Returns (list of Kx3 polylines, number of skeleton edges skipped by\n"
        "max_passes).",
        "vertices"_a,
        "holes"_a,
        "tool_diameter"_a,
        "stepover"_a,
        "pitch"_a,
        "min_trochoid_radius"_a,
        "max_trochoid_radius"_a,
        "mat_scale"_a,
        "radial_clearance"_a,
        "samples_per_cycle"_a,
        "max_passes"_a,
        "climb"_a = true);

    m.def(
        "trochoidal_mat_toolpath_circular",
        &pmp_trochoidal_mat_toolpath_circular,
        "Create linked circular toolpath with path ordering, leads, and links.\n\n"
        "Returns tuple:\n"
        "- meta (Nx4 float): [path_index, type(0=line,1=arc,2=circle), clockwise, operation(0=cut,1=lead_in,2=lead_out,3=link,4=retract,5=plunge)]\n"
        "- start points (Nx3 float)\n"
        "- end points (Nx3 float)\n"
        "- centers (Nx3 float)\n"
        "- radii (Nx1 float)\n"
        "- polyline (Mx3 float): tessellated 3D point sequence\n"
        "- start tangents (Nx3 float): unit tangent at arc/circle start\n"
        "- end tangents (Nx3 float): unit tangent at arc/circle end\n"
        "- skipped (int): skeleton edges dropped once max_passes was reached",
        "vertices"_a,
        "holes"_a,
        "tool_diameter"_a,
        "stepover"_a,
        "pitch"_a,
        "min_trochoid_radius"_a,
        "max_trochoid_radius"_a,
        "mat_scale"_a,
        "radial_clearance"_a,
        "samples_per_cycle"_a,
        "max_passes"_a,
        "climb"_a = true,
        "lead_in"_a = 0.0,
        "lead_out"_a = 0.0,
        "link_paths"_a = true,
        "optimize_order"_a = true,
        "cut_z"_a = 0.0,
        "clearance_z"_a = 0.0,
        "has_clearance_z"_a = false,
        "retract_at_end"_a = true,
        "samples_per_radian"_a = 10.0);
}
