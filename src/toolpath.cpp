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
typedef CGAL::Straight_skeleton_2<K> Ss;

typedef std::shared_ptr<Ss> SsPtr;

namespace
{

typedef K::Vector_2 Vector_2;
typedef K::Segment_2 Segment_2;
typedef K::Circle_2 Circle_2;
typedef K::Direction_2 Direction_2;

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

std::vector<Circle_2>
trochoid_circles(
    const Point_2& p0,
    const Point_2& p1,
    double r0,
    double r1,
    double pitch)
{
    if (p0 == p1) {
        return {};
    }
    const double length = approx_distance(p0, p1);

    const int cycles = std::max(2, static_cast<int>(std::ceil(length / std::max(pitch, 1e-12))));

    std::vector<Circle_2> circles;
    circles.reserve(cycles + 1);
    for (int i = 0; i <= cycles; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(cycles);
        const Point_2 center = CGAL::barycenter(p0, 1.0 - t, p1, t);
        const double radius = std::max(0.0, r0 + t * (r1 - r0));
        circles.emplace_back(center, radius * radius);
    }
    return circles;
}

bool
external_tangents(
    const Circle_2& c0,
    const Circle_2& c1,
    Segment_2& tangent_a,
    Segment_2& tangent_b)
{
    const Vector_2 d = c1.center() - c0.center();
    if (c0.center() == c1.center()) {
        return false;
    }

    const double length = approx_length(d);
    const double r0 = approx_radius(c0);
    const double r1 = approx_radius(c1);

    double delta = r1 - r0;
    if (std::abs(delta) >= length) {
        delta = (delta >= 0.0 ? 1.0 : -1.0) * (length - 1e-9);
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

    tangent_a = Segment_2(
        Point_2(c0x + r0 * n1x, c0y + r0 * n1y),
        Point_2(c1x + r1 * n1x, c1y + r1 * n1y));
    tangent_b = Segment_2(
        Point_2(c0x + r0 * n2x, c0y + r0 * n2y),
        Point_2(c1x + r1 * n2x, c1y + r1 * n2y));
    return true;
}

std::vector<TrochoidArc>
trochoid_chain(
    const std::vector<Circle_2>& circles,
    const Vector_2& edge_direction,
    bool climb_milling)
{
    if (circles.size() < 2) {
        return {};
    }

    // Climb milling: tool rotation matches feed direction.
    const CGAL::Orientation arc_ori = climb_milling ? CGAL::CLOCKWISE : CGAL::COUNTERCLOCKWISE;
    const bool cw = (arc_ori == CGAL::CLOCKWISE);

    std::vector<TrochoidArc> chain;
    chain.reserve(circles.size() * 3);

    Point_2 prev_arrival;
    bool has_prev = false;

    for (std::size_t i = 0; i + 1 < circles.size(); ++i) {
        Segment_2 ta, tb;
        if (!external_tangents(circles[i], circles[i + 1], ta, tb)) {
            continue;
        }

        // Always pick the same tangent side (no alternation).
        const Point_2& ci = circles[i].center();
        const auto orient_a = CGAL::orientation(ci, ci + edge_direction, ta.source());
        const Segment_2& tangent = climb_milling
            ? (orient_a == CGAL::LEFT_TURN ? ta : tb)
            : (orient_a == CGAL::RIGHT_TURN ? ta : tb);

        if (!circles[i].is_degenerate()) {
            const double ri = approx_radius(circles[i]);

            if (has_prev && prev_arrival != tangent.source()) {
                // Varying radius: arrival ≠ departure.  Single arc from arrival
                // to departure, always taking the long way around (>π).
                // The gap side depends on radius gradient direction, so we
                // try the nominal winding first and flip if it takes the short way.
                auto arc = TrochoidArc::make_arc(ci, ri, prev_arrival, tangent.source(), cw);
                if (arc.sweep() < std::numbers::pi) {
                    arc = TrochoidArc::make_arc(ci, ri, prev_arrival, tangent.source(), !cw);
                }
                chain.push_back(arc);
            } else {
                // First circle or constant radius: full 360°.
                const Point_2& cp = has_prev ? prev_arrival : tangent.source();
                chain.push_back(TrochoidArc::make_circle(ci, ri, cp, arc_ori));
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
    if (has_prev && !circles.back().is_degenerate()) {
        const double rlast = approx_radius(circles.back());
        chain.push_back(TrochoidArc::make_circle(circles.back().center(), rlast, prev_arrival, arc_ori));
    }

    return chain;
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

    if (polygon.is_clockwise_oriented()) {
        polygon.reverse_orientation();
    }

    return polygon;
}

K::FT
squared_distance_to_boundary(const Point_2& point, const Polygon_2& boundary)
{
    if (boundary.size() < 2) {
        return K::FT(0);
    }

    K::FT sq_distance = CGAL::squared_distance(point, *boundary.edges_begin());
    for (auto edge_iter = std::next(boundary.edges_begin()); edge_iter != boundary.edges_end(); ++edge_iter) {
        const K::FT sq = CGAL::squared_distance(point, *edge_iter);
        if (sq < sq_distance) {
            sq_distance = sq;
        }
    }
    return sq_distance;
}

double
approx_distance_to_boundary(const Point_2& point, const Polygon_2& boundary)
{
    return std::sqrt(std::max(0.0, CGAL::to_double(squared_distance_to_boundary(point, boundary))));
}

std::tuple<std::vector<Point_2>, std::vector<double>>
polygon_medial_axis_transform_internal(const Polygon_2& polygon)
{
    SsPtr skeleton = CGAL::create_interior_straight_skeleton_2(polygon.vertices_begin(), polygon.vertices_end());
    if (!skeleton) {
        return std::make_tuple(std::vector<Point_2>{}, std::vector<double>{});
    }

    std::vector<Point_2> points;
    std::vector<double> radii;
    points.reserve(skeleton->size_of_vertices());
    radii.reserve(skeleton->size_of_vertices());

    for (auto vertex_iter = skeleton->vertices_begin(); vertex_iter != skeleton->vertices_end(); ++vertex_iter) {
        const auto& point = vertex_iter->point();
        Point_2 point_xy(point.x(), point.y());
        const K::FT sq_dist = squared_distance_to_boundary(point_xy, polygon);

        // Keep only interior MAT samples (exact comparison against zero).
        if (sq_dist > K::FT(0)) {
            points.push_back(point_xy);
            radii.push_back(std::sqrt(std::max(0.0, CGAL::to_double(sq_dist))));
        }
    }

    return std::make_tuple(points, radii);
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

// Walk the straight skeleton, clip edges to the valid cutter-center domain,
// and return trochoid chains for each MAT edge.
std::vector<std::vector<TrochoidArc>>
mat_edge_chains(
    const Polygon_2& boundary,
    const SsPtr& skeleton,
    double tool_radius,
    double radial_clearance,
    double mat_scale,
    double min_trochoid_radius,
    double max_trochoid_radius,
    double pitch,
    int max_passes)
{
    const double min_centerline_distance = tool_radius + radial_clearance;
    const double min_cd_sq = min_centerline_distance * min_centerline_distance;

    auto compute_radius = [&](double boundary_distance) {
        const double available = mat_scale * std::max(0.0, boundary_distance - tool_radius - radial_clearance);
        if (available <= 0.0) return 0.0;
        double r = std::min(max_trochoid_radius, available);
        if (min_trochoid_radius > 0.0) {
            r = std::min(available, std::max(r, min_trochoid_radius));
        }
        return r;
    };

    std::vector<std::vector<TrochoidArc>> chains;
    chains.reserve(static_cast<std::size_t>(max_passes));

    for (auto edge_iter = skeleton->halfedges_begin();
         edge_iter != skeleton->halfedges_end(); ++edge_iter) {
        if (!(edge_iter->is_bisector() || edge_iter->is_inner_bisector())) continue;

        const auto v0 = edge_iter->vertex();
        const auto v1 = edge_iter->opposite()->vertex();
        if (&*v0 >= &*v1) continue;

        Point_2 p0(v0->point().x(), v0->point().y());
        Point_2 p1(v1->point().x(), v1->point().y());

        // Exact squared-distance early rejection (avoids sqrt on discarded edges).
        K::FT sq0 = squared_distance_to_boundary(p0, boundary);
        K::FT sq1 = squared_distance_to_boundary(p1, boundary);
        const K::FT min_cd_sq_ft(min_cd_sq);
        if (sq0 < min_cd_sq_ft && sq1 < min_cd_sq_ft) continue;

        // Edge passes filter — compute distances for clipping and radius.
        double d0 = std::sqrt(std::max(0.0, CGAL::to_double(sq0)));
        double d1 = std::sqrt(std::max(0.0, CGAL::to_double(sq1)));

        // Clip edge to the valid cutter-center domain.
        if (d0 < min_centerline_distance || d1 < min_centerline_distance) {
            const double denom = d1 - d0;
            if (std::abs(denom) < 1e-12) continue;  // division guard
            const double t = std::max(0.0, std::min(1.0, (min_centerline_distance - d0) / denom));
            const Point_2 cp = CGAL::barycenter(p0, 1.0 - t, p1, t);
            // Recompute only the clipped endpoint.
            if (d0 < min_centerline_distance) {
                p0 = cp;
                sq0 = squared_distance_to_boundary(p0, boundary);
                d0 = std::sqrt(std::max(0.0, CGAL::to_double(sq0)));
            } else {
                p1 = cp;
                sq1 = squared_distance_to_boundary(p1, boundary);
                d1 = std::sqrt(std::max(0.0, CGAL::to_double(sq1)));
            }
        }

        // Canonicalize: narrower end first, tie-break via exact predicate.
        if (sq0 > sq1 || (sq0 == sq1 && CGAL::compare_xy(p0, p1) == CGAL::LARGER)) {
            std::swap(p0, p1);
            std::swap(d0, d1);
        }

        const Vector_2 edge_dir = p1 - p0;
        auto circles = trochoid_circles(p0, p1, compute_radius(d0), compute_radius(d1), pitch);
        if (circles.size() < 2) continue;

        auto chain = trochoid_chain(circles, edge_dir, /*climb_milling=*/true);
        if (!chain.empty()) {
            chains.push_back(std::move(chain));
        }

        if (static_cast<int>(chains.size()) >= max_passes) break;
    }

    return chains;
}

} // namespace

std::tuple<compas::RowMatrixXd, compas::RowMatrixXd>
pmp_polygon_medial_axis_transform(
    Eigen::Ref<const compas::RowMatrixXd> vertices)
{
    Polygon_2 polygon = data_to_polygon(vertices);
    auto [mat_points, mat_radii] = polygon_medial_axis_transform_internal(polygon);

    compas::RowMatrixXd points_matrix(mat_points.size(), 3);
    compas::RowMatrixXd radii_matrix(mat_radii.size(), 1);

    for (std::size_t i = 0; i < mat_points.size(); ++i) {
        points_matrix(i, 0) = mat_points[i].x();
        points_matrix(i, 1) = mat_points[i].y();
        points_matrix(i, 2) = 0.0;
        radii_matrix(i, 0) = mat_radii[i];
    }

    return std::make_tuple(points_matrix, radii_matrix);
}

std::vector<compas::RowMatrixXd>
pmp_trochoidal_mat_toolpath(
    Eigen::Ref<const compas::RowMatrixXd> vertices,
    double tool_diameter,
    double stepover,
    double pitch,
    double min_trochoid_radius,
    double max_trochoid_radius,
    double mat_scale,
    double radial_clearance,
    int samples_per_cycle,
    int max_passes)
{
    validate_toolpath_params(tool_diameter, stepover, pitch,
        min_trochoid_radius, max_trochoid_radius, samples_per_cycle, max_passes);

    Polygon_2 boundary = data_to_polygon(vertices);
    SsPtr skeleton = CGAL::create_interior_straight_skeleton_2(
        boundary.vertices_begin(), boundary.vertices_end());
    if (!skeleton) return {};

    const double tool_radius = 0.5 * tool_diameter;
    auto chains = mat_edge_chains(boundary, skeleton, tool_radius,
        radial_clearance, mat_scale, min_trochoid_radius, max_trochoid_radius,
        pitch, max_passes);

    std::vector<compas::RowMatrixXd> toolpaths;
    toolpaths.reserve(chains.size());
    for (auto& chain : chains) {
        auto pts = tessellate_chain(chain, samples_per_cycle);
        deduplicate_consecutive_points(pts, 1e-9);
        if (pts.size() >= 2) {
            toolpaths.push_back(points_to_matrix(pts));
        }
    }
    return toolpaths;
}

static std::tuple<compas::RowMatrixXd, compas::RowMatrixXd, compas::RowMatrixXd, compas::RowMatrixXd, compas::RowMatrixXd>
pmp_trochoidal_mat_toolpath_circular_raw(
    Eigen::Ref<const compas::RowMatrixXd> vertices,
    double tool_diameter,
    double stepover,
    double pitch,
    double min_trochoid_radius,
    double max_trochoid_radius,
    double mat_scale,
    double radial_clearance,
    int samples_per_cycle,
    int max_passes)
{
    validate_toolpath_params(tool_diameter, stepover, pitch,
        min_trochoid_radius, max_trochoid_radius, samples_per_cycle, max_passes);

    Polygon_2 boundary = data_to_polygon(vertices);
    SsPtr skeleton = CGAL::create_interior_straight_skeleton_2(
        boundary.vertices_begin(), boundary.vertices_end());
    if (!skeleton) {
        return std::make_tuple(compas::RowMatrixXd(0, 3), compas::RowMatrixXd(0, 3),
            compas::RowMatrixXd(0, 3), compas::RowMatrixXd(0, 3), compas::RowMatrixXd(0, 1));
    }

    const double tool_radius = 0.5 * tool_diameter;
    auto chains = mat_edge_chains(boundary, skeleton, tool_radius,
        radial_clearance, mat_scale, min_trochoid_radius, max_trochoid_radius,
        pitch, max_passes);

    // Flatten chains into primitives, filter chains without arcs.
    std::vector<TrochoidArc> all_primitives;
    std::vector<int> primitive_path_indices;
    int path_idx = 0;
    for (const auto& chain : chains) {
        bool has_arc = std::any_of(chain.begin(), chain.end(),
            [](const TrochoidArc& a) { return !a.is_line(); });
        if (!has_arc) continue;
        for (const auto& arc : chain) {
            all_primitives.push_back(arc);
            primitive_path_indices.push_back(path_idx);
        }
        path_idx++;
    }

    const int n = static_cast<int>(all_primitives.size());
    compas::RowMatrixXd meta(n, 3);
    compas::RowMatrixXd starts(n, 3);
    compas::RowMatrixXd ends(n, 3);
    compas::RowMatrixXd centers(n, 3);
    compas::RowMatrixXd radii(n, 1);

    for (int i = 0; i < n; ++i) {
        const auto& arc = all_primitives[i];
        meta(i, 0) = static_cast<double>(primitive_path_indices[i]);
        meta(i, 1) = arc.is_line() ? 0.0 : 1.0;
        meta(i, 2) = arc.is_clockwise() ? 1.0 : 0.0;

        starts(i, 0) = CGAL::to_double(arc.start.x());
        starts(i, 1) = CGAL::to_double(arc.start.y());
        starts(i, 2) = 0.0;
        ends(i, 0) = CGAL::to_double(arc.end.x());
        ends(i, 1) = CGAL::to_double(arc.end.y());
        ends(i, 2) = 0.0;
        centers(i, 0) = CGAL::to_double(arc.circle.center().x());
        centers(i, 1) = CGAL::to_double(arc.circle.center().y());
        centers(i, 2) = 0.0;
        radii(i, 0) = arc.radius();
    }

    return std::make_tuple(meta, starts, ends, centers, radii);
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

} // namespace

std::tuple<compas::RowMatrixXd, compas::RowMatrixXd, compas::RowMatrixXd, compas::RowMatrixXd, compas::RowMatrixXd>
pmp_trochoidal_mat_toolpath_circular(
    Eigen::Ref<const compas::RowMatrixXd> vertices,
    double tool_diameter,
    double stepover,
    double pitch,
    double min_trochoid_radius,
    double max_trochoid_radius,
    double mat_scale,
    double radial_clearance,
    int samples_per_cycle,
    int max_passes,
    double lead_in,
    double lead_out,
    bool link_paths,
    bool optimize_order,
    double cut_z,
    double clearance_z,
    bool has_clearance_z,
    bool retract_at_end)
{
    // Validate and build skeleton + chains directly (no serialize/deserialize round-trip)
    validate_toolpath_params(tool_diameter, stepover, pitch,
        min_trochoid_radius, max_trochoid_radius, samples_per_cycle, max_passes);

    Polygon_2 boundary = data_to_polygon(vertices);
    SsPtr skeleton = CGAL::create_interior_straight_skeleton_2(
        boundary.vertices_begin(), boundary.vertices_end());
    if (!skeleton) {
        compas::RowMatrixXd empty_meta(0, 4);
        compas::RowMatrixXd empty3(0, 3);
        compas::RowMatrixXd empty1(0, 1);
        return std::make_tuple(empty_meta, empty3, empty3, empty3, empty1);
    }

    const double tool_radius = 0.5 * tool_diameter;
    auto chains = mat_edge_chains(boundary, skeleton, tool_radius,
        radial_clearance, mat_scale, min_trochoid_radius, max_trochoid_radius,
        pitch, max_passes);

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
        return std::make_tuple(empty_meta, empty3, empty3, empty3, empty1);
    }

    // Greedy nearest-neighbor path ordering
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

        // Compute lead-in start point via start tangent
        const Vector_2 st = path.front().start_tangent();
        const double st_len = approx_length(st);
        Point_2 lead_in_pt = path_start;  // default: no lead-in
        bool has_start_tangent = false;
        if (lead_in > 0.0 && st_len > 1e-12) {  // st_len: division guard
            double inv = 1.0 / st_len;
            lead_in_pt = Point_2(
                CGAL::to_double(path_start.x()) - lead_in * CGAL::to_double(st.x()) * inv,
                CGAL::to_double(path_start.y()) - lead_in * CGAL::to_double(st.y()) * inv);
            has_start_tangent = true;
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
                    operations.push_back(make_tp_line(cur_xy, cut_z, lead_in_pt, cut_z, 3 /*link*/, pidx));
                    cur_xy = lead_in_pt; cur_z = cut_z;
                }
            }
        } else {
            cur_xy = lead_in_pt; cur_z = cut_z;
        }

        // Lead-in
        if (lead_in > 0.0 && has_start_tangent) {
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

        // Lead-out via end tangent
        {
            const Vector_2 et = path.back().end_tangent();
            const double et_len = approx_length(et);
            if (lead_out > 0.0 && et_len > 1e-12) {  // et_len: division guard
                double inv = 1.0 / et_len;
                Point_2 lo_end(
                    CGAL::to_double(path_end.x()) + lead_out * CGAL::to_double(et.x()) * inv,
                    CGAL::to_double(path_end.y()) + lead_out * CGAL::to_double(et.y()) * inv);
                operations.push_back(make_tp_line(path_end, cut_z, lo_end, cut_z, 2 /*lead_out*/, pidx));
                cur_xy = lo_end; cur_z = cut_z;
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
    }

    return std::make_tuple(meta, starts, ends, centers_out, radii_out);
}

NB_MODULE(_toolpath, m)
{
    m.def(
        "polygon_medial_axis_transform",
        &pmp_polygon_medial_axis_transform,
        "Compute MAT samples from polygon interior straight skeleton.\n\n"
        "Parameters\n"
        "----------\n"
        "vertices : array-like\n"
        "    Matrix of polygon vertices (Nx3, float64)\n"
        "\n"
        "Returns\n"
        "-------\n"
        "tuple\n"
        "    - MAT points (Mx3, float64)\n"
        "    - MAT radii (Mx1, float64)",
        "vertices"_a);

    m.def(
        "trochoidal_mat_toolpath",
        &pmp_trochoidal_mat_toolpath,
        "Create MAT-constrained trochoidal toolpaths for a polygon pocket.\n\n"
        "Parameters\n"
        "----------\n"
        "vertices : array-like\n"
        "    Matrix of polygon vertices (Nx3, float64)\n"
        "tool_diameter : float\n"
        "    Cutter diameter\n"
        "stepover : float\n"
        "    Reserved for API compatibility (currently unused)\n"
        "pitch : float\n"
        "    Trochoid pitch\n"
        "min_trochoid_radius : float\n"
        "    Minimum trochoid radius\n"
        "max_trochoid_radius : float\n"
        "    Maximum trochoid radius\n"
        "mat_scale : float\n"
        "    MAT availability scale factor\n"
        "radial_clearance : float\n"
        "    Clearance from available radius\n"
        "samples_per_cycle : int\n"
        "    Samples per trochoid cycle\n"
        "max_passes : int\n"
        "    Maximum number of emitted MAT-edge toolpaths\n"
        "\n"
        "Returns\n"
        "-------\n"
        "list\n"
        "    List of toolpath polylines (each Kx3, float64)",
        "vertices"_a,
        "tool_diameter"_a,
        "stepover"_a,
        "pitch"_a,
        "min_trochoid_radius"_a,
        "max_trochoid_radius"_a,
        "mat_scale"_a,
        "radial_clearance"_a,
        "samples_per_cycle"_a,
        "max_passes"_a);

    m.def(
        "trochoidal_mat_toolpath_circular",
        &pmp_trochoidal_mat_toolpath_circular,
        "Create linked circular toolpath with path ordering, leads, and links.\n\n"
        "Returns tuple:\n"
        "- meta (Nx4 float): [path_index, type(0=line,1=arc,2=circle), clockwise, operation(0=cut,1=lead_in,2=lead_out,3=link,4=retract,5=plunge)]\n"
        "- start points (Nx3 float)\n"
        "- end points (Nx3 float)\n"
        "- centers (Nx3 float)\n"
        "- radii (Nx1 float)",
        "vertices"_a,
        "tool_diameter"_a,
        "stepover"_a,
        "pitch"_a,
        "min_trochoid_radius"_a,
        "max_trochoid_radius"_a,
        "mat_scale"_a,
        "radial_clearance"_a,
        "samples_per_cycle"_a,
        "max_passes"_a,
        "lead_in"_a = 0.0,
        "lead_out"_a = 0.0,
        "link_paths"_a = true,
        "optimize_order"_a = true,
        "cut_z"_a = 0.0,
        "clearance_z"_a = 0.0,
        "has_clearance_z"_a = false,
        "retract_at_end"_a = true);
}
