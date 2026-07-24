#include "stock_2.h"
#include "engagement_2.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>

namespace {

// Convert an Nx3 double matrix (rationals by construction) to a linear
// circle-segment general polygon, validating simplicity via Polygon_2.
GpsPolygon data_to_gps_polygon(Eigen::Ref<const compas::RowMatrixXd> vertices)
{
    if (vertices.rows() < 3) {
        throw std::invalid_argument("Expected at least three polygon vertices.");
    }
    CGAL::Polygon_2<Epeck> simple_check;
    for (int i = 0; i < vertices.rows(); ++i) {
        simple_check.push_back(EPoint(vertices(i, 0), vertices(i, 1)));
    }
    if (!simple_check.is_simple()) {
        throw std::invalid_argument("Polygon boundary must be simple (no self-intersections).");
    }
    if (simple_check.is_clockwise_oriented()) {
        simple_check.reverse_orientation();
    }

    GpsPolygon polygon;
    const std::size_t n = simple_check.size();
    for (std::size_t i = 0; i < n; ++i) {
        const EPoint& a = simple_check[i];
        const EPoint& b = simple_check[(i + 1) % n];
        // Vendored CGAL 6.0.1: X_monotone_curve_2 (_X_monotone_circle_segment_2)
        // has a two-kernel-point constructor but no Segment_2 constructor, so
        // build the linear arc directly from the exact endpoints.
        polygon.push_back(GpsXCurve(a, b));
    }
    return polygon;
}

GpsPolygon exact_linear_polygon(const std::vector<EPoint>& vertices)
{
    GpsPolygon polygon;
    for (std::size_t index = 0; index < vertices.size(); ++index) {
        polygon.push_back(GpsXCurve(
            vertices[index],
            vertices[(index + 1) % vertices.size()]));
    }
    return polygon;
}

Gps exact_disk_union(
    const std::vector<EPoint>& centers,
    const Epeck::FT& radius)
{
    std::vector<GpsPolygon> disks;
    disks.reserve(centers.size());
    for (const EPoint& center : centers) {
        disks.push_back(disk_polygon(center, radius));
    }
    Gps region;
    region.join(disks.begin(), disks.end());
    return region;
}

void validate_depletion_trace(const DepletionTrace& trace)
{
    if (trace.center_count != trace.center_parameters.size()
        || !trace.exact_incidence
        || !trace.exact_parameters_in_range
        || !trace.exact_anchors_present
        || !trace.exact_removal_radius_valid
        || !trace.exact_chord_bound_holds
        || (trace.cyclic && !trace.exact_seam_chord_bound_holds)
        || trace.strategy_version.empty()) {
        throw ExactDepletionConstructionError(
            "exact depletion structural trace failed validation.");
    }
}

bool exact_set_is_subset(const Gps& subset, const Gps& superset)
{
    Gps difference(subset);
    difference.difference(superset);
    return difference.is_empty();
}

} // namespace

// Full disk as a two-arc general polygon (counterclockwise). Split at the two
// x-extreme (vertical-tangency) points; radius comes from a double so it is an
// exact rational, hence center.x() +/- radius are plain rationals that GpsPoint
// accepts directly (Task 1 pattern). This mirrors CGAL's own full-circle
// subdivision in Arr_circle_segment_traits_2::Make_x_monotone_2.
// Declared in stock_2.h with external linkage so engagement_2.cpp can build the
// same exact cutter disk it intersects the stock against.
GpsPolygon disk_polygon(const EPoint& center, const Epeck::FT& radius)
{
    const Epeck::FT r_sq = radius * radius;
    const ECircle circle(center, r_sq, CGAL::COUNTERCLOCKWISE);
    const GpsPoint p_min(center.x() - radius, center.y());  // leftmost point
    const GpsPoint p_max(center.x() + radius, center.y());  // rightmost point
    GpsPolygon polygon;
    // Vendored CGAL 6.0.1: X_monotone_curve_2(Circle_2, source, target,
    // Orientation) — verified against Arr_geometry_traits/Circle_segment_2.h.
    polygon.push_back(GpsXCurve(circle, p_min, p_max, CGAL::COUNTERCLOCKWISE));  // lower half
    polygon.push_back(GpsXCurve(circle, p_max, p_min, CGAL::COUNTERCLOCKWISE));  // upper half
    return polygon;
}

Stock2::Stock2(Eigen::Ref<const compas::RowMatrixXd> boundary,
               const std::vector<compas::RowMatrixXd>& holes)
    : set_(std::make_unique<Gps>())
{
    set_->insert(data_to_gps_polygon(boundary));
    for (const auto& hole : holes) {
        Gps hole_set;
        hole_set.insert(data_to_gps_polygon(hole));
        set_->difference(hole_set);
    }
}

Stock2::Stock2(std::unique_ptr<Gps> set) noexcept
    : set_(std::move(set))
{
}

bool Stock2::contains(double x, double y) const
{
    return set_->oriented_side(GpsPoint(Epeck::FT(x), Epeck::FT(y))) == CGAL::ON_POSITIVE_SIDE;
}

bool Stock2::is_empty() const
{
    return set_->is_empty();
}

Stock2 Stock2::clone() const
{
    return Stock2(std::make_unique<Gps>(*set_));
}

bool Stock2::is_subset_of(const Stock2& other) const
{
    return exact_set_is_subset(*set_, *other.set_);
}

bool Stock2::exactly_equals(const Stock2& other) const
{
    return exact_set_is_subset(*set_, *other.set_)
        && exact_set_is_subset(*other.set_, *set_);
}

// Fraction of the tool radius allowed as chain under-coverage slack; the
// generator's radial_clearance margin (1e-3 * D = 2e-3 * r) dominates it.
constexpr double CHAIN_SLACK_FRACTION = 1e-4;

void Stock2::subtract_disk(double cx, double cy, double radius)
{
    if (radius <= 0.0) throw std::invalid_argument("radius should be positive.");
    Gps region;
    region.insert(disk_polygon(EPoint(cx, cy), Epeck::FT(radius)));
    set_->difference(region);
}

// Subtract the union of exact tool disks centered at the given points. One
// chain implementation shared by the capsule (Task 2) and arc (Task 3) paths:
// callers only choose where the centers sit; exact predicates still decide
// point-in on the constructed disks, so the union under-covers the true sweep.
void Stock2::subtract_point_chain(const std::vector<std::pair<double, double>>& centers,
                                  double radius)
{
    // AGGREGATED union: build every disk, union them in ONE sweep-based pass, then a
    // single difference. The prior per-disk `region.join(disk_set)` was the O(N^2)
    // Boolean_set_operations anti-pattern -- each incremental join re-swept the whole
    // growing union. General_polygon_set_2::join(first, last) unions the range in one
    // arrangement. The result is identical (union is associative/commutative) -- this
    // is purely a cost fix, not a semantic change.
    std::vector<GpsPolygon> disks;
    disks.reserve(centers.size());
    for (const auto& [x, y] : centers) {
        disks.push_back(disk_polygon(EPoint(x, y), Epeck::FT(radius)));
    }
    Gps region;
    region.join(disks.begin(), disks.end());
    set_->difference(region);
}

void Stock2::subtract_capsule(double x0, double y0, double x1, double y1, double radius)
{
    if (radius <= 0.0) throw std::invalid_argument("radius should be positive.");
    const double dx = x1 - x0;
    const double dy = y1 - y0;
    const double len_sq = dx * dx + dy * dy;
    if (len_sq == 0.0) {
        subtract_disk(x0, y0, radius);
        return;
    }
    // Disk chain spacing s with per-disk radius r: the chain covers the true
    // capsule shrunk by delta = r - sqrt(r^2 - (s/2)^2) <= s^2/(4r). Choose s
    // so delta <= CHAIN_SLACK_FRACTION * r  =>  s = 2r*sqrt(fraction). s only
    // sets the (certified, under-covering) construction density.
    const double len = std::sqrt(len_sq);
    const double spacing = 2.0 * radius * std::sqrt(CHAIN_SLACK_FRACTION);
    const int n = std::max(1, static_cast<int>(std::ceil(len / spacing)));
    std::vector<std::pair<double, double>> centers;
    centers.reserve(n + 1);
    for (int i = 0; i <= n; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(n);
        centers.emplace_back(x0 + t * dx, y0 + t * dy);
    }
    subtract_point_chain(centers, radius);
}

void Stock2::subtract_arc_sweep(double cx, double cy, double sx, double sy,
                                double ex, double ey, bool cw, double tool_radius)
{
    if (tool_radius <= 0.0) throw std::invalid_argument("tool_radius should be positive.");
    const double rx = sx - cx, ry = sy - cy;
    const double guide_r = std::hypot(rx, ry);
    if (guide_r == 0.0) { subtract_disk(cx, cy, tool_radius); return; }

    double a0 = std::atan2(ry, rx);
    double a1 = std::atan2(ey - cy, ex - cx);
    double sweep = cw ? a0 - a1 : a1 - a0;
    if (sweep <= 0.0) sweep += 2.0 * std::numbers::pi;
    const bool full = (sx == ex && sy == ey);
    if (full) sweep = 2.0 * std::numbers::pi;

    // Chain spacing along the guide: disks of tool_radius at arc-length step s
    // under-cover the true sweep by delta <= s^2/(4*tool_radius) (chord
    // sagitta bound; the straight-chord bound is conservative for the arc chain
    // because consecutive centers are closer along the chord than the arc).
    const double spacing = 2.0 * tool_radius * std::sqrt(CHAIN_SLACK_FRACTION);
    const double arc_len = guide_r * sweep;
    const int n = std::max(4, static_cast<int>(std::ceil(arc_len / spacing)));

    std::vector<std::pair<double, double>> centers;
    centers.reserve(n + 1);
    for (int i = 0; i <= n; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(n);
        const double a = cw ? a0 - sweep * t : a0 + sweep * t;
        centers.emplace_back(cx + guide_r * std::cos(a), cy + guide_r * std::sin(a));
    }
    subtract_point_chain(centers, tool_radius);
}

DepletionTrace Stock2::subtract_exact_segment(
    const ExactSegmentMotion2& motion,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit)
{
    ExactDepletionConstruction2 construction = construct_exact_segment_depletion(
        motion,
        tool_radius,
        max_chord,
        center_count_limit);
    std::unique_ptr<Gps> trial = std::make_unique<Gps>(*set_);
    Gps removal = exact_disk_union(construction.centers, tool_radius);
    trial->difference(removal);
    validate_depletion_trace(construction.trace);
    set_.swap(trial);
    return std::move(construction.trace);
}

DepletionTrace Stock2::subtract_exact_full_circle(
    const ExactCircleMotion2& motion,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit)
{
    ExactDepletionConstruction2 construction = construct_exact_full_circle_depletion(
        motion,
        tool_radius,
        max_chord,
        center_count_limit);
    std::unique_ptr<Gps> trial = std::make_unique<Gps>(*set_);
    Gps removal = exact_disk_union(construction.centers, tool_radius);
    trial->difference(removal);
    validate_depletion_trace(construction.trace);
    set_.swap(trial);
    return std::move(construction.trace);
}

bool exact_segment_undercover_holds(
    const ExactSegmentMotion2& motion,
    const Epeck::FT& exact_length,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit)
{
    if (CGAL::compare(exact_length, Epeck::FT(0)) != CGAL::LARGER
        || CGAL::compare(
               (motion.end - motion.start).squared_length(),
               exact_length * exact_length)
            != CGAL::EQUAL) {
        throw ExactDepletionConstructionError(
            "segment sweep oracle requires an exact positive Pythagorean length.");
    }
    ExactDepletionConstruction2 construction = construct_exact_segment_depletion(
        motion,
        tool_radius,
        max_chord,
        center_count_limit);
    Gps removal = exact_disk_union(construction.centers, tool_radius);

    const EVector direction = motion.end - motion.start;
    const EVector normal(
        -direction.y() * tool_radius / exact_length,
        direction.x() * tool_radius / exact_length);
    const std::vector<EPoint> rectangle{
        motion.start - normal,
        motion.end - normal,
        motion.end + normal,
        motion.start + normal,
    };
    std::vector<GpsPolygon> sweep_parts{
        exact_linear_polygon(rectangle),
        disk_polygon(motion.start, tool_radius),
        disk_polygon(motion.end, tool_radius),
    };
    Gps sweep;
    sweep.join(sweep_parts.begin(), sweep_parts.end());
    return exact_set_is_subset(removal, sweep);
}

bool exact_full_circle_undercover_holds(
    const ExactCircleMotion2& motion,
    const Epeck::FT& guide_radius,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit)
{
    if (CGAL::compare(guide_radius, Epeck::FT(0)) != CGAL::LARGER
        || CGAL::compare(
               motion.phase_vector.squared_length(),
               guide_radius * guide_radius)
            != CGAL::EQUAL) {
        throw ExactDepletionConstructionError(
            "circle sweep oracle requires an exact positive rational guide radius.");
    }
    ExactDepletionConstruction2 construction = construct_exact_full_circle_depletion(
        motion,
        tool_radius,
        max_chord,
        center_count_limit);
    Gps removal = exact_disk_union(construction.centers, tool_radius);

    Gps sweep;
    sweep.insert(disk_polygon(
        motion.center,
        guide_radius + tool_radius));
    if (CGAL::compare(guide_radius, tool_radius) == CGAL::LARGER) {
        Gps inner;
        inner.insert(disk_polygon(
            motion.center,
            guide_radius - tool_radius));
        sweep.difference(inner);
    }
    return exact_set_is_subset(removal, sweep);
}

NB_MODULE(_stock_2, m)
{
    nb::exception<ExactDepletionConstructionError> construction_error(
        m,
        "ExactDepletionConstructionError");
    nb::exception<ExactDepletionCenterLimitError>(
        m,
        "ExactDepletionCenterLimitError",
        construction_error.ptr());

    nb::class_<DepletionTrace>(m, "DepletionTrace")
        .def_prop_ro("center_count", [](const DepletionTrace& trace) {
            return trace.center_count;
        })
        .def_prop_ro("center_parameters", [](const DepletionTrace& trace) {
            std::vector<std::tuple<int, std::size_t, std::size_t>> result;
            result.reserve(trace.center_parameters.size());
            for (const ExactCenterParameter2& parameter : trace.center_parameters) {
                result.emplace_back(
                    parameter.chart,
                    parameter.numerator,
                    parameter.denominator);
            }
            return result;
        })
        .def_prop_ro("max_chord", [](const DepletionTrace& trace) {
            return CGAL::to_double(trace.max_chord);
        })
        .def_prop_ro("removal_radius", [](const DepletionTrace& trace) {
            return CGAL::to_double(trace.removal_radius);
        })
        .def_prop_ro("strategy_version", [](const DepletionTrace& trace) {
            return nb::bytes(
                trace.strategy_version.data(),
                trace.strategy_version.size());
        })
        .def_ro("cyclic", &DepletionTrace::cyclic)
        .def_ro("exact_incidence", &DepletionTrace::exact_incidence)
        .def_ro("exact_parameters_in_range", &DepletionTrace::exact_parameters_in_range)
        .def_ro("exact_anchors_present", &DepletionTrace::exact_anchors_present)
        .def_ro("exact_removal_radius_valid", &DepletionTrace::exact_removal_radius_valid)
        .def_ro("exact_chord_bound_holds", &DepletionTrace::exact_chord_bound_holds)
        .def_ro("exact_seam_chord_bound_holds", &DepletionTrace::exact_seam_chord_bound_holds);

    nb::class_<Stock2>(m, "Stock2")
        .def(nb::init<Eigen::Ref<const compas::RowMatrixXd>,
                      const std::vector<compas::RowMatrixXd>&>(),
             "boundary"_a, "holes"_a)
        .def("contains", &Stock2::contains, "x"_a, "y"_a)
        .def("is_empty", &Stock2::is_empty)
        .def("clone", &Stock2::clone)
        .def("is_subset_of", &Stock2::is_subset_of, "other"_a)
        .def("exactly_equals", &Stock2::exactly_equals, "other"_a)
        .def("subtract_capsule", &Stock2::subtract_capsule,
             "x0"_a, "y0"_a, "x1"_a, "y1"_a, "radius"_a)
        .def("subtract_arc_sweep", &Stock2::subtract_arc_sweep,
             "cx"_a, "cy"_a, "sx"_a, "sy"_a, "ex"_a, "ey"_a, "cw"_a, "tool_radius"_a)
        .def(
            "subtract_exact_segment",
            [](Stock2& stock,
               double x0,
               double y0,
               double x1,
               double y1,
               double tool_radius,
               double max_chord,
               std::size_t center_count_limit) {
                return stock.subtract_exact_segment(
                    {EPoint(x0, y0), EPoint(x1, y1)},
                    Epeck::FT(tool_radius),
                    Epeck::FT(max_chord),
                    center_count_limit);
            },
            "x0"_a,
            "y0"_a,
            "x1"_a,
            "y1"_a,
            "tool_radius"_a,
            "max_chord"_a,
            "center_count_limit"_a)
        .def(
            "subtract_exact_full_circle",
            [](Stock2& stock,
               double cx,
               double cy,
               double phase_x,
               double phase_y,
               bool clockwise,
               double tool_radius,
               double max_chord,
               std::size_t center_count_limit) {
                return stock.subtract_exact_full_circle(
                    {EPoint(cx, cy), EVector(phase_x, phase_y), clockwise},
                    Epeck::FT(tool_radius),
                    Epeck::FT(max_chord),
                    center_count_limit);
            },
            "cx"_a,
            "cy"_a,
            "phase_x"_a,
            "phase_y"_a,
            "clockwise"_a,
            "tool_radius"_a,
            "max_chord"_a,
            "center_count_limit"_a)
        .def("subtract_disk", &Stock2::subtract_disk, "cx"_a, "cy"_a, "radius"_a);

    m.def(
        "exact_segment_point_is_incident",
        [](double x0, double y0, double x1, double y1, double px, double py) {
            return exact_segment_point_is_incident(
                {EPoint(x0, y0), EPoint(x1, y1)},
                EPoint(px, py));
        },
        "x0"_a,
        "y0"_a,
        "x1"_a,
        "y1"_a,
        "px"_a,
        "py"_a);
    m.def(
        "exact_circle_point_is_incident",
        [](double cx,
           double cy,
           double phase_x,
           double phase_y,
           double px,
           double py) {
            return exact_circle_point_is_incident(
                {EPoint(cx, cy), EVector(phase_x, phase_y), false},
                EPoint(px, py));
        },
        "cx"_a,
        "cy"_a,
        "phase_x"_a,
        "phase_y"_a,
        "px"_a,
        "py"_a);
    m.def(
        "exact_segment_undercover_holds",
        [](double x0,
           double y0,
           double x1,
           double y1,
           double exact_length,
           double tool_radius,
           double max_chord,
           std::size_t center_count_limit) {
            return exact_segment_undercover_holds(
                {EPoint(x0, y0), EPoint(x1, y1)},
                Epeck::FT(exact_length),
                Epeck::FT(tool_radius),
                Epeck::FT(max_chord),
                center_count_limit);
        },
        "x0"_a,
        "y0"_a,
        "x1"_a,
        "y1"_a,
        "exact_length"_a,
        "tool_radius"_a,
        "max_chord"_a,
        "center_count_limit"_a);
    m.def(
        "exact_full_circle_undercover_holds",
        [](double cx,
           double cy,
           double phase_x,
           double phase_y,
           double guide_radius,
           double tool_radius,
           double max_chord,
           std::size_t center_count_limit) {
            return exact_full_circle_undercover_holds(
                {EPoint(cx, cy), EVector(phase_x, phase_y), false},
                Epeck::FT(guide_radius),
                Epeck::FT(tool_radius),
                Epeck::FT(max_chord),
                center_count_limit);
        },
        "cx"_a,
        "cy"_a,
        "phase_x"_a,
        "phase_y"_a,
        "guide_radius"_a,
        "tool_radius"_a,
        "max_chord"_a,
        "center_count_limit"_a);

    register_engagement(m);
}
