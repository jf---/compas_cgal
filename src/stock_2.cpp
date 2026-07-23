#include "stock_2.h"
#include "engagement_2.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

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

// Full disk as a two-arc general polygon (counterclockwise). Split at the two
// x-extreme (vertical-tangency) points; radius comes from a double so it is an
// exact rational, hence center.x() +/- radius are plain rationals that GpsPoint
// accepts directly (Task 1 pattern). This mirrors CGAL's own full-circle
// subdivision in Arr_circle_segment_traits_2::Make_x_monotone_2.
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

} // namespace

Stock2::Stock2(Eigen::Ref<const compas::RowMatrixXd> boundary,
               const std::vector<compas::RowMatrixXd>& holes)
{
    set_.insert(data_to_gps_polygon(boundary));
    for (const auto& hole : holes) {
        Gps hole_set;
        hole_set.insert(data_to_gps_polygon(hole));
        set_.difference(hole_set);
    }
}

bool Stock2::contains(double x, double y) const
{
    return set_.oriented_side(GpsPoint(Epeck::FT(x), Epeck::FT(y))) == CGAL::ON_POSITIVE_SIDE;
}

bool Stock2::is_empty() const
{
    return set_.is_empty();
}

// Fraction of the tool radius allowed as chain under-coverage slack; the
// generator's radial_clearance margin (1e-3 * D = 2e-3 * r) dominates it.
constexpr double CHAIN_SLACK_FRACTION = 1e-4;

void Stock2::subtract_disk(double cx, double cy, double radius)
{
    if (radius <= 0.0) throw std::invalid_argument("radius should be positive.");
    Gps region;
    region.insert(disk_polygon(EPoint(cx, cy), Epeck::FT(radius)));
    set_.difference(region);
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
    // so delta <= CHAIN_SLACK_FRACTION * r  =>  s = 2r*sqrt(fraction). Exact
    // predicates still decide point-in on the constructed disks; s only sets
    // the (certified, under-covering) construction density.
    const double len = std::sqrt(len_sq);
    const double spacing = 2.0 * radius * std::sqrt(CHAIN_SLACK_FRACTION);
    const int n = std::max(1, static_cast<int>(std::ceil(len / spacing)));
    Gps region;
    for (int i = 0; i <= n; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(n);
        Gps disk_set;
        disk_set.insert(disk_polygon(EPoint(x0 + t * dx, y0 + t * dy), Epeck::FT(radius)));
        region.join(disk_set);
    }
    set_.difference(region);
}

NB_MODULE(_stock_2, m)
{
    nb::class_<Stock2>(m, "Stock2")
        .def(nb::init<Eigen::Ref<const compas::RowMatrixXd>,
                      const std::vector<compas::RowMatrixXd>&>(),
             "boundary"_a, "holes"_a)
        .def("contains", &Stock2::contains, "x"_a, "y"_a)
        .def("is_empty", &Stock2::is_empty)
        .def("subtract_capsule", &Stock2::subtract_capsule,
             "x0"_a, "y0"_a, "x1"_a, "y1"_a, "radius"_a)
        .def("subtract_disk", &Stock2::subtract_disk, "cx"_a, "cy"_a, "radius"_a);

    register_engagement(m);
}
