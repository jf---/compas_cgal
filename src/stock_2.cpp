#include "stock_2.h"
#include "engagement_2.h"

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

NB_MODULE(_stock_2, m)
{
    nb::class_<Stock2>(m, "Stock2")
        .def(nb::init<Eigen::Ref<const compas::RowMatrixXd>,
                      const std::vector<compas::RowMatrixXd>&>(),
             "boundary"_a, "holes"_a)
        .def("contains", &Stock2::contains, "x"_a, "y"_a)
        .def("is_empty", &Stock2::is_empty);

    register_engagement(m);
}
