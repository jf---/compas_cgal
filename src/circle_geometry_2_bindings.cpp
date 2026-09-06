// Exact circle decisions for the approximate predecessor-placement consumer.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/squared_distance_2.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/pair.h>

#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace nb = nanobind;
namespace {
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using FT = Kernel::FT;
using Point = Kernel::Point_2;
using XY = std::array<double, 2>;

class InvalidCircleGeometryError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};
class NoCircleIntersectionError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

Point point(const XY& xy)
{
    if (!std::isfinite(xy[0]) || !std::isfinite(xy[1])) {
        throw InvalidCircleGeometryError("Circle coordinates must be finite.");
    }
    return Point(xy[0], xy[1]);
}

FT radius(double value)
{
    if (!std::isfinite(value) || value <= 0) {
        throw InvalidCircleGeometryError("Circle radii must be finite and positive.");
    }
    return FT(value);
}

bool disk_contains_disk(const XY& first, double first_radius,
                        const XY& second, double second_radius)
{
    const Point a = point(first), b = point(second);
    const FT difference = radius(first_radius) - radius(second_radius);
    return CGAL::sign(difference) != CGAL::NEGATIVE
        && CGAL::compare_squared_distance(a, b, CGAL::square(difference))
               != CGAL::LARGER;
}

std::pair<double, double> swept_disk_intersection(
    const XY& first, double first_radius, const XY& second,
    double second_radius, double tool_radius)
{
    const Point a = point(first), b = point(second);
    const FT tool = radius(tool_radius);
    const FT first_squared = CGAL::square(radius(first_radius) + tool);
    const FT second_squared = CGAL::square(radius(second_radius) + tool);
    const FT distance_squared = CGAL::squared_distance(a, b);
    if (CGAL::is_zero(distance_squared)) {
        throw NoCircleIntersectionError("Concentric disks have no unique intersection chord.");
    }
    const FT numerator = first_squared - second_squared - distance_squared;
    const FT height_squared = second_squared
        - CGAL::square(numerator) / (4 * distance_squared);
    if (CGAL::sign(height_squared) == CGAL::NEGATIVE) {
        throw NoCircleIntersectionError("Swept disks have no real intersection.");
    }
    // Reporting coordinates in the successor-centered frame, +x from first
    // center to second. Existence was decided above on exact quantities.
    return {CGAL::to_double(numerator)
                / (2 * std::sqrt(CGAL::to_double(distance_squared))),
            CGAL::to_double(height_squared)};
}
} // namespace

NB_MODULE(_circle_geometry_2, m)
{
    nb::exception<InvalidCircleGeometryError>(m, "InvalidCircleGeometryError");
    nb::exception<NoCircleIntersectionError>(m, "NoCircleIntersectionError");
    m.def("disk_contains_disk", &disk_contains_disk);
    m.def("swept_disk_intersection", &swept_disk_intersection);
    m.def("orientation", [](const XY& a, const XY& b, const XY& c) {
        return static_cast<int>(CGAL::orientation(point(a), point(b), point(c)));
    });
}
