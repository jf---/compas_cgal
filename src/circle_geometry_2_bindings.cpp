// Exact circle decisions for the approximate predecessor-placement consumer.
#include "boundary_contact_projection_2.h"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/squared_distance_2.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/vector.h>

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

FT intersection_numerator(const FT& first_squared, const FT& second_squared,
                          const FT& distance_squared)
{
    if (CGAL::is_zero(distance_squared)) {
        throw NoCircleIntersectionError("Concentric disks have no unique intersection chord.");
    }
    const FT numerator = first_squared - second_squared - distance_squared;
    if (CGAL::sign(4 * distance_squared * second_squared - CGAL::square(numerator)) == CGAL::NEGATIVE) {
        throw NoCircleIntersectionError("Swept disks have no real intersection.");
    }
    return numerator;
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
    const FT numerator = intersection_numerator(first_squared, second_squared, distance_squared);
    const FT height_squared = second_squared
        - CGAL::square(numerator) / (4 * distance_squared);
    // Reporting coordinates in the successor-centered frame, +x from first
    // center to second. Existence was decided above on exact quantities.
    return {CGAL::to_double(numerator)
                / (2 * std::sqrt(CGAL::to_double(distance_squared))),
            CGAL::to_double(height_squared)};
}
double corrected_engagement_cosine_squared(
    const XY& first, double first_radius, const XY& second,
    double second_radius, double tool_radius)
{
    const Point a = point(first), b = point(second);
    const FT tool = radius(tool_radius), current_radius = radius(second_radius);
    const FT outer = current_radius + tool;
    const FT first_squared = CGAL::square(radius(first_radius) + tool);
    const FT distance_squared = CGAL::squared_distance(a, b);
    const FT numerator = intersection_numerator(first_squared, CGAL::square(outer), distance_squared);
    // q is the radial projection of the swept-disk intersection onto the
    // machining circle. Expanding |q-a|² cancels both square roots:
    // rho² + (rho / outer_radius) * numerator + center_distance².
    const FT q_distance_squared = CGAL::square(current_radius)
        + (current_radius / outer) * numerator + distance_squared;
    if (CGAL::sign(q_distance_squared) != CGAL::POSITIVE) {
        throw NoCircleIntersectionError("Corrected contact has no unique chord direction.");
    }
    const FT chord_numerator = q_distance_squared + CGAL::square(tool) - first_squared;
    const FT ratio_squared = CGAL::square(chord_numerator)
        / (4 * CGAL::square(tool) * q_distance_squared);
    if (CGAL::compare(ratio_squared, FT(1)) == CGAL::LARGER) {
        throw NoCircleIntersectionError("Corrected tool circle has no real contour intersection.");
    }
    return CGAL::to_double(ratio_squared);
}
} // namespace

NB_MODULE(_circle_geometry_2, m)
{
    nb::exception<boundary_contact_projection::InvalidBoundaryProjectionError>(m, "InvalidBoundaryProjectionError");
    nb::exception<boundary_contact_projection::AmbiguousBoundaryProjectionError>(m, "AmbiguousBoundaryProjectionError");
    nb::exception<boundary_contact_projection::BoundaryProjectionDistanceError>(m, "BoundaryProjectionDistanceError");
    m.def("project_boundary_contact", &boundary_contact_projection::project);
    nb::exception<InvalidCircleGeometryError>(m, "InvalidCircleGeometryError");
    nb::exception<NoCircleIntersectionError>(m, "NoCircleIntersectionError");
    m.def("disk_contains_disk", &disk_contains_disk);
    m.def("swept_disk_intersection", &swept_disk_intersection);
    m.def("corrected_engagement_cosine_squared", &corrected_engagement_cosine_squared);
    m.def("orientation", [](const XY& a, const XY& b, const XY& c) {
        return static_cast<int>(CGAL::orientation(point(a), point(b), point(c)));
    });
}
