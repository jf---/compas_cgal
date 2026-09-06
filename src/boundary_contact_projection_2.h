#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/squared_distance_2.h>
#include <array>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace boundary_contact_projection {
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using FT = Kernel::FT;
using Point = Kernel::Point_2;
using XY = std::array<double, 2>;

class InvalidBoundaryProjectionError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};
class AmbiguousBoundaryProjectionError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};
class BoundaryProjectionDistanceError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

inline Point point(const XY& xy)
{
    if (!std::isfinite(xy[0]) || !std::isfinite(xy[1])) {
        throw InvalidBoundaryProjectionError("Boundary coordinates must be finite.");
    }
    return Point(xy[0], xy[1]);
}

// Select on exact geometry. The returned parameter is an approximate view for
// the polygon-draft consumer, not a native exact-contact identity.
inline std::pair<std::size_t, double> project(
    const std::vector<XY>& boundary, const XY& contact, double maximum_distance)
{
    if (boundary.size() < 3 || !std::isfinite(maximum_distance) || maximum_distance < 0) {
        throw InvalidBoundaryProjectionError("Projection needs a polygon and nonnegative finite distance bound.");
    }
    const Point query = point(contact);
    std::optional<Point> closest;
    FT minimum, parameter;
    std::size_t side = 0;
    bool ambiguous = false;
    for (std::size_t i = 0; i < boundary.size(); ++i) {
        const Point a = point(boundary[i]);
        const Point b = point(boundary[(i + 1) % boundary.size()]);
        if (a == b) {
            throw InvalidBoundaryProjectionError("Boundary contains a zero-length side.");
        }
        const auto edge = b - a;
        FT t = ((query - a) * edge) / edge.squared_length();
        if (CGAL::sign(t) == CGAL::NEGATIVE) t = FT(0);
        if (CGAL::compare(t, FT(1)) == CGAL::LARGER) t = FT(1);
        const Point projected = a + t * edge;
        const FT distance = CGAL::squared_distance(query, projected);
        if (!closest || CGAL::compare(distance, minimum) == CGAL::SMALLER) {
            closest = projected;
            minimum = distance;
            side = i;
            parameter = t;
            ambiguous = false;
        } else if (CGAL::compare(distance, minimum) == CGAL::EQUAL && projected != *closest) {
            ambiguous = true;
        }
    }
    if (CGAL::compare(minimum, CGAL::square(FT(maximum_distance))) == CGAL::LARGER) {
        throw BoundaryProjectionDistanceError("Closest offset contact exceeds the boundary evidence bound.");
    }
    if (ambiguous) {
        throw AmbiguousBoundaryProjectionError("Distinct offset contacts have equal minimum distance.");
    }
    if (CGAL::compare(parameter, FT(1)) == CGAL::EQUAL) {
        side = (side + 1) % boundary.size();
        parameter = FT(0);
    }
    return {side, CGAL::to_double(parameter)};
}
} // namespace boundary_contact_projection
