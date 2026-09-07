#include "reachable_boundary_sampling_2.h"
#include "reachable_errors_2.h"

#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/number_utils.h>
#include <cmath>

namespace {
ReachFT lift(const ReachPoint::CoordNT& coordinate)
{
    if (!coordinate.is_extended()) return coordinate.a0();
    return coordinate.a0() + coordinate.a1() * CGAL::sqrt(coordinate.root());
}
}

ReachKernelPoint reachable_kernel_point(const ReachPoint& point)
{
    return ReachKernelPoint(lift(point.x()), lift(point.y()));
}

ReachPoint sample_reachable_boundary(const ReachableBoundaryCurve2& primitive, double parameter)
{
    if (!std::isfinite(parameter) || parameter < 0 || parameter > 1) {
        throw ReachableArrangementTopologyError("Boundary sampling parameter must lie in [0,1].");
    }
    const auto& curve = primitive.curve;
    if (parameter == 0) return curve.source();
    if (parameter == 1) return curve.target();
    const ReachFT u(parameter);
    const auto start = reachable_kernel_point(curve.source());
    const auto end = reachable_kernel_point(curve.target());
    if (curve.is_linear()) {
        const auto point = CGAL::barycenter(start, ReachFT(1) - u, end, u);
        return ReachPoint(point.x(), point.y());
    }
    const auto circle = curve.supporting_circle();
    const auto center = circle.center();
    const auto first = start - center, last = end - center;
    auto middle = first + last;
    if (middle == ReachKernel::Vector_2(CGAL::NULL_VECTOR)) {
        // Antipodal endpoints: the directed perpendicular selects the intended
        // semicircle. Every x-monotone circular primitive spans at most pi.
        middle = first.perpendicular(curve.orientation());
    }
    middle = middle * CGAL::sqrt(circle.squared_radius() / middle.squared_length());
    if (CGAL::compare(u, ReachFT(1) / ReachFT(2)) == CGAL::EQUAL) {
        const auto point = center + middle;
        return ReachPoint(point.x(), point.y());
    }
    const bool first_half = CGAL::compare(u, ReachFT(1) / ReachFT(2)) != CGAL::LARGER;
    const ReachFT local = first_half ? ReachFT(2) * u : ReachFT(2) * u - ReachFT(1);
    const auto a = first_half ? first : middle;
    const auto b = first_half ? middle : last;
    const auto chord = (ReachFT(1) - local) * a + local * b;
    // The sample lies on the supporting circle by construction. Re-deciding
    // that incidence would be an identically-zero CORE decision, refined to
    // its root bound on every call; the sampling tests witness it instead.
    const auto point = center + chord * CGAL::sqrt(circle.squared_radius() / chord.squared_length());
    return ReachPoint(point.x(), point.y());
}
