#include "reachable_boundary_sampling_2.h"
#include "reachable_errors_2.h"

#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/number_utils.h>
#include <cmath>
#include <numbers>

namespace {
// Pieces whose angular extent is below this take the exact bisecting sample.
// A chord direction rounded to double places the sample within about 1e-16
// rad of the target angle; 1e-9 rad keeps seven orders of margin for every
// ordinary piece while catching the x-monotone slivers that splitting leaves
// within double rounding of a vertical tangency.
constexpr double SLIVER_SWEEP_RAD = 1e-9;

ReachFT lift(const ReachPoint::CoordNT& coordinate, const ReachFT& root_value)
{
    if (!coordinate.is_extended()) return coordinate.a0();
    return coordinate.a0() + coordinate.a1() * root_value;
}
}

ReachKernelPoint reachable_kernel_point(const ReachPoint& point)
{
    const auto& x = point.x();
    const auto& y = point.y();
    // One radical node per point when both coordinates share the root, so
    // that CORE bounds one square root instead of two independent ones.
    const bool shared = x.is_extended() && y.is_extended() && x.root() == y.root();
    const ReachFT root_x = x.is_extended() ? CGAL::sqrt(x.root()) : ReachFT(0);
    const ReachFT root_y = shared ? root_x : (y.is_extended() ? CGAL::sqrt(y.root()) : ReachFT(0));
    return ReachKernelPoint(lift(x, root_x), lift(y, root_y));
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
    // Rational chord sampling. Normalising a bisecting vector onto the circle
    // needs a division by its squared length, a squared sum of coordinate
    // differences whose value is tiny against the coordinates it is built
    // from; CORE's floating filter charges that quotient at the divisor's
    // condition and every decision downstream inherited the bill (measured
    // 2026-09-08: Monstera queries at 0.3-0.9 s from the midpoint against
    // 3-10 ms from an endpoint). Instead the target angle is chosen on
    // reporting values, and the sample is the second intersection of the
    // circle with the chord from the piece's start in a double direction.
    // That point lies exactly on the exact circle, involves no cancelling
    // divisor, and its membership in the piece is verified exactly below.
    const auto circle = curve.supporting_circle();
    const auto center = circle.center();
    const auto first_radial = start - center, last_radial = end - center;
    const double cx = CGAL::to_double(center.x()), cy = CGAL::to_double(center.y());
    const double start_angle = std::atan2(CGAL::to_double(start.y()) - cy, CGAL::to_double(start.x()) - cx);
    const double end_angle = std::atan2(CGAL::to_double(end.y()) - cy, CGAL::to_double(end.x()) - cx);
    double sweep = end_angle - start_angle;
    // A sliver's raw difference may be zero or carry the wrong sign, so it
    // is recognised before the orientation wrap and again after it: no
    // x-monotone piece spans more than a half turn, so a near-full turn is
    // a wrapped sliver.
    const bool sliver = std::fabs(sweep) < SLIVER_SWEEP_RAD;
    if (curve.orientation() == CGAL::COUNTERCLOCKWISE) {
        if (sweep <= 0) sweep += 2 * std::numbers::pi;
    } else {
        if (sweep >= 0) sweep -= 2 * std::numbers::pi;
    }
    if (sliver || std::fabs(sweep) < SLIVER_SWEEP_RAD || std::fabs(sweep) > 2 * std::numbers::pi - SLIVER_SWEEP_RAD) {
        // Sliver regime: an x-monotone split within double rounding of a
        // vertical tangency. No double chord direction can land inside, so
        // the exact bisecting construction is used; its conditioning cost is
        // confined to these pieces. For parameter 1/2 this is the exact
        // angular midpoint, otherwise a chord interpolation of the two
        // halves, normalised onto the exact circle.
        auto middle = first_radial + last_radial;
        if (middle == ReachKernelVector(CGAL::NULL_VECTOR)) {
            middle = first_radial.perpendicular(curve.orientation());
        }
        middle = middle * CGAL::sqrt(circle.squared_radius() / middle.squared_length());
        if (CGAL::compare(u, ReachFT(1) / ReachFT(2)) == CGAL::EQUAL) {
            const auto point = center + middle;
            return ReachPoint(point.x(), point.y());
        }
        const bool first_half = CGAL::compare(u, ReachFT(1) / ReachFT(2)) != CGAL::LARGER;
        const ReachFT local = first_half ? ReachFT(2) * u : ReachFT(2) * u - ReachFT(1);
        const auto a = first_half ? first_radial : middle;
        const auto b = first_half ? middle : last_radial;
        const auto chord = (ReachFT(1) - local) * a + local * b;
        const auto point = center + chord * CGAL::sqrt(circle.squared_radius() / chord.squared_length());
        return ReachPoint(point.x(), point.y());
    }
    const double target_angle = start_angle + parameter * sweep;
    // The chord from the start to the point at target_angle is perpendicular
    // to the radius bisecting the two angles.
    const double chord_angle = (start_angle + target_angle) / 2 + std::numbers::pi / 2;
    const ReachKernelVector chord(ReachFT(std::cos(chord_angle)), ReachFT(std::sin(chord_angle)));
    const ReachFT along = ReachFT(-2) * CGAL::scalar_product(first_radial, chord) / chord.squared_length();
    const auto point = start + along * chord;
    // Membership of the exact sample in the trimmed piece, decided exactly on
    // directions. Outside the sliver regime the direction rounding is orders
    // of magnitude below the piece's extent, so a failure here is a defect.
    const auto direction = (point - center).direction();
    const auto first = first_radial.direction(), last = last_radial.direction();
    const bool inside = curve.orientation() == CGAL::COUNTERCLOCKWISE
        ? direction.counterclockwise_in_between(first, last)
        : direction.counterclockwise_in_between(last, first);
    if (!inside || direction == first || direction == last) {
        throw ReachableArrangementTopologyError("Boundary sample fell outside its piece.");
    }
    return ReachPoint(point.x(), point.y());
}
