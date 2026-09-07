#include "native_boundary_curve_2.h"

#include <CGAL/Boolean_set_operations_2/Gps_polygon_validation.h>
#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/function_objects.h>
#include <cmath>
#include <iterator>

namespace {
ReachKernelPoint point(const NativeBoundaryCurve2::XY& xy)
{
    if (!std::isfinite(xy[0]) || !std::isfinite(xy[1])) {
        throw InvalidNativeBoundaryCurveError("Boundary coordinates must be finite world-XY millimetres.");
    }
    return ReachKernelPoint(xy[0], xy[1]);
}
void distinct(const ReachKernelPoint& start, const ReachKernelPoint& end)
{
    if (start == end) {
        throw InvalidNativeBoundaryCurveError("Boundary primitive endpoints must be distinct.");
    }
}
}

NativeBoundaryCurve2 NativeBoundaryCurve2::line(const XY& start, const XY& end)
{
    const auto a = point(start), b = point(end);
    distinct(a, b);
    return NativeBoundaryCurve2(ReachCurve(a, b), ReachFT(0));
}

NativeBoundaryCurve2 NativeBoundaryCurve2::arc(
    const XY& start, const XY& end, const XY& center, bool counterclockwise)
{
    const auto a = point(start), b = point(end), authored_center = point(center);
    distinct(a, b);
    // This is the unique closest center satisfying equal endpoint radii.
    // Both endpoints remain the exact injected authored values.
    const auto fitted_center = CGAL::bisector(a, b).projection(authored_center);
    const ReachFT squared_radius = CGAL::squared_distance(fitted_center, a);
    if (CGAL::sign(squared_radius) != CGAL::POSITIVE) {
        throw InvalidNativeBoundaryCurveError("Fitted arc radius must be positive.");
    }
    const ReachKernel::Circle_2 circle(fitted_center, squared_radius,
        counterclockwise ? CGAL::COUNTERCLOCKWISE : CGAL::CLOCKWISE);
    if (!circle.has_on_boundary(a) || !circle.has_on_boundary(b)) {
        throw InvalidNativeBoundaryCurveError("Fitted arc lost exact endpoint incidence.");
    }
    return NativeBoundaryCurve2(ReachCurve(circle, ReachPoint(a.x(), a.y()), ReachPoint(b.x(), b.y())),
        CGAL::squared_distance(authored_center, fitted_center));
}

double NativeBoundaryCurve2::center_adjustment_mm() const
{
    return CGAL::to_double(CGAL::sqrt(adjustment_squared_));
}
NativeBoundaryCurve2::XY NativeBoundaryCurve2::center_mm() const
{
    if (!is_arc()) throw InvalidNativeBoundaryCurveError("A line has no circular center.");
    const auto center = curve_.supporting_circle().center();
    return {CGAL::to_double(center.x()), CGAL::to_double(center.y())};
}
double NativeBoundaryCurve2::radius_mm() const
{
    if (!is_arc()) throw InvalidNativeBoundaryCurveError("A line has no circular radius.");
    return CGAL::to_double(CGAL::sqrt(curve_.supporting_circle().squared_radius()));
}

NativeBoundary2::NativeBoundary2(std::vector<NativeBoundaryCurve2> curves)
    : curves_(std::move(curves))
{
    if (curves_.empty()) throw InvalidNativeBoundaryChainError("Boundary chain must not be empty.");
    ReachTraits traits;
    ReachPolygon polygon;
    for (std::size_t i = 0; i < curves_.size(); ++i) {
        if (!traits.equal_2_object()(curves_[i].end(), curves_[(i + 1) % curves_.size()].start())) {
            throw InvalidNativeBoundaryChainError("Boundary chain must close with exactly shared endpoints.");
        }
        std::vector<ReachXCurve> pieces;
        traits.make_x_monotone_2_object()(curves_[i].curve(),
            CGAL::dispatch_or_drop_output<ReachXCurve>(std::back_inserter(pieces)));
        for (const auto& piece : pieces) {
            polygon.push_back(piece);
            cycle_.curves.push_back({piece, {}});
        }
    }
    if (!CGAL::is_closed_polygon(polygon, traits) || !CGAL::is_simple_polygon(polygon, traits)) {
        throw InvalidNativeBoundaryChainError("Boundary chain must be closed and simple without crossings or retracing.");
    }
    cycle_.orientation = polygon.orientation();
    if (cycle_.orientation == CGAL::COLLINEAR) {
        throw InvalidNativeBoundaryChainError("Boundary chain must enclose nonzero area.");
    }
    if (cycle_.orientation == CGAL::CLOCKWISE) polygon.reverse_orientation();
    design_ = ReachSet(polygon);
}

ExactRegion2 NativeBoundary2::design_region() const
{
    return ExactRegion2::build(design_, ExactRegionRole2::Design, "native-line-arc-design");
}
