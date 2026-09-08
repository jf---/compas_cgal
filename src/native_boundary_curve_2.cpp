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
    const auto a = point(start), b = point(end);
    point(center);  // validates finiteness; the value itself is fitted below
    distinct(a, b);
    // The fit is a rational construction: the unique closest centre with equal
    // endpoint radii is the projection of the supplied centre onto the
    // perpendicular bisector. It is computed in exact rational arithmetic and
    // injected as leaves, so every later expression over this circle (the
    // arrangement, the polygon set, the medial query) stays shallow. Left as
    // an expression DAG, CORE re-bounds and re-approximates it under every
    // decision downstream: measured 2026-09-08, 25 ms per competitor line.
    // Both endpoints remain the exact injected authored values.
    using Rational = CORE::BigRat;
    const Rational ax(start[0]), ay(start[1]), bx(end[0]), by(end[1]), cx(center[0]), cy(center[1]);
    const Rational dx = bx - ax, dy = by - ay;
    const Rational mx = (ax + bx) / 2, my = (ay + by) / 2;
    const Rational along = ((cx - mx) * dx + (cy - my) * dy) / (dx * dx + dy * dy);
    const Rational fx = cx - along * dx, fy = cy - along * dy;
    const Rational squared_radius = (ax - fx) * (ax - fx) + (ay - fy) * (ay - fy);
    if (squared_radius <= 0) {
        throw InvalidNativeBoundaryCurveError("Fitted arc radius must be positive.");
    }
    const ReachKernelPoint fitted_center((ReachFT(fx)), (ReachFT(fy)));
    const ReachKernel::Circle_2 circle(fitted_center, ReachFT(squared_radius),
        counterclockwise ? CGAL::COUNTERCLOCKWISE : CGAL::CLOCKWISE);
    // |a - f|^2 is the squared radius by definition and |b - f| = |a - f|
    // because f lies on the bisector: endpoint incidence holds by construction
    // and is witnessed by the import tests, never re-decided here.
    const Rational adjustment_squared = (cx - fx) * (cx - fx) + (cy - fy) * (cy - fy);
    return NativeBoundaryCurve2(ReachCurve(circle, ReachPoint(a.x(), a.y()), ReachPoint(b.x(), b.y())),
        ReachFT(adjustment_squared));
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
    design_ = std::make_shared<ReachSet>(polygon);
}

ExactRegion2 NativeBoundary2::design_region() const
{
    return ExactRegion2::build(design_, ExactRegionRole2::Design, "native-line-arc-design", {});
}
