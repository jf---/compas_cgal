#include "exact_stock_region_2.h"

#include "exact_sweep_2.h"

#include <iterator>
#include <utility>
#include <vector>

#include <CGAL/number_utils.h>

namespace {

ReachFT lift_rational(const Epeck::FT& value)
{
    return ReachFT(CGAL::exact(value));
}

ReachFT lift_coordinate(
    const GpsPoint::CoordNT& coordinate)
{
    const ReachFT base =
        lift_rational(coordinate.a0());
    if (!coordinate.is_extended()) {
        return base;
    }
    return base
        + lift_rational(coordinate.a1())
            * CGAL::sqrt(
                lift_rational(
                    coordinate.root()));
}

ReachPoint lift_point(const GpsPoint& point)
{
    return ReachPoint(
        lift_coordinate(point.x()),
        lift_coordinate(point.y()));
}

ReachXCurve lift_curve(const GpsXCurve& curve)
{
    const ReachPoint source =
        lift_point(curve.source());
    const ReachPoint target =
        lift_point(curve.target());
    if (curve.is_linear()) {
        const Epeck::Line_2 support =
            curve.supporting_line();
        return ReachXCurve(
            ReachKernel::Line_2(
                lift_rational(support.a()),
                lift_rational(support.b()),
                lift_rational(support.c())),
            source,
            target);
    }
    if (!curve.is_circular()) {
        throw ExactStockRegionLiftError(
            "stock boundary is neither linear nor circular.");
    }
    const Epeck::Circle_2 support =
        curve.supporting_circle();
    return ReachXCurve(
        ReachKernel::Circle_2(
            ReachKernelPoint(
                lift_rational(
                    support.center().x()),
                lift_rational(
                    support.center().y())),
            lift_rational(
                support.squared_radius()),
            support.orientation()),
        source,
        target,
        curve.orientation());
}

ReachPolygon lift_polygon(
    const GpsPolygon& polygon)
{
    ReachPolygon result;
    for (auto curve = polygon.curves_begin();
         curve != polygon.curves_end();
         ++curve) {
        result.push_back(lift_curve(*curve));
    }
    if (result.size() != polygon.size()
        || result.orientation()
            != polygon.orientation()) {
        throw ExactStockRegionLiftError(
            "lifted stock boundary changed exact topology.");
    }
    return result;
}

ReachPolygonWithHoles lift_component(
    const GpsPolygonWithHoles& component)
{
    std::vector<ReachPolygon> holes;
    holes.reserve(
        static_cast<std::size_t>(
            std::distance(
                component.holes_begin(),
                component.holes_end())));
    for (auto hole = component.holes_begin();
         hole != component.holes_end();
         ++hole) {
        holes.push_back(lift_polygon(*hole));
    }
    return ReachPolygonWithHoles(
        lift_polygon(
            component.outer_boundary()),
        holes.begin(),
        holes.end());
}

} // namespace

ReachSet lift_exact_stock_region(
    const Stock2& stock)
{
    std::vector<GpsPolygonWithHoles> components;
    stock.set().polygons_with_holes(
        std::back_inserter(components));
    std::vector<ReachPolygonWithHoles> lifted;
    lifted.reserve(components.size());
    for (const GpsPolygonWithHoles& component :
         components) {
        lifted.push_back(
            lift_component(component));
    }
    return reach_join_parts({}, lifted);
}
