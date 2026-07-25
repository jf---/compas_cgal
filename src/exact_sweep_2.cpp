#include "exact_sweep_2.h"

#include <iterator>
#include <string>
#include <utility>

#include <CGAL/number_utils.h>

namespace {

void require_positive(const ReachFT& value, const std::string& name)
{
    if (CGAL::compare(value, ReachFT(0)) != CGAL::LARGER) {
        throw ReachableDomainConstructionError(name + " must be positive.");
    }
}

ReachPolygon counterclockwise_polygon(ReachPolygon polygon)
{
    const CGAL::Orientation orientation = polygon.orientation();
    if (orientation == CGAL::COLLINEAR) {
        throw ReachableDomainConstructionError(
            "exact sweep polygon must enclose nonzero area.");
    }
    if (orientation == CGAL::CLOCKWISE) {
        polygon.reverse_orientation();
    }
    if (polygon.orientation() != CGAL::COUNTERCLOCKWISE) {
        throw ReachableDomainConstructionError(
            "exact sweep polygon did not normalize counterclockwise.");
    }
    return polygon;
}

ReachPolygon linear_polygon(const std::vector<ReachKernelPoint>& vertices)
{
    ReachPolygon polygon;
    for (std::size_t index = 0; index < vertices.size(); ++index) {
        polygon.push_back(ReachXCurve(
            vertices[index],
            vertices[(index + 1) % vertices.size()]));
    }
    return counterclockwise_polygon(std::move(polygon));
}

ReachFT lift_coordinate(const ReachPoint::CoordNT& coordinate)
{
    if (!coordinate.is_extended()) {
        return coordinate.a0();
    }
    return coordinate.a0()
        + coordinate.a1() * CGAL::sqrt(coordinate.root());
}

ReachKernelPoint kernel_point(const ReachPoint& point)
{
    return ReachKernelPoint(
        lift_coordinate(point.x()),
        lift_coordinate(point.y()));
}

void append_curve_polygon(
    ReachPolygon& polygon,
    const ReachCurve& curve)
{
    ReachTraits traits;
    traits.make_x_monotone_2_object()(
        curve,
        CGAL::dispatch_or_drop_output<ReachXCurve>(
            std::back_inserter(polygon)));
}

void append_sector_center_closure(
    ReachPolygon& sector,
    const ReachKernelPoint& outer_end,
    const ReachKernelPoint& center,
    const ReachKernelPoint& outer_start)
{
    sector.push_back(ReachXCurve(outer_end, center));
    sector.push_back(ReachXCurve(center, outer_start));
}

} // namespace

ReachPolygon reach_disk_polygon(
    const ReachKernelPoint& center,
    const ReachFT& radius)
{
    require_positive(radius, "exact disk radius");
    ReachTraits traits;
    ReachKernel::Circle_2 circle(
        center,
        radius * radius,
        CGAL::COUNTERCLOCKWISE);
    ReachPolygon polygon;
    traits.make_x_monotone_2_object()(
        ReachCurve(circle),
        CGAL::dispatch_or_drop_output<ReachXCurve>(
            std::back_inserter(polygon)));
    if (polygon.size() != 2) {
        throw ReachableDomainConstructionError(
            "exact disk boundary did not split into two x-monotone arcs.");
    }
    return counterclockwise_polygon(std::move(polygon));
}

std::vector<ReachPolygon> reach_capsule_parts(
    const ReachKernelPoint& start,
    const ReachKernelPoint& end,
    const ReachFT& radius)
{
    require_positive(radius, "exact capsule radius");
    const ReachKernelVector direction = end - start;
    if (direction == ReachKernelVector(CGAL::NULL_VECTOR)) {
        return {reach_disk_polygon(start, radius)};
    }
    const ReachFT length = CGAL::sqrt(direction.squared_length());
    const ReachKernelVector normal =
        direction.perpendicular(CGAL::COUNTERCLOCKWISE)
        * (radius / length);
    return {
        linear_polygon(
            {start - normal, end - normal, end + normal, start + normal}),
        reach_disk_polygon(start, radius),
        reach_disk_polygon(end, radius),
    };
}

std::vector<ReachPolygon> reach_arc_sweep_parts(
    const ReachXCurve& guide_arc,
    const ReachFT& tool_radius)
{
    require_positive(tool_radius, "exact arc-sweep tool radius");
    if (!guide_arc.is_circular()) {
        throw ReachableDomainConstructionError(
            "exact arc sweep requires a circular guide curve.");
    }
    const ReachKernel::Circle_2 guide_circle =
        guide_arc.supporting_circle();
    const ReachKernelPoint center = guide_circle.center();
    const ReachFT guide_radius =
        CGAL::sqrt(guide_circle.squared_radius());
    require_positive(guide_radius, "exact arc-sweep guide radius");
    const ReachKernelPoint start = kernel_point(guide_arc.source());
    const ReachKernelPoint end = kernel_point(guide_arc.target());
    const ReachKernelVector start_direction = start - center;
    const ReachKernelVector end_direction = end - center;
    const ReachFT outer_radius = guide_radius + tool_radius;
    const ReachKernelPoint outer_start =
        center + start_direction * (outer_radius / guide_radius);
    const ReachKernelPoint outer_end =
        center + end_direction * (outer_radius / guide_radius);
    const CGAL::Orientation orientation = guide_arc.orientation();
    ReachKernel::Circle_2 outer_circle(
        center,
        outer_radius * outer_radius,
        orientation);
    ReachPolygon sector;
    append_curve_polygon(
        sector,
        ReachCurve(
            outer_circle,
            ReachPoint(outer_start.x(), outer_start.y()),
            ReachPoint(outer_end.x(), outer_end.y())));

    const CGAL::Comparison_result radius_comparison =
        CGAL::compare(guide_radius, tool_radius);
    if (radius_comparison == CGAL::LARGER) {
        const ReachFT inner_radius = guide_radius - tool_radius;
        const ReachKernelPoint inner_start =
            center + start_direction * (inner_radius / guide_radius);
        const ReachKernelPoint inner_end =
            center + end_direction * (inner_radius / guide_radius);
        sector.push_back(ReachXCurve(outer_end, inner_end));
        const CGAL::Orientation reverse_orientation =
            orientation == CGAL::COUNTERCLOCKWISE
            ? CGAL::CLOCKWISE
            : CGAL::COUNTERCLOCKWISE;
        ReachKernel::Circle_2 inner_circle(
            center,
            inner_radius * inner_radius,
            reverse_orientation);
        append_curve_polygon(
            sector,
            ReachCurve(
                inner_circle,
                ReachPoint(inner_end.x(), inner_end.y()),
                ReachPoint(inner_start.x(), inner_start.y())));
        sector.push_back(ReachXCurve(inner_start, outer_start));
    }
    else if (radius_comparison == CGAL::EQUAL) {
        append_sector_center_closure(
            sector,
            outer_end,
            center,
            outer_start);
    }
    else {
        append_sector_center_closure(
            sector,
            outer_end,
            center,
            outer_start);
    }
    return {
        counterclockwise_polygon(std::move(sector)),
        reach_disk_polygon(start, tool_radius),
        reach_disk_polygon(end, tool_radius),
    };
}

ReachSet reach_join_parts(
    const std::vector<ReachPolygon>& polygons,
    const std::vector<ReachPolygonWithHoles>& polygons_with_holes)
{
    ReachSet result;
    result.join(
        polygons.begin(),
        polygons.end(),
        polygons_with_holes.begin(),
        polygons_with_holes.end());
    return result;
}

ReachSet reach_full_circle_sweep(
    const ReachKernelPoint& center,
    const ReachKernelVector& phase_vector,
    const ReachFT& tool_radius)
{
    require_positive(tool_radius, "exact full-circle tool radius");
    if (phase_vector == ReachKernelVector(CGAL::NULL_VECTOR)) {
        throw ReachableDomainConstructionError(
            "exact full-circle sweep requires a nonzero phase vector.");
    }
    const ReachFT guide_radius =
        CGAL::sqrt(phase_vector.squared_length());
    ReachSet result;
    result.insert(reach_disk_polygon(
        center,
        guide_radius + tool_radius));
    const CGAL::Comparison_result radius_comparison =
        CGAL::compare(guide_radius, tool_radius);
    if (radius_comparison == CGAL::LARGER) {
        ReachSet inner;
        inner.insert(reach_disk_polygon(
            center,
            guide_radius - tool_radius));
        result.difference(inner);
    }
    else if (radius_comparison == CGAL::EQUAL) {
        return result;
    }
    else {
        return result;
    }
    return result;
}
