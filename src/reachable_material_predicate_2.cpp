#include "reachable_material_predicate_2.h"

#include "exact_sweep_2.h"
#include "reachable_arrangement_2.h"
#include "reachable_errors_2.h"

#include <cmath>
#include <iterator>
#include <memory>
#include <mutex>
#include <utility>
#include <variant>
#include <vector>

#include <CGAL/Arr_trapezoid_ric_point_location.h>
#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_2.h>

namespace {

using ReachPointLocation2 =
    CGAL::Arr_trapezoid_ric_point_location<ReachSet::Arrangement_2>;

struct ReachCurveDistance2 {
    ReachXCurve curve;
    ReachKernelPoint source;
    ReachKernelPoint target;
    bool circular;
    ReachKernelPoint circle_center;
    ReachFT circle_squared_radius;
};

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

void append_forbidden_capsule_parts(
    const CanonicalReachRing2& ring,
    const ReachFT& radius,
    std::vector<ReachPolygon>& parts)
{
    for (std::size_t index = 0; index < ring.points.size(); ++index) {
        std::vector<ReachPolygon> capsule = reach_capsule_parts(
            ring.points[index],
            ring.points[(index + 1) % ring.points.size()],
            radius);
        parts.insert(
            parts.end(),
            std::make_move_iterator(capsule.begin()),
            std::make_move_iterator(capsule.end()));
    }
}

ReachCurveDistance2 curve_distance_data(const ReachXCurve& curve)
{
    if (!curve.is_linear() && !curve.is_circular()) {
        throw ReachableMaterialPredicateGeometryError(
            "center boundary curve is neither linear nor circular");
    }
    const ReachKernelPoint source = kernel_point(curve.source());
    const ReachKernelPoint target = kernel_point(curve.target());
    if (curve.is_linear()) {
        if (source == target) {
            throw ReachableMaterialPredicateGeometryError(
                "center boundary line curve is degenerate");
        }
        return {
            curve,
            source,
            target,
            false,
            ReachKernelPoint(),
            ReachFT(0),
        };
    }
    const ReachKernel::Circle_2 circle = curve.supporting_circle();
    if (CGAL::compare(circle.squared_radius(), ReachFT(0))
        != CGAL::LARGER) {
        throw ReachableMaterialPredicateGeometryError(
            "center boundary circular curve has non-positive radius");
    }
    return {
        curve,
        source,
        target,
        true,
        circle.center(),
        circle.squared_radius(),
    };
}

bool radial_distance_within_radius(
    const ReachFT& center_distance_squared,
    const ReachFT& circle_squared_radius,
    const ReachFT& radius_squared)
{
    const ReachFT sum = circle_squared_radius + radius_squared;
    const ReachFT product =
        ReachFT(4) * circle_squared_radius * radius_squared;
    if (CGAL::compare(center_distance_squared, sum) == CGAL::LARGER) {
        const ReachFT excess = center_distance_squared - sum;
        return CGAL::compare(excess * excess, product) != CGAL::LARGER;
    }
    if (CGAL::compare(circle_squared_radius, radius_squared)
        != CGAL::LARGER) {
        return true;
    }
    const ReachFT deficit = sum - center_distance_squared;
    return CGAL::compare(deficit * deficit, product) != CGAL::LARGER;
}

bool curve_within_radius(
    const ReachKernelPoint& query,
    const ReachCurveDistance2& data,
    const ReachFT& radius)
{
    const ReachFT squared_radius = radius * radius;
    if (!data.circular) {
        return CGAL::compare(
                   CGAL::squared_distance(
                       query,
                       ReachKernel::Segment_2(data.source, data.target)),
                   squared_radius)
            != CGAL::LARGER;
    }

    const ReachFT center_distance_squared =
        CGAL::squared_distance(query, data.circle_center);
    if (CGAL::is_zero(center_distance_squared)) {
        return CGAL::compare(
                   data.circle_squared_radius,
                   squared_radius)
            != CGAL::LARGER;
    }
    if (!radial_distance_within_radius(
            center_distance_squared,
            data.circle_squared_radius,
            squared_radius)) {
        return false;
    }

    const ReachKernelVector radial = query - data.circle_center;
    const ReachKernelPoint projection =
        data.circle_center
        + radial * CGAL::sqrt(
            data.circle_squared_radius / center_distance_squared);
    const ReachPoint projected(projection.x(), projection.y());
    if (data.curve.is_in_x_range(projected)
        && data.curve.point_position(projected) == CGAL::EQUAL) {
        return true;
    }
    return CGAL::compare(
               CGAL::squared_distance(query, data.source),
               squared_radius)
            != CGAL::LARGER
        || CGAL::compare(
               CGAL::squared_distance(query, data.target),
               squared_radius)
            != CGAL::LARGER;
}

void append_boundary_curves(
    const ReachPolygon& boundary,
    std::vector<ReachCurveDistance2>& curves)
{
    for (auto curve = boundary.curves_begin();
         curve != boundary.curves_end();
         ++curve) {
        curves.push_back(curve_distance_data(*curve));
    }
}

} // namespace

bool reach_curve_within_radius(
    const ReachKernelPoint& query,
    const ReachXCurve& curve,
    const ReachFT& radius)
{
    return curve_within_radius(
        query,
        curve_distance_data(curve),
        radius);
}

struct ReachableMaterialPredicateStorage2 {
    ReachableMaterialPredicateStorage2(
        ReachSet design_value,
        ReachSet center_value,
        std::vector<ReachCurveDistance2> boundary_curve_values,
        ReachFT radius_value,
        ReachableDomainBuildAudit2 audit_value)
        : design(std::move(design_value))
        , center(std::move(center_value))
        , design_point_location(design.arrangement())
        , center_point_location(center.arrangement())
        , boundary_curves(std::move(boundary_curve_values))
        , radius(std::move(radius_value))
        , audit(std::move(audit_value))
    {
    }

    ReachSet design;
    ReachSet center;
    mutable ReachPointLocation2 design_point_location;
    mutable ReachPointLocation2 center_point_location;
    mutable std::mutex design_point_location_mutex;
    mutable std::mutex center_point_location_mutex;
    std::vector<ReachCurveDistance2> boundary_curves;
    ReachFT radius;
    ReachableDomainBuildAudit2 audit;
};

bool region_contains(
    ReachPointLocation2& point_location,
    std::mutex& point_location_mutex,
    const ReachPoint& query)
{
    const std::lock_guard<std::mutex> lock(point_location_mutex);
    using Arrangement = ReachSet::Arrangement_2;
    const auto located = point_location.locate(query);
    const auto* face =
        std::get_if<Arrangement::Face_const_handle>(&located);
    return face == nullptr || (*face)->contained();
}

ReachableMaterialPredicate2 ReachableMaterialPredicate2::build(
    Eigen::Ref<const compas::RowMatrixXd> design_boundary,
    const std::vector<compas::RowMatrixXd>& holes,
    double tool_radius)
{
    CanonicalReachInput2 input = canonical_reach_input(
        design_boundary,
        holes,
        tool_radius);
    validate_canonical_reach_input(input);
    ReachPolygonWithHoles design_polygon =
        reachable_design_polygon(input);

    std::vector<ReachPolygon> forbidden_parts;
    append_forbidden_capsule_parts(
        input.outer,
        input.radius,
        forbidden_parts);
    for (const CanonicalReachRing2& hole : input.holes) {
        append_forbidden_capsule_parts(
            hole,
            input.radius,
            forbidden_parts);
    }
    ReachSet design(design_polygon);
    ReachSet center = design;
    center.difference(reach_join_parts(forbidden_parts, {}));

    std::vector<ReachPolygonWithHoles> center_components;
    center.polygons_with_holes(
        std::back_inserter(center_components));
    if (center_components.empty()) {
        throw PocketNotMachinableError(
            "Phase 1 center domain is empty for the declared tool radius");
    }
    if (center_components.size() != 1) {
        throw PocketNotMachinableError(
            "Phase 1 requires one connected center domain");
    }

    ReachableDomainBuildAudit2 audit;
    ++audit.geometry_passes;
    ++audit.center_extractions;
    ++audit.center_predicate_constructions;
    ++audit.design_point_locators;
    ++audit.center_point_locators;
    std::vector<ReachCurveDistance2> boundary_curves;
    append_boundary_curves(
        center_components.front().outer_boundary(),
        boundary_curves);
    for (auto hole = center_components.front().holes_begin();
         hole != center_components.front().holes_end();
         ++hole) {
        append_boundary_curves(*hole, boundary_curves);
    }
    return ReachableMaterialPredicate2(
        std::make_shared<const ReachableMaterialPredicateStorage2>(
            std::move(design),
            std::move(center),
            std::move(boundary_curves),
            input.radius,
            audit));
}

ReachableMaterialPredicate2::ReachableMaterialPredicate2(
    std::shared_ptr<const ReachableMaterialPredicateStorage2> storage)
    : storage_(std::move(storage))
{
}

bool ReachableMaterialPredicate2::contains(double x, double y) const
{
    if (!std::isfinite(x) || !std::isfinite(y)) {
        throw InvalidReachableDomainInputError(
            "reachable-material query coordinates must be finite binary64 values");
    }
    const ReachPoint query_point{ReachFT(x), ReachFT(y)};
    if (!region_contains(
            storage_->design_point_location,
            storage_->design_point_location_mutex,
            query_point)) {
        return false;
    }
    if (region_contains(
            storage_->center_point_location,
            storage_->center_point_location_mutex,
            query_point)) {
        return true;
    }
    const ReachKernelPoint query{ReachFT(x), ReachFT(y)};
    for (const ReachCurveDistance2& curve : storage_->boundary_curves) {
        if (curve_within_radius(query, curve, storage_->radius)) {
            return true;
        }
    }
    return false;
}

const ReachableDomainBuildAudit2&
ReachableMaterialPredicate2::build_audit_for_native_gate() const
{
    return storage_->audit;
}
