#include "reachable_domain_2.h"

#include "exact_sweep_2.h"
#include "reachable_arrangement_2.h"
#include "reachable_errors_2.h"
#include "reachable_input_2.h"

#include <iterator>
#include <string>
#include <utility>
#include <vector>

#include <CGAL/number_utils.h>

namespace {

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

void append_boundary_sweep_parts(
    const ReachPolygon& boundary,
    const ReachFT& radius,
    std::vector<ReachPolygon>& parts)
{
    for (auto curve = boundary.curves_begin();
         curve != boundary.curves_end();
         ++curve) {
        std::vector<ReachPolygon> curve_parts =
            curve->is_linear()
            ? reach_capsule_parts(
                kernel_point(curve->source()),
                kernel_point(curve->target()),
                radius)
            : reach_arc_sweep_parts(*curve, radius);
        parts.insert(
            parts.end(),
            std::make_move_iterator(curve_parts.begin()),
            std::make_move_iterator(curve_parts.end()));
    }
}

ReachSet build_reachable_material_once(
    const ReachPolygonWithHoles& center,
    const ReachFT& radius,
    ReachableDomainBuildAudit2& audit)
{
    std::vector<ReachPolygon> parts;
    append_boundary_sweep_parts(
        center.outer_boundary(),
        radius,
        parts);
    for (auto hole = center.holes_begin();
         hole != center.holes_end();
         ++hole) {
        append_boundary_sweep_parts(*hole, radius, parts);
    }
    ReachSet material = reach_join_parts(parts, {center});
    ++audit.material_batch_unions;
    return material;
}

} // namespace

ReachableDomain2::ReachableDomain2(
    Eigen::Ref<const compas::RowMatrixXd> design_boundary,
    const std::vector<compas::RowMatrixXd>& holes,
    double tool_radius)
    : ReachableDomain2(
        build_state(design_boundary, holes, tool_radius))
{
}

ReachableDomain2::ReachableDomain2(State state)
    : state_(std::move(state))
{
}

ReachableDomain2::State ReachableDomain2::build_state(
    Eigen::Ref<const compas::RowMatrixXd> design_boundary,
    const std::vector<compas::RowMatrixXd>& holes,
    double tool_radius)
{
    CanonicalReachInput2 input =
        canonical_reach_input(
            design_boundary,
            holes,
            tool_radius);
    ReachableArrangement2 selected =
        build_reachable_arrangement(std::move(input));
    ++selected.audit.geometry_passes;

    ReachSet design(selected.design_polygon);
    ReachSet center(selected.center_polygon);
    ReachSet material = build_reachable_material_once(
        selected.center_polygon,
        selected.input.radius,
        selected.audit);

    ++selected.audit.subset_decisions;
    const bool subset = reach_exact_subset(material, design);
    if (!subset) {
        throw ReachableMaterialContainmentError(
            "exact reachable material is not contained in the design");
    }

    ReachSet residual(design);
    residual.difference(material);
    ++selected.audit.residual_differences;

    ReachableDomainCertificate2 certificate =
        build_reachable_certificate(selected, subset);
    const std::string input_recipe =
        selected.input.recipe_record;
    const std::string design_recipe = reach_tagged_record(
        "exact-region-design-v2",
        {input_recipe});
    const std::string center_recipe = reach_tagged_record(
        "exact-region-center-domain-v2",
        {input_recipe});
    const std::string material_recipe = reach_tagged_record(
        "exact-region-reachable-material-v2",
        {center_recipe});
    const std::string residual_recipe = reach_tagged_record(
        "exact-region-unreachable-residual-v2",
        {design_recipe, material_recipe});

    return State{
        ExactRegion2::build(
            std::move(design),
            ExactRegionRole2::Design,
            design_recipe),
        ExactRegion2::build(
            std::move(center),
            ExactRegionRole2::CenterDomain,
            center_recipe),
        ExactRegion2::build(
            std::move(material),
            ExactRegionRole2::ReachableMaterial,
            material_recipe),
        ExactRegion2::build(
            std::move(residual),
            ExactRegionRole2::UnreachableResidual,
            residual_recipe),
        std::move(certificate),
        selected.audit,
    };
}

ExactRegion2 ReachableDomain2::design_region() const
{
    return state_.design.clone();
}

ExactRegion2 ReachableDomain2::center_domain() const
{
    return state_.center_domain.clone();
}

ExactRegion2 ReachableDomain2::reachable_material() const
{
    return state_.reachable_material.clone();
}

ExactRegion2 ReachableDomain2::unreachable_residual() const
{
    return state_.unreachable_residual.clone();
}

ReachableDomainCertificate2 ReachableDomain2::certificate() const
{
    return state_.certificate;
}

const ReachableDomainBuildAudit2&
ReachableDomain2::build_audit_for_native_gate() const
{
    return state_.audit;
}
