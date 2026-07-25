#include "exact_region_2.h"
#include "exact_sweep_2.h"

#include <iterator>
#include <stdexcept>
#include <vector>

#include <CGAL/iterator.h>
#include <CGAL/number_utils.h>

namespace {

void require(bool condition, const char* message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void exact_region_storage_gate()
{
    ReachSet disk;
    disk.insert(reach_disk_polygon(ReachKernelPoint(0, 0), ReachFT(2)));
    const ExactRegion2 original = ExactRegion2::build(
        std::move(disk),
        ExactRegionRole2::ReachableMaterial,
        "native-gate-region");
    const ExactRegion2 clone = original.clone();
    require(
        original.shares_storage_with_for_audit(clone),
        "read-only clone deep-copied exact storage");
    require(clone.contains(2.0, 0.0), "closed disk lost exact boundary");
}

void require_counterclockwise_parts(
    const std::vector<ReachPolygon>& parts,
    const char* message)
{
    for (const ReachPolygon& part : parts) {
        require(part.orientation() == CGAL::COUNTERCLOCKWISE, message);
    }
}

void irrational_capsule_gate()
{
    const ReachFT root_thirteen = CGAL::sqrt(ReachFT(13));
    const std::vector<ReachPolygon> capsule = reach_capsule_parts(
        ReachKernelPoint(0, 0),
        ReachKernelPoint(2, 3),
        ReachFT(5));
    require_counterclockwise_parts(
        capsule,
        "irrational capsule emitted non-CCW standalone polygon");
    const ReachSet capsule_set = reach_join_parts(capsule, {});
    const ReachPoint irrational_side_boundary(
        ReachFT(1) - ReachFT(15) / root_thirteen,
        ReachFT(3) / ReachFT(2) + ReachFT(10) / root_thirteen);
    require(
        capsule_set.oriented_side(irrational_side_boundary)
            == CGAL::ON_ORIENTED_BOUNDARY,
        "sqrt(13) capsule lost exact irrational side boundary");
}

ReachXCurve quarter_guide_arc(
    CGAL::Orientation orientation,
    const ReachFT& radius)
{
    const ReachKernelPoint center(0, 0);
    const ReachKernelPoint start(radius, ReachFT(0));
    const ReachKernelPoint end = orientation == CGAL::COUNTERCLOCKWISE
        ? ReachKernelPoint(ReachFT(0), radius)
        : ReachKernelPoint(ReachFT(0), -radius);
    const ReachKernel::Circle_2 circle(
        center,
        radius * radius,
        orientation);
    const ReachCurve curve(
        circle,
        ReachPoint(start.x(), start.y()),
        ReachPoint(end.x(), end.y()));
    std::vector<ReachXCurve> pieces;
    ReachTraits().make_x_monotone_2_object()(
        curve,
        CGAL::dispatch_or_drop_output<ReachXCurve>(
            std::back_inserter(pieces)));
    require(
        pieces.size() == 1,
        "quarter guide arc did not remain x-monotone");
    return pieces.front();
}

void arc_sweep_case(
    CGAL::Orientation orientation,
    const ReachFT& guide_radius,
    const ReachFT& tool_radius,
    CGAL::Oriented_side expected_center_side)
{
    const ReachXCurve guide =
        quarter_guide_arc(orientation, guide_radius);
    const std::vector<ReachPolygon> parts =
        reach_arc_sweep_parts(guide, tool_radius);
    require(parts.size() == 3, "arc sweep lost primitive parts");
    require_counterclockwise_parts(
        parts,
        "arc sweep emitted non-CCW standalone polygon");
    const ReachSet sweep = reach_join_parts(parts, {});
    require(
        reach_component_count(sweep) == 1,
        "arc sweep did not produce one connected component");
    require(
        sweep.oriented_side(ReachPoint(ReachFT(0), ReachFT(0)))
            == expected_center_side,
        "arc sweep center topology contradicts radius branch");
    require(
        sweep.oriented_side(
            ReachPoint(guide_radius + tool_radius, ReachFT(0)))
            == CGAL::ON_ORIENTED_BOUNDARY,
        "arc sweep lost exact outer endpoint boundary");
}

void arc_sweep_orientation_and_radius_gate()
{
    arc_sweep_case(
        CGAL::COUNTERCLOCKWISE,
        ReachFT(2),
        ReachFT(1),
        CGAL::ON_NEGATIVE_SIDE);
    arc_sweep_case(
        CGAL::CLOCKWISE,
        ReachFT(2),
        ReachFT(1),
        CGAL::ON_NEGATIVE_SIDE);
    arc_sweep_case(
        CGAL::COUNTERCLOCKWISE,
        ReachFT(1),
        ReachFT(1),
        CGAL::ON_ORIENTED_BOUNDARY);
    arc_sweep_case(
        CGAL::CLOCKWISE,
        ReachFT(1),
        ReachFT(1),
        CGAL::ON_ORIENTED_BOUNDARY);
    arc_sweep_case(
        CGAL::COUNTERCLOCKWISE,
        ReachFT(1),
        ReachFT(2),
        CGAL::ON_POSITIVE_SIDE);
    arc_sweep_case(
        CGAL::CLOCKWISE,
        ReachFT(1),
        ReachFT(2),
        CGAL::ON_POSITIVE_SIDE);
}

void full_circle_case(
    const ReachKernelVector& phase_vector,
    const ReachFT& tool_radius,
    const ReachFT& outer_radius,
    std::size_t expected_hole_count,
    CGAL::Oriented_side expected_center_side)
{
    const ReachSet sweep = reach_full_circle_sweep(
        ReachKernelPoint(0, 0),
        phase_vector,
        tool_radius);
    std::vector<ReachPolygonWithHoles> components;
    sweep.polygons_with_holes(std::back_inserter(components));
    require(
        components.size() == 1,
        "full-circle sweep did not produce one component");
    require(
        static_cast<std::size_t>(std::distance(
            components.front().holes_begin(),
            components.front().holes_end()))
            == expected_hole_count,
        "full-circle sweep hole count contradicts radius branch");
    require(
        sweep.oriented_side(ReachPoint(ReachFT(0), ReachFT(0)))
            == expected_center_side,
        "full-circle sweep center topology contradicts radius branch");
    require(
        sweep.oriented_side(ReachPoint(outer_radius, ReachFT(0)))
            == CGAL::ON_ORIENTED_BOUNDARY,
        "full-circle sweep lost exact outer boundary");
}

void full_circle_radius_gate()
{
    const ReachSet annulus = reach_full_circle_sweep(
        ReachKernelPoint(0, 0),
        ReachKernelVector(2, 0),
        ReachFT(1));
    require(
        annulus.oriented_side(ReachPoint(ReachFT(1), ReachFT(0)))
            == CGAL::ON_ORIENTED_BOUNDARY,
        "full-circle annulus lost exact inner boundary");
    full_circle_case(
        ReachKernelVector(2, 0),
        ReachFT(1),
        ReachFT(3),
        1,
        CGAL::ON_NEGATIVE_SIDE);
    full_circle_case(
        ReachKernelVector(1, 0),
        ReachFT(1),
        ReachFT(2),
        0,
        CGAL::ON_POSITIVE_SIDE);
    full_circle_case(
        ReachKernelVector(1, 0),
        ReachFT(2),
        ReachFT(3),
        0,
        CGAL::ON_POSITIVE_SIDE);
}

} // namespace

int main()
{
    exact_region_storage_gate();
    irrational_capsule_gate();
    arc_sweep_orientation_and_radius_gate();
    full_circle_radius_gate();
}
