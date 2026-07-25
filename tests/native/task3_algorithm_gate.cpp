#include "exact_region_2.h"
#include "exact_sweep_2.h"

#include <stdexcept>

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

void exact_sweep_gate()
{
    const std::vector<ReachPolygon> capsule = reach_capsule_parts(
        ReachKernelPoint(0, 0),
        ReachKernelPoint(3, 4),
        ReachFT(5));
    const ReachSet capsule_set = reach_join_parts(capsule, {});
    require(
        capsule_set.oriented_side(ReachPoint(ReachFT(-2.5), ReachFT(5)))
            != CGAL::ON_NEGATIVE_SIDE,
        "sqrt(13) capsule lost side boundary");
}

} // namespace

int main()
{
    exact_region_storage_gate();
    exact_sweep_gate();
}
