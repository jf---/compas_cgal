#include "exact_disk_region_2.h"

Gps build_exact_disk_union_region_2(
    const std::vector<EPoint>& centers,
    const Epeck::FT& radius)
{
    if (centers.empty()) {
        throw ExactDiskRegionEmptyCentersError(
            "exact disk union requires at least one center");
    }
    if (CGAL::sign(radius) != CGAL::POSITIVE) {
        throw ExactDiskRegionRadiusError(
            "exact disk union radius must be positive");
    }
    std::vector<GpsPolygon> disks;
    disks.reserve(centers.size());
    for (const EPoint& center : centers) {
        disks.push_back(disk_polygon(center, radius));
    }
    Gps region;
    region.join(disks.begin(), disks.end());
    return region;
}
