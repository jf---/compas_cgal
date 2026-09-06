#include "remaining_material_2.h"

#include "coverage_2.h"
#include "exact_sweep_2.h"

#include <cmath>
#include <string>
#include <utility>

ExactRegion2 remaining_material(
    const ExactRegion2& target,
    Eigen::Ref<const compas::RowMatrixXd> circles,
    Eigen::Ref<const compas::RowMatrixXd> segments,
    Eigen::Ref<const compas::RowMatrixXd> disks,
    double tool_radius)
{
    if (target.role() != ExactRegionRole2::Design
        && target.role() != ExactRegionRole2::ReachableMaterial) {
        throw CoverageTransitionError(
            "remaining material requires a design or reachable-material target.");
    }
    if (circles.cols() != 3 || segments.cols() != 4 || disks.cols() != 2) {
        throw InvalidCoverageGeometryError(
            "remaining material requires circles Nx3, segments Nx4, disks Nx2.");
    }
    if (!circles.allFinite() || !segments.allFinite() || !disks.allFinite()
        || !std::isfinite(tool_radius)) {
        throw InvalidCoverageGeometryError(
            "remaining material requires finite motion coordinates and radii.");
    }
    const ReachFT radius(tool_radius);
    if (CGAL::sign(radius) != CGAL::POSITIVE) {
        throw InvalidCoverageGeometryError("cutter radius must be positive.");
    }
    // Validate the complete input before an empty residual can short-circuit
    // geometry. Otherwise an invalid tail could silently qualify as covered.
    for (Eigen::Index row = 0; row < circles.rows(); ++row) {
        if (CGAL::sign(ReachFT(circles(row, 2))) != CGAL::POSITIVE) {
            throw InvalidCoverageGeometryError("circle guide radius must be positive.");
        }
    }
    for (Eigen::Index row = 0; row < segments.rows(); ++row) {
        if (ReachKernelPoint(segments(row, 0), segments(row, 1))
            == ReachKernelPoint(segments(row, 2), segments(row, 3))) {
            throw InvalidCoverageGeometryError(
                "segment endpoints must be distinct; use a disk for a plunge.");
        }
    }
    try {
        ReachSet remaining(target.set());
        for (Eigen::Index row = 0;
             row < circles.rows() && !remaining.is_empty(); ++row) {
            remaining.difference(reach_full_circle_sweep(
                ReachKernelPoint(circles(row, 0), circles(row, 1)),
                ReachKernelVector(circles(row, 2), 0), radius));
        }
        for (Eigen::Index row = 0;
             row < segments.rows() && !remaining.is_empty(); ++row) {
            remaining.difference(reach_join_parts(
                reach_capsule_parts(
                    ReachKernelPoint(segments(row, 0), segments(row, 1)),
                    ReachKernelPoint(segments(row, 2), segments(row, 3)), radius),
                {}));
        }
        for (Eigen::Index row = 0;
             row < disks.rows() && !remaining.is_empty(); ++row) {
            remaining.difference(reach_disk_polygon(
                ReachKernelPoint(disks(row, 0), disks(row, 1)), radius));
        }
        return ExactRegion2::build(
            std::move(remaining), ExactRegionRole2::CoverageResidual,
            "full-motion remaining material");
    }
    catch (const std::exception& error) {
        throw CoverageTransitionError(
            std::string("exact remaining material failed: ") + error.what());
    }
}
