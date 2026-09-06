#include "stock_2.h"

#include <cmath>

bool Stock2::can_remove_circle(
    const std::array<double, 3>& previous,
    const std::array<double, 3>& current,
    const std::array<double, 3>& next,
    Eigen::Ref<const compas::RowMatrixXd> connector_samples,
    double tool_radius) const
{
    if (connector_samples.cols() != 2 || connector_samples.rows() < 2) {
        throw InvalidCircleRemovalInputError(
            "circle removal requires at least two world-XY connector samples.");
    }
    if (!connector_samples.allFinite() || !std::isfinite(tool_radius)) {
        throw InvalidCircleRemovalInputError(
            "circle removal requires finite connector coordinates and cutter radius.");
    }
    if (CGAL::sign(Epeck::FT(tool_radius)) != CGAL::POSITIVE) {
        throw InvalidCircleRemovalInputError("circle removal requires a positive cutter radius.");
    }
    // Validate every input before an empty clipped circle can authorize removal.
    for (const auto* circle : {&previous, &current, &next}) {
        for (double value : *circle) {
            if (!std::isfinite(value)) {
                throw InvalidCircleRemovalInputError("circle removal requires finite circle coordinates and radius.");
            }
        }
        if (CGAL::sign(Epeck::FT((*circle)[2])) == CGAL::NEGATIVE) {
            throw InvalidCircleRemovalInputError("circle removal requires nonnegative guide radii.");
        }
    }

    Stock2 lost = clone();
    lost.intersect_circle_sweep(current[0], current[1], current[2], tool_radius);
    lost.subtract_circle_sweep(previous[0], previous[1], previous[2], tool_radius);
    lost.subtract_circle_sweep(next[0], next[1], next[2], tool_radius);
    for (Eigen::Index row = 1;
         row < connector_samples.rows() && !lost.is_empty(); ++row) {
        const EPoint start(connector_samples(row - 1, 0), connector_samples(row - 1, 1));
        const EPoint end(connector_samples(row, 0), connector_samples(row, 1));
        if (start == end) {
            continue;  // A zero-progress fan contributes no additional motion.
        }
        lost.subtract_capsule_quad(
            connector_samples(row - 1, 0), connector_samples(row - 1, 1),
            connector_samples(row, 0), connector_samples(row, 1), tool_radius);
    }
    // Quad sweeps under-cover the actual connectors. Empty therefore proves
    // preservation; nonempty merely retains the circle, without claiming a gap.
    return lost.is_empty();
}
