#pragma once

#include "reachable_input_2.h"

#include <vector>

class CutterCentreDomain2 {
public:
    static CutterCentreDomain2 build(
        Eigen::Ref<const compas::RowMatrixXd> design_boundary,
        const std::vector<compas::RowMatrixXd>& holes,
        double tool_radius);

    bool contains(double x, double y) const;

    // The closed segment [a, b] lies in the domain: a is contained and the
    // segment keeps squared distance >= r^2 from every outer and hole edge.
    // Exact over rationals of the injected doubles; no sampling, no epsilon.
    // Positive clearance from every edge means the segment never meets a ring,
    // so it stays in a's face of the ring arrangement; the clearance condition
    // is the point predicate's, minimised over the segment.
    bool contains_segment(double ax, double ay, double bx, double by) const;

private:
    explicit CutterCentreDomain2(CanonicalReachInput2 input);

    CanonicalReachInput2 input_;
    ReachFT squared_radius_;
};
