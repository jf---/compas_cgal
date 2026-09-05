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

private:
    explicit CutterCentreDomain2(CanonicalReachInput2 input);

    CanonicalReachInput2 input_;
    ReachFT squared_radius_;
};
