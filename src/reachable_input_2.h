#pragma once

#include "exact_region_2.h"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include <Eigen/Core>

namespace compas {

using RowMatrixXd =
    Eigen::Matrix<
        double,
        Eigen::Dynamic,
        Eigen::Dynamic,
        Eigen::RowMajor>;

} // namespace compas

struct CanonicalReachRing2 {
    bool outer = false;
    std::size_t canonical_ordinal = 0;
    std::vector<ReachKernelPoint> points;
    std::vector<std::array<double, 2>> binary64_points;
    std::string record;
};

struct ReachableArrangement2;

struct CanonicalReachInput2 {
    CanonicalReachRing2 outer;
    std::vector<CanonicalReachRing2> holes;
    ReachFT radius = ReachFT(0);
    double binary64_radius = 0.0;
    std::string recipe_record;

private:
    std::size_t input_vertex_count_ = 0;
    std::size_t ring_rotation_comparisons_ = 0;
    std::string integrity_record_;

    friend CanonicalReachInput2 canonical_reach_input(
        Eigen::Ref<const compas::RowMatrixXd> boundary,
        const std::vector<compas::RowMatrixXd>& holes,
        double tool_radius);
    friend ReachableArrangement2 build_reachable_arrangement(
        CanonicalReachInput2 input);
    friend void validate_canonical_reach_input(
        const CanonicalReachInput2& input);
};

CanonicalReachInput2 canonical_reach_input(
    Eigen::Ref<const compas::RowMatrixXd> boundary,
    const std::vector<compas::RowMatrixXd>& holes,
    double tool_radius);

void validate_canonical_reach_input(
    const CanonicalReachInput2& input);
