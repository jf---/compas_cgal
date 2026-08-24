#pragma once

#include <Eigen/Core>

namespace compas {

using RowMatrixXd = Eigen::Matrix<
    double,
    Eigen::Dynamic,
    Eigen::Dynamic,
    Eigen::RowMajor>;
using RowMatrixXi = Eigen::Matrix<
    int,
    Eigen::Dynamic,
    Eigen::Dynamic,
    Eigen::RowMajor>;

} // namespace compas
