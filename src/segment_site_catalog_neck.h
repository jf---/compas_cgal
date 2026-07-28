#pragma once

#include "reachable_input_2.h"
#include "segment_site_neck_clearance.h"

#include <stdexcept>
#include <vector>

class IncompleteMatClearanceProfileGraphError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class UnsupportedCanonicalMatClearanceProfileError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class MatClearanceProfileGraph2 {
public:
  static MatClearanceProfileGraph2
  build(MatExactGraph2 graph, std::vector<MatClearanceEdgeProfile2> profiles);

  const MatExactGraph2 &graph() const noexcept;
  const std::vector<MatClearanceEdgeProfile2> &profiles() const noexcept;

private:
  MatClearanceProfileGraph2(MatExactGraph2 graph,
                            std::vector<MatClearanceEdgeProfile2> profiles);

  MatExactGraph2 graph_;
  std::vector<MatClearanceEdgeProfile2> profiles_;
};

MatClearanceProfileGraph2
canonical_l_shape_mat_clearance_graph(const CanonicalReachInput2 &input,
                                      const CORE::BigRat &radius_squared);
