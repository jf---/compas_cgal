#pragma once

#include "segment_site_catalog_neck.h"
#include "segment_site_mat_sampling.h"

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

class IncompleteMatProposalSamplingGraphError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class MatSamplingOffsetOverflowError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class UnsupportedCanonicalMatSamplingCurveError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class MatProposalSamplingGraph2 {
public:
  static MatProposalSamplingGraph2
  build(MatClearanceProfileGraph2 profile_graph,
        std::vector<MatProposalSamplingRun2> runs);

  const MatClearanceProfileGraph2 &profile_graph() const noexcept;
  const std::vector<std::int64_t> &sample_offsets() const noexcept;
  const std::vector<MatWorldXYProposalSample2> &samples() const noexcept;

private:
  MatProposalSamplingGraph2(MatClearanceProfileGraph2 profile_graph,
                            std::vector<std::int64_t> sample_offsets,
                            std::vector<MatWorldXYProposalSample2> samples);

  MatClearanceProfileGraph2 profile_graph_;
  std::vector<std::int64_t> sample_offsets_;
  std::vector<MatWorldXYProposalSample2> samples_;
};

MatProposalSamplingGraph2 canonical_l_shape_mat_proposal_graph(
    const CanonicalReachInput2 &input, const CORE::BigRat &radius_squared,
    const MatStationSpacingMm2 &station_spacing,
    const MatSagittaBoundMm2 &max_sagitta, std::size_t max_refinement_depth);
