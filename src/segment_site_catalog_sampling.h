#pragma once

#include "segment_site_catalog_neck.h"
#include "segment_site_mat_sampling.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
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
  const std::optional<CORE::BigRat> &clearance_radius_squared() const noexcept;
  const std::vector<std::array<std::int64_t, 2>> &
  sample_exact_flags() const noexcept;

private:
  MatProposalSamplingGraph2(
      MatClearanceProfileGraph2 profile_graph,
      std::vector<std::int64_t> sample_offsets,
      std::vector<MatWorldXYProposalSample2> samples,
      std::optional<CORE::BigRat> clearance_radius_squared,
      std::vector<std::array<std::int64_t, 2>> sample_exact_flags);

  static MatProposalSamplingGraph2
  build_verified(MatClearanceProfileGraph2 profile_graph,
                 std::vector<MatProposalSamplingRun2> runs,
                 CORE::BigRat clearance_radius_squared);
  static MatProposalSamplingGraph2
  build_impl(MatClearanceProfileGraph2 profile_graph,
             std::vector<MatProposalSamplingRun2> runs,
             std::optional<CORE::BigRat> clearance_radius_squared,
             bool exact_verdicts);

  MatClearanceProfileGraph2 profile_graph_;
  std::vector<std::int64_t> sample_offsets_;
  std::vector<MatWorldXYProposalSample2> samples_;
  std::optional<CORE::BigRat> clearance_radius_squared_;
  std::vector<std::array<std::int64_t, 2>> sample_exact_flags_;

  friend MatProposalSamplingGraph2 canonical_l_shape_mat_proposal_graph(
      const CanonicalReachInput2 &input, const CORE::BigRat &radius_squared,
      const MatStationSpacingMm2 &station_spacing,
      const MatSagittaBoundMm2 &max_sagitta, std::size_t max_refinement_depth);
};

MatProposalSamplingGraph2 canonical_l_shape_mat_proposal_graph(
    const CanonicalReachInput2 &input, const CORE::BigRat &radius_squared,
    const MatStationSpacingMm2 &station_spacing,
    const MatSagittaBoundMm2 &max_sagitta, std::size_t max_refinement_depth);
