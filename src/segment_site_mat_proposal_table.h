#pragma once

#include "segment_site_catalog_sampling.h"
#include "segment_site_graph_csr.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

class InvalidMatProposalToolRadiusError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class UnverifiedMatProposalSamplingGraphError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class MismatchedMatProposalToolRadiusError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class InvalidMatProposalReportingDataError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class MatToolRadiusMm2 {
public:
  static MatToolRadiusMm2 build(double value);

  double value() const noexcept;
  const CORE::BigRat &squared_exact() const noexcept;

private:
  MatToolRadiusMm2(double value, CORE::BigRat squared_exact);

  double value_;
  CORE::BigRat squared_exact_;
};

struct MatNumericSampleTable2 {
  std::vector<std::array<double, 3>> sample_centers;
  std::vector<double> sample_clearance;
  std::vector<double> sample_guide_radius;
  std::vector<std::array<std::int64_t, 2>> sample_flags;
  std::vector<std::int64_t> edge_sample_offsets;
  std::vector<double> sample_parameter;

  bool operator==(const MatNumericSampleTable2 &) const = default;
};

struct MatNumericProposalTable2 {
  MatNumericGraphTable2 graph;
  MatNumericSampleTable2 samples;

  bool operator==(const MatNumericProposalTable2 &) const = default;
};

MatNumericSampleTable2
numeric_sample_table(const MatProposalSamplingGraph2 &sampled,
                     const MatToolRadiusMm2 &tool_radius);

MatNumericProposalTable2 canonical_l_shape_mat_numeric_proposal_table(
    const CanonicalReachInput2 &input,
    const MatStationSpacingMm2 &station_spacing,
    const MatSagittaBoundMm2 &max_sagitta, std::size_t max_refinement_depth);
