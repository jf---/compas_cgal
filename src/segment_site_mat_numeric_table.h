#pragma once

#include "segment_site_mat_proposal_table.h"
#include "segment_site_neck_evidence_bytes.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

class MissingMatNodeReportingCoordinateError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class UnknownMatNumericNeckCutEdgeError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class MatNumericNeckCutOverflowError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

struct MatNumericNeckCutTable2 {
  std::vector<std::string> neck_evidence;
  std::vector<std::int64_t> neck_cut_offsets;
  std::vector<std::int64_t> neck_cut_edge_ids;

  bool operator==(const MatNumericNeckCutTable2 &) const = default;
};

struct MatNumericMatTable2 {
  std::vector<std::array<double, 3>> nodes;
  MatNumericProposalTable2 proposal;
  std::vector<std::string> neck_evidence;
  std::vector<std::int64_t> neck_cut_offsets;
  std::vector<std::int64_t> neck_cut_edge_ids;

  bool operator==(const MatNumericMatTable2 &) const = default;
};

std::vector<std::array<double, 3>>
numeric_node_reporting_coordinates(const MatProposalSamplingGraph2 &sampled);

MatNumericNeckCutTable2
numeric_neck_cut_table(const MatExactGraph2 &graph,
                       const std::vector<MatNeckEvidenceV1> &evidence);

MatNumericMatTable2
canonical_l_shape_mat_numeric_table(const CanonicalReachInput2 &input,
                                    const MatStationSpacingMm2 &station_spacing,
                                    const MatSagittaBoundMm2 &max_sagitta,
                                    std::size_t max_refinement_depth);
