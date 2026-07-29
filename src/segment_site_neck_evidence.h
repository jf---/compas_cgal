#pragma once

#include "segment_site_catalog_neck.h"
#include "segment_site_neck.h"

#include <optional>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

class InvalidMatNeckEvidenceGraphError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class InconsistentMatNeckNodeWidthError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class InvalidMatNeckEvidenceRecordError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class MatStrictEdgeNeckLocation2 {
public:
  MatStrictEdgeNeckLocation2(std::string edge_id,
                             std::string parameter_root_id);

  const std::string &edge_id() const noexcept;
  const std::string &parameter_root_id() const noexcept;

private:
  std::string edge_id_;
  std::string parameter_root_id_;
};

class MatClearanceEndpointNeckLocation2 {
public:
  MatClearanceEndpointNeckLocation2(std::string edge_id, std::string node_id,
                                    MatParameterEndpoint2 endpoint);

  const std::string &edge_id() const noexcept;
  const std::string &node_id() const noexcept;
  const MatParameterEndpoint2 &endpoint() const noexcept;

private:
  std::string edge_id_;
  std::string node_id_;
  MatParameterEndpoint2 endpoint_;
};

class MatSharedVertexNeckLocation2 {
public:
  MatSharedVertexNeckLocation2(std::string node_id,
                               std::vector<std::string> minimizing_edge_ids);

  const std::string &node_id() const noexcept;
  const std::vector<std::string> &minimizing_edge_ids() const noexcept;

private:
  std::string node_id_;
  std::vector<std::string> minimizing_edge_ids_;
};

class MatPlateauNeckLocation2 {
public:
  MatPlateauNeckLocation2(std::vector<std::string> node_ids,
                          std::vector<std::string> edge_ids);

  const std::vector<std::string> &node_ids() const noexcept;
  const std::vector<std::string> &edge_ids() const noexcept;

private:
  std::vector<std::string> node_ids_;
  std::vector<std::string> edge_ids_;
};

using MatExactNeckLocation2 =
    std::variant<MatStrictEdgeNeckLocation2, MatClearanceEndpointNeckLocation2,
                 MatSharedVertexNeckLocation2, MatPlateauNeckLocation2>;

class MatExactNeckEvidence2 {
public:
  static MatExactNeckEvidence2 build(MatExactNeckLocation2 location,
                                     std::string owner_id,
                                     std::vector<std::string> defining_site_ids,
                                     MatExactSquaredWidth2 squared_width,
                                     MatNeckSeparatingCut2 separating_cut);

  const MatExactNeckLocation2 &location() const noexcept;
  const std::string &owner_id() const noexcept;
  const std::vector<std::string> &defining_site_ids() const noexcept;
  const MatExactSquaredWidth2 &squared_width() const noexcept;
  const MatNeckSeparatingCut2 &separating_cut() const noexcept;

private:
  MatExactNeckEvidence2(MatExactNeckLocation2 location, std::string owner_id,
                        std::vector<std::string> defining_site_ids,
                        MatExactSquaredWidth2 squared_width,
                        MatNeckSeparatingCut2 separating_cut);

  MatExactNeckLocation2 location_;
  std::string owner_id_;
  std::vector<std::string> defining_site_ids_;
  MatExactSquaredWidth2 squared_width_;
  MatNeckSeparatingCut2 separating_cut_;
};

std::string mat_neck_location_tag(const MatExactNeckLocation2 &location);

std::vector<std::string>
mat_neck_location_edge_ids(const MatExactNeckLocation2 &location);

std::vector<std::string>
mat_neck_location_node_ids(const MatExactNeckLocation2 &location);

std::optional<std::string>
mat_neck_location_parameter_root_id(const MatExactNeckLocation2 &location);

std::vector<MatExactNeckEvidence2>
exact_neck_evidence(const MatClearanceProfileGraph2 &bundle);
