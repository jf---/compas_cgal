#include "segment_site_neck_evidence_bytes.h"

#include "canonical_encoding.h"
#include "continuous_tea_2/sha256.h"

#include <cstdint>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace {

using CanonicalFields = std::vector<std::pair<std::string, std::string>>;

std::string canonical_integer(const std::int64_t value) {
  return canonical_encode_integer(ExactAlgebraicInteger1(value));
}

std::string
canonical_boundary_feature(const MatEndpointBoundaryFeature2 &feature) {
  return canonical_encode_tagged_union(
      "mat-neck-endpoint-boundary-feature-v1",
      canonical_encode_component_map({
          {
              "domain-kind",
              canonical_integer(static_cast<std::int64_t>(feature.domain_kind)),
          },
          {
              "component",
              canonical_integer(feature.component),
          },
          {
              "curve-kind",
              canonical_integer(static_cast<std::int64_t>(feature.curve_kind)),
          },
          {
              "source-site-or-ring",
              canonical_integer(feature.source_site_or_ring),
          },
          {
              "derived-feature-index",
              canonical_integer(feature.derived_feature_index),
          },
      }));
}

std::string canonical_boundary_features(const MatParameterEndpoint2 &endpoint) {
  std::vector<std::string> records;
  records.reserve(endpoint.exact_evidence.boundary_features.size());
  for (const MatEndpointBoundaryFeature2 &feature :
       endpoint.exact_evidence.boundary_features) {
    records.push_back(canonical_boundary_feature(feature));
  }
  return canonical_encode_sequence(records);
}

std::string canonical_endpoint(const MatParameterEndpoint2 &endpoint) {
  if (!endpoint.parameter.has_value()) {
    throw InvalidMatNeckEvidenceBytesError(
        "neck endpoint evidence has no exact parameter");
  }
  return canonical_encode_tagged_union(
      "mat-neck-endpoint-v1",
      canonical_encode_component_map({
          {
              "parameter-root-id",
              mat_endpoint_root_identity_v1(endpoint),
          },
          {
              "provenance-ids",
              canonical_encode_sequence(endpoint.provenance_ids),
          },
          {
              "original-voronoi-vertex",
              canonical_encode_boolean(
                  endpoint.exact_evidence.original_voronoi_vertex),
          },
          {
              "domain-boundary",
              canonical_encode_boolean(endpoint.exact_evidence.domain_boundary),
          },
          {
              "clearance-root",
              canonical_encode_boolean(endpoint.exact_evidence.clearance_root),
          },
          {
              "boundary-features",
              canonical_boundary_features(endpoint),
          },
      }));
}

std::string canonical_separating_cut(const MatNeckSeparatingCut2 &cut) {
  std::vector<std::string> partitions;
  partitions.reserve(cut.edge_partitions().size());
  for (const std::vector<std::string> &partition : cut.edge_partitions()) {
    partitions.push_back(canonical_encode_sequence(partition));
  }
  return canonical_encode_tagged_union(
      "mat-neck-separating-cut-v1",
      canonical_encode_component_map({
          {
              "edge-partitions",
              canonical_encode_sequence(partitions),
          },
      }));
}

CanonicalFields common_fields(const MatExactNeckEvidence2 &evidence) {
  return {
      {
          "owner-id",
          evidence.owner_id(),
      },
      {
          "defining-site-ids",
          canonical_encode_sequence(evidence.defining_site_ids()),
      },
      {
          "squared-width-root-id",
          evidence.squared_width().root_id(),
      },
      {
          "separating-cut",
          canonical_separating_cut(evidence.separating_cut()),
      },
  };
}

std::string canonical_neck_evidence(const MatExactNeckEvidence2 &evidence) {
  return std::visit(
      [&evidence](const auto &location) {
        using Location = std::decay_t<decltype(location)>;
        CanonicalFields fields = common_fields(evidence);
        const std::string tag = mat_neck_location_tag(evidence.location());
        if constexpr (std::is_same_v<Location, MatStrictEdgeNeckLocation2>) {
          fields.emplace_back("edge-id", location.edge_id());
          fields.emplace_back("parameter-root-id",
                              location.parameter_root_id());
        } else if constexpr (std::is_same_v<
                                 Location, MatClearanceEndpointNeckLocation2>) {
          fields.emplace_back("edge-id", location.edge_id());
          fields.emplace_back("node-id", location.node_id());
          fields.emplace_back("endpoint",
                              canonical_endpoint(location.endpoint()));
        } else if constexpr (std::is_same_v<Location,
                                            MatSharedVertexNeckLocation2>) {
          fields.emplace_back("node-id", location.node_id());
          fields.emplace_back(
              "minimizing-edge-ids",
              canonical_encode_sequence(location.minimizing_edge_ids()));
        } else {
          static_assert(std::is_same_v<Location, MatPlateauNeckLocation2>);
          fields.emplace_back("node-ids",
                              canonical_encode_sequence(location.node_ids()));
          fields.emplace_back("edge-ids",
                              canonical_encode_sequence(location.edge_ids()));
        }
        return canonical_encode_tagged_union(
            tag, canonical_encode_component_map(fields));
      },
      evidence.location());
}

} // namespace

MatNeckEvidenceV1::MatNeckEvidenceV1(MatExactNeckEvidence2 evidence,
                                     std::string canonical_bytes,
                                     std::string canonical_digest)
    : evidence_(std::move(evidence)),
      canonical_bytes_(std::move(canonical_bytes)),
      canonical_digest_(std::move(canonical_digest)) {}

MatNeckEvidenceV1 MatNeckEvidenceV1::build(MatExactNeckEvidence2 evidence) {
  std::string bytes = canonical_neck_evidence(evidence);
  std::string digest = sha256_bytes(bytes);
  return MatNeckEvidenceV1(std::move(evidence), std::move(bytes),
                           std::move(digest));
}

const MatExactNeckEvidence2 &MatNeckEvidenceV1::evidence() const noexcept {
  return evidence_;
}

const std::string &MatNeckEvidenceV1::canonical_bytes() const noexcept {
  return canonical_bytes_;
}

const std::string &MatNeckEvidenceV1::canonical_digest() const noexcept {
  return canonical_digest_;
}

std::vector<MatNeckEvidenceV1>
exact_neck_evidence_v1(const MatClearanceProfileGraph2 &bundle) {
  std::vector<MatExactNeckEvidence2> evidence = exact_neck_evidence(bundle);
  std::vector<MatNeckEvidenceV1> records;
  records.reserve(evidence.size());
  for (MatExactNeckEvidence2 &item : evidence) {
    records.push_back(MatNeckEvidenceV1::build(std::move(item)));
  }
  return records;
}

std::vector<MatNeckEvidenceV1>
verify_neck_evidence_v1(const MatClearanceProfileGraph2 &bundle,
                        const std::vector<std::string> &records) {
  std::vector<MatNeckEvidenceV1> expected = exact_neck_evidence_v1(bundle);
  if (records.size() != expected.size()) {
    throw InvalidMatNeckEvidenceBytesError(
        "neck evidence exact replay expected " +
        std::to_string(expected.size()) + " records, received " +
        std::to_string(records.size()));
  }
  for (std::size_t index = 0; index < records.size(); ++index) {
    if (records[index] != expected[index].canonical_bytes()) {
      throw InvalidMatNeckEvidenceBytesError(
          "neck evidence bytes disagree with exact replay at record " +
          std::to_string(index));
    }
  }
  return expected;
}
