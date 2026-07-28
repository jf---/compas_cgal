#include "segment_site_mat_certificate.h"

#include "canonical_encoding.h"
#include "continuous_tea_2/sha256.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace {

std::string canonical_integer(const std::int64_t value) {
  return canonical_encode_integer(ExactAlgebraicInteger1(value));
}

std::string canonical_cardinality(const std::size_t value) {
  return canonical_encode_integer(ExactAlgebraicInteger1(value));
}

std::string
canonical_integer_sequence(const std::vector<std::int64_t> &values) {
  std::vector<std::string> records;
  records.reserve(values.size());
  for (const std::int64_t value : values) {
    records.push_back(canonical_integer(value));
  }
  return canonical_encode_sequence(records);
}

template <std::size_t Columns>
std::string canonical_integer_rows(
    const std::vector<std::array<std::int64_t, Columns>> &rows) {
  std::vector<std::string> records;
  records.reserve(rows.size());
  for (const auto &row : rows) {
    std::vector<std::string> fields;
    fields.reserve(Columns);
    for (const std::int64_t value : row) {
      fields.push_back(canonical_integer(value));
    }
    records.push_back(canonical_encode_sequence(fields));
  }
  return canonical_encode_sequence(records);
}

std::string
canonical_boundary_feature(const MatEndpointBoundaryFeature2 &feature) {
  return canonical_encode_tagged_union(
      "segment-site-mat-boundary-feature-v1",
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

std::string canonical_endpoint(const MatParameterEndpoint2 &endpoint) {
  if (!endpoint.parameter.has_value()) {
    throw InvalidMatCertificateEncodingError(
        "MAT certificate endpoint has no exact parameter");
  }
  std::vector<std::string> features;
  features.reserve(endpoint.exact_evidence.boundary_features.size());
  for (const MatEndpointBoundaryFeature2 &feature :
       endpoint.exact_evidence.boundary_features) {
    features.push_back(canonical_boundary_feature(feature));
  }
  return canonical_encode_tagged_union(
      "segment-site-mat-endpoint-v1",
      canonical_encode_component_map({
          {
              "parameter-root-id",
              algebraic_root_identity_v1(*endpoint.parameter),
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
              canonical_encode_sequence(features),
          },
      }));
}

std::string canonical_site(const CanonicalMatSite2 &site) {
  std::vector<std::string> endpoint_ids;
  if (site.endpoint_site_ids().has_value()) {
    endpoint_ids.assign(site.endpoint_site_ids()->begin(),
                        site.endpoint_site_ids()->end());
  }
  return canonical_encode_tagged_union(
      "segment-site-mat-site-v1",
      canonical_encode_component_map({
          {
              "stable-id",
              site.stable_id(),
          },
          {
              "kind",
              canonical_integer(
                  static_cast<std::int64_t>(site.provenance().kind)),
          },
          {
              "ring",
              canonical_integer(site.provenance().ring),
          },
          {
              "feature",
              canonical_integer(site.provenance().feature),
          },
          {
              "endpoint-site-ids",
              canonical_encode_sequence(endpoint_ids),
          },
      }));
}

std::string canonical_node(const MatExactGraphNode2 &node) {
  return canonical_encode_tagged_union(
      "segment-site-mat-node-v1",
      canonical_encode_component_map({
          {
              "node-id",
              node.node_id,
          },
          {
              "provenance-ids",
              canonical_encode_sequence(node.provenance_ids),
          },
          {
              "generator-site-ids",
              canonical_encode_sequence(node.generator_site_ids),
          },
          {
              "parent-site-ids",
              canonical_encode_sequence(node.parent_site_ids),
          },
      }));
}

std::string canonical_edge(const MatExactGraphEdge2 &edge) {
  return canonical_encode_tagged_union(
      "segment-site-mat-edge-v1",
      canonical_encode_component_map({
          {
              "edge-id",
              edge.edge_id,
          },
          {
              "primitive-kind",
              edge.primitive_kind,
          },
          {
              "original-dual-id",
              edge.original_dual_id,
          },
          {
              "source-node-id",
              edge.source_node_id,
          },
          {
              "target-node-id",
              edge.target_node_id,
          },
          {
              "source-endpoint",
              canonical_endpoint(edge.source_endpoint),
          },
          {
              "target-endpoint",
              canonical_endpoint(edge.target_endpoint),
          },
          {
              "generator-site-ids",
              canonical_encode_sequence(edge.generator_site_ids),
          },
          {
              "parent-site-ids",
              canonical_encode_sequence(edge.parent_site_ids),
          },
          {
              "clip-component-index",
              canonical_integer(edge.clip_component_index),
          },
          {
              "admissible-center-component",
              canonical_encode_boolean(edge.admissible_center_component),
          },
      }));
}

std::string canonical_profile(const MatClearanceEdgeProfile2 &profile) {
  std::vector<std::string> coefficients;
  coefficients.reserve(profile.squared_clearance().coefficients().size());
  for (const CORE::BigRat &coefficient :
       profile.squared_clearance().coefficients()) {
    coefficients.push_back(canonical_encode_rational(coefficient));
  }
  return canonical_encode_tagged_union(
      "segment-site-mat-clearance-profile-v1",
      canonical_encode_component_map({
          {
              "edge-id",
              profile.edge_id(),
          },
          {
              "defining-site-ids",
              canonical_encode_sequence(profile.defining_site_ids()),
          },
          {
              "lower-endpoint",
              canonical_endpoint(profile.lower()),
          },
          {
              "upper-endpoint",
              canonical_endpoint(profile.upper()),
          },
          {
              "squared-clearance-coefficients",
              canonical_encode_sequence(coefficients),
          },
      }));
}

std::string canonical_numeric_projection(const MatNumericGraphTable2 &graph) {
  return canonical_encode_tagged_union(
      "segment-site-mat-numeric-projection-v1",
      canonical_encode_component_map({
          {
              "node-identities",
              canonical_encode_sequence(graph.node_ids),
          },
          {
              "original-dual-identities",
              canonical_encode_sequence(graph.original_dual_ids),
          },
          {
              "edges",
              canonical_integer_rows(graph.edges),
          },
          {
              "node-site-offsets",
              canonical_integer_sequence(graph.node_site_offsets),
          },
          {
              "node-site-ids",
              canonical_integer_sequence(graph.node_site_ids),
          },
          {
              "site-provenance",
              canonical_integer_rows(graph.site_provenance),
          },
          {
              "edge-endpoint-provenance-flags",
              canonical_integer_rows(graph.edge_endpoint_provenance_flags),
          },
          {
              "endpoint-feature-offsets",
              canonical_integer_sequence(graph.endpoint_feature_offsets),
          },
          {
              "endpoint-features",
              canonical_integer_rows(graph.endpoint_features),
          },
          {
              "edge-exact-flags",
              canonical_integer_rows(graph.edge_exact_flags),
          },
      }));
}

std::string canonical_neck_projection(const MatNumericNeckCutTable2 &necks) {
  return canonical_encode_tagged_union(
      "segment-site-mat-neck-projection-v1",
      canonical_encode_component_map({
          {
              "neck-evidence",
              canonical_encode_sequence(necks.neck_evidence),
          },
          {
              "neck-cut-offsets",
              canonical_integer_sequence(necks.neck_cut_offsets),
          },
          {
              "neck-cut-edge-ids",
              canonical_integer_sequence(necks.neck_cut_edge_ids),
          },
      }));
}

std::string
canonical_certificate_bytes(const MatClearanceProfileGraph2 &exact,
                            const CanonicalMatSiteCatalog2 &catalog,
                            const MatNumericGraphTable2 &graph,
                            const MatNumericNeckCutTable2 &necks,
                            const std::string &center_domain_digest) {
  std::vector<std::string> sites;
  sites.reserve(catalog.sites().size());
  for (const CanonicalMatSite2 &site : catalog.sites()) {
    sites.push_back(canonical_site(site));
  }
  std::vector<std::string> nodes;
  nodes.reserve(exact.graph().nodes.size());
  for (const MatExactGraphNode2 &node : exact.graph().nodes) {
    nodes.push_back(canonical_node(node));
  }
  std::vector<std::string> edges;
  edges.reserve(exact.graph().edges.size());
  for (const MatExactGraphEdge2 &edge : exact.graph().edges) {
    edges.push_back(canonical_edge(edge));
  }
  std::vector<std::string> profiles;
  profiles.reserve(exact.profiles().size());
  for (const MatClearanceEdgeProfile2 &profile : exact.profiles()) {
    profiles.push_back(canonical_profile(profile));
  }

  return canonical_encode_tagged_union(
      "segment-site-mat-certificate-v1",
      canonical_encode_component_map({
          {
              "center-domain-digest",
              canonical_encode_bytes(center_domain_digest),
          },
          {
              "site-records",
              canonical_encode_sequence(sites),
          },
          {
              "node-records",
              canonical_encode_sequence(nodes),
          },
          {
              "edge-records",
              canonical_encode_sequence(edges),
          },
          {
              "clearance-profile-records",
              canonical_encode_sequence(profiles),
          },
          {
              "numeric-projection",
              canonical_numeric_projection(graph),
          },
          {
              "neck-projection",
              canonical_neck_projection(necks),
          },
          {
              "rejected-incident-transitions",
              canonical_cardinality(
                  exact.graph().rejected_incident_transitions),
          },
          {
              "matched-generator-sites",
              canonical_cardinality(exact.graph().matched_generator_sites),
          },
      }));
}

} // namespace

MatCertificateV1::MatCertificateV1(std::string canonical_bytes,
                                   std::string canonical_digest)
    : canonical_bytes_(std::move(canonical_bytes)),
      canonical_digest_(std::move(canonical_digest)) {}

const std::string &MatCertificateV1::canonical_bytes() const noexcept {
  return canonical_bytes_;
}

const std::string &MatCertificateV1::canonical_digest() const noexcept {
  return canonical_digest_;
}

std::string MatCertificateV1::release_canonical_bytes() && noexcept {
  return std::move(canonical_bytes_);
}

MatCertifiedExactProjection2
certified_mat_exact_projection_v1(const MatClearanceProfileGraph2 &exact,
                                  const CanonicalMatSiteCatalog2 &catalog,
                                  const std::string &center_domain_digest) {
  if (center_domain_digest.size() != 32) {
    throw InvalidMatCenterDomainDigestError(
        "MAT certificate center-domain digest must contain 32 bytes");
  }
  MatNumericGraphTable2 graph = numeric_graph_table(exact.graph(), catalog);
  std::vector<MatNeckEvidenceV1> evidence = exact_neck_evidence_v1(exact);
  MatNumericNeckCutTable2 necks =
      numeric_neck_cut_table(exact.graph(), evidence);
  std::string bytes = canonical_certificate_bytes(exact, catalog, graph, necks,
                                                  center_domain_digest);
  std::string digest = sha256_bytes(bytes);
  MatCertificateV1 certificate(std::move(bytes), std::move(digest));
  return {
      std::move(graph),
      std::move(necks),
      std::move(evidence),
      std::move(certificate),
  };
}

MatCertificateV1
replay_mat_certificate_v1(const MatClearanceProfileGraph2 &exact,
                          const CanonicalMatSiteCatalog2 &catalog,
                          const std::string &center_domain_digest,
                          const std::string &record) {
  MatCertifiedExactProjection2 expected =
      certified_mat_exact_projection_v1(exact, catalog, center_domain_digest);
  if (record != expected.certificate.canonical_bytes()) {
    throw InvalidMatCertificateReplayError(
        "MAT certificate bytes disagree with exact replay");
  }
  return std::move(expected.certificate);
}
