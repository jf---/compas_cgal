#include "canonical_encoding.h"
#include "continuous_tea_2/sha256.h"
#include "segment_site_neck_evidence.h"
#include "segment_site_neck_evidence_bytes.h"

#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace {

MatExactGraphNode2 node(const std::string &node_id) {
  return {
      node_id,
      {"provenance-" + node_id},
      {"generator-" + node_id},
      {"parent-" + node_id},
  };
}

MatParameterEndpoint2 endpoint(const CORE::BigRat &parameter,
                               const bool clearance_root = false) {
  ExactAlgebraicKernel1 kernel;
  const auto algebraic = kernel.construct_algebraic_real_1_object()(parameter);
  return {
      algebraic,
      {
          algebraic_root_identity_v1(algebraic),
      },
      {
          !clearance_root,
          false,
          clearance_root,
          {},
      },
  };
}

MatExactGraphEdge2 edge(const std::string &edge_id,
                        const std::string &source_node_id,
                        const std::string &target_node_id,
                        const bool lower_clearance_root = false,
                        const bool upper_clearance_root = false) {
  return {
      edge_id,
      "LINE",
      "dual-" + edge_id,
      source_node_id,
      target_node_id,
      endpoint(0, lower_clearance_root),
      endpoint(1, upper_clearance_root),
      {
          edge_id + "-site-a",
          edge_id + "-site-b",
      },
      {
          edge_id + "-parent-a",
          edge_id + "-parent-b",
      },
      0,
      true,
  };
}

MatClearanceEdgeProfile2 profile(const MatExactGraphEdge2 &edge,
                                 RationalPolynomial coefficients) {
  return MatClearanceEdgeProfile2::build(
      edge, MatSquaredClearanceFunction2::build(std::move(coefficients)));
}

MatClearanceProfileGraph2 strict_bridge_bundle() {
  std::vector<MatExactGraphEdge2> edges{
      edge("e0", "n0", "n1"),
  };
  std::vector<MatClearanceEdgeProfile2> profiles{
      profile(edges[0],
              {
                  CORE::BigRat(5, 4),
                  -1,
                  1,
              }),
  };
  return MatClearanceProfileGraph2::build(
      {
          {
              node("n0"),
              node("n1"),
          },
          std::move(edges),
          0,
          0,
      },
      std::move(profiles));
}

MatClearanceProfileGraph2
shared_vertex_bundle(const bool single_clearance_minimum) {
  std::vector<MatExactGraphEdge2> edges{
      edge("e0", "n0", "n1", false, single_clearance_minimum),
      edge("e1", "n1", "n2"),
      edge("e2", "n1", "n3"),
  };
  std::vector<MatClearanceEdgeProfile2> profiles{
        profile(edges[0], {2, -2, 1}),
        profile(
            edges[1],
            single_clearance_minimum
                ? RationalPolynomial{
                      1,
                      CORE::BigRat(-1, 2),
                  }
                : RationalPolynomial{
                      1,
                      0,
                      1,
                  }),
        profile(
            edges[2],
            {
                1,
                CORE::BigRat(-1, 2),
            }),
    };
  return MatClearanceProfileGraph2::build(
      {
          {
              node("n0"),
              node("n1"),
              node("n2"),
              node("n3"),
          },
          std::move(edges),
          0,
          0,
      },
      std::move(profiles));
}

MatClearanceProfileGraph2 plateau_bundle(const bool descending_exit) {
  std::vector<MatExactGraphEdge2> edges{
      edge("e0", "n0", "n1"),
      edge("e1", "n1", "n2"),
      edge("e2", "n2", "n3"),
      edge("e3", "n3", "n4"),
  };
  std::vector<MatClearanceEdgeProfile2> profiles{
        profile(edges[0], {2, -2, 1}),
        profile(edges[1], {1}),
        profile(edges[2], {1}),
        profile(
            edges[3],
            descending_exit
                ? RationalPolynomial{
                      1,
                      CORE::BigRat(-1, 2),
                  }
                : RationalPolynomial{
                      1,
                      0,
                      1,
                  }),
    };
  return MatClearanceProfileGraph2::build(
      {
          {
              node("n0"),
              node("n1"),
              node("n2"),
              node("n3"),
              node("n4"),
          },
          std::move(edges),
          0,
          0,
      },
      std::move(profiles));
}

MatClearanceProfileGraph2 inconsistent_node_width_bundle() {
  std::vector<MatExactGraphEdge2> edges{
      edge("e0", "n0", "n1"),
      edge("e1", "n1", "n2"),
  };
  std::vector<MatClearanceEdgeProfile2> profiles{
      profile(edges[0], {2, -2, 1}),
      profile(edges[1], {2, 0, 1}),
  };
  return MatClearanceProfileGraph2::build(
      {
          {
              node("n0"),
              node("n1"),
              node("n2"),
          },
          std::move(edges),
          0,
          0,
      },
      std::move(profiles));
}

MatClearanceProfileGraph2 strict_cycle_bundle() {
  std::vector<MatExactGraphEdge2> edges{
      edge("e0", "n0", "n1"),
      edge("e1", "n0", "n2"),
      edge("e2", "n1", "n2"),
  };
  std::vector<MatClearanceEdgeProfile2> profiles{
      profile(edges[0], {CORE::BigRat(5, 4), -1, 1}),
      profile(edges[1], {CORE::BigRat(5, 4), -1, 1}),
      profile(edges[2], {CORE::BigRat(5, 4), -1, 1}),
  };
  return MatClearanceProfileGraph2::build(
      {
          {
              node("n0"),
              node("n1"),
              node("n2"),
          },
          std::move(edges),
          0,
          0,
      },
      std::move(profiles));
}

MatClearanceProfileGraph2 two_strict_bridges_bundle() {
  std::vector<MatExactGraphEdge2> edges{
      edge("e0", "n0", "n1"),
      edge("e1", "n1", "n2"),
  };
  std::vector<MatClearanceEdgeProfile2> profiles{
      profile(edges[0], {CORE::BigRat(5, 4), -1, 1}),
      profile(edges[1], {CORE::BigRat(5, 4), -1, 1}),
  };
  return MatClearanceProfileGraph2::build(
      {
          {
              node("n0"),
              node("n1"),
              node("n2"),
          },
          std::move(edges),
          0,
          0,
      },
      std::move(profiles));
}

bool width_is(const MatExactNeckEvidence2 &evidence,
              const CORE::BigRat &expected) {
  ExactAlgebraicKernel1 kernel;
  return kernel.compare_1_object()(evidence.squared_width().value(),
                                   expected) == CGAL::EQUAL;
}

bool strict_bridge_is_certified() {
  const auto evidence = exact_neck_evidence(strict_bridge_bundle());
  if (evidence.size() != 1 || evidence[0].owner_id() != "e0" ||
      evidence[0].defining_site_ids() !=
          std::vector<std::string>{
              "e0-site-a",
              "e0-site-b",
          } ||
      !width_is(evidence[0], 4) ||
      evidence[0].separating_cut().edge_partitions() !=
          std::vector<std::vector<std::string>>{
              {"e0"},
              {"e0"},
          }) {
    return false;
  }
  const auto *location =
      std::get_if<MatStrictEdgeNeckLocation2>(&evidence[0].location());
  return location != nullptr && location->edge_id() == "e0" &&
         location->parameter_root_id() == algebraic_root_id_v1({-1, 2}, 0);
}

bool incident_endpoint_minima_merge_once() {
  const auto evidence = exact_neck_evidence(shared_vertex_bundle(false));
  if (evidence.size() != 1 || evidence[0].owner_id() != "n1" ||
      !width_is(evidence[0], 4) ||
      evidence[0].separating_cut().edge_partitions() !=
          std::vector<std::vector<std::string>>{
              {"e0"},
              {"e1"},
              {"e2"},
          }) {
    return false;
  }
  const auto *location =
      std::get_if<MatSharedVertexNeckLocation2>(&evidence[0].location());
  return location != nullptr && location->node_id() == "n1" &&
         location->minimizing_edge_ids() == std::vector<std::string>{
                                                "e0",
                                                "e1",
                                            };
}

bool remaining_clearance_endpoint_is_owned() {
  const auto evidence = exact_neck_evidence(shared_vertex_bundle(true));
  if (evidence.size() != 1 || evidence[0].owner_id() != "n1" ||
      !width_is(evidence[0], 4)) {
    return false;
  }
  const auto *location =
      std::get_if<MatClearanceEndpointNeckLocation2>(&evidence[0].location());
  return location != nullptr && location->edge_id() == "e0" &&
         location->node_id() == "n1" &&
         location->endpoint().exact_evidence.clearance_root;
}

bool constant_cells_form_one_maximal_plateau() {
  const auto evidence = exact_neck_evidence(plateau_bundle(false));
  if (evidence.size() != 1 || evidence[0].owner_id() != "e1" ||
      evidence[0].defining_site_ids() !=
          std::vector<std::string>{
              "e1-site-a",
              "e1-site-b",
              "e2-site-a",
              "e2-site-b",
          } ||
      !width_is(evidence[0], 4) ||
      evidence[0].separating_cut().edge_partitions() !=
          std::vector<std::vector<std::string>>{
              {"e0"},
              {"e3"},
          }) {
    return false;
  }
  const auto *location =
      std::get_if<MatPlateauNeckLocation2>(&evidence[0].location());
  return location != nullptr &&
         location->edge_ids() ==
             std::vector<std::string>{
                 "e1",
                 "e2",
             } &&
         location->node_ids() == std::vector<std::string>{
                                     "n1",
                                     "n2",
                                     "n3",
                                 };
}

bool descending_exit_does_not_erase_plateau() {
  const auto evidence = exact_neck_evidence(plateau_bundle(true));
  return evidence.size() == 1 &&
         std::holds_alternative<MatPlateauNeckLocation2>(
             evidence[0].location()) &&
         evidence[0].owner_id() == "e1" && width_is(evidence[0], 4);
}

bool malformed_evidence_inputs_fail_loudly() {
  bool width_rejected = false;
  try {
    static_cast<void>(exact_neck_evidence(inconsistent_node_width_bundle()));
  } catch (const InconsistentMatNeckNodeWidthError &) {
    width_rejected = true;
  }

  bool isolated_rejected = false;
  try {
    const MatClearanceProfileGraph2 valid = strict_bridge_bundle();
    MatExactGraph2 graph = valid.graph();
    graph.nodes.push_back(node("n2"));
    static_cast<void>(exact_neck_evidence(
        MatClearanceProfileGraph2::build(std::move(graph), valid.profiles())));
  } catch (const InvalidMatNeckEvidenceGraphError &) {
    isolated_rejected = true;
  }

  bool record_rejected = false;
  try {
    const auto valid = exact_neck_evidence(strict_bridge_bundle());
    static_cast<void>(MatExactNeckEvidence2::build(
        valid[0].location(), "wrong-owner", valid[0].defining_site_ids(),
        valid[0].squared_width(), valid[0].separating_cut()));
  } catch (const InvalidMatNeckEvidenceRecordError &) {
    record_rejected = true;
  }

  bool endpoint_rejected = false;
  try {
    const auto valid = exact_neck_evidence(shared_vertex_bundle(true));
    const auto *location =
        std::get_if<MatClearanceEndpointNeckLocation2>(&valid[0].location());
    MatParameterEndpoint2 malformed = location->endpoint();
    malformed.provenance_ids.clear();
    static_cast<void>(MatExactNeckEvidence2::build(
        MatClearanceEndpointNeckLocation2(
            location->edge_id(), location->node_id(), std::move(malformed)),
        valid[0].owner_id(), valid[0].defining_site_ids(),
        valid[0].squared_width(), valid[0].separating_cut()));
  } catch (const InvalidMatNeckEvidenceRecordError &) {
    endpoint_rejected = true;
  }

  bool endpoint_projection_rejected = false;
  try {
    static_cast<void>(mat_neck_location_parameter_root_id(
        MatClearanceEndpointNeckLocation2(
            "e0",
            "n0",
            MatParameterEndpoint2{})));
  } catch (const InvalidMatNeckEvidenceRecordError &) {
    endpoint_projection_rejected = true;
  }

  bool site_count_rejected = false;
  try {
    const auto valid = exact_neck_evidence(strict_bridge_bundle());
    std::vector<std::string> malformed_sites = valid[0].defining_site_ids();
    malformed_sites.push_back("zz-site");
    static_cast<void>(MatExactNeckEvidence2::build(
        valid[0].location(), valid[0].owner_id(), std::move(malformed_sites),
        valid[0].squared_width(), valid[0].separating_cut()));
  } catch (const InvalidMatNeckEvidenceRecordError &) {
    site_count_rejected = true;
  }
  return width_rejected && isolated_rejected && record_rejected &&
         endpoint_rejected && endpoint_projection_rejected &&
         site_count_rejected;
}

bool nonseparating_minima_are_filtered() {
  return exact_neck_evidence(strict_cycle_bundle()).empty();
}

bool canonical_encoder_matches_algebraic_root_id() {
  return canonical_encode_tagged_union(
             "algebraic-root-id-v1", canonical_encode_component_map({
                                         {
                                             "coefficients",
                                             canonical_encode_sequence({
                                                 canonical_encode_integer(-2),
                                                 canonical_encode_integer(0),
                                                 canonical_encode_integer(1),
                                             }),
                                         },
                                         {
                                             "root-ordinal",
                                             canonical_encode_integer(0),
                                         },
                                     })) == algebraic_root_id_v1({-2, 0, 1}, 0);
}

bool malformed_canonical_encoding_fails_loudly() {
  bool empty_tag_rejected = false;
  try {
    static_cast<void>(canonical_encode_tagged_union("", "payload"));
  } catch (const InvalidCanonicalEncodingError &) {
    empty_tag_rejected = true;
  }

  bool empty_key_rejected = false;
  try {
    static_cast<void>(canonical_encode_component_map({{"", "value"}}));
  } catch (const InvalidCanonicalEncodingError &) {
    empty_key_rejected = true;
  }

  bool duplicate_key_rejected = false;
  try {
    static_cast<void>(canonical_encode_component_map(
        {{"field", "first"}, {"field", "second"}}));
  } catch (const InvalidCanonicalEncodingError &) {
    duplicate_key_rejected = true;
  }
  return empty_tag_rejected && empty_key_rejected && duplicate_key_rejected;
}

bool strict_record_matches_frozen_schema() {
  const std::vector<MatNeckEvidenceV1> records =
      exact_neck_evidence_v1(strict_bridge_bundle());
  if (records.size() != 1) {
    return false;
  }
  const std::string cut = canonical_encode_tagged_union(
      "mat-neck-separating-cut-v1",
      canonical_encode_component_map({
          {
              "edge-partitions",
              canonical_encode_sequence({
                  canonical_encode_sequence({"e0"}),
                  canonical_encode_sequence({"e0"}),
              }),
          },
      }));
  const std::string expected = canonical_encode_tagged_union(
      "mat-neck-strict-edge-v1", canonical_encode_component_map({
                                     {
                                         "owner-id",
                                         "e0",
                                     },
                                     {
                                         "defining-site-ids",
                                         canonical_encode_sequence({
                                             "e0-site-a",
                                             "e0-site-b",
                                         }),
                                     },
                                     {
                                         "squared-width-root-id",
                                         algebraic_root_id_v1({-4, 1}, 0),
                                     },
                                     {
                                         "separating-cut",
                                         cut,
                                     },
                                     {
                                         "edge-id",
                                         "e0",
                                     },
                                     {
                                         "parameter-root-id",
                                         algebraic_root_id_v1({-1, 2}, 0),
                                     },
                                 }));
  return records[0].canonical_bytes() == expected;
}

bool all_neck_variants_have_canonical_bytes() {
  const std::vector<MatNeckEvidenceV1> strict =
      exact_neck_evidence_v1(strict_bridge_bundle());
  const std::vector<MatNeckEvidenceV1> endpoint =
      exact_neck_evidence_v1(shared_vertex_bundle(true));
  const std::vector<MatNeckEvidenceV1> shared =
      exact_neck_evidence_v1(shared_vertex_bundle(false));
  const std::vector<MatNeckEvidenceV1> plateau =
      exact_neck_evidence_v1(plateau_bundle(false));
  if (strict.size() != 1 || endpoint.size() != 1 || shared.size() != 1 ||
      plateau.size() != 1) {
    return false;
  }
  const std::vector<const MatNeckEvidenceV1 *> records{
      &strict[0],
      &endpoint[0],
      &shared[0],
      &plateau[0],
  };
  const std::vector<std::string> tags{
      "mat-neck-strict-edge-v1",
      "mat-neck-clearance-endpoint-v1",
      "mat-neck-shared-vertex-v1",
      "mat-neck-plateau-v1",
  };
  for (std::size_t index = 0; index < records.size(); ++index) {
    if (records[index]->canonical_bytes().find(tags[index]) ==
            std::string::npos ||
        records[index]->canonical_digest() !=
            sha256_bytes(records[index]->canonical_bytes())) {
      return false;
    }
  }
  return strict[0].canonical_bytes() != endpoint[0].canonical_bytes() &&
         endpoint[0].canonical_bytes() != shared[0].canonical_bytes() &&
         shared[0].canonical_bytes() != plateau[0].canonical_bytes();
}

bool all_neck_variant_projections_are_complete() {
  const std::vector<MatExactNeckEvidence2> strict =
      exact_neck_evidence(strict_bridge_bundle());
  const std::vector<MatExactNeckEvidence2> endpoint =
      exact_neck_evidence(shared_vertex_bundle(true));
  const std::vector<MatExactNeckEvidence2> shared =
      exact_neck_evidence(shared_vertex_bundle(false));
  const std::vector<MatExactNeckEvidence2> plateau =
      exact_neck_evidence(plateau_bundle(false));
  if (strict.size() != 1 || endpoint.size() != 1 || shared.size() != 1 ||
      plateau.size() != 1) {
    return false;
  }
  const std::optional<std::string> endpoint_parameter =
      mat_neck_location_parameter_root_id(endpoint[0].location());
  return mat_neck_location_tag(strict[0].location()) ==
             "mat-neck-strict-edge-v1" &&
         mat_neck_location_edge_ids(strict[0].location()) ==
             std::vector<std::string>{"e0"} &&
         mat_neck_location_node_ids(strict[0].location()).empty() &&
         mat_neck_location_parameter_root_id(strict[0].location()) ==
             algebraic_root_id_v1({-1, 2}, 0) &&
         mat_neck_location_tag(endpoint[0].location()) ==
             "mat-neck-clearance-endpoint-v1" &&
         mat_neck_location_edge_ids(endpoint[0].location()) ==
             std::vector<std::string>{"e0"} &&
         mat_neck_location_node_ids(endpoint[0].location()) ==
             std::vector<std::string>{"n1"} &&
         endpoint_parameter.has_value() &&
         endpoint_parameter ==
             algebraic_root_identity_v1(
                 *std::get<MatClearanceEndpointNeckLocation2>(
                      endpoint[0].location())
                      .endpoint()
                      .parameter) &&
         mat_neck_location_tag(shared[0].location()) ==
             "mat-neck-shared-vertex-v1" &&
         mat_neck_location_edge_ids(shared[0].location()) ==
             std::vector<std::string>({"e0", "e1"}) &&
         mat_neck_location_node_ids(shared[0].location()) ==
             std::vector<std::string>{"n1"} &&
         !mat_neck_location_parameter_root_id(shared[0].location())
              .has_value() &&
         mat_neck_location_tag(plateau[0].location()) ==
             "mat-neck-plateau-v1" &&
         mat_neck_location_edge_ids(plateau[0].location()) ==
             std::vector<std::string>({"e1", "e2"}) &&
         mat_neck_location_node_ids(plateau[0].location()) ==
             std::vector<std::string>({"n1", "n2", "n3"}) &&
         !mat_neck_location_parameter_root_id(plateau[0].location())
              .has_value();
}

bool canonical_neck_bytes_replay_exactly() {
  const MatClearanceProfileGraph2 bundle = two_strict_bridges_bundle();
  const std::vector<MatNeckEvidenceV1> first = exact_neck_evidence_v1(bundle);
  const std::vector<MatNeckEvidenceV1> repeated =
      exact_neck_evidence_v1(bundle);
  if (first.size() != 2 || repeated.size() != first.size()) {
    return false;
  }
  std::vector<std::string> records;
  for (std::size_t index = 0; index < first.size(); ++index) {
    if (first[index].canonical_bytes() != repeated[index].canonical_bytes() ||
        first[index].canonical_digest() != repeated[index].canonical_digest()) {
      return false;
    }
    records.push_back(first[index].canonical_bytes());
  }
  const auto verified = verify_neck_evidence_v1(bundle, records);
  if (verified.size() != first.size()) {
    return false;
  }

  bool missing_rejected = false;
  try {
    std::vector<std::string> missing = records;
    missing.pop_back();
    static_cast<void>(verify_neck_evidence_v1(bundle, missing));
  } catch (const InvalidMatNeckEvidenceBytesError &) {
    missing_rejected = true;
  }

  bool mutation_rejected = false;
  try {
    std::vector<std::string> mutated = records;
    mutated[0].back() ^= 1;
    static_cast<void>(verify_neck_evidence_v1(bundle, mutated));
  } catch (const InvalidMatNeckEvidenceBytesError &) {
    mutation_rejected = true;
  }

  bool order_rejected = false;
  try {
    std::vector<std::string> reordered = records;
    std::swap(reordered[0], reordered[1]);
    static_cast<void>(verify_neck_evidence_v1(bundle, reordered));
  } catch (const InvalidMatNeckEvidenceBytesError &) {
    order_rejected = true;
  }
  return missing_rejected && mutation_rejected && order_rejected;
}

} // namespace

bool neck_evidence_gate() {
  return strict_bridge_is_certified() &&
         incident_endpoint_minima_merge_once() &&
         remaining_clearance_endpoint_is_owned() &&
         constant_cells_form_one_maximal_plateau() &&
         descending_exit_does_not_erase_plateau() &&
         malformed_evidence_inputs_fail_loudly() &&
         nonseparating_minima_are_filtered() &&
         canonical_encoder_matches_algebraic_root_id() &&
         malformed_canonical_encoding_fails_loudly() &&
         strict_record_matches_frozen_schema() &&
         all_neck_variants_have_canonical_bytes() &&
         all_neck_variant_projections_are_complete() &&
         canonical_neck_bytes_replay_exactly();
}
