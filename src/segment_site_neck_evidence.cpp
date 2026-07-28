#include "segment_site_neck_evidence.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <optional>
#include <utility>
#include <variant>
#include <vector>

namespace {

struct MatEdgeEndpointBehaviors2 {
  MatClearanceEndpointBehavior2 lower;
  MatClearanceEndpointBehavior2 upper;
};

struct MatIncidentHalfedge2 {
  std::size_t edge;
  bool lower;
};

bool strictly_ordered(const std::vector<std::string> &identities) {
  return std::all_of(
             identities.begin(), identities.end(),
             [](const std::string &identity) { return !identity.empty(); }) &&
         std::is_sorted(identities.begin(), identities.end()) &&
         std::adjacent_find(identities.begin(), identities.end()) ==
             identities.end();
}

bool strictly_ordered_nonempty(const std::vector<std::string> &identities) {
  return !identities.empty() && strictly_ordered(identities);
}

bool canonical_clearance_endpoint(const MatParameterEndpoint2 &endpoint) {
  return endpoint.parameter.has_value() &&
         endpoint.exact_evidence.clearance_root &&
         endpoint.exact_evidence.domain_boundary ==
             !endpoint.exact_evidence.boundary_features.empty() &&
         strictly_ordered_nonempty(endpoint.provenance_ids) &&
         std::is_sorted(endpoint.exact_evidence.boundary_features.begin(),
                        endpoint.exact_evidence.boundary_features.end()) &&
         std::adjacent_find(endpoint.exact_evidence.boundary_features.begin(),
                            endpoint.exact_evidence.boundary_features.end()) ==
             endpoint.exact_evidence.boundary_features.end() &&
         std::binary_search(endpoint.provenance_ids.begin(),
                            endpoint.provenance_ids.end(),
                            algebraic_root_identity_v1(*endpoint.parameter));
}

std::size_t node_index(const MatExactGraph2 &graph,
                       const std::string &node_id) {
  const auto found = std::lower_bound(
      graph.nodes.begin(), graph.nodes.end(), node_id,
      [](const MatExactGraphNode2 &node, const std::string &identity) {
        return node.node_id < identity;
      });
  if (found == graph.nodes.end() || found->node_id != node_id) {
    throw InvalidMatNeckEvidenceGraphError(
        "MAT neck evidence edge references an unknown node");
  }
  return static_cast<std::size_t>(std::distance(graph.nodes.begin(), found));
}

const MatClearanceEndpointBehavior2 &
behavior(const std::vector<MatEdgeEndpointBehaviors2> &behaviors,
         const MatIncidentHalfedge2 &halfedge) {
  return halfedge.lower ? behaviors[halfedge.edge].lower
                        : behaviors[halfedge.edge].upper;
}

const MatParameterEndpoint2 &
endpoint(const std::vector<MatClearanceEdgeProfile2> &profiles,
         const MatIncidentHalfedge2 &halfedge) {
  return halfedge.lower ? profiles[halfedge.edge].lower()
                        : profiles[halfedge.edge].upper();
}

std::vector<std::string>
defining_sites(const std::vector<MatClearanceEdgeProfile2> &profiles,
               const std::vector<std::size_t> &edges) {
  std::vector<std::string> result;
  for (const std::size_t edge : edges) {
    result.insert(result.end(), profiles[edge].defining_site_ids().begin(),
                  profiles[edge].defining_site_ids().end());
  }
  std::sort(result.begin(), result.end());
  result.erase(std::unique(result.begin(), result.end()), result.end());
  return result;
}

std::vector<std::string> edge_ids(const MatExactGraph2 &graph,
                                  const std::vector<std::size_t> &edges) {
  std::vector<std::string> result;
  result.reserve(edges.size());
  for (const std::size_t edge : edges) {
    result.push_back(graph.edges[edge].edge_id);
  }
  return result;
}

std::optional<MatNeckSeparatingCut2>
strict_edge_cut(const MatNeckCutIndex2 &index, const std::string &edge_id) {
  try {
    return index.strict_edge_cut(edge_id);
  } catch (const NonSeparatingMatNeckCandidateError &) {
    return std::nullopt;
  }
}

std::optional<MatNeckSeparatingCut2> vertex_cut(const MatNeckCutIndex2 &index,
                                                const std::string &node_id) {
  try {
    return index.vertex_cut(node_id);
  } catch (const NonSeparatingMatNeckCandidateError &) {
    return std::nullopt;
  }
}

std::optional<MatNeckSeparatingCut2>
plateau_cut(const MatNeckCutIndex2 &index,
            const std::vector<std::string> &node_ids,
            const std::vector<std::string> &edge_ids) {
  try {
    return index.plateau_cut(node_ids, edge_ids);
  } catch (const NonSeparatingMatNeckCandidateError &) {
    return std::nullopt;
  }
}

std::string plateau_owner(const std::vector<std::string> &node_ids,
                          const std::vector<std::string> &edge_ids) {
  return std::min(node_ids.front(), edge_ids.front());
}

} // namespace

MatStrictEdgeNeckLocation2::MatStrictEdgeNeckLocation2(
    std::string edge_id, std::string parameter_root_id)
    : edge_id_(std::move(edge_id)),
      parameter_root_id_(std::move(parameter_root_id)) {}

const std::string &MatStrictEdgeNeckLocation2::edge_id() const noexcept {
  return edge_id_;
}

const std::string &
MatStrictEdgeNeckLocation2::parameter_root_id() const noexcept {
  return parameter_root_id_;
}

MatClearanceEndpointNeckLocation2::MatClearanceEndpointNeckLocation2(
    std::string edge_id, std::string node_id, MatParameterEndpoint2 endpoint)
    : edge_id_(std::move(edge_id)), node_id_(std::move(node_id)),
      endpoint_(std::move(endpoint)) {}

const std::string &MatClearanceEndpointNeckLocation2::edge_id() const noexcept {
  return edge_id_;
}

const std::string &MatClearanceEndpointNeckLocation2::node_id() const noexcept {
  return node_id_;
}

const MatParameterEndpoint2 &
MatClearanceEndpointNeckLocation2::endpoint() const noexcept {
  return endpoint_;
}

MatSharedVertexNeckLocation2::MatSharedVertexNeckLocation2(
    std::string node_id, std::vector<std::string> minimizing_edge_ids)
    : node_id_(std::move(node_id)),
      minimizing_edge_ids_(std::move(minimizing_edge_ids)) {}

const std::string &MatSharedVertexNeckLocation2::node_id() const noexcept {
  return node_id_;
}

const std::vector<std::string> &
MatSharedVertexNeckLocation2::minimizing_edge_ids() const noexcept {
  return minimizing_edge_ids_;
}

MatPlateauNeckLocation2::MatPlateauNeckLocation2(
    std::vector<std::string> node_ids, std::vector<std::string> edge_ids)
    : node_ids_(std::move(node_ids)), edge_ids_(std::move(edge_ids)) {}

const std::vector<std::string> &
MatPlateauNeckLocation2::node_ids() const noexcept {
  return node_ids_;
}

const std::vector<std::string> &
MatPlateauNeckLocation2::edge_ids() const noexcept {
  return edge_ids_;
}

MatExactNeckEvidence2::MatExactNeckEvidence2(
    MatExactNeckLocation2 location, std::string owner_id,
    std::vector<std::string> defining_site_ids,
    MatExactSquaredWidth2 squared_width, MatNeckSeparatingCut2 separating_cut)
    : location_(std::move(location)), owner_id_(std::move(owner_id)),
      defining_site_ids_(std::move(defining_site_ids)),
      squared_width_(std::move(squared_width)),
      separating_cut_(std::move(separating_cut)) {}

MatExactNeckEvidence2 MatExactNeckEvidence2::build(
    MatExactNeckLocation2 location, std::string owner_id,
    std::vector<std::string> defining_site_ids,
    MatExactSquaredWidth2 squared_width, MatNeckSeparatingCut2 separating_cut) {
  bool location_valid = false;
  bool requires_two_sites = false;
  if (const auto *strict = std::get_if<MatStrictEdgeNeckLocation2>(&location)) {
    requires_two_sites = true;
    location_valid =
        owner_id == strict->edge_id() && !strict->parameter_root_id().empty();
  } else if (const auto *clearance =
                 std::get_if<MatClearanceEndpointNeckLocation2>(&location)) {
    requires_two_sites = true;
    location_valid = owner_id == clearance->node_id() &&
                     !clearance->edge_id().empty() &&
                     canonical_clearance_endpoint(clearance->endpoint());
  } else if (const auto *shared =
                 std::get_if<MatSharedVertexNeckLocation2>(&location)) {
    location_valid = owner_id == shared->node_id() &&
                     shared->minimizing_edge_ids().size() >= 2 &&
                     strictly_ordered_nonempty(shared->minimizing_edge_ids());
  } else if (const auto *plateau =
                 std::get_if<MatPlateauNeckLocation2>(&location)) {
    location_valid =
        strictly_ordered_nonempty(plateau->node_ids()) &&
        strictly_ordered_nonempty(plateau->edge_ids()) &&
        owner_id == plateau_owner(plateau->node_ids(), plateau->edge_ids());
  }
  if (!location_valid || !strictly_ordered_nonempty(defining_site_ids) ||
      (requires_two_sites && defining_site_ids.size() != 2) ||
      squared_width.root_id().empty() ||
      separating_cut.edge_partitions().size() < 2) {
    throw InvalidMatNeckEvidenceRecordError(
        "exact MAT neck evidence record is not canonical");
  }
  return MatExactNeckEvidence2(
      std::move(location), std::move(owner_id), std::move(defining_site_ids),
      std::move(squared_width), std::move(separating_cut));
}

const MatExactNeckLocation2 &MatExactNeckEvidence2::location() const noexcept {
  return location_;
}

const std::string &MatExactNeckEvidence2::owner_id() const noexcept {
  return owner_id_;
}

const std::vector<std::string> &
MatExactNeckEvidence2::defining_site_ids() const noexcept {
  return defining_site_ids_;
}

const MatExactSquaredWidth2 &
MatExactNeckEvidence2::squared_width() const noexcept {
  return squared_width_;
}

const MatNeckSeparatingCut2 &
MatExactNeckEvidence2::separating_cut() const noexcept {
  return separating_cut_;
}

std::vector<MatExactNeckEvidence2>
exact_neck_evidence(const MatClearanceProfileGraph2 &bundle) {
  const MatExactGraph2 &graph = bundle.graph();
  const auto &profiles = bundle.profiles();
  const MatNeckCutIndex2 cut_index = MatNeckCutIndex2::build(graph);

  std::vector<MatEdgeEndpointBehaviors2> behaviors;
  behaviors.reserve(profiles.size());
  for (const auto &profile : profiles) {
    behaviors.push_back({
        lower_endpoint_clearance_behavior(profile),
        upper_endpoint_clearance_behavior(profile),
    });
  }

  std::vector<std::array<std::size_t, 2>> edge_nodes;
  edge_nodes.reserve(graph.edges.size());
  std::vector<std::vector<MatIncidentHalfedge2>> incident(graph.nodes.size());
  for (std::size_t edge = 0; edge < graph.edges.size(); ++edge) {
    const std::size_t source =
        node_index(graph, graph.edges[edge].source_node_id);
    const std::size_t target =
        node_index(graph, graph.edges[edge].target_node_id);
    edge_nodes.push_back({
        source,
        target,
    });
    incident[source].push_back({
        edge,
        true,
    });
    incident[target].push_back({
        edge,
        false,
    });
  }

  ExactAlgebraicKernel1 kernel;
  for (const auto &halfedges : incident) {
    if (halfedges.empty()) {
      throw InvalidMatNeckEvidenceGraphError(
          "MAT neck evidence graph contains an isolated node");
    }
    const auto &reference =
        behavior(behaviors, halfedges.front()).squared_width().value();
    for (const auto &halfedge : halfedges) {
      if (kernel.compare_1_object()(
              reference,
              behavior(behaviors, halfedge).squared_width().value()) !=
          CGAL::EQUAL) {
        throw InconsistentMatNeckNodeWidthError(
            "incident MAT clearance profiles disagree at a node");
      }
    }
  }

  std::vector<MatExactNeckEvidence2> strict_evidence;
  for (std::size_t edge = 0; edge < profiles.size(); ++edge) {
    for (const auto &minimum : strict_edge_clearance_minima(profiles[edge])) {
      auto cut = strict_edge_cut(cut_index, minimum.edge_id());
      if (!cut.has_value()) {
        continue;
      }
      strict_evidence.push_back(MatExactNeckEvidence2::build(
          MatStrictEdgeNeckLocation2(minimum.edge_id(),
                                     minimum.parameter_root_id()),
          minimum.edge_id(), minimum.defining_site_ids(),
          minimum.squared_width(), std::move(*cut)));
    }
  }

  std::vector<bool> constant_edges(graph.edges.size(), false);
  for (std::size_t edge = 0; edge < profiles.size(); ++edge) {
    constant_edges[edge] = profiles[edge].squared_clearance().is_constant();
  }
  std::vector<bool> visited_constant_edges(graph.edges.size(), false);
  std::vector<bool> plateau_owned_nodes(graph.nodes.size(), false);
  std::vector<MatExactNeckEvidence2> plateau_evidence;
  for (std::size_t seed = 0; seed < graph.edges.size(); ++seed) {
    if (!constant_edges[seed] || visited_constant_edges[seed]) {
      continue;
    }
    std::vector<std::size_t> pending{
        seed,
    };
    std::vector<std::size_t> component_edges;
    visited_constant_edges[seed] = true;
    while (!pending.empty()) {
      const std::size_t edge = pending.back();
      pending.pop_back();
      component_edges.push_back(edge);
      for (const std::size_t node : edge_nodes[edge]) {
        for (const auto &halfedge : incident[node]) {
          if (constant_edges[halfedge.edge] &&
              !visited_constant_edges[halfedge.edge]) {
            visited_constant_edges[halfedge.edge] = true;
            pending.push_back(halfedge.edge);
          }
        }
      }
    }
    std::sort(component_edges.begin(), component_edges.end());
    std::vector<bool> component_nodes(graph.nodes.size(), false);
    for (const std::size_t edge : component_edges) {
      component_nodes[edge_nodes[edge][0]] = true;
      component_nodes[edge_nodes[edge][1]] = true;
    }
    std::vector<std::string> component_node_ids;
    for (std::size_t node = 0; node < graph.nodes.size(); ++node) {
      if (!component_nodes[node]) {
        continue;
      }
      plateau_owned_nodes[node] = true;
      component_node_ids.push_back(graph.nodes[node].node_id);
    }

    std::vector<std::string> component_edge_ids =
        edge_ids(graph, component_edges);
    auto cut = plateau_cut(cut_index, component_node_ids, component_edge_ids);
    if (!cut.has_value()) {
      continue;
    }
    plateau_evidence.push_back(MatExactNeckEvidence2::build(
        MatPlateauNeckLocation2(component_node_ids, component_edge_ids),
        plateau_owner(component_node_ids, component_edge_ids),
        defining_sites(profiles, component_edges),
        behaviors[component_edges.front()].lower.squared_width(),
        std::move(*cut)));
  }

  std::vector<MatExactNeckEvidence2> endpoint_evidence;
  std::vector<MatExactNeckEvidence2> shared_evidence;
  for (std::size_t node = 0; node < graph.nodes.size(); ++node) {
    if (plateau_owned_nodes[node]) {
      continue;
    }
    std::vector<MatIncidentHalfedge2> minimizing;
    for (const auto &halfedge : incident[node]) {
      if (behavior(behaviors, halfedge).inward_clearance_sign() ==
          CGAL::POSITIVE) {
        minimizing.push_back(halfedge);
      }
    }
    if (minimizing.size() >= 2) {
      auto cut = vertex_cut(cut_index, graph.nodes[node].node_id);
      if (!cut.has_value()) {
        continue;
      }
      std::vector<std::size_t> minimizing_edges;
      minimizing_edges.reserve(minimizing.size());
      for (const auto &halfedge : minimizing) {
        minimizing_edges.push_back(halfedge.edge);
      }
      shared_evidence.push_back(MatExactNeckEvidence2::build(
          MatSharedVertexNeckLocation2(graph.nodes[node].node_id,
                                       edge_ids(graph, minimizing_edges)),
          graph.nodes[node].node_id, defining_sites(profiles, minimizing_edges),
          behavior(behaviors, minimizing.front()).squared_width(),
          std::move(*cut)));
      continue;
    }
    if (minimizing.size() != 1 ||
        !endpoint(profiles, minimizing.front()).exact_evidence.clearance_root) {
      continue;
    }
    auto cut = vertex_cut(cut_index, graph.nodes[node].node_id);
    if (!cut.has_value()) {
      continue;
    }
    const auto &halfedge = minimizing.front();
    endpoint_evidence.push_back(MatExactNeckEvidence2::build(
        MatClearanceEndpointNeckLocation2(graph.edges[halfedge.edge].edge_id,
                                          graph.nodes[node].node_id,
                                          endpoint(profiles, halfedge)),
        graph.nodes[node].node_id, profiles[halfedge.edge].defining_site_ids(),
        behavior(behaviors, halfedge).squared_width(), std::move(*cut)));
  }

  std::vector<MatExactNeckEvidence2> result;
  result.reserve(strict_evidence.size() + endpoint_evidence.size() +
                 shared_evidence.size() + plateau_evidence.size());
  result.insert(result.end(), std::make_move_iterator(strict_evidence.begin()),
                std::make_move_iterator(strict_evidence.end()));
  result.insert(result.end(),
                std::make_move_iterator(endpoint_evidence.begin()),
                std::make_move_iterator(endpoint_evidence.end()));
  result.insert(result.end(), std::make_move_iterator(shared_evidence.begin()),
                std::make_move_iterator(shared_evidence.end()));
  result.insert(result.end(), std::make_move_iterator(plateau_evidence.begin()),
                std::make_move_iterator(plateau_evidence.end()));
  return result;
}
