#include "segment_site_mat_numeric_table.h"

#include "segment_site_catalog.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace {

std::size_t node_index(const MatExactGraph2 &graph,
                       const std::string &node_id) {
  const auto found = std::lower_bound(
      graph.nodes.begin(), graph.nodes.end(), node_id,
      [](const MatExactGraphNode2 &node, const std::string &identity) {
        return node.node_id < identity;
      });
  if (found == graph.nodes.end() || found->node_id != node_id) {
    throw MissingMatNodeReportingCoordinateError(
        "MAT reporting edge references an unknown node");
  }
  return static_cast<std::size_t>(std::distance(graph.nodes.begin(), found));
}

std::int64_t edge_index(const MatExactGraph2 &graph,
                        const std::string &edge_id) {
  const auto found = std::lower_bound(
      graph.edges.begin(), graph.edges.end(), edge_id,
      [](const MatExactGraphEdge2 &edge, const std::string &identity) {
        return edge.edge_id < identity;
      });
  if (found == graph.edges.end() || found->edge_id != edge_id) {
    throw UnknownMatNumericNeckCutEdgeError(
        "MAT neck cut references an unknown canonical edge");
  }
  const auto index =
      static_cast<std::size_t>(std::distance(graph.edges.begin(), found));
  if (index >
      static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max())) {
    throw MatNumericNeckCutOverflowError(
        "MAT neck cut edge index exceeds int64");
  }
  return static_cast<std::int64_t>(index);
}

} // namespace

std::vector<std::array<double, 3>>
numeric_node_reporting_coordinates(const MatProposalSamplingGraph2 &sampled) {
  const MatExactGraph2 &graph = sampled.profile_graph().graph();
  const std::array<std::int64_t, 2> verified{1, 1};
  if (sampled.sample_exact_flags().size() != sampled.samples().size() ||
      !std::all_of(
          sampled.sample_exact_flags().begin(),
          sampled.sample_exact_flags().end(),
          [&verified](const auto &flags) { return flags == verified; })) {
    throw UnverifiedMatProposalSamplingGraphError(
        "MAT node reporting coordinates require verified samples");
  }

  std::vector<std::optional<std::array<double, 3>>> assigned(
      graph.nodes.size());
  for (std::size_t edge = 0; edge < graph.edges.size(); ++edge) {
    const std::size_t begin =
        static_cast<std::size_t>(sampled.sample_offsets()[edge]);
    const std::size_t end =
        static_cast<std::size_t>(sampled.sample_offsets()[edge + 1]);
    const auto bind = [&assigned](const std::size_t node,
                                  const MatWorldXYProposalSample2 &sample) {
      if (!assigned[node].has_value()) {
        assigned[node] =
            std::array<double, 3>{sample.x_mm(), sample.y_mm(), 0.0};
      }
    };
    bind(node_index(graph, graph.edges[edge].source_node_id),
         sampled.samples()[begin]);
    bind(node_index(graph, graph.edges[edge].target_node_id),
         sampled.samples()[end - 1]);
  }

  std::vector<std::array<double, 3>> nodes;
  nodes.reserve(assigned.size());
  for (const auto &coordinate : assigned) {
    if (!coordinate.has_value()) {
      throw MissingMatNodeReportingCoordinateError(
          "MAT graph node has no incident reporting sample");
    }
    nodes.push_back(*coordinate);
  }
  return nodes;
}

MatNumericNeckCutTable2
numeric_neck_cut_table(const MatExactGraph2 &graph,
                       const std::vector<MatNeckEvidenceV1> &evidence) {
  MatNumericNeckCutTable2 table;
  table.neck_cut_offsets.push_back(0);
  table.neck_evidence.reserve(evidence.size());
  for (const MatNeckEvidenceV1 &record : evidence) {
    table.neck_evidence.push_back(record.canonical_bytes());
    std::vector<std::string> cut_edges;
    for (const auto &partition :
         record.evidence().separating_cut().edge_partitions()) {
      cut_edges.insert(cut_edges.end(), partition.begin(), partition.end());
    }
    std::sort(cut_edges.begin(), cut_edges.end());
    cut_edges.erase(std::unique(cut_edges.begin(), cut_edges.end()),
                    cut_edges.end());
    for (const std::string &edge_id : cut_edges) {
      table.neck_cut_edge_ids.push_back(edge_index(graph, edge_id));
    }
    if (table.neck_cut_edge_ids.size() >
        static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max())) {
      throw MatNumericNeckCutOverflowError("MAT neck cut offset exceeds int64");
    }
    table.neck_cut_offsets.push_back(
        static_cast<std::int64_t>(table.neck_cut_edge_ids.size()));
  }
  return table;
}

MatNumericMatTable2
canonical_l_shape_mat_numeric_table(const CanonicalReachInput2 &input,
                                    const MatStationSpacingMm2 &station_spacing,
                                    const MatSagittaBoundMm2 &max_sagitta,
                                    const std::size_t max_refinement_depth) {
  const MatToolRadiusMm2 tool_radius =
      MatToolRadiusMm2::build(input.binary64_radius);
  const MatProposalSamplingGraph2 sampled =
      canonical_l_shape_mat_proposal_graph(input, tool_radius.squared_exact(),
                                           station_spacing, max_sagitta,
                                           max_refinement_depth);
  MatNumericProposalTable2 proposal{
      numeric_graph_table(sampled.profile_graph().graph(),
                          canonical_mat_site_catalog(input)),
      numeric_sample_table(sampled, tool_radius),
  };
  const MatNumericNeckCutTable2 necks =
      numeric_neck_cut_table(sampled.profile_graph().graph(),
                             exact_neck_evidence_v1(sampled.profile_graph()));
  return {
      numeric_node_reporting_coordinates(sampled),
      std::move(proposal),
      necks.neck_evidence,
      necks.neck_cut_offsets,
      necks.neck_cut_edge_ids,
  };
}
