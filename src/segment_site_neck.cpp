#include "segment_site_neck.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <limits>
#include <utility>
#include <vector>

namespace {

std::size_t validated_node_index(
    const MatExactGraph2& graph,
    const std::string& node_id)
{
    const auto found = std::lower_bound(
        graph.nodes.begin(),
        graph.nodes.end(),
        node_id,
        [](const MatExactGraphNode2& node,
           const std::string& identity) {
            return node.node_id < identity;
        });
    if (found == graph.nodes.end()
        || found->node_id != node_id) {
        throw InvalidMatNeckCutGraphError(
            "MAT neck graph edge references an unknown node");
    }
    return static_cast<std::size_t>(
        std::distance(
            graph.nodes.begin(),
            found));
}

bool strictly_ordered_nonempty(
    const std::vector<std::string>& identities)
{
    return !identities.empty()
        && std::all_of(
            identities.begin(),
            identities.end(),
            [](const std::string& identity) {
                return !identity.empty();
            })
        && std::is_sorted(
            identities.begin(),
            identities.end())
        && std::adjacent_find(
               identities.begin(),
               identities.end())
            == identities.end();
}

} // namespace

MatNeckSeparatingCut2::MatNeckSeparatingCut2(
    std::vector<std::vector<std::string>>
        edge_partitions)
    : edge_partitions_(
          std::move(edge_partitions))
{
}

const std::vector<std::vector<std::string>>&
MatNeckSeparatingCut2::edge_partitions() const noexcept
{
    return edge_partitions_;
}

MatNeckCutIndex2::MatNeckCutIndex2(
    const MatExactGraph2& graph,
    std::vector<std::array<std::size_t, 2>>
        edge_nodes,
    std::vector<std::vector<std::size_t>>
        incident_edges)
    : graph_(&graph),
      edge_nodes_(std::move(edge_nodes)),
      incident_edges_(
          std::move(incident_edges))
{
}

MatNeckCutIndex2 MatNeckCutIndex2::build(
    const MatExactGraph2& graph)
{
    for (std::size_t index = 0;
         index < graph.nodes.size();
         ++index) {
        if (graph.nodes[index].node_id.empty()
            || (index > 0
                && graph.nodes[index].node_id
                    <= graph.nodes[index - 1]
                           .node_id)) {
            throw InvalidMatNeckCutGraphError(
                "MAT neck graph nodes are not strictly ordered");
        }
    }
    std::vector<std::array<std::size_t, 2>>
        edge_nodes;
    edge_nodes.reserve(graph.edges.size());
    std::vector<std::vector<std::size_t>>
        incident_edges(graph.nodes.size());
    for (std::size_t index = 0;
         index < graph.edges.size();
         ++index) {
        const MatExactGraphEdge2& edge =
            graph.edges[index];
        if (edge.edge_id.empty()
            || (index > 0
                && edge.edge_id
                    <= graph.edges[index - 1]
                           .edge_id)
            || edge.source_node_id
                == edge.target_node_id) {
            throw InvalidMatNeckCutGraphError(
                "MAT neck graph edges are not canonical");
        }
        const std::size_t source =
            validated_node_index(
                graph,
                edge.source_node_id);
        const std::size_t target =
            validated_node_index(
                graph,
                edge.target_node_id);
        edge_nodes.push_back({source, target});
        incident_edges[source].push_back(index);
        incident_edges[target].push_back(index);
    }
    return MatNeckCutIndex2(
        graph,
        std::move(edge_nodes),
        std::move(incident_edges));
}

std::size_t MatNeckCutIndex2::edge_index(
    const std::string& edge_id) const
{
    const auto found = std::lower_bound(
        graph_->edges.begin(),
        graph_->edges.end(),
        edge_id,
        [](const MatExactGraphEdge2& edge,
           const std::string& identity) {
            return edge.edge_id < identity;
        });
    if (found == graph_->edges.end()
        || found->edge_id != edge_id) {
        throw UnknownMatNeckCutTargetError(
            "MAT neck cut references an unknown edge");
    }
    return static_cast<std::size_t>(
        std::distance(
            graph_->edges.begin(),
            found));
}

std::size_t MatNeckCutIndex2::node_index(
    const std::string& node_id) const
{
    const auto found = std::lower_bound(
        graph_->nodes.begin(),
        graph_->nodes.end(),
        node_id,
        [](const MatExactGraphNode2& node,
           const std::string& identity) {
            return node.node_id < identity;
        });
    if (found == graph_->nodes.end()
        || found->node_id != node_id) {
        throw UnknownMatNeckCutTargetError(
            "MAT neck cut references an unknown node");
    }
    return static_cast<std::size_t>(
        std::distance(
            graph_->nodes.begin(),
            found));
}

std::size_t MatNeckCutIndex2::other_node(
    const std::size_t edge,
    const std::size_t node) const
{
    if (edge_nodes_[edge][0] == node) {
        return edge_nodes_[edge][1];
    }
    if (edge_nodes_[edge][1] == node) {
        return edge_nodes_[edge][0];
    }
    throw InvalidMatNeckCutGraphError(
        "MAT neck graph incidence index is inconsistent");
}

std::vector<bool>
MatNeckCutIndex2::connected_component(
    const std::size_t start) const
{
    std::vector<bool> included(
        graph_->nodes.size(),
        false);
    std::vector<std::size_t> pending{start};
    included[start] = true;
    while (!pending.empty()) {
        const std::size_t node =
            pending.back();
        pending.pop_back();
        for (const std::size_t edge :
             incident_edges_[node]) {
            const std::size_t neighbor =
                other_node(edge, node);
            if (!included[neighbor]) {
                included[neighbor] = true;
                pending.push_back(neighbor);
            }
        }
    }
    return included;
}

MatNeckSeparatingCut2
MatNeckCutIndex2::separating_cut(
    const std::vector<bool>& active_nodes,
    const std::vector<bool>& removed_nodes,
    const std::vector<bool>& removed_edges) const
{
    constexpr std::size_t no_component =
        std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> component_of(
        graph_->nodes.size(),
        no_component);
    std::size_t component_count = 0;
    for (std::size_t seed = 0;
         seed < graph_->nodes.size();
         ++seed) {
        if (!active_nodes[seed]
            || removed_nodes[seed]
            || component_of[seed]
                != no_component) {
            continue;
        }
        std::vector<std::size_t> pending{seed};
        component_of[seed] = component_count;
        while (!pending.empty()) {
            const std::size_t node =
                pending.back();
            pending.pop_back();
            for (const std::size_t edge :
                 incident_edges_[node]) {
                if (removed_edges[edge]) {
                    continue;
                }
                const std::size_t neighbor =
                    other_node(edge, node);
                if (!active_nodes[neighbor]
                    || removed_nodes[neighbor]
                    || component_of[neighbor]
                        != no_component) {
                    continue;
                }
                component_of[neighbor] =
                    component_count;
                pending.push_back(neighbor);
            }
        }
        ++component_count;
    }

    std::vector<std::vector<std::string>>
        partitions(component_count);
    for (std::size_t edge = 0;
         edge < graph_->edges.size();
         ++edge) {
        const auto endpoints =
            edge_nodes_[edge];
        const std::size_t source_component =
            component_of[endpoints[0]];
        const std::size_t target_component =
            component_of[endpoints[1]];
        const std::string& edge_id =
            graph_->edges[edge].edge_id;
        if (source_component != no_component
            && source_component
                == target_component
            && !removed_edges[edge]) {
            partitions[source_component]
                .push_back(edge_id);
            continue;
        }
        if (source_component != no_component
            && ((target_component
                     != no_component
                 && source_component
                     != target_component
                 && removed_edges[edge])
                || removed_nodes[
                    endpoints[1]])) {
            partitions[source_component]
                .push_back(edge_id);
        }
        if (target_component != no_component
            && ((source_component
                     != no_component
                 && source_component
                     != target_component
                 && removed_edges[edge])
                || removed_nodes[
                    endpoints[0]])) {
            partitions[target_component]
                .push_back(edge_id);
        }
    }
    partitions.erase(
        std::remove_if(
            partitions.begin(),
            partitions.end(),
            [](const auto& partition) {
                return partition.empty();
            }),
        partitions.end());
    std::sort(
        partitions.begin(),
        partitions.end());
    if (partitions.size() < 2) {
        throw NonSeparatingMatNeckCandidateError(
            "MAT neck candidate does not separate two traversal sides");
    }
    return MatNeckSeparatingCut2(
        std::move(partitions));
}

MatNeckSeparatingCut2
MatNeckCutIndex2::strict_edge_cut(
    const std::string& edge_id) const
{
    const std::size_t target =
        edge_index(edge_id);
    const std::vector<bool> active =
        connected_component(
            edge_nodes_[target][0]);
    std::vector<bool> removed_nodes(
        graph_->nodes.size(),
        false);
    std::vector<bool> removed_edges(
        graph_->edges.size(),
        false);
    removed_edges[target] = true;
    return separating_cut(
        active,
        removed_nodes,
        removed_edges);
}

MatNeckSeparatingCut2
MatNeckCutIndex2::vertex_cut(
    const std::string& node_id) const
{
    const std::size_t target =
        node_index(node_id);
    const std::vector<bool> active =
        connected_component(target);
    std::vector<bool> removed_nodes(
        graph_->nodes.size(),
        false);
    std::vector<bool> removed_edges(
        graph_->edges.size(),
        false);
    removed_nodes[target] = true;
    return separating_cut(
        active,
        removed_nodes,
        removed_edges);
}

MatNeckSeparatingCut2
MatNeckCutIndex2::plateau_cut(
    const std::vector<std::string>&
        plateau_node_ids,
    const std::vector<std::string>&
        plateau_edge_ids) const
{
    if (!strictly_ordered_nonempty(
            plateau_node_ids)
        || !strictly_ordered_nonempty(
            plateau_edge_ids)) {
        throw InvalidMatNeckPlateauError(
            "MAT neck plateau identities are not canonical");
    }
    std::vector<bool> plateau_nodes(
        graph_->nodes.size(),
        false);
    for (const std::string& node_id :
         plateau_node_ids) {
        plateau_nodes[node_index(node_id)] =
            true;
    }
    std::vector<bool> plateau_edges(
        graph_->edges.size(),
        false);
    for (const std::string& edge_id :
         plateau_edge_ids) {
        const std::size_t edge =
            edge_index(edge_id);
        const auto endpoints =
            edge_nodes_[edge];
        if (!plateau_nodes[endpoints[0]]
            || !plateau_nodes[endpoints[1]]) {
            throw InvalidMatNeckPlateauError(
                "MAT neck plateau edge leaves its node set");
        }
        plateau_edges[edge] = true;
    }
    for (std::size_t edge = 0;
         edge < graph_->edges.size();
         ++edge) {
        const auto endpoints =
            edge_nodes_[edge];
        if (plateau_nodes[endpoints[0]]
            && plateau_nodes[endpoints[1]]
            && !plateau_edges[edge]) {
            throw InvalidMatNeckPlateauError(
                "MAT neck plateau omits an induced edge");
        }
    }

    const std::size_t first =
        node_index(plateau_node_ids.front());
    std::vector<bool> reached(
        graph_->nodes.size(),
        false);
    std::vector<std::size_t> pending{first};
    reached[first] = true;
    while (!pending.empty()) {
        const std::size_t node =
            pending.back();
        pending.pop_back();
        for (const std::size_t edge :
             incident_edges_[node]) {
            if (!plateau_edges[edge]) {
                continue;
            }
            const std::size_t neighbor =
                other_node(edge, node);
            if (!reached[neighbor]) {
                reached[neighbor] = true;
                pending.push_back(neighbor);
            }
        }
    }
    for (std::size_t node = 0;
         node < graph_->nodes.size();
         ++node) {
        if (plateau_nodes[node]
            && !reached[node]) {
            throw InvalidMatNeckPlateauError(
                "MAT neck plateau subgraph is disconnected");
        }
    }
    const std::vector<bool> active =
        connected_component(first);
    return separating_cut(
        active,
        plateau_nodes,
        plateau_edges);
}

MatNeckSeparatingCut2
strict_edge_neck_separating_cut(
    const MatExactGraph2& graph,
    const std::string& edge_id)
{
    return MatNeckCutIndex2::build(graph)
        .strict_edge_cut(edge_id);
}

MatNeckSeparatingCut2
vertex_neck_separating_cut(
    const MatExactGraph2& graph,
    const std::string& node_id)
{
    return MatNeckCutIndex2::build(graph)
        .vertex_cut(node_id);
}

MatNeckSeparatingCut2
plateau_neck_separating_cut(
    const MatExactGraph2& graph,
    const std::vector<std::string>&
        plateau_node_ids,
    const std::vector<std::string>&
        plateau_edge_ids)
{
    return MatNeckCutIndex2::build(graph)
        .plateau_cut(
            plateau_node_ids,
            plateau_edge_ids);
}
