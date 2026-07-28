#include "segment_site_neck.h"
#include "segment_site_catalog_graph.h"

#include <algorithm>
#include <array>
#include <string>
#include <utility>
#include <vector>

namespace {

template <typename Graph>
constexpr bool mat_neck_index_buildable =
    requires(Graph&& graph) {
        MatNeckCutIndex2::build(
            std::forward<Graph>(graph));
    };

static_assert(
    mat_neck_index_buildable<MatExactGraph2&>);
static_assert(
    !mat_neck_index_buildable<MatExactGraph2>);

compas::RowMatrixXd matrix(
    const std::vector<std::array<double, 2>>&
        points)
{
    compas::RowMatrixXd result(
        static_cast<Eigen::Index>(
            points.size()),
        2);
    for (std::size_t index = 0;
         index < points.size();
         ++index) {
        result(
            static_cast<Eigen::Index>(index),
            0) = points[index][0];
        result(
            static_cast<Eigen::Index>(index),
            1) = points[index][1];
    }
    return result;
}

CanonicalReachInput2 l_shape_input()
{
    return canonical_reach_input(
        matrix(
            {
                {0.0, 0.0},
                {6.0, 0.0},
                {6.0, 2.0},
                {2.0, 2.0},
                {2.0, 6.0},
                {0.0, 6.0},
            }),
        {},
        1.0);
}

MatExactGraphNode2 node(std::string node_id)
{
    return {
        std::move(node_id),
        {},
        {},
        {},
    };
}

MatExactGraphEdge2 edge(
    std::string edge_id,
    std::string source,
    std::string target)
{
    const std::string dual_id =
        "dual/" + edge_id;
    return {
        std::move(edge_id),
        "LINE",
        dual_id,
        std::move(source),
        std::move(target),
        {},
        {},
        {},
        {},
        0,
        true,
    };
}

MatExactGraph2 chain_graph()
{
    return {
        {
            node("n0"),
            node("n1"),
            node("n2"),
        },
        {
            edge("e0", "n0", "n1"),
            edge("e1", "n1", "n2"),
        },
        0,
        0,
    };
}

MatExactGraph2 single_edge_graph()
{
    return {
        {
            node("n0"),
            node("n1"),
        },
        {
            edge("e0", "n0", "n1"),
        },
        0,
        0,
    };
}

MatExactGraph2 y_graph()
{
    return {
        {
            node("center"),
            node("leaf-a"),
            node("leaf-b"),
            node("leaf-c"),
        },
        {
            edge("edge-a", "center", "leaf-a"),
            edge("edge-b", "center", "leaf-b"),
            edge("edge-c", "center", "leaf-c"),
        },
        0,
        0,
    };
}

MatExactGraph2 cycle_graph()
{
    return {
        {
            node("n0"),
            node("n1"),
            node("n2"),
            node("n3"),
        },
        {
            edge("e0", "n0", "n1"),
            edge("e1", "n1", "n2"),
            edge("e2", "n2", "n3"),
            edge("e3", "n0", "n3"),
        },
        0,
        0,
    };
}

bool strict_edge_cut_is_canonical()
{
    const MatNeckSeparatingCut2 first =
        strict_edge_neck_separating_cut(
            chain_graph(),
            "e0");
    const MatNeckSeparatingCut2 repeated =
        strict_edge_neck_separating_cut(
            chain_graph(),
            "e0");
    const MatNeckSeparatingCut2 open_halves =
        strict_edge_neck_separating_cut(
            single_edge_graph(),
            "e0");
    return first == repeated
        && first.edge_partitions()
            == std::vector<
                std::vector<std::string>>{
                {"e0"},
                {"e0", "e1"},
            }
        && open_halves.edge_partitions()
            == std::vector<
                std::vector<std::string>>{
                {"e0"},
                {"e0"},
            };
}

bool vertex_cut_is_canonical()
{
    const MatNeckSeparatingCut2 chain =
        vertex_neck_separating_cut(
            chain_graph(),
            "n1");
    const MatNeckSeparatingCut2 branch =
        vertex_neck_separating_cut(
            y_graph(),
            "center");
    return chain.edge_partitions()
            == std::vector<
                std::vector<std::string>>{
                {"e0"},
                {"e1"},
            }
        && branch.edge_partitions()
            == std::vector<
                std::vector<std::string>>{
                {"edge-a"},
                {"edge-b"},
                {"edge-c"},
            };
}

bool plateau_cut_is_canonical()
{
    const MatExactGraph2 graph{
        {
            node("n0"),
            node("n1"),
            node("n2"),
            node("n3"),
        },
        {
            edge("e0", "n0", "n1"),
            edge("e1", "n1", "n2"),
            edge("e2", "n2", "n3"),
        },
        0,
        0,
    };
    return plateau_neck_separating_cut(
               graph,
               {"n1", "n2"},
               {"e1"})
            .edge_partitions()
        == std::vector<
            std::vector<std::string>>{
            {"e0"},
            {"e2"},
        };
}

bool l_shape_batch_index_matches_tree_topology()
{
    const MatExactGraph2 graph =
        canonical_l_shape_mat_graph(
            l_shape_input(),
            CORE::BigRat(0));
    if (graph.nodes.size() != 10
        || graph.edges.size() != 9) {
        return false;
    }
    const MatNeckCutIndex2 index =
        MatNeckCutIndex2::build(graph);
    std::vector<std::size_t> degrees(
        graph.nodes.size(),
        0);
    for (const MatExactGraphEdge2& edge :
         graph.edges) {
        if (index.strict_edge_cut(edge.edge_id)
                .edge_partitions()
                .size()
            != 2) {
            return false;
        }
        const auto source = std::lower_bound(
            graph.nodes.begin(),
            graph.nodes.end(),
            edge.source_node_id,
            [](const MatExactGraphNode2& node,
               const std::string& identity) {
                return node.node_id < identity;
            });
        const auto target = std::lower_bound(
            graph.nodes.begin(),
            graph.nodes.end(),
            edge.target_node_id,
            [](const MatExactGraphNode2& node,
               const std::string& identity) {
                return node.node_id < identity;
            });
        if (source == graph.nodes.end()
            || target == graph.nodes.end()) {
            return false;
        }
        ++degrees[static_cast<std::size_t>(
            std::distance(
                graph.nodes.begin(),
                source))];
        ++degrees[static_cast<std::size_t>(
            std::distance(
                graph.nodes.begin(),
                target))];
    }

    std::size_t leaves = 0;
    std::size_t degree_two = 0;
    std::size_t degree_three = 0;
    for (std::size_t node_index = 0;
         node_index < graph.nodes.size();
         ++node_index) {
        const std::size_t degree =
            degrees[node_index];
        if (degree == 1) {
            ++leaves;
            try {
                static_cast<void>(
                    index.vertex_cut(
                        graph.nodes[node_index]
                            .node_id));
            } catch (
                const NonSeparatingMatNeckCandidateError&) {
                continue;
            }
            return false;
        }
        if (index.vertex_cut(
                     graph.nodes[node_index]
                         .node_id)
                .edge_partitions()
                .size()
            != degree) {
            return false;
        }
        degree_two += degree == 2;
        degree_three += degree == 3;
    }
    return leaves == 5
        && degree_two == 2
        && degree_three == 3;
}

bool nonseparating_candidates_are_rejected()
{
    bool leaf_rejected = false;
    try {
        static_cast<void>(
            vertex_neck_separating_cut(
                chain_graph(),
                "n0"));
    } catch (
        const NonSeparatingMatNeckCandidateError&) {
        leaf_rejected = true;
    }

    bool cycle_edge_rejected = false;
    try {
        static_cast<void>(
            strict_edge_neck_separating_cut(
                cycle_graph(),
                "e0"));
    } catch (
        const NonSeparatingMatNeckCandidateError&) {
        cycle_edge_rejected = true;
    }

    bool cycle_vertex_rejected = false;
    try {
        static_cast<void>(
            vertex_neck_separating_cut(
                cycle_graph(),
                "n0"));
    } catch (
        const NonSeparatingMatNeckCandidateError&) {
        cycle_vertex_rejected = true;
    }
    return leaf_rejected
        && cycle_edge_rejected
        && cycle_vertex_rejected;
}

bool unrelated_component_does_not_create_a_cut()
{
    MatExactGraph2 disconnected{
        {
            node("a0"),
            node("a1"),
            node("b0"),
            node("b1"),
        },
        {
            edge("e0", "a0", "a1"),
            edge("e1", "b0", "b1"),
        },
        0,
        0,
    };
    try {
        static_cast<void>(
            vertex_neck_separating_cut(
                disconnected,
                "a0"));
    } catch (
        const NonSeparatingMatNeckCandidateError&) {
        return true;
    }
    return false;
}

bool malformed_cut_inputs_are_rejected()
{
    bool graph_order_rejected = false;
    try {
        MatExactGraph2 malformed =
            chain_graph();
        std::swap(
            malformed.nodes[0],
            malformed.nodes[1]);
        static_cast<void>(
            vertex_neck_separating_cut(
                malformed,
                "n1"));
    } catch (
        const InvalidMatNeckCutGraphError&) {
        graph_order_rejected = true;
    }

    bool unknown_target_rejected = false;
    try {
        static_cast<void>(
            strict_edge_neck_separating_cut(
                chain_graph(),
                "unknown"));
    } catch (
        const UnknownMatNeckCutTargetError&) {
        unknown_target_rejected = true;
    }

    bool malformed_plateau_rejected = false;
    try {
        static_cast<void>(
            plateau_neck_separating_cut(
                chain_graph(),
                {"n2", "n1"},
                {"e1"}));
    } catch (
        const InvalidMatNeckPlateauError&) {
        malformed_plateau_rejected = true;
    }

    bool disconnected_plateau_rejected = false;
    try {
        const MatExactGraph2 disconnected{
            {
                node("a0"),
                node("a1"),
                node("b0"),
                node("b1"),
            },
            {
                edge("e0", "a0", "a1"),
                edge("e1", "b0", "b1"),
            },
            0,
            0,
        };
        static_cast<void>(
            plateau_neck_separating_cut(
                disconnected,
                {"a0", "b0"},
                {}));
    } catch (
        const InvalidMatNeckPlateauError&) {
        disconnected_plateau_rejected =
            true;
    }
    bool omitted_induced_edge_rejected =
        false;
    try {
        static_cast<void>(
            plateau_neck_separating_cut(
                chain_graph(),
                {"n0", "n1", "n2"},
                {"e0"}));
    } catch (
        const InvalidMatNeckPlateauError&) {
        omitted_induced_edge_rejected =
            true;
    }
    bool escaping_plateau_edge_rejected =
        false;
    try {
        static_cast<void>(
            plateau_neck_separating_cut(
                chain_graph(),
                {"n0", "n1"},
                {"e1"}));
    } catch (
        const InvalidMatNeckPlateauError&) {
        escaping_plateau_edge_rejected =
            true;
    }
    return graph_order_rejected
        && unknown_target_rejected
        && malformed_plateau_rejected
        && disconnected_plateau_rejected
        && omitted_induced_edge_rejected
        && escaping_plateau_edge_rejected;
}

} // namespace

bool neck_cut_gate()
{
    return strict_edge_cut_is_canonical()
        && vertex_cut_is_canonical()
        && plateau_cut_is_canonical()
        && l_shape_batch_index_matches_tree_topology()
        && nonseparating_candidates_are_rejected()
        && unrelated_component_does_not_create_a_cut()
        && malformed_cut_inputs_are_rejected();
}
