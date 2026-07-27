#include "segment_site_graph_emission.h"

void append_exact_graph_components(
    const std::string& dual_id,
    const std::string& primitive_kind,
    const std::vector<std::string>& generator_ids,
    const std::vector<MatAdmissibleComponent2>& components,
    MatExactGraph2& graph,
    std::map<std::string, std::size_t>& node_indices)
{
    const auto append_node =
        [&graph, &node_indices, &dual_id, &generator_ids](
            const MatParameterEndpoint2& bound_endpoint)
        {
            const std::string node_id =
                stable_endpoint_node_identity_v1(
                    dual_id,
                    bound_endpoint);
            const auto existing = node_indices.find(node_id);
            if (existing != node_indices.end())
            {
                union_stable_ids(
                    graph.nodes[existing->second]
                        .generator_site_ids,
                    generator_ids);
                union_stable_ids(
                    graph.nodes[existing->second]
                        .provenance_ids,
                    bound_endpoint.provenance_ids);
                return node_id;
            }
            node_indices.emplace(
                node_id,
                graph.nodes.size());
            graph.nodes.push_back(
                {
                    node_id,
                    bound_endpoint.provenance_ids,
                    generator_ids,
                });
            return node_id;
        };

    for (const MatAdmissibleComponent2& component : components)
    {
        const MatParameterEndpoint2 source_endpoint =
            exact_graph_endpoint_binding(component.lower);
        const MatParameterEndpoint2 target_endpoint =
            exact_graph_endpoint_binding(component.upper);
        const std::string source =
            append_node(source_endpoint);
        const std::string target =
            append_node(target_endpoint);
        graph.edges.push_back(
            {
                component.component_id,
                primitive_kind,
                source,
                target,
                source_endpoint,
                target_endpoint,
                generator_ids,
            });
    }
}
