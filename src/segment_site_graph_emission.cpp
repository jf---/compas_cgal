#include "segment_site_graph_emission.h"

std::vector<MatAdmissibleComponent2>
one_dimensional_graph_components(
    const std::vector<MatAdmissibleComponent2>& components)
{
    ExactAlgebraicKernel1 kernel;
    const auto compare =
        kernel.compare_1_object();
    std::vector<MatAdmissibleComponent2> retained;
    retained.reserve(components.size());
    for (const MatAdmissibleComponent2& component :
         components) {
        if (!component.lower.parameter.has_value()
            || !component.upper.parameter.has_value()) {
            throw InvalidSegmentSiteGraphComponentError(
                "segment-site graph component is unbounded");
        }
        const CGAL::Comparison_result order =
            compare(
                *component.lower.parameter,
                *component.upper.parameter);
        if (order == CGAL::LARGER) {
            throw InvalidSegmentSiteGraphComponentError(
                "segment-site graph component bounds are reversed");
        }
        if (order == CGAL::SMALLER) {
            retained.push_back(component);
        }
    }
    return retained;
}

void append_exact_graph_components(
    const std::string& dual_id,
    const std::string& primitive_kind,
    const std::vector<std::string>& generator_ids,
    const std::vector<MatAdmissibleComponent2>& components,
    MatExactGraph2& graph,
    std::map<std::string, std::size_t>& node_indices)
{
    append_exact_graph_components(
        dual_id,
        dual_id,
        primitive_kind,
        generator_ids,
        generator_ids,
        components,
        graph,
        node_indices);
}

void append_exact_graph_components(
    const std::string& piece_dual_id,
    const std::string& original_dual_id,
    const std::string& primitive_kind,
    const std::vector<std::string>& generator_ids,
    const std::vector<std::string>& parent_site_ids,
    const std::vector<MatAdmissibleComponent2>& components,
    MatExactGraph2& graph,
    std::map<std::string, std::size_t>& node_indices)
{
    const auto append_node =
        [&graph,
         &node_indices,
         &piece_dual_id,
         &generator_ids,
         &parent_site_ids](
            const MatParameterEndpoint2& bound_endpoint)
        {
            const std::string node_id =
                stable_endpoint_node_identity_v1(
                    piece_dual_id,
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
                        .parent_site_ids,
                    parent_site_ids);
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
                    parent_site_ids,
                });
            return node_id;
        };

    for (const MatAdmissibleComponent2& component :
         components)
    {
        if (component.component_index < 0) {
            throw InvalidSegmentSiteGraphComponentError(
                "MAT graph component has no exact clip index");
        }
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
                original_dual_id,
                source,
                target,
                source_endpoint,
                target_endpoint,
                generator_ids,
                parent_site_ids,
                component.component_index,
                true,
            });
    }
}
