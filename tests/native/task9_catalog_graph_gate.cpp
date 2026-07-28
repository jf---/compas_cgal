#include "segment_site_catalog_graph.h"
#include "segment_site_catalog.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

using EdgeSignature = std::tuple<
    std::string,
    std::vector<std::string>,
    std::vector<std::string>,
    std::string,
    std::string,
    std::string,
    std::string,
    std::vector<std::string>,
    std::vector<std::string>>;

using NodeSignature = std::tuple<
    std::string,
    std::vector<std::string>,
    std::vector<std::string>,
    std::vector<std::string>,
    std::size_t>;

compas::RowMatrixXd rectangle(
    const bool transformed)
{
    const std::array<std::array<double, 2>, 4>
        canonical{{
            {-4.0, -2.0},
            {4.0, -2.0},
            {4.0, 2.0},
            {-4.0, 2.0},
        }};
    const std::array<std::array<double, 2>, 4>
        reversed{{
            {4.0, 2.0},
            {4.0, -2.0},
            {-4.0, -2.0},
            {-4.0, 2.0},
        }};
    const auto& points =
        transformed ? reversed : canonical;
    compas::RowMatrixXd result(
        static_cast<Eigen::Index>(points.size()),
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

CanonicalReachInput2 rectangle_input(
    const bool transformed)
{
    return canonical_reach_input(
        rectangle(transformed),
        {},
        1.0);
}

template <std::size_t PointCount>
CanonicalReachInput2 polygon_input(
    const std::array<
        std::array<double, 2>,
        PointCount>& points)
{
    compas::RowMatrixXd polygon(
        static_cast<Eigen::Index>(
            points.size()),
        2);
    for (std::size_t index = 0;
         index < points.size();
         ++index) {
        polygon(
            static_cast<Eigen::Index>(index),
            0) = points[index][0];
        polygon(
            static_cast<Eigen::Index>(index),
            1) = points[index][1];
    }
    return canonical_reach_input(
        polygon,
        {},
        1.0);
}

bool endpoints_equal(
    const MatParameterEndpoint2& left,
    const MatParameterEndpoint2& right)
{
    if (!left.parameter.has_value()
        || !right.parameter.has_value()
        || left.provenance_ids
            != right.provenance_ids) {
        return false;
    }
    ExactAlgebraicKernel1 kernel;
    return kernel.compare_1_object()(
               *left.parameter,
               *right.parameter)
        == CGAL::EQUAL;
}

bool graphs_equal(
    const MatExactGraph2& left,
    const MatExactGraph2& right)
{
    if (left.nodes.size() != right.nodes.size()
        || left.edges.size() != right.edges.size()
        || left.rejected_incident_transitions
            != right.rejected_incident_transitions
        || left.matched_generator_sites
            != right.matched_generator_sites) {
        return false;
    }
    for (std::size_t index = 0;
         index < left.nodes.size();
         ++index) {
        const MatExactGraphNode2& first =
            left.nodes[index];
        const MatExactGraphNode2& second =
            right.nodes[index];
        if (first.node_id != second.node_id
            || first.provenance_ids
                != second.provenance_ids
            || first.generator_site_ids
                != second.generator_site_ids
            || first.parent_site_ids
                != second.parent_site_ids) {
            return false;
        }
    }
    for (std::size_t index = 0;
         index < left.edges.size();
         ++index) {
        const MatExactGraphEdge2& first =
            left.edges[index];
        const MatExactGraphEdge2& second =
            right.edges[index];
        if (first.edge_id != second.edge_id
            || first.primitive_kind
                != second.primitive_kind
            || first.original_dual_id
                != second.original_dual_id
            || first.source_node_id
                != second.source_node_id
            || first.target_node_id
                != second.target_node_id
            || first.generator_site_ids
                != second.generator_site_ids
            || first.parent_site_ids
                != second.parent_site_ids
            || !endpoints_equal(
                first.source_endpoint,
                second.source_endpoint)
            || !endpoints_equal(
                first.target_endpoint,
                second.target_endpoint)) {
            return false;
        }
    }
    return true;
}

std::map<std::string, std::string>
friendly_site_names(
    const CanonicalMatSiteCatalog2& catalog)
{
    const std::array<std::string, 4> point_names{
        "lower-left",
        "lower-right",
        "upper-right",
        "upper-left",
    };
    const std::array<std::string, 4> segment_names{
        "bottom-segment",
        "right-segment",
        "top-segment",
        "left-segment",
    };
    std::map<std::string, std::string> names;
    for (const CanonicalMatSite2& site :
         catalog.sites()) {
        const MatSiteProvenance2& provenance =
            site.provenance();
        if (provenance.ring != 0
            || provenance.feature < 0
            || provenance.feature
                >= static_cast<std::int64_t>(
                    point_names.size())) {
            return {};
        }
        const std::size_t feature =
            static_cast<std::size_t>(
                provenance.feature);
        names.emplace(
            site.stable_id(),
            provenance.kind == MatSiteKind2::Point
                ? point_names[feature]
                : segment_names[feature]);
    }
    return names;
}

std::vector<std::string> normalized_ids(
    const std::vector<std::string>& ids,
    const std::map<std::string, std::string>&
        names)
{
    std::vector<std::string> normalized;
    normalized.reserve(ids.size());
    for (const std::string& id : ids) {
        const auto found = names.find(id);
        if (found != names.end()) {
            normalized.push_back(found->second);
            continue;
        }
        const auto prefix = std::find_if(
            names.begin(),
            names.end(),
            [&id](const auto& candidate) {
                return id.size()
                        > candidate.first.size()
                    && id.compare(
                           0,
                           candidate.first.size(),
                           candidate.first)
                        == 0
                    && id[candidate.first.size()]
                        == '/';
            });
        normalized.push_back(
            prefix == names.end()
                ? id
                : prefix->second
                    + id.substr(
                        prefix->first.size()));
    }
    std::sort(
        normalized.begin(),
        normalized.end());
    return normalized;
}

std::string projected_dual_identity(
    const MatExactGraphEdge2& edge,
    const std::map<std::string, std::string>&
        site_names)
{
    if (edge.generator_site_ids.size() != 2) {
        return {};
    }
    const std::array<std::string, 3> kinds{
        "segment-segment",
        "segment-segment/branch-negative",
        "segment-segment/branch-positive",
    };
    const std::vector<std::string>
        friendly_generators =
            normalized_ids(
                edge.generator_site_ids,
                site_names);
    for (const std::string& kind : kinds) {
        if (edge.original_dual_id
            == stable_dual_identity_v1(
                kind,
                edge.generator_site_ids)) {
            return stable_dual_identity_v1(
                kind,
                friendly_generators);
        }
    }
    return {};
}

std::string projected_node_identity(
    const MatExactGraphNode2& node,
    const MatExactGraph2& graph,
    const std::map<std::string, std::string>&
        site_names)
{
    if (node.generator_site_ids.size() == 3) {
        if (node.node_id
            != stable_voronoi_node_identity_v1(
                node.generator_site_ids)) {
            return {};
        }
        return stable_voronoi_node_identity_v1(
            normalized_ids(
                node.generator_site_ids,
                site_names));
    }
    if (node.generator_site_ids.size() != 2) {
        return {};
    }

    const MatExactGraphEdge2* owner = nullptr;
    const MatParameterEndpoint2* endpoint =
        nullptr;
    for (const MatExactGraphEdge2& edge :
         graph.edges) {
        const auto register_incidence =
            [&owner, &endpoint, &edge](
                const MatParameterEndpoint2&
                    candidate) {
                if (owner != nullptr) {
                    return false;
                }
                owner = &edge;
                endpoint = &candidate;
                return true;
            };
        if (edge.source_node_id == node.node_id
            && !register_incidence(
                edge.source_endpoint)) {
            return {};
        }
        if (edge.target_node_id == node.node_id
            && !register_incidence(
                edge.target_endpoint)) {
            return {};
        }
    }
    if (owner == nullptr || endpoint == nullptr
        || !endpoint->parameter.has_value()
        || node.generator_site_ids
            != owner->generator_site_ids
        || node.parent_site_ids
            != owner->parent_site_ids
        || node.node_id
            != stable_endpoint_node_identity_v1(
                owner->original_dual_id,
                *endpoint)) {
        return {};
    }
    const std::string projected_dual =
        projected_dual_identity(
            *owner,
            site_names);
    if (projected_dual.empty()) {
        return {};
    }
    return stable_endpoint_node_identity_v1(
        projected_dual,
        *endpoint);
}

std::map<std::string, std::string>
certificate_namespace(
    const MatExactGraph2& graph,
    const std::map<std::string, std::string>&
        site_names)
{
    std::map<std::string, std::string> names =
        site_names;
    for (const MatExactGraphEdge2& edge :
         graph.edges) {
        const std::string projected_dual =
            projected_dual_identity(
                edge,
                site_names);
        if (projected_dual.empty()
            || edge.edge_id
                != edge.original_dual_id
                    + "/component-0") {
            return {};
        }
        names.emplace(
            edge.original_dual_id,
            projected_dual);
        names.emplace(
            edge.edge_id,
            projected_dual + "/component-0");
    }
    for (const MatExactGraphNode2& node :
         graph.nodes) {
        const std::string projected_node =
            projected_node_identity(
                node,
                graph,
                site_names);
        if (projected_node.empty()) {
            return {};
        }
        names.emplace(
            node.node_id,
            projected_node);
    }
    return names;
}

bool derived_identities_are_valid(
    const MatExactGraph2& graph)
{
    std::set<std::string> node_ids;
    for (const MatExactGraphNode2& node :
         graph.nodes) {
        if (!node_ids.insert(node.node_id).second
            || !std::is_sorted(
                node.generator_site_ids.begin(),
                node.generator_site_ids.end())
            || !std::is_sorted(
                node.parent_site_ids.begin(),
                node.parent_site_ids.end())
            || projected_node_identity(
                   node,
                   graph,
                   {})
                .empty()) {
            return false;
        }
    }
    if (!std::is_sorted(
            graph.nodes.begin(),
            graph.nodes.end(),
            [](const MatExactGraphNode2& left,
               const MatExactGraphNode2& right) {
                return left.node_id
                    < right.node_id;
            })
        || !std::is_sorted(
            graph.edges.begin(),
            graph.edges.end(),
            [](const MatExactGraphEdge2& left,
               const MatExactGraphEdge2& right) {
                return left.edge_id
                    < right.edge_id;
            })) {
        return false;
    }

    for (const MatExactGraphEdge2& edge :
         graph.edges) {
        if (projected_dual_identity(edge, {})
                .empty()
            || edge.edge_id
                != edge.original_dual_id
                    + "/component-0"
            || edge.primitive_kind != "LINE"
            || edge.generator_site_ids
                != edge.parent_site_ids
            || node_ids.count(edge.source_node_id)
                != 1
            || node_ids.count(edge.target_node_id)
                != 1
            || !edge.source_endpoint.parameter.has_value()
            || !edge.target_endpoint.parameter.has_value()) {
            return false;
        }
        const std::string source_root =
            algebraic_root_identity_v1(
                *edge.source_endpoint.parameter);
        const std::string target_root =
            algebraic_root_identity_v1(
                *edge.target_endpoint.parameter);
        if (std::find(
                edge.source_endpoint.provenance_ids.begin(),
                edge.source_endpoint.provenance_ids.end(),
                source_root)
                == edge.source_endpoint.provenance_ids.end()
            || std::find(
                   edge.target_endpoint.provenance_ids.begin(),
                   edge.target_endpoint.provenance_ids.end(),
                   target_root)
                == edge.target_endpoint.provenance_ids.end()) {
            return false;
        }
    }
    return true;
}

std::vector<EdgeSignature> edge_signatures(
    const MatExactGraph2& graph,
    const std::map<std::string, std::string>&
        names)
{
    std::vector<EdgeSignature> signatures;
    signatures.reserve(graph.edges.size());
    for (const MatExactGraphEdge2& edge :
         graph.edges) {
        signatures.emplace_back(
            edge.primitive_kind,
            normalized_ids(
                edge.generator_site_ids,
                names),
            normalized_ids(
                edge.parent_site_ids,
                names),
            normalized_ids(
                {edge.original_dual_id},
                names)
                .front(),
            normalized_ids(
                {edge.edge_id},
                names)
                .front(),
            algebraic_root_identity_v1(
                *edge.source_endpoint.parameter),
            algebraic_root_identity_v1(
                *edge.target_endpoint.parameter),
            normalized_ids(
                edge.source_endpoint.provenance_ids,
                names),
            normalized_ids(
                edge.target_endpoint.provenance_ids,
                names));
    }
    std::sort(
        signatures.begin(),
        signatures.end());
    return signatures;
}

std::vector<NodeSignature> node_signatures(
    const MatExactGraph2& graph,
    const std::map<std::string, std::string>&
        names)
{
    std::map<std::string, std::size_t> degree;
    for (const MatExactGraphEdge2& edge :
         graph.edges) {
        ++degree[edge.source_node_id];
        ++degree[edge.target_node_id];
    }
    std::vector<NodeSignature> signatures;
    signatures.reserve(graph.nodes.size());
    for (const MatExactGraphNode2& node :
         graph.nodes) {
        signatures.emplace_back(
            normalized_ids(
                {node.node_id},
                names)
                .front(),
            normalized_ids(
                node.generator_site_ids,
                names),
            normalized_ids(
                node.parent_site_ids,
                names),
            normalized_ids(
                node.provenance_ids,
                names),
            degree[node.node_id]);
    }
    std::sort(
        signatures.begin(),
        signatures.end());
    return signatures;
}

bool equivalent_to_validated_fixture(
    const MatExactGraph2& catalog_graph,
    const MatExactGraph2& fixture_graph,
    const std::map<std::string, std::string>&
        names)
{
    const auto namespace_names =
        certificate_namespace(
            catalog_graph,
            names);
    const bool namespace_valid =
        !namespace_names.empty();
    const bool rejected_equal =
        catalog_graph.rejected_incident_transitions
        == fixture_graph.rejected_incident_transitions;
    const bool matched_equal =
        catalog_graph.matched_generator_sites
        == fixture_graph.matched_generator_sites;
    const bool edges_equal =
        edge_signatures(
            catalog_graph,
            namespace_names)
        == edge_signatures(
            fixture_graph,
            {});
    const bool nodes_equal =
        node_signatures(
            catalog_graph,
            namespace_names)
        == node_signatures(
            fixture_graph,
            {});
    return namespace_valid && rejected_equal
        && matched_equal && edges_equal
        && nodes_equal;
}

bool catalog_graph_is_exact_and_invariant()
{
    const CanonicalReachInput2 input =
        rectangle_input(false);
    const CanonicalReachInput2 transformed =
        rectangle_input(true);
    const CanonicalMatSiteCatalog2 catalog =
        canonical_mat_site_catalog(input);
    const auto names =
        friendly_site_names(catalog);
    if (names.size() != 8) {
        return false;
    }
    const MatExactGraph2 radius_zero =
        canonical_rectangle_mat_graph(
            input,
            0);
    const MatExactGraph2 repeated =
        canonical_rectangle_mat_graph(
            input,
            0);
    const MatExactGraph2 reversed =
        canonical_rectangle_mat_graph(
            transformed,
            0);
    const MatExactGraph2 radius_one =
        canonical_rectangle_mat_graph(
            input,
            1);
    const MatExactGraph2 plateau =
        canonical_rectangle_mat_graph(
            input,
            4);
    const MatExactGraph2 rejected =
        canonical_rectangle_mat_graph(
            input,
            5);

    bool negative_rejected = false;
    try {
        static_cast<void>(
            canonical_rectangle_mat_graph(
                input,
                -1));
    } catch (
        const NegativeClearanceRadiusSquaredError&) {
        negative_rejected = true;
    }

    const bool repeat_equal =
        graphs_equal(radius_zero, repeated);
    const bool symmetry_equal =
        graphs_equal(radius_zero, reversed);
    const bool ids_zero =
        derived_identities_are_valid(radius_zero);
    const bool ids_one =
        derived_identities_are_valid(radius_one);
    const bool ids_plateau =
        derived_identities_are_valid(plateau);
    const bool ids_rejected =
        derived_identities_are_valid(rejected);
    const bool fixture_zero =
        equivalent_to_validated_fixture(
            radius_zero,
            segment_site_rectangle_graph_spike(0),
            names);
    const bool fixture_one =
        equivalent_to_validated_fixture(
            radius_one,
            segment_site_rectangle_graph_spike(1),
            names);
    const bool fixture_plateau =
        equivalent_to_validated_fixture(
            plateau,
            segment_site_rectangle_graph_spike(4),
            names);
    const bool fixture_rejected =
        equivalent_to_validated_fixture(
            rejected,
            segment_site_rectangle_graph_spike(5),
            names);
    return negative_rejected
        && radius_zero.edges.size() == 5
        && radius_zero.nodes.size() == 6
        && radius_one.edges.size() == 5
        && radius_one.nodes.size() == 6
        && plateau.edges.size() == 1
        && plateau.nodes.size() == 2
        && rejected.edges.empty()
        && rejected.nodes.empty()
        && repeat_equal
        && symmetry_equal
        && ids_zero
        && ids_one
        && ids_plateau
        && ids_rejected
        && fixture_zero
        && fixture_one
        && fixture_plateau
        && fixture_rejected;
}

bool unsupported_rectangle_inputs_are_rejected()
{
    const CanonicalReachInput2 triangle =
        polygon_input(
            std::array<
                std::array<double, 2>,
                3>{{
                {-4.0, -2.0},
                {4.0, -2.0},
                {0.0, 2.0},
            }});
    const CanonicalReachInput2 trapezoid =
        polygon_input(
            std::array<
                std::array<double, 2>,
                4>{{
                {-4.0, -2.0},
                {4.0, -2.0},
                {3.0, 2.0},
                {-3.0, 2.0},
            }});
    const auto is_rejected =
        [](const CanonicalReachInput2& input) {
            try {
                static_cast<void>(
                    canonical_rectangle_mat_graph(
                        input,
                        0));
            } catch (
                const UnsupportedCanonicalMatRectangleGraphError&) {
                return true;
            }
            return false;
        };
    return is_rejected(triangle)
        && is_rejected(trapezoid);
}

} // namespace

bool catalog_graph_gate()
{
    const bool exact_and_invariant =
        catalog_graph_is_exact_and_invariant();
    const bool unsupported_rejected =
        unsupported_rectangle_inputs_are_rejected();
    return exact_and_invariant
        && unsupported_rejected;
}
