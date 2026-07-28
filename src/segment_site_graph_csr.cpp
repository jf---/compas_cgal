#include "segment_site_graph_csr.h"

#include "segment_site_catalog.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <string_view>
#include <utility>

namespace {

void require_canonical_node_order(
    const std::vector<MatExactGraphNode2>& nodes)
{
    for (std::size_t index = 0;
         index < nodes.size();
         ++index) {
        if (nodes[index].node_id.empty()
            || (index > 0
                && nodes[index].node_id
                    <= nodes[index - 1].node_id)) {
            throw NonCanonicalMatGraphNodeOrderError(
                "MAT graph nodes are not strictly ordered by stable identity");
        }
    }
}

void require_canonical_site_order(
    const MatExactGraphNode2& node)
{
    if (node.generator_site_ids.empty()
        || std::any_of(
            node.generator_site_ids.begin(),
            node.generator_site_ids.end(),
            [](const std::string& site_id) {
                return site_id.empty();
            })
        || !std::is_sorted(
            node.generator_site_ids.begin(),
            node.generator_site_ids.end())
        || std::adjacent_find(
               node.generator_site_ids.begin(),
               node.generator_site_ids.end())
            != node.generator_site_ids.end()) {
        throw NonCanonicalMatNodeSiteOrderError(
            "MAT node sites are not nonempty, unique, and ordered");
    }
}

std::int64_t checked_table_index(
    const std::size_t index,
    const std::string_view context)
{
    if (static_cast<std::uintmax_t>(index)
        > static_cast<std::uintmax_t>(
            std::numeric_limits<std::int64_t>::max())) {
        throw MatNumericGraphTableOverflowError(
            std::string(context)
            + " exceeds int64 range");
    }
    return static_cast<std::int64_t>(index);
}

void require_canonical_edge_order(
    const std::vector<MatExactGraphEdge2>& edges)
{
    for (std::size_t index = 0;
         index < edges.size();
         ++index) {
        if (edges[index].edge_id.empty()
            || (index > 0
                && edges[index].edge_id
                    <= edges[index - 1].edge_id)) {
            throw NonCanonicalMatGraphEdgeOrderError(
                "MAT graph edges are not strictly ordered by stable identity");
        }
    }
}

std::int64_t node_index(
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
        throw UnknownMatGraphNodeIdentityError(
            "MAT edge references a node outside the canonical graph");
    }
    return checked_table_index(
        static_cast<std::size_t>(
            std::distance(
                graph.nodes.begin(),
                found)),
        "MAT node index");
}

std::int64_t curve_kind(
    const std::string& primitive_kind)
{
    if (primitive_kind == "LINE") {
        return MAT_CURVE_LINE;
    }
    if (primitive_kind == "PARABOLA") {
        return MAT_CURVE_PARABOLA;
    }
    throw UnsupportedMatGraphCurveKindError(
        "MAT graph edge has an unsupported curve kind");
}

std::int64_t original_dual_kind(
    const std::array<std::int64_t, 2>& site_indices,
    const CanonicalMatSiteCatalog2& catalog)
{
    const MatSiteKind2 first =
        catalog.sites()[
            static_cast<std::size_t>(
                site_indices[0])]
            .provenance()
            .kind;
    const MatSiteKind2 second =
        catalog.sites()[
            static_cast<std::size_t>(
                site_indices[1])]
            .provenance()
            .kind;
    if (first == MatSiteKind2::Point
        && second == MatSiteKind2::Point) {
        return MAT_DUAL_POINT_POINT;
    }
    if (first == MatSiteKind2::OpenSegment
        && second == MatSiteKind2::OpenSegment) {
        return MAT_DUAL_SEGMENT_SEGMENT;
    }
    return MAT_DUAL_POINT_SEGMENT;
}

std::int64_t original_dual_index(
    const std::vector<std::string>& identities,
    const std::string& identity)
{
    const auto found = std::lower_bound(
        identities.begin(),
        identities.end(),
        identity);
    if (found == identities.end()
        || *found != identity) {
        throw InvalidMatGraphEdgeRecordError(
            "MAT edge has no canonical original-dual identity");
    }
    return checked_table_index(
        static_cast<std::size_t>(
            std::distance(
                identities.begin(),
                found)),
        "MAT original-dual index");
}

void require_endpoint_provenance(
    const MatParameterEndpoint2& endpoint,
    const MatExactGraphNode2& node)
{
    if (!endpoint.parameter.has_value()
        || endpoint.provenance_ids.empty()
        || !std::is_sorted(
            endpoint.provenance_ids.begin(),
            endpoint.provenance_ids.end())
        || std::adjacent_find(
               endpoint.provenance_ids.begin(),
               endpoint.provenance_ids.end())
            != endpoint.provenance_ids.end()
        || !std::is_sorted(
            node.provenance_ids.begin(),
            node.provenance_ids.end())
        || std::adjacent_find(
               node.provenance_ids.begin(),
               node.provenance_ids.end())
            != node.provenance_ids.end()
        || !std::includes(
            node.provenance_ids.begin(),
            node.provenance_ids.end(),
            endpoint.provenance_ids.begin(),
            endpoint.provenance_ids.end())) {
        throw InvalidMatEndpointEvidenceError(
            "MAT endpoint provenance is not contained in its exact node");
    }
    const MatEndpointExactEvidence2& evidence =
        endpoint.exact_evidence;
    const std::string root_id =
        algebraic_root_identity_v1(
            *endpoint.parameter);
    if ((!evidence.original_voronoi_vertex
         && !evidence.domain_boundary
         && !evidence.clearance_root)
        || (evidence.domain_boundary
            != !evidence.boundary_features.empty())
        || !std::is_sorted(
            evidence.boundary_features.begin(),
            evidence.boundary_features.end())
        || std::adjacent_find(
               evidence.boundary_features.begin(),
               evidence.boundary_features.end())
            != evidence.boundary_features.end()
        || !std::binary_search(
            endpoint.provenance_ids.begin(),
            endpoint.provenance_ids.end(),
            root_id)
        || (evidence.original_voronoi_vertex
            && node.generator_site_ids.size() < 3)) {
        throw InvalidMatEndpointEvidenceError(
            "MAT endpoint exact-event evidence is incomplete");
    }
}

std::int64_t endpoint_flags(
    const MatParameterEndpoint2& endpoint)
{
    std::int64_t flags = 0;
    if (endpoint.exact_evidence
            .original_voronoi_vertex) {
        flags |=
            MAT_ENDPOINT_ORIGINAL_VORONOI_VERTEX;
    }
    if (endpoint.exact_evidence.domain_boundary) {
        flags |= MAT_ENDPOINT_DOMAIN_BOUNDARY;
    }
    if (endpoint.exact_evidence.clearance_root) {
        flags |= MAT_ENDPOINT_CLEARANCE_ROOT;
    }
    return flags;
}

std::array<std::int64_t, 5>
numeric_boundary_feature(
    const MatEndpointBoundaryFeature2& feature,
    const CanonicalMatSiteCatalog2& catalog)
{
    if (feature.domain_kind
            != MatEndpointDomainKind2::Design
        || feature.component != 0
        || feature.curve_kind
            != MatExactCurveKind2::Line
        || feature.source_site_or_ring < 0
        || feature.derived_feature_index < 0
        || std::none_of(
            catalog.sites().begin(),
            catalog.sites().end(),
            [&feature](
                const CanonicalMatSite2& site) {
                const MatSiteProvenance2&
                    provenance =
                        site.provenance();
                return provenance.kind
                        == MatSiteKind2::OpenSegment
                    && provenance.ring
                        == feature
                               .source_site_or_ring
                    && provenance.feature
                        == feature
                               .derived_feature_index;
            })) {
        throw InvalidMatEndpointEvidenceError(
            "MAT design-boundary feature has an invalid numeric contract");
    }
    return {
        static_cast<std::int64_t>(
            feature.domain_kind),
        feature.component,
        static_cast<std::int64_t>(
            feature.curve_kind),
        feature.source_site_or_ring,
        feature.derived_feature_index,
    };
}

std::vector<std::array<std::int64_t, 5>>
numeric_endpoint_features(
    const MatParameterEndpoint2& endpoint,
    const MatExactGraphEdge2& edge,
    const MatExactGraphNode2& node,
    const CanonicalMatSiteCatalog2& catalog)
{
    std::vector<std::array<std::int64_t, 5>>
        features;
    if (endpoint.exact_evidence
            .boundary_features.size()
        > std::numeric_limits<std::size_t>::max()
            - node.generator_site_ids.size()) {
        throw MatNumericGraphTableOverflowError(
            "MAT endpoint-feature reserve size overflows");
    }
    features.reserve(
        endpoint.exact_evidence
            .boundary_features.size()
        + node.generator_site_ids.size());
    for (const MatEndpointBoundaryFeature2& feature :
         endpoint.exact_evidence.boundary_features) {
        features.push_back(
            numeric_boundary_feature(
                feature,
                catalog));
    }
    if (endpoint.exact_evidence.clearance_root) {
        const auto& contact_site_ids =
            endpoint.exact_evidence
                    .original_voronoi_vertex
            ? node.generator_site_ids
            : edge.generator_site_ids;
        for (const std::string& site_id :
             contact_site_ids) {
            const std::int64_t site_index =
                catalog.index_of(site_id);
            const MatSiteKind2 site_kind =
                catalog.sites()[
                    static_cast<std::size_t>(
                        site_index)]
                    .provenance()
                    .kind;
            features.push_back(
                {
                    MAT_ENDPOINT_DOMAIN_CLEARANCE,
                    edge.clip_component_index,
                    site_kind == MatSiteKind2::Point
                        ? MAT_CURVE_CIRCLE
                        : MAT_CURVE_LINE,
                    site_index,
                    0,
                });
        }
    }
    std::sort(features.begin(), features.end());
    features.erase(
        std::unique(
            features.begin(),
            features.end()),
        features.end());
    return features;
}

} // namespace

MatNodeSiteCsr2
node_site_csr(const MatExactGraph2& graph)
{
    require_canonical_node_order(graph.nodes);
    MatNodeSiteCsr2 csr{{0}, {}};
    const std::uintmax_t max_offset =
        static_cast<std::uintmax_t>(
            std::numeric_limits<std::int64_t>::max());
    for (const MatExactGraphNode2& node :
         graph.nodes) {
        require_canonical_site_order(node);
        const std::uintmax_t retained_size =
            static_cast<std::uintmax_t>(
                csr.node_site_ids.size());
        const std::uintmax_t appended_size =
            static_cast<std::uintmax_t>(
                node.generator_site_ids.size());
        if (retained_size > max_offset
            || appended_size
                > max_offset - retained_size) {
            throw MatNodeSiteCsrOverflowError(
                "MAT node-site CSR exceeds int64 offset range");
        }
        csr.node_site_ids.insert(
            csr.node_site_ids.end(),
            node.generator_site_ids.begin(),
            node.generator_site_ids.end());
        csr.node_site_offsets.push_back(
            static_cast<std::int64_t>(
                csr.node_site_ids.size()));
    }
    return csr;
}

MatNumericNodeSiteCsr2
numeric_node_site_csr(
    const MatExactGraph2& graph,
    const CanonicalMatSiteCatalog2& catalog)
{
    MatNodeSiteCsr2 stable =
        node_site_csr(graph);
    MatNumericNodeSiteCsr2 numeric{
        std::move(stable.node_site_offsets),
        {},
    };
    numeric.node_site_ids.reserve(
        stable.node_site_ids.size());
    for (const std::string& stable_id :
         stable.node_site_ids) {
        numeric.node_site_ids.push_back(
            catalog.index_of(stable_id));
    }
    return numeric;
}

MatNumericGraphTable2 numeric_graph_table(
    const MatExactGraph2& graph,
    const CanonicalMatSiteCatalog2& catalog)
{
    const MatNumericNodeSiteCsr2 node_sites =
        numeric_node_site_csr(
            graph,
            catalog);
    require_canonical_edge_order(graph.edges);

    MatNumericGraphTable2 table;
    table.node_site_offsets =
        node_sites.node_site_offsets;
    table.node_site_ids =
        node_sites.node_site_ids;
    table.site_provenance =
        catalog.site_provenance();
    table.node_ids.reserve(graph.nodes.size());
    for (const MatExactGraphNode2& node :
         graph.nodes) {
        table.node_ids.push_back(node.node_id);
    }
    table.original_dual_ids.reserve(
        graph.edges.size());
    for (const MatExactGraphEdge2& edge :
         graph.edges) {
        if (edge.original_dual_id.empty()) {
            throw InvalidMatGraphEdgeRecordError(
                "MAT edge original-dual identity is empty");
        }
        table.original_dual_ids.push_back(
            edge.original_dual_id);
    }
    std::sort(
        table.original_dual_ids.begin(),
        table.original_dual_ids.end());
    table.original_dual_ids.erase(
        std::unique(
            table.original_dual_ids.begin(),
            table.original_dual_ids.end()),
        table.original_dual_ids.end());

    table.edges.reserve(graph.edges.size());
    table.edge_endpoint_provenance_flags.reserve(
        graph.edges.size());
    if (graph.edges.size()
        > (std::numeric_limits<std::size_t>::max()
               - 1)
            / 2) {
        throw MatNumericGraphTableOverflowError(
            "MAT endpoint-feature offset count overflows");
    }
    table.endpoint_feature_offsets.reserve(
        2 * graph.edges.size() + 1);
    table.endpoint_feature_offsets.push_back(0);
    table.edge_exact_flags.reserve(
        graph.edges.size());

    for (const MatExactGraphEdge2& edge :
         graph.edges) {
        if (edge.source_node_id
                == edge.target_node_id
            || edge.clip_component_index < 0
            || !edge.admissible_center_component
            || edge.generator_site_ids.size() != 2
            || edge.generator_site_ids[1]
                <= edge.generator_site_ids[0]) {
            throw InvalidMatGraphEdgeRecordError(
                "MAT edge topology or generator contract is invalid");
        }
        const std::int64_t source =
            node_index(
                graph,
                edge.source_node_id);
        const std::int64_t target =
            node_index(
                graph,
                edge.target_node_id);
        const auto& source_node =
            graph.nodes[
                static_cast<std::size_t>(source)];
        const auto& target_node =
            graph.nodes[
                static_cast<std::size_t>(target)];
        require_endpoint_provenance(
            edge.source_endpoint,
            source_node);
        require_endpoint_provenance(
            edge.target_endpoint,
            target_node);

        const std::array<std::int64_t, 2>
            site_indices{
                catalog.index_of(
                    edge.generator_site_ids[0]),
                catalog.index_of(
                    edge.generator_site_ids[1]),
            };
        table.edges.push_back(
            {
                source,
                target,
                curve_kind(
                    edge.primitive_kind),
                site_indices[0],
                site_indices[1],
                original_dual_kind(
                    site_indices,
                    catalog),
                original_dual_index(
                    table.original_dual_ids,
                    edge.original_dual_id),
                edge.clip_component_index,
            });
        table.edge_endpoint_provenance_flags
            .push_back(
                {
                    endpoint_flags(
                        edge.source_endpoint),
                    endpoint_flags(
                        edge.target_endpoint),
                });

        const auto append_features =
            [&table](
                std::vector<
                    std::array<std::int64_t, 5>>
                    features) {
                table.endpoint_features.insert(
                    table.endpoint_features.end(),
                    std::make_move_iterator(
                        features.begin()),
                    std::make_move_iterator(
                        features.end()));
                table.endpoint_feature_offsets
                    .push_back(
                        checked_table_index(
                            table.endpoint_features
                                .size(),
                            "MAT endpoint-feature offset"));
            };
        append_features(
            numeric_endpoint_features(
                edge.source_endpoint,
                edge,
                source_node,
                catalog));
        append_features(
            numeric_endpoint_features(
                edge.target_endpoint,
                edge,
                target_node,
                catalog));
        table.edge_exact_flags.push_back(
            {1, 1, 1});
    }
    return table;
}
