#include "segment_site_graph_csr.h"
#include "segment_site_catalog.h"
#include "segment_site_catalog_graph.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace {

compas::RowMatrixXd matrix(
    const std::vector<std::array<double, 2>>& points)
{
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
        transformed
            ? matrix(
                  {
                      {4.0, 2.0},
                      {4.0, -2.0},
                      {-4.0, -2.0},
                      {-4.0, 2.0},
                  })
            : matrix(
                  {
                      {-4.0, -2.0},
                      {4.0, -2.0},
                      {4.0, 2.0},
                      {-4.0, 2.0},
                  }),
        {},
        1.0);
}

CanonicalReachInput2 l_shape_input(
    const bool transformed)
{
    return canonical_reach_input(
        transformed
            ? matrix(
                  {
                      {0.0, 6.0},
                      {2.0, 6.0},
                      {2.0, 2.0},
                      {6.0, 2.0},
                      {6.0, 0.0},
                      {0.0, 0.0},
                  })
            : matrix(
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

std::size_t endpoint_flag_count(
    const MatNumericGraphTable2& table,
    const std::int64_t flag)
{
    std::size_t count = 0;
    for (const auto& edge_flags :
         table.edge_endpoint_provenance_flags) {
        count +=
            (edge_flags[0] & flag) != 0;
        count +=
            (edge_flags[1] & flag) != 0;
    }
    return count;
}

std::size_t feature_row_count(
    const MatNumericGraphTable2& table,
    const std::int64_t domain_kind,
    const std::int64_t curve_kind)
{
    return static_cast<std::size_t>(
        std::count_if(
            table.endpoint_features.begin(),
            table.endpoint_features.end(),
            [domain_kind, curve_kind](
                const auto& row) {
                return row[0] == domain_kind
                    && row[2] == curve_kind;
            }));
}

bool endpoint_feature_csr_is_canonical(
    const MatNumericGraphTable2& table)
{
    if (table.endpoint_feature_offsets.size()
            != 2 * table.edges.size() + 1
        || table.endpoint_feature_offsets.empty()
        || table.endpoint_feature_offsets.front() != 0
        || table.endpoint_feature_offsets.back()
            != static_cast<std::int64_t>(
                table.endpoint_features.size())) {
        return false;
    }
    for (std::size_t endpoint = 0;
         endpoint + 1
             < table.endpoint_feature_offsets.size();
         ++endpoint) {
        const auto first =
            table.endpoint_features.begin()
            + table.endpoint_feature_offsets[endpoint];
        const auto last =
            table.endpoint_features.begin()
            + table.endpoint_feature_offsets[
                endpoint + 1];
        if (!std::is_sorted(first, last)
            || std::adjacent_find(first, last)
                != last) {
            return false;
        }
    }
    return true;
}

bool rectangle_public_table_is_exact()
{
    const CanonicalReachInput2 input =
        rectangle_input(false);
    const CanonicalMatSiteCatalog2 catalog =
        canonical_mat_site_catalog(input);
    const MatNumericGraphTable2 zero =
        numeric_graph_table(
            canonical_rectangle_mat_graph(
                input,
                CORE::BigRat(0)),
            catalog);
    const MatNumericGraphTable2 positive =
        numeric_graph_table(
            canonical_rectangle_mat_graph(
                input,
                CORE::BigRat(1)),
            catalog);

    const std::array<std::int64_t, 3>
        verified{1, 1, 1};
    return zero.node_ids.size() == 6
        && zero.edges.size() == 5
        && zero.original_dual_ids.size() == 5
        && zero.node_site_offsets.size() == 7
        && zero.node_site_ids.size() == 18
        && zero.site_provenance
            == catalog.site_provenance()
        && zero.edge_endpoint_provenance_flags.size()
            == 5
        && endpoint_flag_count(
               zero,
               MAT_ENDPOINT_ORIGINAL_VORONOI_VERTEX)
            == 10
        && endpoint_flag_count(
               zero,
               MAT_ENDPOINT_DOMAIN_BOUNDARY)
            == 4
        && endpoint_flag_count(
               zero,
               MAT_ENDPOINT_CLEARANCE_ROOT)
            == 4
        && zero.endpoint_features.size() == 20
        && feature_row_count(
               zero,
               MAT_ENDPOINT_DOMAIN_DESIGN,
               MAT_CURVE_LINE)
            == 8
        && feature_row_count(
               zero,
               MAT_ENDPOINT_DOMAIN_CLEARANCE,
               MAT_CURVE_LINE)
            == 8
        && feature_row_count(
               zero,
               MAT_ENDPOINT_DOMAIN_CLEARANCE,
               MAT_CURVE_CIRCLE)
            == 4
        && endpoint_feature_csr_is_canonical(zero)
        && std::all_of(
               zero.edge_exact_flags.begin(),
               zero.edge_exact_flags.end(),
               [&verified](const auto& flags) {
                   return flags == verified;
               })
        && positive.node_ids.size() == 6
        && positive.edges.size() == 5
        && positive.node_site_ids.size() == 14
        && endpoint_flag_count(
               positive,
               MAT_ENDPOINT_ORIGINAL_VORONOI_VERTEX)
            == 6
        && endpoint_flag_count(
               positive,
               MAT_ENDPOINT_DOMAIN_BOUNDARY)
            == 0
        && endpoint_flag_count(
               positive,
               MAT_ENDPOINT_CLEARANCE_ROOT)
            == 4
        && positive.endpoint_features.size() == 8
        && feature_row_count(
               positive,
               MAT_ENDPOINT_DOMAIN_CLEARANCE,
               MAT_CURVE_LINE)
            == 8
        && endpoint_feature_csr_is_canonical(positive);
}

bool l_shape_public_table_is_canonical()
{
    const CanonicalReachInput2 canonical =
        l_shape_input(false);
    const CanonicalReachInput2 transformed =
        l_shape_input(true);
    const MatNumericGraphTable2 first =
        numeric_graph_table(
            canonical_l_shape_mat_graph(
                canonical,
                CORE::BigRat(0)),
            canonical_mat_site_catalog(canonical));
    const MatNumericGraphTable2 second =
        numeric_graph_table(
            canonical_l_shape_mat_graph(
                transformed,
                CORE::BigRat(0)),
            canonical_mat_site_catalog(transformed));
    std::size_t lines = 0;
    std::size_t parabolas = 0;
    std::size_t point_segment = 0;
    std::size_t segment_segment = 0;
    for (const auto& edge : first.edges) {
        lines += edge[2] == MAT_CURVE_LINE;
        parabolas +=
            edge[2] == MAT_CURVE_PARABOLA;
        point_segment +=
            edge[5]
            == MAT_DUAL_POINT_SEGMENT;
        segment_segment +=
            edge[5]
            == MAT_DUAL_SEGMENT_SEGMENT;
    }
    return first == second
        && first.node_ids.size() == 10
        && first.edges.size() == 9
        && first.original_dual_ids.size() == 9
        && first.site_provenance.size() == 12
        && lines == 7
        && parabolas == 2
        && point_segment == 2
        && segment_segment == 7
        && endpoint_flag_count(
               first,
               MAT_ENDPOINT_ORIGINAL_VORONOI_VERTEX)
            == 18
        && endpoint_feature_csr_is_canonical(first);
}

bool malformed_public_graph_is_rejected()
{
    const CanonicalReachInput2 input =
        rectangle_input(false);
    const CanonicalMatSiteCatalog2 catalog =
        canonical_mat_site_catalog(input);
    const MatExactGraph2 graph =
        canonical_rectangle_mat_graph(
            input,
            CORE::BigRat(1));

    bool edge_order_rejected = false;
    try {
        MatExactGraph2 reversed = graph;
        std::swap(
            reversed.edges[0],
            reversed.edges[1]);
        static_cast<void>(
            numeric_graph_table(
                reversed,
                catalog));
    } catch (
        const NonCanonicalMatGraphEdgeOrderError&) {
        edge_order_rejected = true;
    }

    bool node_reference_rejected = false;
    try {
        MatExactGraph2 unknown = graph;
        unknown.edges[0].source_node_id =
            "unknown-node";
        static_cast<void>(
            numeric_graph_table(
                unknown,
                catalog));
    } catch (
        const UnknownMatGraphNodeIdentityError&) {
        node_reference_rejected = true;
    }

    bool endpoint_evidence_rejected = false;
    try {
        MatExactGraph2 unevidenced = graph;
        unevidenced.edges[0]
            .source_endpoint.exact_evidence = {};
        static_cast<void>(
            numeric_graph_table(
                unevidenced,
                catalog));
    } catch (
        const InvalidMatEndpointEvidenceError&) {
        endpoint_evidence_rejected = true;
    }

    bool component_rejected = false;
    try {
        MatExactGraph2 invalid = graph;
        invalid.edges[0].clip_component_index = -1;
        static_cast<void>(
            numeric_graph_table(
                invalid,
                catalog));
    } catch (
        const InvalidMatGraphEdgeRecordError&) {
        component_rejected = true;
    }

    bool admissibility_rejected = false;
    try {
        MatExactGraph2 unverified = graph;
        unverified.edges[0]
            .admissible_center_component = false;
        static_cast<void>(
            numeric_graph_table(
                unverified,
                catalog));
    } catch (
        const InvalidMatGraphEdgeRecordError&) {
        admissibility_rejected = true;
    }

    return edge_order_rejected
        && node_reference_rejected
        && endpoint_evidence_rejected
        && component_rejected
        && admissibility_rejected;
}

} // namespace

bool graph_csr_gate()
{
    const MatExactGraph2 graph{
        {
            {
                "node-a",
                {},
                {"corner", "left-segment"},
                {"left-segment"},
            },
            {
                "node-b",
                {},
                {
                    "bottom-segment",
                    "right-segment",
                    "top-segment",
                },
                {
                    "bottom-segment",
                    "right-segment",
                    "top-segment",
                },
            },
        },
        {},
        0,
        5,
    };
    const MatNodeSiteCsr2 csr =
        node_site_csr(graph);
    const MatNodeSiteCsr2 empty =
        node_site_csr({{}, {}, 0, 0});
    const MatNodeSiteCsr2 rectangle =
        node_site_csr(
            segment_site_rectangle_graph_spike());

    bool node_order_rejected = false;
    try {
        MatExactGraph2 reversed = graph;
        std::swap(
            reversed.nodes[0],
            reversed.nodes[1]);
        static_cast<void>(
            node_site_csr(reversed));
    } catch (
        const NonCanonicalMatGraphNodeOrderError&) {
        node_order_rejected = true;
    }
    bool empty_node_id_rejected = false;
    try {
        MatExactGraph2 empty_node_id = graph;
        empty_node_id.nodes[0].node_id.clear();
        static_cast<void>(
            node_site_csr(empty_node_id));
    } catch (
        const NonCanonicalMatGraphNodeOrderError&) {
        empty_node_id_rejected = true;
    }
    bool site_order_rejected = false;
    try {
        MatExactGraph2 duplicated = graph;
        duplicated.nodes[0].generator_site_ids = {
            "corner",
            "corner",
        };
        static_cast<void>(
            node_site_csr(duplicated));
    } catch (
        const NonCanonicalMatNodeSiteOrderError&) {
        site_order_rejected = true;
    }
    bool unsorted_sites_rejected = false;
    try {
        MatExactGraph2 unsorted = graph;
        unsorted.nodes[0].generator_site_ids = {
            "left-segment",
            "corner",
        };
        static_cast<void>(
            node_site_csr(unsorted));
    } catch (
        const NonCanonicalMatNodeSiteOrderError&) {
        unsorted_sites_rejected = true;
    }
    bool empty_sites_rejected = false;
    try {
        MatExactGraph2 empty_sites = graph;
        empty_sites.nodes[0]
            .generator_site_ids.clear();
        static_cast<void>(
            node_site_csr(empty_sites));
    } catch (
        const NonCanonicalMatNodeSiteOrderError&) {
        empty_sites_rejected = true;
    }

    return csr.node_site_offsets
            == std::vector<std::int64_t>{
                0,
                2,
                5,
            }
        && csr.node_site_ids
            == std::vector<std::string>{
                "corner",
                "left-segment",
                "bottom-segment",
                "right-segment",
                "top-segment",
            }
        && empty.node_site_offsets
            == std::vector<std::int64_t>{0}
        && empty.node_site_ids.empty()
        && rectangle.node_site_offsets
            == std::vector<std::int64_t>{
                0,
                3,
                6,
                9,
                12,
                15,
                18,
            }
        && rectangle.node_site_ids
            == std::vector<std::string>{
                "left-segment",
                "top-segment",
                "upper-left",
                "right-segment",
                "top-segment",
                "upper-right",
                "bottom-segment",
                "lower-right",
                "right-segment",
                "bottom-segment",
                "left-segment",
                "lower-left",
                "bottom-segment",
                "left-segment",
                "top-segment",
                "bottom-segment",
                "right-segment",
                "top-segment",
            }
        && node_order_rejected
        && empty_node_id_rejected
        && site_order_rejected
        && unsorted_sites_rejected
        && empty_sites_rejected
        && rectangle_public_table_is_exact()
        && l_shape_public_table_is_canonical()
        && malformed_public_graph_is_rejected();
}
