#include "segment_site_graph_csr.h"

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

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
        && empty_sites_rejected;
}
