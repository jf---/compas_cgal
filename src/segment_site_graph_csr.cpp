#include "segment_site_graph_csr.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>

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
