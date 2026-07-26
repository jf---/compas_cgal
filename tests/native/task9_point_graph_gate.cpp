#include "segment_site_mat.h"

#include <algorithm>
#include <string>
#include <vector>

namespace
{

template <typename Error, typename Callable>
bool throws_named(Callable&& callable)
{
    try
    {
        callable();
    }
    catch (const Error&)
    {
        return true;
    }
    return false;
}

MatDomainPolygonWithHoles2 production_domain()
{
    MatDomainPolygon2 outer;
    outer.push_back({-8, -8});
    outer.push_back({8, -8});
    outer.push_back({8, 8});
    outer.push_back({-8, 8});
    return MatDomainPolygonWithHoles2(outer);
}

std::vector<NormalizedPointSource2> production_points()
{
    return {
        {"site-delta", -3, -2},
        {"site-alpha", 4, -1},
        {"site-charlie", 1, 5},
        {"site-bravo", -1, 1},
    };
}

std::string graph_signature(const MatExactGraph2& graph)
{
    std::string signature;
    for (const MatExactGraphNode2& node : graph.nodes)
    {
        signature += node.node_id;
        for (const std::string& site_id : node.generator_site_ids)
        {
            signature += "/" + site_id;
        }
        signature += ";";
    }
    signature += "|";
    for (const MatExactGraphEdge2& edge : graph.edges)
    {
        signature += edge.edge_id + "/" + edge.primitive_kind
            + "/" + edge.source_node_id + "/" + edge.target_node_id;
        for (const std::string& site_id : edge.generator_site_ids)
        {
            signature += "/" + site_id;
        }
        signature += ";";
    }
    return signature;
}

} // namespace

bool point_graph_production_gate()
{
    const MatDomainPolygonWithHoles2 domain = production_domain();
    const std::vector<NormalizedPointSource2> points =
        production_points();
    const MatExactGraph2 first =
        exact_point_site_graph(points, domain, 0);
    const MatExactGraph2 second =
        exact_point_site_graph(points, domain, 0);
    std::vector<NormalizedPointSource2> reversed_points = points;
    std::reverse(
        reversed_points.begin(),
        reversed_points.end());
    const MatExactGraph2 reversed =
        exact_point_site_graph(reversed_points, domain, 0);

    // The square-domain diameter squared is 512, so 1024 excludes every
    // admissible-center point regardless of which in-domain site defines it.
    const CORE::BigRat outside_domain_radius_squared = 1024;
    const MatExactGraph2 excluded =
        exact_point_site_graph(
            points,
            domain,
            outside_domain_radius_squared);

    std::vector<NormalizedPointSource2> duplicate_identity = points;
    duplicate_identity.back().stable_site_id =
        duplicate_identity.front().stable_site_id;
    std::vector<NormalizedPointSource2> coincident = points;
    coincident.back().x = coincident.front().x;
    coincident.back().y = coincident.front().y;
    const std::vector<NormalizedPointSource2> cocircular{
        {"square-a", -1, -1},
        {"square-b", 1, -1},
        {"square-c", 1, 1},
        {"square-d", -1, 1},
    };

    return !first.nodes.empty()
        && !first.edges.empty()
        && graph_signature(first) == graph_signature(second)
        && graph_signature(first) == graph_signature(reversed)
        && excluded.nodes.empty()
        && excluded.edges.empty()
        && throws_named<DuplicatePointSiteIdentityError>(
            [&duplicate_identity, &domain]()
            {
                exact_point_site_graph(
                    duplicate_identity,
                    domain,
                    0);
            })
        && throws_named<CoincidentPointSitesError>(
            [&coincident, &domain]()
            {
                exact_point_site_graph(coincident, domain, 0);
            })
        && throws_named<EmptyPointSiteIdentityError>(
            [&points, &domain]()
            {
                std::vector<NormalizedPointSource2> empty_identity =
                    points;
                empty_identity.front().stable_site_id.clear();
                exact_point_site_graph(
                    empty_identity,
                    domain,
                    0);
            })
        && throws_named<InsufficientPointSitesError>(
            [&points, &domain]()
            {
                exact_point_site_graph(
                    {points.front()},
                    domain,
                    0);
            })
        && throws_named<NegativeSquaredRadiusError>(
            [&points, &domain]()
            {
                exact_point_site_graph(points, domain, -1);
            })
        && throws_named<DegeneratePointSiteTopologyError>(
            [&cocircular, &domain]()
            {
                exact_point_site_graph(cocircular, domain, 0);
            });
}
