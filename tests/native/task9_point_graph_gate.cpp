#include "segment_site_mat.h"

#include <algorithm>
#include <iostream>
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

bool require_contract(
    const bool condition,
    const char* contract)
{
    if (!condition)
    {
        std::cerr
            << "point graph contract failed: "
            << contract << '\n';
    }
    return condition;
}

void append_framed(
    std::string& target,
    const std::string& value)
{
    target += std::to_string(value.size())
        + ":" + value;
}

void append_count(
    std::string& target,
    const std::size_t value)
{
    append_framed(target, std::to_string(value));
}

void append_endpoint_signature(
    std::string& signature,
    const MatParameterEndpoint2& endpoint)
{
    if (endpoint.parameter.has_value())
    {
        append_framed(signature, "bounded");
        append_framed(
            signature,
            algebraic_root_identity_v1(
                *endpoint.parameter));
    }
    else
    {
        append_framed(signature, "unbounded");
    }
    append_count(
        signature,
        endpoint.provenance_ids.size());
    for (const std::string& provenance :
         endpoint.provenance_ids)
    {
        append_framed(signature, provenance);
    }
}

void append_node_provenance_signature(
    std::string& signature,
    const std::vector<std::string>& provenance_ids)
{
    append_count(
        signature,
        provenance_ids.size());
    for (const std::string& provenance :
         provenance_ids)
    {
        append_framed(signature, provenance);
    }
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

MatDomainPolygonWithHoles2 line_domain()
{
    MatDomainPolygon2 outer;
    outer.push_back({-3, -4});
    outer.push_back({3, -4});
    outer.push_back({3, 4});
    outer.push_back({-3, 4});
    return MatDomainPolygonWithHoles2(outer);
}

MatDomainPolygonWithHoles2 line_domain_with_hole()
{
    MatDomainPolygon2 outer;
    outer.push_back({-3, -4});
    outer.push_back({3, -4});
    outer.push_back({3, 4});
    outer.push_back({-3, 4});
    MatDomainPolygon2 hole;
    hole.push_back({-1, -1});
    hole.push_back({-1, 1});
    hole.push_back({1, 1});
    hole.push_back({1, -1});
    const std::vector<MatDomainPolygon2> holes{hole};
    return MatDomainPolygonWithHoles2(
        outer,
        holes.begin(),
        holes.end());
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
    append_framed(signature, "nodes");
    append_count(signature, graph.nodes.size());
    for (const MatExactGraphNode2& node : graph.nodes)
    {
        append_framed(signature, node.node_id);
        append_node_provenance_signature(
            signature,
            node.provenance_ids);
        append_count(
            signature,
            node.generator_site_ids.size());
        for (const std::string& site_id : node.generator_site_ids)
        {
            append_framed(signature, site_id);
        }
    }
    append_framed(signature, "edges");
    append_count(signature, graph.edges.size());
    for (const MatExactGraphEdge2& edge : graph.edges)
    {
        append_framed(signature, edge.edge_id);
        append_framed(signature, edge.primitive_kind);
        append_framed(signature, edge.source_node_id);
        append_framed(signature, edge.target_node_id);
        append_endpoint_signature(
            signature,
            edge.source_endpoint);
        append_endpoint_signature(
            signature,
            edge.target_endpoint);
        append_count(
            signature,
            edge.generator_site_ids.size());
        for (const std::string& site_id : edge.generator_site_ids)
        {
            append_framed(signature, site_id);
        }
    }
    append_framed(
        signature,
        "rejected-incident-transitions");
    append_count(
        signature,
        graph.rejected_incident_transitions);
    append_framed(
        signature,
        "matched-generator-sites");
    append_count(
        signature,
        graph.matched_generator_sites);
    return signature;
}

std::size_t primitive_count(
    const MatExactGraph2& graph,
    const std::string& primitive_kind)
{
    return static_cast<std::size_t>(std::count_if(
        graph.edges.begin(),
        graph.edges.end(),
        [&primitive_kind](
            const MatExactGraphEdge2& edge)
        {
            return edge.primitive_kind
                == primitive_kind;
        }));
}

bool exact_topology(
    const MatExactGraph2& graph,
    const std::size_t lines,
    const std::size_t rays,
    const std::size_t segments)
{
    const std::size_t actual_lines =
        primitive_count(graph, "LINE");
    const std::size_t actual_rays =
        primitive_count(graph, "RAY");
    const std::size_t actual_segments =
        primitive_count(graph, "SEGMENT");
    const bool matches =
        graph.edges.size() == lines + rays + segments
        && actual_lines == lines
        && actual_rays == rays
        && actual_segments == segments;
    if (!matches)
    {
        std::cerr
            << "topology expected L/R/S "
            << lines << "/" << rays << "/"
            << segments << ", observed "
            << actual_lines << "/" << actual_rays
            << "/" << actual_segments
            << ", total " << graph.edges.size()
            << '\n';
    }
    return matches;
}

bool has_hole_endpoint_provenance(
    const MatExactGraph2& graph)
{
    std::size_t count = 0;
    for (const MatExactGraphEdge2& edge : graph.edges)
    {
        for (const MatParameterEndpoint2* endpoint :
             {&edge.source_endpoint, &edge.target_endpoint})
        {
            if (std::any_of(
                    endpoint->provenance_ids.begin(),
                    endpoint->provenance_ids.end(),
                    [](const std::string& provenance)
                    {
                        return provenance.find(
                                   "/D-hole-0/")
                            != std::string::npos;
                    }))
            {
                ++count;
            }
        }
    }
    return count == 2;
}

bool edge_endpoint_is_bound(
    const MatParameterEndpoint2& endpoint)
{
    if (!endpoint.parameter.has_value())
    {
        return false;
    }
    const std::string root_id =
        algebraic_root_identity_v1(
            *endpoint.parameter);
    return std::find(
               endpoint.provenance_ids.begin(),
               endpoint.provenance_ids.end(),
               root_id)
        != endpoint.provenance_ids.end();
}

bool edge_bindings_are_exact(
    const MatExactGraph2& graph)
{
    return std::all_of(
        graph.edges.begin(),
        graph.edges.end(),
        [](const MatExactGraphEdge2& edge)
        {
            return edge_endpoint_is_bound(
                       edge.source_endpoint)
                && edge_endpoint_is_bound(
                       edge.target_endpoint);
        });
}

bool bijection_audit_rejects_unmatched_site()
{
    SegmentSiteDelaunay2 delaunay;
    const MatTraits::Point_2 left(-2, 0);
    const MatTraits::Point_2 right(2, 0);
    delaunay.insert(left);
    delaunay.insert(right);
    const std::vector<GeneratorSite2> incomplete{
        {
            "left",
            MatTraits::Site_2::construct_site_2(
                left),
        },
    };
    return throws_named<GeneratorSiteBijectionError>(
        [&delaunay, &incomplete]()
        {
            require_generator_site_bijection(
                delaunay,
                incomplete);
        });
}

bool dual_identity_is_injective()
{
    const MatDomainPolygonWithHoles2 domain =
        line_domain();
    const MatExactGraph2 first =
        exact_point_site_graph(
            {
                {"a", -2, 0},
                {"b/c", 2, 0},
            },
            domain,
            0);
    const MatExactGraph2 second =
        exact_point_site_graph(
            {
                {"a/b", -2, 0},
                {"c", 2, 0},
            },
            domain,
            0);
    if (first.edges.size() != 1
        || second.edges.size() != 1)
    {
        return false;
    }
    return first.edges.front().edge_id
            == "mat-dual/v1/5:point/1:a/3:b/c/component-0"
        && second.edges.front().edge_id
            == "mat-dual/v1/5:point/3:a/b/1:c/component-0"
        && first.edges.front().edge_id
            != second.edges.front().edge_id
        && first.edges.front().source_node_id
            != second.edges.front().source_node_id
        && first.edges.front().target_node_id
            != second.edges.front().target_node_id;
}

bool identity_helpers_reject_unordered_sites()
{
    return throws_named<InvalidDualIdentityError>(
               []()
               {
                   stable_dual_identity_v1(
                       "point",
                       {"right", "left"});
               })
        && throws_named<InvalidDualIdentityError>(
               []()
               {
                   stable_voronoi_node_identity_v1(
                       {"b", "a", "c"});
               });
}

bool analytic_two_site_line()
{
    return exact_topology(
        exact_point_site_graph(
            {
                {"left", -2, 0},
                {"right", 2, 0},
            },
            line_domain(),
            0),
        1,
        0,
        0);
}

bool analytic_three_site_hull()
{
    return exact_topology(
        exact_point_site_graph(
            {
                {"lower-left", -3, -2},
                {"lower-right", 3, -2},
                {"top", 0, 3},
            },
            production_domain(),
            0),
        0,
        3,
        0);
}

bool analytic_interior_site()
{
    return exact_topology(
        exact_point_site_graph(
            {
                {"lower-left", -4, -3},
                {"lower-right", 4, -3},
                {"top", 0, 5},
                {"interior", 0, 0},
            },
            production_domain(),
            0),
        0,
        3,
        3);
}

bool analytic_collinear_sites()
{
    return exact_topology(
        exact_point_site_graph(
            {
                {"left", -3, 0},
                {"middle", 0, 0},
                {"right", 4, 0},
            },
            line_domain(),
            0),
        2,
        0,
        0);
}

bool analytic_hole_splits_unbounded_line()
{
    const MatExactGraph2 graph =
        exact_point_site_graph(
            {
                {"left", -2, 0},
                {"right", 2, 0},
            },
            line_domain_with_hole(),
            0);
    return exact_topology(graph, 2, 0, 0)
        && graph.nodes.size() == 4
        && has_hole_endpoint_provenance(graph);
}

bool analytic_radius_clipping()
{
    const std::vector<NormalizedPointSource2> points{
        {"left", -2, 0},
        {"right", 2, 0},
    };
    const MatExactGraph2 split =
        exact_point_site_graph(
            points,
            line_domain(),
            5);
    const MatExactGraph2 tangent =
        exact_point_site_graph(
            points,
            line_domain(),
            4);
    return exact_topology(split, 2, 0, 0)
        && split.nodes.size() == 4
        && exact_topology(tangent, 1, 0, 0)
        && tangent.nodes.size() == 2;
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

    bool valid = true;
    valid = require_contract(
                !first.nodes.empty(),
                "production graph has nodes")
        && valid;
    valid = require_contract(
                !first.edges.empty(),
                "production graph has edges")
        && valid;
    valid = require_contract(
                graph_signature(first)
                    == graph_signature(second),
                "canonical signature is repeatable")
        && valid;
    valid = require_contract(
                graph_signature(first)
                    == graph_signature(reversed),
                "canonical signature ignores insertion order")
        && valid;
    valid = require_contract(
                first.matched_generator_sites
                        == points.size()
                    && second.matched_generator_sites
                        == points.size()
                    && reversed.matched_generator_sites
                        == points.size(),
                "live sites map bijectively to caller identities")
        && valid;
    valid = require_contract(
                edge_bindings_are_exact(first)
                    && edge_bindings_are_exact(second)
                    && edge_bindings_are_exact(reversed),
                "edge endpoints own exact bound parameters")
        && valid;
    valid = require_contract(
                bijection_audit_rejects_unmatched_site(),
                "bijection audit rejects unmatched live sites")
        && valid;
    valid = require_contract(
                excluded.nodes.empty()
                    && excluded.edges.empty(),
                "radius exclusion yields empty graph")
        && valid;
    valid = require_contract(
                dual_identity_is_injective(),
                "dual identity is injective")
        && valid;
    valid = require_contract(
                identity_helpers_reject_unordered_sites(),
                "derived identities require strict site order")
        && valid;
    valid = require_contract(
                analytic_two_site_line(),
                "two sites yield one line")
        && valid;
    valid = require_contract(
                analytic_three_site_hull(),
                "three hull sites yield three rays")
        && valid;
    valid = require_contract(
                analytic_interior_site(),
                "interior site yields three rays and three segments")
        && valid;
    valid = require_contract(
                analytic_collinear_sites(),
                "three collinear sites yield two lines")
        && valid;
    valid = require_contract(
                analytic_hole_splits_unbounded_line(),
                "hole splits an unbounded line")
        && valid;
    valid = require_contract(
                analytic_radius_clipping(),
                "exact radius clipping preserves tangent equality")
        && valid;
    valid = require_contract(
                throws_named<DuplicatePointSiteIdentityError>(
            [&duplicate_identity, &domain]()
            {
                exact_point_site_graph(
                    duplicate_identity,
                    domain,
                    0);
            }),
                "duplicate identities fail loudly")
        && valid;
    valid = require_contract(
                throws_named<CoincidentPointSitesError>(
            [&coincident, &domain]()
            {
                exact_point_site_graph(coincident, domain, 0);
            }),
                "coincident sites fail loudly")
        && valid;
    valid = require_contract(
                throws_named<EmptyPointSiteIdentityError>(
            [&points, &domain]()
            {
                std::vector<NormalizedPointSource2> empty_identity =
                    points;
                empty_identity.front().stable_site_id.clear();
                exact_point_site_graph(
                    empty_identity,
                    domain,
                    0);
            }),
                "empty identities fail loudly")
        && valid;
    valid = require_contract(
                throws_named<InsufficientPointSitesError>(
            [&points, &domain]()
            {
                exact_point_site_graph(
                    {points.front()},
                    domain,
                    0);
            }),
                "insufficient sites fail loudly")
        && valid;
    valid = require_contract(
                throws_named<NegativeSquaredRadiusError>(
            [&points, &domain]()
            {
                exact_point_site_graph(points, domain, -1);
            }),
                "negative squared radius fails loudly")
        && valid;
    valid = require_contract(
                throws_named<DegeneratePointSiteTopologyError>(
            [&cocircular, &domain]()
            {
                exact_point_site_graph(cocircular, domain, 0);
            }),
                "cocircular topology fails loudly")
        && valid;
    return valid;
}
