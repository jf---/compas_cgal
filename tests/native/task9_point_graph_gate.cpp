#include "segment_site_mat.h"

#include <algorithm>
#include <iostream>
#include <set>
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

struct ExactEdgeGolden2
{
    std::string edge_id;
    std::string primitive_kind;
    std::vector<std::string> generator_site_ids;
    CORE::BigRat source_parameter;
    CORE::BigRat target_parameter;
    std::string source_node_id;
    std::string target_node_id;
    std::vector<std::string> source_required_provenance_ids;
    std::vector<std::string> target_required_provenance_ids;
};

std::string rational_root_identity(
    const CORE::BigRat& parameter)
{
    ExactAlgebraicKernel1 kernel;
    return algebraic_root_identity_v1(
        kernel.construct_algebraic_real_1_object()(
            parameter));
}

std::string parameter_node_identity(
    const std::string& dual_id,
    const CORE::BigRat& parameter)
{
    return dual_id + "/node/"
        + rational_root_identity(parameter);
}

bool endpoint_matches_golden(
    const MatParameterEndpoint2& endpoint,
    const CORE::BigRat& parameter,
    const std::vector<std::string>&
        required_provenance_ids)
{
    if (!endpoint.parameter.has_value())
    {
        return false;
    }
    ExactAlgebraicKernel1 kernel;
    const auto expected =
        kernel.construct_algebraic_real_1_object()(
            parameter);
    const std::string root_id =
        rational_root_identity(parameter);
    const bool has_required_provenance =
        std::all_of(
            required_provenance_ids.begin(),
            required_provenance_ids.end(),
            [&endpoint](const std::string& provenance)
            {
                return std::find(
                           endpoint.provenance_ids.begin(),
                           endpoint.provenance_ids.end(),
                           provenance)
                    != endpoint.provenance_ids.end();
            });
    return kernel.compare_1_object()(
               *endpoint.parameter,
               expected)
            == CGAL::EQUAL
        && algebraic_root_identity_v1(
               *endpoint.parameter)
            == root_id
        && std::find(
               endpoint.provenance_ids.begin(),
               endpoint.provenance_ids.end(),
               root_id)
            != endpoint.provenance_ids.end()
        && has_required_provenance;
}

bool exact_edges_match_golden(
    const MatExactGraph2& graph,
    const std::size_t matched_generator_sites,
    const std::vector<ExactEdgeGolden2>& expected)
{
    bool matches =
        graph.edges.size() == expected.size()
        && graph.matched_generator_sites
            == matched_generator_sites
        && graph.rejected_incident_transitions == 0;
    std::set<std::string> expected_node_ids;
    for (const ExactEdgeGolden2& golden : expected)
    {
        expected_node_ids.insert(golden.source_node_id);
        expected_node_ids.insert(golden.target_node_id);
        const auto edge = std::find_if(
            graph.edges.begin(),
            graph.edges.end(),
            [&golden](const MatExactGraphEdge2& candidate)
            {
                return candidate.edge_id
                    == golden.edge_id;
            });
        if (edge == graph.edges.end())
        {
            std::cerr
                << "missing exact edge "
                << golden.edge_id << '\n';
            matches = false;
            continue;
        }
        const bool edge_matches =
            edge->primitive_kind
                == golden.primitive_kind
            && edge->generator_site_ids
                == golden.generator_site_ids
            && edge->source_node_id
                == golden.source_node_id
            && edge->target_node_id
                == golden.target_node_id
            && endpoint_matches_golden(
                edge->source_endpoint,
                golden.source_parameter,
                golden.source_required_provenance_ids)
            && endpoint_matches_golden(
                edge->target_endpoint,
                golden.target_parameter,
                golden.target_required_provenance_ids);
        if (!edge_matches)
        {
            std::cerr
                << "exact edge mismatch "
                << golden.edge_id << '\n'
                << "  primitive: "
                << edge->primitive_kind << '\n'
                << "  source node: "
                << edge->source_node_id << '\n'
                << "  target node: "
                << edge->target_node_id << '\n';
            matches = false;
        }
    }
    std::set<std::string> actual_node_ids;
    for (const MatExactGraphNode2& node : graph.nodes)
    {
        actual_node_ids.insert(node.node_id);
    }
    if (actual_node_ids != expected_node_ids
        || actual_node_ids.size() != graph.nodes.size())
    {
        std::cerr
            << "exact node connectivity mismatch\n";
        matches = false;
    }
    return matches;
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
    const MatExactGraph2 graph =
        exact_point_site_graph(
            {
                {"left", -2, 0},
                {"right", 2, 0},
            },
            line_domain(),
            0);
    const std::string dual_id =
        "mat-dual/v1/5:point/4:left/5:right";
    return exact_edges_match_golden(
        graph,
        2,
        {
            {
                dual_id + "/component-0",
                "LINE",
                {"left", "right"},
                CORE::BigRat(-1),
                CORE::BigRat(1),
                parameter_node_identity(
                    dual_id,
                    CORE::BigRat(-1)),
                parameter_node_identity(
                    dual_id,
                    CORE::BigRat(1)),
            },
        });
}

bool analytic_three_site_hull()
{
    const MatExactGraph2 graph =
        exact_point_site_graph(
            {
                {"lower-left", -3, -2},
                {"lower-right", 3, -2},
                {"top", 0, 3},
            },
            production_domain(),
            0);
    const std::string lower_dual =
        "mat-dual/v1/5:point/10:lower-left/11:lower-right";
    const std::string left_top_dual =
        "mat-dual/v1/5:point/10:lower-left/3:top";
    const std::string right_top_dual =
        "mat-dual/v1/5:point/11:lower-right/3:top";
    const std::string circumcenter =
        "voronoi-node/v1/10:lower-left/11:lower-right/3:top";
    return exact_edges_match_golden(
        graph,
        3,
        {
            {
                lower_dual + "/component-0",
                "RAY",
                {"lower-left", "lower-right"},
                CORE::BigRat(-1),
                CORE::BigRat(4, 15),
                parameter_node_identity(
                    lower_dual,
                    CORE::BigRat(-1)),
                circumcenter,
            },
            {
                left_top_dual + "/component-0",
                "RAY",
                {"lower-left", "top"},
                CORE::BigRat(-3, 10),
                CORE::BigRat(13, 10),
                circumcenter,
                parameter_node_identity(
                    left_top_dual,
                    CORE::BigRat(13, 10)),
            },
            {
                right_top_dual + "/component-0",
                "RAY",
                {"lower-right", "top"},
                CORE::BigRat(-13, 10),
                CORE::BigRat(3, 10),
                parameter_node_identity(
                    right_top_dual,
                    CORE::BigRat(-13, 10)),
                circumcenter,
            },
        });
}

bool analytic_interior_site()
{
    const MatExactGraph2 graph =
        exact_point_site_graph(
            {
                {"lower-left", -4, -3},
                {"lower-right", 4, -3},
                {"top", 0, 5},
                {"interior", 0, 0},
            },
            production_domain(),
            0);
    const std::string lower_dual =
        "mat-dual/v1/5:point/10:lower-left/11:lower-right";
    const std::string left_top_dual =
        "mat-dual/v1/5:point/10:lower-left/3:top";
    const std::string right_top_dual =
        "mat-dual/v1/5:point/11:lower-right/3:top";
    const std::string interior_left_dual =
        "mat-dual/v1/5:point/8:interior/10:lower-left";
    const std::string interior_right_dual =
        "mat-dual/v1/5:point/8:interior/11:lower-right";
    const std::string interior_top_dual =
        "mat-dual/v1/5:point/8:interior/3:top";
    const std::string lower_node =
        "voronoi-node/v1/8:interior/10:lower-left/11:lower-right";
    const std::string left_node =
        "voronoi-node/v1/8:interior/10:lower-left/3:top";
    const std::string right_node =
        "voronoi-node/v1/8:interior/11:lower-right/3:top";
    return exact_edges_match_golden(
        graph,
        4,
        {
            {
                lower_dual + "/component-0",
                "RAY",
                {"lower-left", "lower-right"},
                CORE::BigRat(-5, 8),
                CORE::BigRat(-7, 48),
                parameter_node_identity(
                    lower_dual,
                    CORE::BigRat(-5, 8)),
                lower_node,
            },
            {
                left_top_dual + "/component-0",
                "RAY",
                {"lower-left", "top"},
                CORE::BigRat(3, 8),
                CORE::BigRat(3, 4),
                left_node,
                parameter_node_identity(
                    left_top_dual,
                    CORE::BigRat(3, 4)),
            },
            {
                right_top_dual + "/component-0",
                "RAY",
                {"lower-right", "top"},
                CORE::BigRat(-3, 4),
                CORE::BigRat(-3, 8),
                parameter_node_identity(
                    right_top_dual,
                    CORE::BigRat(-3, 4)),
                right_node,
            },
            {
                interior_left_dual + "/component-0",
                "SEGMENT",
                {"interior", "lower-left"},
                CORE::BigRat(-1),
                CORE::BigRat(2, 3),
                left_node,
                lower_node,
            },
            {
                interior_right_dual + "/component-0",
                "SEGMENT",
                {"interior", "lower-right"},
                CORE::BigRat(-2, 3),
                CORE::BigRat(1),
                lower_node,
                right_node,
            },
            {
                interior_top_dual + "/component-0",
                "SEGMENT",
                {"interior", "top"},
                CORE::BigRat(-1),
                CORE::BigRat(1),
                right_node,
                left_node,
            },
        });
}

bool analytic_collinear_sites()
{
    const MatExactGraph2 graph =
        exact_point_site_graph(
            {
                {"left", -3, 0},
                {"middle", 0, 0},
                {"right", 4, 0},
            },
            line_domain(),
            0);
    const std::string left_dual =
        "mat-dual/v1/5:point/4:left/6:middle";
    const std::string right_dual =
        "mat-dual/v1/5:point/6:middle/5:right";
    return exact_edges_match_golden(
        graph,
        3,
        {
            {
                left_dual + "/component-0",
                "LINE",
                {"left", "middle"},
                CORE::BigRat(-4, 3),
                CORE::BigRat(4, 3),
                parameter_node_identity(
                    left_dual,
                    CORE::BigRat(-4, 3)),
                parameter_node_identity(
                    left_dual,
                    CORE::BigRat(4, 3)),
            },
            {
                right_dual + "/component-0",
                "LINE",
                {"middle", "right"},
                CORE::BigRat(-1),
                CORE::BigRat(1),
                parameter_node_identity(
                    right_dual,
                    CORE::BigRat(-1)),
                parameter_node_identity(
                    right_dual,
                    CORE::BigRat(1)),
            },
        });
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
    const std::string dual_id =
        "mat-dual/v1/5:point/4:left/5:right";
    return exact_edges_match_golden(
        graph,
        2,
        {
            {
                dual_id + "/component-0",
                "LINE",
                {"left", "right"},
                CORE::BigRat(-1),
                CORE::BigRat(-1, 4),
                parameter_node_identity(
                    dual_id,
                    CORE::BigRat(-1)),
                parameter_node_identity(
                    dual_id,
                    CORE::BigRat(-1, 4)),
                {dual_id + "/D-outer/edge-0"},
                {dual_id + "/D-hole-0/edge-3"},
            },
            {
                dual_id + "/component-1",
                "LINE",
                {"left", "right"},
                CORE::BigRat(1, 4),
                CORE::BigRat(1),
                parameter_node_identity(
                    dual_id,
                    CORE::BigRat(1, 4)),
                parameter_node_identity(
                    dual_id,
                    CORE::BigRat(1)),
                {dual_id + "/D-hole-0/edge-1"},
                {dual_id + "/D-outer/edge-2"},
            },
        });
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
    const std::string dual_id =
        "mat-dual/v1/5:point/4:left/5:right";
    return exact_edges_match_golden(
               split,
               2,
               {
                   {
                       dual_id + "/component-0",
                       "LINE",
                       {"left", "right"},
                       CORE::BigRat(-1),
                       CORE::BigRat(-1, 4),
                       parameter_node_identity(
                           dual_id,
                           CORE::BigRat(-1)),
                       parameter_node_identity(
                           dual_id,
                           CORE::BigRat(-1, 4)),
                       {dual_id + "/D-outer/edge-0"},
                       {},
                   },
                   {
                       dual_id + "/component-1",
                       "LINE",
                       {"left", "right"},
                       CORE::BigRat(1, 4),
                       CORE::BigRat(1),
                       parameter_node_identity(
                           dual_id,
                           CORE::BigRat(1, 4)),
                       parameter_node_identity(
                           dual_id,
                           CORE::BigRat(1)),
                       {},
                       {dual_id + "/D-outer/edge-2"},
                   },
               })
        && exact_edges_match_golden(
               tangent,
               2,
               {
                   {
                       dual_id + "/component-0",
                       "LINE",
                       {"left", "right"},
                       CORE::BigRat(-1),
                       CORE::BigRat(1),
                       parameter_node_identity(
                           dual_id,
                           CORE::BigRat(-1)),
                       parameter_node_identity(
                           dual_id,
                           CORE::BigRat(1)),
                       {dual_id + "/D-outer/edge-0"},
                       {dual_id + "/D-outer/edge-2"},
                   },
               });
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
                "two-site exact line golden")
        && valid;
    valid = require_contract(
                analytic_three_site_hull(),
                "three-site ray and circumcenter goldens")
        && valid;
    valid = require_contract(
                analytic_interior_site(),
                "interior-site exact connectivity golden")
        && valid;
    valid = require_contract(
                analytic_collinear_sites(),
                "collinear dimension-one exact goldens")
        && valid;
    valid = require_contract(
                analytic_hole_splits_unbounded_line(),
                "hole split roots and provenance goldens")
        && valid;
    valid = require_contract(
                analytic_radius_clipping(),
                "clearance split and tangent goldens")
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
