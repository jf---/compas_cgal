#include "segment_site_mat.h"

#include <algorithm>
#include <cstdlib>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

bool segment_limiter_gate();

namespace {

bool has_reconstructible_root_id(
    const MatParameterEndpoint2& endpoint)
{
    if (!endpoint.parameter.has_value()) {
        return false;
    }
    const std::string root_id =
        algebraic_root_identity_v1(*endpoint.parameter);
    return std::find(
               endpoint.provenance_ids.begin(),
               endpoint.provenance_ids.end(),
               root_id)
        != endpoint.provenance_ids.end();
}

bool has_local_solution_identity(
    const MatParameterEndpoint2& endpoint)
{
    return std::any_of(
        endpoint.provenance_ids.begin(),
        endpoint.provenance_ids.end(),
        [](const std::string& provenance) {
            return provenance.find("solution-")
                != std::string::npos;
        });
}

bool parameter_domains_are_exact()
{
    using Point = MatTraits::Point_2;
    const MatTraits::Line_2 line(Point(0, 0), Point(1, 0));
    const MatTraits::Ray_2 ray(Point(0, 0), Point(1, 0));
    const MatTraits::Segment_2 segment(
        Point(0, 0),
        Point(1, 0));
    const SegmentSiteParabola2 parabola(
        Point(0, 1),
        MatTraits::Line_2(0, 1, 1),
        Point(-2, 1),
        Point(2, 1));

    const MatParameterDomain2 line_domain =
        exact_parameter_domain(line);
    const MatParameterDomain2 ray_domain =
        exact_parameter_domain(ray);
    const MatParameterDomain2 segment_domain =
        exact_parameter_domain(segment);
    const MatParameterDomain2 parabola_domain =
        exact_parameter_domain(parabola);
    return !line_domain.lower.has_value()
        && !line_domain.upper.has_value()
        && ray_domain.lower == MatTraits::FT(0)
        && !ray_domain.upper.has_value()
        && segment_domain.lower == MatTraits::FT(0)
        && segment_domain.upper == MatTraits::FT(1)
        && parabola_domain.lower.has_value()
        && parabola_domain.upper.has_value()
        && CGAL::compare(
               *parabola_domain.lower,
               *parabola_domain.upper)
            != CGAL::LARGER;
}

bool clearance_boundaries_are_exact()
{
    const RationalPrimitiveParameterization2 line{
        {CORE::BigRat(0), CORE::BigRat(1)},
        {CORE::BigRat(0)},
        std::nullopt,
        std::nullopt,
    };
    const RationalPrimitiveParameterization2 ray{
        line.x_coefficients,
        line.y_coefficients,
        CORE::BigRat(0),
        std::nullopt,
    };
    const RationalPrimitiveParameterization2 segment{
        {CORE::BigRat(0), CORE::BigRat(2)},
        {CORE::BigRat(0)},
        CORE::BigRat(0),
        CORE::BigRat(1),
    };
    const RationalPrimitiveParameterization2 parabola{
        {CORE::BigRat(0), CORE::BigRat(1)},
        {CORE::BigRat(0), CORE::BigRat(0), CORE::BigRat(1)},
        CORE::BigRat(-1),
        CORE::BigRat(1),
    };

    const ClearanceRootBoundary2 line_boundary =
        point_clearance_boundary(line, 0, 0, 1);
    const ClearanceRootBoundary2 ray_boundary =
        point_clearance_boundary(ray, 0, 0, 1);
    const ClearanceRootBoundary2 segment_boundary =
        point_clearance_boundary(segment, 0, 0, 1);
    const ClearanceRootBoundary2 parabola_boundary =
        point_clearance_boundary(parabola, 0, 0, 1);
    const ClearanceRootBoundary2 zero =
        point_clearance_boundary(
            {{CORE::BigRat(0)}, {CORE::BigRat(0)}, 0, 1},
            0,
            0,
            0);
    const ClearanceRootBoundary2 positive =
        point_clearance_boundary(
            {{CORE::BigRat(1)}, {CORE::BigRat(0)}, 0, 1},
            0,
            0,
            0);
    const ClearanceRootBoundary2 negative =
        point_clearance_boundary(
            {{CORE::BigRat(0)}, {CORE::BigRat(0)}, 0, 1},
            0,
            0,
            1);
    const RationalPrimitiveParameterization2 clipped_line{
        line.x_coefficients,
        line.y_coefficients,
        CORE::BigRat(-2),
        CORE::BigRat(2),
    };
    const ClearanceRootBoundary2 clipped_line_boundary =
        point_clearance_boundary(clipped_line, 0, 0, 1);

    const std::vector<MatAdmissibleComponent2> line_components =
        maximal_clearance_components(
            "line-dual",
            line,
            line_boundary);
    const std::vector<MatAdmissibleComponent2> ray_components =
        maximal_clearance_components(
            "ray-dual",
            ray,
            ray_boundary);
    const std::vector<MatAdmissibleComponent2> segment_components =
        maximal_clearance_components(
            "segment-dual",
            segment,
            segment_boundary);
    const std::vector<MatAdmissibleComponent2> parabola_components =
        maximal_clearance_components(
            "parabola-dual",
            parabola,
            parabola_boundary);
    const std::vector<MatAdmissibleComponent2>
        clipped_line_components =
            maximal_clearance_components(
                "clipped-line-dual",
                clipped_line,
                clipped_line_boundary);
    const std::vector<MatAdmissibleComponent2> zero_components =
        maximal_clearance_components(
            "zero-dual",
            {{CORE::BigRat(0)}, {CORE::BigRat(0)}, 0, 1},
            zero);
    const std::vector<MatAdmissibleComponent2> positive_components =
        maximal_clearance_components(
            "positive-dual",
            {{CORE::BigRat(1)}, {CORE::BigRat(0)}, 0, 1},
            positive);
    const std::vector<MatAdmissibleComponent2> negative_components =
        maximal_clearance_components(
            "negative-dual",
            {{CORE::BigRat(0)}, {CORE::BigRat(0)}, 0, 1},
            negative);

    MatDomainPolygon2 outer;
    outer.push_back({-3, -2});
    outer.push_back({3, -2});
    outer.push_back({3, 2});
    outer.push_back({-3, 2});
    MatDomainPolygon2 hole;
    hole.push_back({-1, -1});
    hole.push_back({-1, 1});
    hole.push_back({1, 1});
    hole.push_back({1, -1});
    const std::vector<MatDomainPolygon2> holes{hole};
    const MatDomainPolygonWithHoles2 domain(
        outer,
        holes.begin(),
        holes.end());

    const std::vector<MatAdmissibleComponent2>
        clipped_domain_line =
            clip_linear_clearance_components(
                "D-line-dual",
                line,
                line_boundary,
                domain);
    const std::vector<MatAdmissibleComponent2>
        clipped_domain_ray =
            clip_linear_clearance_components(
                "D-ray-dual",
                ray,
                ray_boundary,
                domain);
    const RationalPrimitiveParameterization2
        crossing_segment{
            {CORE::BigRat(-4), CORE::BigRat(8)},
            {CORE::BigRat(0)},
            CORE::BigRat(0),
            CORE::BigRat(1),
        };
    const ClearanceRootBoundary2
        crossing_segment_boundary =
            point_clearance_boundary(
                crossing_segment,
                0,
                10,
                0);
    const std::vector<MatAdmissibleComponent2>
        clipped_domain_segment =
            clip_linear_clearance_components(
                "D-segment-dual",
                crossing_segment,
                crossing_segment_boundary,
                domain);
    const RationalPrimitiveParameterization2
        crossing_parabola{
            {CORE::BigRat(0), CORE::BigRat(1)},
            {
                CORE::BigRat(0),
                CORE::BigRat(0),
                CORE::BigRat(1),
            },
            CORE::BigRat(-2),
            CORE::BigRat(2),
        };
    const ClearanceRootBoundary2
        crossing_parabola_boundary =
            point_clearance_boundary(
                crossing_parabola,
                0,
                0,
                0);
    const std::vector<MatAdmissibleComponent2>
        clipped_domain_parabola =
            clip_parabola_clearance_components(
                "D-parabola-dual",
                crossing_parabola,
                crossing_parabola_boundary,
                domain);
    MatDomainPolygon2 source_outer;
    source_outer.push_back({-4, -1});
    source_outer.push_back({4, -1});
    source_outer.push_back({4, 2});
    source_outer.push_back({-4, 2});
    const MatDomainPolygonWithHoles2 source_domain(
        source_outer,
        holes.begin(),
        holes.end());
    const std::vector<MatAdmissibleComponent2>
        source_parabola_components =
            clip_source_parabola_clearance_components(
                "source-parabola-dual",
                {
                    "irrational-focus",
                    {0, 1},
                    {1, 0},
                    2,
                },
                {
                    "open-directrix",
                    0,
                    1,
                    1,
                },
                CORE::BigRat(-4),
                CORE::BigRat(4),
                {
                    CGAL::POSITIVE,
                    {},
                    {},
                },
                source_domain);

    return line_boundary.roots.size() == 2
        && ray_boundary.roots.size() == 1
        && segment_boundary.roots.size() == 1
        && parabola_boundary.roots.size() == 2
        && zero.constant_sign == CGAL::ZERO
        && positive.constant_sign == CGAL::POSITIVE
        && negative.constant_sign == CGAL::NEGATIVE
        && line_components.size() == 2
        && !line_components.front().lower.parameter.has_value()
        && !line_components.back().upper.parameter.has_value()
        && ray_components.size() == 1
        && segment_components.size() == 1
        && parabola_components.size() == 2
        && clipped_line_components.size() == 2
        && clipped_line_components.front().lower.parameter.has_value()
        && clipped_line_components.back().upper.parameter.has_value()
        && clipped_line_components.front().component_id
            == "clipped-line-dual/component-0"
        && has_reconstructible_root_id(
            clipped_line_components.front().upper)
        && has_reconstructible_root_id(
            clipped_line_components.back().lower)
        && zero_components.size() == 1
        && positive_components.size() == 1
        && negative_components.empty()
        && clipped_domain_line.size() == 2
        && clipped_domain_line.front().lower.parameter.has_value()
        && clipped_domain_line.front().upper.parameter.has_value()
        && clipped_domain_line.back().lower.parameter.has_value()
        && clipped_domain_line.back().upper.parameter.has_value()
        && clipped_domain_ray.size() == 1
        && clipped_domain_segment.size() == 2
        && clipped_domain_segment.front()
            .lower.provenance_ids.front()
            == "D-segment-dual/D-outer/edge-3"
        && clipped_domain_segment.back()
            .upper.provenance_ids.front()
            == "D-segment-dual/D-outer/edge-1"
        && clipped_domain_parabola.size() == 2
        && has_reconstructible_root_id(
            clipped_domain_parabola.front().lower)
        && has_reconstructible_root_id(
            clipped_domain_parabola.back().upper)
        && source_parabola_components.size() == 2
        && has_reconstructible_root_id(
            source_parabola_components.front().lower)
        && has_reconstructible_root_id(
            source_parabola_components.back().upper)
        && !has_local_solution_identity(
            source_parabola_components.front().lower)
        && !has_local_solution_identity(
            source_parabola_components.back().upper);
}

bool live_graph_is_exact()
{
    const MatExactGraph2 graph =
        segment_site_live_graph_spike();
    const std::size_t lines = static_cast<std::size_t>(
        std::count_if(
            graph.edges.begin(),
            graph.edges.end(),
            [](const MatExactGraphEdge2& edge) {
                return edge.primitive_kind == "LINE";
            }));
    const std::size_t parabolas = static_cast<std::size_t>(
        std::count_if(
            graph.edges.begin(),
            graph.edges.end(),
            [](const MatExactGraphEdge2& edge) {
                return edge.primitive_kind == "PARABOLA";
            }));
    std::set<std::string> node_ids;
    std::set<std::string> edge_ids;
    for (const MatExactGraphNode2& node : graph.nodes) {
        node_ids.insert(node.node_id);
    }
    for (const MatExactGraphEdge2& edge : graph.edges) {
        edge_ids.insert(edge.edge_id);
    }
    return lines == 2
        && parabolas == 2
        && graph.rejected_incident_transitions > 0
        && node_ids.size() == graph.nodes.size()
        && edge_ids.size() == graph.edges.size()
        && std::all_of(
            graph.nodes.begin(),
            graph.nodes.end(),
            [](const MatExactGraphNode2& node) {
                return has_reconstructible_root_id(
                           node.endpoint)
                    && node.generator_site_ids.size() >= 2;
            })
        && std::all_of(
            graph.edges.begin(),
            graph.edges.end(),
            [](const MatExactGraphEdge2& edge) {
                return edge.generator_site_ids.size() == 2
                    && edge.source_node_id
                        != edge.target_node_id;
            });
}

const MatExactGraph2& generic_live_graph()
{
    static const MatExactGraph2 graph =
        segment_site_generic_graph_spike();
    return graph;
}

bool generic_graph_extracts_exact_linear_topology()
{
    const MatExactGraph2& graph = generic_live_graph();
    const bool has_ray = std::any_of(
        graph.edges.begin(),
        graph.edges.end(),
        [](const MatExactGraphEdge2& edge) {
            return edge.primitive_kind == "RAY";
        });
    const bool has_segment = std::any_of(
        graph.edges.begin(),
        graph.edges.end(),
        [](const MatExactGraphEdge2& edge) {
            return edge.primitive_kind == "SEGMENT";
        });
    return has_ray && has_segment;
}

bool generic_graph_unions_shared_voronoi_nodes()
{
    const MatExactGraph2& graph = generic_live_graph();
    return std::any_of(
        graph.nodes.begin(),
        graph.nodes.end(),
        [&graph](const MatExactGraphNode2& node) {
            const std::size_t incident_edges =
                static_cast<std::size_t>(std::count_if(
                    graph.edges.begin(),
                    graph.edges.end(),
                    [&node](const MatExactGraphEdge2& edge) {
                        return edge.source_node_id
                                == node.node_id
                            || edge.target_node_id
                                == node.node_id;
                    }));
            return incident_edges >= 2
                && node.generator_site_ids.size() >= 3;
        });
}

bool generic_graph_derives_open_segment_bounds()
{
    const MatExactGraph2& graph = generic_live_graph();
    return std::any_of(
        graph.edges.begin(),
        graph.edges.end(),
        [&graph](const MatExactGraphEdge2& edge) {
            if (edge.primitive_kind != "PARABOLA") {
                return false;
            }
            const auto is_node = [&edge](
                                     const MatExactGraphNode2& node) {
                return node.node_id == edge.source_node_id
                    || node.node_id == edge.target_node_id;
            };
            return std::count_if(
                       graph.nodes.begin(),
                       graph.nodes.end(),
                       [&is_node](const MatExactGraphNode2& node) {
                           return is_node(node)
                               && std::any_of(
                                   node.endpoint
                                       .provenance_ids.begin(),
                                   node.endpoint
                                       .provenance_ids.end(),
                                   [](const std::string& id) {
                                       return id
                                               == "segment-source"
                                           || id
                                               == "segment-target";
                                   });
                       })
                == 2;
        });
}

struct PointLimiterFixture2 {
    MatTraits::Point_2 segment_source;
    MatTraits::Point_2 segment_target;
    MatTraits::Point_2 focus;
    MatTraits::Point_2 limiter;
    MatExactPointSiteSource2 focus_record;
    MatExactOpenSegmentSource2 segment_record;
    MatExactPointSiteSource2 segment_source_record;
    MatExactPointSiteSource2 segment_target_record;
    MatExactPointSiteSource2 limiter_record;
};

SegmentSiteDelaunay2 point_limiter_delaunay(
    const PointLimiterFixture2& source)
{
    SegmentSiteDelaunay2 delaunay;
    delaunay.insert(
        source.segment_source,
        source.segment_target);
    delaunay.insert(source.focus);
    delaunay.insert(source.limiter);
    return delaunay;
}

struct LivePointLimiterFixture2 {
    explicit LivePointLimiterFixture2(
        PointLimiterFixture2 source)
        : source(std::move(source)),
          delaunay(point_limiter_delaunay(this->source)),
          voronoi(delaunay)
    {
    }

    PointLimiterFixture2 source;
    SegmentSiteDelaunay2 delaunay;
    SegmentSiteVoronoi2 voronoi;
};

struct BoundEndpointObservation2 {
    MatParameterEndpoint2 endpoint;
    SegmentSiteVoronoi2::Halfedge_handle halfedge;
    MatTraits::Point_2 point;
};

bool exact_site_equal(
    const MatTraits::Site_2& lhs,
    const MatTraits::Site_2& rhs)
{
    if (lhs.is_point() != rhs.is_point()) {
        return false;
    }
    if (lhs.is_point()) {
        return lhs.point() == rhs.point();
    }
    return (lhs.source() == rhs.source()
            && lhs.target() == rhs.target())
        || (lhs.source() == rhs.target()
            && lhs.target() == rhs.source());
}

std::vector<BoundEndpointObservation2>
bind_point_limiter_fixture(
    LivePointLimiterFixture2& fixture)
{
    const MatTraits::Site_2 focus_site =
        MatTraits::Site_2::construct_site_2(
            fixture.source.focus);
    const MatTraits::Site_2 segment_site =
        MatTraits::Site_2::construct_site_2(
            fixture.source.segment_source,
            fixture.source.segment_target);
    const MatTraits::Site_2 limiter_site =
        MatTraits::Site_2::construct_site_2(
            fixture.source.limiter);
    std::vector<BoundEndpointObservation2> result;
    for (auto halfedge = fixture.voronoi.halfedges_begin();
         halfedge != fixture.voronoi.halfedges_end();
         ++halfedge)
    {
        const MatTraits::Site_2 up = halfedge->up()->site();
        const MatTraits::Site_2 down = halfedge->down()->site();
        if (!((exact_site_equal(up, focus_site)
               && exact_site_equal(down, segment_site))
              || (exact_site_equal(up, segment_site)
                  && exact_site_equal(down, focus_site))))
        {
            continue;
        }
        const bool limiter_is_left =
            halfedge->has_source()
            && exact_site_equal(
                halfedge->left()->site(),
                limiter_site);
        const bool limiter_is_right =
            halfedge->has_target()
            && exact_site_equal(
                halfedge->right()->site(),
                limiter_site);
        if (limiter_is_left == limiter_is_right) {
            continue;
        }
        const MatTraits::Point_2 point =
            limiter_is_left
            ? halfedge->source()->point()
            : halfedge->target()->point();
        SegmentSiteParabola2 parabola;
        if (!CGAL::assign(
                parabola,
                fixture.voronoi.dual().primal(
                    halfedge->dual()))
            || (point != parabola.p1
                && point != parabola.p2))
        {
            continue;
        }
        result.push_back(
            {
                bind_point_limiter_parabola_endpoint(
                    fixture.source.focus_record,
                    fixture.source.segment_record,
                    fixture.source.segment_source_record,
                    fixture.source.segment_target_record,
                    fixture.source.limiter_record,
                    fixture.voronoi,
                    halfedge),
                halfedge,
                point,
            });
    }
    return result;
}

std::string canonical_endpoint_bytes(
    const MatParameterEndpoint2& endpoint)
{
    if (!endpoint.parameter.has_value()) {
        return "unbounded";
    }
    std::string result =
        algebraic_root_identity_v1(*endpoint.parameter);
    for (const std::string& provenance : endpoint.provenance_ids) {
        result += ":";
        result += std::to_string(provenance.size());
        result += ":";
        result += provenance;
    }
    return result;
}

PointLimiterFixture2 irrational_point_limiter_fixture()
{
    const CORE::Expr radical = CORE::sqrt(CORE::Expr(2));
    return {
        {-100, 0},
        {100, 0},
        {radical, 3},
        {0, CORE::BigRat(2, 3)},
        {
            "focus",
            {0, 1},
            {3, 0},
            2,
        },
        {
            "open-segment",
            0,
            1,
            0,
        },
        {
            "segment-source",
            {-100, 0},
            {0, 0},
            2,
        },
        {
            "segment-target",
            {100, 0},
            {0, 0},
            2,
        },
        {
            "limiter",
            {0, 0},
            {CORE::BigRat(2, 3), 0},
            2,
        },
    };
}

PointLimiterFixture2 rational_repeated_factor_fixture()
{
    return {
        {-100, 0},
        {100, 0},
        {1, 3},
        {0, 3},
        {
            "focus-rational",
            {1, 0},
            {3, 0},
            1,
        },
        {
            "open-segment-rational",
            0,
            1,
            0,
        },
        {
            "segment-source-rational",
            {-100, 0},
            {0, 0},
            1,
        },
        {
            "segment-target-rational",
            {100, 0},
            {0, 0},
            1,
        },
        {
            "limiter-rational",
            {0, 0},
            {3, 0},
            1,
        },
    };
}

bool point_limiter_binds_live_parabola_endpoint()
{
    LivePointLimiterFixture2 fixture(
        irrational_point_limiter_fixture());
    const std::vector<BoundEndpointObservation2> observations =
        bind_point_limiter_fixture(fixture);
    if (observations.size() != 2) {
        return false;
    }
    ExactAlgebraicKernel1 kernel;
    return observations[0].endpoint.parameter.has_value()
        && observations[1].endpoint.parameter.has_value()
        && kernel.compare_1_object()(
               *observations[0].endpoint.parameter,
               *observations[1].endpoint.parameter)
            == CGAL::EQUAL
        && observations[0].point == observations[1].point
        && observations[0].halfedge->opposite()
            == observations[1].halfedge
        && canonical_endpoint_bytes(observations[0].endpoint)
            == canonical_endpoint_bytes(observations[1].endpoint)
        && std::none_of(
               observations[0].endpoint.provenance_ids.begin(),
               observations[0].endpoint.provenance_ids.end(),
               [](const std::string& provenance) {
                   return provenance == "live-parabola/p1"
                       || provenance == "live-parabola/p2";
               });
}

bool rational_repeated_factor_binds_without_solver_precondition()
{
    LivePointLimiterFixture2 fixture(
        rational_repeated_factor_fixture());
    const std::vector<BoundEndpointObservation2> observations =
        bind_point_limiter_fixture(fixture);
    if (observations.size() != 2
        || !observations[0].endpoint.parameter.has_value()) {
        return false;
    }
    const std::vector<ExactAlgebraicInteger1>
        governing_factor{-1, 2};
    const std::string expected_root_id =
        algebraic_root_id_v1(governing_factor, 0);
    ExactAlgebraicKernel1 kernel;
    return kernel.compare_1_object()(
               *observations[0].endpoint.parameter,
               CORE::BigRat(1, 2))
            == CGAL::EQUAL
        && std::find(
               observations[0].endpoint.provenance_ids.begin(),
               observations[0].endpoint.provenance_ids.end(),
               expected_root_id)
            != observations[0].endpoint.provenance_ids.end()
        && std::find(
               observations[0].endpoint.provenance_ids.begin(),
               observations[0].endpoint.provenance_ids.end(),
               "point-limiter/norm-factor-multiplicity/2")
            != observations[0].endpoint.provenance_ids.end()
        && canonical_endpoint_bytes(observations[0].endpoint)
            == canonical_endpoint_bytes(observations[1].endpoint);
}

bool different_quadratic_field_is_rejected()
{
    LivePointLimiterFixture2 fixture(
        irrational_point_limiter_fixture());
    const std::vector<BoundEndpointObservation2> observations =
        bind_point_limiter_fixture(fixture);
    if (observations.empty()) {
        return false;
    }
    fixture.source.limiter_record.radicand = 3;
    try {
        bind_point_limiter_parabola_endpoint(
            fixture.source.focus_record,
            fixture.source.segment_record,
            fixture.source.segment_source_record,
            fixture.source.segment_target_record,
            fixture.source.limiter_record,
            fixture.voronoi,
            observations[0].halfedge);
    }
    catch (const InvalidRationalPrimitiveError&) {
        return true;
    }
    return false;
}

bool unrelated_live_generator_is_rejected()
{
    LivePointLimiterFixture2 fixture(
        irrational_point_limiter_fixture());
    const std::vector<BoundEndpointObservation2> observations =
        bind_point_limiter_fixture(fixture);
    if (observations.empty()) {
        return false;
    }
    MatExactPointSiteSource2 unrelated_focus =
        fixture.source.focus_record;
    unrelated_focus.x = {99, 0};
    unrelated_focus.y = {99, 0};
    try {
        bind_point_limiter_parabola_endpoint(
            unrelated_focus,
            fixture.source.segment_record,
            fixture.source.segment_source_record,
            fixture.source.segment_target_record,
            fixture.source.limiter_record,
            fixture.voronoi,
            observations[0].halfedge);
    }
    catch (const MismatchedLiveParabolaBridgeError&) {
        return true;
    }
    return false;
}

bool strict_open_segment_feature_rejects_decoy()
{
    const PointLimiterFixture2 fixture =
        rational_repeated_factor_fixture();
    const MatExactPointSiteSource2 decoy_source{
        "decoy-source",
        {2, 0},
        {0, 0},
        1,
    };
    const MatExactPointSiteSource2 decoy_target{
        "decoy-target",
        {3, 0},
        {0, 0},
        1,
    };
    const SourceParabolaParameterization2 parabola =
        source_parameterization(
            fixture.focus_record,
            fixture.segment_record);
    ExactAlgebraicKernel1 kernel;
    const auto parameter =
        kernel.construct_algebraic_real_1_object()(
            CORE::BigRat(1, 2));
    return exact_open_segment_feature_contains(
               parabola,
               fixture.segment_source_record,
               fixture.segment_target_record,
               parameter)
        && !exact_open_segment_feature_contains(
               parabola,
               decoy_source,
               decoy_target,
               parameter);
}

bool collapsed_live_parabola_is_rejected()
{
    SegmentSiteParabola2 collapsed(
        MatTraits::Point_2(0, 1),
        MatTraits::Line_2(0, 1, 0),
        MatTraits::Point_2(0, 1),
        MatTraits::Point_2(0, 1));
    try {
        require_distinct_live_parabola_endpoints(
            collapsed);
    }
    catch (const AmbiguousLiveParabolaEndpointError&) {
        return true;
    }
    return false;
}

} // namespace

int main()
{
    const SegmentSiteMatCompileEvidence2 evidence =
        segment_site_mat_compile_spike();
    return evidence.delaunay_valid
            && evidence.assigned_dual_primitives > 0
            && evidence.exact_clearance_roots == 2
            && evidence.matched_generator_sites
                == evidence.delaunay_vertices
            && parameter_domains_are_exact()
            && clearance_boundaries_are_exact()
            && live_graph_is_exact()
            && generic_graph_extracts_exact_linear_topology()
            && generic_graph_unions_shared_voronoi_nodes()
            && generic_graph_derives_open_segment_bounds()
            && point_limiter_binds_live_parabola_endpoint()
            && rational_repeated_factor_binds_without_solver_precondition()
            && different_quadratic_field_is_rejected()
            && unrelated_live_generator_is_rejected()
            && strict_open_segment_feature_rejects_decoy()
            && collapsed_live_parabola_is_rejected()
            && segment_limiter_gate()
        ? EXIT_SUCCESS
        : EXIT_FAILURE;
}
