#include "segment_site_mat.h"

#include <algorithm>
#include <iterator>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <CGAL/Object.h>
#include <CGAL/Polynomial_traits_d.h>

namespace {

std::size_t exact_clearance_root_count()
{
    using Polynomial =
        ExactAlgebraicKernel1::Polynomial_1;
    const std::vector<ExactAlgebraicInteger1>
        coefficients{-1, 0, 1};
    const Polynomial clearance =
        typename CGAL::Polynomial_traits_d<
            Polynomial>::Construct_polynomial()(
            coefficients.begin(),
            coefficients.end());
    std::vector<
        ExactAlgebraicKernel1::Algebraic_real_1>
        roots;
    ExactAlgebraicKernel1 kernel;
    kernel.solve_1_object()(
        clearance,
        true,
        std::back_inserter(roots));
    return roots.size();
}

std::size_t assign_dual_primitive(
    const CGAL::Object& dual)
{
    MatTraits::Line_2 line;
    MatTraits::Ray_2 ray;
    MatTraits::Segment_2 segment;
    SegmentSiteParabola2 parabola;
    if (CGAL::assign(line, dual)
        || CGAL::assign(ray, dual)
        || CGAL::assign(segment, dual)) {
        return 1;
    }
    if (!CGAL::assign(parabola, dual)) {
        return 0;
    }

    // `t` and `f` are the exact algebraic parameterization API. Drawing
    // helpers (`compute_k`, `generate_points`, streaming) are forbidden.
    const MatTraits::FT first_parameter =
        parabola.t(parabola.p1);
    const MatTraits::Point_2 reconstructed =
        parabola.f(first_parameter);
    return reconstructed == parabola.p1 ? 1 : 0;
}

std::vector<std::string> ordered_generator_ids(
    std::string first,
    std::string second)
{
    if (second < first) {
        std::swap(first, second);
    }
    return {std::move(first), std::move(second)};
}

std::string endpoint_node_id(
    const std::string& dual_id,
    const MatParameterEndpoint2& endpoint)
{
    if (!endpoint.parameter.has_value()) {
        throw InvalidRationalPrimitiveError(
            "emitted exact graph endpoint is unbounded");
    }
    return dual_id + "/node/"
        + algebraic_root_identity_v1(*endpoint.parameter);
}

void union_generator_ids(
    std::vector<std::string>& target,
    const std::vector<std::string>& source)
{
    target.insert(
        target.end(),
        source.begin(),
        source.end());
    std::sort(target.begin(), target.end());
    target.erase(
        std::unique(target.begin(), target.end()),
        target.end());
}

void append_components(
    const std::string& dual_id,
    const std::string& primitive_kind,
    const std::vector<std::string>& generator_ids,
    const std::vector<MatAdmissibleComponent2>& components,
    MatExactGraph2& graph,
    std::map<std::string, std::size_t>& node_indices)
{
    const auto append_node =
        [&graph, &node_indices, &dual_id, &generator_ids](
            const MatParameterEndpoint2& endpoint) {
            const std::string node_id =
                endpoint_node_id(dual_id, endpoint);
            MatParameterEndpoint2 bound_endpoint = endpoint;
            const std::string root_id =
                algebraic_root_identity_v1(
                    *endpoint.parameter);
            if (std::find(
                    bound_endpoint.provenance_ids.begin(),
                    bound_endpoint.provenance_ids.end(),
                    root_id)
                == bound_endpoint.provenance_ids.end()) {
                bound_endpoint.provenance_ids.push_back(
                    root_id);
            }
            const auto existing = node_indices.find(node_id);
            if (existing != node_indices.end()) {
                union_generator_ids(
                    graph.nodes[existing->second]
                        .generator_site_ids,
                    generator_ids);
                return node_id;
            }
            node_indices.emplace(
                node_id,
                graph.nodes.size());
            graph.nodes.push_back(
                {
                    node_id,
                    std::move(bound_endpoint),
                    generator_ids,
                });
            return node_id;
        };

    for (const MatAdmissibleComponent2& component : components) {
        const std::string source =
            append_node(component.lower);
        const std::string target =
            append_node(component.upper);
        graph.edges.push_back(
            {
                component.component_id,
                primitive_kind,
                source,
                target,
                generator_ids,
            });
    }
}

} // namespace

SegmentSiteMatCompileEvidence2
segment_site_mat_compile_spike()
{
    using Point = MatTraits::Point_2;
    using Site = MatTraits::Site_2;
    const Point lower_left(0, 0);
    const Point lower_right(4, 0);
    const Point upper_left(0, 3);
    const Point upper_right(4, 3);
    const Point isolated_point(6, 1);
    const std::vector<GeneratorSite2> generators{
        {"lower-left", Site::construct_site_2(lower_left)},
        {"lower-right", Site::construct_site_2(lower_right)},
        {"upper-left", Site::construct_site_2(upper_left)},
        {"upper-right", Site::construct_site_2(upper_right)},
        {"isolated", Site::construct_site_2(isolated_point)},
        {"lower-segment",
         Site::construct_site_2(lower_left, lower_right)},
        {"upper-segment",
         Site::construct_site_2(upper_left, upper_right)},
    };

    SegmentSiteDelaunay2 delaunay;
    delaunay.insert(lower_left, lower_right);
    delaunay.insert(upper_left, upper_right);
    delaunay.insert(isolated_point);

    SegmentSiteVoronoi2 voronoi(delaunay);
    static_cast<void>(voronoi);

    std::size_t assigned = 0;
    for (auto edge = delaunay.finite_edges_begin();
         edge != delaunay.finite_edges_end();
         ++edge) {
        assigned += assign_dual_primitive(
            delaunay.primal(edge));
    }
    const std::size_t delaunay_vertices = static_cast<std::size_t>(
        std::distance(
            delaunay.finite_vertices_begin(),
            delaunay.finite_vertices_end()));
    return {
        delaunay.is_valid(),
        assigned,
        exact_clearance_root_count(),
        matched_generator_site_count(delaunay, generators),
        delaunay_vertices,
    };
}

MatExactGraph2 segment_site_live_graph_spike()
{
    using Point = MatTraits::Point_2;
    using Site = MatTraits::Site_2;
    MatExactGraph2 graph{{}, {}, 0};
    std::map<std::string, std::size_t> node_indices;

    const Point line_left(-2, 0);
    const Point line_right(2, 0);
    const std::vector<GeneratorSite2> line_generators{
        {"line-left", Site::construct_site_2(line_left)},
        {"line-right", Site::construct_site_2(line_right)},
    };
    SegmentSiteDelaunay2 line_delaunay;
    line_delaunay.insert(line_left);
    line_delaunay.insert(line_right);
    for (auto edge = line_delaunay.finite_edges_begin();
         edge != line_delaunay.finite_edges_end();
         ++edge) {
        const auto value = *edge;
        const MatTraits::Site_2 first =
            value.first->vertex(
                CGAL::Triangulation_cw_ccw_2::ccw(
                    value.second))->site();
        const MatTraits::Site_2 second =
            value.first->vertex(
                CGAL::Triangulation_cw_ccw_2::cw(
                    value.second))->site();
        MatTraits::Line_2 live_line;
        if (!CGAL::assign(
                live_line,
                line_delaunay.primal(edge))) {
            continue;
        }
        const std::vector<std::string> generator_ids =
            ordered_generator_ids(
                stable_generator_site_id(
                    first,
                    line_generators),
                stable_generator_site_id(
                    second,
                    line_generators));
        MatDomainPolygon2 outer;
        outer.push_back({-1, -3});
        outer.push_back({1, -3});
        outer.push_back({1, 3});
        outer.push_back({-1, 3});
        const MatDomainPolygonWithHoles2 domain(outer);
        const RationalPrimitiveParameterization2 parameterization{
            {CORE::BigRat(0)},
            {CORE::BigRat(0), CORE::BigRat(1)},
            std::nullopt,
            std::nullopt,
        };
        const std::string dual_id =
            "live-dual/line-left/line-right";
        append_components(
            dual_id,
            "LINE",
            generator_ids,
            clip_linear_clearance_components(
                dual_id,
                parameterization,
                point_clearance_boundary(
                    parameterization,
                    -2,
                    0,
                    5),
                domain),
            graph,
            node_indices);
    }

    const Point segment_source(0, 0);
    const Point segment_target(4, 0);
    const Point focus(2, 3);
    const std::vector<GeneratorSite2> parabola_generators{
        {"segment-source",
         Site::construct_site_2(segment_source)},
        {"segment-target",
         Site::construct_site_2(segment_target)},
        {"open-segment",
         Site::construct_site_2(
             segment_source,
             segment_target)},
        {"focus", Site::construct_site_2(focus)},
    };
    SegmentSiteDelaunay2 parabola_delaunay;
    parabola_delaunay.insert(
        segment_source,
        segment_target);
    parabola_delaunay.insert(focus);
    for (auto edge = parabola_delaunay.finite_edges_begin();
         edge != parabola_delaunay.finite_edges_end();
         ++edge) {
        const auto value = *edge;
        const MatTraits::Site_2 first =
            value.first->vertex(
                CGAL::Triangulation_cw_ccw_2::ccw(
                    value.second))->site();
        const MatTraits::Site_2 second =
            value.first->vertex(
                CGAL::Triangulation_cw_ccw_2::cw(
                    value.second))->site();
        if (first.is_point() == second.is_point()) {
            continue;
        }
        const MatTraits::Site_2& point =
            first.is_point() ? first : second;
        const MatTraits::Site_2& segment =
            first.is_segment() ? first : second;
        if (point.point() == segment.source()
            || point.point() == segment.target()) {
            ++graph.rejected_incident_transitions;
            continue;
        }
        SegmentSiteParabola2 live_parabola;
        if (!CGAL::assign(
                live_parabola,
                parabola_delaunay.primal(edge))) {
            continue;
        }
        const std::vector<std::string> generator_ids =
            ordered_generator_ids(
                stable_generator_site_id(
                    point,
                    parabola_generators),
                stable_generator_site_id(
                    segment,
                    parabola_generators));
        MatDomainPolygon2 outer;
        outer.push_back({0, 1});
        outer.push_back({4, 1});
        outer.push_back({4, CORE::BigRat(13, 6)});
        outer.push_back({0, CORE::BigRat(13, 6)});
        MatDomainPolygon2 hole;
        hole.push_back(
            {CORE::BigRat(3, 2), 1});
        hole.push_back(
            {CORE::BigRat(3, 2), 2});
        hole.push_back(
            {CORE::BigRat(5, 2), 2});
        hole.push_back(
            {CORE::BigRat(5, 2), 1});
        const std::vector<MatDomainPolygon2> holes{hole};
        const MatDomainPolygonWithHoles2 domain(
            outer,
            holes.begin(),
            holes.end());
        const std::string dual_id =
            "live-dual/focus/open-segment";
        append_components(
            dual_id,
            "PARABOLA",
            generator_ids,
            clip_source_parabola_clearance_components(
                dual_id,
                {
                    "focus",
                    {2, 0},
                    {3, 0},
                    1,
                },
                {
                    "open-segment",
                    0,
                    1,
                    0,
                },
                CORE::BigRat(-2),
                CORE::BigRat(2),
                {CGAL::POSITIVE, {}, {}},
                domain),
            graph,
            node_indices);
    }
    return graph;
}

