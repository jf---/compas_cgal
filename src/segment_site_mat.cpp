#include "segment_site_mat.h"
#include "segment_site_graph_emission.h"

#include <algorithm>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <CGAL/Object.h>
#include <CGAL/Kernel/global_functions_2.h>

namespace {

struct NormalizedOpenSegmentSource2 {
    std::string stable_site_id;
    std::string source_point_id;
    std::string target_point_id;
};

void require_valid_point_graph_inputs(
    const std::vector<NormalizedPointSource2>& points,
    const CORE::BigRat& radius_squared)
{
    if (points.size() < 2) {
        throw InsufficientPointSitesError(
            "point-site graph requires at least two sites");
    }
    if (radius_squared < 0) {
        throw NegativeSquaredRadiusError(
            "point-site graph squared radius is negative");
    }
    std::set<std::string> stable_ids;
    std::set<std::pair<CORE::BigRat, CORE::BigRat>> coordinates;
    for (const NormalizedPointSource2& point : points) {
        if (point.stable_site_id.empty()) {
            throw EmptyPointSiteIdentityError(
                "point-site graph identity is empty");
        }
        if (!stable_ids.insert(point.stable_site_id).second) {
            throw DuplicatePointSiteIdentityError(
                "point-site graph identity is duplicated");
        }
        if (!coordinates.emplace(point.x, point.y).second) {
            throw CoincidentPointSitesError(
                "point-site graph sites coincide");
        }
    }
}

const NormalizedPointSource2& point_source(
    const std::vector<NormalizedPointSource2>& points,
    const std::string& stable_site_id)
{
    const auto found = std::find_if(
        points.begin(),
        points.end(),
        [&stable_site_id](
            const NormalizedPointSource2& point) {
            return point.stable_site_id == stable_site_id;
        });
    if (found == points.end()) {
        throw InvalidRationalPrimitiveError(
            "point site has no normalized source record");
    }
    return *found;
}

std::vector<GeneratorSite2> point_generators(
    const std::vector<NormalizedPointSource2>& points)
{
    std::vector<GeneratorSite2> generators;
    generators.reserve(points.size());
    for (const NormalizedPointSource2& point : points) {
        if (point.stable_site_id.empty()) {
            throw InvalidRationalPrimitiveError(
                "normalized point identity is empty");
        }
        generators.push_back(
            {
                point.stable_site_id,
                MatTraits::Site_2::construct_site_2(
                    MatTraits::Point_2(point.x, point.y)),
            });
    }
    return generators;
}

RationalPrimitiveParameterization2 point_bisector(
    const NormalizedPointSource2& first,
    const NormalizedPointSource2& second,
    std::optional<CORE::BigRat> lower,
    std::optional<CORE::BigRat> upper)
{
    const CORE::BigRat direction_x =
        first.y - second.y;
    const CORE::BigRat direction_y =
        second.x - first.x;
    if (direction_x == 0 && direction_y == 0) {
        throw InvalidRationalPrimitiveError(
            "point bisector sources coincide");
    }
    return {
        {
            (first.x + second.x) / 2,
            direction_x,
        },
        {
            (first.y + second.y) / 2,
            direction_y,
        },
        std::move(lower),
        std::move(upper),
    };
}

CORE::BigRat parameter_at(
    const RationalPrimitiveParameterization2& primitive,
    const CORE::BigRat& x,
    const CORE::BigRat& y)
{
    if (primitive.x_coefficients[1] != 0) {
        return (x - primitive.x_coefficients[0])
            / primitive.x_coefficients[1];
    }
    return (y - primitive.y_coefficients[0])
        / primitive.y_coefficients[1];
}

struct RationalVoronoiNode2 {
    CORE::BigRat x;
    CORE::BigRat y;
    std::string node_id;
    std::vector<std::string> generator_site_ids;
};

RationalVoronoiNode2 rational_voronoi_node(
    const SegmentSiteDelaunay2::Face_handle& face,
    const std::vector<GeneratorSite2>& generators,
    const std::vector<NormalizedPointSource2>& points)
{
    std::vector<std::string> ids;
    ids.reserve(3);
    for (int index = 0; index < 3; ++index) {
        const MatTraits::Site_2 site =
            face->vertex(index)->site();
        if (!site.is_point()) {
            throw InvalidRationalPrimitiveError(
                "linear Voronoi node has nonpoint source");
        }
        ids.push_back(
            stable_generator_site_id(site, generators));
    }
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    if (ids.size() != 3) {
        throw InvalidRationalPrimitiveError(
            "Voronoi node source triple is degenerate");
    }

    const NormalizedPointSource2& a =
        point_source(points, ids[0]);
    const NormalizedPointSource2& b =
        point_source(points, ids[1]);
    const NormalizedPointSource2& c =
        point_source(points, ids[2]);
    const MatDomainKernel2::Point_2 point_a(a.x, a.y);
    const MatDomainKernel2::Point_2 point_b(b.x, b.y);
    const MatDomainKernel2::Point_2 point_c(c.x, c.y);
    if (CGAL::collinear(point_a, point_b, point_c)) {
        throw InvalidRationalPrimitiveError(
            "Voronoi source triple is collinear");
    }
    const MatDomainKernel2::Point_2 circumcenter =
        CGAL::circumcenter(point_a, point_b, point_c);
    return {
        circumcenter.x(),
        circumcenter.y(),
        stable_voronoi_node_identity_v1(ids),
        std::move(ids),
    };
}

void register_voronoi_node_location(
    const RationalVoronoiNode2& node,
    std::map<
        std::pair<CORE::BigRat, CORE::BigRat>,
        std::string>& node_ids_by_location)
{
    const auto [existing, inserted] =
        node_ids_by_location.emplace(
            std::make_pair(node.x, node.y),
            node.node_id);
    if (!inserted && existing->second != node.node_id) {
        throw DegeneratePointSiteTopologyError(
            "distinct point-site triples share an exact Voronoi node");
    }
}

struct CanonicalNodeAlias2 {
    std::string node_id;
    std::vector<std::string> generator_site_ids;
};

void register_node_alias(
    const std::string& dual_id,
    const CORE::BigRat& parameter,
    const RationalVoronoiNode2& node,
    std::map<std::string, CanonicalNodeAlias2>& aliases)
{
    ExactAlgebraicKernel1 kernel;
    const MatParameterEndpoint2 endpoint{
        kernel.construct_algebraic_real_1_object()(parameter),
        {},
    };
    aliases.emplace(
        stable_endpoint_node_identity_v1(
            dual_id,
            endpoint),
        CanonicalNodeAlias2{
            node.node_id,
            node.generator_site_ids,
        });
}

void canonicalize_original_nodes(
    MatExactGraph2& graph,
    const std::map<std::string, CanonicalNodeAlias2>& aliases)
{
    for (MatExactGraphEdge2& edge : graph.edges) {
        const auto source = aliases.find(edge.source_node_id);
        if (source != aliases.end()) {
            edge.source_node_id = source->second.node_id;
        }
        const auto target = aliases.find(edge.target_node_id);
        if (target != aliases.end()) {
            edge.target_node_id = target->second.node_id;
        }
    }

    std::vector<MatExactGraphNode2> canonical;
    std::map<std::string, std::size_t> indices;
    for (MatExactGraphNode2& node : graph.nodes) {
        const auto alias = aliases.find(node.node_id);
        if (alias != aliases.end()) {
            node.node_id = alias->second.node_id;
            union_stable_ids(
                node.generator_site_ids,
                alias->second.generator_site_ids);
        }
        const auto existing = indices.find(node.node_id);
        if (existing == indices.end()) {
            indices.emplace(node.node_id, canonical.size());
            canonical.push_back(std::move(node));
            continue;
        }
        MatExactGraphNode2& retained =
            canonical[existing->second];
        union_stable_ids(
            retained.generator_site_ids,
            node.generator_site_ids);
        union_stable_ids(
            retained.provenance_ids,
            node.provenance_ids);
    }
    graph.nodes = std::move(canonical);
}

void append_dimension_one_point_site_graph(
    const SegmentSiteDelaunay2& delaunay,
    const std::vector<GeneratorSite2>& generators,
    const std::vector<NormalizedPointSource2>& points,
    const MatDomainPolygonWithHoles2& domain,
    const CORE::BigRat& radius_squared,
    MatExactGraph2& graph,
    std::map<std::string, std::size_t>& node_indices)
{
    for (auto edge = delaunay.finite_edges_begin();
         edge != delaunay.finite_edges_end();
         ++edge) {
        const auto value = *edge;
        const std::vector<std::string> generator_ids =
            ordered_generator_site_ids(
                stable_generator_site_id(
                    value.first->vertex(
                        CGAL::Triangulation_cw_ccw_2::ccw(
                            value.second))->site(),
                    generators),
                stable_generator_site_id(
                    value.first->vertex(
                        CGAL::Triangulation_cw_ccw_2::cw(
                            value.second))->site(),
                    generators));
        MatTraits::Line_2 live_line;
        if (!CGAL::assign(
                live_line,
                delaunay.primal(edge))) {
            throw InvalidRationalPrimitiveError(
                "dimension-one point dual is not a line");
        }
        const NormalizedPointSource2& first =
            point_source(points, generator_ids[0]);
        const RationalPrimitiveParameterization2 primitive =
            point_bisector(
                first,
                point_source(points, generator_ids[1]),
                std::nullopt,
                std::nullopt);
        const std::string dual_id =
            stable_dual_identity_v1("point", generator_ids);
        append_exact_graph_components(
            dual_id,
            "LINE",
            generator_ids,
            clip_linear_clearance_components(
                dual_id,
                primitive,
                point_clearance_boundary(
                    primitive,
                    first.x,
                    first.y,
                    radius_squared),
                domain),
            graph,
            node_indices);
    }
}

void append_point_site_graph(
    const std::vector<NormalizedPointSource2>& points,
    const MatDomainPolygonWithHoles2& domain,
    const CORE::BigRat& radius_squared,
    MatExactGraph2& graph,
    std::map<std::string, std::size_t>& node_indices,
    std::map<std::string, CanonicalNodeAlias2>& aliases)
{
    const std::vector<GeneratorSite2> generators =
        point_generators(points);
    SegmentSiteDelaunay2 delaunay;
    for (const NormalizedPointSource2& point : points) {
        delaunay.insert(MatTraits::Point_2(point.x, point.y));
    }
    graph.matched_generator_sites =
        require_generator_site_bijection(delaunay, generators);
    if (delaunay.dimension() == 1) {
        append_dimension_one_point_site_graph(
            delaunay,
            generators,
            points,
            domain,
            radius_squared,
            graph,
            node_indices);
        return;
    }
    std::map<
        std::pair<CORE::BigRat, CORE::BigRat>,
        std::string> node_ids_by_location;

    for (auto edge = delaunay.finite_edges_begin();
         edge != delaunay.finite_edges_end();
         ++edge) {
        const auto value = *edge;
        const int opposite = value.second;
        const MatTraits::Site_2 first_site =
            value.first->vertex(
                CGAL::Triangulation_cw_ccw_2::ccw(
                    opposite))->site();
        const MatTraits::Site_2 second_site =
            value.first->vertex(
                CGAL::Triangulation_cw_ccw_2::cw(
                    opposite))->site();
        const std::vector<std::string> generator_ids =
            ordered_generator_site_ids(
                stable_generator_site_id(
                    first_site,
                    generators),
                stable_generator_site_id(
                    second_site,
                    generators));
        const NormalizedPointSource2& first =
            point_source(points, generator_ids[0]);
        const NormalizedPointSource2& second =
            point_source(points, generator_ids[1]);
        const std::string dual_id =
            stable_dual_identity_v1(
                "point",
                generator_ids);

        const auto neighbor =
            value.first->neighbor(opposite);
        const bool first_finite =
            !delaunay.is_infinite(value.first);
        const bool second_finite =
            !delaunay.is_infinite(neighbor);
        std::optional<RationalVoronoiNode2> first_node;
        std::optional<RationalVoronoiNode2> second_node;
        if (first_finite) {
            first_node = rational_voronoi_node(
                value.first,
                generators,
                points);
            register_voronoi_node_location(
                *first_node,
                node_ids_by_location);
        }
        if (second_finite) {
            second_node = rational_voronoi_node(
                neighbor,
                generators,
                points);
            register_voronoi_node_location(
                *second_node,
                node_ids_by_location);
        }

        const RationalPrimitiveParameterization2 unbounded =
            point_bisector(
                first,
                second,
                std::nullopt,
                std::nullopt);
        std::optional<CORE::BigRat> first_parameter;
        std::optional<CORE::BigRat> second_parameter;
        if (first_node.has_value()) {
            first_parameter = parameter_at(
                unbounded,
                first_node->x,
                first_node->y);
        }
        if (second_node.has_value()) {
            second_parameter = parameter_at(
                unbounded,
                second_node->x,
                second_node->y);
        }

        std::string primitive_kind;
        std::optional<CORE::BigRat> lower;
        std::optional<CORE::BigRat> upper;
        const CGAL::Object dual = delaunay.primal(edge);
        MatTraits::Segment_2 live_segment;
        MatTraits::Ray_2 live_ray;
        MatTraits::Line_2 live_line;
        if (CGAL::assign(live_segment, dual)) {
            if (!first_parameter.has_value()
                || !second_parameter.has_value()) {
                throw InvalidRationalPrimitiveError(
                    "live segment lacks two exact source nodes");
            }
            lower = std::min(
                *first_parameter,
                *second_parameter);
            upper = std::max(
                *first_parameter,
                *second_parameter);
            primitive_kind = "SEGMENT";
        } else if (CGAL::assign(live_ray, dual)) {
            const CORE::BigRat finite_parameter =
                first_parameter.has_value()
                ? *first_parameter
                : *second_parameter;
            const auto finite_face = first_finite
                ? value.first
                : neighbor;
            int finite_opposite = opposite;
            if (!first_finite) {
                finite_opposite = -1;
                for (int index = 0; index < 3; ++index) {
                    if (finite_face->neighbor(index)
                        == value.first) {
                        finite_opposite = index;
                        break;
                    }
                }
                if (finite_opposite < 0) {
                    throw InvalidRationalPrimitiveError(
                        "ray faces are not reciprocal");
                }
            }
            const std::string limiter_id =
                stable_generator_site_id(
                    finite_face->vertex(
                        finite_opposite)->site(),
                    generators);
            const NormalizedPointSource2& limiter =
                point_source(points, limiter_id);
            const CORE::BigRat direction_x =
                unbounded.x_coefficients[1];
            const CORE::BigRat direction_y =
                unbounded.y_coefficients[1];
            const CORE::BigRat derivative =
                2
                * (direction_x
                       * (limiter.x - first.x)
                   + direction_y
                       * (limiter.y - first.y));
            if (derivative == 0) {
                throw InvalidRationalPrimitiveError(
                    "ray limiter does not orient dual");
            }
            if (CGAL::sign(derivative)
                == CGAL::POSITIVE) {
                upper = finite_parameter;
            } else {
                lower = finite_parameter;
            }
            primitive_kind = "RAY";
        } else if (CGAL::assign(live_line, dual)) {
            primitive_kind = "LINE";
        } else {
            continue;
        }

        const RationalPrimitiveParameterization2 primitive =
            point_bisector(
                first,
                second,
                lower,
                upper);
        append_exact_graph_components(
            dual_id,
            primitive_kind,
            generator_ids,
            clip_linear_clearance_components(
                dual_id,
                primitive,
                point_clearance_boundary(
                    primitive,
                    first.x,
                    first.y,
                    radius_squared),
                domain),
            graph,
            node_indices);
        if (first_node.has_value()) {
            register_node_alias(
                dual_id,
                *first_parameter,
                *first_node,
                aliases);
        }
        if (second_node.has_value()) {
            register_node_alias(
                dual_id,
                *second_parameter,
                *second_node,
                aliases);
        }
    }
}

MatExactOpenSegmentSource2 exact_segment_source(
    const NormalizedOpenSegmentSource2& segment,
    const std::vector<NormalizedPointSource2>& points)
{
    const NormalizedPointSource2& source =
        point_source(points, segment.source_point_id);
    const NormalizedPointSource2& target =
        point_source(points, segment.target_point_id);
    const CORE::BigRat line_a = source.y - target.y;
    const CORE::BigRat line_b = target.x - source.x;
    if (line_a == 0 && line_b == 0) {
        throw InvalidRationalPrimitiveError(
            "normalized open segment is degenerate");
    }
    return {
        segment.stable_site_id,
        line_a,
        line_b,
        source.x * target.y - target.x * source.y,
    };
}

std::pair<
    std::pair<CORE::BigRat, std::string>,
    std::pair<CORE::BigRat, std::string>>
source_parameter_bounds(
    const MatExactPointSiteSource2& focus,
    const MatExactOpenSegmentSource2& segment,
    const NormalizedOpenSegmentSource2& normalized,
    const std::vector<NormalizedPointSource2>& points)
{
    if (focus.x.radical != 0 || focus.y.radical != 0) {
        throw InvalidRationalPrimitiveError(
            "source endpoint bounds require rational focus");
    }
    const CORE::BigRat line_norm =
        segment.line_a * segment.line_a
        + segment.line_b * segment.line_b;
    const CORE::BigRat signed_distance =
        segment.line_a * focus.x.rational
        + segment.line_b * focus.y.rational
        + segment.line_c;
    const CORE::BigRat projection_x =
        focus.x.rational
        - signed_distance * segment.line_a / line_norm;
    const CORE::BigRat projection_y =
        focus.y.rational
        - signed_distance * segment.line_b / line_norm;
    const auto endpoint_parameter =
        [&segment,
         &points,
         &projection_x,
         &projection_y,
         &line_norm](const std::string& point_id) {
            const NormalizedPointSource2& endpoint =
                point_source(points, point_id);
            return (
                (endpoint.x - projection_x)
                    * -segment.line_b
                + (endpoint.y - projection_y)
                    * segment.line_a)
                / line_norm;
        };
    std::pair<CORE::BigRat, std::string> source{
        endpoint_parameter(normalized.source_point_id),
        normalized.source_point_id,
    };
    std::pair<CORE::BigRat, std::string> target{
        endpoint_parameter(normalized.target_point_id),
        normalized.target_point_id,
    };
    if (target.first < source.first) {
        std::swap(source, target);
    }
    return {std::move(source), std::move(target)};
}

void add_bound_provenance(
    std::vector<MatAdmissibleComponent2>& components,
    const CORE::BigRat& lower,
    const std::string& lower_id,
    const CORE::BigRat& upper,
    const std::string& upper_id)
{
    ExactAlgebraicKernel1 kernel;
    const auto compare = kernel.compare_1_object();
    for (MatAdmissibleComponent2& component : components) {
        for (MatParameterEndpoint2* endpoint :
             {&component.lower, &component.upper}) {
            if (compare(*endpoint->parameter, lower)
                == CGAL::EQUAL) {
                endpoint->provenance_ids.push_back(lower_id);
            }
            if (compare(*endpoint->parameter, upper)
                == CGAL::EQUAL) {
                endpoint->provenance_ids.push_back(upper_id);
            }
            std::sort(
                endpoint->provenance_ids.begin(),
                endpoint->provenance_ids.end());
            endpoint->provenance_ids.erase(
                std::unique(
                    endpoint->provenance_ids.begin(),
                    endpoint->provenance_ids.end()),
                endpoint->provenance_ids.end());
        }
    }
}

void append_point_segment_graph(
    const std::vector<NormalizedPointSource2>& points,
    const NormalizedOpenSegmentSource2& segment_record,
    const std::string& focus_id,
    const MatDomainPolygonWithHoles2& domain,
    MatExactGraph2& graph,
    std::map<std::string, std::size_t>& node_indices)
{
    const NormalizedPointSource2& source =
        point_source(points, segment_record.source_point_id);
    const NormalizedPointSource2& target =
        point_source(points, segment_record.target_point_id);
    const NormalizedPointSource2& focus =
        point_source(points, focus_id);
    const MatTraits::Point_2 source_point(source.x, source.y);
    const MatTraits::Point_2 target_point(target.x, target.y);
    const MatTraits::Point_2 focus_point(focus.x, focus.y);
    const std::vector<GeneratorSite2> generators{
        {
            source.stable_site_id,
            MatTraits::Site_2::construct_site_2(
                source_point),
        },
        {
            target.stable_site_id,
            MatTraits::Site_2::construct_site_2(
                target_point),
        },
        {
            segment_record.stable_site_id,
            MatTraits::Site_2::construct_site_2(
                source_point,
                target_point),
        },
        {
            focus.stable_site_id,
            MatTraits::Site_2::construct_site_2(
                focus_point),
        },
    };
    SegmentSiteDelaunay2 delaunay;
    delaunay.insert(source_point, target_point);
    delaunay.insert(focus_point);
    const MatExactPointSiteSource2 exact_focus{
        focus.stable_site_id,
        {focus.x, 0},
        {focus.y, 0},
        1,
    };
    const MatExactOpenSegmentSource2 exact_segment =
        exact_segment_source(segment_record, points);
    const auto bounds = source_parameter_bounds(
        exact_focus,
        exact_segment,
        segment_record,
        points);

    for (auto edge = delaunay.finite_edges_begin();
         edge != delaunay.finite_edges_end();
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
                delaunay.primal(edge))) {
            continue;
        }
        const std::vector<std::string> generator_ids =
            ordered_generator_site_ids(
                stable_generator_site_id(
                    point,
                    generators),
                stable_generator_site_id(
                    segment,
                    generators));
        const std::string dual_id =
            stable_dual_identity_v1(
                "point-segment",
                generator_ids);
        std::vector<MatAdmissibleComponent2> components =
            clip_source_parabola_clearance_components(
                dual_id,
                exact_focus,
                exact_segment,
                bounds.first.first,
                bounds.second.first,
                {CGAL::POSITIVE, {}, {}},
                domain);
        add_bound_provenance(
            components,
            bounds.first.first,
            bounds.first.second,
            bounds.second.first,
            bounds.second.second);
        append_exact_graph_components(
            dual_id,
            "PARABOLA",
            generator_ids,
            components,
            graph,
            node_indices);
    }
}

} // namespace

MatExactGraph2 exact_point_site_graph(
    const std::vector<NormalizedPointSource2>& points,
    const MatDomainPolygonWithHoles2& domain,
    const CORE::BigRat& radius_squared)
{
    require_valid_point_graph_inputs(points, radius_squared);
    MatExactGraph2 graph{{}, {}, 0, 0};
    std::map<std::string, std::size_t> node_indices;
    std::map<std::string, CanonicalNodeAlias2> aliases;
    append_point_site_graph(
        points,
        domain,
        radius_squared,
        graph,
        node_indices,
        aliases);
    canonicalize_original_nodes(graph, aliases);
    const auto collapsed_edge = std::find_if(
        graph.edges.begin(),
        graph.edges.end(),
        [](const MatExactGraphEdge2& edge)
        {
            return edge.source_node_id == edge.target_node_id;
        });
    if (collapsed_edge != graph.edges.end()) {
        throw DegeneratePointSiteTopologyError(
            "point-site graph contains a collapsed exact dual");
    }
    std::sort(
        graph.nodes.begin(),
        graph.nodes.end(),
        [](const MatExactGraphNode2& lhs,
           const MatExactGraphNode2& rhs)
        {
            return lhs.node_id < rhs.node_id;
        });
    std::sort(
        graph.edges.begin(),
        graph.edges.end(),
        [](const MatExactGraphEdge2& lhs,
           const MatExactGraphEdge2& rhs)
        {
            return lhs.edge_id < rhs.edge_id;
        });
    return graph;
}

MatExactGraph2 segment_site_generic_graph_spike()
{
    MatDomainPolygon2 linear_outer;
    linear_outer.push_back({-10, -10});
    linear_outer.push_back({10, -10});
    linear_outer.push_back({10, 10});
    linear_outer.push_back({-10, 10});
    const MatDomainPolygonWithHoles2 linear_domain(
        linear_outer);
    MatExactGraph2 graph = exact_point_site_graph(
        {
            {"linear-a", 0, 0},
            {"linear-b", 8, 0},
            {"linear-c", 0, 8},
            {"linear-interior", 2, 2},
        },
        linear_domain,
        0);
    std::map<std::string, std::size_t> node_indices;
    for (std::size_t index = 0; index < graph.nodes.size(); ++index) {
        node_indices.emplace(graph.nodes[index].node_id, index);
    }

    MatDomainPolygon2 parabola_outer;
    parabola_outer.push_back({0, 1});
    parabola_outer.push_back({4, 1});
    parabola_outer.push_back(
        {4, CORE::BigRat(13, 6)});
    parabola_outer.push_back(
        {0, CORE::BigRat(13, 6)});
    const MatDomainPolygonWithHoles2 parabola_domain(
        parabola_outer);
    append_point_segment_graph(
        {
            {"segment-source", 0, 0},
            {"segment-target", 4, 0},
            {"focus", 2, 3},
        },
        {
            "open-segment",
            "segment-source",
            "segment-target",
        },
        "focus",
        parabola_domain,
        graph,
        node_indices);

    return graph;
}
