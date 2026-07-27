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

const NormalizedOpenSegmentSource2& normalized_segment_source(
    const std::vector<NormalizedOpenSegmentSource2>& segments,
    const std::string& stable_site_id)
{
    const auto found = std::find_if(
        segments.begin(),
        segments.end(),
        [&stable_site_id](
            const NormalizedOpenSegmentSource2& segment) {
            return segment.stable_site_id == stable_site_id;
        });
    if (found == segments.end()) {
        throw InvalidRationalPrimitiveError(
            "segment site has no normalized source record");
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
    std::vector<std::string> parent_site_ids;
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
            union_stable_ids(
                node.parent_site_ids,
                alias->second.parent_site_ids);
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
            retained.parent_site_ids,
            node.parent_site_ids);
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
    return canonical_open_segment_source(
        segment.stable_site_id,
        source.x,
        source.y,
        target.x,
        target.y);
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

MatExactPointSiteSource2 exact_point_site_source(
    const NormalizedPointSource2& point)
{
    return {
        point.stable_site_id,
        {point.x, 0},
        {point.y, 0},
        1,
    };
}

MatParameterEndpoint2 rational_endpoint(
    const std::pair<CORE::BigRat, std::string>& bound)
{
    ExactAlgebraicKernel1 kernel;
    return {
        kernel.construct_algebraic_real_1_object()(
            bound.first),
        {bound.second},
    };
}

MatParameterEndpoint2 bind_point_segment_endpoint_owner(
    const MatTraits::Site_2& owner,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge,
    const MatExactPointSiteSource2& focus,
    const MatExactOpenSegmentSource2& segment,
    const MatExactPointSiteSource2& segment_source,
    const MatExactPointSiteSource2& segment_target,
    const std::pair<
        std::pair<CORE::BigRat, std::string>,
        std::pair<CORE::BigRat, std::string>>&
        feature_bounds,
    const std::vector<NormalizedPointSource2>& points,
    const std::vector<NormalizedOpenSegmentSource2>& segments,
    const std::vector<GeneratorSite2>& generators)
{
    const std::string owner_id =
        stable_generator_site_id(owner, generators);
    if (owner_id == feature_bounds.first.second) {
        return rational_endpoint(feature_bounds.first);
    }
    if (owner_id == feature_bounds.second.second) {
        return rational_endpoint(feature_bounds.second);
    }
    if (owner.is_point()) {
        return bind_point_limiter_parabola_endpoint(
            focus,
            segment,
            segment_source,
            segment_target,
            exact_point_site_source(
                point_source(points, owner_id)),
            voronoi,
            halfedge);
    }
    const NormalizedOpenSegmentSource2& limiter =
        normalized_segment_source(segments, owner_id);
    const NormalizedPointSource2& limiter_source =
        point_source(points, limiter.source_point_id);
    const NormalizedPointSource2& limiter_target =
        point_source(points, limiter.target_point_id);
    return bind_segment_limiter_parabola_endpoint(
        focus,
        segment,
        segment_source,
        segment_target,
        exact_segment_source(limiter, points),
        exact_point_site_source(limiter_source),
        exact_point_site_source(limiter_target),
        voronoi,
        halfedge);
}

std::pair<MatParameterEndpoint2, MatParameterEndpoint2>
bind_point_segment_cell_endpoints(
    const SegmentSiteVoronoi2& voronoi,
    const MatExactPointSiteSource2& focus,
    const MatExactOpenSegmentSource2& segment,
    const MatExactPointSiteSource2& segment_source,
    const MatExactPointSiteSource2& segment_target,
    const std::pair<
        std::pair<CORE::BigRat, std::string>,
        std::pair<CORE::BigRat, std::string>>&
        feature_bounds,
    const std::vector<NormalizedPointSource2>& points,
    const std::vector<NormalizedOpenSegmentSource2>& segments,
    const std::vector<GeneratorSite2>& generators)
{
    const std::vector<std::string> expected_generators =
        ordered_generator_site_ids(
            focus.stable_site_id,
            segment.stable_site_id);
    for (auto halfedge = voronoi.halfedges_begin();
         halfedge != voronoi.halfedges_end();
         ++halfedge) {
        const std::string up_id =
            stable_generator_site_id(
                halfedge->up()->site(),
                generators);
        const std::string down_id =
            stable_generator_site_id(
                halfedge->down()->site(),
                generators);
        if (ordered_generator_site_ids(up_id, down_id)
                != expected_generators
            || up_id != expected_generators.front()) {
            continue;
        }
        if (!halfedge->has_source()
            || !halfedge->has_target()) {
            throw UnboundLiveParabolaEndpointError(
                "P-S graph halfedge is not bounded by two exact owners");
        }

        MatParameterEndpoint2 source =
            bind_point_segment_endpoint_owner(
                halfedge->left()->site(),
                voronoi,
                halfedge,
                focus,
                segment,
                segment_source,
                segment_target,
                feature_bounds,
                points,
                segments,
                generators);
        MatParameterEndpoint2 target =
            bind_point_segment_endpoint_owner(
                halfedge->right()->site(),
                voronoi,
                halfedge,
                focus,
                segment,
                segment_source,
                segment_target,
                feature_bounds,
                points,
                segments,
                generators);
        ExactAlgebraicKernel1 kernel;
        if (kernel.compare_1_object()(
                *target.parameter,
                *source.parameter)
            == CGAL::SMALLER) {
            std::swap(source, target);
        }
        if (kernel.compare_1_object()(
                *source.parameter,
                *target.parameter)
            != CGAL::SMALLER) {
            throw AmbiguousLiveParabolaEndpointError(
                "P-S graph halfedge endpoints are not strictly ordered");
        }
        return {
            std::move(source),
            std::move(target),
        };
    }
    throw UnboundLiveParabolaEndpointError(
        "P-S graph has no canonical live halfedge");
}

void append_point_segment_graph(
    const std::vector<NormalizedPointSource2>& points,
    const std::vector<NormalizedOpenSegmentSource2>& segments,
    const std::string& source_segment_id,
    const std::string& focus_id,
    const MatDomainPolygonWithHoles2& domain,
    const CORE::BigRat& radius_squared,
    MatExactGraph2& graph,
    std::map<std::string, std::size_t>& node_indices)
{
    const NormalizedOpenSegmentSource2& segment_record =
        normalized_segment_source(
            segments,
            source_segment_id);
    const NormalizedPointSource2& source =
        point_source(points, segment_record.source_point_id);
    const NormalizedPointSource2& target =
        point_source(points, segment_record.target_point_id);
    const NormalizedPointSource2& focus =
        point_source(points, focus_id);
    std::vector<GeneratorSite2> generators;
    generators.reserve(points.size() + segments.size());
    for (const NormalizedPointSource2& point : points) {
        generators.push_back(
            {
                point.stable_site_id,
                MatTraits::Site_2::construct_site_2(
                    MatTraits::Point_2(
                        point.x,
                        point.y)),
            });
    }
    for (const NormalizedOpenSegmentSource2& segment :
         segments) {
        const NormalizedPointSource2& segment_source_point =
            point_source(points, segment.source_point_id);
        const NormalizedPointSource2& segment_target_point =
            point_source(points, segment.target_point_id);
        generators.push_back(
            {
                segment.stable_site_id,
                MatTraits::Site_2::construct_site_2(
                    MatTraits::Point_2(
                        segment_source_point.x,
                        segment_source_point.y),
                    MatTraits::Point_2(
                        segment_target_point.x,
                        segment_target_point.y)),
            });
    }
    SegmentSiteDelaunay2 delaunay;
    for (const NormalizedOpenSegmentSource2& segment :
         segments) {
        const NormalizedPointSource2& segment_source_point =
            point_source(points, segment.source_point_id);
        const NormalizedPointSource2& segment_target_point =
            point_source(points, segment.target_point_id);
        delaunay.insert(
            MatTraits::Point_2(
                segment_source_point.x,
                segment_source_point.y),
            MatTraits::Point_2(
                segment_target_point.x,
                segment_target_point.y));
    }
    for (const NormalizedPointSource2& point : points) {
        const bool is_segment_endpoint =
            std::any_of(
                segments.begin(),
                segments.end(),
                [&point](
                    const NormalizedOpenSegmentSource2& segment) {
                    return point.stable_site_id
                            == segment.source_point_id
                        || point.stable_site_id
                            == segment.target_point_id;
                });
        if (is_segment_endpoint) {
            continue;
        }
        delaunay.insert(
            MatTraits::Point_2(point.x, point.y));
    }
    graph.matched_generator_sites +=
        require_generator_site_bijection(
            delaunay,
            generators);
    SegmentSiteVoronoi2 voronoi(delaunay);
    const MatExactPointSiteSource2 exact_focus =
        exact_point_site_source(focus);
    const MatExactOpenSegmentSource2 exact_segment =
        exact_segment_source(segment_record, points);
    const auto bounds = source_parameter_bounds(
        exact_focus,
        exact_segment,
        segment_record,
        points);
    const ClearanceRootBoundary2 clearance_boundary =
        source_parabola_clearance_boundary(
            exact_focus,
            exact_segment,
            radius_squared);

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
        const std::vector<std::string> generator_ids =
            ordered_generator_site_ids(
                stable_generator_site_id(
                    point,
                    generators),
                stable_generator_site_id(
                    segment,
                    generators));
        if (generator_ids
            != ordered_generator_site_ids(
                focus.stable_site_id,
                segment_record.stable_site_id)) {
            continue;
        }
        SegmentSiteParabola2 live_parabola;
        if (!CGAL::assign(
                live_parabola,
                delaunay.primal(edge))) {
            throw MismatchedLiveParabolaBridgeError(
                "nonincident P-S dual is not a parabola");
        }
        const std::string dual_id =
            stable_dual_identity_v1(
                "point-segment",
                generator_ids);
        const auto owned_endpoints =
            bind_point_segment_cell_endpoints(
                voronoi,
                exact_focus,
                exact_segment,
                exact_point_site_source(source),
                exact_point_site_source(target),
                bounds,
                points,
                segments,
                generators);
        const std::vector<MatAdmissibleComponent2> components =
            clip_source_parabola_clearance_components(
                dual_id,
                exact_focus,
                exact_segment,
                owned_endpoints.first,
                owned_endpoints.second,
                clearance_boundary,
                domain);
        append_exact_graph_components(
            dual_id,
            "PARABOLA",
            generator_ids,
            components,
            graph,
            node_indices);
    }
}

struct ParallelSegmentFeatureCell2 {
    RationalPrimitiveParameterization2 primitive;
    RationalDomainRoot2 lower;
    RationalDomainRoot2 upper;
};

ParallelSegmentFeatureCell2 parallel_segment_feature_cell(
    const MatExactOpenSegmentSource2& first,
    const NormalizedOpenSegmentSource2& first_record,
    const MatExactOpenSegmentSource2& second,
    const NormalizedOpenSegmentSource2& second_record,
    const std::vector<NormalizedPointSource2>& points)
{
    RationalPrimitiveParameterization2 primitive =
        parallel_segment_bisector_parameterization(
            first,
            second);
    struct SourceEndpointParameter2 {
        CORE::BigRat parameter;
        std::string stable_id;
        std::string segment_id;
    };
    const auto endpoint_parameter =
        [&primitive, &points](
            const std::string& endpoint_id,
            const std::string& segment_id) {
            const NormalizedPointSource2& endpoint =
                point_source(points, endpoint_id);
            return SourceEndpointParameter2{
                parallel_segment_tangent_parameter(
                    primitive,
                    endpoint.x,
                    endpoint.y),
                endpoint.stable_site_id,
                segment_id,
            };
        };
    std::vector<SourceEndpointParameter2> endpoints{
        endpoint_parameter(
            first_record.source_point_id,
            first_record.stable_site_id),
        endpoint_parameter(
            first_record.target_point_id,
            first_record.stable_site_id),
        endpoint_parameter(
            second_record.source_point_id,
            second_record.stable_site_id),
        endpoint_parameter(
            second_record.target_point_id,
            second_record.stable_site_id),
    };
    const auto segment_bounds =
        [&endpoints](const std::string& segment_id) {
            std::vector<CORE::BigRat> values;
            for (const SourceEndpointParameter2& endpoint :
                 endpoints) {
                if (endpoint.segment_id == segment_id) {
                    values.push_back(endpoint.parameter);
                }
            }
            if (values.size() != 2) {
                throw InvalidRationalPrimitiveError(
                    "parallel S-S source has incomplete endpoints");
            }
            if (values[1] < values[0]) {
                std::swap(values[0], values[1]);
            }
            return std::pair<CORE::BigRat, CORE::BigRat>{
                values[0],
                values[1],
            };
        };
    const auto first_bounds =
        segment_bounds(first_record.stable_site_id);
    const auto second_bounds =
        segment_bounds(second_record.stable_site_id);
    const CORE::BigRat lower =
        std::max(first_bounds.first, second_bounds.first);
    const CORE::BigRat upper =
        std::min(first_bounds.second, second_bounds.second);
    if (lower >= upper) {
        throw EmptyParallelSegmentFeatureDomainError(
            "parallel open segments have no positive-length overlap");
    }
    const auto provenance_at =
        [&endpoints](const CORE::BigRat& parameter) {
            std::vector<std::string> provenance;
            for (const SourceEndpointParameter2& endpoint :
                 endpoints) {
                if (endpoint.parameter == parameter) {
                    provenance.push_back(endpoint.stable_id);
                }
            }
            union_stable_ids(provenance, {});
            return provenance;
        };
    primitive.domain_lower = lower;
    primitive.domain_upper = upper;
    return {
        std::move(primitive),
        {
            lower,
            provenance_at(lower),
        },
        {
            upper,
            provenance_at(upper),
        },
    };
}

void append_parallel_segment_segment_graph(
    const std::vector<NormalizedPointSource2>& points,
    const std::vector<NormalizedOpenSegmentSource2>& segments,
    const std::string& first_segment_id,
    const std::string& second_segment_id,
    const MatDomainPolygonWithHoles2& domain,
    const CORE::BigRat& radius_squared,
    const std::vector<MatExactPointSiteSource2>& point_limiters,
    const std::vector<MatExactOpenSegmentSource2>& segment_limiters,
    MatExactGraph2& graph,
    std::map<std::string, std::size_t>& node_indices)
{
    const std::vector<std::string> generator_ids =
        ordered_generator_site_ids(
            first_segment_id,
            second_segment_id);
    const NormalizedOpenSegmentSource2& first_record =
        normalized_segment_source(
            segments,
            generator_ids[0]);
    const NormalizedOpenSegmentSource2& second_record =
        normalized_segment_source(
            segments,
            generator_ids[1]);
    const MatExactOpenSegmentSource2 first =
        exact_segment_source(first_record, points);
    const MatExactOpenSegmentSource2 second =
        exact_segment_source(second_record, points);
    const ParallelSegmentFeatureCell2 feature =
        parallel_segment_feature_cell(
            first,
            first_record,
            second,
            second_record,
            points);

    std::vector<GeneratorSite2> generators;
    generators.reserve(points.size() + segments.size());
    for (const NormalizedPointSource2& point : points) {
        generators.push_back(
            {
                point.stable_site_id,
                MatTraits::Site_2::construct_site_2(
                    MatTraits::Point_2(point.x, point.y)),
            });
    }
    for (const NormalizedOpenSegmentSource2& segment :
         segments) {
        const NormalizedPointSource2& source =
            point_source(points, segment.source_point_id);
        const NormalizedPointSource2& target =
            point_source(points, segment.target_point_id);
        generators.push_back(
            {
                segment.stable_site_id,
                MatTraits::Site_2::construct_site_2(
                    MatTraits::Point_2(source.x, source.y),
                    MatTraits::Point_2(target.x, target.y)),
            });
    }

    SegmentSiteDelaunay2 delaunay;
    for (const NormalizedOpenSegmentSource2& segment :
         segments) {
        const NormalizedPointSource2& source =
            point_source(points, segment.source_point_id);
        const NormalizedPointSource2& target =
            point_source(points, segment.target_point_id);
        delaunay.insert(
            MatTraits::Point_2(source.x, source.y),
            MatTraits::Point_2(target.x, target.y));
    }
    for (const NormalizedPointSource2& point : points) {
        const bool is_segment_endpoint =
            std::any_of(
                segments.begin(),
                segments.end(),
                [&point](
                    const NormalizedOpenSegmentSource2& segment) {
                    return point.stable_site_id
                            == segment.source_point_id
                        || point.stable_site_id
                            == segment.target_point_id;
                });
        if (!is_segment_endpoint) {
            delaunay.insert(
                MatTraits::Point_2(point.x, point.y));
        }
    }
    const std::size_t matched_generator_sites =
        require_generator_site_bijection(
            delaunay,
            generators);
    SegmentSiteVoronoi2 voronoi(delaunay);
    const ClearanceRootBoundary2 boundary =
        parallel_segment_clearance_boundary(
            feature.primitive,
            first,
            second,
            radius_squared);
    const std::string dual_id =
        stable_dual_identity_v1(
            "segment-segment",
            generator_ids);

    std::vector<SegmentSiteVoronoi2::Halfedge_handle>
        matching_halfedges;
    for (auto halfedge = voronoi.halfedges_begin();
         halfedge != voronoi.halfedges_end();
         ++halfedge) {
        const MatTraits::Site_2 up =
            halfedge->up()->site();
        const MatTraits::Site_2 down =
            halfedge->down()->site();
        if (!up.is_segment() || !down.is_segment()) {
            continue;
        }
        const std::string up_id =
            stable_generator_site_id(up, generators);
        const std::string down_id =
            stable_generator_site_id(down, generators);
        if (ordered_generator_site_ids(up_id, down_id)
                != generator_ids
            || up_id != generator_ids.front()) {
            continue;
        }
        matching_halfedges.push_back(halfedge);
    }
    if (matching_halfedges.size() != 1) {
        throw UnboundLiveSegmentSegmentEndpointError(
            "parallel S-S graph has no unique canonical halfedge");
    }
    const auto owned_endpoints =
        bind_parallel_segment_segment_cell_endpoints(
            feature.primitive,
            feature.lower,
            feature.upper,
            first,
            second,
            point_limiters,
            segment_limiters,
            generator_ids,
            generators,
            voronoi,
            matching_halfedges.front());
    std::vector<MatAdmissibleComponent2> components =
        clip_bounded_linear_clearance_components(
            dual_id,
            feature.primitive,
            owned_endpoints.first,
            owned_endpoints.second,
            boundary,
            domain);
    graph.matched_generator_sites +=
        matched_generator_sites;
    append_exact_graph_components(
        dual_id,
        "LINE",
        generator_ids,
        components,
        graph,
        node_indices);
}

bool same_unoriented_segment(
    const MatTraits::Segment_2& lhs,
    const MatTraits::Segment_2& rhs)
{
    return (lhs.source() == rhs.source()
            && lhs.target() == rhs.target())
        || (lhs.source() == rhs.target()
            && lhs.target() == rhs.source());
}

std::vector<GeneratorSite2> segment_site_generators(
    const std::vector<NormalizedPointSource2>& points,
    const std::vector<NormalizedOpenSegmentSource2>& segments)
{
    std::vector<GeneratorSite2> generators;
    generators.reserve(points.size() + segments.size());
    for (const NormalizedPointSource2& point : points) {
        generators.push_back(
            {
                point.stable_site_id,
                MatTraits::Site_2::construct_site_2(
                    MatTraits::Point_2(point.x, point.y)),
            });
    }
    for (const NormalizedOpenSegmentSource2& segment :
         segments) {
        const NormalizedPointSource2& source =
            point_source(points, segment.source_point_id);
        const NormalizedPointSource2& target =
            point_source(points, segment.target_point_id);
        generators.push_back(
            {
                segment.stable_site_id,
                MatTraits::Site_2::construct_site_2(
                    MatTraits::Point_2(source.x, source.y),
                    MatTraits::Point_2(target.x, target.y)),
            });
    }
    return generators;
}

void insert_segment_site_sources(
    const std::vector<NormalizedPointSource2>& points,
    const std::vector<NormalizedOpenSegmentSource2>& segments,
    SegmentSiteDelaunay2& delaunay)
{
    for (const NormalizedOpenSegmentSource2& segment :
         segments) {
        const NormalizedPointSource2& source =
            point_source(points, segment.source_point_id);
        const NormalizedPointSource2& target =
            point_source(points, segment.target_point_id);
        delaunay.insert(
            MatTraits::Point_2(source.x, source.y),
            MatTraits::Point_2(target.x, target.y));
    }
    for (const NormalizedPointSource2& point : points) {
        const bool is_segment_endpoint = std::any_of(
            segments.begin(),
            segments.end(),
            [&point](
                const NormalizedOpenSegmentSource2& segment) {
                return point.stable_site_id
                        == segment.source_point_id
                    || point.stable_site_id
                        == segment.target_point_id;
            });
        if (!is_segment_endpoint) {
            delaunay.insert(
                MatTraits::Point_2(point.x, point.y));
        }
    }
}

std::vector<std::string> raw_feature_parent_ids(
    const std::string& feature_id,
    const std::vector<NormalizedOpenSegmentSource2>& segments)
{
    const auto segment = std::find_if(
        segments.begin(),
        segments.end(),
        [&feature_id](
            const NormalizedOpenSegmentSource2& candidate) {
            return candidate.stable_site_id == feature_id;
        });
    if (segment != segments.end()) {
        return {segment->stable_site_id};
    }

    std::vector<std::string> parents;
    for (const NormalizedOpenSegmentSource2& candidate :
         segments) {
        if (candidate.source_point_id == feature_id
            || candidate.target_point_id == feature_id) {
            parents.push_back(candidate.stable_site_id);
        }
    }
    if (parents.empty()) {
        return {feature_id};
    }
    std::sort(parents.begin(), parents.end());
    parents.erase(
        std::unique(parents.begin(), parents.end()),
        parents.end());
    return parents;
}

struct RawFeatureParentDisposition2 {
    bool maps_to_target;
    bool self_transition;
};

RawFeatureParentDisposition2 classify_raw_feature_pair(
    const std::string& first_feature_id,
    const std::string& second_feature_id,
    const std::vector<std::string>& target_parent_ids,
    const std::vector<NormalizedOpenSegmentSource2>& segments)
{
    const std::vector<std::string> first_parents =
        raw_feature_parent_ids(first_feature_id, segments);
    const std::vector<std::string> second_parents =
        raw_feature_parent_ids(second_feature_id, segments);
    std::set<std::vector<std::string>> distinct_pairs;
    std::set<std::string> self_parents;
    for (const std::string& first : first_parents) {
        for (const std::string& second : second_parents) {
            if (first == second) {
                self_parents.insert(first);
            } else {
                distinct_pairs.insert(
                    ordered_generator_site_ids(first, second));
            }
        }
    }
    if (!self_parents.empty() && !distinct_pairs.empty()) {
        throw AmbiguousCompositeSiteOwnerError(
            "raw feature pair maps to both self and distinct parents");
    }
    if (distinct_pairs.size() > 1) {
        throw AmbiguousCompositeSiteOwnerError(
            "raw feature pair maps to multiple parent pairs");
    }
    if (!distinct_pairs.empty()) {
        return {
            *distinct_pairs.begin() == target_parent_ids,
            false,
        };
    }
    const bool target_self = std::any_of(
        self_parents.begin(),
        self_parents.end(),
        [&target_parent_ids](const std::string& parent) {
            return std::find(
                       target_parent_ids.begin(),
                       target_parent_ids.end(),
                       parent)
                != target_parent_ids.end();
        });
    return {false, target_self};
}

std::string nonparallel_original_dual_id(
    const int branch_sign,
    const std::vector<std::string>& parent_site_ids)
{
    if (branch_sign != -1 && branch_sign != 1) {
        throw IncompleteCompositeSegmentChainError(
            "nonparallel normalized branch sign is invalid");
    }
    return stable_dual_identity_v1(
        branch_sign < 0
            ? "segment-segment/branch-negative"
            : "segment-segment/branch-positive",
        parent_site_ids);
}

NonparallelSegmentFeatureDomain2
nonparallel_segment_feature_domain(
    const NonparallelSegmentBisectorParameterization2& primitive,
    const std::vector<NormalizedPointSource2>& points,
    const std::vector<NormalizedOpenSegmentSource2>& segments)
{
    const NormalizedOpenSegmentSource2& first_record =
        normalized_segment_source(
            segments,
            primitive.first_segment_id);
    const NormalizedOpenSegmentSource2& second_record =
        normalized_segment_source(
            segments,
            primitive.second_segment_id);
    const MatExactOpenSegmentSource2 first =
        exact_segment_source(first_record, points);
    const MatExactOpenSegmentSource2 second =
        exact_segment_source(second_record, points);
    const auto endpoint =
        [&primitive, &points](
            const MatExactOpenSegmentSource2& segment,
            const std::string& point_id) {
            const NormalizedPointSource2& point =
                point_source(points, point_id);
            return MatQuadraticFieldDomainBoundary2{
                nonparallel_segment_tangent_parameter(
                    primitive,
                    segment,
                    point.x,
                    point.y),
                {point_id},
            };
        };
    return intersect_nonparallel_segment_feature_domains(
        endpoint(first, first_record.source_point_id),
        endpoint(first, first_record.target_point_id),
        endpoint(second, second_record.source_point_id),
        endpoint(second, second_record.target_point_id),
        primitive.radicand);
}

MatQuadraticFieldValue2 multiply_quadratic_field_values(
    const MatQuadraticFieldValue2& lhs,
    const MatQuadraticFieldValue2& rhs,
    const CORE::BigRat& radicand)
{
    return {
        lhs.rational * rhs.rational
            + radicand * lhs.radical * rhs.radical,
        lhs.rational * rhs.radical
            + lhs.radical * rhs.rational,
    };
}

MatQuadraticFieldValue2 evaluate_affine_field_coordinate(
    const std::vector<CORE::BigRat>& rational,
    const std::vector<CORE::BigRat>& radical,
    const MatQuadraticFieldValue2& parameter,
    const CORE::BigRat& radicand)
{
    if (rational.size() != 2 || radical.size() != 2) {
        throw IncompleteCompositeSegmentChainError(
            "nonparallel segment chart is not affine");
    }
    const MatQuadraticFieldValue2 product =
        multiply_quadratic_field_values(
            {rational[1], radical[1]},
            parameter,
            radicand);
    return {
        rational[0] + product.rational,
        radical[0] + product.radical,
    };
}

CORE::Expr composite_core_field_value(
    const MatQuadraticFieldValue2& value,
    const CORE::BigRat& radicand)
{
    const CORE::Expr root =
        CORE::sqrt(CORE::Expr(radicand));
    return CORE::Expr(value.rational)
        + CORE::Expr(value.radical) * root;
}

MatTraits::Point_2 nonparallel_segment_point_at(
    const NonparallelSegmentBisectorParameterization2& primitive,
    const MatQuadraticFieldValue2& parameter)
{
    const MatQuadraticFieldValue2 x =
        evaluate_affine_field_coordinate(
            primitive.x_rational,
            primitive.x_radical,
            parameter,
            primitive.radicand);
    const MatQuadraticFieldValue2 y =
        evaluate_affine_field_coordinate(
            primitive.y_rational,
            primitive.y_radical,
            parameter,
            primitive.radicand);
    return {
        composite_core_field_value(x, primitive.radicand),
        composite_core_field_value(y, primitive.radicand),
    };
}

struct NormalizedNonparallelSegmentCell2 {
    MatTraits::Segment_2 representative;
    MatTraits::Segment_2 adaptor_span;
    NonparallelSegmentBisectorParameterization2 primitive;
    std::string original_dual_id;
};

struct CompositeParabolaPiece2 {
    SegmentSiteParabola2 parabola;
    std::string point_feature_id;
    std::string segment_feature_id;
    std::vector<std::string> generator_site_ids;
    std::size_t normalized_cell_index;
};

struct EmittedExactPoint2 {
    MatTraits::Point_2 point;
    std::string node_id;
};

std::size_t composite_cell_for_parabola(
    const SegmentSiteParabola2& parabola,
    const std::vector<NormalizedNonparallelSegmentCell2>& cells)
{
    std::optional<std::size_t> match;
    for (std::size_t index = 0;
         index < cells.size();
         ++index) {
        const MatTraits::Segment_2& representative =
            cells[index].representative;
        std::vector<MatTraits::Point_2> shared;
        for (const MatTraits::Point_2& segment_point :
             {representative.source(),
              representative.target()}) {
            if (segment_point == parabola.p1
                || segment_point == parabola.p2) {
                shared.push_back(segment_point);
            }
        }
        if (shared.empty()) {
            continue;
        }
        if (shared.size() != 1) {
            throw IncompleteCompositeSegmentChainError(
                "raw parabola shares multiple endpoints with one S-S piece");
        }
        const MatTraits::Point_2 segment_terminal =
            representative.source() == shared.front()
            ? representative.target()
            : representative.source();
        const MatTraits::Point_2 parabola_terminal =
            parabola.p1 == shared.front()
            ? parabola.p2
            : parabola.p1;
        if (!same_unoriented_segment(
                MatTraits::Segment_2(
                    segment_terminal,
                    parabola_terminal),
                cells[index].adaptor_span)) {
            continue;
        }
        if (match.has_value()) {
            throw IncompleteCompositeSegmentChainError(
                "raw parabola completes multiple normalized cells");
        }
        match = index;
    }
    if (!match.has_value()) {
        throw IncompleteCompositeSegmentChainError(
            "raw parabola completes no normalized cell");
    }
    return *match;
}

void register_composite_transition_aliases(
    const std::vector<EmittedExactPoint2>& segment_endpoints,
    const std::vector<EmittedExactPoint2>& parabola_endpoints,
    const std::vector<std::string>& parent_site_ids,
    const std::vector<std::string>& transition_feature_ids,
    std::map<std::string, CanonicalNodeAlias2>& aliases)
{
    for (const EmittedExactPoint2& parabola_endpoint :
         parabola_endpoints) {
        std::vector<const EmittedExactPoint2*> matches;
        for (const EmittedExactPoint2& segment_endpoint :
             segment_endpoints) {
            if (segment_endpoint.point
                == parabola_endpoint.point) {
                matches.push_back(&segment_endpoint);
            }
        }
        if (matches.size() != 1) {
            throw IncompleteCompositeSegmentChainError(
                "composite transition has no unique S-S endpoint");
        }
        const auto [existing, inserted] = aliases.emplace(
            matches.front()->node_id,
            CanonicalNodeAlias2{
                parabola_endpoint.node_id,
                transition_feature_ids,
                parent_site_ids,
            });
        if (!inserted
            && existing->second.node_id
                != parabola_endpoint.node_id) {
            throw IncompleteCompositeSegmentChainError(
                "composite transition aliases disagree");
        }
    }
}

void append_nonparallel_segment_segment_graph(
    const std::vector<NormalizedPointSource2>& points,
    const std::vector<NormalizedOpenSegmentSource2>& segments,
    const std::string& first_segment_id,
    const std::string& second_segment_id,
    const MatDomainPolygonWithHoles2& domain,
    const CORE::BigRat& radius_squared,
    MatExactGraph2& graph,
    std::map<std::string, std::size_t>& node_indices)
{
    if (radius_squared < 0) {
        throw NegativeClearanceRadiusSquaredError(
            "nonparallel S-S graph squared clearance radius is negative");
    }
    const std::vector<std::string> parent_site_ids =
        ordered_generator_site_ids(
            first_segment_id,
            second_segment_id);
    const MatExactOpenSegmentSource2 first_segment =
        exact_segment_source(
            normalized_segment_source(
                segments,
                parent_site_ids[0]),
            points);
    const MatExactOpenSegmentSource2 second_segment =
        exact_segment_source(
            normalized_segment_source(
                segments,
                parent_site_ids[1]),
            points);
    const std::vector<GeneratorSite2> generators =
        segment_site_generators(points, segments);
    SegmentSiteDelaunay2 delaunay;
    insert_segment_site_sources(
        points,
        segments,
        delaunay);
    graph.matched_generator_sites +=
        require_generator_site_bijection(
            delaunay,
            generators);
    SegmentSiteVoronoi2 voronoi(delaunay);

    std::vector<NormalizedNonparallelSegmentCell2> cells;
    for (auto halfedge = voronoi.halfedges_begin();
         halfedge != voronoi.halfedges_end();
         ++halfedge) {
        const MatTraits::Site_2 up = halfedge->up()->site();
        const MatTraits::Site_2 down =
            halfedge->down()->site();
        if (!up.is_segment() || !down.is_segment()) {
            continue;
        }
        const std::string up_id =
            stable_generator_site_id(up, generators);
        const std::string down_id =
            stable_generator_site_id(down, generators);
        if (ordered_generator_site_ids(up_id, down_id)
                != parent_site_ids
            || up_id != parent_site_ids.front()) {
            continue;
        }
        if (!halfedge->has_source()
            || !halfedge->has_target()) {
            throw IncompleteCompositeSegmentChainError(
                "normalized nonparallel S-S cell is unbounded");
        }
        MatTraits::Segment_2 representative;
        if (!CGAL::assign(
                representative,
                voronoi.dual().primal(
                    halfedge->dual()))) {
            throw UnsupportedCompositeSegmentPrimitiveError(
                "normalized nonparallel representative is not a segment");
        }
        NonparallelSegmentBisectorParameterization2 primitive =
            nonparallel_segment_bisector_parameterization(
                first_segment,
                second_segment,
                representative);
        cells.push_back(
            {
                representative,
                MatTraits::Segment_2(
                    halfedge->source()->point(),
                    halfedge->target()->point()),
                primitive,
                nonparallel_original_dual_id(
                    primitive.branch_sign,
                    parent_site_ids),
            });
    }
    if (cells.size() != 2
        || cells[0].original_dual_id
            == cells[1].original_dual_id) {
        throw IncompleteCompositeSegmentChainError(
            "nonparallel S-S pair has no two stable normalized cells");
    }

    std::vector<CompositeParabolaPiece2> parabolas;
    for (auto edge = delaunay.finite_edges_begin();
         edge != delaunay.finite_edges_end();
         ++edge) {
        const auto raw = *edge;
        const auto face = raw.first;
        const int index = raw.second;
        const MatTraits::Site_2 first =
            face->vertex(delaunay.ccw(index))->site();
        const MatTraits::Site_2 second =
            face->vertex(delaunay.cw(index))->site();
        const std::string first_id =
            stable_generator_site_id(first, generators);
        const std::string second_id =
            stable_generator_site_id(second, generators);
        const RawFeatureParentDisposition2 disposition =
            classify_raw_feature_pair(
                first_id,
                second_id,
                parent_site_ids,
                segments);
        const bool rejected =
            voronoi.edge_rejector()(delaunay, raw);
        if (disposition.self_transition) {
            if (rejected) {
                ++graph.rejected_incident_transitions;
            }
            continue;
        }
        if (!disposition.maps_to_target
            || !rejected) {
            continue;
        }
        if (first.is_point() == second.is_point()) {
            throw UnsupportedCompositeSegmentPrimitiveError(
                "retained composite transition is not point-segment");
        }
        SegmentSiteParabola2 parabola;
        if (!CGAL::assign(
                parabola,
                delaunay.primal(raw))) {
            throw UnsupportedCompositeSegmentPrimitiveError(
                "retained point-segment transition is not a parabola");
        }
        const std::string point_id =
            first.is_point() ? first_id : second_id;
        const std::string segment_id =
            first.is_segment() ? first_id : second_id;
        parabolas.push_back(
            {
                parabola,
                point_id,
                segment_id,
                ordered_generator_site_ids(
                    point_id,
                    segment_id),
                composite_cell_for_parabola(
                    parabola,
                    cells),
            });
    }
    if (parabolas.size() != 1
        || graph.rejected_incident_transitions != 1) {
        throw IncompleteCompositeSegmentChainError(
            "nonparallel fixture has no unique retained transition");
    }

    std::vector<EmittedExactPoint2> segment_endpoints;
    for (const NormalizedNonparallelSegmentCell2& cell :
         cells) {
        const NonparallelSegmentFeatureDomain2 feature_domain =
            nonparallel_segment_feature_domain(
                cell.primitive,
                points,
                segments);
        const std::vector<MatAdmissibleComponent2> components =
            clip_nonparallel_segment_clearance_components(
                cell.original_dual_id,
                cell.primitive,
                feature_domain,
                nonparallel_segment_clearance_boundary(
                    cell.primitive,
                    first_segment,
                    second_segment,
                    radius_squared),
                domain);
        if (components.size() != 1) {
            throw IncompleteCompositeSegmentChainError(
                "nonparallel fixture S-S cell did not retain one component");
        }
        const MatTraits::Point_2 lower_point =
            nonparallel_segment_point_at(
                cell.primitive,
                feature_domain.lower.parameter);
        const MatTraits::Point_2 upper_point =
            nonparallel_segment_point_at(
                cell.primitive,
                feature_domain.upper.parameter);
        if (!same_unoriented_segment(
                MatTraits::Segment_2(
                    lower_point,
                    upper_point),
                cell.representative)) {
            throw IncompleteCompositeSegmentChainError(
                "S-S feature domain and raw representative disagree");
        }
        append_exact_graph_components(
            cell.original_dual_id,
            cell.original_dual_id,
            "LINE",
            parent_site_ids,
            parent_site_ids,
            components,
            graph,
            node_indices);
        segment_endpoints.push_back(
            {
                lower_point,
                stable_endpoint_node_identity_v1(
                    cell.original_dual_id,
                    components.front().lower),
            });
        segment_endpoints.push_back(
            {
                upper_point,
                stable_endpoint_node_identity_v1(
                    cell.original_dual_id,
                    components.front().upper),
            });
    }

    std::map<std::string, CanonicalNodeAlias2> aliases;
    for (const CompositeParabolaPiece2& piece :
         parabolas) {
        const NormalizedNonparallelSegmentCell2& cell =
            cells[piece.normalized_cell_index];
        const NormalizedPointSource2& focus_record =
            point_source(points, piece.point_feature_id);
        const NormalizedOpenSegmentSource2& source_record =
            normalized_segment_source(
                segments,
                piece.segment_feature_id);
        const std::string limiter_id =
            parent_site_ids[0] == piece.segment_feature_id
            ? parent_site_ids[1]
            : parent_site_ids[0];
        const NormalizedOpenSegmentSource2& limiter_record =
            normalized_segment_source(
                segments,
                limiter_id);
        const MatExactPointSiteSource2 focus =
            exact_point_site_source(focus_record);
        const MatExactOpenSegmentSource2 source =
            exact_segment_source(source_record, points);
        const MatExactOpenSegmentSource2 limiter =
            exact_segment_source(limiter_record, points);
        const MatExactPointSiteSource2 source_point =
            exact_point_site_source(
                point_source(
                    points,
                    source_record.source_point_id));
        const MatExactPointSiteSource2 source_target =
            exact_point_site_source(
                point_source(
                    points,
                    source_record.target_point_id));
        const MatExactPointSiteSource2 limiter_source =
            exact_point_site_source(
                point_source(
                    points,
                    limiter_record.source_point_id));
        const MatExactPointSiteSource2 limiter_target =
            exact_point_site_source(
                point_source(
                    points,
                    limiter_record.target_point_id));
        EmittedExactPoint2 first_endpoint{
            piece.parabola.p1,
            {},
        };
        EmittedExactPoint2 second_endpoint{
            piece.parabola.p2,
            {},
        };
        MatParameterEndpoint2 first_bound =
            bind_segment_limiter_parabola_endpoint(
                focus,
                source,
                source_point,
                source_target,
                limiter,
                limiter_source,
                limiter_target,
                piece.parabola,
                piece.parabola.p1);
        MatParameterEndpoint2 second_bound =
            bind_segment_limiter_parabola_endpoint(
                focus,
                source,
                source_point,
                source_target,
                limiter,
                limiter_source,
                limiter_target,
                piece.parabola,
                piece.parabola.p2);
        ExactAlgebraicKernel1 kernel;
        if (kernel.compare_1_object()(
                *second_bound.parameter,
                *first_bound.parameter)
            == CGAL::SMALLER) {
            std::swap(first_bound, second_bound);
            std::swap(first_endpoint, second_endpoint);
        }
        const std::string piece_dual_id =
            stable_dual_identity_v1(
                cell.primitive.branch_sign < 0
                    ? "point-segment/composite-branch-negative"
                    : "point-segment/composite-branch-positive",
                piece.generator_site_ids);
        const std::vector<MatAdmissibleComponent2> components =
            clip_source_parabola_clearance_components(
                piece_dual_id,
                focus,
                source,
                first_bound,
                second_bound,
                source_parabola_clearance_boundary(
                    focus,
                    source,
                    radius_squared),
                domain);
        if (components.size() != 1) {
            throw IncompleteCompositeSegmentChainError(
                "nonparallel fixture parabola did not retain one component");
        }
        append_exact_graph_components(
            piece_dual_id,
            cell.original_dual_id,
            "PARABOLA",
            piece.generator_site_ids,
            parent_site_ids,
            components,
            graph,
            node_indices);
        first_endpoint.node_id =
            stable_endpoint_node_identity_v1(
                piece_dual_id,
                first_bound);
        second_endpoint.node_id =
            stable_endpoint_node_identity_v1(
                piece_dual_id,
                second_bound);
        std::vector<std::string> transition_feature_ids =
            piece.generator_site_ids;
        union_stable_ids(
            transition_feature_ids,
            parent_site_ids);
        register_composite_transition_aliases(
            segment_endpoints,
            {first_endpoint, second_endpoint},
            parent_site_ids,
            transition_feature_ids,
            aliases);
    }
    canonicalize_original_nodes(graph, aliases);
    std::sort(
        graph.nodes.begin(),
        graph.nodes.end(),
        [](const MatExactGraphNode2& lhs,
           const MatExactGraphNode2& rhs) {
            return lhs.node_id < rhs.node_id;
        });
    std::sort(
        graph.edges.begin(),
        graph.edges.end(),
        [](const MatExactGraphEdge2& lhs,
           const MatExactGraphEdge2& rhs) {
            return lhs.edge_id < rhs.edge_id;
        });
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

MatExactGraph2 segment_site_live_graph_spike()
{
    MatDomainPolygon2 linear_outer;
    linear_outer.push_back({-1, -3});
    linear_outer.push_back({1, -3});
    linear_outer.push_back({1, 3});
    linear_outer.push_back({-1, 3});
    const MatDomainPolygonWithHoles2 linear_domain(
        linear_outer);
    MatExactGraph2 graph = exact_point_site_graph(
        {
            {"line-left", -2, 0},
            {"line-right", 2, 0},
        },
        linear_domain,
        5);
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
    const MatDomainPolygonWithHoles2 parabola_domain(
        parabola_outer,
        holes.begin(),
        holes.end());
    append_point_segment_graph(
        {
            {"segment-source", 0, 0},
            {"segment-target", 4, 0},
            {"focus", 2, 3},
        },
        {
            {
                "open-segment",
                "segment-source",
                "segment-target",
            },
        },
        "open-segment",
        "focus",
        parabola_domain,
        0,
        graph,
        node_indices);
    return graph;
}

MatExactGraph2
segment_site_true_radius_graph_spike()
{
    MatExactGraph2 graph{{}, {}, 0, 0};
    std::map<std::string, std::size_t> node_indices;
    MatDomainPolygon2 outer;
    outer.push_back({-5, 0});
    outer.push_back({5, 0});
    outer.push_back({5, 6});
    outer.push_back({-5, 6});
    const MatDomainPolygonWithHoles2 domain(outer);
    append_point_segment_graph(
        {
            {"radius-segment-source", -4, 0},
            {"radius-segment-target", 4, 0},
            {"radius-focus", 0, 2},
        },
        {
            {
                "radius-segment",
                "radius-segment-source",
                "radius-segment-target",
            },
        },
        "radius-segment",
        "radius-focus",
        domain,
        CORE::BigRat(9, 4),
        graph,
        node_indices);
    return graph;
}

namespace {

MatExactGraph2 segment_site_generic_graph_spike_impl(
    const bool reverse_segment_endpoints)
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
    parabola_outer.push_back({-5, 0});
    parabola_outer.push_back({5, 0});
    parabola_outer.push_back({5, 6});
    parabola_outer.push_back({-5, 6});
    const MatDomainPolygonWithHoles2 parabola_domain(
        parabola_outer);
    append_point_segment_graph(
        {
            {"segment-source", -4, 0},
            {"segment-target", 4, 0},
            {"focus", 0, 2},
            {"limiter", 2, 2},
        },
        {
            {
                "open-segment",
                reverse_segment_endpoints
                    ? "segment-target"
                    : "segment-source",
                reverse_segment_endpoints
                    ? "segment-source"
                    : "segment-target",
            },
        },
        "open-segment",
        "focus",
        parabola_domain,
        0,
        graph,
        node_indices);

    MatDomainPolygon2 segment_limiter_outer;
    segment_limiter_outer.push_back({-110, -1});
    segment_limiter_outer.push_back({110, -1});
    segment_limiter_outer.push_back({110, 1800});
    segment_limiter_outer.push_back({-110, 1800});
    const MatDomainPolygonWithHoles2
        segment_limiter_domain(
            segment_limiter_outer);
    append_point_segment_graph(
        {
            {"source-segment-source", -100, 0},
            {"source-segment-target", 100, 0},
            {"segment-focus", 1, 3},
            {"segment-limiter-source", 0, 1},
            {"segment-limiter-target", 0, 5},
        },
        {
            {
                "source-open-segment",
                reverse_segment_endpoints
                    ? "source-segment-target"
                    : "source-segment-source",
                reverse_segment_endpoints
                    ? "source-segment-source"
                    : "source-segment-target",
            },
            {
                "segment-limiter",
                reverse_segment_endpoints
                    ? "segment-limiter-target"
                    : "segment-limiter-source",
                reverse_segment_endpoints
                    ? "segment-limiter-source"
                    : "segment-limiter-target",
            },
        },
        "source-open-segment",
        "segment-focus",
        segment_limiter_domain,
        0,
        graph,
        node_indices);

    return graph;
}

} // namespace

MatExactGraph2 segment_site_generic_graph_spike()
{
    return segment_site_generic_graph_spike_impl(false);
}

MatExactGraph2
segment_site_reversed_source_graph_spike()
{
    return segment_site_generic_graph_spike_impl(true);
}

namespace {

MatDomainPolygonWithHoles2 segment_segment_spike_domain(
    const CORE::BigRat& lower_x,
    const CORE::BigRat& upper_x)
{
    MatDomainPolygon2 outer;
    outer.push_back({lower_x, -1});
    outer.push_back({upper_x, -1});
    outer.push_back({upper_x, 4});
    outer.push_back({lower_x, 4});
    return MatDomainPolygonWithHoles2(outer);
}

MatExactGraph2 segment_segment_graph_spike_impl(
    const bool reverse_segment_endpoints,
    const CORE::BigRat& radius_squared,
    const bool disjoint_features,
    const std::vector<NormalizedPointSource2>& limiter_points,
    const std::vector<NormalizedOpenSegmentSource2>& limiter_segments,
    const std::vector<MatExactPointSiteSource2>& point_limiters,
    const std::vector<MatExactOpenSegmentSource2>& segment_limiters,
    const MatDomainPolygonWithHoles2& domain)
{
    MatExactGraph2 graph{{}, {}, 0, 0};
    std::map<std::string, std::size_t> node_indices;
    std::vector<NormalizedPointSource2> points{
        {
            "lower-left",
            disjoint_features ? -4 : -2,
            0,
        },
        {
            "lower-right",
            disjoint_features ? -2 : 6,
            0,
        },
        {"upper-left", 0, 3},
        {"upper-right", 4, 3},
    };
    points.insert(
        points.end(),
        limiter_points.begin(),
        limiter_points.end());
    std::vector<NormalizedOpenSegmentSource2> segments{
        {
            "lower-segment",
            reverse_segment_endpoints
                ? "lower-right"
                : "lower-left",
            reverse_segment_endpoints
                ? "lower-left"
                : "lower-right",
        },
        {
            "upper-segment",
            reverse_segment_endpoints
                ? "upper-right"
                : "upper-left",
            reverse_segment_endpoints
                ? "upper-left"
                : "upper-right",
        },
    };
    segments.insert(
        segments.end(),
        limiter_segments.begin(),
        limiter_segments.end());
    append_parallel_segment_segment_graph(
        points,
        segments,
        "lower-segment",
        "upper-segment",
        domain,
        radius_squared,
        point_limiters,
        segment_limiters,
        graph,
        node_indices);
    return graph;
}

MatExactGraph2 point_limited_parallel_segment_graph_spike_impl(
    const bool reverse_segment_endpoints,
    const CORE::BigRat& radius_squared)
{
    return segment_segment_graph_spike_impl(
        reverse_segment_endpoints,
        radius_squared,
        false,
        {
            {
                "external-limiter",
                5,
                CORE::BigRat(3, 2),
            },
        },
        {},
        {
            {
                "external-limiter",
                {5, 0},
                {
                    CORE::BigRat(3, 2),
                    0,
                },
                1,
            },
        },
        {},
        segment_segment_spike_domain(-1, 5));
}

MatExactGraph2 segment_limited_parallel_segment_graph_spike_impl(
    const bool reverse_segment_endpoints,
    const CORE::BigRat& radius_squared)
{
    return segment_segment_graph_spike_impl(
        reverse_segment_endpoints,
        radius_squared,
        false,
        {
            {"external-limiter-lower", 5, 1},
            {"external-limiter-upper", 5, 2},
        },
        {
            {
                "external-segment-limiter",
                reverse_segment_endpoints
                    ? "external-limiter-upper"
                    : "external-limiter-lower",
                reverse_segment_endpoints
                    ? "external-limiter-lower"
                    : "external-limiter-upper",
            },
        },
        {},
        {
            canonical_open_segment_source(
                "external-segment-limiter",
                5,
                1,
                5,
                2),
        },
        segment_segment_spike_domain(-1, 5));
}

std::vector<NormalizedPointSource2> rectangle_points()
{
    return {
        {"lower-left", -4, -2},
        {"lower-right", 4, -2},
        {"upper-right", 4, 2},
        {"upper-left", -4, 2},
    };
}

std::vector<NormalizedOpenSegmentSource2> rectangle_segments(
    const bool reverse_segment_endpoints)
{
    return {
        {
            "bottom-segment",
            reverse_segment_endpoints
                ? "lower-right"
                : "lower-left",
            reverse_segment_endpoints
                ? "lower-left"
                : "lower-right",
        },
        {
            "right-segment",
            reverse_segment_endpoints
                ? "upper-right"
                : "lower-right",
            reverse_segment_endpoints
                ? "lower-right"
                : "upper-right",
        },
        {
            "top-segment",
            reverse_segment_endpoints
                ? "upper-left"
                : "upper-right",
            reverse_segment_endpoints
                ? "upper-right"
                : "upper-left",
        },
        {
            "left-segment",
            reverse_segment_endpoints
                ? "lower-left"
                : "upper-left",
            reverse_segment_endpoints
                ? "upper-left"
                : "lower-left",
        },
    };
}

MatDomainPolygonWithHoles2 rectangle_domain()
{
    MatDomainPolygon2 outer;
    outer.push_back({-4, -2});
    outer.push_back({4, -2});
    outer.push_back({4, 2});
    outer.push_back({-4, 2});
    return MatDomainPolygonWithHoles2(outer);
}

std::vector<MatExactOpenSegmentSource2>
exact_segment_sources(
    const std::vector<NormalizedPointSource2>& points,
    const std::vector<NormalizedOpenSegmentSource2>& segments)
{
    std::vector<MatExactOpenSegmentSource2> sources;
    sources.reserve(segments.size());
    for (const NormalizedOpenSegmentSource2& segment :
         segments) {
        sources.push_back(
            exact_segment_source(segment, points));
    }
    return sources;
}

CanonicalNodeAlias2 rectangle_vertex_alias(
    const SegmentSiteVoronoi2::Vertex_handle& vertex,
    const std::vector<GeneratorSite2>& generators,
    const std::vector<NormalizedOpenSegmentSource2>& segments)
{
    std::vector<std::string> feature_ids;
    auto first = vertex->incident_halfedges();
    auto halfedge = first;
    do {
        union_stable_ids(
            feature_ids,
            {
                stable_generator_site_id(
                    halfedge->up()->site(),
                    generators),
                stable_generator_site_id(
                    halfedge->down()->site(),
                    generators),
            });
    } while (++halfedge != first);
    if (feature_ids.size() != 3) {
        throw InvalidSegmentSiteNodeProvenanceError(
            "rectangle adaptor vertex does not have three exact features");
    }

    std::vector<std::string> parent_ids;
    for (const std::string& feature_id : feature_ids) {
        union_stable_ids(
            parent_ids,
            raw_feature_parent_ids(
                feature_id,
                segments));
    }
    return {
        stable_voronoi_node_identity_v1(
            feature_ids),
        std::move(feature_ids),
        std::move(parent_ids),
    };
}

const MatParameterEndpoint2& endpoint_owned_by(
    const std::pair<
        MatParameterEndpoint2,
        MatParameterEndpoint2>& endpoints,
    const std::string& owner_id)
{
    const auto has_owner =
        [&owner_id](
            const MatParameterEndpoint2& endpoint) {
            return std::find(
                       endpoint.provenance_ids.begin(),
                       endpoint.provenance_ids.end(),
                       owner_id)
                != endpoint.provenance_ids.end();
        };
    const bool first_matches =
        has_owner(endpoints.first);
    const bool second_matches =
        has_owner(endpoints.second);
    if (first_matches == second_matches) {
        throw InvalidSegmentSiteNodeProvenanceError(
            "rectangle endpoint does not have one exact adaptor owner");
    }
    return first_matches
        ? endpoints.first
        : endpoints.second;
}

void register_rectangle_vertex_alias(
    const std::string& dual_id,
    const MatParameterEndpoint2& endpoint,
    const CanonicalNodeAlias2& alias,
    const std::string& owner_id,
    std::map<std::string, CanonicalNodeAlias2>& aliases)
{
    if (std::find(
            alias.generator_site_ids.begin(),
            alias.generator_site_ids.end(),
            owner_id)
        == alias.generator_site_ids.end()) {
        throw InvalidSegmentSiteNodeProvenanceError(
            "rectangle endpoint owner is absent from its adaptor vertex");
    }
    const std::string endpoint_node_id =
        stable_endpoint_node_identity_v1(
            dual_id,
            endpoint);
    const auto [existing, inserted] =
        aliases.emplace(
            endpoint_node_id,
            alias);
    if (!inserted
        && (existing->second.node_id
                != alias.node_id
            || existing->second.generator_site_ids
                != alias.generator_site_ids
            || existing->second.parent_site_ids
                != alias.parent_site_ids)) {
        throw InvalidSegmentSiteNodeProvenanceError(
            "rectangle exact endpoint aliases disagree");
    }
}

void register_rectangle_halfedge_aliases(
    const std::string& dual_id,
    const std::pair<
        MatParameterEndpoint2,
        MatParameterEndpoint2>& endpoints,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge,
    const std::vector<GeneratorSite2>& generators,
    const std::vector<NormalizedOpenSegmentSource2>& segments,
    std::map<std::string, CanonicalNodeAlias2>& aliases)
{
    const std::string source_owner_id =
        stable_generator_site_id(
            halfedge->left()->site(),
            generators);
    const std::string target_owner_id =
        stable_generator_site_id(
            halfedge->right()->site(),
            generators);
    const CanonicalNodeAlias2 source_alias =
        rectangle_vertex_alias(
            halfedge->source(),
            generators,
            segments);
    const CanonicalNodeAlias2 target_alias =
        rectangle_vertex_alias(
            halfedge->target(),
            generators,
            segments);
    if (source_alias.node_id == target_alias.node_id) {
        throw InvalidSegmentSiteNodeProvenanceError(
            "rectangle live endpoints share a feature triple; left="
            + source_owner_id
            + "; right="
            + target_owner_id
            + "; node="
            + source_alias.node_id);
    }
    register_rectangle_vertex_alias(
        dual_id,
        endpoint_owned_by(
            endpoints,
            source_owner_id),
        source_alias,
        source_owner_id,
        aliases);
    register_rectangle_vertex_alias(
        dual_id,
        endpoint_owned_by(
            endpoints,
            target_owner_id),
        target_alias,
        target_owner_id,
        aliases);
}

bool is_normalized_segment_endpoint(
    const NormalizedOpenSegmentSource2& segment,
    const std::string& point_id)
{
    return segment.source_point_id == point_id
        || segment.target_point_id == point_id;
}

MatExactGraph2 rectangle_graph_spike_impl(
    const bool reverse_segment_endpoints,
    const CORE::BigRat& radius_squared)
{
    if (radius_squared < 0) {
        throw NegativeClearanceRadiusSquaredError(
            "rectangle graph squared clearance radius is negative");
    }
    const std::vector<NormalizedPointSource2> points =
        rectangle_points();
    const std::vector<NormalizedOpenSegmentSource2> segments =
        rectangle_segments(
            reverse_segment_endpoints);
    const std::vector<GeneratorSite2> generators =
        segment_site_generators(points, segments);
    const std::vector<MatExactOpenSegmentSource2>
        exact_segments =
            exact_segment_sources(
                points,
                segments);
    SegmentSiteDelaunay2 delaunay;
    insert_segment_site_sources(
        points,
        segments,
        delaunay);

    MatExactGraph2 graph{{}, {}, 0, 0};
    graph.matched_generator_sites =
        require_generator_site_bijection(
            delaunay,
            generators);
    SegmentSiteVoronoi2 voronoi(delaunay);
    const MatDomainPolygonWithHoles2 domain =
        rectangle_domain();
    std::map<std::string, std::size_t> node_indices;
    std::map<std::string, CanonicalNodeAlias2> aliases;
    std::set<std::vector<std::string>>
        segment_pairs;
    std::set<std::vector<std::string>>
        incident_pairs;

    for (auto halfedge = voronoi.halfedges_begin();
         halfedge != voronoi.halfedges_end();
         ++halfedge) {
        const MatTraits::Site_2 up =
            halfedge->up()->site();
        const MatTraits::Site_2 down =
            halfedge->down()->site();
        const std::string up_id =
            stable_generator_site_id(
                up,
                generators);
        const std::string down_id =
            stable_generator_site_id(
                down,
                generators);
        const std::vector<std::string> pair =
            ordered_generator_site_ids(
                up_id,
                down_id);
        if (up_id != pair.front()) {
            continue;
        }

        if (up.is_segment() && down.is_segment()) {
            if (!halfedge->has_source()
                || !halfedge->has_target()) {
                throw IncompleteCompositeSegmentChainError(
                    "rectangle S-S halfedge is unbounded");
            }
            MatTraits::Segment_2 representative;
            if (!CGAL::assign(
                    representative,
                    voronoi.dual().primal(
                        halfedge->dual()))) {
                throw UnsupportedCompositeSegmentPrimitiveError(
                    "rectangle S-S dual is not a segment");
            }
            if (!segment_pairs.insert(pair).second) {
                throw IncompleteCompositeSegmentChainError(
                    "rectangle has duplicate canonical S-S cells");
            }

            const NormalizedOpenSegmentSource2&
                first_record =
                    normalized_segment_source(
                        segments,
                        pair[0]);
            const NormalizedOpenSegmentSource2&
                second_record =
                    normalized_segment_source(
                        segments,
                        pair[1]);
            const MatExactOpenSegmentSource2 first =
                exact_segment_source(
                    first_record,
                    points);
            const MatExactOpenSegmentSource2 second =
                exact_segment_source(
                    second_record,
                    points);
            const CORE::BigRat support_determinant =
                first.line_a * second.line_b
                - first.line_b * second.line_a;

            std::string dual_id;
            std::pair<
                MatParameterEndpoint2,
                MatParameterEndpoint2> endpoints;
            std::vector<MatAdmissibleComponent2>
                components;
            if (support_determinant == 0) {
                const ParallelSegmentFeatureCell2
                    feature =
                        parallel_segment_feature_cell(
                            first,
                            first_record,
                            second,
                            second_record,
                            points);
                dual_id =
                    stable_dual_identity_v1(
                        "segment-segment",
                        pair);
                endpoints =
                    bind_parallel_segment_segment_cell_endpoints(
                        feature.primitive,
                        feature.lower,
                        feature.upper,
                        first,
                        second,
                        {},
                        exact_segments,
                        pair,
                        generators,
                        voronoi,
                        halfedge);
                components =
                    clip_bounded_linear_clearance_components(
                        dual_id,
                        feature.primitive,
                        endpoints.first,
                        endpoints.second,
                        parallel_segment_clearance_boundary(
                            feature.primitive,
                            first,
                            second,
                            radius_squared),
                        domain);
            } else {
                const NonparallelSegmentBisectorParameterization2
                    primitive =
                        nonparallel_segment_bisector_parameterization(
                            first,
                            second,
                            representative);
                const NonparallelSegmentFeatureDomain2
                    feature_domain =
                        nonparallel_segment_feature_domain(
                            primitive,
                            points,
                            segments);
                dual_id =
                    nonparallel_original_dual_id(
                        primitive.branch_sign,
                        pair);
                endpoints =
                    bind_nonparallel_segment_segment_cell_endpoints(
                        primitive,
                        feature_domain,
                        first,
                        second,
                        exact_segments,
                        pair,
                        generators,
                        voronoi,
                        halfedge);
                ExactAlgebraicKernel1 kernel;
                const CORE::BigRat live_cell_witness =
                    kernel.bound_between_1_object()(
                        *endpoints.first.parameter,
                        *endpoints.second.parameter);
                if (!nonparallel_segment_domain_contains(
                        domain,
                        primitive,
                        live_cell_witness)) {
                    throw IncompleteCompositeSegmentChainError(
                        "rectangle live S-S open cell is outside its domain");
                }
                components =
                    clip_bounded_nonparallel_segment_clearance_components(
                        dual_id,
                        primitive,
                        feature_domain,
                        endpoints.first,
                        endpoints.second,
                        nonparallel_segment_clearance_boundary(
                            primitive,
                            first,
                            second,
                            radius_squared),
                        domain);
            }
            components =
                one_dimensional_graph_components(
                    components);
            append_exact_graph_components(
                dual_id,
                "LINE",
                pair,
                components,
                graph,
                node_indices);
            register_rectangle_halfedge_aliases(
                dual_id,
                endpoints,
                halfedge,
                generators,
                segments,
                aliases);
            continue;
        }

        if (up.is_point() == down.is_point()) {
            throw UnsupportedCompositeSegmentPrimitiveError(
                "rectangle raw transition is not point-segment");
        }
        const std::string point_id =
            up.is_point() ? up_id : down_id;
        const std::string segment_id =
            up.is_segment() ? up_id : down_id;
        if (!is_normalized_segment_endpoint(
                normalized_segment_source(
                    segments,
                    segment_id),
                point_id)) {
            throw IncompleteCompositeSegmentChainError(
                "rectangle has a nonincident point-segment transition");
        }
        MatTraits::Ray_2 ray;
        if (halfedge->has_source()
                == halfedge->has_target()
            || !CGAL::assign(
                ray,
                voronoi.dual().primal(
                    halfedge->dual()))) {
            throw UnsupportedCompositeSegmentPrimitiveError(
                "rectangle incident point-segment transition is not a ray");
        }
        if (!incident_pairs.insert(pair).second) {
            throw IncompleteCompositeSegmentChainError(
                "rectangle has duplicate incident point-segment transitions");
        }
        ++graph.rejected_incident_transitions;
    }

    const std::set<std::vector<std::string>>
        expected_segment_pairs{
            {"bottom-segment", "left-segment"},
            {"bottom-segment", "right-segment"},
            {"bottom-segment", "top-segment"},
            {"left-segment", "top-segment"},
            {"right-segment", "top-segment"},
        };
    const std::set<std::vector<std::string>>
        expected_incident_pairs{
            {"bottom-segment", "lower-left"},
            {"bottom-segment", "lower-right"},
            {"left-segment", "lower-left"},
            {"left-segment", "upper-left"},
            {"lower-right", "right-segment"},
            {"right-segment", "upper-right"},
            {"top-segment", "upper-left"},
            {"top-segment", "upper-right"},
        };
    if (segment_pairs != expected_segment_pairs
        || incident_pairs != expected_incident_pairs
        || graph.rejected_incident_transitions
            != expected_incident_pairs.size()) {
        throw IncompleteCompositeSegmentChainError(
            "rectangle adaptor traversal is incomplete");
    }

    for (const MatExactGraphEdge2& edge : graph.edges) {
        if (edge.source_node_id == edge.target_node_id) {
            throw InvalidSegmentSiteNodeProvenanceError(
                "rectangle graph emitted a zero-dimensional edge; dual="
                + edge.original_dual_id);
        }
        const auto source_alias =
            aliases.find(edge.source_node_id);
        const auto target_alias =
            aliases.find(edge.target_node_id);
        if (source_alias != aliases.end()
            && target_alias != aliases.end()
            && source_alias->second.node_id
                == target_alias->second.node_id) {
            throw InvalidSegmentSiteNodeProvenanceError(
                std::string(
                    "rectangle endpoint aliases collapse; reversed=")
                + (reverse_segment_endpoints
                       ? "true"
                       : "false")
                + "; dual="
                + edge.original_dual_id
                + "; node="
                + source_alias->second.node_id);
        }
    }
    canonicalize_original_nodes(
        graph,
        aliases);
    const auto collapsed_edge =
        std::find_if(
            graph.edges.begin(),
            graph.edges.end(),
            [](const MatExactGraphEdge2& edge) {
                return edge.source_node_id
                    == edge.target_node_id;
            });
    if (collapsed_edge != graph.edges.end()) {
        throw InvalidSegmentSiteNodeProvenanceError(
            "rectangle graph collapsed dual "
            + collapsed_edge->original_dual_id
            + " at "
            + collapsed_edge->source_node_id);
    }
    std::sort(
        graph.nodes.begin(),
        graph.nodes.end(),
        [](const MatExactGraphNode2& lhs,
           const MatExactGraphNode2& rhs) {
            return lhs.node_id < rhs.node_id;
        });
    std::sort(
        graph.edges.begin(),
        graph.edges.end(),
        [](const MatExactGraphEdge2& lhs,
           const MatExactGraphEdge2& rhs) {
            return lhs.edge_id < rhs.edge_id;
        });
    return graph;
}

MatExactGraph2 rectangle_central_parallel_graph_spike_impl(
    const bool reverse_segment_endpoints,
    const CORE::BigRat& radius_squared)
{
    const std::vector<NormalizedPointSource2> points =
        rectangle_points();
    const std::vector<NormalizedOpenSegmentSource2> segments =
        rectangle_segments(
            reverse_segment_endpoints);
    MatExactGraph2 graph{{}, {}, 0, 0};
    std::map<std::string, std::size_t> node_indices;
    append_parallel_segment_segment_graph(
        points,
        segments,
        "bottom-segment",
        "top-segment",
        rectangle_domain(),
        radius_squared,
        {},
        {
            exact_segment_source(
                normalized_segment_source(
                    segments,
                    "left-segment"),
                points),
            exact_segment_source(
                normalized_segment_source(
                    segments,
                    "right-segment"),
                points),
        },
        graph,
        node_indices);
    return graph;
}

MatExactGraph2 rectangle_lower_left_graph_spike_impl(
    const bool reverse_segment_endpoints,
    const CORE::BigRat& radius_squared)
{
    if (radius_squared < 0)
    {
        throw NegativeClearanceRadiusSquaredError(
            "rectangle corner squared clearance radius is negative");
    }
    const std::vector<NormalizedPointSource2> points =
        rectangle_points();
    const std::vector<NormalizedOpenSegmentSource2> segments =
        rectangle_segments(
            reverse_segment_endpoints);
    const std::vector<std::string> parent_site_ids{
        "bottom-segment",
        "left-segment",
    };
    const MatExactOpenSegmentSource2 first_segment =
        exact_segment_source(
            normalized_segment_source(
                segments,
                parent_site_ids[0]),
            points);
    const MatExactOpenSegmentSource2 second_segment =
        exact_segment_source(
            normalized_segment_source(
                segments,
                parent_site_ids[1]),
            points);
    const std::vector<GeneratorSite2> generators =
        segment_site_generators(
            points,
            segments);
    SegmentSiteDelaunay2 delaunay;
    insert_segment_site_sources(
        points,
        segments,
        delaunay);
    MatExactGraph2 graph{{}, {}, 0, 0};
    graph.matched_generator_sites =
        require_generator_site_bijection(
            delaunay,
            generators);
    SegmentSiteVoronoi2 voronoi(delaunay);
    std::vector<SegmentSiteVoronoi2::Halfedge_handle>
        matching_halfedges;
    for (auto halfedge = voronoi.halfedges_begin();
         halfedge != voronoi.halfedges_end();
         ++halfedge)
    {
        const MatTraits::Site_2 up =
            halfedge->up()->site();
        const MatTraits::Site_2 down =
            halfedge->down()->site();
        if (!up.is_segment()
            || !down.is_segment())
        {
            continue;
        }
        const std::string up_id =
            stable_generator_site_id(
                up,
                generators);
        const std::string down_id =
            stable_generator_site_id(
                down,
                generators);
        if (ordered_generator_site_ids(
                up_id,
                down_id)
                == parent_site_ids
            && up_id == parent_site_ids.front())
        {
            matching_halfedges.push_back(
                halfedge);
        }
    }
    if (matching_halfedges.size() != 1)
    {
        throw IncompleteCompositeSegmentChainError(
            "rectangle lower-left corner has no unique S-S halfedge");
    }
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge =
        matching_halfedges.front();
    MatTraits::Segment_2 representative;
    if (!CGAL::assign(
            representative,
            voronoi.dual().primal(
                halfedge->dual())))
    {
        throw UnsupportedCompositeSegmentPrimitiveError(
            "rectangle lower-left S-S dual is not a segment");
    }
    const NonparallelSegmentBisectorParameterization2 primitive =
        nonparallel_segment_bisector_parameterization(
            first_segment,
            second_segment,
            representative);
    const NonparallelSegmentFeatureDomain2 feature_domain =
        nonparallel_segment_feature_domain(
            primitive,
            points,
            segments);
    const auto endpoints =
        bind_nonparallel_segment_segment_cell_endpoints(
            primitive,
            feature_domain,
            first_segment,
            second_segment,
            {
                exact_segment_source(
                    normalized_segment_source(
                        segments,
                        "top-segment"),
                    points),
            },
            parent_site_ids,
            generators,
            voronoi,
            halfedge);
    const std::string dual_id =
        nonparallel_original_dual_id(
            primitive.branch_sign,
            parent_site_ids);
    const std::vector<MatAdmissibleComponent2> components =
        clip_bounded_nonparallel_segment_clearance_components(
            dual_id,
            primitive,
            feature_domain,
            endpoints.first,
            endpoints.second,
            nonparallel_segment_clearance_boundary(
                primitive,
                first_segment,
                second_segment,
                radius_squared),
            rectangle_domain());
    std::map<std::string, std::size_t> node_indices;
    append_exact_graph_components(
        dual_id,
        "LINE",
        parent_site_ids,
        components,
        graph,
        node_indices);
    return graph;
}

MatExactGraph2 nonparallel_segment_graph_spike_impl(
    const bool reverse_segment_endpoints,
    const CORE::BigRat& radius_squared)
{
    MatDomainPolygon2 outer;
    outer.push_back({-100, -100});
    outer.push_back({100, -100});
    outer.push_back({100, 100});
    outer.push_back({-100, 100});
    MatExactGraph2 graph{{}, {}, 0, 0};
    std::map<std::string, std::size_t> node_indices;
    append_nonparallel_segment_segment_graph(
        {
            {"lower-left", -20, 0},
            {"lower-right", 8, 0},
            {"diagonal-lower", 5, -4},
            {"diagonal-upper", 20, 11},
        },
        {
            {
                "lower-segment",
                reverse_segment_endpoints
                    ? "lower-right"
                    : "lower-left",
                reverse_segment_endpoints
                    ? "lower-left"
                    : "lower-right",
            },
            {
                "diagonal-segment",
                reverse_segment_endpoints
                    ? "diagonal-upper"
                    : "diagonal-lower",
                reverse_segment_endpoints
                    ? "diagonal-lower"
                    : "diagonal-upper",
            },
        },
        "lower-segment",
        "diagonal-segment",
        MatDomainPolygonWithHoles2(outer),
        radius_squared,
        graph,
        node_indices);
    return graph;
}

} // namespace

MatExactGraph2
segment_site_segment_segment_graph_spike()
{
    return segment_site_segment_segment_graph_spike(
        CORE::BigRat(9, 4));
}

MatExactGraph2
segment_site_segment_segment_graph_spike(
    const CORE::BigRat& radius_squared)
{
    return segment_segment_graph_spike_impl(
        false,
        radius_squared,
        false,
        {{"limiter", 6, 1}},
        {},
        {},
        {},
        segment_segment_spike_domain(-1, 5));
}

MatExactGraph2
segment_site_reversed_segment_segment_graph_spike()
{
    return segment_segment_graph_spike_impl(
        true,
        CORE::BigRat(9, 4),
        false,
        {{"limiter", 6, 1}},
        {},
        {},
        {},
        segment_segment_spike_domain(-1, 5));
}

MatExactGraph2
segment_site_point_limited_parallel_segment_graph_spike()
{
    return segment_site_point_limited_parallel_segment_graph_spike(
        CORE::BigRat(9, 4));
}

MatExactGraph2
segment_site_point_limited_parallel_segment_graph_spike(
    const CORE::BigRat& radius_squared)
{
    return point_limited_parallel_segment_graph_spike_impl(
        false,
        radius_squared);
}

MatExactGraph2
segment_site_reversed_point_limited_parallel_segment_graph_spike()
{
    return point_limited_parallel_segment_graph_spike_impl(
        true,
        CORE::BigRat(9, 4));
}

MatExactGraph2
segment_site_segment_limited_parallel_segment_graph_spike()
{
    return segment_site_segment_limited_parallel_segment_graph_spike(
        CORE::BigRat(9, 4));
}

MatExactGraph2
segment_site_segment_limited_parallel_segment_graph_spike(
    const CORE::BigRat& radius_squared)
{
    return segment_limited_parallel_segment_graph_spike_impl(
        false,
        radius_squared);
}

MatExactGraph2
segment_site_reversed_segment_limited_parallel_segment_graph_spike()
{
    return segment_limited_parallel_segment_graph_spike_impl(
        true,
        CORE::BigRat(9, 4));
}

MatExactGraph2
segment_site_rectangle_central_parallel_graph_spike()
{
    return segment_site_rectangle_central_parallel_graph_spike(4);
}

MatExactGraph2
segment_site_rectangle_central_parallel_graph_spike(
    const CORE::BigRat& radius_squared)
{
    return rectangle_central_parallel_graph_spike_impl(
        false,
        radius_squared);
}

MatExactGraph2
segment_site_reversed_rectangle_central_parallel_graph_spike()
{
    return rectangle_central_parallel_graph_spike_impl(
        true,
        4);
}

MatExactGraph2
segment_site_rectangle_lower_left_graph_spike()
{
    return segment_site_rectangle_lower_left_graph_spike(0);
}

MatExactGraph2
segment_site_rectangle_lower_left_graph_spike(
    const CORE::BigRat& radius_squared)
{
    return rectangle_lower_left_graph_spike_impl(
        false,
        radius_squared);
}

MatExactGraph2
segment_site_reversed_rectangle_lower_left_graph_spike()
{
    return rectangle_lower_left_graph_spike_impl(
        true,
        0);
}

MatExactGraph2
segment_site_rectangle_graph_spike()
{
    return segment_site_rectangle_graph_spike(0);
}

MatExactGraph2
segment_site_rectangle_graph_spike(
    const CORE::BigRat& radius_squared)
{
    return rectangle_graph_spike_impl(
        false,
        radius_squared);
}

MatExactGraph2
segment_site_reversed_rectangle_graph_spike()
{
    return rectangle_graph_spike_impl(
        true,
        0);
}

MatExactGraph2
segment_site_nonparallel_segment_segment_graph_spike()
{
    return segment_site_nonparallel_segment_segment_graph_spike(0);
}

MatExactGraph2
segment_site_nonparallel_segment_segment_graph_spike(
    const CORE::BigRat& radius_squared)
{
    return nonparallel_segment_graph_spike_impl(
        false,
        radius_squared);
}

MatExactGraph2
segment_site_reversed_nonparallel_segment_segment_graph_spike()
{
    return nonparallel_segment_graph_spike_impl(true, 0);
}

MatExactGraph2
segment_site_disjoint_parallel_segment_graph_spike()
{
    return segment_segment_graph_spike_impl(
        false,
        CORE::BigRat(9, 4),
        true,
        {{"limiter", 6, 1}},
        {},
        {},
        {},
        segment_segment_spike_domain(-1, 5));
}

MatExactGraph2
segment_site_external_limited_parallel_segment_graph_spike()
{
    return segment_segment_graph_spike_impl(
        false,
        CORE::BigRat(9, 4),
        false,
        {
            {
                "external-limiter",
                5,
                CORE::BigRat(3, 2),
            },
        },
        {},
        {},
        {},
        segment_segment_spike_domain(-1, 5));
}

MatExactGraph2
segment_site_domain_coincident_parallel_segment_graph_spike()
{
    return segment_segment_graph_spike_impl(
        false,
        CORE::BigRat(9, 4),
        false,
        {{"limiter", 6, 1}},
        {},
        {},
        {},
        segment_segment_spike_domain(0, 4));
}
