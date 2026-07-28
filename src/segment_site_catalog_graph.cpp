#include "segment_site_catalog_graph.h"

#include "segment_site_catalog.h"
#include "segment_site_graph_emission.h"
#include "segment_site_rational_sources.h"
#include "segment_site_voronoi.h"

#include <algorithm>
#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <CGAL/Object.h>

namespace {

struct ParallelSegmentFeatureCell2 {
    RationalPrimitiveParameterization2 primitive;
    RationalDomainRoot2 lower;
    RationalDomainRoot2 upper;
};

struct CanonicalNodeAlias2 {
    std::string node_id;
    std::vector<std::string> generator_site_ids;
    std::vector<std::string> parent_site_ids;
};

const CanonicalMatRationalPointSource2&
point_source(
    const CanonicalMatRationalSources2& sources,
    const std::string& stable_id)
{
    const auto found = std::lower_bound(
        sources.points().begin(),
        sources.points().end(),
        stable_id,
        [](const CanonicalMatRationalPointSource2&
               point,
           const std::string& identity) {
            return point.stable_site_id
                < identity;
        });
    if (found == sources.points().end()
        || found->stable_site_id != stable_id) {
        throw UnknownCanonicalMatGraphSourceError(
            "canonical MAT graph point has no rational source");
    }
    return *found;
}

const CanonicalMatRationalOpenSegmentSource2&
segment_source(
    const CanonicalMatRationalSources2& sources,
    const std::string& stable_id)
{
    const auto found = std::lower_bound(
        sources.segments().begin(),
        sources.segments().end(),
        stable_id,
        [](const CanonicalMatRationalOpenSegmentSource2&
               segment,
           const std::string& identity) {
            return segment.stable_site_id
                < identity;
        });
    if (found == sources.segments().end()
        || found->stable_site_id != stable_id) {
        throw UnknownCanonicalMatGraphSourceError(
            "canonical MAT graph segment has no rational source");
    }
    return *found;
}

std::vector<MatExactOpenSegmentSource2>
segment_supports(
    const CanonicalMatRationalSources2& sources)
{
    std::vector<MatExactOpenSegmentSource2> supports;
    supports.reserve(sources.segments().size());
    for (const auto& segment :
         sources.segments()) {
        supports.push_back(segment.support);
    }
    return supports;
}

CORE::BigRat support_determinant(
    const MatExactOpenSegmentSource2& first,
    const MatExactOpenSegmentSource2& second)
{
    return first.line_a * second.line_b
        - first.line_b * second.line_a;
}

void require_rectangle_input(
    const CanonicalReachInput2& input,
    const CanonicalMatSiteCatalog2& catalog,
    const CanonicalMatRationalSources2& sources)
{
    if (!input.holes.empty()
        || input.outer.points.size() != 4
        || catalog.sites().size() != 8
        || sources.points().size() != 4
        || sources.segments().size() != 4) {
        throw UnsupportedCanonicalMatRectangleGraphError(
            "canonical rectangle MAT graph requires one four-edge outer ring");
    }
    const auto& first =
        sources.segments()[0].support;
    const auto& second =
        sources.segments()[1].support;
    const auto& third =
        sources.segments()[2].support;
    const auto& fourth =
        sources.segments()[3].support;
    if (support_determinant(first, third) != 0
        || support_determinant(second, fourth)
            != 0
        || support_determinant(first, second)
            == 0
        || first.line_a * second.line_a
                + first.line_b * second.line_b
            != 0) {
        throw UnsupportedCanonicalMatRectangleGraphError(
            "canonical four-edge ring is not an exact rectangle");
    }
}

MatDomainPolygonWithHoles2 exact_domain(
    const CanonicalReachInput2& input)
{
    const auto exact_ring =
        [](const CanonicalReachRing2& ring) {
            MatDomainPolygon2 polygon;
            for (const ReachKernelPoint& point :
                 ring.points) {
                polygon.push_back(
                    {
                        exact_mat_input_rational(
                            point.x()),
                        exact_mat_input_rational(
                            point.y()),
                    });
            }
            return polygon;
        };
    MatDomainPolygon2 outer =
        exact_ring(input.outer);
    std::vector<MatDomainPolygon2> holes;
    holes.reserve(input.holes.size());
    for (const CanonicalReachRing2& hole :
         input.holes) {
        holes.push_back(exact_ring(hole));
    }
    return MatDomainPolygonWithHoles2(
        outer,
        holes.begin(),
        holes.end());
}

ParallelSegmentFeatureCell2
parallel_segment_feature_cell(
    const CanonicalMatRationalOpenSegmentSource2&
        first,
    const CanonicalMatRationalOpenSegmentSource2&
        second,
    const CanonicalMatRationalSources2& sources)
{
    RationalPrimitiveParameterization2 primitive =
        parallel_segment_bisector_parameterization(
            first.support,
            second.support);
    struct EndpointParameter2 {
        CORE::BigRat parameter;
        std::string point_id;
        std::string segment_id;
    };
    const auto endpoint =
        [&primitive, &sources](
            const std::string& point_id,
            const std::string& segment_id) {
            const auto& point =
                point_source(
                    sources,
                    point_id);
            return EndpointParameter2{
                parallel_segment_tangent_parameter(
                    primitive,
                    point.x,
                    point.y),
                point.stable_site_id,
                segment_id,
            };
        };
    std::vector<EndpointParameter2> endpoints{
        endpoint(
            first.source_point_id,
            first.stable_site_id),
        endpoint(
            first.target_point_id,
            first.stable_site_id),
        endpoint(
            second.source_point_id,
            second.stable_site_id),
        endpoint(
            second.target_point_id,
            second.stable_site_id),
    };
    const auto segment_bounds =
        [&endpoints](
            const std::string& segment_id) {
            std::vector<CORE::BigRat> values;
            for (const EndpointParameter2& endpoint :
                 endpoints) {
                if (endpoint.segment_id
                    == segment_id) {
                    values.push_back(
                        endpoint.parameter);
                }
            }
            if (values.size() != 2) {
                throw IncompleteCanonicalMatRectangleGraphError(
                    "parallel canonical segment has incomplete endpoints");
            }
            if (values[1] < values[0]) {
                std::swap(
                    values[0],
                    values[1]);
            }
            return std::pair<
                CORE::BigRat,
                CORE::BigRat>{
                values[0],
                values[1],
            };
        };
    const auto first_bounds =
        segment_bounds(first.stable_site_id);
    const auto second_bounds =
        segment_bounds(second.stable_site_id);
    const CORE::BigRat lower =
        std::max(
            first_bounds.first,
            second_bounds.first);
    const CORE::BigRat upper =
        std::min(
            first_bounds.second,
            second_bounds.second);
    if (lower >= upper) {
        throw EmptyParallelSegmentFeatureDomainError(
            "canonical parallel segments have no positive overlap");
    }
    const auto provenance_at =
        [&endpoints](
            const CORE::BigRat& parameter) {
            std::vector<std::string> provenance;
            for (const EndpointParameter2& endpoint :
                 endpoints) {
                if (endpoint.parameter
                    == parameter) {
                    provenance.push_back(
                        endpoint.point_id);
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

NonparallelSegmentFeatureDomain2
nonparallel_segment_feature_domain(
    const NonparallelSegmentBisectorParameterization2&
        primitive,
    const CanonicalMatRationalSources2& sources)
{
    const auto& first =
        segment_source(
            sources,
            primitive.first_segment_id);
    const auto& second =
        segment_source(
            sources,
            primitive.second_segment_id);
    const auto endpoint =
        [&primitive, &sources](
            const CanonicalMatRationalOpenSegmentSource2&
                segment,
            const std::string& point_id) {
            const auto& point =
                point_source(
                    sources,
                    point_id);
            return MatQuadraticFieldDomainBoundary2{
                nonparallel_segment_tangent_parameter(
                    primitive,
                    segment.support,
                    point.x,
                    point.y),
                {point.stable_site_id},
            };
        };
    return intersect_nonparallel_segment_feature_domains(
        endpoint(
            first,
            first.source_point_id),
        endpoint(
            first,
            first.target_point_id),
        endpoint(
            second,
            second.source_point_id),
        endpoint(
            second,
            second.target_point_id),
        primitive.radicand);
}

std::string nonparallel_dual_id(
    const NonparallelSegmentBisectorParameterization2&
        primitive,
    const std::vector<std::string>& generator_ids)
{
    if (primitive.branch_sign != -1
        && primitive.branch_sign != 1) {
        throw IncompleteCanonicalMatRectangleGraphError(
            "canonical nonparallel branch sign is invalid");
    }
    return stable_dual_identity_v1(
        primitive.branch_sign < 0
            ? "segment-segment/branch-negative"
            : "segment-segment/branch-positive",
        generator_ids);
}

std::map<
    std::string,
    std::vector<std::string>>
parent_ids_by_feature(
    const CanonicalMatSiteCatalog2& catalog)
{
    std::map<
        std::string,
        std::vector<std::string>>
        parents;
    for (const CanonicalMatSite2& site :
         catalog.sites()) {
        parents.emplace(
            site.stable_id(),
            std::vector<std::string>{});
    }
    for (const CanonicalMatSite2& site :
         catalog.sites()) {
        if (site.site().is_point()) {
            continue;
        }
        auto segment = parents.find(
            site.stable_id());
        if (segment == parents.end()
            || !site.endpoint_site_ids()
                    .has_value()) {
            throw InvalidCanonicalMatRectangleNodeError(
                "canonical segment parent ownership is incomplete");
        }
        segment->second.push_back(
            site.stable_id());
        for (const std::string& endpoint_id :
             *site.endpoint_site_ids()) {
            auto endpoint =
                parents.find(endpoint_id);
            if (endpoint == parents.end()) {
                throw InvalidCanonicalMatRectangleNodeError(
                    "canonical segment endpoint is outside the catalog");
            }
            endpoint->second.push_back(
                site.stable_id());
        }
    }
    for (auto& [feature_id, feature_parents] :
         parents) {
        static_cast<void>(feature_id);
        std::sort(
            feature_parents.begin(),
            feature_parents.end());
        feature_parents.erase(
            std::unique(
                feature_parents.begin(),
                feature_parents.end()),
            feature_parents.end());
        if (feature_parents.empty()
            || feature_parents.size() > 2) {
            throw InvalidCanonicalMatRectangleNodeError(
                "canonical feature has invalid parent incidence");
        }
    }
    return parents;
}

CanonicalNodeAlias2 canonical_vertex_alias(
    const SegmentSiteVoronoi2::Vertex_handle&
        vertex,
    const CanonicalMatSiteGeometryIndex2&
        site_index,
    const std::map<
        std::string,
        std::vector<std::string>>&
        parent_ids)
{
    std::vector<std::string> feature_ids;
    auto first = vertex->incident_halfedges();
    auto halfedge = first;
    do {
        union_stable_ids(
            feature_ids,
            {
                site_index.stable_id(
                    halfedge->up()->site()),
                site_index.stable_id(
                    halfedge->down()->site()),
            });
    } while (++halfedge != first);
    if (feature_ids.size() != 3) {
        throw InvalidCanonicalMatRectangleNodeError(
            "canonical rectangle vertex does not have three exact features");
    }

    std::vector<std::string> parent_site_ids;
    for (const std::string& feature_id :
         feature_ids) {
        const auto found =
            parent_ids.find(feature_id);
        if (found == parent_ids.end()) {
            throw InvalidCanonicalMatRectangleNodeError(
                "canonical vertex feature has no parent record");
        }
        union_stable_ids(
            parent_site_ids,
            found->second);
    }
    return {
        stable_voronoi_node_identity_v1(
            feature_ids),
        std::move(feature_ids),
        std::move(parent_site_ids),
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
    const bool first = has_owner(endpoints.first);
    const bool second = has_owner(endpoints.second);
    if (first == second) {
        throw InvalidCanonicalMatRectangleNodeError(
            "canonical endpoint does not have one exact owner");
    }
    return first
        ? endpoints.first
        : endpoints.second;
}

void register_vertex_alias(
    const std::string& dual_id,
    const MatParameterEndpoint2& endpoint,
    const CanonicalNodeAlias2& alias,
    const std::string& owner_id,
    std::map<std::string, CanonicalNodeAlias2>&
        aliases)
{
    if (std::find(
            alias.generator_site_ids.begin(),
            alias.generator_site_ids.end(),
            owner_id)
        == alias.generator_site_ids.end()) {
        throw InvalidCanonicalMatRectangleNodeError(
            "canonical endpoint owner is absent from its adaptor vertex");
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
        throw InvalidCanonicalMatRectangleNodeError(
            "canonical endpoint aliases disagree");
    }
}

void register_halfedge_aliases(
    const std::string& dual_id,
    const std::pair<
        MatParameterEndpoint2,
        MatParameterEndpoint2>& endpoints,
    const SegmentSiteVoronoi2::Halfedge_handle&
        halfedge,
    const CanonicalMatSiteGeometryIndex2&
        site_index,
    const std::map<
        std::string,
        std::vector<std::string>>&
        parent_ids,
    std::map<std::string, CanonicalNodeAlias2>&
        aliases)
{
    const std::string source_owner_id =
        site_index.stable_id(
            halfedge->left()->site());
    const std::string target_owner_id =
        site_index.stable_id(
            halfedge->right()->site());
    const CanonicalNodeAlias2 source_alias =
        canonical_vertex_alias(
            halfedge->source(),
            site_index,
            parent_ids);
    const CanonicalNodeAlias2 target_alias =
        canonical_vertex_alias(
            halfedge->target(),
            site_index,
            parent_ids);
    if (source_alias.node_id
        == target_alias.node_id) {
        throw InvalidCanonicalMatRectangleNodeError(
            "canonical halfedge endpoints share one feature triple");
    }
    register_vertex_alias(
        dual_id,
        endpoint_owned_by(
            endpoints,
            source_owner_id),
        source_alias,
        source_owner_id,
        aliases);
    register_vertex_alias(
        dual_id,
        endpoint_owned_by(
            endpoints,
            target_owner_id),
        target_alias,
        target_owner_id,
        aliases);
}

void canonicalize_nodes(
    MatExactGraph2& graph,
    const std::map<std::string, CanonicalNodeAlias2>&
        aliases)
{
    for (MatExactGraphEdge2& edge :
         graph.edges) {
        const auto source =
            aliases.find(edge.source_node_id);
        if (source != aliases.end()) {
            edge.source_node_id =
                source->second.node_id;
        }
        const auto target =
            aliases.find(edge.target_node_id);
        if (target != aliases.end()) {
            edge.target_node_id =
                target->second.node_id;
        }
    }

    std::vector<MatExactGraphNode2> canonical;
    std::map<std::string, std::size_t> indices;
    for (MatExactGraphNode2& node :
         graph.nodes) {
        const auto alias =
            aliases.find(node.node_id);
        if (alias != aliases.end()) {
            node.node_id = alias->second.node_id;
            union_stable_ids(
                node.generator_site_ids,
                alias->second
                    .generator_site_ids);
            union_stable_ids(
                node.parent_site_ids,
                alias->second.parent_site_ids);
        }
        const auto existing =
            indices.find(node.node_id);
        if (existing == indices.end()) {
            indices.emplace(
                node.node_id,
                canonical.size());
            canonical.push_back(
                std::move(node));
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

bool is_segment_endpoint(
    const CanonicalMatRationalOpenSegmentSource2&
        segment,
    const std::string& point_id)
{
    return segment.source_point_id == point_id
        || segment.target_point_id == point_id;
}

} // namespace

MatExactGraph2 canonical_rectangle_mat_graph(
    const CanonicalReachInput2& input,
    const CORE::BigRat& radius_squared)
{
    if (radius_squared < 0) {
        throw NegativeClearanceRadiusSquaredError(
            "canonical rectangle MAT squared clearance radius is negative");
    }
    const CanonicalMatSiteCatalog2 catalog =
        canonical_mat_site_catalog(input);
    const CanonicalMatRationalSources2 sources =
        CanonicalMatRationalSources2::build(
            catalog);
    require_rectangle_input(
        input,
        catalog,
        sources);
    CanonicalMatDelaunaySource2 delaunay =
        CanonicalMatDelaunaySource2::build(
            catalog);
    const CanonicalMatVoronoiSource2 source =
        CanonicalMatVoronoiSource2::build(
            std::move(delaunay));
    const MatDomainPolygonWithHoles2 domain =
        exact_domain(input);
    const std::vector<MatExactOpenSegmentSource2>
        limiters =
            segment_supports(sources);
    const auto parent_ids =
        parent_ids_by_feature(catalog);

    MatExactGraph2 graph{
        {},
        {},
        0,
        source.site_index().size(),
    };
    std::map<std::string, std::size_t>
        node_indices;
    std::map<std::string, CanonicalNodeAlias2>
        aliases;
    std::set<std::vector<std::string>>
        segment_pairs;
    std::set<std::vector<std::string>>
        incident_pairs;

    for (auto halfedge =
             source.voronoi().halfedges_begin();
         halfedge
         != source.voronoi().halfedges_end();
         ++halfedge) {
        const MatTraits::Site_2 up =
            halfedge->up()->site();
        const MatTraits::Site_2 down =
            halfedge->down()->site();
        const std::string up_id =
            source.site_index().stable_id(up);
        const std::string down_id =
            source.site_index().stable_id(down);
        const std::vector<std::string> pair =
            ordered_generator_site_ids(
                up_id,
                down_id);
        if (up_id != pair.front()) {
            continue;
        }

        if (up.is_segment()
            && down.is_segment()) {
            if (!halfedge->has_source()
                || !halfedge->has_target()) {
                throw IncompleteCanonicalMatRectangleGraphError(
                    "canonical rectangle S-S halfedge is unbounded");
            }
            MatTraits::Segment_2 representative;
            if (!CGAL::assign(
                    representative,
                    source.voronoi()
                        .dual()
                        .primal(
                            halfedge->dual()))) {
                throw UnsupportedCanonicalMatRectangleGraphError(
                    "canonical rectangle S-S dual is not a segment");
            }
            if (!segment_pairs.insert(pair).second) {
                throw IncompleteCanonicalMatRectangleGraphError(
                    "canonical rectangle has duplicate S-S cells");
            }
            const auto& first =
                segment_source(
                    sources,
                    pair[0]);
            const auto& second =
                segment_source(
                    sources,
                    pair[1]);
            const CORE::BigRat determinant =
                support_determinant(
                    first.support,
                    second.support);

            std::string dual_id;
            std::pair<
                MatParameterEndpoint2,
                MatParameterEndpoint2> endpoints;
            std::vector<MatAdmissibleComponent2>
                components;
            if (determinant == 0) {
                const ParallelSegmentFeatureCell2
                    feature =
                        parallel_segment_feature_cell(
                            first,
                            second,
                            sources);
                dual_id = stable_dual_identity_v1(
                    "segment-segment",
                    pair);
                endpoints =
                    bind_parallel_segment_segment_cell_endpoints(
                        feature.primitive,
                        feature.lower,
                        feature.upper,
                        first.support,
                        second.support,
                        {},
                        limiters,
                        pair,
                        source.site_index(),
                        source.voronoi(),
                        halfedge);
                components =
                    clip_bounded_linear_clearance_components(
                        dual_id,
                        feature.primitive,
                        endpoints.first,
                        endpoints.second,
                        parallel_segment_clearance_boundary(
                            feature.primitive,
                            first.support,
                            second.support,
                            radius_squared),
                        domain);
            } else {
                const NonparallelSegmentBisectorParameterization2
                    primitive =
                        nonparallel_segment_bisector_parameterization(
                            first.support,
                            second.support,
                            representative);
                const NonparallelSegmentFeatureDomain2
                    feature_domain =
                        nonparallel_segment_feature_domain(
                            primitive,
                            sources);
                dual_id =
                    nonparallel_dual_id(
                        primitive,
                        pair);
                endpoints =
                    bind_nonparallel_segment_segment_cell_endpoints(
                        primitive,
                        feature_domain,
                        first.support,
                        second.support,
                        limiters,
                        pair,
                        source.site_index(),
                        source.voronoi(),
                        halfedge);
                ExactAlgebraicKernel1 kernel;
                const CORE::BigRat witness =
                    kernel.bound_between_1_object()(
                        *endpoints.first.parameter,
                        *endpoints.second.parameter);
                if (!nonparallel_segment_domain_contains(
                        domain,
                        primitive,
                        witness)) {
                    throw IncompleteCanonicalMatRectangleGraphError(
                        "canonical rectangle S-S open cell is outside its domain");
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
                            first.support,
                            second.support,
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
            register_halfedge_aliases(
                dual_id,
                endpoints,
                halfedge,
                source.site_index(),
                parent_ids,
                aliases);
            continue;
        }

        if (up.is_point() == down.is_point()) {
            throw UnsupportedCanonicalMatRectangleGraphError(
                "canonical rectangle transition is not point-segment");
        }
        const std::string point_id =
            up.is_point() ? up_id : down_id;
        const std::string segment_id =
            up.is_segment() ? up_id : down_id;
        if (!is_segment_endpoint(
                segment_source(
                    sources,
                    segment_id),
                point_id)) {
            throw UnsupportedCanonicalMatRectangleGraphError(
                "canonical rectangle has a nonincident P-S transition");
        }
        MatTraits::Ray_2 ray;
        if (halfedge->has_source()
                == halfedge->has_target()
            || !CGAL::assign(
                ray,
                source.voronoi()
                    .dual()
                    .primal(
                        halfedge->dual()))) {
            throw UnsupportedCanonicalMatRectangleGraphError(
                "canonical incident P-S transition is not a ray");
        }
        if (!incident_pairs.insert(pair).second) {
            throw IncompleteCanonicalMatRectangleGraphError(
                "canonical rectangle has duplicate incident transitions");
        }
        ++graph.rejected_incident_transitions;
    }

    if (segment_pairs.size() != 5
        || incident_pairs.size() != 8
        || graph.rejected_incident_transitions
            != incident_pairs.size()) {
        throw IncompleteCanonicalMatRectangleGraphError(
            "canonical rectangle adaptor traversal is incomplete");
    }

    for (const MatExactGraphEdge2& edge :
         graph.edges) {
        if (edge.source_node_id
            == edge.target_node_id) {
            throw InvalidCanonicalMatRectangleNodeError(
                "canonical rectangle emitted a zero-dimensional edge");
        }
    }
    canonicalize_nodes(
        graph,
        aliases);
    const auto collapsed =
        std::find_if(
            graph.edges.begin(),
            graph.edges.end(),
            [](const MatExactGraphEdge2& edge) {
                return edge.source_node_id
                    == edge.target_node_id;
            });
    if (collapsed != graph.edges.end()) {
        throw InvalidCanonicalMatRectangleNodeError(
            "canonical rectangle graph collapsed an exact dual");
    }
    std::sort(
        graph.nodes.begin(),
        graph.nodes.end(),
        [](const MatExactGraphNode2& left,
           const MatExactGraphNode2& right) {
            return left.node_id < right.node_id;
        });
    std::sort(
        graph.edges.begin(),
        graph.edges.end(),
        [](const MatExactGraphEdge2& left,
           const MatExactGraphEdge2& right) {
            return left.edge_id < right.edge_id;
        });
    return graph;
}
