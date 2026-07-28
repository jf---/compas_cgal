#include "segment_site_endpoint_binding.h"
#include "segment_site_clipping.h"
#include "segment_site_rational_sources.h"
#include "segment_site_voronoi.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

compas::RowMatrixXd rectangle()
{
    const std::array<std::array<double, 2>, 4> points{{
        {-4.0, -2.0},
        {4.0, -2.0},
        {4.0, 2.0},
        {-4.0, 2.0},
    }};
    compas::RowMatrixXd result(
        static_cast<Eigen::Index>(points.size()),
        2);
    for (std::size_t index = 0;
         index < points.size();
         ++index) {
        result(
            static_cast<Eigen::Index>(index),
            0) = points[index][0];
        result(
            static_cast<Eigen::Index>(index),
            1) = points[index][1];
    }
    return result;
}

CanonicalMatSiteCatalog2 rectangle_catalog()
{
    return canonical_mat_site_catalog(
        canonical_reach_input(
            rectangle(),
            {},
            1.0));
}

std::vector<GeneratorSite2> linear_generators(
    const CanonicalMatSiteCatalog2& catalog)
{
    std::vector<GeneratorSite2> generators;
    generators.reserve(catalog.sites().size());
    for (const CanonicalMatSite2& site :
         catalog.sites()) {
        generators.push_back(
            {
                site.stable_id(),
                site.site(),
            });
    }
    return generators;
}

std::vector<MatExactOpenSegmentSource2>
segment_limiters(
    const CanonicalMatRationalSources2& sources)
{
    std::vector<MatExactOpenSegmentSource2> segments;
    segments.reserve(sources.segments().size());
    for (const auto& segment :
         sources.segments()) {
        segments.push_back(segment.support);
    }
    return segments;
}

bool endpoints_equal(
    const std::pair<
        MatParameterEndpoint2,
        MatParameterEndpoint2>& left,
    const std::pair<
        MatParameterEndpoint2,
        MatParameterEndpoint2>& right)
{
    ExactAlgebraicKernel1 kernel;
    const auto endpoint_equal =
        [&kernel](
            const MatParameterEndpoint2& first,
            const MatParameterEndpoint2& second) {
            if (!first.parameter.has_value()
                || !second.parameter.has_value()
                || first.provenance_ids
                    != second.provenance_ids) {
                return false;
            }
            return kernel.compare_1_object()(
                       *first.parameter,
                       *second.parameter)
                == CGAL::EQUAL;
        };
    return endpoint_equal(left.first, right.first)
        && endpoint_equal(left.second, right.second);
}

RationalDomainRoot2 rational_bound(
    const RationalPrimitiveParameterization2&
        primitive,
    const CORE::BigRat& parameter,
    const CanonicalMatRationalSources2& sources)
{
    std::vector<std::string> provenance;
    for (const auto& point : sources.points()) {
        if (parallel_segment_tangent_parameter(
                primitive,
                point.x,
                point.y)
            == parameter) {
            provenance.push_back(
                point.stable_site_id);
        }
    }
    std::sort(
        provenance.begin(),
        provenance.end());
    return {
        parameter,
        std::move(provenance),
    };
}

bool parallel_binding_uses_same_exact_records()
{
    const CanonicalMatSiteCatalog2 catalog =
        rectangle_catalog();
    const CanonicalMatRationalSources2 sources =
        CanonicalMatRationalSources2::build(
            catalog);
    CanonicalMatDelaunaySource2 delaunay =
        CanonicalMatDelaunaySource2::build(
            catalog);
    const CanonicalMatVoronoiSource2 source =
        CanonicalMatVoronoiSource2::build(
            std::move(delaunay));
    const std::vector<GeneratorSite2> generators =
        linear_generators(catalog);
    const std::vector<MatExactOpenSegmentSource2>
        limiters =
            segment_limiters(sources);
    const auto& bottom = sources.segments()[0];
    const auto& top = sources.segments()[2];
    const std::vector<std::string> generator_ids =
        ordered_generator_site_ids(
            bottom.stable_site_id,
            top.stable_site_id);
    RationalPrimitiveParameterization2 primitive =
        parallel_segment_bisector_parameterization(
            bottom.support,
            top.support);
    primitive.domain_lower = -4;
    primitive.domain_upper = 4;
    const RationalDomainRoot2 lower =
        rational_bound(
            primitive,
            -4,
            sources);
    const RationalDomainRoot2 upper =
        rational_bound(
            primitive,
            4,
            sources);

    std::size_t matches = 0;
    for (auto halfedge =
             source.voronoi().halfedges_begin();
         halfedge
         != source.voronoi().halfedges_end();
         ++halfedge) {
        const MatTraits::Site_2 up =
            halfedge->up()->site();
        const MatTraits::Site_2 down =
            halfedge->down()->site();
        if (!up.is_segment()
            || !down.is_segment()) {
            continue;
        }
        const std::string up_id =
            source.site_index().stable_id(up);
        if (up_id != generator_ids.front()
            || ordered_generator_site_ids(
                   up_id,
                   source.site_index().stable_id(
                       down))
                != generator_ids) {
            continue;
        }
        const auto linear =
            bind_parallel_segment_segment_cell_endpoints(
                primitive,
                lower,
                upper,
                bottom.support,
                top.support,
                {},
                limiters,
                generator_ids,
                generators,
                source.voronoi(),
                halfedge);
        const auto indexed =
            bind_parallel_segment_segment_cell_endpoints(
                primitive,
                lower,
                upper,
                bottom.support,
                top.support,
                {},
                limiters,
                generator_ids,
                source.site_index(),
                source.voronoi(),
                halfedge);
        if (!endpoints_equal(linear, indexed)) {
            return false;
        }
        ++matches;
    }
    return matches == 1;
}

MatQuadraticFieldDomainBoundary2
nonparallel_boundary(
    const NonparallelSegmentBisectorParameterization2&
        primitive,
    const MatExactOpenSegmentSource2& segment,
    const CanonicalMatRationalPointSource2& point)
{
    return {
        nonparallel_segment_tangent_parameter(
            primitive,
            segment,
            point.x,
            point.y),
        {point.stable_site_id},
    };
}

bool nonparallel_binding_uses_same_exact_records()
{
    const CanonicalMatSiteCatalog2 catalog =
        rectangle_catalog();
    const CanonicalMatRationalSources2 sources =
        CanonicalMatRationalSources2::build(
            catalog);
    CanonicalMatDelaunaySource2 delaunay =
        CanonicalMatDelaunaySource2::build(
            catalog);
    const CanonicalMatVoronoiSource2 source =
        CanonicalMatVoronoiSource2::build(
            std::move(delaunay));
    const std::vector<GeneratorSite2> generators =
        linear_generators(catalog);
    const std::vector<MatExactOpenSegmentSource2>
        limiters =
            segment_limiters(sources);
    const auto& bottom = sources.segments()[0];
    const auto& left = sources.segments()[3];
    const std::vector<std::string> generator_ids =
        ordered_generator_site_ids(
            bottom.stable_site_id,
            left.stable_site_id);

    std::size_t matches = 0;
    for (auto halfedge =
             source.voronoi().halfedges_begin();
         halfedge
         != source.voronoi().halfedges_end();
         ++halfedge) {
        const MatTraits::Site_2 up =
            halfedge->up()->site();
        const MatTraits::Site_2 down =
            halfedge->down()->site();
        if (!up.is_segment()
            || !down.is_segment()) {
            continue;
        }
        const std::string up_id =
            source.site_index().stable_id(up);
        if (up_id != generator_ids.front()
            || ordered_generator_site_ids(
                   up_id,
                   source.site_index().stable_id(
                       down))
                != generator_ids) {
            continue;
        }
        MatTraits::Segment_2 representative;
        if (!CGAL::assign(
                representative,
                source.voronoi().dual().primal(
                    halfedge->dual()))) {
            return false;
        }
        const NonparallelSegmentBisectorParameterization2
            primitive =
                nonparallel_segment_bisector_parameterization(
                    bottom.support,
                    left.support,
                    representative);
        const NonparallelSegmentFeatureDomain2
            feature_domain =
                intersect_nonparallel_segment_feature_domains(
                    nonparallel_boundary(
                        primitive,
                        bottom.support,
                        sources.points()[0]),
                    nonparallel_boundary(
                        primitive,
                        bottom.support,
                        sources.points()[1]),
                    nonparallel_boundary(
                        primitive,
                        left.support,
                        sources.points()[3]),
                    nonparallel_boundary(
                        primitive,
                        left.support,
                        sources.points()[0]),
                    primitive.radicand);
        const auto linear =
            bind_nonparallel_segment_segment_cell_endpoints(
                primitive,
                feature_domain,
                bottom.support,
                left.support,
                limiters,
                generator_ids,
                generators,
                source.voronoi(),
                halfedge);
        const auto indexed =
            bind_nonparallel_segment_segment_cell_endpoints(
                primitive,
                feature_domain,
                bottom.support,
                left.support,
                limiters,
                generator_ids,
                source.site_index(),
                source.voronoi(),
                halfedge);
        if (!endpoints_equal(linear, indexed)) {
            return false;
        }
        ++matches;
    }
    return matches == 1;
}

compas::RowMatrixXd l_shape(
    const bool transformed)
{
    const std::array<std::array<double, 2>, 6>
        canonical{{
            {0.0, 0.0},
            {6.0, 0.0},
            {6.0, 2.0},
            {2.0, 2.0},
            {2.0, 6.0},
            {0.0, 6.0},
        }};
    const std::array<std::array<double, 2>, 6>
        reversed{{
            {0.0, 6.0},
            {2.0, 6.0},
            {2.0, 2.0},
            {6.0, 2.0},
            {6.0, 0.0},
            {0.0, 0.0},
        }};
    const auto& points =
        transformed ? reversed : canonical;
    compas::RowMatrixXd result(
        static_cast<Eigen::Index>(points.size()),
        2);
    for (std::size_t index = 0;
         index < points.size();
         ++index) {
        result(
            static_cast<Eigen::Index>(index),
            0) = points[index][0];
        result(
            static_cast<Eigen::Index>(index),
            1) = points[index][1];
    }
    return result;
}

MatDomainPolygonWithHoles2 l_shape_domain()
{
    MatDomainPolygon2 outer;
    outer.push_back({0, 0});
    outer.push_back({6, 0});
    outer.push_back({6, 2});
    outer.push_back({2, 2});
    outer.push_back({2, 6});
    outer.push_back({0, 6});
    return MatDomainPolygonWithHoles2(outer);
}

MatExactPointSiteSource2 exact_point_source(
    const CanonicalMatRationalPointSource2&
        point)
{
    return {
        point.stable_site_id,
        {point.x, 0},
        {point.y, 0},
        1,
    };
}

const CanonicalMatRationalPointSource2&
rational_point(
    const CanonicalMatRationalSources2& sources,
    const std::string& stable_id)
{
    const auto found = std::find_if(
        sources.points().begin(),
        sources.points().end(),
        [&stable_id](const auto& point) {
            return point.stable_site_id
                == stable_id;
        });
    if (found == sources.points().end()) {
        throw UnknownCanonicalMatParabolaSourceError(
            "P-S gate point is outside the rational sources");
    }
    return *found;
}

const CanonicalMatRationalOpenSegmentSource2&
rational_segment(
    const CanonicalMatRationalSources2& sources,
    const std::string& stable_id)
{
    const auto found = std::find_if(
        sources.segments().begin(),
        sources.segments().end(),
        [&stable_id](const auto& segment) {
            return segment.stable_site_id
                == stable_id;
        });
    if (found == sources.segments().end()) {
        throw UnknownCanonicalMatParabolaSourceError(
            "P-S gate segment is outside the rational sources");
    }
    return *found;
}

std::pair<
    MatParameterEndpoint2,
    MatParameterEndpoint2>
ordered_endpoints(
    MatParameterEndpoint2 first,
    MatParameterEndpoint2 second)
{
    if (!first.parameter.has_value()
        || !second.parameter.has_value()) {
        return {
            std::move(first),
            std::move(second),
        };
    }
    ExactAlgebraicKernel1 kernel;
    if (kernel.compare_1_object()(
            *second.parameter,
            *first.parameter)
        == CGAL::SMALLER) {
        std::swap(first, second);
    }
    return {
        std::move(first),
        std::move(second),
    };
}

using EndpointSignature = std::tuple<
    std::string,
    std::vector<std::string>>;

using ParabolaBindingSignature = std::tuple<
    std::vector<std::string>,
    EndpointSignature,
    EndpointSignature>;

EndpointSignature endpoint_signature(
    const MatParameterEndpoint2& endpoint)
{
    if (!endpoint.parameter.has_value()) {
        return {};
    }
    return {
        algebraic_root_identity_v1(
            *endpoint.parameter),
        endpoint.provenance_ids,
    };
}

std::vector<ParabolaBindingSignature>
l_shape_parabola_bindings(
    const bool transformed)
{
    const CanonicalMatSiteCatalog2 catalog =
        canonical_mat_site_catalog(
            canonical_reach_input(
                l_shape(transformed),
                {},
                1.0));
    const CanonicalMatRationalSources2 sources =
        CanonicalMatRationalSources2::build(
            catalog);
    CanonicalMatDelaunaySource2 delaunay =
        CanonicalMatDelaunaySource2::build(
            catalog);
    const CanonicalMatVoronoiSource2 source =
        CanonicalMatVoronoiSource2::build(
            std::move(delaunay));
    const auto& focus_record =
        sources.points()[3];
    const MatExactPointSiteSource2 focus =
        exact_point_source(focus_record);
    std::vector<ParabolaBindingSignature>
        signatures;

    for (auto halfedge =
             source.voronoi().halfedges_begin();
         halfedge
         != source.voronoi().halfedges_end();
         ++halfedge) {
        const MatTraits::Site_2 up =
            halfedge->up()->site();
        const MatTraits::Site_2 down =
            halfedge->down()->site();
        if (up.is_point() == down.is_point()) {
            continue;
        }
        const std::string up_id =
            source.site_index().stable_id(up);
        const std::string down_id =
            source.site_index().stable_id(down);
        const std::vector<std::string>
            generator_ids =
                ordered_generator_site_ids(
                    up_id,
                    down_id);
        if (up_id != generator_ids.front()) {
            continue;
        }
        const std::string point_id =
            up.is_point() ? up_id : down_id;
        const std::string segment_id =
            up.is_segment() ? up_id : down_id;
        if (point_id
                != focus_record.stable_site_id
            || (segment_id
                    != sources.segments()[0]
                           .stable_site_id
                && segment_id
                    != sources.segments()[5]
                           .stable_site_id)) {
            continue;
        }
        if (!halfedge->has_source()
            || !halfedge->has_target()) {
            return {};
        }
        SegmentSiteParabola2 live;
        if (!CGAL::assign(
                live,
                source.voronoi()
                    .dual()
                    .primal(
                        halfedge->dual()))) {
            return {};
        }
        require_distinct_live_parabola_endpoints(
            live);

        const auto indexed =
            bind_point_segment_cell_endpoints(
                sources,
                point_id,
                segment_id,
                source.site_index(),
                source.voronoi(),
                halfedge);
        const auto& source_segment =
            rational_segment(
                sources,
                segment_id);
        const MatExactPointSiteSource2
            segment_source =
                exact_point_source(
                    rational_point(
                        sources,
                        source_segment
                            .source_point_id));
        const MatExactPointSiteSource2
            segment_target =
                exact_point_source(
                    rational_point(
                        sources,
                        source_segment
                            .target_point_id));
        const auto bind_owner =
            [&sources,
             &focus,
             &source_segment,
             &segment_source,
             &segment_target,
             &source,
             &halfedge](
                const MatTraits::Site_2& owner) {
                const auto& limiter =
                    rational_segment(
                        sources,
                        source.site_index()
                            .stable_id(owner));
                return bind_segment_limiter_parabola_endpoint(
                    focus,
                    source_segment.support,
                    segment_source,
                    segment_target,
                    limiter.support,
                    exact_point_source(
                        rational_point(
                            sources,
                            limiter
                                .source_point_id)),
                    exact_point_source(
                        rational_point(
                            sources,
                            limiter
                                .target_point_id)),
                    source.voronoi(),
                    halfedge);
            };
        const auto direct =
            ordered_endpoints(
                bind_owner(
                    halfedge->left()
                        ->site()),
                bind_owner(
                    halfedge->right()
                        ->site()));
        if (!endpoints_equal(
                indexed,
                direct)) {
            return {};
        }
        signatures.emplace_back(
            generator_ids,
            endpoint_signature(
                indexed.first),
            endpoint_signature(
                indexed.second));
    }
    std::sort(
        signatures.begin(),
        signatures.end());
    return signatures;
}

bool point_segment_binding_uses_catalog_index()
{
    const auto canonical =
        l_shape_parabola_bindings(false);
    return canonical.size() == 2
        && canonical
            == l_shape_parabola_bindings(false)
        && canonical
            == l_shape_parabola_bindings(true);
}

bool point_segment_binding_rejects_invalid_sources()
{
    const CanonicalMatSiteCatalog2 catalog =
        canonical_mat_site_catalog(
            canonical_reach_input(
                l_shape(false),
                {},
                1.0));
    const CanonicalMatRationalSources2 sources =
        CanonicalMatRationalSources2::build(
            catalog);
    CanonicalMatDelaunaySource2 delaunay =
        CanonicalMatDelaunaySource2::build(
            catalog);
    const CanonicalMatVoronoiSource2 source =
        CanonicalMatVoronoiSource2::build(
            std::move(delaunay));
    const std::string focus_id =
        sources.points()[3].stable_site_id;
    const std::string segment_id =
        sources.segments()[0].stable_site_id;

    for (auto halfedge =
             source.voronoi().halfedges_begin();
         halfedge
         != source.voronoi().halfedges_end();
         ++halfedge) {
        const std::string up_id =
            source.site_index().stable_id(
                halfedge->up()->site());
        const std::string down_id =
            source.site_index().stable_id(
                halfedge->down()->site());
        if (up_id
                != ordered_generator_site_ids(
                       up_id,
                       down_id)
                       .front()
            || ordered_generator_site_ids(
                   up_id,
                   down_id)
                != ordered_generator_site_ids(
                    focus_id,
                    segment_id)) {
            continue;
        }
        bool unknown_point = false;
        try {
            static_cast<void>(
                bind_point_segment_cell_endpoints(
                    sources,
                    "unknown-point",
                    segment_id,
                    source.site_index(),
                    source.voronoi(),
                    halfedge));
        } catch (
            const UnknownCanonicalMatParabolaSourceError&) {
            unknown_point = true;
        }
        bool unknown_segment = false;
        try {
            static_cast<void>(
                bind_point_segment_cell_endpoints(
                    sources,
                    focus_id,
                    "unknown-segment",
                    source.site_index(),
                    source.voronoi(),
                    halfedge));
        } catch (
            const UnknownCanonicalMatParabolaSourceError&) {
            unknown_segment = true;
        }
        bool incident = false;
        try {
            static_cast<void>(
                bind_point_segment_cell_endpoints(
                    sources,
                    sources.points()[0]
                        .stable_site_id,
                    segment_id,
                    source.site_index(),
                    source.voronoi(),
                    halfedge));
        } catch (
            const IncidentCanonicalMatParabolaSourceError&) {
            incident = true;
        }
        bool mismatched = false;
        try {
            static_cast<void>(
                bind_point_segment_cell_endpoints(
                    sources,
                    sources.points()[4]
                        .stable_site_id,
                    segment_id,
                    source.site_index(),
                    source.voronoi(),
                    halfedge));
        } catch (
            const MismatchedLiveParabolaBridgeError&) {
            mismatched = true;
        }
        return unknown_point
            && unknown_segment
            && incident
            && mismatched;
    }
    return false;
}

bool primitives_equal(
    const RationalPrimitiveParameterization2& left,
    const RationalPrimitiveParameterization2& right)
{
    return left.x_coefficients
            == right.x_coefficients
        && left.y_coefficients
            == right.y_coefficients
        && left.domain_lower
            == right.domain_lower
        && left.domain_upper
            == right.domain_upper;
}

std::vector<RationalPrimitiveParameterization2>
l_shape_point_point_ray_bindings(
    const bool transformed)
{
    const CanonicalMatSiteCatalog2 catalog =
        canonical_mat_site_catalog(
            canonical_reach_input(
                l_shape(transformed),
                {},
                1.0));
    const CanonicalMatRationalSources2 sources =
        CanonicalMatRationalSources2::build(
            catalog);
    CanonicalMatDelaunaySource2 delaunay =
        CanonicalMatDelaunaySource2::build(
            catalog);
    const CanonicalMatVoronoiSource2 source =
        CanonicalMatVoronoiSource2::build(
            std::move(delaunay));
    const auto& first = sources.points()[2];
    const auto& second = sources.points()[4];
    const std::vector<std::string> generator_ids =
        ordered_generator_site_ids(
            first.stable_site_id,
            second.stable_site_id);
    std::vector<RationalPrimitiveParameterization2>
        bindings;

    for (auto halfedge =
             source.voronoi().halfedges_begin();
         halfedge
         != source.voronoi().halfedges_end();
         ++halfedge) {
        const MatTraits::Site_2 up =
            halfedge->up()->site();
        const MatTraits::Site_2 down =
            halfedge->down()->site();
        if (!up.is_point()
            || !down.is_point()) {
            continue;
        }
        const std::string up_id =
            source.site_index().stable_id(up);
        const std::string down_id =
            source.site_index().stable_id(
                down);
        if (up_id != generator_ids.front()
            || ordered_generator_site_ids(
                   up_id,
                   down_id)
                != generator_ids) {
            continue;
        }
        RationalPrimitiveParameterization2
            primitive =
                bind_point_point_ray_parameterization(
                    sources,
                    generator_ids,
                    source.site_index(),
                    source.voronoi(),
                    halfedge);
        const std::string dual_id =
            stable_dual_identity_v1(
                "point",
                generator_ids);
        if (!clip_linear_clearance_components(
                 dual_id,
                 primitive,
                 point_clearance_boundary(
                     primitive,
                     first.x,
                     first.y,
                     0),
                 l_shape_domain())
                 .empty()) {
            return {};
        }
        bindings.push_back(
            std::move(primitive));
    }
    return bindings;
}

bool point_point_ray_binding_is_exact_and_invariant()
{
    const auto canonical =
        l_shape_point_point_ray_bindings(
            false);
    const auto repeated =
        l_shape_point_point_ray_bindings(
            false);
    const auto reversed =
        l_shape_point_point_ray_bindings(
            true);
    if (canonical.size() != 1
        || repeated.size() != 1
        || reversed.size() != 1) {
        return false;
    }
    const RationalPrimitiveParameterization2&
        primitive = canonical.front();
    return primitive.x_coefficients
            == std::vector<CORE::BigRat>{
                4,
                -4,
            }
        && primitive.y_coefficients
            == std::vector<CORE::BigRat>{
                4,
                -4,
            }
        && !primitive.domain_lower.has_value()
        && primitive.domain_upper
            == std::optional<CORE::BigRat>(
                CORE::BigRat(-1, 2))
        && primitives_equal(
            primitive,
            repeated.front())
        && primitives_equal(
            primitive,
            reversed.front());
}

bool point_point_ray_binding_rejects_invalid_sources()
{
    const CanonicalMatSiteCatalog2 catalog =
        canonical_mat_site_catalog(
            canonical_reach_input(
                l_shape(false),
                {},
                1.0));
    const CanonicalMatRationalSources2 sources =
        CanonicalMatRationalSources2::build(
            catalog);
    CanonicalMatDelaunaySource2 delaunay =
        CanonicalMatDelaunaySource2::build(
            catalog);
    const CanonicalMatVoronoiSource2 source =
        CanonicalMatVoronoiSource2::build(
            std::move(delaunay));
    const std::vector<std::string> generator_ids =
        ordered_generator_site_ids(
            sources.points()[2]
                .stable_site_id,
            sources.points()[4]
                .stable_site_id);

    for (auto halfedge =
             source.voronoi().halfedges_begin();
         halfedge
         != source.voronoi().halfedges_end();
         ++halfedge) {
        const std::string up_id =
            source.site_index().stable_id(
                halfedge->up()->site());
        const std::string down_id =
            source.site_index().stable_id(
                halfedge->down()->site());
        if (up_id != generator_ids.front()
            || ordered_generator_site_ids(
                   up_id,
                   down_id)
                != generator_ids) {
            continue;
        }
        bool unknown = false;
        try {
            static_cast<void>(
                bind_point_point_ray_parameterization(
                    sources,
                    ordered_generator_site_ids(
                        "unknown-point",
                        generator_ids.back()),
                    source.site_index(),
                    source.voronoi(),
                    halfedge));
        } catch (
            const UnknownCanonicalMatPointRaySourceError&) {
            unknown = true;
        }
        bool mismatched = false;
        try {
            static_cast<void>(
                bind_point_point_ray_parameterization(
                    sources,
                    ordered_generator_site_ids(
                        sources.points()[0]
                            .stable_site_id,
                        generator_ids.back()),
                    source.site_index(),
                    source.voronoi(),
                    halfedge));
        } catch (
            const MismatchedLivePointPointRayError&) {
            mismatched = true;
        }
        bool duplicate = false;
        try {
            static_cast<void>(
                bind_point_point_ray_parameterization(
                    sources,
                    {
                        generator_ids.front(),
                        generator_ids.front(),
                    },
                    source.site_index(),
                    source.voronoi(),
                    halfedge));
        } catch (
            const InvalidDualIdentityError&) {
            duplicate = true;
        }
        return unknown
            && mismatched
            && duplicate;
    }
    return false;
}

} // namespace

bool indexed_endpoint_binding_gate()
{
    return parallel_binding_uses_same_exact_records()
        && nonparallel_binding_uses_same_exact_records()
        && point_segment_binding_uses_catalog_index()
        && point_segment_binding_rejects_invalid_sources()
        && point_point_ray_binding_is_exact_and_invariant()
        && point_point_ray_binding_rejects_invalid_sources();
}
