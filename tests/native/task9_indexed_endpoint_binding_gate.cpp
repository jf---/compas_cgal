#include "segment_site_endpoint_binding.h"
#include "segment_site_rational_sources.h"
#include "segment_site_voronoi.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>
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

} // namespace

bool indexed_endpoint_binding_gate()
{
    return parallel_binding_uses_same_exact_records()
        && nonparallel_binding_uses_same_exact_records();
}
