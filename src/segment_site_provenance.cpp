#include "segment_site_provenance.h"

#include <algorithm>
#include <iterator>
#include <limits>
#include <set>
#include <utility>

#include <CGAL/Polynomial_traits_d.h>

namespace {

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

std::string length_framed(const std::string& value)
{
    return std::to_string(value.size()) + ":" + value;
}

} // namespace

std::string algebraic_root_identity_v1(
    const ExactAlgebraicKernel1::Algebraic_real_1& root)
{
    using Polynomial = ExactAlgebraicKernel1::Polynomial_1;
    ExactAlgebraicKernel1 kernel;
    const Polynomial factor =
        typename CGAL::Polynomial_traits_d<
            Polynomial>::Canonicalize()(
            kernel.compute_polynomial_1_object()(root));
    const int degree = CGAL::degree(factor);
    std::vector<ExactAlgebraicInteger1> coefficients;
    coefficients.reserve(
        static_cast<std::size_t>(degree + 1));
    for (int index = 0; index <= degree; ++index) {
        coefficients.push_back(factor[index]);
    }

    std::vector<ExactAlgebraicKernel1::Algebraic_real_1> roots;
    kernel.solve_1_object()(
        factor,
        true,
        std::back_inserter(roots));
    for (std::size_t ordinal = 0;
         ordinal < roots.size();
         ++ordinal) {
        if (kernel.compare_1_object()(roots[ordinal], root)
            == CGAL::EQUAL) {
            return algebraic_root_id_v1(
                coefficients,
                ordinal);
        }
    }
    throw InvalidAlgebraicPolynomialError(
        "canonical factor does not contain clipping root");
}

std::size_t matched_generator_site_count(
    const SegmentSiteDelaunay2& delaunay,
    const std::vector<GeneratorSite2>& generators)
{
    std::set<std::string> matched_ids;
    for (auto vertex = delaunay.finite_vertices_begin();
         vertex != delaunay.finite_vertices_end();
         ++vertex) {
        const MatTraits::Site_2& site = vertex->site();
        const GeneratorSite2* match = nullptr;
        for (const GeneratorSite2& generator : generators) {
            if (!exact_site_equal(site, generator.site)) {
                continue;
            }
            if (match != nullptr) {
                return 0;
            }
            match = &generator;
        }
        if (match == nullptr
            || !matched_ids.insert(match->stable_id).second) {
            return 0;
        }
    }
    return matched_ids.size();
}

std::size_t require_generator_site_bijection(
    const SegmentSiteDelaunay2& delaunay,
    const std::vector<GeneratorSite2>& generators)
{
    std::set<std::string> matched_ids;
    for (auto vertex = delaunay.finite_vertices_begin();
         vertex != delaunay.finite_vertices_end();
         ++vertex) {
        const MatTraits::Site_2& site = vertex->site();
        const GeneratorSite2* match = nullptr;
        for (const GeneratorSite2& generator : generators) {
            if (!exact_site_equal(site, generator.site)) {
                continue;
            }
            if (match != nullptr) {
                throw GeneratorSiteBijectionError(
                    "live site maps to multiple caller identities");
            }
            match = &generator;
        }
        if (match == nullptr) {
            throw GeneratorSiteBijectionError(
                "live site has no caller identity");
        }
        if (!matched_ids.insert(match->stable_id).second) {
            throw GeneratorSiteBijectionError(
                "caller identity maps to multiple live sites");
        }
    }
    if (matched_ids.size() != generators.size()) {
        throw GeneratorSiteBijectionError(
            "caller identity has no live site");
    }
    return matched_ids.size();
}

std::string stable_generator_site_id(
    const MatTraits::Site_2& site,
    const std::vector<GeneratorSite2>& generators)
{
    const GeneratorSite2* match = nullptr;
    for (const GeneratorSite2& generator : generators) {
        if (!exact_site_equal(site, generator.site)) {
            continue;
        }
        if (match != nullptr) {
            throw InvalidRationalPrimitiveError(
                "generator-site provenance is not bijective");
        }
        match = &generator;
    }
    if (match == nullptr) {
        throw InvalidRationalPrimitiveError(
            "generator site has no stable source record");
    }
    return match->stable_id;
}

std::vector<std::string> ordered_generator_site_ids(
    std::string first,
    std::string second)
{
    if (second < first) {
        std::swap(first, second);
    }
    return {std::move(first), std::move(second)};
}

void union_stable_ids(
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

MatEndpointBoundaryFeature2 design_line_endpoint_feature(
    const std::size_t ring,
    const std::size_t feature)
{
    const auto max_index =
        static_cast<std::uintmax_t>(
            std::numeric_limits<std::int64_t>::max());
    if (static_cast<std::uintmax_t>(ring)
            > max_index
        || static_cast<std::uintmax_t>(feature)
            > max_index) {
        throw MatEndpointFeatureOverflowError(
            "MAT endpoint boundary feature exceeds int64 range");
    }
    return {
        MatEndpointDomainKind2::Design,
        0,
        MatExactCurveKind2::Line,
        static_cast<std::int64_t>(ring),
        static_cast<std::int64_t>(feature),
    };
}

void union_endpoint_boundary_features(
    std::vector<MatEndpointBoundaryFeature2>& target,
    const std::vector<MatEndpointBoundaryFeature2>& source)
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

void union_endpoint_evidence(
    MatParameterEndpoint2& target,
    const MatParameterEndpoint2& source)
{
    union_stable_ids(
        target.provenance_ids,
        source.provenance_ids);
    target.exact_evidence.original_voronoi_vertex =
        target.exact_evidence.original_voronoi_vertex
        || source.exact_evidence
               .original_voronoi_vertex;
    target.exact_evidence.domain_boundary =
        target.exact_evidence.domain_boundary
        || source.exact_evidence.domain_boundary;
    target.exact_evidence.clearance_root =
        target.exact_evidence.clearance_root
        || source.exact_evidence.clearance_root;
    auto& features =
        target.exact_evidence.boundary_features;
    union_endpoint_boundary_features(
        features,
        source.exact_evidence.boundary_features);
}

std::string stable_dual_identity_v1(
    const std::string& dual_kind,
    const std::vector<std::string>& ordered_generator_ids)
{
    if (dual_kind.empty()) {
        throw InvalidDualIdentityError(
            "dual kind is empty");
    }
    if (ordered_generator_ids.size() != 2) {
        throw InvalidDualIdentityError(
            "dual identity requires two generators");
    }
    if (ordered_generator_ids[0].empty()
        || ordered_generator_ids[1].empty()) {
        throw InvalidDualIdentityError(
            "dual generator identity is empty");
    }
    if (ordered_generator_ids[1]
        <= ordered_generator_ids[0]) {
        throw InvalidDualIdentityError(
            "dual generator identities are not strictly ordered");
    }
    return "mat-dual/v1/"
        + length_framed(dual_kind)
        + "/" + length_framed(ordered_generator_ids[0])
        + "/" + length_framed(ordered_generator_ids[1]);
}

std::string stable_voronoi_node_identity_v1(
    const std::vector<std::string>& ordered_generator_ids)
{
    if (ordered_generator_ids.size() != 3) {
        throw InvalidDualIdentityError(
            "Voronoi node identity requires three generators");
    }
    std::string identity = "voronoi-node/v1";
    for (std::size_t index = 0;
         index < ordered_generator_ids.size();
         ++index) {
        const std::string& generator_id =
            ordered_generator_ids[index];
        if (generator_id.empty()) {
            throw InvalidDualIdentityError(
                "Voronoi node generator identity is empty");
        }
        if (index > 0
            && generator_id
                <= ordered_generator_ids[index - 1]) {
            throw InvalidDualIdentityError(
                "Voronoi node generators are not strictly ordered");
        }
        identity += "/" + length_framed(generator_id);
    }
    return identity;
}

std::string stable_normalized_voronoi_node_identity_v1(
    const std::vector<std::string>& ordered_generator_ids)
{
    if (ordered_generator_ids.size() < 3) {
        throw InvalidDualIdentityError(
            "normalized Voronoi node identity requires at least three generators");
    }
    if (ordered_generator_ids.size() == 3) {
        return stable_voronoi_node_identity_v1(
            ordered_generator_ids);
    }

    std::string identity =
        "normalized-voronoi-node/v1";
    for (std::size_t index = 0;
         index < ordered_generator_ids.size();
         ++index) {
        const std::string& generator_id =
            ordered_generator_ids[index];
        if (generator_id.empty()) {
            throw InvalidDualIdentityError(
                "normalized Voronoi node generator identity is empty");
        }
        if (index > 0
            && generator_id
                <= ordered_generator_ids[index - 1]) {
            throw InvalidDualIdentityError(
                "normalized Voronoi node generators are not strictly ordered");
        }
        identity += "/"
            + length_framed(generator_id);
    }
    return identity;
}

MatParameterEndpoint2 exact_graph_endpoint_binding(
    const MatParameterEndpoint2& endpoint)
{
    if (!endpoint.parameter.has_value()) {
        throw InvalidRationalPrimitiveError(
            "exact graph endpoint is unbounded");
    }
    MatParameterEndpoint2 bound = endpoint;
    union_stable_ids(
        bound.provenance_ids,
        {algebraic_root_identity_v1(*bound.parameter)});
    return bound;
}

std::string stable_endpoint_node_identity_v1(
    const std::string& dual_id,
    const MatParameterEndpoint2& endpoint)
{
    if (!endpoint.parameter.has_value()) {
        throw InvalidRationalPrimitiveError(
            "exact graph endpoint is unbounded");
    }
    return dual_id + "/node/"
        + algebraic_root_identity_v1(*endpoint.parameter);
}

MatParameterEndpoint2 domain_endpoint(
    const std::string& dual_id,
    const char* side,
    const std::optional<CORE::BigRat>& value,
    const ExactAlgebraicKernel1& kernel)
{
    const std::string provenance =
        dual_id + "/domain-" + side;
    if (!value.has_value()) {
        return {
            std::nullopt,
            {provenance + "-unbounded"},
        };
    }
    return {
        kernel.construct_algebraic_real_1_object()(*value),
        {provenance},
    };
}


using AlgebraicParameter =
    ExactAlgebraicKernel1::Algebraic_real_1;

void append_root_endpoint(
    std::vector<TaggedEndpoint2>& endpoints,
    TaggedEndpoint2& upper,
    const AlgebraicParameter& root,
    const std::string& root_id,
    const ExactAlgebraicKernel1& kernel)
{
    const auto compare = kernel.compare_1_object();
    if (endpoints.front().endpoint.parameter.has_value()
        && compare(
               root,
               *endpoints.front().endpoint.parameter)
            == CGAL::EQUAL) {
        endpoints.front().clearance_root = true;
        endpoints.front().endpoint.provenance_ids.push_back(
            root_id);
        endpoints.front().endpoint.exact_evidence
            .clearance_root = true;
        return;
    }
    if (upper.endpoint.parameter.has_value()
        && compare(root, *upper.endpoint.parameter)
            == CGAL::EQUAL) {
        upper.clearance_root = true;
        upper.endpoint.provenance_ids.push_back(root_id);
        upper.endpoint.exact_evidence.clearance_root =
            true;
        return;
    }
    endpoints.push_back(
        {
            {
                root,
                {root_id},
                {
                    false,
                    false,
                    true,
                    {},
                },
            },
            true,
        });
}
