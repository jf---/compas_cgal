#include "segment_site_provenance.h"

#include <iterator>
#include <set>

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
        return;
    }
    if (upper.endpoint.parameter.has_value()
        && compare(root, *upper.endpoint.parameter)
            == CGAL::EQUAL) {
        upper.clearance_root = true;
        upper.endpoint.provenance_ids.push_back(root_id);
        return;
    }
    endpoints.push_back(
        {
            {root, {root_id}},
            true,
        });
}
