#include "segment_site_mat.h"

#include <algorithm>
#include <iterator>
#include <set>
#include <string>
#include <vector>

#include <CGAL/Object.h>
#include <CGAL/Polynomial_traits_d.h>

namespace {

using RationalPolynomial = std::vector<CORE::BigRat>;

void trim(RationalPolynomial& polynomial)
{
    while (polynomial.size() > 1
           && polynomial.back() == 0) {
        polynomial.pop_back();
    }
}

RationalPolynomial multiply(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs)
{
    RationalPolynomial product(
        lhs.size() + rhs.size() - 1,
        CORE::BigRat(0));
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            product[i + j] += lhs[i] * rhs[j];
        }
    }
    trim(product);
    return product;
}

void add_in_place(
    RationalPolynomial& target,
    const RationalPolynomial& source)
{
    target.resize(
        std::max(target.size(), source.size()),
        CORE::BigRat(0));
    for (std::size_t i = 0; i < source.size(); ++i) {
        target[i] += source[i];
    }
    trim(target);
}

std::vector<ExactAlgebraicInteger1> primitive_integer_coefficients(
    const RationalPolynomial& polynomial)
{
    ExactAlgebraicInteger1 common_denominator = 1;
    for (const CORE::BigRat& coefficient : polynomial) {
        common_denominator *= CORE::denominator(coefficient);
    }

    std::vector<ExactAlgebraicInteger1> primitive;
    primitive.reserve(polynomial.size());
    for (const CORE::BigRat& coefficient : polynomial) {
        primitive.push_back(
            CORE::numerator(coefficient)
            * CORE::div_exact(
                common_denominator,
                CORE::denominator(coefficient)));
    }

    ExactAlgebraicInteger1 divisor = 0;
    for (const ExactAlgebraicInteger1& coefficient : primitive) {
        const ExactAlgebraicInteger1 magnitude =
            coefficient < 0 ? -coefficient : coefficient;
        divisor = CORE::gcd(divisor, magnitude);
    }
    if (divisor != 0 && divisor != 1) {
        for (ExactAlgebraicInteger1& coefficient : primitive) {
            coefficient = CORE::div_exact(
                coefficient,
                divisor);
        }
    }
    if (primitive.back() < 0) {
        for (ExactAlgebraicInteger1& coefficient : primitive) {
            coefficient = -coefficient;
        }
    }
    return primitive;
}

bool root_is_in_domain(
    const ExactAlgebraicKernel1::Algebraic_real_1& root,
    const RationalPrimitiveParameterization2& primitive,
    const ExactAlgebraicKernel1& kernel)
{
    const auto compare = kernel.compare_1_object();
    if (primitive.domain_lower.has_value()
        && compare(root, *primitive.domain_lower)
            == CGAL::SMALLER) {
        return false;
    }
    return !primitive.domain_upper.has_value()
        || compare(root, *primitive.domain_upper)
            != CGAL::LARGER;
}

struct GeneratorSite2 {
    std::string stable_id;
    MatTraits::Site_2 site;
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

} // namespace

MatParameterDomain2
exact_parameter_domain(const MatTraits::Line_2&)
{
    return {std::nullopt, std::nullopt};
}

MatParameterDomain2
exact_parameter_domain(const MatTraits::Ray_2&)
{
    return {MatTraits::FT(0), std::nullopt};
}

MatParameterDomain2
exact_parameter_domain(const MatTraits::Segment_2&)
{
    return {MatTraits::FT(0), MatTraits::FT(1)};
}

MatParameterDomain2
exact_parameter_domain(const SegmentSiteParabola2& parabola)
{
    MatTraits::FT lower = parabola.t(parabola.p1);
    MatTraits::FT upper = parabola.t(parabola.p2);
    if (CGAL::compare(lower, upper) == CGAL::LARGER) {
        std::swap(lower, upper);
    }
    return {lower, upper};
}

ClearanceRootBoundary2 point_clearance_boundary(
    const RationalPrimitiveParameterization2& primitive,
    const CORE::BigRat& site_x,
    const CORE::BigRat& site_y,
    const CORE::BigRat& radius_squared)
{
    if (primitive.x_coefficients.empty()
        || primitive.y_coefficients.empty()) {
        throw InvalidRationalPrimitiveError(
            "primitive coordinate polynomial is empty");
    }
    if (primitive.domain_lower.has_value()
        && primitive.domain_upper.has_value()
        && *primitive.domain_lower > *primitive.domain_upper) {
        throw InvalidRationalPrimitiveError(
            "primitive parameter domain is reversed");
    }

    RationalPolynomial dx = primitive.x_coefficients;
    RationalPolynomial dy = primitive.y_coefficients;
    dx.front() -= site_x;
    dy.front() -= site_y;
    trim(dx);
    trim(dy);

    RationalPolynomial clearance = multiply(dx, dx);
    add_in_place(clearance, multiply(dy, dy));
    clearance.front() -= radius_squared;
    trim(clearance);
    if (clearance.size() == 1) {
        return {
            CGAL::sign(clearance.front()),
            {},
            {},
        };
    }

    const std::vector<ExactAlgebraicInteger1> coefficients =
        primitive_integer_coefficients(clearance);
    using Polynomial = ExactAlgebraicKernel1::Polynomial_1;
    const Polynomial polynomial =
        typename CGAL::Polynomial_traits_d<
            Polynomial>::Construct_polynomial()(
            coefficients.begin(),
            coefficients.end());
    ExactAlgebraicKernel1 kernel;
    std::vector<ExactAlgebraicKernel1::Algebraic_real_1>
        isolated_roots;
    kernel.solve_1_object()(
        polynomial,
        true,
        std::back_inserter(isolated_roots));

    std::vector<ExactAlgebraicKernel1::Algebraic_real_1> roots;
    std::copy_if(
        isolated_roots.begin(),
        isolated_roots.end(),
        std::back_inserter(roots),
        [&primitive, &kernel](const auto& root) {
            return root_is_in_domain(
                root,
                primitive,
                kernel);
        });
    return {
        std::nullopt,
        coefficients,
        roots,
    };
}

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
