#include "segment_site_endpoint_binding.h"
#include "segment_site_delaunay.h"
#include "segment_site_rational_sources.h"

#include <algorithm>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

#include <CGAL/CORE/poly/Poly.h>
#include <CGAL/Fraction_traits.h>
#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/Object.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/number_utils.h>

namespace
{

using Algebraic = ExactAlgebraicKernel1::Algebraic_real_1;
using Polynomial = ExactAlgebraicKernel1::Polynomial_1;

struct ExactFactorRootWitness2
{
    std::vector<ExactAlgebraicInteger1> primitive_factor;
    std::size_t factor_local_ordinal;
    Algebraic root;
    std::string root_id;
    std::size_t source_factor_multiplicity;
};

struct LiveParabolaEndpointBridge2
{
    MatTraits::Point_2 point;
    SegmentSiteParabola2 parabola;
};

struct FieldPolynomial2
{
    RationalPolynomial rational;
    RationalPolynomial radical;
};

void normalize(RationalPolynomial& polynomial)
{
    trim(polynomial);
    if (polynomial.empty())
    {
        polynomial.push_back(0);
    }
}

RationalPolynomial scaled(RationalPolynomial polynomial,
                          const CORE::BigRat& factor)
{
    for (CORE::BigRat& coefficient : polynomial)
    {
        coefficient *= factor;
    }
    normalize(polynomial);
    return polynomial;
}

RationalPolynomial difference(RationalPolynomial lhs,
                              const RationalPolynomial& rhs)
{
    add_in_place(lhs, scaled(rhs, -1));
    normalize(lhs);
    return lhs;
}

FieldPolynomial2 field_constant(const MatQuadraticFieldValue2& value)
{
    return {{value.rational}, {value.radical}};
}

FieldPolynomial2 field_coordinate(const std::vector<CORE::BigRat>& rational,
                                  const std::vector<CORE::BigRat>& radical)
{
    return {rational, radical};
}

FieldPolynomial2 field_add(FieldPolynomial2 lhs, const FieldPolynomial2& rhs)
{
    add_in_place(lhs.rational, rhs.rational);
    add_in_place(lhs.radical, rhs.radical);
    normalize(lhs.rational);
    normalize(lhs.radical);
    return lhs;
}

FieldPolynomial2 field_subtract(FieldPolynomial2 lhs,
                                const FieldPolynomial2& rhs)
{
    lhs.rational = difference(lhs.rational, rhs.rational);
    lhs.radical = difference(lhs.radical, rhs.radical);
    return lhs;
}

FieldPolynomial2 field_scale(FieldPolynomial2 value, const CORE::BigRat& factor)
{
    value.rational = scaled(std::move(value.rational), factor);
    value.radical = scaled(std::move(value.radical), factor);
    return value;
}

FieldPolynomial2 field_multiply(const FieldPolynomial2& lhs,
                                const FieldPolynomial2& rhs,
                                const CORE::BigRat& radicand)
{
    RationalPolynomial rational = multiply(lhs.rational, rhs.rational);
    add_in_place(rational,
                 scaled(multiply(lhs.radical, rhs.radical), radicand));
    RationalPolynomial radical = multiply(lhs.rational, rhs.radical);
    add_in_place(radical, multiply(lhs.radical, rhs.rational));
    normalize(rational);
    normalize(radical);
    return {
        std::move(rational),
        std::move(radical),
    };
}

FieldPolynomial2 field_square(const FieldPolynomial2& value,
                              const CORE::BigRat& radicand)
{
    return field_multiply(value, value, radicand);
}

bool is_zero(const RationalPolynomial& polynomial)
{
    RationalPolynomial normalized = polynomial;
    normalize(normalized);
    return normalized.size() == 1 && normalized.front() == 0;
}

Polynomial integer_polynomial(const RationalPolynomial& polynomial)
{
    const std::vector<ExactAlgebraicInteger1> coefficients =
        primitive_integer_coefficients(polynomial);
    return
        typename CGAL::Polynomial_traits_d<Polynomial>::Construct_polynomial()(
            coefficients.begin(), coefficients.end());
}

CGAL::Sign rational_sign_at(const RationalPolynomial& polynomial,
                            const Algebraic& root,
                            const ExactAlgebraicKernel1& kernel)
{
    if (is_zero(polynomial))
    {
        return CGAL::ZERO;
    }
    const CGAL::Sign normalized_sign =
        kernel.sign_at_1_object()(integer_polynomial(polynomial), root);
    return polynomial.back() < 0 ? CGAL::opposite(normalized_sign)
                                 : normalized_sign;
}

CGAL::Sign field_sign_at(const FieldPolynomial2& value,
                         const CORE::BigRat& radicand, const Algebraic& root,
                         const ExactAlgebraicKernel1& kernel)
{
    const CGAL::Sign rational_sign =
        rational_sign_at(value.rational, root, kernel);
    const CGAL::Sign radical_sign =
        rational_sign_at(value.radical, root, kernel);
    if (radical_sign == CGAL::ZERO)
    {
        return rational_sign;
    }
    if (rational_sign == CGAL::ZERO)
    {
        return radical_sign;
    }
    if (rational_sign == radical_sign)
    {
        return rational_sign;
    }
    RationalPolynomial magnitude = multiply(value.rational, value.rational);
    add_in_place(magnitude,
                 scaled(multiply(value.radical, value.radical), -radicand));
    const CGAL::Sign magnitude_sign = rational_sign_at(magnitude, root, kernel);
    if (magnitude_sign == CGAL::ZERO)
    {
        return CGAL::ZERO;
    }
    return magnitude_sign == CGAL::POSITIVE ? rational_sign : radical_sign;
}

FieldPolynomial2 squared_distance(const FieldPolynomial2& x,
                                  const FieldPolynomial2& y,
                                  const MatExactPointSiteSource2& point,
                                  const CORE::BigRat& radicand)
{
    return field_add(
        field_square(field_subtract(x, field_constant(point.x)), radicand),
        field_square(field_subtract(y, field_constant(point.y)), radicand));
}

FieldPolynomial2 dot(const FieldPolynomial2& ax, const FieldPolynomial2& ay,
                     const FieldPolynomial2& bx, const FieldPolynomial2& by,
                     const CORE::BigRat& radicand)
{
    return field_add(field_multiply(ax, bx, radicand),
                     field_multiply(ay, by, radicand));
}

CORE::Expr core_field_value(const MatQuadraticFieldValue2& value,
                            const CORE::BigRat& radicand)
{
    return CORE::Expr(value.rational) +
           CORE::sqrt(CORE::Expr(radicand)) * CORE::Expr(value.radical);
}

CORE::Expr core_evaluate(const FieldPolynomial2& value,
                         const CORE::BigRat& radicand,
                         const CORE::Expr& parameter)
{
    const auto evaluate = [&parameter](const RationalPolynomial& polynomial)
    {
        CORE::Expr result(0);
        for (auto coefficient = polynomial.rbegin();
             coefficient != polynomial.rend(); ++coefficient)
        {
            result *= parameter;
            result += CORE::Expr(*coefficient);
        }
        return result;
    };
    return evaluate(value.rational) +
           CORE::sqrt(CORE::Expr(radicand)) * evaluate(value.radical);
}

CORE::Expr exact_core_root(const ExactFactorRootWitness2& witness)
{
    return CORE::rootOf(
        CORE::Polynomial<CORE::BigInt>(witness.primitive_factor),
        static_cast<int>(witness.factor_local_ordinal + 1));
}

std::vector<ExactAlgebraicInteger1>
factor_coefficients(const Polynomial& factor)
{
    const int degree = CGAL::degree(factor);
    std::vector<ExactAlgebraicInteger1> coefficients;
    coefficients.reserve(static_cast<std::size_t>(degree + 1));
    for (int index = 0; index <= degree; ++index)
    {
        coefficients.push_back(factor[index]);
    }
    return coefficients;
}

std::vector<ExactFactorRootWitness2>
factor_root_witnesses(const Polynomial& polynomial,
                      const ExactAlgebraicKernel1& kernel)
{
    std::vector<std::pair<Polynomial, ExactAlgebraicKernel1::Multiplicity_type>>
        raw_factors;
    kernel.square_free_factorize_1_object()(polynomial,
                                            std::back_inserter(raw_factors));
    std::vector<ExactFactorRootWitness2> witnesses;
    for (const auto& [raw_factor, multiplicity] : raw_factors)
    {
        const Polynomial factor =
            typename CGAL::Polynomial_traits_d<Polynomial>::Canonicalize()(
                raw_factor);
        const std::vector<ExactAlgebraicInteger1> coefficients =
            factor_coefficients(factor);
        std::vector<Algebraic> roots;
        kernel.solve_1_object()(factor, true, std::back_inserter(roots));
        for (std::size_t ordinal = 0; ordinal < roots.size(); ++ordinal)
        {
            witnesses.push_back({
                coefficients,
                ordinal,
                roots[ordinal],
                algebraic_root_id_v1(coefficients, ordinal),
                static_cast<std::size_t>(multiplicity),
            });
        }
    }
    return witnesses;
}

void require_same_field(const MatExactPointSiteSource2& point,
                        const CORE::BigRat& radicand)
{
    if (point.radicand != radicand)
    {
        throw InvalidRationalPrimitiveError(
            "endpoint sources use different quadratic fields");
    }
}

bool has_radical_coefficient(
    const MatExactPointSiteSource2& point)
{
    return point.x.radical != 0
        || point.y.radical != 0;
}

bool rational_is_square(
    const CORE::BigRat& value)
{
    using FractionTraits =
        CGAL::Fraction_traits<CORE::BigRat>;
    ExactAlgebraicInteger1 numerator;
    ExactAlgebraicInteger1 denominator;
    typename FractionTraits::Decompose()(
        value,
        numerator,
        denominator);
    ExactAlgebraicInteger1 numerator_root;
    ExactAlgebraicInteger1 denominator_root;
    return CGAL::is_square(
               numerator,
               numerator_root)
        && CGAL::is_square(
               denominator,
               denominator_root);
}

template <typename... Points>
void require_canonical_quadratic_field(
    const MatExactPointSiteSource2& focus,
    const Points&... points)
{
    const CORE::BigRat radicand = focus.radicand;
    if (radicand <= 0)
    {
        throw NonCanonicalQuadraticFieldError(
            "quadratic-field radicand must be positive");
    }
    (require_same_field(points, radicand), ...);
    if (!rational_is_square(radicand))
    {
        return;
    }
    if (has_radical_coefficient(focus)
        || (has_radical_coefficient(points) || ...))
    {
        throw NonCanonicalQuadraticFieldError(
            "square radicand requires zero radical coefficients");
    }
}

bool exact_point_records_equal(
    const MatExactPointSiteSource2& lhs,
    const MatExactPointSiteSource2& rhs)
{
    return lhs.radicand == rhs.radicand
        && lhs.x.rational == rhs.x.rational
        && lhs.x.radical == rhs.x.radical
        && lhs.y.rational == rhs.y.rational
        && lhs.y.radical == rhs.y.radical;
}

void require_distinct_segment_endpoints(
    const MatExactPointSiteSource2& source,
    const MatExactPointSiteSource2& target,
    const char* role)
{
    if (exact_point_records_equal(source, target))
    {
        throw CoincidentSegmentEndpointsError(
            std::string(role)
            + " segment endpoints coincide");
    }
}

void require_nonzero_line_normal(
    const MatExactOpenSegmentSource2& segment,
    const char* role)
{
    if (segment.line_a == 0
        && segment.line_b == 0)
    {
        throw ZeroSegmentLineNormalError(
            std::string(role)
            + " segment has a zero line normal");
    }
}

FieldPolynomial2 line_value(const MatExactPointSiteSource2& point,
                            const MatExactOpenSegmentSource2& segment)
{
    return field_add(
        field_add(field_scale(field_constant(point.x), segment.line_a),
                  field_scale(field_constant(point.y), segment.line_b)),
        field_constant({segment.line_c, 0}));
}

FieldPolynomial2 line_value(const FieldPolynomial2& x,
                            const FieldPolynomial2& y,
                            const MatExactOpenSegmentSource2& segment)
{
    return field_add(
        field_add(field_scale(x, segment.line_a),
                  field_scale(y, segment.line_b)),
        field_constant({segment.line_c, 0}));
}

RationalPolynomial field_norm(const FieldPolynomial2& value,
                              const CORE::BigRat& radicand)
{
    RationalPolynomial norm = multiply(value.rational, value.rational);
    add_in_place(
        norm,
        scaled(multiply(value.radical, value.radical), -radicand));
    normalize(norm);
    return norm;
}

struct FieldZeroSetPolynomial2
{
    RationalPolynomial polynomial;
    std::string factor_kind;
};

FieldZeroSetPolynomial2 field_zero_set_polynomial(
    const FieldPolynomial2& value,
    const CORE::BigRat& radicand)
{
    if (is_zero(value.radical))
    {
        RationalPolynomial rational =
            value.rational;
        normalize(rational);
        return {
            std::move(rational),
            "equation-factor-multiplicity",
        };
    }
    if (is_zero(value.rational))
    {
        RationalPolynomial radical =
            value.radical;
        normalize(radical);
        return {
            std::move(radical),
            "equation-factor-multiplicity",
        };
    }
    return {
        field_norm(
            value,
            radicand),
        "norm-factor-multiplicity",
    };
}

std::size_t polynomial_degree(
    const RationalPolynomial& polynomial)
{
    RationalPolynomial normalized = polynomial;
    normalize(normalized);
    return normalized.size() - 1;
}

FieldPolynomial2 segment_limiter_equation(
    const MatExactPointSiteSource2& focus,
    const MatExactOpenSegmentSource2& segment,
    const MatExactOpenSegmentSource2& limiter)
{
    const SourceParabolaParameterization2 source =
        source_parameterization(focus, segment);
    const FieldPolynomial2 x =
        field_coordinate(source.x_rational, source.x_radical);
    const FieldPolynomial2 y =
        field_coordinate(source.y_rational, source.y_radical);
    const CORE::BigRat limiter_line_norm =
        limiter.line_a * limiter.line_a
        + limiter.line_b * limiter.line_b;
    require_nonzero_line_normal(limiter, "limiter");
    const FieldPolynomial2 equation =
        field_subtract(
            field_scale(
                squared_distance(x, y, focus, source.radicand),
                limiter_line_norm),
            field_square(
                line_value(x, y, limiter),
                source.radicand));
    if (polynomial_degree(equation.rational) > 4
        || polynomial_degree(equation.radical) > 4)
    {
        throw InvalidRationalPrimitiveError(
            "segment limiter equation exceeds quartic degree");
    }
    return equation;
}

MatTraits::Point_2 exact_live_point(const MatExactPointSiteSource2& point,
                                    const CORE::BigRat& radicand)
{
    return {
        core_field_value(point.x, radicand),
        core_field_value(point.y, radicand),
    };
}

bool exact_site_equal(const MatTraits::Site_2& lhs,
                      const MatTraits::Site_2& rhs)
{
    if (lhs.is_point() != rhs.is_point())
    {
        return false;
    }
    if (lhs.is_point())
    {
        return lhs.point() == rhs.point();
    }
    return (lhs.source() == rhs.source() && lhs.target() == rhs.target()) ||
           (lhs.source() == rhs.target() && lhs.target() == rhs.source());
}

LiveParabolaEndpointBridge2
live_endpoint_bridge(const MatExactPointSiteSource2& focus,
                     const MatExactPointSiteSource2& segment_source,
                     const MatExactPointSiteSource2& segment_target,
                     const MatTraits::Site_2& expected_limiter,
                     const CORE::BigRat& radicand,
                     const SegmentSiteVoronoi2& voronoi,
                     const SegmentSiteVoronoi2::Halfedge_handle& halfedge)
{
    const MatTraits::Site_2 expected_focus =
        MatTraits::Site_2::construct_site_2(exact_live_point(focus, radicand));
    const MatTraits::Site_2 expected_segment =
        MatTraits::Site_2::construct_site_2(
            exact_live_point(segment_source, radicand),
            exact_live_point(segment_target, radicand));
    const MatTraits::Site_2 live_up_generator = halfedge->up()->site();
    const MatTraits::Site_2 live_down_generator = halfedge->down()->site();
    const bool generators_match =
        (exact_site_equal(live_up_generator, expected_focus) &&
         exact_site_equal(live_down_generator, expected_segment)) ||
        (exact_site_equal(live_up_generator, expected_segment) &&
         exact_site_equal(live_down_generator, expected_focus));
    if (!generators_match)
    {
        throw MismatchedLiveParabolaBridgeError(
            "source records do not match live up/down generators");
    }
    const bool limiter_is_left =
        halfedge->has_source()
        && exact_site_equal(
            halfedge->left()->site(),
            expected_limiter);
    const bool limiter_is_right =
        halfedge->has_target()
        && exact_site_equal(
            halfedge->right()->site(),
            expected_limiter);
    if (!limiter_is_left && !limiter_is_right)
    {
        throw MismatchedLiveParabolaBridgeError(
            "limiter record matches neither live left nor right site");
    }
    if (limiter_is_left && limiter_is_right)
    {
        throw AmbiguousLiveParabolaEndpointError(
            "limiter record matches both live left and right sites");
    }
    if ((limiter_is_left && !halfedge->has_source()) ||
        (limiter_is_right && !halfedge->has_target()))
    {
        throw MismatchedLiveParabolaBridgeError(
            "limiter-owned endpoint is unbounded");
    }
    SegmentSiteParabola2 parabola;
    if (!CGAL::assign(parabola, voronoi.dual().primal(halfedge->dual())))
    {
        throw MismatchedLiveParabolaBridgeError(
            "live halfedge dual is not a parabola");
    }
    require_distinct_live_parabola_endpoints(parabola);
    const MatTraits::Point_2 point = limiter_is_left
                                         ? halfedge->source()->point()
                                         : halfedge->target()->point();
    if (point != parabola.p1 && point != parabola.p2)
    {
        throw MismatchedLiveParabolaBridgeError(
            "limiter-owned vertex is not a live parabola endpoint");
    }
    return {
        point,
        std::move(parabola),
    };
}

bool exact_segment_feature_contains(
    const SourceParabolaParameterization2& parabola,
    const MatExactPointSiteSource2& segment_source,
    const MatExactPointSiteSource2& segment_target,
    const Algebraic& parameter,
    bool include_endpoints)
{
    require_same_field(segment_source, parabola.radicand);
    require_same_field(segment_target, parabola.radicand);
    const FieldPolynomial2 x =
        field_coordinate(parabola.x_rational, parabola.x_radical);
    const FieldPolynomial2 y =
        field_coordinate(parabola.y_rational, parabola.y_radical);
    const FieldPolynomial2 segment_dx = field_subtract(
        field_constant(segment_target.x),
        field_constant(segment_source.x));
    const FieldPolynomial2 segment_dy = field_subtract(
        field_constant(segment_target.y),
        field_constant(segment_source.y));
    const FieldPolynomial2 segment_length_squared =
        dot(
            segment_dx,
            segment_dy,
            segment_dx,
            segment_dy,
            parabola.radicand);
    ExactAlgebraicKernel1 kernel;
    if (field_sign_at(
            segment_length_squared,
            parabola.radicand,
            parameter,
            kernel) != CGAL::POSITIVE)
    {
        throw InvalidRationalPrimitiveError(
            "segment feature endpoints coincide");
    }
    const FieldPolynomial2 projection =
        dot(
            field_subtract(
                x,
                field_constant(segment_source.x)),
            field_subtract(
                y,
                field_constant(segment_source.y)),
            segment_dx,
            segment_dy,
            parabola.radicand);
    const FieldPolynomial2 projection_remainder =
        field_subtract(segment_length_squared, projection);
    const CGAL::Sign projection_sign =
        field_sign_at(
            projection,
            parabola.radicand,
            parameter,
            kernel);
    const CGAL::Sign remainder_sign =
        field_sign_at(
            projection_remainder,
            parabola.radicand,
            parameter,
            kernel);
    if (include_endpoints)
    {
        return projection_sign != CGAL::NEGATIVE
            && remainder_sign != CGAL::NEGATIVE;
    }
    return projection_sign == CGAL::POSITIVE
        && remainder_sign == CGAL::POSITIVE;
}

} // namespace

bool exact_open_segment_feature_contains(
    const SourceParabolaParameterization2& parabola,
    const MatExactPointSiteSource2& segment_source,
    const MatExactPointSiteSource2& segment_target, const Algebraic& parameter)
{
    return exact_segment_feature_contains(
        parabola,
        segment_source,
        segment_target,
        parameter,
        false);
}

void require_distinct_live_parabola_endpoints(
    const SegmentSiteParabola2& parabola)
{
    if (parabola.p1 == parabola.p2)
    {
        throw AmbiguousLiveParabolaEndpointError(
            "live parabola endpoints collapse to one point");
    }
}

MatParameterEndpoint2 bind_point_limiter_parabola_endpoint(
    const MatExactPointSiteSource2& focus,
    const MatExactOpenSegmentSource2& segment,
    const MatExactPointSiteSource2& segment_source,
    const MatExactPointSiteSource2& segment_target,
    const MatExactPointSiteSource2& limiter,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge)
{
    require_canonical_quadratic_field(
        focus,
        segment_source,
        segment_target,
        limiter);
    const CORE::BigRat radicand = focus.radicand;
    require_nonzero_line_normal(segment, "source");
    require_distinct_segment_endpoints(
        segment_source,
        segment_target,
        "source");
    const MatTraits::Site_2 expected_limiter =
        MatTraits::Site_2::construct_site_2(
            exact_live_point(limiter, radicand));
    const LiveParabolaEndpointBridge2 live =
        live_endpoint_bridge(
            focus,
            segment_source,
            segment_target,
            expected_limiter,
            radicand,
            voronoi,
            halfedge);
    if (quadratic_field_sign(
            {
                line_value(segment_source, segment).rational.front(),
                line_value(segment_source, segment).radical.front(),
            },
            radicand) != CGAL::ZERO ||
        quadratic_field_sign(
            {
                line_value(segment_target, segment).rational.front(),
                line_value(segment_target, segment).radical.front(),
            },
            radicand) != CGAL::ZERO)
    {
        throw InvalidRationalPrimitiveError(
            "open-segment endpoint is off directrix");
    }

    const SourceParabolaParameterization2 source =
        source_parameterization(focus, segment);
    const FieldPolynomial2 x =
        field_coordinate(source.x_rational, source.x_radical);
    const FieldPolynomial2 y =
        field_coordinate(source.y_rational, source.y_radical);
    const FieldPolynomial2 equation =
        field_subtract(squared_distance(x, y, focus, radicand),
                       squared_distance(x, y, limiter, radicand));
    RationalPolynomial norm = field_norm(equation, radicand);
    if (is_zero(norm))
    {
        throw UnboundLiveParabolaEndpointError(
            "point limiter equality is not zero-dimensional");
    }
    const Polynomial norm_polynomial = integer_polynomial(norm);
    ExactAlgebraicKernel1 kernel;
    const std::vector<ExactFactorRootWitness2> witnesses =
        factor_root_witnesses(norm_polynomial, kernel);

    const CORE::Expr vertex_x = core_field_value(
        {
            source.x_rational.front(),
            source.x_radical.front(),
        },
        radicand);
    const CORE::Expr vertex_y = core_field_value(
        {
            source.y_rational.front(),
            source.y_radical.front(),
        },
        radicand);
    const CORE::BigRat line_norm =
        segment.line_a * segment.line_a + segment.line_b * segment.line_b;
    const CORE::Expr live_parameter =
        ((live.point.x() - vertex_x) * CORE::Expr(-segment.line_b) +
         (live.point.y() - vertex_y) * CORE::Expr(segment.line_a)) /
        CORE::Expr(line_norm);

    std::vector<MatParameterEndpoint2> matches;
    std::size_t equation_roots = 0;
    std::size_t feature_roots = 0;
    std::size_t parameter_matches = 0;
    std::size_t point_matches = 0;
    for (const ExactFactorRootWitness2& witness : witnesses)
    {
        const Algebraic& root = witness.root;
        if (field_sign_at(equation, radicand, root, kernel) != CGAL::ZERO)
        {
            continue;
        }
        ++equation_roots;
        if (!exact_open_segment_feature_contains(source, segment_source,
                                                 segment_target, root))
        {
            continue;
        }
        ++feature_roots;
        const CORE::Expr parameter = exact_core_root(witness);
        if (parameter != live_parameter)
        {
            continue;
        }
        ++parameter_matches;
        const MatTraits::Point_2 candidate(
            core_evaluate(x, radicand, parameter),
            core_evaluate(y, radicand, parameter));
        if (candidate != live.point)
        {
            continue;
        }
        ++point_matches;
        matches.push_back({
            root,
            {
                witness.root_id,
                limiter.stable_site_id,
                "point-limiter/norm-factor-multiplicity/" +
                    std::to_string(witness.source_factor_multiplicity),
            },
        });
    }
    if (matches.size() != 1)
    {
        throw UnboundLiveParabolaEndpointError(
            "live point binds " + std::to_string(matches.size()) +
            " roots from " + std::to_string(witnesses.size()) +
            " square-free factor roots, " + std::to_string(equation_roots) +
            " radical, " + std::to_string(feature_roots) + " feature, " +
            std::to_string(parameter_matches) + " parameter, " +
            std::to_string(point_matches) + " point matches");
    }
    return std::move(matches.front());
}

namespace
{

void require_valid_segment_limiter_binding_inputs(
    const MatExactPointSiteSource2& focus,
    const MatExactOpenSegmentSource2& segment,
    const MatExactPointSiteSource2& segment_source,
    const MatExactPointSiteSource2& segment_target,
    const MatExactOpenSegmentSource2& limiter,
    const MatExactPointSiteSource2& limiter_source,
    const MatExactPointSiteSource2& limiter_target)
{
    require_canonical_quadratic_field(
        focus,
        segment_source,
        segment_target,
        limiter_source,
        limiter_target);
    const CORE::BigRat radicand = focus.radicand;
    require_nonzero_line_normal(segment, "source");
    require_nonzero_line_normal(limiter, "limiter");
    require_distinct_segment_endpoints(
        segment_source,
        segment_target,
        "source");
    require_distinct_segment_endpoints(
        limiter_source,
        limiter_target,
        "limiter");
    if (quadratic_field_sign(
            {
                line_value(segment_source, segment).rational.front(),
                line_value(segment_source, segment).radical.front(),
            },
            radicand) != CGAL::ZERO
        || quadratic_field_sign(
               {
                   line_value(segment_target, segment).rational.front(),
                   line_value(segment_target, segment).radical.front(),
               },
               radicand) != CGAL::ZERO)
    {
        throw InvalidRationalPrimitiveError(
            "open-segment endpoint is off directrix");
    }
    if (quadratic_field_sign(
            {
                line_value(limiter_source, limiter).rational.front(),
                line_value(limiter_source, limiter).radical.front(),
            },
            radicand) != CGAL::ZERO
        || quadratic_field_sign(
               {
                   line_value(limiter_target, limiter).rational.front(),
                   line_value(limiter_target, limiter).radical.front(),
               },
               radicand) != CGAL::ZERO)
    {
        throw InvalidRationalPrimitiveError(
            "segment limiter endpoint is off supporting line");
    }
}

MatParameterEndpoint2
bind_segment_limiter_parabola_endpoint_from_live(
    const MatExactPointSiteSource2& focus,
    const MatExactOpenSegmentSource2& segment,
    const MatExactPointSiteSource2& segment_source,
    const MatExactPointSiteSource2& segment_target,
    const MatExactOpenSegmentSource2& limiter,
    const MatExactPointSiteSource2& limiter_source,
    const MatExactPointSiteSource2& limiter_target,
    const LiveParabolaEndpointBridge2& live)
{
    const CORE::BigRat radicand = focus.radicand;
    const SourceParabolaParameterization2 source =
        source_parameterization(focus, segment);
    const FieldPolynomial2 x =
        field_coordinate(source.x_rational, source.x_radical);
    const FieldPolynomial2 y =
        field_coordinate(source.y_rational, source.y_radical);
    const FieldPolynomial2 equation =
        segment_limiter_equation(
            focus,
            segment,
            limiter);
    const RationalPolynomial norm =
        field_norm(equation, radicand);
    if (is_zero(norm))
    {
        // Unreachable for validated adaptor edges; protects producer drift.
        throw IdenticallyZeroSegmentLimiterEquationError(
            "segment limiter equality is identically zero");
    }
    const Polynomial norm_polynomial =
        integer_polynomial(norm);
    ExactAlgebraicKernel1 kernel;
    const std::vector<ExactFactorRootWitness2> witnesses =
        factor_root_witnesses(norm_polynomial, kernel);

    const CORE::Expr vertex_x = core_field_value(
        {
            source.x_rational.front(),
            source.x_radical.front(),
        },
        radicand);
    const CORE::Expr vertex_y = core_field_value(
        {
            source.y_rational.front(),
            source.y_radical.front(),
        },
        radicand);
    const CORE::BigRat line_norm =
        segment.line_a * segment.line_a
        + segment.line_b * segment.line_b;
    const CORE::Expr live_parameter =
        ((live.point.x() - vertex_x)
             * CORE::Expr(-segment.line_b)
         + (live.point.y() - vertex_y)
             * CORE::Expr(segment.line_a))
        / CORE::Expr(line_norm);

    std::vector<MatParameterEndpoint2> matches;
    std::size_t equation_roots = 0;
    std::size_t source_feature_roots = 0;
    std::size_t limiter_feature_roots = 0;
    std::size_t parameter_matches = 0;
    std::size_t point_matches = 0;
    for (const ExactFactorRootWitness2& witness : witnesses)
    {
        const Algebraic& root = witness.root;
        if (field_sign_at(
                equation,
                radicand,
                root,
                kernel) != CGAL::ZERO)
        {
            continue;
        }
        ++equation_roots;
        if (!exact_segment_feature_contains(
                source,
                segment_source,
                segment_target,
                root,
                false))
        {
            continue;
        }
        ++source_feature_roots;
        if (!exact_segment_feature_contains(
                source,
                limiter_source,
                limiter_target,
                root,
                true))
        {
            continue;
        }
        ++limiter_feature_roots;
        const CORE::Expr parameter =
            exact_core_root(witness);
        if (parameter != live_parameter)
        {
            continue;
        }
        ++parameter_matches;
        const MatTraits::Point_2 candidate(
            core_evaluate(x, radicand, parameter),
            core_evaluate(y, radicand, parameter));
        if (candidate != live.point)
        {
            continue;
        }
        ++point_matches;
        matches.push_back({
            root,
            {
                witness.root_id,
                limiter.stable_site_id,
                "segment-limiter/norm-factor-multiplicity/"
                    + std::to_string(
                        witness.source_factor_multiplicity),
            },
        });
    }
    if (matches.size() != 1)
    {
        throw UnboundLiveParabolaEndpointError(
            "live segment-limiter point binds "
            + std::to_string(matches.size())
            + " roots from "
            + std::to_string(witnesses.size())
            + " square-free factor roots, "
            + std::to_string(equation_roots)
            + " equation, "
            + std::to_string(source_feature_roots)
            + " source feature, "
            + std::to_string(limiter_feature_roots)
            + " limiter feature, "
            + std::to_string(parameter_matches)
            + " parameter, "
            + std::to_string(point_matches)
            + " point matches");
    }
    return std::move(matches.front());
}

} // namespace

MatParameterEndpoint2 bind_segment_limiter_parabola_endpoint(
    const MatExactPointSiteSource2& focus,
    const MatExactOpenSegmentSource2& segment,
    const MatExactPointSiteSource2& segment_source,
    const MatExactPointSiteSource2& segment_target,
    const MatExactOpenSegmentSource2& limiter,
    const MatExactPointSiteSource2& limiter_source,
    const MatExactPointSiteSource2& limiter_target,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge)
{
    require_valid_segment_limiter_binding_inputs(
        focus,
        segment,
        segment_source,
        segment_target,
        limiter,
        limiter_source,
        limiter_target);
    const CORE::BigRat radicand = focus.radicand;
    const MatTraits::Site_2 expected_limiter =
        MatTraits::Site_2::construct_site_2(
            exact_live_point(limiter_source, radicand),
            exact_live_point(limiter_target, radicand));
    const LiveParabolaEndpointBridge2 live =
        live_endpoint_bridge(
            focus,
            segment_source,
            segment_target,
            expected_limiter,
            radicand,
            voronoi,
            halfedge);
    return bind_segment_limiter_parabola_endpoint_from_live(
        focus,
        segment,
        segment_source,
        segment_target,
        limiter,
        limiter_source,
        limiter_target,
        live);
}

MatParameterEndpoint2 bind_segment_limiter_parabola_endpoint(
    const MatExactPointSiteSource2& focus,
    const MatExactOpenSegmentSource2& segment,
    const MatExactPointSiteSource2& segment_source,
    const MatExactPointSiteSource2& segment_target,
    const MatExactOpenSegmentSource2& limiter,
    const MatExactPointSiteSource2& limiter_source,
    const MatExactPointSiteSource2& limiter_target,
    const SegmentSiteParabola2& live_parabola,
    const MatTraits::Point_2& live_point)
{
    require_valid_segment_limiter_binding_inputs(
        focus,
        segment,
        segment_source,
        segment_target,
        limiter,
        limiter_source,
        limiter_target);
    require_distinct_live_parabola_endpoints(
        live_parabola);
    if (live_point != live_parabola.p1
        && live_point != live_parabola.p2)
    {
        throw MismatchedLiveParabolaBridgeError(
            "explicit live point is not a parabola endpoint");
    }
    return bind_segment_limiter_parabola_endpoint_from_live(
        focus,
        segment,
        segment_source,
        segment_target,
        limiter,
        limiter_source,
        limiter_target,
        {live_point, live_parabola});
}

namespace
{

const CanonicalMatRationalPointSource2&
canonical_parabola_point_source(
    const CanonicalMatRationalSources2& sources,
    const std::string& stable_id)
{
    const auto found = std::lower_bound(
        sources.points().begin(),
        sources.points().end(),
        stable_id,
        [](const CanonicalMatRationalPointSource2& point,
           const std::string& identity)
        {
            return point.stable_site_id < identity;
        });
    if (found == sources.points().end()
        || found->stable_site_id != stable_id)
    {
        throw UnknownCanonicalMatParabolaSourceError(
            "canonical P-S point has no rational source");
    }
    return *found;
}

const CanonicalMatRationalOpenSegmentSource2&
canonical_parabola_segment_source(
    const CanonicalMatRationalSources2& sources,
    const std::string& stable_id)
{
    const auto found = std::lower_bound(
        sources.segments().begin(),
        sources.segments().end(),
        stable_id,
        [](const CanonicalMatRationalOpenSegmentSource2& segment,
           const std::string& identity)
        {
            return segment.stable_site_id < identity;
        });
    if (found == sources.segments().end()
        || found->stable_site_id != stable_id)
    {
        throw UnknownCanonicalMatParabolaSourceError(
            "canonical P-S segment has no rational source");
    }
    return *found;
}

MatExactPointSiteSource2 exact_canonical_parabola_point(
    const CanonicalMatRationalPointSource2& point)
{
    return {
        point.stable_site_id,
        {point.x, 0},
        {point.y, 0},
        1,
    };
}

std::pair<RationalDomainRoot2, RationalDomainRoot2>
canonical_parabola_feature_bounds(
    const MatExactPointSiteSource2& focus,
    const CanonicalMatRationalOpenSegmentSource2& segment,
    const CanonicalMatRationalSources2& sources)
{
    if (focus.x.radical != 0
        || focus.y.radical != 0)
    {
        throw InvalidCanonicalMatParabolaFeatureDomainError(
            "canonical P-S focus is not rational");
    }
    const CORE::BigRat line_norm =
        segment.support.line_a * segment.support.line_a
        + segment.support.line_b * segment.support.line_b;
    if (line_norm == 0)
    {
        throw InvalidCanonicalMatParabolaFeatureDomainError(
            "canonical P-S source segment has zero line normal");
    }
    const CORE::BigRat signed_distance =
        segment.support.line_a * focus.x.rational
        + segment.support.line_b * focus.y.rational
        + segment.support.line_c;
    if (signed_distance == 0)
    {
        throw InvalidCanonicalMatParabolaFeatureDomainError(
            "canonical P-S focus lies on its directrix");
    }
    const CORE::BigRat projection_x =
        focus.x.rational
        - signed_distance * segment.support.line_a
            / line_norm;
    const CORE::BigRat projection_y =
        focus.y.rational
        - signed_distance * segment.support.line_b
            / line_norm;
    const auto endpoint =
        [&segment,
         &sources,
         &projection_x,
         &projection_y,
         &line_norm](
            const std::string& point_id)
        {
            const auto& point =
                canonical_parabola_point_source(
                    sources,
                    point_id);
            return RationalDomainRoot2{
                (
                    (point.x - projection_x)
                        * -segment.support.line_b
                    + (point.y - projection_y)
                        * segment.support.line_a)
                    / line_norm,
                {point.stable_site_id},
            };
        };
    RationalDomainRoot2 source =
        endpoint(segment.source_point_id);
    RationalDomainRoot2 target =
        endpoint(segment.target_point_id);
    if (source.parameter == target.parameter)
    {
        throw InvalidCanonicalMatParabolaFeatureDomainError(
            "canonical P-S feature interval is empty");
    }
    if (target.parameter < source.parameter)
    {
        std::swap(source, target);
    }
    return {
        std::move(source),
        std::move(target),
    };
}

MatTraits::Point_2 canonical_parabola_chart_point(
    const SourceParabolaParameterization2& primitive,
    const CORE::BigRat& parameter)
{
    return {
        core_evaluate(
            {
                primitive.x_rational,
                primitive.x_radical,
            },
            primitive.radicand,
            CORE::Expr(parameter)),
        core_evaluate(
            {
                primitive.y_rational,
                primitive.y_radical,
            },
            primitive.radicand,
            CORE::Expr(parameter)),
    };
}

MatParameterEndpoint2 bind_canonical_parabola_feature_endpoint(
    const SourceParabolaParameterization2& primitive,
    const RationalDomainRoot2& bound,
    const std::string& owner_id,
    const MatTraits::Point_2& live_point)
{
    if (std::find(
            bound.provenance_ids.begin(),
            bound.provenance_ids.end(),
            owner_id)
            == bound.provenance_ids.end()
        || canonical_parabola_chart_point(
               primitive,
               bound.parameter)
            != live_point)
    {
        throw MismatchedLiveParabolaBridgeError(
            "canonical P-S feature owner differs from live endpoint");
    }
    ExactAlgebraicKernel1 kernel;
    return exact_graph_endpoint_binding(
        {
            kernel.construct_algebraic_real_1_object()(
                bound.parameter),
            bound.provenance_ids,
        });
}

MatParameterEndpoint2 bind_canonical_parabola_owner(
    const MatTraits::Site_2& owner,
    const MatTraits::Point_2& live_point,
    const SourceParabolaParameterization2& primitive,
    const RationalDomainRoot2& lower,
    const RationalDomainRoot2& upper,
    const MatExactPointSiteSource2& focus,
    const CanonicalMatRationalOpenSegmentSource2& segment,
    const MatExactPointSiteSource2& segment_source,
    const MatExactPointSiteSource2& segment_target,
    const CanonicalMatRationalSources2& sources,
    const CanonicalMatSiteGeometryIndex2& site_index,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge)
{
    const std::string owner_id =
        site_index.stable_id(owner);
    if (std::find(
            lower.provenance_ids.begin(),
            lower.provenance_ids.end(),
            owner_id)
        != lower.provenance_ids.end())
    {
        return bind_canonical_parabola_feature_endpoint(
            primitive,
            lower,
            owner_id,
            live_point);
    }
    if (std::find(
            upper.provenance_ids.begin(),
            upper.provenance_ids.end(),
            owner_id)
        != upper.provenance_ids.end())
    {
        return bind_canonical_parabola_feature_endpoint(
            primitive,
            upper,
            owner_id,
            live_point);
    }
    if (owner.is_point())
    {
        return bind_point_limiter_parabola_endpoint(
            focus,
            segment.support,
            segment_source,
            segment_target,
            exact_canonical_parabola_point(
                canonical_parabola_point_source(
                    sources,
                    owner_id)),
            voronoi,
            halfedge);
    }
    const auto& limiter =
        canonical_parabola_segment_source(
            sources,
            owner_id);
    return bind_segment_limiter_parabola_endpoint(
        focus,
        segment.support,
        segment_source,
        segment_target,
        limiter.support,
        exact_canonical_parabola_point(
            canonical_parabola_point_source(
                sources,
                limiter.source_point_id)),
        exact_canonical_parabola_point(
            canonical_parabola_point_source(
                sources,
                limiter.target_point_id)),
        voronoi,
        halfedge);
}

} // namespace

std::pair<MatParameterEndpoint2, MatParameterEndpoint2>
bind_point_segment_cell_endpoints(
    const CanonicalMatRationalSources2& sources,
    const std::string& focus_id,
    const std::string& segment_id,
    const CanonicalMatSiteGeometryIndex2& site_index,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge)
{
    const auto& focus_record =
        canonical_parabola_point_source(
            sources,
            focus_id);
    const auto& segment =
        canonical_parabola_segment_source(
            sources,
            segment_id);
    if (focus_id == segment.source_point_id
        || focus_id == segment.target_point_id)
    {
        throw IncidentCanonicalMatParabolaSourceError(
            "canonical P-S focus is a source-segment endpoint");
    }
    const MatExactPointSiteSource2 focus =
        exact_canonical_parabola_point(
            focus_record);
    const MatExactPointSiteSource2 segment_source =
        exact_canonical_parabola_point(
            canonical_parabola_point_source(
                sources,
                segment.source_point_id));
    const MatExactPointSiteSource2 segment_target =
        exact_canonical_parabola_point(
            canonical_parabola_point_source(
                sources,
                segment.target_point_id));
    const std::vector<std::string> generators =
        ordered_generator_site_ids(
            focus_id,
            segment_id);
    if (ordered_generator_site_ids(
            site_index.stable_id(
                halfedge->up()->site()),
            site_index.stable_id(
                halfedge->down()->site()))
            != generators)
    {
        throw MismatchedLiveParabolaBridgeError(
            "canonical P-S sources differ from live generators");
    }
    if (!halfedge->has_source()
        || !halfedge->has_target())
    {
        throw UnboundLiveParabolaEndpointError(
            "canonical P-S halfedge is not bounded");
    }
    SegmentSiteParabola2 live;
    if (!CGAL::assign(
            live,
            voronoi.dual().primal(
                halfedge->dual())))
    {
        throw MismatchedLiveParabolaBridgeError(
            "canonical P-S dual is not a parabola");
    }
    require_distinct_live_parabola_endpoints(
        live);
    const MatTraits::Point_2 live_source =
        halfedge->source()->point();
    const MatTraits::Point_2 live_target =
        halfedge->target()->point();
    if ((live_source != live.p1
            && live_source != live.p2)
        || (live_target != live.p1
            && live_target != live.p2))
    {
        throw MismatchedLiveParabolaBridgeError(
            "canonical P-S adaptor endpoints differ from the live parabola");
    }
    const SourceParabolaParameterization2 primitive =
        source_parameterization(
            focus,
            segment.support);
    const auto bounds =
        canonical_parabola_feature_bounds(
            focus,
            segment,
            sources);
    MatParameterEndpoint2 source =
        bind_canonical_parabola_owner(
            halfedge->left()->site(),
            live_source,
            primitive,
            bounds.first,
            bounds.second,
            focus,
            segment,
            segment_source,
            segment_target,
            sources,
            site_index,
            voronoi,
            halfedge);
    MatParameterEndpoint2 target =
        bind_canonical_parabola_owner(
            halfedge->right()->site(),
            live_target,
            primitive,
            bounds.first,
            bounds.second,
            focus,
            segment,
            segment_source,
            segment_target,
            sources,
            site_index,
            voronoi,
            halfedge);
    ExactAlgebraicKernel1 kernel;
    const CGAL::Comparison_result order =
        kernel.compare_1_object()(
            *source.parameter,
            *target.parameter);
    if (order == CGAL::LARGER)
    {
        std::swap(source, target);
    }
    else if (order != CGAL::SMALLER)
    {
        throw AmbiguousLiveParabolaEndpointError(
            "canonical P-S endpoints are not strictly ordered");
    }
    return {
        std::move(source),
        std::move(target),
    };
}

RationalPrimitiveParameterization2
bind_point_point_ray_parameterization(
    const CanonicalMatRationalSources2& sources,
    const std::vector<std::string>& generator_ids,
    const CanonicalMatSiteGeometryIndex2& site_index,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge)
{
    static_cast<void>(
        stable_dual_identity_v1(
            "point",
            generator_ids));
    const auto point_source =
        [&sources](
            const std::string& stable_id)
            -> const CanonicalMatRationalPointSource2&
        {
            const auto found = std::lower_bound(
                sources.points().begin(),
                sources.points().end(),
                stable_id,
                [](const CanonicalMatRationalPointSource2&
                       point,
                   const std::string& identity)
                {
                    return point.stable_site_id
                        < identity;
                });
            if (found == sources.points().end()
                || found->stable_site_id
                    != stable_id)
            {
                throw UnknownCanonicalMatPointRaySourceError(
                    "canonical P-P ray point has no rational source");
            }
            return *found;
        };
    const auto& first =
        point_source(generator_ids[0]);
    const auto& second =
        point_source(generator_ids[1]);
    const MatTraits::Site_2 up =
        halfedge->up()->site();
    const MatTraits::Site_2 down =
        halfedge->down()->site();
    if (!up.is_point()
        || !down.is_point()
        || ordered_generator_site_ids(
               site_index.stable_id(up),
               site_index.stable_id(down))
            != generator_ids)
    {
        throw MismatchedLivePointPointRayError(
            "canonical P-P ray sources differ from live generators");
    }
    if (halfedge->has_source()
            == halfedge->has_target())
    {
        throw UnboundLivePointPointRayEndpointError(
            "canonical P-P ray does not have exactly one finite endpoint");
    }

    MatTraits::Ray_2 live_ray;
    if (!CGAL::assign(
            live_ray,
            voronoi.dual().primal(
                halfedge->dual())))
    {
        throw MismatchedLivePointPointRayError(
            "canonical P-P live dual is not a ray");
    }
    const MatTraits::Point_2 live_endpoint =
        halfedge->has_source()
        ? halfedge->source()->point()
        : halfedge->target()->point();
    if (live_endpoint != live_ray.source())
    {
        throw MismatchedLivePointPointRayError(
            "canonical P-P adaptor endpoint differs from the live ray source");
    }
    const MatTraits::Point_2 first_point(
        CORE::Expr(first.x),
        CORE::Expr(first.y));
    const MatTraits::Point_2 second_point(
        CORE::Expr(second.x),
        CORE::Expr(second.y));
    if (CGAL::squared_distance(
            live_endpoint,
            first_point)
        != CGAL::squared_distance(
            live_endpoint,
            second_point))
    {
        throw MismatchedLivePointPointRayError(
            "canonical P-P ray source is not equidistant from its generators");
    }

    const auto rational_coordinate =
        [](const MatTraits::FT& coordinate)
        {
            try
            {
                return exact_mat_input_rational(
                    coordinate);
            }
            catch (
                const NonRationalCanonicalMatCoordinateError&)
            {
                throw NonRationalLivePointPointRayEndpointError(
                    "canonical P-P ray endpoint is not rational");
            }
        };
    const CORE::BigRat endpoint_x =
        rational_coordinate(
            live_endpoint.x());
    const CORE::BigRat endpoint_y =
        rational_coordinate(
            live_endpoint.y());
    const CORE::BigRat direction_x =
        first.y - second.y;
    const CORE::BigRat direction_y =
        second.x - first.x;
    if (direction_x == 0
        && direction_y == 0)
    {
        throw MismatchedLivePointPointRayError(
            "canonical P-P ray generators coincide");
    }
    const CORE::BigRat origin_x =
        (first.x + second.x) / 2;
    const CORE::BigRat origin_y =
        (first.y + second.y) / 2;
    const CORE::BigRat endpoint_parameter =
        direction_x != 0
        ? (endpoint_x - origin_x)
            / direction_x
        : (endpoint_y - origin_y)
            / direction_y;
    if (origin_x
                + direction_x
                    * endpoint_parameter
            != endpoint_x
        || origin_y
                + direction_y
                    * endpoint_parameter
            != endpoint_y)
    {
        throw MismatchedLivePointPointRayError(
            "canonical P-P ray endpoint is outside its bisector chart");
    }

    const MatKernel::Vector_2 live_direction =
        live_ray.to_vector();
    const MatTraits::FT cross =
        CORE::Expr(direction_x)
            * live_direction.y()
        - CORE::Expr(direction_y)
            * live_direction.x();
    const MatTraits::FT dot =
        CORE::Expr(direction_x)
            * live_direction.x()
        + CORE::Expr(direction_y)
            * live_direction.y();
    if (CGAL::sign(cross) != CGAL::ZERO
        || CGAL::sign(dot) == CGAL::ZERO)
    {
        throw MismatchedLivePointPointRayError(
            "canonical P-P chart and live ray direction disagree");
    }

    std::optional<CORE::BigRat> lower;
    std::optional<CORE::BigRat> upper;
    if (CGAL::sign(dot) == CGAL::POSITIVE)
    {
        lower = endpoint_parameter;
    }
    else
    {
        upper = endpoint_parameter;
    }
    return {
        {
            origin_x,
            direction_x,
        },
        {
            origin_y,
            direction_y,
        },
        std::move(lower),
        std::move(upper),
    };
}

namespace
{

void require_canonical_parallel_segment_chart(
    const RationalPrimitiveParameterization2& primitive,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const std::vector<std::string>& generator_ids)
{
    if (ordered_generator_site_ids(
            first_segment.stable_site_id,
            second_segment.stable_site_id)
            != generator_ids)
    {
        throw MismatchedLiveSegmentSegmentBridgeError(
            "parallel S-S source records do not match generator identities");
    }
    const RationalPrimitiveParameterization2 canonical =
        parallel_segment_bisector_parameterization(
            first_segment,
            second_segment);
    if (primitive.x_coefficients != canonical.x_coefficients
        || primitive.y_coefficients != canonical.y_coefficients)
    {
        throw MismatchedLiveSegmentSegmentBridgeError(
            "parallel S-S limiter chart is not canonical");
    }
}

RationalPolynomial parallel_squared_support_distance(
    const RationalPrimitiveParameterization2& primitive,
    const MatExactOpenSegmentSource2& segment)
{
    if (primitive.x_coefficients.size() != 2
        || primitive.y_coefficients.size() != 2)
    {
        throw InvalidRationalPrimitiveError(
            "parallel S-S point limiter requires an affine chart");
    }
    const CORE::BigRat line_norm =
        segment.line_a * segment.line_a
        + segment.line_b * segment.line_b;
    if (line_norm == 0)
    {
        throw ZeroSegmentLineNormalError(
            "parallel S-S point limiter source has zero line normal");
    }
    RationalPolynomial line_value =
        scaled(
            primitive.x_coefficients,
            segment.line_a);
    add_in_place(
        line_value,
        scaled(
            primitive.y_coefficients,
            segment.line_b));
    line_value.front() += segment.line_c;
    normalize(line_value);
    return scaled(
        multiply(line_value, line_value),
        CORE::BigRat(1) / line_norm);
}

RationalPolynomial parallel_support_clearance(
    const RationalPrimitiveParameterization2& primitive,
    const MatExactOpenSegmentSource2& segment)
{
    RationalPolynomial clearance =
        parallel_squared_support_distance(
            primitive,
            segment);
    normalize(clearance);
    if (clearance.size() != 1)
    {
        throw MismatchedLiveSegmentSegmentBridgeError(
            "parallel S-S chart is not parallel to its source supports");
    }
    return clearance;
}

void require_rational_parallel_point_limiter(
    const MatExactPointSiteSource2& limiter)
{
    if (limiter.radicand <= 0
        || limiter.x.radical != 0
        || limiter.y.radical != 0)
    {
        throw NonRationalParallelSegmentPointLimiterError(
            "parallel S-S point limiter is not rational");
    }
}

RationalPolynomial parallel_point_limiter_equation(
    const RationalPrimitiveParameterization2& primitive,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const MatExactPointSiteSource2& limiter)
{
    require_rational_parallel_point_limiter(
        limiter);
    RationalPolynomial first_clearance =
        parallel_support_clearance(
            primitive,
            first_segment);
    RationalPolynomial second_clearance =
        parallel_support_clearance(
            primitive,
            second_segment);
    normalize(first_clearance);
    normalize(second_clearance);
    if (first_clearance != second_clearance)
    {
        throw MismatchedLiveSegmentSegmentBridgeError(
            "parallel S-S point limiter sources have unequal clearance");
    }
    RationalPolynomial dx = primitive.x_coefficients;
    RationalPolynomial dy = primitive.y_coefficients;
    dx.front() -= limiter.x.rational;
    dy.front() -= limiter.y.rational;
    RationalPolynomial equation = multiply(dx, dx);
    add_in_place(equation, multiply(dy, dy));
    add_in_place(
        equation,
        scaled(first_clearance, -1));
    normalize(equation);
    if (polynomial_degree(equation) > 2)
    {
        throw InvalidRationalPrimitiveError(
            "parallel S-S point limiter equation exceeds quadratic degree");
    }
    return equation;
}

RationalPolynomial parallel_segment_limiter_equation(
    const RationalPrimitiveParameterization2& primitive,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const MatExactOpenSegmentSource2& limiter)
{
    RationalPolynomial first_clearance =
        parallel_support_clearance(
            primitive,
            first_segment);
    RationalPolynomial second_clearance =
        parallel_support_clearance(
            primitive,
            second_segment);
    normalize(first_clearance);
    normalize(second_clearance);
    if (first_clearance != second_clearance)
    {
        throw MismatchedLiveSegmentSegmentBridgeError(
            "parallel S-S segment limiter sources have unequal clearance");
    }
    RationalPolynomial equation =
        parallel_squared_support_distance(
            primitive,
            limiter);
    add_in_place(
        equation,
        scaled(first_clearance, -1));
    normalize(equation);
    if (polynomial_degree(equation) > 2)
    {
        throw InvalidRationalPrimitiveError(
            "parallel S-S segment limiter equation exceeds quadratic degree");
    }
    return equation;
}

template <typename CandidatePredicate>
MatParameterEndpoint2 bind_parallel_limiter_equation_endpoint(
    const RationalPrimitiveParameterization2& primitive,
    const RationalDomainRoot2& lower,
    const RationalDomainRoot2& upper,
    const RationalPolynomial& equation,
    const std::string& limiter_id,
    const std::string& provenance_kind,
    const MatTraits::Point_2& live_point,
    CandidatePredicate&& candidate_is_owned)
{
    if (is_zero(equation))
    {
        throw UnboundLiveSegmentSegmentEndpointError(
            "parallel S-S limiter equality is not zero-dimensional");
    }
    ExactAlgebraicKernel1 kernel;
    const Polynomial polynomial =
        integer_polynomial(equation);
    const std::vector<ExactFactorRootWitness2> witnesses =
        factor_root_witnesses(
            polynomial,
            kernel);
    const auto construct =
        kernel.construct_algebraic_real_1_object();
    const auto compare =
        kernel.compare_1_object();
    const Algebraic feature_lower =
        construct(lower.parameter);
    const Algebraic feature_upper =
        construct(upper.parameter);
    const FieldPolynomial2 x{
        primitive.x_coefficients,
        {0},
    };
    const FieldPolynomial2 y{
        primitive.y_coefficients,
        {0},
    };
    std::vector<MatParameterEndpoint2> matches;
    for (const ExactFactorRootWitness2& witness :
         witnesses)
    {
        if (compare(witness.root, feature_lower)
                == CGAL::SMALLER
            || compare(witness.root, feature_upper)
                == CGAL::LARGER)
        {
            continue;
        }
        const CORE::Expr parameter =
            exact_core_root(witness);
        const MatTraits::Point_2 candidate(
            core_evaluate(x, 1, parameter),
            core_evaluate(y, 1, parameter));
        if (candidate != live_point
            || !candidate_is_owned(candidate))
        {
            continue;
        }
        MatParameterEndpoint2 endpoint{
            witness.root,
            {
                witness.root_id,
                limiter_id,
                provenance_kind
                    + "/equation-factor-multiplicity/"
                    + std::to_string(
                        witness.source_factor_multiplicity),
            },
        };
        if (compare(witness.root, feature_lower)
            == CGAL::EQUAL)
        {
            union_stable_ids(
                endpoint.provenance_ids,
                lower.provenance_ids);
        }
        if (compare(witness.root, feature_upper)
            == CGAL::EQUAL)
        {
            union_stable_ids(
                endpoint.provenance_ids,
                upper.provenance_ids);
        }
        matches.push_back(std::move(endpoint));
    }
    if (matches.size() != 1)
    {
        throw UnboundLiveSegmentSegmentEndpointError(
            "parallel S-S live point binds "
            + std::to_string(matches.size())
            + " limiter roots");
    }
    return exact_graph_endpoint_binding(
        matches.front());
}

MatParameterEndpoint2 bind_parallel_point_limiter_endpoint(
    const RationalPrimitiveParameterization2& primitive,
    const RationalDomainRoot2& lower,
    const RationalDomainRoot2& upper,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const MatExactPointSiteSource2& limiter,
    const MatTraits::Site_2& live_owner,
    const MatTraits::Point_2& live_point)
{
    require_rational_parallel_point_limiter(
        limiter);
    if (!live_owner.is_point()
        || live_owner.point()
            != exact_live_point(
                limiter,
                limiter.radicand))
    {
        throw MismatchedLiveSegmentSegmentBridgeError(
            "parallel S-S point limiter record differs from live owner");
    }
    const RationalPolynomial equation =
        parallel_point_limiter_equation(
            primitive,
            first_segment,
            second_segment,
            limiter);
    return bind_parallel_limiter_equation_endpoint(
        primitive,
        lower,
        upper,
        equation,
        limiter.stable_site_id,
        "parallel-point-limiter",
        live_point,
        [](const MatTraits::Point_2&) {
            return true;
        });
}

MatParameterEndpoint2 bind_parallel_segment_limiter_endpoint(
    const RationalPrimitiveParameterization2& primitive,
    const RationalDomainRoot2& lower,
    const RationalDomainRoot2& upper,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const MatExactOpenSegmentSource2& limiter,
    const MatTraits::Site_2& live_owner,
    const MatTraits::Point_2& live_point)
{
    if (!live_owner.is_segment())
    {
        throw MismatchedLiveSegmentSegmentBridgeError(
            "parallel S-S segment limiter record differs from live owner");
    }
    const MatTraits::Point_2 limiter_source =
        live_owner.source();
    const MatTraits::Point_2 limiter_target =
        live_owner.target();
    if (limiter_source == limiter_target)
    {
        throw UnboundLiveSegmentSegmentEndpointError(
            "parallel S-S segment limiter has coincident endpoints");
    }
    const auto line_sign =
        [&limiter](const MatTraits::Point_2& point) {
            return CGAL::sign(
                CORE::Expr(limiter.line_a) * point.x()
                + CORE::Expr(limiter.line_b) * point.y()
                + CORE::Expr(limiter.line_c));
        };
    if (line_sign(limiter_source) != CGAL::ZERO
        || line_sign(limiter_target) != CGAL::ZERO)
    {
        throw MismatchedLiveSegmentSegmentBridgeError(
            "parallel S-S segment limiter support differs from live owner");
    }
    const RationalPolynomial equation =
        parallel_segment_limiter_equation(
            primitive,
            first_segment,
            second_segment,
            limiter);
    const CORE::Expr direction_x =
        limiter_target.x() - limiter_source.x();
    const CORE::Expr direction_y =
        limiter_target.y() - limiter_source.y();
    const CORE::Expr direction_norm =
        direction_x * direction_x
        + direction_y * direction_y;
    return bind_parallel_limiter_equation_endpoint(
        primitive,
        lower,
        upper,
        equation,
        limiter.stable_site_id,
        "parallel-segment-limiter",
        live_point,
        [&limiter_source,
         &direction_x,
         &direction_y,
         &direction_norm](
            const MatTraits::Point_2& candidate) {
            const CORE::Expr projection =
                (candidate.x() - limiter_source.x())
                    * direction_x
                + (candidate.y() - limiter_source.y())
                    * direction_y;
            return CGAL::sign(projection)
                    == CGAL::POSITIVE
                && CGAL::sign(
                       direction_norm - projection)
                    == CGAL::POSITIVE;
        });
}

bool same_unoriented_live_segment(
    const MatTraits::Segment_2& lhs,
    const MatTraits::Segment_2& rhs)
{
    return (lhs.source() == rhs.source()
            && lhs.target() == rhs.target())
        || (lhs.source() == rhs.target()
            && lhs.target() == rhs.source());
}

bool same_nonparallel_segment_chart(
    const NonparallelSegmentBisectorParameterization2& lhs,
    const NonparallelSegmentBisectorParameterization2& rhs)
{
    return lhs.first_segment_id == rhs.first_segment_id
        && lhs.second_segment_id == rhs.second_segment_id
        && lhs.branch_sign == rhs.branch_sign
        && lhs.x_rational == rhs.x_rational
        && lhs.x_radical == rhs.x_radical
        && lhs.y_rational == rhs.y_rational
        && lhs.y_radical == rhs.y_radical
        && lhs.radicand == rhs.radicand;
}

void require_rational_nonparallel_point_limiter(
    const MatExactPointSiteSource2& limiter)
{
    if (limiter.radicand <= 0
        || limiter.x.radical != 0
        || limiter.y.radical != 0)
    {
        throw NonRationalNonparallelSegmentPointLimiterError(
            "nonparallel S-S point limiter is not rational");
    }
}

FieldPolynomial2 nonparallel_point_limiter_equation(
    const NonparallelSegmentBisectorParameterization2& primitive,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const MatExactPointSiteSource2& limiter)
{
    require_nonzero_line_normal(
        first_segment,
        "first source");
    require_nonzero_line_normal(
        second_segment,
        "second source");
    require_rational_nonparallel_point_limiter(
        limiter);
    const FieldPolynomial2 x =
        field_coordinate(
            primitive.x_rational,
            primitive.x_radical);
    const FieldPolynomial2 y =
        field_coordinate(
            primitive.y_rational,
            primitive.y_radical);
    const CORE::BigRat first_norm =
        first_segment.line_a * first_segment.line_a
        + first_segment.line_b * first_segment.line_b;
    const CORE::BigRat second_norm =
        second_segment.line_a * second_segment.line_a
        + second_segment.line_b * second_segment.line_b;
    const FieldPolynomial2 first_squared =
        field_square(
            line_value(x, y, first_segment),
            primitive.radicand);
    const FieldPolynomial2 second_squared =
        field_square(
            line_value(x, y, second_segment),
            primitive.radicand);
    const FieldPolynomial2 first_equidistance =
        field_scale(
            first_squared,
            second_norm);
    const FieldPolynomial2 second_equidistance =
        field_scale(
            second_squared,
            first_norm);
    if (first_equidistance.rational
            != second_equidistance.rational
        || first_equidistance.radical
            != second_equidistance.radical)
    {
        throw MismatchedLiveNonparallelSegmentBridgeError(
            "nonparallel S-S chart is not equidistant from both sources");
    }
    const FieldPolynomial2 equation =
        field_subtract(
            field_scale(
                squared_distance(
                    x,
                    y,
                    limiter,
                    primitive.radicand),
                first_norm),
            first_squared);
    if (polynomial_degree(equation.rational) > 2
        || polynomial_degree(equation.radical) > 2)
    {
        throw InvalidRationalPrimitiveError(
            "nonparallel S-S point limiter equation exceeds quadratic degree");
    }
    return equation;
}

FieldPolynomial2 nonparallel_segment_limiter_equation(
    const NonparallelSegmentBisectorParameterization2& primitive,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const MatExactOpenSegmentSource2& limiter)
{
    require_nonzero_line_normal(
        first_segment,
        "first source");
    require_nonzero_line_normal(
        second_segment,
        "second source");
    require_nonzero_line_normal(
        limiter,
        "limiter");
    const FieldPolynomial2 x =
        field_coordinate(
            primitive.x_rational,
            primitive.x_radical);
    const FieldPolynomial2 y =
        field_coordinate(
            primitive.y_rational,
            primitive.y_radical);
    const CORE::BigRat first_norm =
        first_segment.line_a * first_segment.line_a
        + first_segment.line_b * first_segment.line_b;
    const CORE::BigRat second_norm =
        second_segment.line_a * second_segment.line_a
        + second_segment.line_b * second_segment.line_b;
    const CORE::BigRat limiter_norm =
        limiter.line_a * limiter.line_a
        + limiter.line_b * limiter.line_b;
    const FieldPolynomial2 first_squared =
        field_square(
            line_value(x, y, first_segment),
            primitive.radicand);
    const FieldPolynomial2 second_squared =
        field_square(
            line_value(x, y, second_segment),
            primitive.radicand);
    const FieldPolynomial2 first_equidistance =
        field_scale(
            first_squared,
            second_norm);
    const FieldPolynomial2 second_equidistance =
        field_scale(
            second_squared,
            first_norm);
    if (first_equidistance.rational
            != second_equidistance.rational
        || first_equidistance.radical
            != second_equidistance.radical)
    {
        throw MismatchedLiveNonparallelSegmentBridgeError(
            "nonparallel S-S chart is not equidistant from both sources");
    }
    const FieldPolynomial2 equation =
        field_subtract(
            field_scale(
                field_square(
                    line_value(x, y, limiter),
                    primitive.radicand),
                first_norm),
            field_scale(
                first_squared,
                limiter_norm));
    if (polynomial_degree(equation.rational) > 2
        || polynomial_degree(equation.radical) > 2)
    {
        throw InvalidRationalPrimitiveError(
            "nonparallel S-S segment limiter equation exceeds quadratic degree");
    }
    return equation;
}

MatTraits::Point_2 nonparallel_segment_chart_point(
    const NonparallelSegmentBisectorParameterization2& primitive,
    const CORE::Expr& parameter)
{
    return {
        core_evaluate(
            {
                primitive.x_rational,
                primitive.x_radical,
            },
            primitive.radicand,
            parameter),
        core_evaluate(
            {
                primitive.y_rational,
                primitive.y_radical,
            },
            primitive.radicand,
            parameter),
    };
}

MatParameterEndpoint2 bind_nonparallel_feature_endpoint(
    const NonparallelSegmentBisectorParameterization2& primitive,
    const NonparallelSegmentFeatureDomain2& feature_domain,
    const MatTraits::Point_2& live_point,
    const std::string& owner_id)
{
    const auto owned =
        [&owner_id](
            const MatQuadraticFieldDomainBoundary2& bound) {
            return std::find(
                       bound.provenance_ids.begin(),
                       bound.provenance_ids.end(),
                       owner_id)
                != bound.provenance_ids.end();
        };
    const auto point =
        [&primitive](
            const MatQuadraticFieldDomainBoundary2& bound) {
            return nonparallel_segment_chart_point(
                primitive,
                core_field_value(
                    bound.parameter,
                    primitive.radicand));
        };
    const bool is_lower =
        owned(feature_domain.lower)
        && live_point == point(feature_domain.lower);
    const bool is_upper =
        owned(feature_domain.upper)
        && live_point == point(feature_domain.upper);
    if (is_lower == is_upper)
    {
        throw MismatchedLiveNonparallelSegmentBridgeError(
            "nonparallel S-S owner differs from matched feature bound");
    }
    const MatQuadraticFieldDomainBoundary2& bound =
        is_lower
        ? feature_domain.lower
        : feature_domain.upper;
    MatParameterEndpoint2 endpoint{
        quadratic_field_algebraic_real(
            bound.parameter,
            feature_domain.radicand),
        bound.provenance_ids,
    };
    union_stable_ids(
        endpoint.provenance_ids,
        {owner_id});
    return exact_graph_endpoint_binding(
        endpoint);
}

MatParameterEndpoint2 bind_nonparallel_point_limiter_endpoint(
    const NonparallelSegmentBisectorParameterization2& primitive,
    const NonparallelSegmentFeatureDomain2& feature_domain,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const MatExactPointSiteSource2& limiter,
    const MatTraits::Site_2& live_owner,
    const MatTraits::Point_2& live_point)
{
    require_rational_nonparallel_point_limiter(
        limiter);
    if (!live_owner.is_point()
        || live_owner.point()
            != exact_live_point(
                limiter,
                limiter.radicand))
    {
        throw MismatchedLiveNonparallelSegmentBridgeError(
            "nonparallel S-S point limiter record differs from live owner");
    }
    const FieldPolynomial2 equation =
        nonparallel_point_limiter_equation(
            primitive,
            first_segment,
            second_segment,
            limiter);
    const FieldZeroSetPolynomial2 zero_set =
        field_zero_set_polynomial(
            equation,
            primitive.radicand);
    if (is_zero(zero_set.polynomial))
    {
        throw IdenticallyZeroNonparallelPointLimiterEquationError(
            "nonparallel S-S point limiter equality is not zero-dimensional");
    }
    ExactAlgebraicKernel1 kernel;
    const std::vector<ExactFactorRootWitness2>
        witnesses =
            factor_root_witnesses(
                integer_polynomial(
                    zero_set.polynomial),
                kernel);
    const Algebraic feature_lower =
        quadratic_field_algebraic_real(
            feature_domain.lower.parameter,
            feature_domain.radicand);
    const Algebraic feature_upper =
        quadratic_field_algebraic_real(
            feature_domain.upper.parameter,
            feature_domain.radicand);
    const auto compare =
        kernel.compare_1_object();
    std::vector<MatParameterEndpoint2> matches;
    for (const ExactFactorRootWitness2& witness :
         witnesses)
    {
        if (compare(witness.root, feature_lower)
                == CGAL::SMALLER
            || compare(witness.root, feature_upper)
                == CGAL::LARGER
            || field_sign_at(
                   equation,
                   primitive.radicand,
                   witness.root,
                   kernel)
                != CGAL::ZERO)
        {
            continue;
        }
        const CORE::Expr parameter =
            exact_core_root(witness);
        const MatTraits::Point_2 candidate =
            nonparallel_segment_chart_point(
                primitive,
                parameter);
        if (candidate != live_point)
        {
            continue;
        }
        MatParameterEndpoint2 endpoint{
            witness.root,
            {
                witness.root_id,
                limiter.stable_site_id,
                "nonparallel-point-limiter/"
                    + zero_set.factor_kind
                    + "/"
                    + std::to_string(
                        witness.source_factor_multiplicity),
            },
        };
        if (compare(witness.root, feature_lower)
            == CGAL::EQUAL)
        {
            union_stable_ids(
                endpoint.provenance_ids,
                feature_domain.lower.provenance_ids);
        }
        if (compare(witness.root, feature_upper)
            == CGAL::EQUAL)
        {
            union_stable_ids(
                endpoint.provenance_ids,
                feature_domain.upper.provenance_ids);
        }
        matches.push_back(std::move(endpoint));
    }
    if (matches.size() != 1)
    {
        throw UnboundLiveNonparallelSegmentEndpointError(
            "nonparallel S-S live point binds "
            + std::to_string(matches.size())
            + " point-limiter roots");
    }
    return exact_graph_endpoint_binding(
        matches.front());
}

MatParameterEndpoint2 bind_nonparallel_segment_limiter_endpoint(
    const NonparallelSegmentBisectorParameterization2& primitive,
    const NonparallelSegmentFeatureDomain2& feature_domain,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const MatExactOpenSegmentSource2& limiter,
    const MatTraits::Site_2& live_owner,
    const MatTraits::Point_2& live_point)
{
    if (!live_owner.is_segment())
    {
        throw MismatchedLiveNonparallelSegmentBridgeError(
            "nonparallel S-S segment limiter record differs from live owner");
    }
    const MatTraits::Point_2 limiter_source =
        live_owner.source();
    const MatTraits::Point_2 limiter_target =
        live_owner.target();
    if (limiter_source == limiter_target)
    {
        throw UnboundLiveNonparallelSegmentEndpointError(
            "nonparallel S-S segment limiter has coincident endpoints");
    }
    const auto line_sign =
        [&limiter](const MatTraits::Point_2& point) {
            return CGAL::sign(
                CORE::Expr(limiter.line_a) * point.x()
                + CORE::Expr(limiter.line_b) * point.y()
                + CORE::Expr(limiter.line_c));
        };
    if (line_sign(limiter_source) != CGAL::ZERO
        || line_sign(limiter_target) != CGAL::ZERO)
    {
        throw MismatchedLiveNonparallelSegmentBridgeError(
            "nonparallel S-S segment limiter support differs from live owner");
    }
    const FieldPolynomial2 equation =
        nonparallel_segment_limiter_equation(
            primitive,
            first_segment,
            second_segment,
            limiter);
    const FieldZeroSetPolynomial2 zero_set =
        field_zero_set_polynomial(
            equation,
            primitive.radicand);
    if (is_zero(zero_set.polynomial))
    {
        throw IdenticallyZeroNonparallelSegmentLimiterEquationError(
            "nonparallel S-S segment limiter equality is not zero-dimensional");
    }
    ExactAlgebraicKernel1 kernel;
    const std::vector<ExactFactorRootWitness2> witnesses =
        factor_root_witnesses(
            integer_polynomial(
                zero_set.polynomial),
            kernel);
    const Algebraic feature_lower =
        quadratic_field_algebraic_real(
            feature_domain.lower.parameter,
            feature_domain.radicand);
    const Algebraic feature_upper =
        quadratic_field_algebraic_real(
            feature_domain.upper.parameter,
            feature_domain.radicand);
    const auto compare =
        kernel.compare_1_object();
    const CORE::Expr direction_x =
        limiter_target.x() - limiter_source.x();
    const CORE::Expr direction_y =
        limiter_target.y() - limiter_source.y();
    const CORE::Expr direction_norm =
        direction_x * direction_x
        + direction_y * direction_y;
    std::vector<MatParameterEndpoint2> matches;
    for (const ExactFactorRootWitness2& witness :
         witnesses)
    {
        if (compare(witness.root, feature_lower)
                == CGAL::SMALLER
            || compare(witness.root, feature_upper)
                == CGAL::LARGER
            || field_sign_at(
                   equation,
                   primitive.radicand,
                   witness.root,
                   kernel)
                != CGAL::ZERO)
        {
            continue;
        }
        const CORE::Expr parameter =
            exact_core_root(witness);
        const MatTraits::Point_2 candidate =
            nonparallel_segment_chart_point(
                primitive,
                parameter);
        const CORE::Expr projection =
            (candidate.x() - limiter_source.x())
                * direction_x
            + (candidate.y() - limiter_source.y())
                * direction_y;
        if (candidate != live_point
            || CGAL::sign(projection)
                != CGAL::POSITIVE
            || CGAL::sign(
                   direction_norm - projection)
                != CGAL::POSITIVE)
        {
            continue;
        }
        MatParameterEndpoint2 endpoint{
            witness.root,
            {
                witness.root_id,
                limiter.stable_site_id,
                "nonparallel-segment-limiter/"
                    + zero_set.factor_kind
                    + "/"
                    + std::to_string(
                        witness.source_factor_multiplicity),
            },
        };
        if (compare(witness.root, feature_lower)
            == CGAL::EQUAL)
        {
            union_stable_ids(
                endpoint.provenance_ids,
                feature_domain.lower.provenance_ids);
        }
        if (compare(witness.root, feature_upper)
            == CGAL::EQUAL)
        {
            union_stable_ids(
                endpoint.provenance_ids,
                feature_domain.upper.provenance_ids);
        }
        matches.push_back(std::move(endpoint));
    }
    if (matches.size() != 1)
    {
        throw UnboundLiveNonparallelSegmentEndpointError(
            "nonparallel S-S live point binds "
            + std::to_string(matches.size())
            + " segment-limiter roots");
    }
    return exact_graph_endpoint_binding(
        matches.front());
}

template <typename StableSiteIdentity>
std::pair<MatParameterEndpoint2, MatParameterEndpoint2>
bind_parallel_segment_segment_cell_endpoints_impl(
    const RationalPrimitiveParameterization2& primitive,
    const RationalDomainRoot2& lower,
    const RationalDomainRoot2& upper,
    const MatExactOpenSegmentSource2* first_segment,
    const MatExactOpenSegmentSource2* second_segment,
    const std::vector<MatExactPointSiteSource2>& point_limiters,
    const std::vector<MatExactOpenSegmentSource2>& segment_limiters,
    const std::vector<std::string>& generator_ids,
    const StableSiteIdentity& stable_site_identity,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge)
{
    if (generator_ids.size() != 2
        || ordered_generator_site_ids(
               stable_site_identity(
                   halfedge->up()->site()),
               stable_site_identity(
                   halfedge->down()->site()))
            != generator_ids)
    {
        throw MismatchedLiveSegmentSegmentBridgeError(
            "parallel S-S source records do not match live generators");
    }
    if (!halfedge->has_source()
        || !halfedge->has_target())
    {
        throw UnboundLiveSegmentSegmentEndpointError(
            "parallel S-S halfedge is not bounded");
    }
    if (lower.parameter >= upper.parameter
        || primitive.domain_lower != lower.parameter
        || primitive.domain_upper != upper.parameter)
    {
        throw MismatchedLiveSegmentSegmentBridgeError(
            "parallel S-S feature domain is inconsistent");
    }
    if (lower.provenance_ids.empty()
        || upper.provenance_ids.empty())
    {
        throw UnboundLiveSegmentSegmentEndpointError(
            "parallel S-S feature endpoint provenance is empty");
    }

    const auto evaluate =
        [&primitive](const CORE::BigRat& parameter) {
            return MatTraits::Point_2(
                evaluate_rational_polynomial(
                    primitive.x_coefficients,
                    parameter),
                evaluate_rational_polynomial(
                    primitive.y_coefficients,
                    parameter));
        };
    const MatTraits::Point_2 expected_lower =
        evaluate(lower.parameter);
    const MatTraits::Point_2 expected_upper =
        evaluate(upper.parameter);
    const MatTraits::Point_2 live_source =
        halfedge->source()->point();
    const MatTraits::Point_2 live_target =
        halfedge->target()->point();
    const MatTraits::Site_2 source_owner =
        halfedge->left()->site();
    const MatTraits::Site_2 target_owner =
        halfedge->right()->site();
    const std::string source_owner_id =
        stable_site_identity(source_owner);
    const std::string target_owner_id =
        stable_site_identity(target_owner);

    MatTraits::Segment_2 primal;
    if (!CGAL::assign(
            primal,
            voronoi.dual().primal(
                halfedge->dual())))
    {
        throw MismatchedLiveSegmentSegmentBridgeError(
            "parallel S-S live dual is not a segment");
    }
    const bool primal_matches =
        (primal.source() == live_source
         && primal.target() == live_target)
        || (primal.source() == live_target
            && primal.target() == live_source);
    if (!primal_matches)
    {
        throw MismatchedLiveSegmentSegmentBridgeError(
            "parallel S-S primal and adaptor endpoints differ");
    }

    const auto feature_owner =
        [&lower, &upper](
            const std::string& owner_id) {
            return std::find(
                       lower.provenance_ids.begin(),
                       lower.provenance_ids.end(),
                       owner_id)
                    != lower.provenance_ids.end()
                || std::find(
                       upper.provenance_ids.begin(),
                       upper.provenance_ids.end(),
                       owner_id)
                    != upper.provenance_ids.end();
        };
    const auto bind_feature =
        [&lower,
         &upper,
         &expected_lower,
         &expected_upper](
            const MatTraits::Point_2& live_point,
            const std::string& owner_id) {
            const bool is_lower =
                live_point == expected_lower
                && std::find(
                       lower.provenance_ids.begin(),
                       lower.provenance_ids.end(),
                       owner_id)
                    != lower.provenance_ids.end();
            const bool is_upper =
                live_point == expected_upper
                && std::find(
                       upper.provenance_ids.begin(),
                       upper.provenance_ids.end(),
                       owner_id)
                    != upper.provenance_ids.end();
            if (is_lower == is_upper)
            {
                throw MismatchedLiveSegmentSegmentBridgeError(
                    "parallel S-S owner differs from matched feature bound");
            }
            const RationalDomainRoot2& bound =
                is_lower ? lower : upper;
            ExactAlgebraicKernel1 kernel;
            MatParameterEndpoint2 endpoint{
                kernel.construct_algebraic_real_1_object()(
                    bound.parameter),
                bound.provenance_ids,
            };
            union_stable_ids(
                endpoint.provenance_ids,
                {owner_id});
            return exact_graph_endpoint_binding(
                endpoint);
        };
    const auto bind =
        [&primitive,
         &lower,
         &upper,
         first_segment,
         second_segment,
         &point_limiters,
         &segment_limiters,
         &feature_owner,
         &bind_feature](
            const MatTraits::Point_2& live_point,
            const MatTraits::Site_2& live_owner,
            const std::string& owner_id) {
            if (feature_owner(owner_id))
            {
                return bind_feature(
                    live_point,
                    owner_id);
            }
            if (first_segment == nullptr
                || second_segment == nullptr)
            {
                throw UnsupportedSegmentSegmentLimiterError(
                    "parallel S-S endpoint has no exact external limiter binder");
            }
            if (live_owner.is_point())
            {
                std::vector<const MatExactPointSiteSource2*> matches;
                for (const MatExactPointSiteSource2& limiter :
                     point_limiters)
                {
                    if (limiter.stable_site_id == owner_id)
                    {
                        matches.push_back(&limiter);
                    }
                }
                if (matches.empty())
                {
                    throw UnsupportedSegmentSegmentLimiterError(
                        "parallel S-S point limiter has no source record");
                }
                if (matches.size() != 1)
                {
                    throw AmbiguousParallelSegmentPointLimiterError(
                        "parallel S-S point limiter identity is duplicated");
                }
                return bind_parallel_point_limiter_endpoint(
                    primitive,
                    lower,
                    upper,
                    *first_segment,
                    *second_segment,
                    *matches.front(),
                    live_owner,
                    live_point);
            }
            if (!live_owner.is_segment())
            {
                throw UnsupportedSegmentSegmentLimiterError(
                    "parallel S-S external limiter has no supported site kind");
            }
            std::vector<const MatExactOpenSegmentSource2*> matches;
            for (const MatExactOpenSegmentSource2& limiter :
                 segment_limiters)
            {
                if (limiter.stable_site_id == owner_id)
                {
                    matches.push_back(&limiter);
                }
            }
            if (matches.empty())
            {
                throw UnsupportedSegmentSegmentLimiterError(
                    "parallel S-S segment limiter has no source record");
            }
            if (matches.size() != 1)
            {
                throw AmbiguousParallelSegmentOpenLimiterError(
                    "parallel S-S segment limiter identity is duplicated");
            }
            return bind_parallel_segment_limiter_endpoint(
                primitive,
                lower,
                upper,
                *first_segment,
                *second_segment,
                *matches.front(),
                live_owner,
                live_point);
        };
    MatParameterEndpoint2 source =
        bind(
            live_source,
            source_owner,
            source_owner_id);
    MatParameterEndpoint2 target =
        bind(
            live_target,
            target_owner,
            target_owner_id);
    ExactAlgebraicKernel1 kernel;
    const CGAL::Comparison_result order =
        kernel.compare_1_object()(
            *source.parameter,
            *target.parameter);
    if (order == CGAL::EQUAL)
    {
        throw MismatchedLiveSegmentSegmentBridgeError(
            "parallel S-S live endpoints collapse");
    }
    if (order == CGAL::SMALLER)
    {
        return {
            std::move(source),
            std::move(target),
        };
    }
    return {
        std::move(target),
        std::move(source),
    };
}

} // namespace

std::pair<MatParameterEndpoint2, MatParameterEndpoint2>
bind_parallel_segment_segment_cell_endpoints(
    const RationalPrimitiveParameterization2& primitive,
    const RationalDomainRoot2& lower,
    const RationalDomainRoot2& upper,
    const std::vector<std::string>& generator_ids,
    const std::vector<GeneratorSite2>& generators,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge)
{
    const auto stable_site_identity =
        [&generators](
            const MatTraits::Site_2& site) {
            return stable_generator_site_id(
                site,
                generators);
        };
    return bind_parallel_segment_segment_cell_endpoints_impl(
        primitive,
        lower,
        upper,
        nullptr,
        nullptr,
        {},
        {},
        generator_ids,
        stable_site_identity,
        voronoi,
        halfedge);
}

namespace {

template <typename StableSiteIdentity>
std::pair<MatParameterEndpoint2, MatParameterEndpoint2>
bind_nonparallel_segment_segment_cell_endpoints_impl(
    const NonparallelSegmentBisectorParameterization2& primitive,
    const NonparallelSegmentFeatureDomain2& feature_domain,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const std::vector<MatExactPointSiteSource2>& point_limiters,
    const std::vector<MatExactOpenSegmentSource2>& segment_limiters,
    const std::vector<std::string>& generator_ids,
    const StableSiteIdentity& stable_site_identity,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge)
{
    if (generator_ids.size() != 2
        || ordered_generator_site_ids(
               stable_site_identity(
                   halfedge->up()->site()),
               stable_site_identity(
                   halfedge->down()->site()))
            != generator_ids
        || primitive.first_segment_id
            != generator_ids[0]
        || primitive.second_segment_id
            != generator_ids[1])
    {
        throw MismatchedLiveNonparallelSegmentBridgeError(
            "nonparallel S-S source records do not match live generators");
    }
    if (!halfedge->has_source()
        || !halfedge->has_target())
    {
        throw UnboundLiveNonparallelSegmentEndpointError(
            "nonparallel S-S halfedge is not bounded");
    }
    if (feature_domain.radicand
            != primitive.radicand
        || feature_domain.lower.provenance_ids.empty()
        || feature_domain.upper.provenance_ids.empty())
    {
        throw MismatchedLiveNonparallelSegmentBridgeError(
            "nonparallel S-S feature domain is inconsistent");
    }

    MatTraits::Segment_2 primal;
    if (!CGAL::assign(
            primal,
            voronoi.dual().primal(
                halfedge->dual())))
    {
        throw MismatchedLiveNonparallelSegmentBridgeError(
            "nonparallel S-S live dual is not a segment");
    }
    const MatTraits::Segment_2 adaptor_span(
        halfedge->source()->point(),
        halfedge->target()->point());
    if (!same_unoriented_live_segment(
            primal,
            adaptor_span))
    {
        throw MismatchedLiveNonparallelSegmentBridgeError(
            "nonparallel S-S primal and adaptor endpoints differ");
    }
    const NonparallelSegmentBisectorParameterization2
        canonical =
            nonparallel_segment_bisector_parameterization(
                first_segment,
                second_segment,
                primal);
    if (!same_nonparallel_segment_chart(
            primitive,
            canonical))
    {
        throw MismatchedLiveNonparallelSegmentBridgeError(
            "nonparallel S-S chart is not canonical for the live dual");
    }

    const MatTraits::Point_2 live_source =
        halfedge->source()->point();
    const MatTraits::Point_2 live_target =
        halfedge->target()->point();
    const MatTraits::Site_2 source_owner =
        halfedge->left()->site();
    const MatTraits::Site_2 target_owner =
        halfedge->right()->site();
    const std::string source_owner_id =
        stable_site_identity(source_owner);
    const std::string target_owner_id =
        stable_site_identity(target_owner);
    const auto feature_owner =
        [&feature_domain](
            const std::string& owner_id) {
            return std::find(
                       feature_domain.lower.provenance_ids.begin(),
                       feature_domain.lower.provenance_ids.end(),
                       owner_id)
                    != feature_domain.lower.provenance_ids.end()
                || std::find(
                       feature_domain.upper.provenance_ids.begin(),
                       feature_domain.upper.provenance_ids.end(),
                       owner_id)
                    != feature_domain.upper.provenance_ids.end();
        };
    const auto bind =
        [&primitive,
         &feature_domain,
         &first_segment,
         &second_segment,
         &point_limiters,
         &segment_limiters,
         &feature_owner](
            const MatTraits::Point_2& live_point,
            const MatTraits::Site_2& live_owner,
            const std::string& owner_id) {
            if (feature_owner(owner_id))
            {
                return bind_nonparallel_feature_endpoint(
                    primitive,
                    feature_domain,
                    live_point,
                    owner_id);
            }
            if (live_owner.is_point())
            {
                std::vector<
                    const MatExactPointSiteSource2*>
                    matches;
                for (const MatExactPointSiteSource2&
                         limiter : point_limiters)
                {
                    if (limiter.stable_site_id
                        == owner_id)
                    {
                        matches.push_back(&limiter);
                    }
                }
                if (matches.empty())
                {
                    throw UnsupportedNonparallelSegmentLimiterError(
                        "nonparallel S-S point limiter has no source record");
                }
                if (matches.size() != 1)
                {
                    throw AmbiguousNonparallelSegmentPointLimiterError(
                        "nonparallel S-S point limiter identity is duplicated");
                }
                return bind_nonparallel_point_limiter_endpoint(
                    primitive,
                    feature_domain,
                    first_segment,
                    second_segment,
                    *matches.front(),
                    live_owner,
                    live_point);
            }
            if (!live_owner.is_segment())
            {
                throw UnsupportedNonparallelSegmentLimiterError(
                    "nonparallel S-S external limiter is not an open segment");
            }
            std::vector<const MatExactOpenSegmentSource2*> matches;
            for (const MatExactOpenSegmentSource2& limiter :
                 segment_limiters)
            {
                if (limiter.stable_site_id == owner_id)
                {
                    matches.push_back(&limiter);
                }
            }
            if (matches.empty())
            {
                throw UnsupportedNonparallelSegmentLimiterError(
                    "nonparallel S-S segment limiter has no source record");
            }
            if (matches.size() != 1)
            {
                throw AmbiguousNonparallelSegmentOpenLimiterError(
                    "nonparallel S-S segment limiter identity is duplicated");
            }
            return bind_nonparallel_segment_limiter_endpoint(
                primitive,
                feature_domain,
                first_segment,
                second_segment,
                *matches.front(),
                live_owner,
                live_point);
        };
    MatParameterEndpoint2 source =
        bind(
            live_source,
            source_owner,
            source_owner_id);
    MatParameterEndpoint2 target =
        bind(
            live_target,
            target_owner,
            target_owner_id);
    ExactAlgebraicKernel1 kernel;
    const CGAL::Comparison_result order =
        kernel.compare_1_object()(
            *source.parameter,
            *target.parameter);
    if (order == CGAL::EQUAL)
    {
        throw MismatchedLiveNonparallelSegmentBridgeError(
            "nonparallel S-S live endpoints collapse");
    }
    if (order == CGAL::SMALLER)
    {
        return {
            std::move(source),
            std::move(target),
        };
    }
    return {
        std::move(target),
        std::move(source),
    };
}

} // namespace

std::pair<MatParameterEndpoint2, MatParameterEndpoint2>
bind_nonparallel_segment_segment_cell_endpoints(
    const NonparallelSegmentBisectorParameterization2& primitive,
    const NonparallelSegmentFeatureDomain2& feature_domain,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const std::vector<MatExactPointSiteSource2>& point_limiters,
    const std::vector<MatExactOpenSegmentSource2>& segment_limiters,
    const std::vector<std::string>& generator_ids,
    const std::vector<GeneratorSite2>& generators,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge)
{
    const auto stable_site_identity =
        [&generators](
            const MatTraits::Site_2& site) {
            return stable_generator_site_id(
                site,
                generators);
        };
    return bind_nonparallel_segment_segment_cell_endpoints_impl(
        primitive,
        feature_domain,
        first_segment,
        second_segment,
        point_limiters,
        segment_limiters,
        generator_ids,
        stable_site_identity,
        voronoi,
        halfedge);
}

std::pair<MatParameterEndpoint2, MatParameterEndpoint2>
bind_nonparallel_segment_segment_cell_endpoints(
    const NonparallelSegmentBisectorParameterization2& primitive,
    const NonparallelSegmentFeatureDomain2& feature_domain,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const std::vector<MatExactOpenSegmentSource2>& segment_limiters,
    const std::vector<std::string>& generator_ids,
    const std::vector<GeneratorSite2>& generators,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge)
{
    const auto stable_site_identity =
        [&generators](
            const MatTraits::Site_2& site) {
            return stable_generator_site_id(
                site,
                generators);
        };
    return bind_nonparallel_segment_segment_cell_endpoints_impl(
        primitive,
        feature_domain,
        first_segment,
        second_segment,
        {},
        segment_limiters,
        generator_ids,
        stable_site_identity,
        voronoi,
        halfedge);
}

std::pair<MatParameterEndpoint2, MatParameterEndpoint2>
bind_nonparallel_segment_segment_cell_endpoints(
    const NonparallelSegmentBisectorParameterization2& primitive,
    const NonparallelSegmentFeatureDomain2& feature_domain,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const std::vector<MatExactPointSiteSource2>& point_limiters,
    const std::vector<MatExactOpenSegmentSource2>& segment_limiters,
    const std::vector<std::string>& generator_ids,
    const CanonicalMatSiteGeometryIndex2& site_index,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge)
{
    const auto stable_site_identity =
        [&site_index](
            const MatTraits::Site_2& site)
            -> const std::string& {
            return site_index.stable_id(site);
        };
    return bind_nonparallel_segment_segment_cell_endpoints_impl(
        primitive,
        feature_domain,
        first_segment,
        second_segment,
        point_limiters,
        segment_limiters,
        generator_ids,
        stable_site_identity,
        voronoi,
        halfedge);
}

std::pair<MatParameterEndpoint2, MatParameterEndpoint2>
bind_nonparallel_segment_segment_cell_endpoints(
    const NonparallelSegmentBisectorParameterization2& primitive,
    const NonparallelSegmentFeatureDomain2& feature_domain,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const std::vector<MatExactOpenSegmentSource2>& segment_limiters,
    const std::vector<std::string>& generator_ids,
    const CanonicalMatSiteGeometryIndex2& site_index,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge)
{
    const auto stable_site_identity =
        [&site_index](
            const MatTraits::Site_2& site)
            -> const std::string& {
            return site_index.stable_id(site);
        };
    return bind_nonparallel_segment_segment_cell_endpoints_impl(
        primitive,
        feature_domain,
        first_segment,
        second_segment,
        {},
        segment_limiters,
        generator_ids,
        stable_site_identity,
        voronoi,
        halfedge);
}

std::pair<MatParameterEndpoint2, MatParameterEndpoint2>
bind_parallel_segment_segment_cell_endpoints(
    const RationalPrimitiveParameterization2& primitive,
    const RationalDomainRoot2& lower,
    const RationalDomainRoot2& upper,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const std::vector<MatExactPointSiteSource2>& point_limiters,
    const std::vector<std::string>& generator_ids,
    const std::vector<GeneratorSite2>& generators,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge)
{
    require_canonical_parallel_segment_chart(
        primitive,
        first_segment,
        second_segment,
        generator_ids);
    const auto stable_site_identity =
        [&generators](
            const MatTraits::Site_2& site) {
            return stable_generator_site_id(
                site,
                generators);
        };
    return bind_parallel_segment_segment_cell_endpoints_impl(
        primitive,
        lower,
        upper,
        &first_segment,
        &second_segment,
        point_limiters,
        {},
        generator_ids,
        stable_site_identity,
        voronoi,
        halfedge);
}

std::pair<MatParameterEndpoint2, MatParameterEndpoint2>
bind_parallel_segment_segment_cell_endpoints(
    const RationalPrimitiveParameterization2& primitive,
    const RationalDomainRoot2& lower,
    const RationalDomainRoot2& upper,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const std::vector<MatExactPointSiteSource2>& point_limiters,
    const std::vector<MatExactOpenSegmentSource2>& segment_limiters,
    const std::vector<std::string>& generator_ids,
    const std::vector<GeneratorSite2>& generators,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge)
{
    require_canonical_parallel_segment_chart(
        primitive,
        first_segment,
        second_segment,
        generator_ids);
    const auto stable_site_identity =
        [&generators](
            const MatTraits::Site_2& site) {
            return stable_generator_site_id(
                site,
                generators);
        };
    return bind_parallel_segment_segment_cell_endpoints_impl(
        primitive,
        lower,
        upper,
        &first_segment,
        &second_segment,
        point_limiters,
        segment_limiters,
        generator_ids,
        stable_site_identity,
        voronoi,
        halfedge);
}

std::pair<MatParameterEndpoint2, MatParameterEndpoint2>
bind_parallel_segment_segment_cell_endpoints(
    const RationalPrimitiveParameterization2& primitive,
    const RationalDomainRoot2& lower,
    const RationalDomainRoot2& upper,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const std::vector<MatExactPointSiteSource2>& point_limiters,
    const std::vector<MatExactOpenSegmentSource2>& segment_limiters,
    const std::vector<std::string>& generator_ids,
    const CanonicalMatSiteGeometryIndex2& site_index,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge)
{
    require_canonical_parallel_segment_chart(
        primitive,
        first_segment,
        second_segment,
        generator_ids);
    const auto stable_site_identity =
        [&site_index](
            const MatTraits::Site_2& site)
            -> const std::string& {
            return site_index.stable_id(site);
        };
    return bind_parallel_segment_segment_cell_endpoints_impl(
        primitive,
        lower,
        upper,
        &first_segment,
        &second_segment,
        point_limiters,
        segment_limiters,
        generator_ids,
        stable_site_identity,
        voronoi,
        halfedge);
}
