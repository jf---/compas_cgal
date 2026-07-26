#include "segment_site_endpoint_binding.h"

#include <algorithm>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

#include <CGAL/CORE/poly/Poly.h>
#include <CGAL/Polynomial_traits_d.h>

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
    if (limiter_line_norm == 0)
    {
        throw InvalidRationalPrimitiveError(
            "segment limiter has a zero supporting-line normal");
    }
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
    const CORE::BigRat radicand = focus.radicand;
    require_same_field(segment_source, radicand);
    require_same_field(segment_target, radicand);
    require_same_field(limiter, radicand);
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
    const CORE::BigRat radicand = focus.radicand;
    require_same_field(segment_source, radicand);
    require_same_field(segment_target, radicand);
    require_same_field(limiter_source, radicand);
    require_same_field(limiter_target, radicand);
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
