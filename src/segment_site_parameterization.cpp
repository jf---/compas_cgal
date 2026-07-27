#include "segment_site_parameterization.h"
#include "segment_site_provenance.h"

#include <algorithm>
#include <optional>
#include <utility>

#include <CGAL/Fraction_traits.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/number_utils.h>

namespace
{

std::optional<CORE::BigRat> rational_square_root(
    const CORE::BigRat& value)
{
    if (value < 0) {
        return std::nullopt;
    }
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
    if (!CGAL::is_square(
            numerator,
            numerator_root)
        || !CGAL::is_square(
            denominator,
            denominator_root)) {
        return std::nullopt;
    }
    return typename FractionTraits::Compose()(
        numerator_root,
        denominator_root);
}

}  // namespace

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

MatExactOpenSegmentSource2 canonical_open_segment_source(
    std::string stable_site_id,
    const CORE::BigRat& source_x,
    const CORE::BigRat& source_y,
    const CORE::BigRat& target_x,
    const CORE::BigRat& target_y)
{
    if (stable_site_id.empty()) {
        throw EmptyOpenSegmentSourceIdentityError(
            "open segment source identity is empty");
    }

    const CORE::BigRat line_a = source_y - target_y;
    const CORE::BigRat line_b = target_x - source_x;
    if (line_a == 0 && line_b == 0) {
        throw DegenerateOpenSegmentSourceError(
            "open segment source endpoints coincide");
    }

    std::vector<ExactAlgebraicInteger1> line =
        primitive_integer_coefficients({
            line_a,
            line_b,
            source_x * target_y - target_x * source_y,
        });
    const auto first_nonzero = std::find_if(
        line.begin(),
        line.end(),
        [](const ExactAlgebraicInteger1& coefficient) {
            return coefficient != 0;
        });
    if (*first_nonzero < 0) {
        for (ExactAlgebraicInteger1& coefficient : line) {
            coefficient = -coefficient;
        }
    }

    return {
        std::move(stable_site_id),
        CORE::BigRat(line[0]),
        CORE::BigRat(line[1]),
        CORE::BigRat(line[2]),
    };
}

RationalPrimitiveParameterization2
parallel_segment_bisector_parameterization(
    const MatExactOpenSegmentSource2& first,
    const MatExactOpenSegmentSource2& second)
{
    if (first.stable_site_id == second.stable_site_id) {
        throw DuplicateOpenSegmentSourceIdentityError(
            "parallel segment sources share one stable identity");
    }
    const MatExactOpenSegmentSource2* ordered_first =
        &first;
    const MatExactOpenSegmentSource2* ordered_second =
        &second;
    if (ordered_second->stable_site_id
        < ordered_first->stable_site_id) {
        std::swap(ordered_first, ordered_second);
    }

    const CORE::BigRat determinant =
        ordered_first->line_a * ordered_second->line_b
        - ordered_second->line_a * ordered_first->line_b;
    if (determinant != 0) {
        throw NonparallelSegmentSupportsError(
            "parallel S-S chart received nonparallel supports");
    }
    const CORE::BigRat normal_scale =
        ordered_first->line_a != 0
        ? ordered_second->line_a / ordered_first->line_a
        : ordered_second->line_b / ordered_first->line_b;
    if (normal_scale <= 0) {
        throw InvalidRationalPrimitiveError(
            "canonical parallel support normals disagree");
    }
    const CORE::BigRat second_constant =
        ordered_second->line_c / normal_scale;
    if (second_constant == ordered_first->line_c) {
        throw CoincidentSegmentSupportsError(
            "parallel segment supports coincide");
    }

    const CORE::BigRat middle_constant =
        (ordered_first->line_c + second_constant) / 2;
    const CORE::BigRat normal_squared =
        ordered_first->line_a * ordered_first->line_a
        + ordered_first->line_b * ordered_first->line_b;
    CORE::BigRat direction_x = ordered_first->line_b;
    CORE::BigRat direction_y = -ordered_first->line_a;
    if (direction_x < 0
        || (direction_x == 0 && direction_y < 0)) {
        direction_x = -direction_x;
        direction_y = -direction_y;
    }
    return {
        {
            -ordered_first->line_a
                * middle_constant / normal_squared,
            direction_x,
        },
        {
            -ordered_first->line_b
                * middle_constant / normal_squared,
            direction_y,
        },
        std::nullopt,
        std::nullopt,
    };
}

CORE::BigRat parallel_segment_tangent_parameter(
    const RationalPrimitiveParameterization2& primitive,
    const CORE::BigRat& x,
    const CORE::BigRat& y)
{
    if (primitive.x_coefficients.size() != 2
        || primitive.y_coefficients.size() != 2) {
        throw InvalidRationalPrimitiveError(
            "parallel S-S chart must be affine linear");
    }
    const CORE::BigRat direction_x =
        primitive.x_coefficients[1];
    const CORE::BigRat direction_y =
        primitive.y_coefficients[1];
    const CORE::BigRat direction_squared =
        direction_x * direction_x
        + direction_y * direction_y;
    if (direction_squared == 0) {
        throw InvalidRationalPrimitiveError(
            "parallel S-S chart direction is zero");
    }
    return (
        (x - primitive.x_coefficients[0]) * direction_x
        + (y - primitive.y_coefficients[0]) * direction_y)
        / direction_squared;
}

NonparallelSegmentBisectorParameterization2
nonparallel_segment_bisector_parameterization(
    const MatExactOpenSegmentSource2& first,
    const MatExactOpenSegmentSource2& second,
    const MatTraits::Segment_2& live_primitive)
{
    if (first.stable_site_id == second.stable_site_id) {
        throw DuplicateOpenSegmentSourceIdentityError(
            "nonparallel segment sources share one stable identity");
    }
    const MatExactOpenSegmentSource2* ordered_first =
        &first;
    const MatExactOpenSegmentSource2* ordered_second =
        &second;
    if (ordered_second->stable_site_id
        < ordered_first->stable_site_id) {
        std::swap(ordered_first, ordered_second);
    }

    const CORE::BigRat determinant =
        ordered_first->line_a * ordered_second->line_b
        - ordered_second->line_a * ordered_first->line_b;
    if (determinant == 0) {
        throw ParallelSegmentSupportsError(
            "nonparallel S-S chart received parallel supports");
    }
    if (live_primitive.source()
        == live_primitive.target()) {
        throw DegenerateLiveSegmentPrimitiveError(
            "nonparallel S-S live primitive has zero length");
    }

    const CORE::BigRat first_norm =
        ordered_first->line_a * ordered_first->line_a
        + ordered_first->line_b * ordered_first->line_b;
    const CORE::BigRat second_norm =
        ordered_second->line_a * ordered_second->line_a
        + ordered_second->line_b * ordered_second->line_b;
    const CORE::BigRat radicand =
        first_norm * second_norm;
    const CORE::Expr radical =
        CORE::sqrt(CORE::Expr(radicand));

    int branch_sign = 0;
    for (const int candidate_sign : {1, -1}) {
        const CORE::Expr candidate_a =
            CORE::Expr(
                second_norm * ordered_first->line_a)
            + CORE::Expr(
                  candidate_sign
                  * ordered_second->line_a)
                * radical;
        const CORE::Expr candidate_b =
            CORE::Expr(
                second_norm * ordered_first->line_b)
            + CORE::Expr(
                  candidate_sign
                  * ordered_second->line_b)
                * radical;
        const CORE::Expr candidate_c =
            CORE::Expr(
                second_norm * ordered_first->line_c)
            + CORE::Expr(
                  candidate_sign
                  * ordered_second->line_c)
                * radical;
        const MatTraits::Line_2 candidate(
            candidate_a,
            candidate_b,
            candidate_c);
        if (candidate.has_on(live_primitive.source())
            && candidate.has_on(
                live_primitive.target())) {
            branch_sign = candidate_sign;
            break;
        }
    }
    if (branch_sign == 0) {
        throw UnboundNonparallelSegmentBranchError(
            "live S-S primitive lies on neither exact bisector branch");
    }

    const CORE::BigRat intersection_x =
        (
            ordered_first->line_b
                * ordered_second->line_c
            - ordered_second->line_b
                * ordered_first->line_c)
        / determinant;
    const CORE::BigRat intersection_y =
        (
            ordered_first->line_c
                * ordered_second->line_a
            - ordered_second->line_c
                * ordered_first->line_a)
        / determinant;
    CORE::BigRat direction_x_rational =
        second_norm * ordered_first->line_b;
    CORE::BigRat direction_x_radical =
        branch_sign * ordered_second->line_b;
    CORE::BigRat direction_y_rational =
        -second_norm * ordered_first->line_a;
    CORE::BigRat direction_y_radical =
        -branch_sign * ordered_second->line_a;
    const CGAL::Sign direction_x_sign =
        quadratic_field_sign(
            {
                direction_x_rational,
                direction_x_radical,
            },
            radicand);
    const CGAL::Sign leading_sign =
        direction_x_sign != CGAL::ZERO
        ? direction_x_sign
        : quadratic_field_sign(
              {
                  direction_y_rational,
                  direction_y_radical,
              },
              radicand);
    if (leading_sign == CGAL::ZERO) {
        throw InvalidRationalPrimitiveError(
            "nonparallel S-S chart direction is zero");
    }
    if (leading_sign == CGAL::NEGATIVE) {
        direction_x_rational =
            -direction_x_rational;
        direction_x_radical =
            -direction_x_radical;
        direction_y_rational =
            -direction_y_rational;
        direction_y_radical =
            -direction_y_radical;
    }

    CORE::BigRat canonical_radicand = radicand;
    if (const std::optional<CORE::BigRat> square_root =
            rational_square_root(radicand);
        square_root.has_value()) {
        direction_x_rational +=
            direction_x_radical * *square_root;
        direction_y_rational +=
            direction_y_radical * *square_root;
        direction_x_radical = 0;
        direction_y_radical = 0;
        canonical_radicand = 1;
    }

    return {
        ordered_first->stable_site_id,
        ordered_second->stable_site_id,
        branch_sign,
        {intersection_x, direction_x_rational},
        {0, direction_x_radical},
        {intersection_y, direction_y_rational},
        {0, direction_y_radical},
        canonical_radicand,
    };
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

CORE::BigRat evaluate_rational_polynomial(
    const std::vector<CORE::BigRat>& coefficients,
    const CORE::BigRat& parameter)
{
    CORE::BigRat value = 0;
    for (auto coefficient = coefficients.rbegin();
         coefficient != coefficients.rend();
         ++coefficient) {
        value *= parameter;
        value += *coefficient;
    }
    return value;
}

bool domain_contains(
    const MatDomainPolygonWithHoles2& domain,
    const CORE::BigRat& x,
    const CORE::BigRat& y)
{
    const MatDomainKernel2::Point_2 point(x, y);
    if (domain.outer_boundary().bounded_side(point)
        == CGAL::ON_UNBOUNDED_SIDE) {
        return false;
    }
    for (auto hole = domain.holes_begin();
         hole != domain.holes_end();
         ++hole) {
        if (hole->bounded_side(point)
            == CGAL::ON_BOUNDED_SIDE) {
            return false;
        }
    }
    return true;
}

ExactAlgebraicKernel2::Polynomial_2
parabola_edge_equation(
    const std::vector<CORE::BigRat>& coordinate,
    const CORE::BigRat& edge_origin,
    const CORE::BigRat& edge_direction)
{
    const CORE::BigRat constant =
        coordinate.empty()
        ? CORE::BigRat(0) - edge_origin
        : coordinate[0] - edge_origin;
    std::vector<CORE::BigRat> coefficients{
        constant,
        coordinate.size() > 1
            ? coordinate[1]
            : CORE::BigRat(0),
        coordinate.size() > 2
            ? coordinate[2]
            : CORE::BigRat(0),
        -edge_direction,
    };
    const std::vector<ExactAlgebraicInteger1> integer =
        primitive_integer_coefficients(coefficients);
    using Polynomial = ExactAlgebraicKernel2::Polynomial_2;
    const Polynomial parameter =
        CGAL::shift(Polynomial(ExactAlgebraicInteger1(1)), 1, 0);
    const Polynomial edge_parameter =
        CGAL::shift(Polynomial(ExactAlgebraicInteger1(1)), 1, 1);
    return Polynomial(integer[0])
        + integer[1] * parameter
        + integer[2] * parameter * parameter
        + integer[3] * edge_parameter;
}

void append_parabola_polygon_intersections(
    const MatDomainPolygon2& polygon,
    const std::string& ring_id,
    const RationalPrimitiveParameterization2& primitive,
    std::vector<AlgebraicDomainRoot2>& algebraic_roots)
{
    ExactAlgebraicKernel1 kernel1;
    ExactAlgebraicKernel2 kernel2;
    std::size_t edge_index = 0;
    for (auto edge = polygon.edges_begin();
         edge != polygon.edges_end();
         ++edge, ++edge_index) {
        const CORE::BigRat ax = edge->source().x();
        const CORE::BigRat ay = edge->source().y();
        const CORE::BigRat ex = edge->target().x() - ax;
        const CORE::BigRat ey = edge->target().y() - ay;
        const auto x_equation = parabola_edge_equation(
            primitive.x_coefficients,
            ax,
            ex);
        const auto y_equation = parabola_edge_equation(
            primitive.y_coefficients,
            ay,
            ey);
        if (CGAL::is_zero(x_equation)
            || CGAL::is_zero(y_equation)
            || !kernel2.is_coprime_2_object()(
                x_equation,
                y_equation)) {
            throw OverlappingDomainBoundaryError(
                "parabola has a positive-dimensional boundary intersection");
        }

        std::vector<
            std::pair<
                ExactAlgebraicKernel2::Algebraic_real_2,
                ExactAlgebraicKernel2::Multiplicity_type>>
            solutions;
        kernel2.solve_2_object()(
            x_equation,
            y_equation,
            std::back_inserter(solutions));
        for (const auto& [solution, multiplicity] : solutions) {
            static_cast<void>(multiplicity);
            const auto parameter =
                kernel2.compute_x_2_object()(solution);
            const auto edge_parameter =
                kernel2.compute_y_2_object()(solution);
            const auto compare = kernel1.compare_1_object();
            if (compare(edge_parameter, CORE::BigRat(0))
                    == CGAL::SMALLER
                || compare(edge_parameter, CORE::BigRat(1))
                    == CGAL::LARGER
                || (primitive.domain_lower.has_value()
                    && compare(
                           parameter,
                           *primitive.domain_lower)
                        == CGAL::SMALLER)
                || (primitive.domain_upper.has_value()
                    && compare(
                           parameter,
                           *primitive.domain_upper)
                        == CGAL::LARGER)) {
                continue;
            }
            algebraic_roots.push_back(
                {
                    parameter,
                    {
                        ring_id + "/edge-"
                        + std::to_string(edge_index),
                        algebraic_root_identity_v1(parameter),
                    },
                });
        }
    }
}

std::vector<AlgebraicDomainRoot2> parabola_domain_roots(
    const std::string& dual_id,
    const RationalPrimitiveParameterization2& primitive,
    const MatDomainPolygonWithHoles2& domain)
{
    if (primitive.x_coefficients.size() > 3
        || primitive.y_coefficients.size() > 3
        || (primitive.x_coefficients.size() < 3
            && primitive.y_coefficients.size() < 3)) {
        throw InvalidRationalPrimitiveError(
            "parabola D clipping requires a quadratic parameterization");
    }
    if (domain.is_unbounded()) {
        throw InvalidRationalPrimitiveError(
            "D clipping requires a bounded outer polygon");
    }

    std::vector<AlgebraicDomainRoot2> roots;
    append_parabola_polygon_intersections(
        domain.outer_boundary(),
        dual_id + "/D-outer",
        primitive,
        roots);
    std::size_t hole_index = 0;
    for (auto hole = domain.holes_begin();
         hole != domain.holes_end();
         ++hole, ++hole_index) {
        append_parabola_polygon_intersections(
            *hole,
            dual_id + "/D-hole-"
                + std::to_string(hole_index),
            primitive,
            roots);
    }
    ExactAlgebraicKernel1 kernel;
    const auto compare = kernel.compare_1_object();
    std::sort(
        roots.begin(),
        roots.end(),
        [&compare](const AlgebraicDomainRoot2& lhs,
                   const AlgebraicDomainRoot2& rhs) {
            return compare(lhs.parameter, rhs.parameter)
                == CGAL::SMALLER;
        });
    std::vector<AlgebraicDomainRoot2> unique;
    for (AlgebraicDomainRoot2& root : roots) {
        if (!unique.empty()
            && compare(
                   unique.back().parameter,
                   root.parameter)
                == CGAL::EQUAL) {
            unique.back().provenance_ids.insert(
                unique.back().provenance_ids.end(),
                root.provenance_ids.begin(),
                root.provenance_ids.end());
        } else {
            unique.push_back(std::move(root));
        }
    }
    return unique;
}

MatQuadraticFieldValue2 field_add(
    const MatQuadraticFieldValue2& lhs,
    const MatQuadraticFieldValue2& rhs)
{
    return {
        lhs.rational + rhs.rational,
        lhs.radical + rhs.radical,
    };
}

MatQuadraticFieldValue2 field_scale(
    const MatQuadraticFieldValue2& value,
    const CORE::BigRat& scale)
{
    return {
        value.rational * scale,
        value.radical * scale,
    };
}

MatQuadraticFieldValue2 field_multiply(
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

MatQuadraticFieldValue2 field_divide(
    const MatQuadraticFieldValue2& numerator,
    const MatQuadraticFieldValue2& denominator,
    const CORE::BigRat& radicand)
{
    const CORE::BigRat norm =
        denominator.rational * denominator.rational
        - radicand
            * denominator.radical
            * denominator.radical;
    if (norm == 0) {
        throw InvalidRationalPrimitiveError(
            "source quadratic-field denominator is zero");
    }
    return field_scale(
        field_multiply(
            numerator,
            {
                denominator.rational,
                -denominator.radical,
            },
            radicand),
        CORE::BigRat(1) / norm);
}

CGAL::Sign quadratic_field_sign(
    const MatQuadraticFieldValue2& value,
    const CORE::BigRat& radicand)
{
    if (value.radical == 0) {
        return CGAL::sign(value.rational);
    }
    if (value.rational == 0) {
        return CGAL::sign(value.radical);
    }
    const CGAL::Sign rational_sign =
        CGAL::sign(value.rational);
    const CGAL::Sign radical_sign =
        CGAL::sign(value.radical);
    if (rational_sign == radical_sign) {
        return rational_sign;
    }
    const CORE::BigRat comparison =
        value.rational * value.rational
        - radicand * value.radical * value.radical;
    const CGAL::Sign magnitude = CGAL::sign(comparison);
    if (magnitude == CGAL::ZERO) {
        return CGAL::ZERO;
    }
    return magnitude == CGAL::POSITIVE
        ? rational_sign
        : radical_sign;
}

SourceParabolaParameterization2 source_parameterization(
    const MatExactPointSiteSource2& point,
    const MatExactOpenSegmentSource2& segment)
{
    if (point.stable_site_id.empty()
        || segment.stable_site_id.empty()) {
        throw InvalidRationalPrimitiveError(
            "source site identity is empty");
    }
    if (point.radicand <= 0) {
        throw InvalidRationalPrimitiveError(
            "source quadratic radicand is not positive");
    }
    const CORE::BigRat line_norm =
        segment.line_a * segment.line_a
        + segment.line_b * segment.line_b;
    if (line_norm == 0) {
        throw InvalidRationalPrimitiveError(
            "source directrix is degenerate");
    }
    const MatQuadraticFieldValue2 signed_distance_numerator =
        field_add(
            field_add(
                field_scale(point.x, segment.line_a),
                field_scale(point.y, segment.line_b)),
            {segment.line_c, 0});
    if (quadratic_field_sign(
            signed_distance_numerator,
            point.radicand)
        == CGAL::ZERO) {
        throw InvalidRationalPrimitiveError(
            "source focus lies on directrix");
    }

    const MatQuadraticFieldValue2 vertex_x =
        field_add(
            point.x,
            field_scale(
                signed_distance_numerator,
                -segment.line_a
                    / (CORE::BigRat(2) * line_norm)));
    const MatQuadraticFieldValue2 vertex_y =
        field_add(
            point.y,
            field_scale(
                signed_distance_numerator,
                -segment.line_b
                    / (CORE::BigRat(2) * line_norm)));
    const MatQuadraticFieldValue2 quadratic_scale =
        field_divide(
            {line_norm, 0},
            field_scale(
                signed_distance_numerator,
                CORE::BigRat(2)),
            point.radicand);
    const MatQuadraticFieldValue2 quadratic_x =
        field_scale(quadratic_scale, segment.line_a);
    const MatQuadraticFieldValue2 quadratic_y =
        field_scale(quadratic_scale, segment.line_b);
    return {
        {
            vertex_x.rational,
            -segment.line_b,
            quadratic_x.rational,
        },
        {
            vertex_x.radical,
            0,
            quadratic_x.radical,
        },
        {
            vertex_y.rational,
            segment.line_a,
            quadratic_y.rational,
        },
        {
            vertex_y.radical,
            0,
            quadratic_y.radical,
        },
        point.radicand,
    };
}

RadicalEquation2 radical_equation(
    const std::vector<CORE::BigRat>& rational,
    const std::vector<CORE::BigRat>& radical,
    const CORE::BigRat& edge_origin,
    const CORE::BigRat& edge_direction)
{
    std::vector<CORE::BigRat> values{
        rational[0] - edge_origin,
        rational[1],
        rational[2],
        -edge_direction,
        radical[0],
        radical[1],
        radical[2],
    };
    const std::vector<ExactAlgebraicInteger1> integer =
        primitive_integer_coefficients(values);
    using Polynomial = ExactAlgebraicKernel2::Polynomial_2;
    const Polynomial parameter =
        CGAL::shift(Polynomial(ExactAlgebraicInteger1(1)), 1, 0);
    const Polynomial edge_parameter =
        CGAL::shift(Polynomial(ExactAlgebraicInteger1(1)), 1, 1);
    return {
        Polynomial(integer[0])
            + integer[1] * parameter
            + integer[2] * parameter * parameter
            + integer[3] * edge_parameter,
        Polynomial(integer[4])
            + integer[5] * parameter
            + integer[6] * parameter * parameter,
    };
}

ExactAlgebraicKernel2::Polynomial_2 radical_norm(
    const RadicalEquation2& equation,
    const CORE::BigRat& radicand)
{
    const std::vector<ExactAlgebraicInteger1> integer =
        primitive_integer_coefficients(
            {CORE::BigRat(1), -radicand});
    return integer[0]
            * equation.rational
            * equation.rational
        + integer[1]
            * equation.radical
            * equation.radical;
}

bool radical_equation_holds(
    const RadicalEquation2& equation,
    const ExactAlgebraicKernel2::Algebraic_real_2& point,
    const ExactAlgebraicKernel2& kernel)
{
    const CGAL::Sign rational =
        kernel.sign_at_2_object()(
            equation.rational,
            point);
    const CGAL::Sign radical =
        kernel.sign_at_2_object()(
            equation.radical,
            point);
    return (rational == CGAL::ZERO
            && radical == CGAL::ZERO)
        || rational == -radical;
}

MatQuadraticFieldValue2 evaluate_source_coordinate(
    const std::vector<CORE::BigRat>& rational,
    const std::vector<CORE::BigRat>& radical,
    const CORE::BigRat& parameter)
{
    return {
        evaluate_rational_polynomial(rational, parameter),
        evaluate_rational_polynomial(radical, parameter),
    };
}

bool quadratic_point_in_polygon(
    const MatDomainPolygon2& polygon,
    const MatQuadraticFieldValue2& x,
    const MatQuadraticFieldValue2& y,
    const CORE::BigRat& radicand)
{
    int winding = 0;
    for (auto edge = polygon.edges_begin();
         edge != polygon.edges_end();
         ++edge) {
        const MatQuadraticFieldValue2 source_y{
            edge->source().y() - y.rational,
            -y.radical,
        };
        const MatQuadraticFieldValue2 target_y{
            edge->target().y() - y.rational,
            -y.radical,
        };
        const MatQuadraticFieldValue2 orientation =
            field_add(
                field_scale(
                    {
                        edge->target().y() - y.rational,
                        -y.radical,
                    },
                    edge->source().x() - x.rational),
                field_scale(
                    {
                        edge->source().y() - y.rational,
                        -y.radical,
                    },
                    x.rational - edge->target().x()));
        const CGAL::Sign source_side =
            quadratic_field_sign(source_y, radicand);
        const CGAL::Sign target_side =
            quadratic_field_sign(target_y, radicand);
        const CGAL::Sign turn =
            quadratic_field_sign(orientation, radicand);
        if (source_side != CGAL::POSITIVE
            && target_side == CGAL::POSITIVE
            && turn == CGAL::POSITIVE) {
            ++winding;
        } else if (
            source_side == CGAL::POSITIVE
            && target_side != CGAL::POSITIVE
            && turn == CGAL::NEGATIVE) {
            --winding;
        }
    }
    return winding != 0;
}

bool source_domain_contains(
    const MatDomainPolygonWithHoles2& domain,
    const SourceParabolaParameterization2& primitive,
    const CORE::BigRat& parameter)
{
    const MatQuadraticFieldValue2 x =
        evaluate_source_coordinate(
            primitive.x_rational,
            primitive.x_radical,
            parameter);
    const MatQuadraticFieldValue2 y =
        evaluate_source_coordinate(
            primitive.y_rational,
            primitive.y_radical,
            parameter);
    if (!quadratic_point_in_polygon(
            domain.outer_boundary(),
            x,
            y,
            primitive.radicand)) {
        return false;
    }
    for (auto hole = domain.holes_begin();
         hole != domain.holes_end();
         ++hole) {
        if (quadratic_point_in_polygon(
                *hole,
                x,
                y,
                primitive.radicand)) {
            return false;
        }
    }
    return true;
}

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
