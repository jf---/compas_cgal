#include "segment_site_endpoint_binding.h"

#include <algorithm>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

#include <CGAL/CORE/poly/Poly.h>
#include <CGAL/Polynomial_traits_d.h>

namespace {

using Algebraic = ExactAlgebraicKernel1::Algebraic_real_1;
using Polynomial = ExactAlgebraicKernel1::Polynomial_1;

struct FieldPolynomial2 {
  RationalPolynomial rational;
  RationalPolynomial radical;
};

void normalize(RationalPolynomial &polynomial) {
  trim(polynomial);
  if (polynomial.empty()) {
    polynomial.push_back(0);
  }
}

RationalPolynomial scaled(RationalPolynomial polynomial,
                          const CORE::BigRat &factor) {
  for (CORE::BigRat &coefficient : polynomial) {
    coefficient *= factor;
  }
  normalize(polynomial);
  return polynomial;
}

RationalPolynomial difference(RationalPolynomial lhs,
                              const RationalPolynomial &rhs) {
  add_in_place(lhs, scaled(rhs, -1));
  normalize(lhs);
  return lhs;
}

FieldPolynomial2 field_constant(const MatQuadraticFieldValue2 &value) {
  return {{value.rational}, {value.radical}};
}

FieldPolynomial2 field_coordinate(const std::vector<CORE::BigRat> &rational,
                                  const std::vector<CORE::BigRat> &radical) {
  return {rational, radical};
}

FieldPolynomial2 field_add(FieldPolynomial2 lhs, const FieldPolynomial2 &rhs) {
  add_in_place(lhs.rational, rhs.rational);
  add_in_place(lhs.radical, rhs.radical);
  normalize(lhs.rational);
  normalize(lhs.radical);
  return lhs;
}

FieldPolynomial2 field_subtract(FieldPolynomial2 lhs,
                                const FieldPolynomial2 &rhs) {
  lhs.rational = difference(lhs.rational, rhs.rational);
  lhs.radical = difference(lhs.radical, rhs.radical);
  return lhs;
}

FieldPolynomial2 field_scale(FieldPolynomial2 value,
                             const CORE::BigRat &factor) {
  value.rational = scaled(std::move(value.rational), factor);
  value.radical = scaled(std::move(value.radical), factor);
  return value;
}

FieldPolynomial2 field_multiply(const FieldPolynomial2 &lhs,
                                const FieldPolynomial2 &rhs,
                                const CORE::BigRat &radicand) {
  RationalPolynomial rational = multiply(lhs.rational, rhs.rational);
  add_in_place(rational, scaled(multiply(lhs.radical, rhs.radical), radicand));
  RationalPolynomial radical = multiply(lhs.rational, rhs.radical);
  add_in_place(radical, multiply(lhs.radical, rhs.rational));
  normalize(rational);
  normalize(radical);
  return {
      std::move(rational),
      std::move(radical),
  };
}

FieldPolynomial2 field_square(const FieldPolynomial2 &value,
                              const CORE::BigRat &radicand) {
  return field_multiply(value, value, radicand);
}

bool is_zero(const RationalPolynomial &polynomial) {
  RationalPolynomial normalized = polynomial;
  normalize(normalized);
  return normalized.size() == 1 && normalized.front() == 0;
}

Polynomial integer_polynomial(const RationalPolynomial &polynomial) {
  const std::vector<ExactAlgebraicInteger1> coefficients =
      primitive_integer_coefficients(polynomial);
  return typename CGAL::Polynomial_traits_d<Polynomial>::Construct_polynomial()(
      coefficients.begin(), coefficients.end());
}

CGAL::Sign rational_sign_at(const RationalPolynomial &polynomial,
                            const Algebraic &root,
                            const ExactAlgebraicKernel1 &kernel) {
  if (is_zero(polynomial)) {
    return CGAL::ZERO;
  }
  const CGAL::Sign normalized_sign =
      kernel.sign_at_1_object()(integer_polynomial(polynomial), root);
  return polynomial.back() < 0 ? CGAL::opposite(normalized_sign)
                               : normalized_sign;
}

CGAL::Sign field_sign_at(const FieldPolynomial2 &value,
                         const CORE::BigRat &radicand, const Algebraic &root,
                         const ExactAlgebraicKernel1 &kernel) {
  const CGAL::Sign rational_sign =
      rational_sign_at(value.rational, root, kernel);
  const CGAL::Sign radical_sign = rational_sign_at(value.radical, root, kernel);
  if (radical_sign == CGAL::ZERO) {
    return rational_sign;
  }
  if (rational_sign == CGAL::ZERO) {
    return radical_sign;
  }
  if (rational_sign == radical_sign) {
    return rational_sign;
  }
  RationalPolynomial magnitude = multiply(value.rational, value.rational);
  add_in_place(magnitude,
               scaled(multiply(value.radical, value.radical), -radicand));
  const CGAL::Sign magnitude_sign = rational_sign_at(magnitude, root, kernel);
  if (magnitude_sign == CGAL::ZERO) {
    return CGAL::ZERO;
  }
  return magnitude_sign == CGAL::POSITIVE ? rational_sign : radical_sign;
}

FieldPolynomial2 squared_distance(const FieldPolynomial2 &x,
                                  const FieldPolynomial2 &y,
                                  const MatExactPointSiteSource2 &point,
                                  const CORE::BigRat &radicand) {
  return field_add(
      field_square(field_subtract(x, field_constant(point.x)), radicand),
      field_square(field_subtract(y, field_constant(point.y)), radicand));
}

FieldPolynomial2 dot(const FieldPolynomial2 &ax, const FieldPolynomial2 &ay,
                     const FieldPolynomial2 &bx, const FieldPolynomial2 &by,
                     const CORE::BigRat &radicand) {
  return field_add(field_multiply(ax, bx, radicand),
                   field_multiply(ay, by, radicand));
}

CORE::Expr core_field_value(const MatQuadraticFieldValue2 &value,
                            const CORE::BigRat &radicand) {
  return CORE::Expr(value.rational) +
         CORE::sqrt(CORE::Expr(radicand)) * CORE::Expr(value.radical);
}

CORE::Expr core_evaluate(const FieldPolynomial2 &value,
                         const CORE::BigRat &radicand,
                         const CORE::Expr &parameter) {
  const auto evaluate = [&parameter](const RationalPolynomial &polynomial) {
    CORE::Expr result(0);
    for (auto coefficient = polynomial.rbegin();
         coefficient != polynomial.rend(); ++coefficient) {
      result *= parameter;
      result += CORE::Expr(*coefficient);
    }
    return result;
  };
  return evaluate(value.rational) +
         CORE::sqrt(CORE::Expr(radicand)) * evaluate(value.radical);
}

CORE::Expr
exact_core_root(const std::vector<ExactAlgebraicInteger1> &coefficients,
                const Algebraic &root, std::size_t ordinal) {
  if (root.is_rational()) {
    return CORE::Expr(root.rational());
  }
  return CORE::rootOf(CORE::Polynomial<CORE::BigInt>(coefficients),
                      static_cast<int>(ordinal + 1));
}

void require_same_field(const MatExactPointSiteSource2 &point,
                        const CORE::BigRat &radicand) {
  if (point.radicand != radicand) {
    throw InvalidRationalPrimitiveError(
        "endpoint sources use different quadratic fields");
  }
}

FieldPolynomial2 line_value(const MatExactPointSiteSource2 &point,
                            const MatExactOpenSegmentSource2 &segment) {
  return field_add(
      field_add(field_scale(field_constant(point.x), segment.line_a),
                field_scale(field_constant(point.y), segment.line_b)),
      field_constant({segment.line_c, 0}));
}

} // namespace

MatLiveParabolaEndpoint2 bind_point_limiter_parabola_endpoint(
    const MatExactPointSiteSource2 &focus,
    const MatExactOpenSegmentSource2 &segment,
    const MatExactPointSiteSource2 &segment_source,
    const MatExactPointSiteSource2 &segment_target,
    const MatExactPointSiteSource2 &limiter,
    const MatTraits::Point_2 &live_point,
    const SegmentSiteParabola2 &live_parabola) {
  const CORE::BigRat radicand = focus.radicand;
  require_same_field(segment_source, radicand);
  require_same_field(segment_target, radicand);
  require_same_field(limiter, radicand);
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
          radicand) != CGAL::ZERO) {
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
  RationalPolynomial norm = multiply(equation.rational, equation.rational);
  add_in_place(norm,
               scaled(multiply(equation.radical, equation.radical), -radicand));
  normalize(norm);
  if (is_zero(norm)) {
    throw UnboundLiveParabolaEndpointError(
        "point limiter equality is not zero-dimensional");
  }
  const std::vector<ExactAlgebraicInteger1> norm_coefficients =
      primitive_integer_coefficients(norm);
  const Polynomial norm_polynomial = integer_polynomial(norm);
  ExactAlgebraicKernel1 kernel;
  std::vector<Algebraic> roots;
  kernel.solve_1_object()(norm_polynomial, true, std::back_inserter(roots));

  const FieldPolynomial2 segment_dx = field_subtract(
      field_constant(segment_target.x), field_constant(segment_source.x));
  const FieldPolynomial2 segment_dy = field_subtract(
      field_constant(segment_target.y), field_constant(segment_source.y));
  const FieldPolynomial2 projection =
      dot(field_subtract(x, field_constant(segment_source.x)),
          field_subtract(y, field_constant(segment_source.y)), segment_dx,
          segment_dy, radicand);
  const FieldPolynomial2 segment_length_squared =
      dot(segment_dx, segment_dy, segment_dx, segment_dy, radicand);
  const FieldPolynomial2 projection_remainder =
      field_subtract(segment_length_squared, projection);

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
      ((live_point.x() - vertex_x) * CORE::Expr(-segment.line_b) +
       (live_point.y() - vertex_y) * CORE::Expr(segment.line_a)) /
      CORE::Expr(line_norm);

  std::vector<MatLiveParabolaEndpoint2> matches;
  std::size_t equation_roots = 0;
  std::size_t feature_roots = 0;
  std::size_t parameter_matches = 0;
  std::size_t point_matches = 0;
  for (std::size_t ordinal = 0; ordinal < roots.size(); ++ordinal) {
    const Algebraic &root = roots[ordinal];
    if (field_sign_at(equation, radicand, root, kernel) != CGAL::ZERO) {
      continue;
    }
    ++equation_roots;
    if (field_sign_at(projection, radicand, root, kernel) != CGAL::POSITIVE ||
        field_sign_at(projection_remainder, radicand, root, kernel) !=
            CGAL::POSITIVE) {
      continue;
    }
    ++feature_roots;
    const CORE::Expr parameter =
        exact_core_root(norm_coefficients, root, ordinal);
    if (parameter != live_parameter) {
      continue;
    }
    ++parameter_matches;
    const MatTraits::Point_2 candidate(core_evaluate(x, radicand, parameter),
                                       core_evaluate(y, radicand, parameter));
    if (candidate != live_point) {
      continue;
    }
    ++point_matches;
    const bool p1 = live_point == live_parabola.p1;
    const bool p2 = live_point == live_parabola.p2;
    if (!p1 && !p2) {
      continue;
    }
    matches.push_back({
        {
            root,
            {
                algebraic_root_identity_v1(root),
                limiter.stable_site_id,
                p1 ? "live-parabola/p1" : "live-parabola/p2",
            },
        },
        true,
    });
  }
  if (matches.size() != 1) {
    throw UnboundLiveParabolaEndpointError(
        "live point binds " + std::to_string(matches.size()) + " roots from " +
        std::to_string(roots.size()) + " norm, " +
        std::to_string(equation_roots) + " radical, " +
        std::to_string(feature_roots) + " feature, " +
        std::to_string(parameter_matches) + " parameter, " +
        std::to_string(point_matches) + " point matches");
  }
  return std::move(matches.front());
}
