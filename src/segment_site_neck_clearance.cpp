#include "segment_site_neck_clearance.h"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <utility>
#include <vector>

#include <CGAL/Polynomial_traits_d.h>

namespace {

constexpr std::size_t
    MAX_MAT_SQUARED_CLEARANCE_DEGREE = 4;
constexpr int
    SQUARED_WIDTH_FROM_CLEARANCE_FACTOR = 4;

using Algebraic =
    ExactAlgebraicKernel1::Algebraic_real_1;
using Polynomial1 =
    ExactAlgebraicKernel1::Polynomial_1;
using Polynomial2 =
    ExactAlgebraicKernel2::Polynomial_2;

struct ScaledIntegerPolynomial2 {
    std::vector<ExactAlgebraicInteger1>
        coefficients;
    ExactAlgebraicInteger1 denominator;
};

Polynomial1 polynomial_1(
    const RationalPolynomial& coefficients)
{
    const auto integer =
        primitive_integer_coefficients(
            coefficients);
    return typename CGAL::Polynomial_traits_d<
        Polynomial1>::Construct_polynomial()(
        integer.begin(),
        integer.end());
}

Polynomial1 sign_preserving_polynomial_1(
    const RationalPolynomial& coefficients)
{
    ExactAlgebraicInteger1 denominator = 1;
    for (const CORE::BigRat& coefficient :
         coefficients) {
        denominator *=
            CORE::denominator(coefficient);
    }
    std::vector<ExactAlgebraicInteger1>
        integer;
    integer.reserve(coefficients.size());
    for (const CORE::BigRat& coefficient :
         coefficients) {
        integer.push_back(
            CORE::numerator(coefficient)
            * CORE::div_exact(
                denominator,
                CORE::denominator(
                    coefficient)));
    }
    return typename CGAL::Polynomial_traits_d<
        Polynomial1>::Construct_polynomial()(
        integer.begin(),
        integer.end());
}

std::vector<Algebraic> real_roots(
    const RationalPolynomial& coefficients)
{
    if (coefficients.size() <= 1) {
        return {};
    }
    ExactAlgebraicKernel1 kernel;
    const Polynomial1 square_free =
        typename CGAL::Polynomial_traits_d<
            Polynomial1>::Canonicalize()(
            kernel.make_square_free_1_object()(
                polynomial_1(coefficients)));
    std::vector<Algebraic> roots;
    kernel.solve_1_object()(
        square_free,
        true,
        std::back_inserter(roots));
    return roots;
}

RationalPolynomial derivative(
    const RationalPolynomial& coefficients)
{
    if (coefficients.size() <= 1) {
        return {0};
    }
    RationalPolynomial result(
        coefficients.size() - 1,
        CORE::BigRat(0));
    for (std::size_t degree = 1;
         degree < coefficients.size();
         ++degree) {
        result[degree - 1] =
            CORE::BigRat(
                static_cast<unsigned long>(
                    degree))
            * coefficients[degree];
    }
    trim(result);
    return result;
}

CGAL::Sign opposite_sign(const CGAL::Sign sign)
{
    if (sign == CGAL::POSITIVE) {
        return CGAL::NEGATIVE;
    }
    if (sign == CGAL::NEGATIVE) {
        return CGAL::POSITIVE;
    }
    return CGAL::ZERO;
}

std::vector<Algebraic> relative_interior_roots(
    const RationalPolynomial& coefficients,
    const Algebraic& lower,
    const Algebraic& upper)
{
    ExactAlgebraicKernel1 kernel;
    const auto compare =
        kernel.compare_1_object();
    std::vector<Algebraic> interior;
    for (const Algebraic& root :
         real_roots(coefficients)) {
        if (compare(lower, root)
                == CGAL::SMALLER
            && compare(root, upper)
                == CGAL::SMALLER) {
            interior.push_back(root);
        }
    }
    return interior;
}

void require_nonnegative_on_profile(
    const RationalPolynomial& coefficients,
    const Algebraic& lower,
    const Algebraic& upper)
{
    ExactAlgebraicKernel1 kernel;
    const Polynomial1 exact_polynomial =
        sign_preserving_polynomial_1(
            coefficients);
    if (kernel.sign_at_1_object()(
            exact_polynomial,
            lower)
            == CGAL::NEGATIVE
        || kernel.sign_at_1_object()(
               exact_polynomial,
               upper)
            == CGAL::NEGATIVE) {
        throw NegativeMatSquaredClearanceError(
            "MAT squared clearance is negative at an edge endpoint");
    }

    std::vector<Algebraic> boundaries{lower};
    const std::vector<Algebraic> roots =
        relative_interior_roots(
            coefficients,
            lower,
            upper);
    boundaries.insert(
        boundaries.end(),
        roots.begin(),
        roots.end());
    boundaries.push_back(upper);
    for (std::size_t index = 0;
         index + 1 < boundaries.size();
         ++index) {
        const CORE::BigRat witness =
            kernel.bound_between_1_object()(
                boundaries[index],
                boundaries[index + 1]);
        if (evaluate_rational_polynomial(
                coefficients,
                witness)
            < 0) {
            throw NegativeMatSquaredClearanceError(
                "MAT squared clearance is negative inside an edge");
        }
    }
}

void require_profile_endpoint(
    const MatParameterEndpoint2& endpoint)
{
    if (!endpoint.parameter.has_value()
        || endpoint.provenance_ids.empty()
        || !std::is_sorted(
            endpoint.provenance_ids.begin(),
            endpoint.provenance_ids.end())
        || std::adjacent_find(
               endpoint.provenance_ids.begin(),
               endpoint.provenance_ids.end())
            != endpoint.provenance_ids.end()) {
        throw InvalidMatClearanceEdgeProfileError(
            "MAT clearance profile endpoint provenance is not canonical");
    }
    const MatEndpointExactEvidence2& evidence =
        endpoint.exact_evidence;
    if ((!evidence.original_voronoi_vertex
         && !evidence.domain_boundary
         && !evidence.clearance_root)
        || evidence.domain_boundary
            != !evidence.boundary_features.empty()
        || !std::is_sorted(
            evidence.boundary_features.begin(),
            evidence.boundary_features.end())
        || std::adjacent_find(
               evidence.boundary_features.begin(),
               evidence.boundary_features.end())
            != evidence.boundary_features.end()
        || !std::binary_search(
            endpoint.provenance_ids.begin(),
            endpoint.provenance_ids.end(),
            algebraic_root_identity_v1(
                *endpoint.parameter))) {
        throw InvalidMatClearanceEdgeProfileError(
            "MAT clearance profile endpoint evidence is incomplete");
    }
}

CGAL::Sign endpoint_inward_clearance_sign(
    const MatClearanceEdgeProfile2& profile,
    const bool lower)
{
    const RationalPolynomial derivative_function =
        derivative(
            profile.squared_clearance()
                .coefficients());
    const Algebraic& lower_parameter =
        *profile.lower().parameter;
    const Algebraic& upper_parameter =
        *profile.upper().parameter;
    const std::vector<Algebraic> critical =
        relative_interior_roots(
            derivative_function,
            lower_parameter,
            upper_parameter);
    ExactAlgebraicKernel1 kernel;
    const CORE::BigRat witness = lower
        ? kernel.bound_between_1_object()(
              lower_parameter,
              critical.empty()
                  ? upper_parameter
                  : critical.front())
        : kernel.bound_between_1_object()(
              critical.empty()
                  ? lower_parameter
                  : critical.back(),
              upper_parameter);
    const CGAL::Sign parameter_sign = CGAL::sign(
        evaluate_rational_polynomial(
            derivative_function,
            witness));
    return lower
        ? parameter_sign
        : opposite_sign(parameter_sign);
}

ScaledIntegerPolynomial2 scaled_integer_polynomial(
    const RationalPolynomial& coefficients)
{
    ExactAlgebraicInteger1 denominator = 1;
    for (const CORE::BigRat& coefficient :
         coefficients) {
        denominator *=
            CORE::denominator(coefficient);
    }
    std::vector<ExactAlgebraicInteger1>
        integer;
    integer.reserve(coefficients.size());
    for (const CORE::BigRat& coefficient :
         coefficients) {
        integer.push_back(
            CORE::numerator(coefficient)
            * CORE::div_exact(
                denominator,
                CORE::denominator(
                    coefficient)));
    }
    return {
        std::move(integer),
        std::move(denominator),
    };
}

Polynomial2 polynomial_2_in_variable(
    const Polynomial1& polynomial,
    const Polynomial2& variable)
{
    Polynomial2 result(
        ExactAlgebraicInteger1(0));
    Polynomial2 power(
        ExactAlgebraicInteger1(1));
    const int degree = CGAL::degree(polynomial);
    for (int index = 0;
         index <= degree;
         ++index) {
        result += polynomial[index] * power;
        power *= variable;
    }
    return result;
}

Polynomial2 polynomial_2_in_variable(
    const std::vector<
        ExactAlgebraicInteger1>& coefficients,
    const Polynomial2& variable)
{
    Polynomial2 result(
        ExactAlgebraicInteger1(0));
    Polynomial2 power(
        ExactAlgebraicInteger1(1));
    for (const ExactAlgebraicInteger1&
             coefficient :
         coefficients) {
        result += coefficient * power;
        power *= variable;
    }
    return result;
}

} // namespace

MatExactSquaredWidth2
MatExactSquaredWidth2::from_clearance(
    const MatSquaredClearanceFunction2&
        clearance,
    const Algebraic& parameter)
{
    RationalPolynomial width =
        clearance.coefficients();
    for (CORE::BigRat& coefficient : width) {
        coefficient *=
            SQUARED_WIDTH_FROM_CLEARANCE_FACTOR;
    }
    const ScaledIntegerPolynomial2 scaled =
        scaled_integer_polynomial(width);
    const Polynomial2 parameter_variable =
        CGAL::shift(
            Polynomial2(
                ExactAlgebraicInteger1(1)),
            1,
            0);
    const Polynomial2 width_variable =
        CGAL::shift(
            Polynomial2(
                ExactAlgebraicInteger1(1)),
            1,
            1);
    ExactAlgebraicKernel1 kernel1;
    const Polynomial2 parameter_equation =
        polynomial_2_in_variable(
            kernel1.compute_polynomial_1_object()(
                parameter),
            parameter_variable);
    const Polynomial2 width_equation =
        polynomial_2_in_variable(
            scaled.coefficients,
            parameter_variable)
        - scaled.denominator * width_variable;

    ExactAlgebraicKernel2 kernel2;
    std::vector<
        std::pair<
            ExactAlgebraicKernel2::Algebraic_real_2,
            ExactAlgebraicKernel2::Multiplicity_type>>
        solutions;
    kernel2.solve_2_object()(
        parameter_equation,
        width_equation,
        std::back_inserter(solutions));
    std::vector<Algebraic> matching;
    for (const auto& [solution, multiplicity] :
         solutions) {
        static_cast<void>(multiplicity);
        if (kernel1.compare_1_object()(
                kernel2.compute_x_2_object()(
                    solution),
                parameter)
            != CGAL::EQUAL) {
            continue;
        }
        const Algebraic value =
            kernel2.compute_y_2_object()(
                solution);
        if (matching.empty()
            || kernel1.compare_1_object()(
                   matching.back(),
                   value)
                != CGAL::EQUAL) {
            matching.push_back(value);
        }
    }
    if (matching.size() != 1) {
        throw IncompleteMatSquaredWidthError(
            "MAT critical parameter does not map to one exact squared width");
    }
    if (kernel1.compare_1_object()(
            matching.front(),
            CORE::BigRat(0))
        == CGAL::SMALLER) {
        throw NegativeMatSquaredClearanceError(
            "MAT exact squared width is negative");
    }
    return MatExactSquaredWidth2(
        matching.front(),
        algebraic_root_identity_v1(
            matching.front()));
}

MatSquaredClearanceFunction2::
    MatSquaredClearanceFunction2(
        RationalPolynomial coefficients)
    : coefficients_(
          std::move(coefficients))
{
}

MatSquaredClearanceFunction2
MatSquaredClearanceFunction2::build(
    RationalPolynomial coefficients)
{
    if (coefficients.empty()) {
        throw InvalidMatSquaredClearanceFunctionError(
            "MAT squared-clearance function has no coefficients");
    }
    trim(coefficients);
    if (coefficients.size() - 1
        > MAX_MAT_SQUARED_CLEARANCE_DEGREE) {
        throw UnsupportedMatSquaredClearanceDegreeError(
            "MAT squared-clearance polynomial exceeds degree four");
    }
    return MatSquaredClearanceFunction2(
        std::move(coefficients));
}

const RationalPolynomial&
MatSquaredClearanceFunction2::coefficients()
    const noexcept
{
    return coefficients_;
}

bool MatSquaredClearanceFunction2::is_constant()
    const noexcept
{
    return coefficients_.size() == 1;
}

MatClearanceEdgeProfile2::
    MatClearanceEdgeProfile2(
        std::string edge_id,
        std::vector<std::string>
            defining_site_ids,
        MatParameterEndpoint2 lower,
        MatParameterEndpoint2 upper,
        MatSquaredClearanceFunction2
            squared_clearance)
    : edge_id_(std::move(edge_id)),
      defining_site_ids_(
          std::move(defining_site_ids)),
      lower_(std::move(lower)),
      upper_(std::move(upper)),
      squared_clearance_(
          std::move(squared_clearance))
{
}

MatClearanceEdgeProfile2
MatClearanceEdgeProfile2::build(
    const MatExactGraphEdge2& edge,
    MatSquaredClearanceFunction2
        squared_clearance)
{
    if (edge.edge_id.empty()
        || edge.original_dual_id.empty()
        || edge.source_node_id.empty()
        || edge.target_node_id.empty()
        || edge.source_node_id
            == edge.target_node_id
        || edge.generator_site_ids.size()
            != 2
        || edge.generator_site_ids[0].empty()
        || edge.generator_site_ids[1].empty()
        || edge.generator_site_ids[0]
            >= edge.generator_site_ids[1]
        || edge.clip_component_index < 0
        || !edge.admissible_center_component) {
        throw InvalidMatClearanceEdgeProfileError(
            "MAT clearance profile graph edge is not canonical");
    }
    require_profile_endpoint(
        edge.source_endpoint);
    require_profile_endpoint(
        edge.target_endpoint);
    ExactAlgebraicKernel1 kernel;
    if (kernel.compare_1_object()(
            *edge.source_endpoint.parameter,
            *edge.target_endpoint.parameter)
        != CGAL::SMALLER) {
        throw InvalidMatClearanceEdgeProfileError(
            "MAT clearance profile parameter interval is not increasing");
    }
    require_nonnegative_on_profile(
        squared_clearance.coefficients(),
        *edge.source_endpoint.parameter,
        *edge.target_endpoint.parameter);
    return MatClearanceEdgeProfile2(
        edge.edge_id,
        edge.generator_site_ids,
        edge.source_endpoint,
        edge.target_endpoint,
        std::move(squared_clearance));
}

const std::string&
MatClearanceEdgeProfile2::edge_id()
    const noexcept
{
    return edge_id_;
}

const std::vector<std::string>&
MatClearanceEdgeProfile2::defining_site_ids()
    const noexcept
{
    return defining_site_ids_;
}

const MatParameterEndpoint2&
MatClearanceEdgeProfile2::lower()
    const noexcept
{
    return lower_;
}

const MatParameterEndpoint2&
MatClearanceEdgeProfile2::upper()
    const noexcept
{
    return upper_;
}

const MatSquaredClearanceFunction2&
MatClearanceEdgeProfile2::squared_clearance()
    const noexcept
{
    return squared_clearance_;
}

MatExactSquaredWidth2::MatExactSquaredWidth2(
    Algebraic value,
    std::string root_id)
    : value_(std::move(value)),
      root_id_(std::move(root_id))
{
}

const Algebraic&
MatExactSquaredWidth2::value() const noexcept
{
    return value_;
}

const std::string&
MatExactSquaredWidth2::root_id() const noexcept
{
    return root_id_;
}

MatClearanceEndpointBehavior2::
    MatClearanceEndpointBehavior2(
        const CGAL::Sign inward_clearance_sign,
        MatExactSquaredWidth2 squared_width)
    : inward_clearance_sign_(
          inward_clearance_sign),
      squared_width_(
          std::move(squared_width))
{
}

CGAL::Sign MatClearanceEndpointBehavior2::
    inward_clearance_sign() const noexcept
{
    return inward_clearance_sign_;
}

const MatExactSquaredWidth2&
MatClearanceEndpointBehavior2::
    squared_width() const noexcept
{
    return squared_width_;
}

MatClearanceEndpointBehavior2
lower_endpoint_clearance_behavior(
    const MatClearanceEdgeProfile2& profile)
{
    return MatClearanceEndpointBehavior2(
        endpoint_inward_clearance_sign(
            profile,
            true),
        MatExactSquaredWidth2::from_clearance(
            profile.squared_clearance(),
            *profile.lower().parameter));
}

MatClearanceEndpointBehavior2
upper_endpoint_clearance_behavior(
    const MatClearanceEdgeProfile2& profile)
{
    return MatClearanceEndpointBehavior2(
        endpoint_inward_clearance_sign(
            profile,
            false),
        MatExactSquaredWidth2::from_clearance(
            profile.squared_clearance(),
            *profile.upper().parameter));
}

MatStrictEdgeClearanceMinimum2::
    MatStrictEdgeClearanceMinimum2(
        std::string edge_id,
        std::vector<std::string>
            defining_site_ids,
        Algebraic parameter,
        std::string parameter_root_id,
        MatExactSquaredWidth2 squared_width)
    : edge_id_(std::move(edge_id)),
      defining_site_ids_(
          std::move(defining_site_ids)),
      parameter_(std::move(parameter)),
      parameter_root_id_(
          std::move(parameter_root_id)),
      squared_width_(
          std::move(squared_width))
{
}

const std::string&
MatStrictEdgeClearanceMinimum2::edge_id()
    const noexcept
{
    return edge_id_;
}

const std::vector<std::string>&
MatStrictEdgeClearanceMinimum2::
    defining_site_ids() const noexcept
{
    return defining_site_ids_;
}

const Algebraic&
MatStrictEdgeClearanceMinimum2::parameter()
    const noexcept
{
    return parameter_;
}

const std::string&
MatStrictEdgeClearanceMinimum2::
    parameter_root_id() const noexcept
{
    return parameter_root_id_;
}

CGAL::Sign MatStrictEdgeClearanceMinimum2::
    left_derivative_sign() const noexcept
{
    return CGAL::NEGATIVE;
}

CGAL::Sign MatStrictEdgeClearanceMinimum2::
    right_derivative_sign() const noexcept
{
    return CGAL::POSITIVE;
}

const MatExactSquaredWidth2&
MatStrictEdgeClearanceMinimum2::
    squared_width() const noexcept
{
    return squared_width_;
}

std::vector<MatStrictEdgeClearanceMinimum2>
strict_edge_clearance_minima(
    const MatClearanceEdgeProfile2& profile)
{
    const RationalPolynomial derivative_function =
        derivative(
            profile.squared_clearance()
                .coefficients());
    if (derivative_function.size() == 1) {
        return {};
    }
    const std::vector<Algebraic> critical =
        relative_interior_roots(
            derivative_function,
            *profile.lower().parameter,
            *profile.upper().parameter);
    ExactAlgebraicKernel1 kernel;
    std::vector<MatStrictEdgeClearanceMinimum2>
        minima;
    for (std::size_t index = 0;
         index < critical.size();
         ++index) {
        const Algebraic& left =
            index == 0
            ? *profile.lower().parameter
            : critical[index - 1];
        const Algebraic& right =
            index + 1 == critical.size()
            ? *profile.upper().parameter
            : critical[index + 1];
        const CORE::BigRat left_witness =
            kernel.bound_between_1_object()(
                left,
                critical[index]);
        const CORE::BigRat right_witness =
            kernel.bound_between_1_object()(
                critical[index],
                right);
        const CGAL::Sign left_sign = CGAL::sign(
            evaluate_rational_polynomial(
                derivative_function,
                left_witness));
        const CGAL::Sign right_sign = CGAL::sign(
            evaluate_rational_polynomial(
                derivative_function,
                right_witness));
        if (left_sign != CGAL::NEGATIVE
            || right_sign != CGAL::POSITIVE) {
            continue;
        }
        minima.push_back(
            MatStrictEdgeClearanceMinimum2(
                profile.edge_id(),
                profile.defining_site_ids(),
                critical[index],
                algebraic_root_identity_v1(
                    critical[index]),
                MatExactSquaredWidth2::
                    from_clearance(
                        profile.squared_clearance(),
                        critical[index])));
    }
    return minima;
}
