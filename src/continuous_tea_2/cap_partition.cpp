#include "cap_partition.h"

#include "../exact_algebraic_1.h"
#include "event_certificate.h"
#include "event_partition.h"
#include "parameter_charts.h"

#include <algorithm>
#include <sstream>
#include <string_view>
#include <utility>

#include <CGAL/CORE/BigRat.h>
#include <CGAL/Polynomial_traits_d.h>

namespace {

using Integer = CORE::BigInt;
using Rational = CORE::BigRat;
using Polynomial = ExactAlgebraicKernel1::Polynomial_1;
using AlgebraicReal = ExactAlgebraicKernel1::Algebraic_real_1;
using RationalPolynomial = std::vector<Rational>;

Integer greatest_common_divisor(Integer left, Integer right)
{
    left = left < 0 ? -left : left;
    right = right < 0 ? -right : right;
    while (right != 0) {
        const Integer remainder = left % right;
        left = right;
        right = remainder;
    }
    return left;
}

Integer least_common_multiple(
    const Integer& left,
    const Integer& right)
{
    return left / greatest_common_divisor(left, right)
        * right;
}

Rational parse_rational(
    const std::string& text,
    std::string_view role)
{
    const std::size_t separator = text.find('/');
    try {
        if (separator == std::string::npos) {
            return Rational(Integer(text));
        }
        if (text.find('/', separator + 1)
            != std::string::npos) {
            throw InvalidAlgebraicPolynomialError(
                std::string(role)
                + " is not an exact rational");
        }
        const Integer numerator(
            text.substr(0, separator));
        const Integer denominator(
            text.substr(separator + 1));
        if (denominator == 0) {
            throw InvalidAlgebraicPolynomialError(
                std::string(role)
                + " has zero denominator");
        }
        return Rational(numerator, denominator);
    } catch (const EventSubstrateError&) {
        throw;
    } catch (const std::exception&) {
        throw InvalidAlgebraicPolynomialError(
            std::string(role)
            + " is not an exact rational");
    }
}

void trim(RationalPolynomial& polynomial)
{
    while (polynomial.size() > 1
           && polynomial.back() == 0) {
        polynomial.pop_back();
    }
}

RationalPolynomial parse_polynomial(
    const std::vector<std::string>& coefficients,
    std::string_view role)
{
    if (coefficients.empty()) {
        throw InvalidAlgebraicPolynomialError(
            std::string(role)
            + " requires coefficients");
    }
    RationalPolynomial result;
    result.reserve(coefficients.size());
    for (const std::string& coefficient : coefficients) {
        result.push_back(
            parse_rational(coefficient, role));
    }
    trim(result);
    return result;
}

Rational coefficient(
    const RationalPolynomial& polynomial,
    std::size_t index)
{
    return index < polynomial.size()
        ? polynomial[index]
        : Rational(0);
}

RationalPolynomial add(
    const RationalPolynomial& left,
    const RationalPolynomial& right)
{
    RationalPolynomial result(
        std::max(left.size(), right.size()),
        Rational(0));
    for (std::size_t index = 0;
         index < result.size();
         ++index) {
        result[index] =
            coefficient(left, index)
            + coefficient(right, index);
    }
    trim(result);
    return result;
}

RationalPolynomial scale(
    const RationalPolynomial& polynomial,
    const Rational& factor)
{
    RationalPolynomial result = polynomial;
    for (Rational& coefficient_value : result) {
        coefficient_value *= factor;
    }
    trim(result);
    return result;
}

RationalPolynomial subtract(
    const RationalPolynomial& left,
    const RationalPolynomial& right)
{
    return add(left, scale(right, -1));
}

RationalPolynomial multiply(
    const RationalPolynomial& left,
    const RationalPolynomial& right)
{
    RationalPolynomial result(
        left.size() + right.size() - 1,
        Rational(0));
    for (std::size_t left_index = 0;
         left_index < left.size();
         ++left_index) {
        for (std::size_t right_index = 0;
             right_index < right.size();
             ++right_index) {
            result[left_index + right_index] +=
                left[left_index] * right[right_index];
        }
    }
    trim(result);
    return result;
}

Rational evaluate(
    const RationalPolynomial& polynomial,
    const Rational& parameter)
{
    Rational result = 0;
    for (std::size_t index = polynomial.size();
         index-- > 0;) {
        result *= parameter;
        result += polynomial[index];
    }
    return result;
}

struct RationalCircleVector {
    RationalPolynomial x_numerator;
    RationalPolynomial y_numerator;
    RationalPolynomial denominator;
};

RationalCircleVector circle_vector(
    const RationalPolynomial& numerator,
    const RationalPolynomial& denominator)
{
    const RationalPolynomial numerator_squared =
        multiply(numerator, numerator);
    const RationalPolynomial denominator_squared =
        multiply(denominator, denominator);
    return {
        subtract(
            denominator_squared,
            numerator_squared),
        scale(
            multiply(numerator, denominator),
            2),
        add(
            denominator_squared,
            numerator_squared),
    };
}

struct PairPolynomials {
    RationalPolynomial cap;
    RationalPolynomial determinant;
    RationalPolynomial dot;
};

PairPolynomials pair_polynomials(
    const RationalPolynomial& first_numerator,
    const RationalPolynomial& first_denominator,
    const RationalPolynomial& second_numerator,
    const RationalPolynomial& second_denominator,
    const Rational& cap)
{
    const RationalCircleVector first =
        circle_vector(
            first_numerator,
            first_denominator);
    const RationalCircleVector second =
        circle_vector(
            second_numerator,
            second_denominator);
    const RationalPolynomial difference =
        subtract(
            multiply(
                first_numerator,
                second_denominator),
            multiply(
                second_numerator,
                first_denominator));
    return {
        subtract(
            scale(
                multiply(difference, difference),
                4),
            scale(
                multiply(
                    first.denominator,
                    second.denominator),
                cap)),
        subtract(
            multiply(
                first.x_numerator,
                second.y_numerator),
            multiply(
                first.y_numerator,
                second.x_numerator)),
        add(
            multiply(
                first.x_numerator,
                second.x_numerator),
            multiply(
                first.y_numerator,
                second.y_numerator)),
    };
}

std::vector<Integer> primitive_integers(
    const RationalPolynomial& polynomial,
    bool positive_leading = true)
{
    Integer denominator_lcm = 1;
    for (const Rational& value : polynomial) {
        denominator_lcm = least_common_multiple(
            denominator_lcm,
            CORE::denominator(value));
    }
    std::vector<Integer> result;
    result.reserve(polynomial.size());
    Integer divisor = 0;
    for (const Rational& value : polynomial) {
        const Integer integer = CORE::numerator(
            value * Rational(denominator_lcm));
        result.push_back(integer);
        divisor = greatest_common_divisor(
            divisor,
            integer);
    }
    if (divisor == 0) {
        return {0};
    }
    for (Integer& value : result) {
        value /= divisor;
    }
    if (positive_leading && result.back() < 0) {
        for (Integer& value : result) {
            value = -value;
        }
    }
    return result;
}

std::vector<std::string> coefficient_text(
    const std::vector<Integer>& coefficients)
{
    std::vector<std::string> result;
    result.reserve(coefficients.size());
    for (const Integer& coefficient_value : coefficients) {
        result.push_back(
            coefficient_value.convert_to<std::string>());
    }
    return result;
}

Polynomial polynomial_from_integers(
    const std::vector<Integer>& coefficients)
{
    return typename CGAL::Polynomial_traits_d<
        Polynomial>::Construct_polynomial()(
        coefficients.begin(),
        coefficients.end());
}

CGAL::Sign sign_at_root(
    const RationalPolynomial& polynomial,
    const AlgebraicRootRecord2& root)
{
    const Polynomial tested =
        polynomial_from_integers(
            primitive_integers(polynomial, false));
    std::vector<Integer> factor_coefficients;
    factor_coefficients.reserve(
        root.factor_coefficients.size());
    for (const std::string& coefficient_value :
         root.factor_coefficients) {
        factor_coefficients.emplace_back(
            coefficient_value);
    }
    const Polynomial factor =
        polynomial_from_integers(
            factor_coefficients);
    ExactAlgebraicKernel1 kernel;
    std::vector<AlgebraicReal> roots;
    kernel.solve_1_object()(
        factor,
        true,
        std::back_inserter(roots));
    if (root.root_ordinal >= roots.size()) {
        throw AlgebraicRootIsolationError(
            "root ordinal exceeds irreducible factor roots");
    }
    return kernel.sign_at_1_object()(
        tested,
        roots[root.root_ordinal]);
}

std::string sign_name(const Rational& value)
{
    return value > 0
        ? "positive"
        : value < 0
        ? "negative"
        : "zero";
}

std::string disposition(
    const Rational& determinant,
    const Rational& dot)
{
    if (determinant > 0) {
        return "ccw";
    }
    if (determinant < 0) {
        return "clockwise-complement";
    }
    if (dot > 0) {
        return "zero";
    }
    if (dot < 0) {
        return "pi";
    }
    return "undefined";
}

std::vector<ChartSeam2> chart_seams(
    const std::vector<ParameterChart2>& charts)
{
    std::vector<ChartSeam2> seams;
    seams.reserve(charts.size());
    for (const ParameterChart2& chart : charts) {
        seams.push_back(
            {
                chart.start_seam_id,
                chart.chart_id,
            });
    }
    return seams;
}

std::string source_payload(
    const std::vector<std::string>& first_numerator,
    const std::vector<std::string>& first_denominator,
    const std::vector<std::string>& second_numerator,
    const std::vector<std::string>& second_denominator,
    const std::string& cap_numerator,
    const std::string& cap_denominator,
    const PartitionEvent2& event)
{
    return encode_string_sequence(
        {
            encode_string_sequence(first_numerator),
            encode_string_sequence(first_denominator),
            encode_string_sequence(second_numerator),
            encode_string_sequence(second_denominator),
            cap_numerator,
            cap_denominator,
            event.kind,
            event.feature_id,
            event.support_id,
            event.trim_id,
            event.vertex_id,
            event.branch_id,
            event.disposition,
        });
}

void assign_cap_projection_contract(
    ProjectionRecord2& projection,
    int actual_degree,
    const std::vector<std::string>& coefficients)
{
    projection.actual_motion_degree = actual_degree;
    projection.actual_rim_degree = 0;
    projection.bound_motion_degree = 4;
    projection.bound_rim_degree = 0;
    projection.degree_bound_id =
        "cap-resultant-(4,0)-v1";
    projection.normalized_coefficient_bytes =
        encode_string_sequence(coefficients);
    std::sort(
        projection.factor_coefficients.begin(),
        projection.factor_coefficients.end(),
        [](const auto& left, const auto& right) {
            const Integer left_leading(left.back());
            const Integer right_leading(right.back());
            if (left_leading != right_leading) {
                return left_leading < right_leading;
            }
            return Integer(left.front())
                < Integer(right.front());
        });
}

void assign_rational_cap_witnesses(
    EventPartitionCertificate2& certificate)
{
    std::vector<Rational> boundaries{Rational(0)};
    for (const AlgebraicRootRecord2& root :
         certificate.roots) {
        if (root.factor_coefficients.size() != 2) {
            return;
        }
        boundaries.emplace_back(
            -Integer(root.factor_coefficients[0]),
            Integer(root.factor_coefficients[1]));
    }
    boundaries.emplace_back(1);
    if (certificate.cells.size() + 1
        != boundaries.size()) {
        throw EventPartitionVerificationError(
            "cap cell count does not match rational roots");
    }
    for (std::size_t index = 0;
         index < certificate.cells.size();
         ++index) {
        const Rational witness =
            (boundaries[index]
             + boundaries[index + 1])
            / Rational(2);
        certificate.cells[index].witness_numerator =
            CORE::numerator(witness)
                .convert_to<std::string>();
        certificate.cells[index].witness_denominator =
            CORE::denominator(witness)
                .convert_to<std::string>();
    }
}

} // namespace

CcwOrientation2 classify_ccw_orientation(
    const std::vector<std::string>& first_numerator,
    const std::vector<std::string>& first_denominator,
    const std::vector<std::string>& second_numerator,
    const std::vector<std::string>& second_denominator,
    const std::string& motion_parameter)
{
    const RationalPolynomial first_num =
        parse_polynomial(
            first_numerator,
            "first branch numerator");
    const RationalPolynomial first_den =
        parse_polynomial(
            first_denominator,
            "first branch denominator");
    const RationalPolynomial second_num =
        parse_polynomial(
            second_numerator,
            "second branch numerator");
    const RationalPolynomial second_den =
        parse_polynomial(
            second_denominator,
            "second branch denominator");
    const Rational parameter =
        parse_rational(
            motion_parameter,
            "motion parameter");
    const Rational first_denominator_value =
        evaluate(first_den, parameter);
    const Rational second_denominator_value =
        evaluate(second_den, parameter);
    if (first_denominator_value == 0
        || second_denominator_value == 0) {
        throw InvalidAlgebraicPolynomialError(
            "branch denominator vanishes at motion parameter");
    }
    const PairPolynomials pair = pair_polynomials(
        first_num,
        first_den,
        second_num,
        second_den,
        Rational(0));
    const Rational determinant =
        evaluate(pair.determinant, parameter);
    const Rational dot = evaluate(pair.dot, parameter);
    const Rational first_orientation =
        first_denominator_value > 0
        ? Rational(1)
        : Rational(-1);
    const Rational second_orientation =
        second_denominator_value > 0
        ? Rational(1)
        : Rational(-1);
    const Rational oriented_determinant =
        determinant * first_orientation
        * second_orientation;
    const Rational oriented_dot =
        dot * first_orientation
        * second_orientation;
    return {
        disposition(
            oriented_determinant,
            oriented_dot),
        sign_name(oriented_determinant),
        sign_name(oriented_dot),
    };
}

EventPartitionCertificate2 partition_cap_crossings(
    const std::vector<std::string>& first_numerator,
    const std::vector<std::string>& first_denominator,
    const std::vector<std::string>& second_numerator,
    const std::vector<std::string>& second_denominator,
    const std::string& cap_numerator,
    const std::string& cap_denominator,
    const PartitionEvent2& event)
{
    const Rational cap =
        parse_rational(cap_numerator, "cap numerator")
        / parse_rational(cap_denominator, "cap denominator");
    if (cap <= 0) {
        throw InvalidAlgebraicPolynomialError(
            "squared chord cap must be positive");
    }
    const RationalPolynomial first_num =
        parse_polynomial(
            first_numerator,
            "first branch numerator");
    const RationalPolynomial first_den =
        parse_polynomial(
            first_denominator,
            "first branch denominator");
    const RationalPolynomial second_num =
        parse_polynomial(
            second_numerator,
            "second branch numerator");
    const RationalPolynomial second_den =
        parse_polynomial(
            second_denominator,
            "second branch denominator");
    const PairPolynomials pair = pair_polynomials(
        first_num,
        first_den,
        second_num,
        second_den,
        cap);
    const std::vector<Integer> cap_integers =
        primitive_integers(pair.cap);
    const std::vector<std::string> cap_coefficients =
        coefficient_text(cap_integers);

    EventPartitionCertificate2 certificate;
    if (cap_integers.size() == 1
        && cap_integers.front() == 0) {
        certificate.build_evidence =
            exact_algebraic_backend_evidence();
        certificate.charts = parameter_charts();
        certificate.projections = {
            {
                "cap-resultant",
                {cap_coefficients},
                {},
                -1,
                -1,
                4,
                0,
                "cap-resultant-(4,0)-v1",
                encode_string_sequence(cap_coefficients),
            },
        };
        const CcwOrientation2 orientation =
            classify_ccw_orientation(
                first_numerator,
                first_denominator,
                second_numerator,
                second_denominator,
                "1/2");
        certificate.overlaps.push_back(
            {
                orientation.disposition == "ccw"
                    ? "identically-equal-cap-interval"
                    : "filtered-cap-complement",
                "0",
                "1",
                "1",
                "2",
                orientation.disposition,
            });
    } else {
        certificate = partition_projections(
            {
                {
                    "cap-resultant",
                    cap_coefficients,
                    {event},
                },
            });
        assign_cap_projection_contract(
            certificate.projections.front(),
            static_cast<int>(cap_integers.size()) - 1,
            cap_coefficients);

        bool any_accepted = false;
        bool any_rejected = false;
        for (const AlgebraicRootRecord2& root :
             certificate.roots) {
            const CGAL::Sign determinant =
                sign_at_root(pair.determinant, root);
            const CGAL::Sign dot =
                sign_at_root(pair.dot, root);
            const bool accepted =
                determinant == CGAL::POSITIVE
                || (determinant == CGAL::ZERO
                    && dot != CGAL::NEGATIVE);
            any_accepted = any_accepted || accepted;
            any_rejected = any_rejected || !accepted;
        }
        if (any_accepted && any_rejected) {
            throw EventPartitionVerificationError(
                "cap roots have mixed orientation dispositions");
        }
        if (any_rejected) {
            certificate.roots.clear();
            certificate.cells.clear();
            certificate.fibres.clear();
            certificate.overlaps.push_back(
                {
                    "filtered-cap-complement",
                    "0",
                    "1",
                    "1",
                    "2",
                    "clockwise-complement",
                });
        } else {
            assign_rational_cap_witnesses(certificate);
            for (EventFibre2& fibre :
                 certificate.fibres) {
                for (PartitionEvent2& accepted_event :
                     fibre.events) {
                    accepted_event
                        .original_equations_rechecked = true;
                    accepted_event
                        .orientation_rechecked = true;
                    accepted_event.trim_disposition =
                        "accepted";
                }
            }
            for (ParameterCell2& cell : certificate.cells) {
                const Rational witness(
                    Integer(cell.witness_numerator),
                    Integer(cell.witness_denominator));
                cell.disposition =
                    evaluate(pair.cap, witness) > 0
                    ? "cap-exceeds"
                    : "below-cap";
            }
        }
    }

    certificate.seams = chart_seams(certificate.charts);
    certificate.source_kind = "cap-crossings-v1";
    certificate.source_payload = source_payload(
        first_numerator,
        first_denominator,
        second_numerator,
        second_denominator,
        cap_numerator,
        cap_denominator,
        event);
    finalize_event_partition(certificate);
    return certificate;
}
