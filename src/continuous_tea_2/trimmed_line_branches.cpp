#include "event_partition.h"

#include "event_partition_internal.h"
#include "segment_source.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <iterator>
#include <optional>
#include <string_view>
#include <utility>

#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Fraction_traits.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/number_utils.h>

namespace {

using continuous_tea_2::event_partition_internal::
    Integer;
using continuous_tea_2::event_partition_internal::
    Polynomial;
using continuous_tea_2::event_partition_internal::
    Rational;
using continuous_tea_2::event_partition_internal::
    coefficient_text;
using continuous_tea_2::event_partition_internal::
    parse_rational;
using continuous_tea_2::event_partition_internal::
    parse_values;
using continuous_tea_2::event_partition_internal::
    polynomial_from_integers;
using continuous_tea_2::event_partition_internal::
    primitive_coefficients;
using continuous_tea_2::event_partition_internal::
    rational_text;

using Kernel = ExactAlgebraicKernel1;
using AlgebraicReal = Kernel::Algebraic_real_1;
using IntegerRows = std::vector<std::vector<Integer>>;
using RationalRows = std::vector<std::vector<Rational>>;
using RationalPolynomial =
    CGAL::Polynomial<Rational>;
using RationalRowFamily =
    CGAL::Polynomial<RationalPolynomial>;
using QuadraticCoefficients =
    std::array<Rational, 3>;

struct TrimRootCandidate {
    AlgebraicReal value;
    AlgebraicRootRecord2 record;
    std::optional<Rational> rational_value;
};

std::string stable_record(
    std::string_view tag,
    const std::vector<std::string>& fields)
{
    std::string result(tag);
    result.push_back('\0');
    for (const std::string& field : fields) {
        const std::uint64_t size = field.size();
        for (int shift = 56; shift >= 0; shift -= 8) {
            result.push_back(
                static_cast<char>(
                    (size >> shift) & 0xffU));
        }
        result += field;
    }
    return result;
}

template <typename Coefficients>
Polynomial polynomial_from_rationals(
    const Coefficients& coefficients)
{
    using FractionTraits =
        CGAL::Fraction_traits<RationalPolynomial>;
    const RationalPolynomial rational_polynomial(
        coefficients.begin(),
        coefficients.end());
    typename FractionTraits::Numerator_type numerator;
    typename FractionTraits::Denominator_type denominator;
    typename FractionTraits::Decompose()(
        rational_polynomial,
        numerator,
        denominator);
    if (denominator <= 0) {
        throw AlgebraicBackendError(
            "polynomial fraction decomposition returned a non-positive denominator");
    }
    return numerator;
}

IntegerRows signed_primitive_rows(
    const RationalRows& rows)
{
    std::vector<RationalPolynomial>
        rational_row_polynomials;
    rational_row_polynomials.reserve(rows.size());
    for (const auto& row : rows) {
        rational_row_polynomials.emplace_back(
            row.begin(),
            row.end());
    }
    const RationalRowFamily rational_family(
        rational_row_polynomials.begin(),
        rational_row_polynomials.end());
    using FamilyFractionTraits =
        CGAL::Fraction_traits<RationalRowFamily>;
    typename FamilyFractionTraits::Numerator_type
        integer_family;
    typename FamilyFractionTraits::Denominator_type
        common_denominator;
    typename FamilyFractionTraits::Decompose()(
        rational_family,
        integer_family,
        common_denominator);
    if (common_denominator <= 0) {
        throw AlgebraicBackendError(
            "predicate-family fraction decomposition returned a non-positive denominator");
    }

    IntegerRows integers;
    integers.reserve(rows.size());
    Integer divisor = 0;
    using IntegerTraits =
        CGAL::Algebraic_structure_traits<Integer>;
    const typename IntegerTraits::Gcd gcd;
    const typename IntegerTraits::Integral_division
        exact_divide;
    const int family_degree =
        CGAL::degree(integer_family);
    for (std::size_t row_index = 0;
         row_index < rows.size();
         ++row_index) {
        std::vector<Integer> integer_row(
            rows[row_index].size(),
            Integer(0));
        if (family_degree >= 0
            && row_index
            <= static_cast<std::size_t>(
                family_degree)) {
            const Polynomial& integer_row_polynomial =
                integer_family[
                    static_cast<unsigned int>(
                        row_index)];
            const int row_degree =
                CGAL::degree(
                    integer_row_polynomial);
            for (std::size_t column_index = 0;
                 column_index < integer_row.size();
                 ++column_index) {
                if (row_degree >= 0
                    && column_index
                    <= static_cast<std::size_t>(
                        row_degree)) {
                    integer_row[column_index] =
                        integer_row_polynomial[
                            static_cast<unsigned int>(
                                column_index)];
                }
            }
        }
        for (const Integer& value : integer_row) {
            divisor = gcd(
                divisor,
                value);
        }
        integers.push_back(
            std::move(integer_row));
    }
    if (divisor != 0) {
        for (auto& row : integers) {
            for (Integer& coefficient : row) {
                coefficient =
                    exact_divide(
                        coefficient,
                        divisor);
            }
        }
    }
    return integers;
}

std::vector<Integer> primitive_rational_vector(
    const std::vector<Rational>& values)
{
    return primitive_coefficients(
        polynomial_from_rationals(values));
}

std::vector<std::vector<std::string>>
coefficient_rows_text(const IntegerRows& rows)
{
    std::vector<std::vector<std::string>> result;
    result.reserve(rows.size());
    for (const auto& row : rows) {
        result.push_back(coefficient_text(row));
    }
    return result;
}

Rational compose_rational(
    Integer numerator,
    Integer denominator)
{
    if (denominator < 0) {
        numerator = -numerator;
        denominator = -denominator;
    }
    if (denominator <= 0) {
        throw AlgebraicBackendError(
            "rational candidate requires a positive denominator");
    }
    return typename CGAL::Fraction_traits<
        Rational>::Compose()(
        numerator,
        denominator);
}

std::vector<Rational> bounded_rational_candidates(
    const Polynomial& factor)
{
    const int degree = CGAL::degree(factor);
    if (degree == 1) {
        return {
            compose_rational(
                -factor[0],
                factor[1]),
        };
    }
    if (degree != 2) {
        return {};
    }
    const Integer discriminant =
        factor[1] * factor[1]
        - Integer(4) * factor[2] * factor[0];
    if (discriminant < 0) {
        return {};
    }
    Integer square_root;
    if (!CGAL::is_square(
            discriminant,
            square_root)) {
        return {};
    }
    std::vector<Rational> candidates = {
        compose_rational(
            -factor[1] - square_root,
            Integer(2) * factor[2]),
    };
    if (square_root != 0) {
        candidates.push_back(
            compose_rational(
                -factor[1] + square_root,
                Integer(2) * factor[2]));
    }
    return candidates;
}

std::pair<Rational, Rational> strict_root_interval(
    const std::vector<AlgebraicReal>& roots,
    std::size_t ordinal,
    const Polynomial& factor)
{
    Kernel kernel;
    auto interval =
        kernel.isolate_1_object()(
            roots[ordinal],
            factor);
    if (interval.first == interval.second) {
        interval.first = ordinal == 0
            ? interval.first - Rational(1)
            : kernel.bound_between_1_object()(
                  roots[ordinal - 1],
                  roots[ordinal]);
        interval.second = ordinal + 1 == roots.size()
            ? interval.second + Rational(1)
            : kernel.bound_between_1_object()(
                  roots[ordinal],
                  roots[ordinal + 1]);
    }
    if (kernel.compare_1_object()(
            interval.first,
            roots[ordinal])
            != CGAL::SMALLER
        || kernel.compare_1_object()(
            roots[ordinal],
            interval.second)
            != CGAL::SMALLER) {
        throw AlgebraicRootIsolationError(
            "trim root interval is not strict");
    }
    for (std::size_t other = 0;
         other < roots.size();
         ++other) {
        if (other == ordinal) {
            continue;
        }
        if (kernel.compare_1_object()(
                roots[other],
                interval.first)
                == CGAL::LARGER
            && kernel.compare_1_object()(
                roots[other],
                interval.second)
                == CGAL::SMALLER) {
            throw AlgebraicRootIsolationError(
                "trim root interval contains another factor root");
        }
    }
    return interval;
}

std::vector<TrimRootCandidate> exact_rim_roots(
    const Polynomial& equation)
{
    Kernel kernel;
    std::vector<
        std::pair<
            Polynomial,
            Kernel::Multiplicity_type>>
        raw_factors;
    kernel.square_free_factorize_1_object()(
        equation,
        std::back_inserter(raw_factors));

    std::vector<TrimRootCandidate> result;
    for (const auto& [raw_factor, multiplicity] :
         raw_factors) {
        const std::vector<Integer> primitive =
            primitive_coefficients(raw_factor);
        const Polynomial factor =
            polynomial_from_integers(primitive);
        std::vector<AlgebraicReal> roots;
        kernel.solve_1_object()(
            factor,
            true,
            std::back_inserter(roots));
        const std::vector<Rational> rational_candidates =
            bounded_rational_candidates(factor);
        for (std::size_t ordinal = 0;
             ordinal < roots.size();
             ++ordinal) {
            const AlgebraicReal& root = roots[ordinal];
            if (kernel.compare_1_object()(
                    root,
                    Rational(-1))
                    == CGAL::SMALLER
                || kernel.compare_1_object()(
                    root,
                    Rational(1))
                    == CGAL::LARGER) {
                continue;
            }
            std::optional<Rational> rational_value;
            for (const Rational& candidate :
                 rational_candidates) {
                if (kernel.compare_1_object()(
                        root,
                        candidate)
                    == CGAL::EQUAL) {
                    rational_value = candidate;
                    break;
                }
            }
            const auto interval =
                strict_root_interval(
                    roots,
                    ordinal,
                    factor);
            result.push_back(
                {
                    root,
                    {
                        algebraic_root_id_v1(
                            primitive,
                            ordinal),
                        coefficient_text(primitive),
                        ordinal,
                        static_cast<unsigned int>(
                            multiplicity),
                        rational_text(interval.first),
                        rational_text(interval.second),
                    },
                    rational_value,
                });
        }
    }
    const auto compare = kernel.compare_1_object();
    std::sort(
        result.begin(),
        result.end(),
        [&compare](
            const TrimRootCandidate& left,
            const TrimRootCandidate& right) {
            return compare(left.value, right.value)
                == CGAL::LARGER;
        });
    return result;
}

CGAL::Sign predicate_endpoint_sign(
    const IntegerRows& rows,
    bool at_motion_end,
    const AlgebraicReal& root)
{
    std::vector<Integer> coefficients =
        rows.front();
    if (at_motion_end) {
        for (std::size_t index = 0;
             index < rows[1].size();
             ++index) {
            coefficients[index] += rows[1][index];
        }
    }
    return Kernel().sign_at_1_object()(
        polynomial_from_integers(coefficients),
        root);
}

bool predicate_has_nonnegative_endpoint(
    const IntegerRows& rows,
    const AlgebraicReal& root)
{
    return predicate_endpoint_sign(
               rows,
               false,
               root)
            != CGAL::NEGATIVE
        || predicate_endpoint_sign(
               rows,
               true,
               root)
            != CGAL::NEGATIVE;
}

CGAL::Sign predicate_row_sum_sign(
    const IntegerRows& lower,
    const IntegerRows& upper,
    std::size_t row,
    const AlgebraicReal& root)
{
    std::vector<Integer> coefficients =
        lower[row];
    for (std::size_t index = 0;
         index < upper[row].size();
         ++index) {
        coefficients[index] += upper[row][index];
    }
    return Kernel().sign_at_1_object()(
        polynomial_from_integers(coefficients),
        root);
}

std::pair<Rational, Rational> clipped_unit_interval(
    const Rational& lambda_zero,
    const Rational& lambda_velocity)
{
    if (lambda_velocity == 0) {
        if (lambda_zero < 0 || lambda_zero > 1) {
            throw TrimFilterError(
                "branch never enters the closed trim");
        }
        return {Rational(0), Rational(1)};
    }
    Rational low = -lambda_zero / lambda_velocity;
    Rational high =
        (Rational(1) - lambda_zero)
        / lambda_velocity;
    if (low > high) {
        std::swap(low, high);
    }
    low = std::max(low, Rational(0));
    high = std::min(high, Rational(1));
    if (low > high) {
        throw TrimFilterError(
            "branch misses the closed trim domain");
    }
    return {low, high};
}

} // namespace

std::vector<TrimmedLineBranch2> solve_trimmed_line_branches(
    const std::vector<std::string>& line_support,
    const std::vector<std::string>& trim_start,
    const std::vector<std::string>& trim_end,
    const std::vector<std::string>& segment_motion,
    const std::string& cutter_radius,
    const std::string& rim_chart)
{
    const std::vector<Rational> line =
        parse_values(line_support, 3, "line support");
    const std::vector<Rational> start =
        parse_values(trim_start, 2, "trim start");
    const std::vector<Rational> end =
        parse_values(trim_end, 2, "trim end");
    const std::vector<Rational> motion =
        parse_values(
            segment_motion,
            4,
            "segment motion");
    const Rational radius =
        parse_rational(cutter_radius, "cutter radius");
    if (line[0] == 0 && line[1] == 0) {
        throw DegenerateBoundarySupportError(
            "line support requires a nonzero normal");
    }
    if (radius <= 0) {
        throw NonPositiveToolRadiusError(
            "cutter radius must be positive");
    }
    const Rational trim_start_evaluation =
        line[0] * start[0]
        + line[1] * start[1] + line[2];
    const Rational trim_end_evaluation =
        line[0] * end[0]
        + line[1] * end[1] + line[2];
    if (trim_start_evaluation != 0
        || trim_end_evaluation != 0) {
        throw TrimEndpointOffSupportError(
            "each trim endpoint must lie on the exact line support");
    }
    if (motion[0] == motion[2]
        && motion[1] == motion[3]) {
        throw ZeroLengthSegmentMotionError(
            "segment motion endpoints must be distinct");
    }
    const Rational trim_dx = end[0] - start[0];
    const Rational trim_dy = end[1] - start[1];
    if (trim_dx == 0 && trim_dy == 0) {
        throw DegenerateTrimError(
            "trim endpoints must be distinct");
    }
    const Rational start_evaluation =
        line[0] * motion[0]
        + line[1] * motion[1] + line[2];
    const Rational end_evaluation =
        line[0] * motion[2]
        + line[1] * motion[3] + line[2];
    if (start_evaluation != end_evaluation) {
        throw UnsupportedLineMotionError(
            "native line branch solver requires motion parallel to support");
    }
    const Rational velocity_x =
        motion[2] - motion[0];
    const Rational velocity_y =
        motion[3] - motion[1];
    if (velocity_x * trim_dy
        != velocity_y * trim_dx) {
        throw EventPartitionVerificationError(
            "on-support line motion is not parallel to its trim");
    }
    Rational orientation = 1;
    if (rim_chart == "rim-half-1-v1") {
        orientation = -1;
    } else if (rim_chart != "rim-half-0-v1") {
        throw ChartCoverageError(
            "trim branch requires a frozen rim chart");
    }
    const QuadraticCoefficients equation = {
        start_evaluation
            + radius * orientation * line[0],
        radius * orientation * Rational(2) * line[1],
        start_evaluation
            - radius * orientation * line[0],
    };
    const Polynomial raw_equation =
        polynomial_from_rationals(equation);
    if (CGAL::is_zero(raw_equation)) {
        throw EventPartitionVerificationError(
            "valid line-contact inputs produced an identically zero support equation");
    }
    if (CGAL::degree(raw_equation) == 0) {
        return {};
    }
    const Polynomial exact_equation =
        polynomial_from_integers(
            primitive_coefficients(
                raw_equation));
    const std::vector<TrimRootCandidate> roots =
        exact_rim_roots(exact_equation);

    const QuadraticCoefficients denominator = {
        Rational(1),
        Rational(0),
        Rational(1),
    };
    const QuadraticCoefficients rim_x = {
        orientation,
        Rational(0),
        -orientation,
    };
    const QuadraticCoefficients rim_y = {
        Rational(0),
        Rational(2) * orientation,
        Rational(0),
    };
    const QuadraticCoefficients offset_x = {
        motion[0] - start[0]
            + radius * rim_x[0],
        radius * rim_x[1],
        motion[0] - start[0]
            + radius * rim_x[2],
    };
    const QuadraticCoefficients offset_y = {
        motion[1] - start[1]
            + radius * rim_y[0],
        radius * rim_y[1],
        motion[1] - start[1]
            + radius * rim_y[2],
    };
    QuadraticCoefficients collinearity{};
    for (std::size_t index = 0;
         index < collinearity.size();
         ++index) {
        collinearity[index] =
            offset_x[index] * trim_dy
            - offset_y[index] * trim_dx;
    }

    const bool use_x = trim_dx != 0;
    const Rational trim_component =
        use_x ? trim_dx : trim_dy;
    const Rational orientation_sign =
        trim_component > 0
        ? Rational(1)
        : Rational(-1);
    const Rational motion_start_component =
        use_x ? motion[0] : motion[1];
    const Rational motion_end_component =
        use_x ? motion[2] : motion[3];
    const Rational trim_start_component =
        use_x ? start[0] : start[1];
    const Rational trim_end_component =
        use_x ? end[0] : end[1];
    const Rational velocity_component =
        motion_end_component
        - motion_start_component;
    const QuadraticCoefficients& rim_component =
        use_x ? rim_x : rim_y;
    RationalRows lower_rows(
        2,
        std::vector<Rational>(3));
    RationalRows upper_rows(
        2,
        std::vector<Rational>(3));
    for (std::size_t index = 0;
         index < denominator.size();
         ++index) {
        lower_rows[0][index] =
            orientation_sign
            * ((motion_start_component
                    - trim_start_component)
                   * denominator[index]
               + radius * rim_component[index]);
        lower_rows[1][index] =
            orientation_sign
            * velocity_component
            * denominator[index];
        upper_rows[0][index] =
            orientation_sign
            * ((trim_end_component
                    - motion_start_component)
                   * denominator[index]
               - radius * rim_component[index]);
        upper_rows[1][index] =
            -orientation_sign
            * velocity_component
            * denominator[index];
    }
    RationalRows predicate_rows = lower_rows;
    predicate_rows.insert(
        predicate_rows.end(),
        upper_rows.begin(),
        upper_rows.end());
    const IntegerRows signed_predicates =
        signed_primitive_rows(
            predicate_rows);
    const IntegerRows lower_predicate = {
        signed_predicates[0],
        signed_predicates[1],
    };
    const IntegerRows upper_predicate = {
        signed_predicates[2],
        signed_predicates[3],
    };
    const std::string local_support_id =
        stable_record(
            "trimmed-line-feature-v1",
            coefficient_text(
                primitive_rational_vector(line)));
    const std::string trim_id = stable_record(
        "closed-line-trim-v1",
        {
            rational_text(start[0]),
            rational_text(start[1]),
            rational_text(end[0]),
            rational_text(end[1]),
        });
    const std::string local_trimmed_feature_id =
        stable_record(
            "trimmed-line-producer-feature-v2",
            {
                local_support_id,
                trim_id,
                "left",
            });
    const std::string motion_id = stable_record(
        "segment-motion-v1",
        {
            rational_text(motion[0]),
            rational_text(motion[1]),
            rational_text(motion[2]),
            rational_text(motion[3]),
        });
    const std::string radius_id = stable_record(
        "cutter-radius-v1",
        {rational_text(radius)});

    Kernel kernel;
    std::vector<TrimmedLineBranch2> branches;
    for (const TrimRootCandidate& root : roots) {
        if (kernel.sign_at_1_object()(
                polynomial_from_rationals(
                    collinearity),
                root.value)
            != CGAL::ZERO) {
            throw EventPartitionVerificationError(
                "rim root does not lie on the exact trim support");
        }
        if (predicate_row_sum_sign(
                lower_predicate,
                upper_predicate,
                1,
                root.value)
            != CGAL::ZERO) {
            throw EventPartitionVerificationError(
                "trim predicate motion coefficients are not complementary");
        }
        if (predicate_row_sum_sign(
                lower_predicate,
                upper_predicate,
                0,
                root.value)
            != CGAL::POSITIVE) {
            throw EventPartitionVerificationError(
                "trim predicate width is not positive");
        }
        if (!predicate_has_nonnegative_endpoint(
                lower_predicate,
                root.value)
            || !predicate_has_nonnegative_endpoint(
                upper_predicate,
                root.value)) {
            continue;
        }

        std::optional<RationalTrimInterval2>
            rational_convenience;
        if (root.rational_value.has_value()) {
            const Rational rational_root =
                *root.rational_value;
            const Rational chart_denominator =
                Rational(1)
                + rational_root * rational_root;
            const Rational rational_rim_x =
                orientation
                * (Rational(1)
                   - rational_root * rational_root)
                / chart_denominator;
            const Rational rational_rim_y =
                orientation
                * Rational(2) * rational_root
                / chart_denominator;
            const Rational point_component =
                motion_start_component
                + radius
                    * (use_x
                        ? rational_rim_x
                        : rational_rim_y);
            const Rational lambda_zero =
                (point_component
                 - trim_start_component)
                / trim_component;
            const Rational lambda_velocity =
                velocity_component
                / trim_component;
            const auto [domain_low, domain_high] =
                clipped_unit_interval(
                    lambda_zero,
                    lambda_velocity);
            rational_convenience =
                RationalTrimInterval2{
                    rational_text(rational_root),
                    rational_text(domain_low),
                    rational_text(domain_high),
                };
        }
        const RationalTrimInterval2 legacy_rational =
            rational_convenience.value_or(
                RationalTrimInterval2{});
        branches.push_back(
            {
                legacy_rational.rim_parameter,
                legacy_rational.motion_domain_low,
                legacy_rational.motion_domain_high,
                rational_convenience,
                root.record,
                coefficient_rows_text(
                    lower_predicate),
                coefficient_rows_text(
                    upper_predicate),
                rational_convenience.has_value(),
                true,
                true,
                "accepted",
                true,
                local_support_id,
                local_support_id,
                local_trimmed_feature_id,
                trim_id,
                stable_record(
                    "trimmed-line-branch-v3",
                    {
                        local_trimmed_feature_id,
                        motion_id,
                        radius_id,
                        rim_chart,
                        root.record.root_id,
                    }),
            });
    }
    return branches;
}
