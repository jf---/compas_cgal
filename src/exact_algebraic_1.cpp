#include "exact_algebraic_1.h"

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <optional>
#include <sstream>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>

#include <boost/multiprecision/cpp_int/import_export.hpp>
#include <boost/multiprecision/integer.hpp>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/version.h>

#ifndef CGAL_USE_CORE
#error "exact_algebraic_1 requires CGAL_USE_CORE=1"
#endif

#ifndef CGAL_CORE_USE_BOOST_BACKEND
#error "exact_algebraic_1 requires CGAL_CORE_USE_BOOST_BACKEND=1"
#endif

#ifndef CGAL_USE_BOOST_MP
#error "exact_algebraic_1 requires CGAL_USE_BOOST_MP=1"
#endif

#ifndef CGAL_DISABLE_GMP
#error "exact_algebraic_1 requires CGAL_DISABLE_GMP=1"
#endif

#ifdef CGAL_CORE_USE_GMP_BACKEND
#error "exact_algebraic_1 forbids the undeclared GMP backend"
#endif

namespace {

using Integer = ExactAlgebraicInteger1;
using Kernel = ExactAlgebraicKernel1;
using Polynomial = Kernel::Polynomial_1;
using AlgebraicReal = Kernel::Algebraic_real_1;
using Rational = Kernel::Bound;

constexpr int MAX_ALGEBRAIC_DEGREE = 4;

struct RootCandidate {
    AlgebraicReal value;
    AlgebraicRootRecord2 record;
    std::vector<PartitionEvent2> events;
};

std::string u64_record(std::size_t value)
{
    std::string result;
    for (int shift = 56; shift >= 0; shift -= 8) {
        result.push_back(
            static_cast<char>(
                (static_cast<std::uint64_t>(value) >> shift)
                & 0xffU));
    }
    return result;
}

std::string ccan_node(char kind, const std::string& payload)
{
    std::string result("CCAN\0\1", 6);
    result.push_back(kind);
    result += u64_record(payload.size());
    result += payload;
    return result;
}

std::string ccan_bytes(const std::string& value)
{
    return ccan_node('B', value);
}

std::string ccan_integer(const Integer& value)
{
    const Integer magnitude = value < 0 ? -value : value;
    std::vector<unsigned char> bytes;
    if (magnitude != 0) {
        export_bits(
            magnitude,
            std::back_inserter(bytes),
            8,
            true);
    }
    std::string payload(
        1,
        value < 0
            ? static_cast<char>(1)
            : static_cast<char>(0));
    payload.append(
        reinterpret_cast<const char*>(bytes.data()),
        bytes.size());
    return ccan_node('I', payload);
}

std::string ccan_sequence(
    const std::vector<std::string>& values)
{
    std::string payload = u64_record(values.size());
    for (const std::string& value : values) {
        payload += ccan_bytes(value);
    }
    return ccan_node('S', payload);
}

std::string ccan_map(
    const std::vector<std::pair<std::string, std::string>>& fields)
{
    std::vector<std::pair<std::string, std::string>> encoded;
    encoded.reserve(fields.size());
    for (const auto& [key, value] : fields) {
        encoded.emplace_back(
            ccan_bytes(key),
            ccan_bytes(value));
    }
    std::sort(encoded.begin(), encoded.end());
    std::string payload = u64_record(encoded.size());
    for (const auto& [key, value] : encoded) {
        payload += key;
        payload += value;
    }
    return ccan_node('M', payload);
}

std::string ccan_tagged(
    const std::string& tag,
    const std::string& payload)
{
    return ccan_node(
        'T',
        ccan_bytes(tag) + ccan_bytes(payload));
}

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

std::string integer_text(const Integer& value)
{
    return value.convert_to<std::string>();
}

std::string rational_text(const Rational& value)
{
    std::ostringstream stream;
    stream << value;
    return stream.str();
}

std::pair<std::string, std::string> rational_parts(
    const Rational& value)
{
    return {
        integer_text(boost::multiprecision::numerator(value)),
        integer_text(boost::multiprecision::denominator(value)),
    };
}

Polynomial polynomial_from_strings(
    const std::vector<std::string>& coefficients)
{
    if (coefficients.empty()) {
        throw InvalidAlgebraicPolynomialError(
            "algebraic polynomial requires coefficients");
    }
    std::vector<Integer> parsed;
    parsed.reserve(coefficients.size());
    try {
        for (const std::string& coefficient : coefficients) {
            parsed.emplace_back(coefficient);
        }
    } catch (const std::exception& error) {
        throw InvalidAlgebraicPolynomialError(
            std::string("invalid integer coefficient: ")
            + error.what());
    }
    while (parsed.size() > 1 && parsed.back() == 0) {
        parsed.pop_back();
    }
    return typename CGAL::Polynomial_traits_d<
        Polynomial>::Construct_polynomial()(
            parsed.begin(),
            parsed.end());
}

std::vector<Integer> primitive_coefficients(
    const Polynomial& polynomial)
{
    if (CGAL::is_zero(polynomial)) {
        throw InvalidAlgebraicPolynomialError(
            "zero polynomial has no algebraic root identity");
    }
    const int degree = CGAL::degree(polynomial);
    if (degree < 1) {
        throw InvalidAlgebraicPolynomialError(
            "constant polynomial has no algebraic root identity");
    }
    std::vector<Integer> result;
    result.reserve(static_cast<std::size_t>(degree + 1));
    Integer common_factor = 0;
    for (int index = 0; index <= degree; ++index) {
        const Integer coefficient = polynomial[index];
        result.push_back(coefficient);
        common_factor = greatest_common_divisor(
            common_factor,
            coefficient);
    }
    if (common_factor == 0) {
        throw InvalidAlgebraicPolynomialError(
            "zero polynomial has no primitive normalization");
    }
    for (Integer& coefficient : result) {
        coefficient /= common_factor;
    }
    if (result.back() < 0) {
        for (Integer& coefficient : result) {
            coefficient = -coefficient;
        }
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

std::vector<Integer> positive_divisors(Integer value)
{
    value = value < 0 ? -value : value;
    if (value == 0) {
        throw InvalidAlgebraicPolynomialError(
            "zero has no finite divisor set");
    }
    std::vector<Integer> result;
    for (Integer divisor = 1;
         divisor * divisor <= value;
         ++divisor) {
        if (value % divisor != 0) {
            continue;
        }
        result.push_back(divisor);
        const Integer complement = value / divisor;
        if (complement != divisor) {
            result.push_back(complement);
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}

Integer evaluate_rational_root_numerator(
    const std::vector<Integer>& coefficients,
    const Integer& numerator,
    const Integer& denominator)
{
    const std::size_t degree = coefficients.size() - 1;
    Integer result = 0;
    Integer numerator_power = 1;
    Integer denominator_power = 1;
    for (std::size_t index = 0; index < degree; ++index) {
        denominator_power *= denominator;
    }
    for (std::size_t index = 0;
         index <= degree;
         ++index) {
        result += coefficients[index]
            * numerator_power * denominator_power;
        numerator_power *= numerator;
        if (index < degree && denominator != 0) {
            denominator_power /= denominator;
        }
    }
    return result;
}

std::optional<std::vector<Integer>> rational_linear_factor(
    const std::vector<Integer>& coefficients)
{
    if (coefficients.front() == 0) {
        return std::vector<Integer>{0, 1};
    }
    const std::vector<Integer> numerators =
        positive_divisors(coefficients.front());
    const std::vector<Integer> denominators =
        positive_divisors(coefficients.back());
    for (const Integer& numerator_magnitude : numerators) {
        for (const Integer& denominator : denominators) {
            if (greatest_common_divisor(
                    numerator_magnitude,
                    denominator)
                != 1) {
                continue;
            }
            for (const int sign : {-1, 1}) {
                const Integer numerator =
                    Integer(sign) * numerator_magnitude;
                if (evaluate_rational_root_numerator(
                        coefficients,
                        numerator,
                        denominator)
                    == 0) {
                    return std::vector<Integer>{
                        -numerator,
                        denominator,
                    };
                }
            }
        }
    }
    return std::nullopt;
}

std::optional<std::vector<Integer>> divide_exact(
    const std::vector<Integer>& dividend,
    const std::vector<Integer>& divisor)
{
    if (divisor.empty() || divisor.back() == 0
        || dividend.size() < divisor.size()) {
        return std::nullopt;
    }
    std::vector<Integer> remainder = dividend;
    std::vector<Integer> quotient(
        dividend.size() - divisor.size() + 1,
        Integer(0));
    for (std::size_t offset = quotient.size();
         offset-- > 0;) {
        const std::size_t leading_index =
            offset + divisor.size() - 1;
        if (remainder[leading_index]
                % divisor.back()
            != 0) {
            return std::nullopt;
        }
        const Integer coefficient =
            remainder[leading_index]
            / divisor.back();
        quotient[offset] = coefficient;
        for (std::size_t index = 0;
             index < divisor.size();
             ++index) {
            remainder[offset + index] -=
                coefficient * divisor[index];
        }
    }
    if (std::any_of(
            remainder.begin(),
            remainder.end(),
            [](const Integer& value) {
                return value != 0;
            })) {
        return std::nullopt;
    }
    return primitive_coefficients(
        polynomial_from_integers(quotient));
}

Integer integer_square_root(const Integer& value)
{
    if (value < 0) {
        throw InvalidAlgebraicPolynomialError(
            "negative integer has no real square root");
    }
    if (value < 2) {
        return value;
    }
    const unsigned int bit_count =
        boost::multiprecision::msb(value) + 1;
    Integer estimate =
        Integer(1) << ((bit_count + 1) / 2);
    while (true) {
        const Integer next =
            (estimate + value / estimate) / 2;
        if (next >= estimate) {
            return estimate;
        }
        estimate = next;
    }
}

std::vector<Integer> quadratic_integer_solutions(
    const Integer& quadratic,
    const Integer& linear,
    const Integer& constant)
{
    if (quadratic == 0) {
        if (linear == 0 || (-constant) % linear != 0) {
            return {};
        }
        return {-constant / linear};
    }
    const Integer discriminant =
        linear * linear
        - Integer(4) * quadratic * constant;
    if (discriminant < 0) {
        return {};
    }
    const Integer square_root =
        integer_square_root(discriminant);
    if (square_root * square_root != discriminant) {
        return {};
    }
    std::vector<Integer> result;
    const Integer denominator = Integer(2) * quadratic;
    for (const int sign : {-1, 1}) {
        const Integer signed_root =
            Integer(sign) * square_root;
        const Integer numerator = -linear + signed_root;
        if (numerator % denominator == 0) {
            result.push_back(numerator / denominator);
        }
    }
    std::sort(result.begin(), result.end());
    result.erase(
        std::unique(result.begin(), result.end()),
        result.end());
    return result;
}

std::optional<std::vector<Integer>>
quadratic_factor_of_quartic(
    const std::vector<Integer>& coefficients)
{
    const Integer& constant = coefficients[0];
    const Integer& linear = coefficients[1];
    const Integer& quadratic = coefficients[2];
    const Integer& cubic = coefficients[3];
    const Integer& quartic = coefficients[4];
    const std::vector<Integer> leading_divisors =
        positive_divisors(quartic);
    const std::vector<Integer> constant_divisors =
        positive_divisors(constant);
    for (const Integer& leading_magnitude :
         leading_divisors) {
        for (const int leading_sign : {-1, 1}) {
            const Integer factor_leading =
                Integer(leading_sign)
                * leading_magnitude;
            const Integer quotient_leading =
                quartic / factor_leading;
            for (const Integer& constant_magnitude :
                 constant_divisors) {
                for (const int constant_sign :
                     {-1, 1}) {
                    const Integer factor_constant =
                        Integer(constant_sign)
                        * constant_magnitude;
                    const Integer quotient_constant =
                        constant / factor_constant;
                    const Integer determinant =
                        quotient_leading
                            * factor_constant
                        - factor_leading
                            * quotient_constant;
                    std::vector<Integer>
                        factor_linear_candidates;
                    if (determinant != 0) {
                        const Integer numerator =
                            cubic * factor_constant
                            - factor_leading * linear;
                        if (numerator % determinant != 0) {
                            continue;
                        }
                        factor_linear_candidates.push_back(
                            numerator / determinant);
                    } else {
                        factor_linear_candidates =
                            quadratic_integer_solutions(
                                -quotient_leading,
                                cubic,
                                factor_leading
                                    * (factor_leading
                                           * quotient_constant
                                       + factor_constant
                                           * quotient_leading
                                       - quadratic));
                    }
                    for (const Integer& factor_linear :
                         factor_linear_candidates) {
                        const Integer remainder =
                            cubic
                            - quotient_leading
                                * factor_linear;
                        if (remainder % factor_leading
                            != 0) {
                            continue;
                        }
                        const Integer quotient_linear =
                            remainder / factor_leading;
                        const std::vector<Integer> factor = {
                            factor_constant,
                            factor_linear,
                            factor_leading,
                        };
                        const std::vector<Integer> quotient = {
                            quotient_constant,
                            quotient_linear,
                            quotient_leading,
                        };
                        if (factor_linear
                                    * quotient_constant
                                + factor_constant
                                    * quotient_linear
                                != linear
                            || factor_linear
                                        * quotient_linear
                                    + factor_leading
                                        * quotient_constant
                                    + factor_constant
                                        * quotient_leading
                                != quadratic) {
                            continue;
                        }
                        const std::vector<Integer> primitive =
                            primitive_coefficients(
                                polynomial_from_integers(
                                    factor));
                        if (divide_exact(
                                coefficients,
                                primitive)
                            .has_value()) {
                            return primitive;
                        }
                    }
                }
            }
        }
    }
    return std::nullopt;
}

std::vector<std::vector<Integer>>
irreducible_factors_degree_four(
    const std::vector<Integer>& coefficients)
{
    const std::vector<Integer> primitive =
        primitive_coefficients(
            polynomial_from_integers(coefficients));
    const int degree =
        static_cast<int>(primitive.size()) - 1;
    if (degree > MAX_ALGEBRAIC_DEGREE) {
        throw UnsupportedAlgebraicDegreeError(
            "exact factorization supports degree 4 or lower");
    }
    if (degree <= 1) {
        return {primitive};
    }
    if (const auto linear =
            rational_linear_factor(primitive)) {
        const auto quotient =
            divide_exact(primitive, *linear);
        if (!quotient) {
            throw InvalidAlgebraicPolynomialError(
                "proven rational root did not divide polynomial");
        }
        std::vector<std::vector<Integer>> result = {
            primitive_coefficients(
                polynomial_from_integers(*linear)),
        };
        std::vector<std::vector<Integer>> remainder =
            irreducible_factors_degree_four(*quotient);
        result.insert(
            result.end(),
            remainder.begin(),
            remainder.end());
        std::sort(result.begin(), result.end());
        return result;
    }
    if (degree < MAX_ALGEBRAIC_DEGREE) {
        return {primitive};
    }
    if (const auto quadratic =
            quadratic_factor_of_quartic(primitive)) {
        const auto quotient =
            divide_exact(primitive, *quadratic);
        if (!quotient) {
            throw InvalidAlgebraicPolynomialError(
                "proven quadratic factor did not divide polynomial");
        }
        std::vector<std::vector<Integer>> result =
            irreducible_factors_degree_four(*quadratic);
        std::vector<std::vector<Integer>> remainder =
            irreducible_factors_degree_four(*quotient);
        result.insert(
            result.end(),
            remainder.begin(),
            remainder.end());
        std::sort(result.begin(), result.end());
        return result;
    }
    return {primitive};
}

std::vector<std::string> coefficient_text(
    const std::vector<Integer>& coefficients)
{
    std::vector<std::string> result;
    result.reserve(coefficients.size());
    for (const Integer& coefficient : coefficients) {
        result.push_back(integer_text(coefficient));
    }
    return result;
}

std::pair<Rational, Rational> strict_interval(
    const AlgebraicReal& root,
    const Polynomial& factor)
{
    Kernel kernel;
    auto interval =
        kernel.isolate_1_object()(root, factor);
    if (interval.first == interval.second) {
        const Rational offset(Integer(1), Integer(16));
        interval.first -= offset;
        interval.second += offset;
    }
    if (kernel.compare_1_object()(interval.first, root)
            != CGAL::SMALLER
        || kernel.compare_1_object()(root, interval.second)
            != CGAL::SMALLER) {
        throw AlgebraicRootIsolationError(
            "isolating interval is not strict");
    }
    return interval;
}

Rational rational_between(
    const AlgebraicReal& left,
    const AlgebraicReal& right)
{
    Kernel kernel;
    if (left.low() == left.high()
        && right.low() == right.high()) {
        return (left.low() + right.low()) / Rational(2);
    }
    return kernel.bound_between_1_object()(left, right);
}

bool event_less(
    const PartitionEvent2& left,
    const PartitionEvent2& right)
{
    return std::tie(
               left.feature_id,
               left.support_id,
               left.trim_id,
               left.vertex_id,
               left.branch_id,
               left.kind,
               left.disposition)
        < std::tie(
               right.feature_id,
               right.support_id,
               right.trim_id,
               right.vertex_id,
               right.branch_id,
               right.kind,
               right.disposition);
}

} // namespace

AlgebraicBackendEvidence2 exact_algebraic_backend_evidence()
{
    static_assert(
        std::is_same_v<
            ExactAlgebraicInteger1,
            boost::multiprecision::cpp_int>);
    return {
        CGAL_VERSION_STR,
        "CORE::BigInt(boost::multiprecision::cpp_int)",
        "CGAL::Algebraic_kernel_d_1<CORE::BigInt>",
        "CGAL::Algebraic_kernel_d_2<CORE::BigInt>",
        "CGAL::Arr_algebraic_segment_traits_2<CORE::BigInt>",
        {
            "CGAL_DISABLE_GMP=1",
            "CGAL_USE_BOOST_MP=1",
            "CGAL_USE_CORE=1",
            "CGAL_CORE_USE_BOOST_BACKEND=1",
        },
    };
}

std::string algebraic_root_id_v1(
    const std::vector<ExactAlgebraicInteger1>& primitive_factor,
    std::size_t root_ordinal)
{
    if (primitive_factor.size() < 2) {
        throw InvalidAlgebraicPolynomialError(
            "algebraic root identity requires a nonconstant factor");
    }
    std::vector<std::string> encoded_coefficients;
    encoded_coefficients.reserve(primitive_factor.size());
    for (const Integer& coefficient : primitive_factor) {
        encoded_coefficients.push_back(
            ccan_integer(coefficient));
    }
    return ccan_tagged(
        "algebraic-root-id-v1",
        ccan_map(
            {
                {
                    "coefficients",
                    ccan_sequence(encoded_coefficients),
                },
                {
                    "root-ordinal",
                    ccan_integer(Integer(root_ordinal)),
                },
            }));
}

EventPartitionCertificate2 partition_integer_projections(
    const std::vector<ProjectionInput2>& projections)
{
    if (projections.empty()) {
        throw InvalidAlgebraicPolynomialError(
            "event partition requires at least one projection");
    }

    Kernel kernel;
    std::vector<RootCandidate> candidates;
    std::vector<ProjectionRecord2> projection_records;
    projection_records.reserve(projections.size());
    for (const ProjectionInput2& input : projections) {
        const Polynomial polynomial =
            polynomial_from_strings(input.coefficients);
        if (CGAL::is_zero(polynomial)
            || CGAL::degree(polynomial) < 1) {
            throw InvalidAlgebraicPolynomialError(
                "projection polynomial must be nonzero and nonconstant");
        }
        if (CGAL::degree(polynomial)
            > MAX_ALGEBRAIC_DEGREE) {
            throw UnsupportedAlgebraicDegreeError(
                "exact factorization supports degree 4 or lower");
        }
        std::vector<std::pair<Polynomial, int>> factors;
        kernel.square_free_factorize_1_object()(
            polynomial,
            std::back_inserter(factors));

        ProjectionRecord2 projection_record;
        projection_record.projection_id = input.projection_id;
        projection_record.coefficient_rows = {
            coefficient_text(
                primitive_coefficients(polynomial)),
        };

        for (const auto& [raw_factor, multiplicity] : factors) {
            const auto irreducible_factors =
                irreducible_factors_degree_four(
                    primitive_coefficients(raw_factor));
            for (const std::vector<Integer>& primitive :
                 irreducible_factors) {
                const Polynomial factor =
                    polynomial_from_integers(primitive);
                const std::vector<std::string> factor_text =
                    coefficient_text(primitive);
                projection_record.factor_coefficients.push_back(
                    factor_text);

                std::vector<AlgebraicReal> roots;
                kernel.solve_1_object()(
                    factor,
                    true,
                    std::back_inserter(roots));
                for (std::size_t ordinal = 0;
                     ordinal < roots.size();
                     ++ordinal) {
                    const AlgebraicReal& root = roots[ordinal];
                    if (kernel.compare_1_object()(
                            root,
                            Rational(0))
                            == CGAL::SMALLER
                        || kernel.compare_1_object()(
                            root,
                            Rational(1))
                            == CGAL::LARGER) {
                        continue;
                    }
                    const auto interval =
                        strict_interval(root, factor);
                    candidates.push_back(
                        {
                            root,
                            {
                                algebraic_root_id_v1(
                                    primitive,
                                    ordinal),
                                factor_text,
                                ordinal,
                                static_cast<unsigned int>(
                                    multiplicity),
                                rational_text(interval.first),
                                rational_text(interval.second),
                            },
                            input.events,
                        });
                }
            }
        }
        std::sort(
            projection_record.factor_coefficients.begin(),
            projection_record.factor_coefficients.end());
        projection_records.push_back(
            std::move(projection_record));
    }

    std::sort(
        candidates.begin(),
        candidates.end(),
        [&kernel](
            const RootCandidate& left,
            const RootCandidate& right) {
            return kernel.compare_1_object()(
                       left.value,
                       right.value)
                == CGAL::SMALLER;
        });

    std::vector<RootCandidate> unique;
    for (RootCandidate& candidate : candidates) {
        if (!unique.empty()
            && kernel.compare_1_object()(
                   unique.back().value,
                   candidate.value)
                == CGAL::EQUAL) {
            if (unique.back().record.root_id
                != candidate.record.root_id) {
                throw InvalidAlgebraicPolynomialError(
                    "equal roots have different irreducible identities");
            }
            unique.back().record.multiplicity = std::max(
                unique.back().record.multiplicity,
                candidate.record.multiplicity);
            unique.back().events.insert(
                unique.back().events.end(),
                candidate.events.begin(),
                candidate.events.end());
            continue;
        }
        unique.push_back(std::move(candidate));
    }

    EventPartitionCertificate2 certificate;
    certificate.build_evidence =
        exact_algebraic_backend_evidence();
    certificate.projections =
        std::move(projection_records);
    certificate.source_kind = "integer-projections-v1";

    for (RootCandidate& root : unique) {
        std::sort(
            root.events.begin(),
            root.events.end(),
            event_less);
        root.events.erase(
            std::unique(
                root.events.begin(),
                root.events.end(),
                [](const PartitionEvent2& left,
                   const PartitionEvent2& right) {
                    return !event_less(left, right)
                        && !event_less(right, left);
                }),
            root.events.end());
        certificate.roots.push_back(root.record);
        certificate.fibres.push_back(
            {
                root.record.root_id,
                root.events,
            });
    }

    std::vector<AlgebraicReal> boundaries;
    boundaries.reserve(unique.size() + 2);
    boundaries.push_back(
        kernel.construct_algebraic_real_1_object()(
            Rational(0)));
    for (const RootCandidate& root : unique) {
        boundaries.push_back(root.value);
    }
    boundaries.push_back(
        kernel.construct_algebraic_real_1_object()(
            Rational(1)));
    for (std::size_t index = 0;
         index + 1 < boundaries.size();
         ++index) {
        if (kernel.compare_1_object()(
                boundaries[index],
                boundaries[index + 1])
            != CGAL::SMALLER) {
            continue;
        }
        const Rational witness = rational_between(
            boundaries[index],
            boundaries[index + 1]);
        const auto [numerator, denominator] =
            rational_parts(witness);
        certificate.cells.push_back(
            {
                index == 0
                    ? std::string()
                    : unique[index - 1].record.root_id,
                index == unique.size()
                    ? std::string()
                    : unique[index].record.root_id,
                numerator,
                denominator,
                {},
            });
    }
    return certificate;
}
