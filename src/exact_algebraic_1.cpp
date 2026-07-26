#include "exact_algebraic_1.h"

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <sstream>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>

#include <boost/multiprecision/cpp_int/import_export.hpp>
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
using PolynomialTraits = CGAL::Polynomial_traits_d<Polynomial>;

struct RootCandidate {
    AlgebraicReal value;
    Polynomial governing_factor;
    AlgebraicRootRecord2 record;
    struct SourceEventEvidence {
        PartitionEvent2 event;
        std::string projection_id;
        std::string factor_id;
        unsigned int multiplicity;
    };
    std::vector<SourceEventEvidence> event_evidence;
};

struct FactorRootCacheEntry {
    Polynomial factor;
    std::vector<AlgebraicReal> roots;
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
    PolynomialTraits::Construct_polynomial construct;
    return construct(parsed.begin(), parsed.end());
}

Polynomial normalized_polynomial(
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
    PolynomialTraits::Canonicalize canonicalize;
    return canonicalize(polynomial);
}

std::vector<Integer> primitive_coefficients(
    const Polynomial& polynomial)
{
    const Polynomial normalized =
        normalized_polynomial(polynomial);
    const int degree = CGAL::degree(normalized);
    std::vector<Integer> result;
    result.reserve(static_cast<std::size_t>(degree + 1));
    for (int index = 0; index <= degree; ++index) {
        result.push_back(normalized[index]);
    }
    return result;
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

Rational parse_rational_interval_bound(
    const std::string& text)
{
    const std::size_t separator = text.find('/');
    try {
        if (separator == std::string::npos) {
            return Rational(Integer(text));
        }
        if (text.find('/', separator + 1)
            != std::string::npos) {
            throw AlgebraicRootIsolationError(
                "root interval bound is not an exact rational");
        }
        const Integer numerator(
            text.substr(0, separator));
        const Integer denominator(
            text.substr(separator + 1));
        if (denominator == 0) {
            throw AlgebraicRootIsolationError(
                "root interval bound has zero denominator");
        }
        return Rational(numerator, denominator);
    } catch (const EventSubstrateError&) {
        throw;
    } catch (const std::exception&) {
        throw AlgebraicRootIsolationError(
            "root interval bound is not an exact rational");
    }
}

void require_strict_isolating_interval(
    const AlgebraicReal& root,
    const std::vector<AlgebraicReal>& factor_roots,
    std::size_t ordinal,
    const std::pair<Rational, Rational>& interval,
    Kernel& kernel)
{
    if (ordinal >= factor_roots.size()
        || kernel.compare_1_object()(
               factor_roots[ordinal],
               root)
            != CGAL::EQUAL) {
        throw AlgebraicRootIsolationError(
            "governing factor does not contain represented root");
    }
    if (kernel.compare_1_object()(interval.first, root)
            != CGAL::SMALLER
        || kernel.compare_1_object()(root, interval.second)
            != CGAL::SMALLER) {
        throw AlgebraicRootIsolationError(
            "isolating interval is not strict");
    }
    for (std::size_t other = 0;
         other < factor_roots.size();
         ++other) {
        if (other == ordinal) {
            continue;
        }
        if (kernel.compare_1_object()(
                interval.first,
                factor_roots[other])
                == CGAL::SMALLER
            && kernel.compare_1_object()(
                factor_roots[other],
                interval.second)
                == CGAL::SMALLER) {
            throw AlgebraicRootIsolationError(
                "isolating interval contains another factor root");
        }
    }
}

std::pair<Rational, Rational> strict_interval(
    const AlgebraicReal& root,
    const Polynomial& factor,
    const std::vector<AlgebraicReal>& factor_roots,
    std::size_t ordinal,
    Kernel& kernel)
{
    auto interval =
        kernel.isolate_1_object()(root, factor);
    if (interval.first == interval.second) {
        interval.first = ordinal == 0
            ? interval.first - Rational(1)
            : kernel.bound_between_1_object()(
                  factor_roots[ordinal - 1],
                  factor_roots[ordinal]);
        interval.second = ordinal + 1 == factor_roots.size()
            ? interval.second + Rational(1)
            : kernel.bound_between_1_object()(
                  factor_roots[ordinal],
                  factor_roots[ordinal + 1]);
    }
    require_strict_isolating_interval(
        root,
        factor_roots,
        ordinal,
        interval,
        kernel);
    return interval;
}

const std::vector<AlgebraicReal>& cached_factor_roots(
    const Polynomial& factor,
    std::vector<FactorRootCacheEntry>& cache,
    Kernel& kernel)
{
    const auto cached = std::find_if(
        cache.begin(),
        cache.end(),
        [&factor](const FactorRootCacheEntry& entry) {
            return entry.factor == factor;
        });
    if (cached != cache.end()) {
        return cached->roots;
    }
    FactorRootCacheEntry entry;
    entry.factor = factor;
    kernel.solve_1_object()(
        factor,
        true,
        std::back_inserter(entry.roots));
    cache.push_back(std::move(entry));
    return cache.back().roots;
}

std::size_t root_ordinal(
    const AlgebraicReal& root,
    const std::vector<AlgebraicReal>& factor_roots,
    Kernel& kernel)
{
    for (std::size_t ordinal = 0;
         ordinal < factor_roots.size();
         ++ordinal) {
        if (kernel.compare_1_object()(
                factor_roots[ordinal],
                root)
            == CGAL::EQUAL) {
            return ordinal;
        }
    }
    throw InvalidAlgebraicPolynomialError(
        "governing polynomial does not contain merged root");
}

void finalize_root_records(
    std::vector<RootCandidate>& roots)
{
    Kernel kernel;
    std::vector<Rational> separators;
    separators.reserve(
        roots.empty() ? 0 : roots.size() - 1);
    for (std::size_t index = 0;
         index + 1 < roots.size();
         ++index) {
        separators.push_back(
            kernel.bound_between_1_object()(
                roots[index].value,
                roots[index + 1].value));
    }
    std::vector<FactorRootCacheEntry> factor_root_cache;
    std::vector<std::pair<Rational, Rational>>
        certified_intervals;
    certified_intervals.reserve(roots.size());
    for (std::size_t index = 0;
         index < roots.size();
         ++index) {
        RootCandidate& candidate = roots[index];
        const std::vector<AlgebraicReal>& factor_roots =
            cached_factor_roots(
                candidate.governing_factor,
                factor_root_cache,
                kernel);
        const std::size_t ordinal = root_ordinal(
            candidate.value,
            factor_roots,
            kernel);
        auto interval = strict_interval(
            candidate.value,
            candidate.governing_factor,
            factor_roots,
            ordinal,
            kernel);
        if (index > 0) {
            interval.first =
                std::max(
                    interval.first,
                    separators[index - 1]);
        }
        if (index + 1 < roots.size()) {
            interval.second =
                std::min(
                    interval.second,
                    separators[index]);
        }
        if (kernel.compare_1_object()(
                interval.first,
                candidate.value)
                != CGAL::SMALLER
            || kernel.compare_1_object()(
                candidate.value,
                interval.second)
                != CGAL::SMALLER) {
            throw AlgebraicRootIsolationError(
                "global root interval is not strict");
        }
        const std::vector<Integer> primitive =
            primitive_coefficients(
                candidate.governing_factor);
        candidate.record = {
            algebraic_root_id_v1(primitive, ordinal),
            coefficient_text(primitive),
            ordinal,
            candidate.record.multiplicity,
            rational_text(interval.first),
            rational_text(interval.second),
        };
        certified_intervals.push_back(interval);
    }
    for (std::size_t index = 0;
         index + 1 < certified_intervals.size();
         ++index) {
        if (certified_intervals[index].second
            > certified_intervals[index + 1].first) {
            throw AlgebraicRootIsolationError(
                "adjacent global root intervals overlap");
        }
    }
}

Polynomial common_governing_factor(
    const Polynomial& left,
    const Polynomial& right)
{
    PolynomialTraits::Gcd_up_to_constant_factor gcd;
    const Polynomial common = gcd(left, right);
    if (CGAL::degree(common) < 1) {
        throw InvalidAlgebraicPolynomialError(
            "equal roots have no common governing polynomial");
    }
    return normalized_polynomial(common);
}

Rational simplest_dyadic_between(
    const AlgebraicReal& left,
    const AlgebraicReal& right)
{
    Kernel kernel;
    Integer denominator = 1;
    while (true) {
        denominator *= 2;
        const Rational scaled_upper =
            left.high() * Rational(denominator);
        Integer numerator =
            CORE::numerator(scaled_upper)
                / CORE::denominator(scaled_upper)
            + 1;
        while (numerator < denominator) {
            const Rational candidate(
                numerator,
                denominator);
            if (kernel.compare_1_object()(
                    left,
                    candidate)
                    == CGAL::SMALLER
                && kernel.compare_1_object()(
                    candidate,
                    right)
                    == CGAL::SMALLER) {
                return candidate;
            }
            if (kernel.compare_1_object()(
                    candidate,
                    right)
                != CGAL::SMALLER) {
                break;
            }
            ++numerator;
        }
    }
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
               left.disposition,
               left.endpoint_role,
               left.left_active_count,
               left.right_active_count,
               left.incidence_permutation_rechecked,
               left.original_equations_rechecked,
               left.orientation_rechecked,
               left.trim_disposition,
               left.pair_sheet_id,
               left.first_feature_id,
               left.second_feature_id,
               left.first_chart_id,
               left.second_chart_id,
               left.first_branch_id,
               left.second_branch_id)
        < std::tie(
               right.feature_id,
               right.support_id,
               right.trim_id,
               right.vertex_id,
               right.branch_id,
               right.kind,
               right.disposition,
               right.endpoint_role,
               right.left_active_count,
               right.right_active_count,
               right.incidence_permutation_rechecked,
               right.original_equations_rechecked,
               right.orientation_rechecked,
               right.trim_disposition,
               right.pair_sheet_id,
               right.first_feature_id,
               right.second_feature_id,
               right.first_chart_id,
               right.second_chart_id,
               right.first_branch_id,
               right.second_branch_id);
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

void validate_algebraic_root_intervals(
    const std::vector<AlgebraicRootRecord2>& roots)
{
    Kernel kernel;
    std::vector<FactorRootCacheEntry> factor_root_cache;
    std::vector<std::pair<Rational, Rational>> intervals;
    intervals.reserve(roots.size());
    for (const AlgebraicRootRecord2& record : roots) {
        const Polynomial factor = normalized_polynomial(
            polynomial_from_strings(
                record.factor_coefficients));
        const std::vector<Integer> primitive =
            primitive_coefficients(factor);
        if (coefficient_text(primitive)
                != record.factor_coefficients
            || algebraic_root_id_v1(
                   primitive,
                   record.root_ordinal)
                != record.root_id) {
            throw AlgebraicRootIsolationError(
                "root identity does not match governing factor");
        }
        const std::vector<AlgebraicReal>& factor_roots =
            cached_factor_roots(
                factor,
                factor_root_cache,
                kernel);
        if (record.root_ordinal
            >= factor_roots.size()) {
            throw AlgebraicRootIsolationError(
                "root ordinal exceeds governing factor roots");
        }
        const AlgebraicReal& root =
            factor_roots[record.root_ordinal];
        const std::pair<Rational, Rational> interval = {
            parse_rational_interval_bound(
                record.interval_low),
            parse_rational_interval_bound(
                record.interval_high),
        };
        require_strict_isolating_interval(
            root,
            factor_roots,
            record.root_ordinal,
            interval,
            kernel);
        intervals.push_back(interval);
    }
    for (std::size_t index = 0;
         index + 1 < intervals.size();
         ++index) {
        if (intervals[index].second
            > intervals[index + 1].first) {
            throw AlgebraicRootIsolationError(
                "adjacent certified root intervals overlap");
        }
    }
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
        Polynomial signed_predicate;
        if (!input.signed_predicate_coefficients.empty()) {
            signed_predicate = polynomial_from_strings(
                input.signed_predicate_coefficients);
            if (normalized_polynomial(signed_predicate)
                != normalized_polynomial(polynomial)) {
                throw InvalidAlgebraicPolynomialError(
                    "signed predicate does not govern the projection zero set");
            }
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
        const int actual_degree =
            CGAL::degree(polynomial);
        projection_record.actual_motion_degree =
            actual_degree;
        projection_record.actual_rim_degree = 0;
        projection_record.bound_motion_degree =
            actual_degree;
        projection_record.bound_rim_degree = 0;
        projection_record.degree_bound_id =
            "exact-univariate-("
            + std::to_string(actual_degree)
            + ",0)-v1";
        projection_record.normalized_coefficient_bytes =
            ccan_sequence(
                projection_record.coefficient_rows.front());
        projection_record.signed_predicate_coefficients =
            input.signed_predicate_coefficients;

        for (const auto& [raw_factor, multiplicity] : factors) {
            const Polynomial factor =
                normalized_polynomial(raw_factor);
            const std::vector<Integer> primitive =
                primitive_coefficients(factor);
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
                std::vector<
                    RootCandidate::SourceEventEvidence>
                    event_evidence;
                event_evidence.reserve(
                    input.events.size());
                for (const PartitionEvent2& event :
                     input.events) {
                    event_evidence.push_back(
                        {
                            event,
                            input.projection_id,
                            ccan_sequence(factor_text),
                            static_cast<unsigned int>(
                                multiplicity),
                        });
                }
                candidates.push_back(
                    {
                        root,
                        factor,
                        {
                            algebraic_root_id_v1(
                                primitive,
                                ordinal),
                            factor_text,
                            ordinal,
                            static_cast<unsigned int>(
                                multiplicity),
                            {},
                            {},
                        },
                        std::move(event_evidence),
                    });
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
            const unsigned int multiplicity = std::max(
                unique.back().record.multiplicity,
                candidate.record.multiplicity);
            unique.back().governing_factor =
                common_governing_factor(
                    unique.back().governing_factor,
                    candidate.governing_factor);
            unique.back().record.multiplicity =
                multiplicity;
            unique.back().event_evidence.insert(
                unique.back().event_evidence.end(),
                candidate.event_evidence.begin(),
                candidate.event_evidence.end());
            continue;
        }
        unique.push_back(std::move(candidate));
    }
    finalize_root_records(unique);

    EventPartitionCertificate2 certificate;
    certificate.build_evidence =
        exact_algebraic_backend_evidence();
    certificate.projections =
        std::move(projection_records);
    certificate.source_kind = "integer-projections-v1";

    for (RootCandidate& root : unique) {
        std::vector<PartitionEvent2> events;
        events.reserve(root.event_evidence.size());
        for (const auto& evidence :
             root.event_evidence) {
            events.push_back(evidence.event);
        }
        std::sort(
            events.begin(),
            events.end(),
            event_less);
        events.erase(
            std::unique(
                events.begin(),
                events.end(),
                [](const PartitionEvent2& left,
                   const PartitionEvent2& right) {
                    return !event_less(left, right)
                        && !event_less(right, left);
                }),
            events.end());
        std::vector<LocalEventWitness2>
            local_event_witnesses;
        local_event_witnesses.reserve(
            root.event_evidence.size());
        for (const auto& evidence :
             root.event_evidence) {
            local_event_witnesses.push_back(
                {
                    evidence.event.kind,
                    evidence.event.feature_id,
                    evidence.event.support_id,
                    evidence.event.trim_id,
                    evidence.event.vertex_id,
                    evidence.event.endpoint_role,
                    evidence.event.disposition,
                    evidence.event.branch_id,
                    root.record.root_id,
                    evidence.projection_id,
                    evidence.factor_id,
                    evidence.multiplicity,
                });
        }
        std::sort(
            local_event_witnesses.begin(),
            local_event_witnesses.end(),
            [](const LocalEventWitness2& left,
               const LocalEventWitness2& right) {
                return std::tie(
                           left.kind,
                           left.feature_id,
                           left.support_id,
                           left.trim_id,
                           left.vertex_id,
                           left.endpoint_role,
                           left.disposition,
                           left.local_branch_id,
                           left.source_projection_id,
                           left.source_factor_id,
                           left.multiplicity)
                    < std::tie(
                           right.kind,
                           right.feature_id,
                           right.support_id,
                           right.trim_id,
                           right.vertex_id,
                           right.endpoint_role,
                           right.disposition,
                           right.local_branch_id,
                           right.source_projection_id,
                           right.source_factor_id,
                           right.multiplicity);
            });
        certificate.roots.push_back(root.record);
        EventFibre2 fibre;
        fibre.root_id = root.record.root_id;
        fibre.events = std::move(events);
        fibre.local_event_witnesses =
            std::move(local_event_witnesses);
        certificate.fibres.push_back(
            std::move(fibre));
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
        const Rational witness = simplest_dyadic_between(
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
