#include "circle_pair_projection.h"

#include "event_certificate.h"

#include <algorithm>
#include <array>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include <CGAL/CORE/BigRat.h>
#include <CGAL/Fraction_traits.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/number_utils.h>

namespace {

using Integer = CORE::BigInt;
using Rational = CORE::BigRat;
using Polynomial1Q = CGAL::Polynomial<Rational>;
using Polynomial2Q = CGAL::Polynomial<Polynomial1Q>;
using Polynomial3Q = CGAL::Polynomial<Polynomial2Q>;
using Polynomial1Z = CGAL::Polynomial<Integer>;

constexpr std::array<const char*, 4> CENTER_CHART_IDS{
    "center-quarter-0-v1",
    "center-quarter-1-v1",
    "center-quarter-2-v1",
    "center-quarter-3-v1",
};

constexpr std::array<const char*, 2> RIM_CHART_IDS{
    "rim-half-0-v1",
    "rim-half-1-v1",
};

struct FullCircleFeatureSource {
    std::string feature_id;
    std::string support_id;
    std::string trim_id;
    std::string support_kind;
};

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
        const Integer numerator(text.substr(0, separator));
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

Polynomial3Q variable(
    int index)
{
    return CGAL::shift(Polynomial3Q(Rational(1)), 1, index);
}

Polynomial3Q pullback(
    const ProjectionRecord2& projection,
    int rim_variable)
{
    const Polynomial3Q motion = variable(0);
    const Polynomial3Q rim = variable(rim_variable);
    Polynomial3Q result(Rational(0));
    Polynomial3Q motion_power(Rational(1));
    for (const auto& row : projection.coefficient_rows) {
        Polynomial3Q rim_power(Rational(1));
        for (const std::string& coefficient : row) {
            result +=
                parse_rational(coefficient,
                               "full-circle pair pullback")
                * motion_power * rim_power;
            rim_power *= rim;
        }
        motion_power *= motion;
    }
    return result;
}

struct RimCoordinates {
    Polynomial3Q x;
    Polynomial3Q y;
    Polynomial3Q denominator;
};

RimCoordinates rim_coordinates(
    int variable_index,
    const std::string& chart)
{
    const Rational sign = chart == "rim-half-0-v1"
                              ? Rational(1)
                              : Rational(-1);
    if (chart != "rim-half-0-v1"
        && chart != "rim-half-1-v1") {
        throw ChartCoverageError(
            "full-circle pair projection uses an unknown "
            "rim chart");
    }
    const Polynomial3Q parameter = variable(variable_index);
    const Polynomial3Q square = parameter * parameter;
    return {
        sign * (Polynomial3Q(Rational(1)) - square),
        Rational(2) * sign * parameter,
        Polynomial3Q(Rational(1)) + square,
    };
}

Polynomial3Q orientation_predicate(
    const RimCoordinates& first,
    const RimCoordinates& second)
{
    return first.x * second.y - first.y * second.x;
}

Polynomial3Q cap_predicate(
    const RimCoordinates& first,
    const RimCoordinates& second,
    const Rational& cap_ratio)
{
    return (Rational(2) - cap_ratio) * first.denominator
               * second.denominator
           - Rational(2)
                 * (first.x * second.x
                    + first.y * second.y);
}

Polynomial1Z eliminate(
    const Polynomial3Q& first_pullback,
    const Polynomial3Q& second_pullback,
    const Polynomial3Q& predicate)
{
    const Polynomial2Q inner =
        typename CGAL::Polynomial_traits_d<
            Polynomial3Q>::Resultant()(second_pullback,
                                       predicate);
    const Polynomial1Q result =
        typename CGAL::Polynomial_traits_d<
            Polynomial2Q>::Resultant()(first_pullback[0],
                                       inner);
    using FractionTraits =
        CGAL::Fraction_traits<Polynomial1Q>;
    typename FractionTraits::Numerator_type numerator;
    typename FractionTraits::Denominator_type denominator;
    typename FractionTraits::Decompose()(result, numerator,
                                         denominator);
    if (denominator <= 0) {
        throw InvalidAlgebraicPolynomialError(
            "full-circle pair resultant denominator is not "
            "positive");
    }
    if (CGAL::is_zero(numerator)) {
        return numerator;
    }
    return typename CGAL::Polynomial_traits_d<
        Polynomial1Z>::Canonicalize()(numerator);
}

Polynomial1Z compose_global_chart(
    const Polynomial1Z& local,
    std::size_t chart)
{
    if (chart >= CENTER_CHART_IDS.size()) {
        throw ChartCoverageError(
            "full-circle pair projection uses an unknown "
            "center chart");
    }
    if (CGAL::is_zero(local)) {
        return local;
    }
    const std::array<Integer, 2> affine_coefficients{
        -Integer(chart),
        Integer(4),
    };
    const Polynomial1Z affine(affine_coefficients.begin(),
                              affine_coefficients.end());
    Polynomial1Z result(Integer(0));
    for (int degree = CGAL::degree(local); degree >= 0;
         --degree) {
        result =
            result * affine + Polynomial1Z(local[degree]);
    }
    return typename CGAL::Polynomial_traits_d<
        Polynomial1Z>::Canonicalize()(result);
}

std::vector<std::string> coefficients(
    const Polynomial1Z& polynomial)
{
    if (CGAL::is_zero(polynomial)) {
        return {"0"};
    }
    std::vector<std::string> result;
    result.reserve(static_cast<std::size_t>(
        CGAL::degree(polynomial) + 1));
    for (int index = 0; index <= CGAL::degree(polynomial);
         ++index) {
        result.push_back(
            polynomial[index].convert_to<std::string>());
    }
    return result;
}

FullCircleFeatureSource decode_source(
    const std::string& encoded,
    const std::string& support_kind)
{
    const std::vector<std::string> fields =
        decode_string_sequence(encoded);
    if (fields.size() != 10) {
        throw EventPartitionVerificationError(
            "full-circle pair source is malformed");
    }
    return {
        fields[0],
        fields[1],
        fields[2],
        support_kind,
    };
}

std::vector<FullCircleFeatureSource> feature_sources(
    const std::vector<std::string>& line_sources,
    const std::vector<std::string>& circle_sources)
{
    std::vector<FullCircleFeatureSource> result;
    result.reserve(line_sources.size()
                   + circle_sources.size());
    for (const std::string& source : line_sources) {
        result.push_back(decode_source(source, "line"));
    }
    for (const std::string& source : circle_sources) {
        result.push_back(decode_source(source, "circle"));
    }
    std::sort(result.begin(), result.end(),
              [](const FullCircleFeatureSource& first,
                 const FullCircleFeatureSource& second) {
                  return std::tie(first.feature_id,
                                  first.support_id,
                                  first.trim_id,
                                  first.support_kind)
                         < std::tie(second.feature_id,
                                    second.support_id,
                                    second.trim_id,
                                    second.support_kind);
              });
    const auto duplicate = std::adjacent_find(
        result.begin(), result.end(),
        [](const FullCircleFeatureSource& first,
           const FullCircleFeatureSource& second) {
            return first.feature_id == second.feature_id;
        });
    if (duplicate != result.end()) {
        throw EventPartitionVerificationError(
            "full-circle pair source repeats one feature");
    }
    return result;
}

const ProjectionRecord2& pullback_for(
    const std::vector<ProjectionRecord2>& pullbacks,
    const FullCircleFeatureSource& source,
    const std::string& center_chart,
    const std::string& rim_chart)
{
    const std::string projection_id =
        encode_string_sequence({
            source.support_kind == "line"
                ? "full-circle-line-pullback-v1"
                : "full-circle-circle-pullback-v1",
            source.feature_id,
            center_chart,
            rim_chart,
        });
    const auto found = std::find_if(
        pullbacks.begin(), pullbacks.end(),
        [&projection_id](
            const ProjectionRecord2& projection) {
            return projection.projection_id
                   == projection_id;
        });
    if (found == pullbacks.end()) {
        throw InvalidAlgebraicPolynomialError(
            "full-circle pair projection is missing its "
            "pullback");
    }
    return *found;
}

ProjectionRecord2 constant_record(
    const std::string& projection_id,
    const std::vector<std::string>& values)
{
    return {
        projection_id,
        {values},
        {},
        0,
        0,
        0,
        0,
        "full-circle-pair-invariant-constant-(0,0)-v1",
        encode_string_sequence(values),
    };
}

void append_pair_projection(
    FullCirclePairProjectionBundle2& result,
    const std::string& projection_tag,
    const FullCircleFeatureSource& first,
    const FullCircleFeatureSource& second,
    const std::string& center_chart,
    std::size_t center_chart_index,
    const std::string& first_rim_chart,
    const std::string& second_rim_chart,
    const Polynomial1Z& local)
{
    const std::string projection_id =
        encode_string_sequence({
            projection_tag,
            first.feature_id,
            second.feature_id,
            first.support_id,
            second.support_id,
            center_chart,
            first_rim_chart,
            second_rim_chart,
        });
    const Polynomial1Z global =
        compose_global_chart(local, center_chart_index);
    const std::vector<std::string> values =
        coefficients(global);
    if (CGAL::is_zero(global)) {
        return;
    }
    if (CGAL::degree(global) == 0) {
        result.constants.push_back(
            constant_record(projection_id, values));
        return;
    }
    result.projections.emplace_back(
        projection_id, values,
        std::vector<PartitionEvent2>{});
}

} // namespace

FullCirclePairProjectionBundle2
derive_full_circle_pair_cap_projections(
    const std::vector<std::string>& line_sources,
    const std::vector<std::string>& circle_sources,
    const std::string& cap_chord_ratio,
    const std::vector<ProjectionRecord2>& pullbacks)
{
    const Rational cap_ratio =
        parse_rational(cap_chord_ratio,
                       "full-circle pair cap chord ratio");
    if (cap_ratio <= 0 || cap_ratio > 4) {
        throw InvalidAlgebraicPolynomialError(
            "full-circle pair cap chord ratio lies outside "
            "(0, 4]");
    }
    const std::vector<FullCircleFeatureSource> sources =
        feature_sources(line_sources, circle_sources);
    FullCirclePairProjectionBundle2 result;
    for (std::size_t first_index = 0;
         first_index < sources.size(); ++first_index) {
        for (std::size_t second_index = first_index + 1;
             second_index < sources.size();
             ++second_index) {
            const FullCircleFeatureSource& first =
                sources[first_index];
            const FullCircleFeatureSource& second =
                sources[second_index];
            if (first.support_kind == second.support_kind
                && first.support_id == second.support_id) {
                continue;
            }
            for (std::size_t center_index = 0;
                 center_index < CENTER_CHART_IDS.size();
                 ++center_index) {
                const std::string center_chart =
                    CENTER_CHART_IDS[center_index];
                for (const std::string first_rim_chart :
                     RIM_CHART_IDS) {
                    for (const std::string
                             second_rim_chart :
                         RIM_CHART_IDS) {
                        const Polynomial3Q first_pullback =
                            pullback(pullback_for(
                                         pullbacks, first,
                                         center_chart,
                                         first_rim_chart),
                                     1);
                        const Polynomial3Q second_pullback =
                            pullback(pullback_for(
                                         pullbacks, second,
                                         center_chart,
                                         second_rim_chart),
                                     2);
                        append_pair_projection(
                            result,
                            "full-circle-pair-orientation-"
                            "v1",
                            first, second, center_chart,
                            center_index, first_rim_chart,
                            second_rim_chart,
                            eliminate(
                                first_pullback,
                                second_pullback,
                                orientation_predicate(
                                    rim_coordinates(
                                        1, first_rim_chart),
                                    rim_coordinates(
                                        2,
                                        second_rim_chart))));
                        if (cap_ratio < Rational(4)) {
                            append_pair_projection(
                                result,
                                "full-circle-pair-cap-v2",
                                first, second, center_chart,
                                center_index,
                                first_rim_chart,
                                second_rim_chart,
                                eliminate(
                                    first_pullback,
                                    second_pullback,
                                    cap_predicate(
                                        rim_coordinates(
                                            1,
                                            first_rim_chart),
                                        rim_coordinates(
                                            2,
                                            second_rim_chart),
                                        cap_ratio)));
                        }
                    }
                }
            }
        }
    }
    return result;
}
