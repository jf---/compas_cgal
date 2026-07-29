#include "circle_pair_projection.h"

#include "event_certificate.h"

#include <algorithm>
#include <array>
#include <map>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include <CGAL/CORE/BigRat.h>
#include <CGAL/Algebraic_structure_traits.h>
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
    int rim_variable,
    const std::string& support_kind)
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
    if (support_kind == "line") {
        return result;
    }
    if (support_kind != "circle") {
        throw InvalidAlgebraicPolynomialError(
            "full-circle pair pullback uses an unknown support kind");
    }
    const Polynomial3Q rim_denominator =
        Polynomial3Q(Rational(1)) + rim * rim;
    Polynomial3Q saturated;
    using PolynomialTraits =
        CGAL::Algebraic_structure_traits<
            Polynomial3Q>;
    if (!typename PolynomialTraits::Divides()(
            rim_denominator,
            result,
            saturated)) {
        throw InvalidAlgebraicPolynomialError(
            "circle pair pullback lacks its exact rim denominator factor");
    }
    return saturated;
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

FullCirclePairRequest2::FullCirclePairRequest2(
    std::string center_chart_id,
    std::string first_feature_id,
    std::string first_support_id,
    std::string first_support_kind,
    std::string first_rim_chart_id,
    std::string second_feature_id,
    std::string second_support_id,
    std::string second_support_kind,
    std::string second_rim_chart_id)
    : center_chart_id_(std::move(center_chart_id)),
      first_feature_id_(std::move(first_feature_id)),
      first_support_id_(std::move(first_support_id)),
      first_support_kind_(std::move(first_support_kind)),
      first_rim_chart_id_(std::move(first_rim_chart_id)),
      second_feature_id_(std::move(second_feature_id)),
      second_support_id_(std::move(second_support_id)),
      second_support_kind_(std::move(second_support_kind)),
      second_rim_chart_id_(std::move(second_rim_chart_id))
{
}

FullCirclePairRequest2
FullCirclePairRequest2::build(
    std::string center_chart_id,
    std::string first_feature_id,
    std::string first_support_id,
    std::string first_support_kind,
    std::string first_rim_chart_id,
    std::string second_feature_id,
    std::string second_support_id,
    std::string second_support_kind,
    std::string second_rim_chart_id)
{
    const auto known_center_chart =
        std::find(
            CENTER_CHART_IDS.begin(),
            CENTER_CHART_IDS.end(),
            center_chart_id)
        != CENTER_CHART_IDS.end();
    const auto known_rim_chart =
        [](const std::string& chart_id) {
            return std::find(
                       RIM_CHART_IDS.begin(),
                       RIM_CHART_IDS.end(),
                       chart_id)
                != RIM_CHART_IDS.end();
        };
    const auto known_support_kind =
        [](const std::string& support_kind) {
            return support_kind == "line"
                || support_kind == "circle";
        };
    if (!known_center_chart
        || !known_rim_chart(first_rim_chart_id)
        || !known_rim_chart(second_rim_chart_id)) {
        throw InvalidFullCirclePairRequestError(
            "full-circle pair request uses an unknown chart");
    }
    if (first_feature_id.empty()
        || first_support_id.empty()
        || second_feature_id.empty()
        || second_support_id.empty()
        || !known_support_kind(first_support_kind)
        || !known_support_kind(second_support_kind)) {
        throw InvalidFullCirclePairRequestError(
            "full-circle pair request has malformed feature identity");
    }
    if (first_feature_id == second_feature_id) {
        throw InvalidFullCirclePairRequestError(
            "full-circle pair request repeats one feature");
    }
    if (first_support_kind == second_support_kind
        && first_support_id == second_support_id) {
        throw InvalidFullCirclePairRequestError(
            "full-circle pair request requires distinct supports");
    }
    const auto first_key =
        std::tie(
            first_feature_id,
            first_support_kind,
            first_support_id,
            first_rim_chart_id);
    const auto second_key =
        std::tie(
            second_feature_id,
            second_support_kind,
            second_support_id,
            second_rim_chart_id);
    if (second_key < first_key) {
        std::swap(first_feature_id, second_feature_id);
        std::swap(first_support_id, second_support_id);
        std::swap(first_support_kind, second_support_kind);
        std::swap(first_rim_chart_id, second_rim_chart_id);
    }
    return FullCirclePairRequest2(
        std::move(center_chart_id),
        std::move(first_feature_id),
        std::move(first_support_id),
        std::move(first_support_kind),
        std::move(first_rim_chart_id),
        std::move(second_feature_id),
        std::move(second_support_id),
        std::move(second_support_kind),
        std::move(second_rim_chart_id));
}

FullCirclePairRequest2
FullCirclePairRequest2::decode(
    const std::string& canonical_source)
{
    std::vector<std::string> fields;
    try {
        fields = decode_string_sequence(
            canonical_source);
    } catch (const EventSubstrateError&) {
        throw InvalidFullCirclePairRequestError(
            "full-circle pair request is not one encoded sequence");
    }
    if (fields.size() != 10
        || fields[0]
            != "full-circle-active-pair-request-v1") {
        throw InvalidFullCirclePairRequestError(
            "full-circle pair request has malformed source fields");
    }
    FullCirclePairRequest2 result =
        FullCirclePairRequest2::build(
            fields[1],
            fields[2],
            fields[3],
            fields[4],
            fields[5],
            fields[6],
            fields[7],
            fields[8],
            fields[9]);
    if (result.canonical_source()
        != canonical_source) {
        throw InvalidFullCirclePairRequestError(
            "full-circle pair request is not canonical");
    }
    return result;
}

const std::string&
FullCirclePairRequest2::center_chart_id() const
{
    return center_chart_id_;
}

const std::string&
FullCirclePairRequest2::first_feature_id() const
{
    return first_feature_id_;
}

const std::string&
FullCirclePairRequest2::first_support_id() const
{
    return first_support_id_;
}

const std::string&
FullCirclePairRequest2::first_support_kind() const
{
    return first_support_kind_;
}

const std::string&
FullCirclePairRequest2::first_rim_chart_id() const
{
    return first_rim_chart_id_;
}

const std::string&
FullCirclePairRequest2::second_feature_id() const
{
    return second_feature_id_;
}

const std::string&
FullCirclePairRequest2::second_support_id() const
{
    return second_support_id_;
}

const std::string&
FullCirclePairRequest2::second_support_kind() const
{
    return second_support_kind_;
}

const std::string&
FullCirclePairRequest2::second_rim_chart_id() const
{
    return second_rim_chart_id_;
}

std::string
FullCirclePairRequest2::canonical_source() const
{
    return encode_string_sequence(
        {
            "full-circle-active-pair-request-v1",
            center_chart_id_,
            first_feature_id_,
            first_support_id_,
            first_support_kind_,
            first_rim_chart_id_,
            second_feature_id_,
            second_support_id_,
            second_support_kind_,
            second_rim_chart_id_,
        });
}

std::vector<std::string>
full_circle_exhaustive_pair_requests(
    const std::vector<std::string>& line_sources,
    const std::vector<std::string>& circle_sources)
{
    const std::vector<FullCircleFeatureSource> sources =
        feature_sources(line_sources, circle_sources);
    std::vector<std::string> requests;
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
            for (const std::string center_chart :
                 CENTER_CHART_IDS) {
                for (const std::string first_rim_chart :
                     RIM_CHART_IDS) {
                    for (const std::string
                             second_rim_chart :
                         RIM_CHART_IDS) {
                        requests.push_back(
                            FullCirclePairRequest2::build(
                                center_chart,
                                first.feature_id,
                                first.support_id,
                                first.support_kind,
                                first_rim_chart,
                                second.feature_id,
                                second.support_id,
                                second.support_kind,
                                second_rim_chart)
                                .canonical_source());
                    }
                }
            }
        }
    }
    std::sort(requests.begin(), requests.end());
    requests.erase(
        std::unique(requests.begin(), requests.end()),
        requests.end());
    return requests;
}

FullCirclePairProjectionBundle2
derive_full_circle_pair_cap_projections(
    const std::vector<std::string>& line_sources,
    const std::vector<std::string>& circle_sources,
    const std::string& cap_chord_ratio,
    const std::vector<ProjectionRecord2>& pullbacks,
    const std::vector<std::string>& pair_requests)
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
    if (!std::is_sorted(
            pair_requests.begin(),
            pair_requests.end())
        || std::adjacent_find(
               pair_requests.begin(),
               pair_requests.end())
            != pair_requests.end()) {
        throw InvalidFullCirclePairRequestError(
            "full-circle pair requests are not a canonical set");
    }
    std::map<std::string, const FullCircleFeatureSource*>
        source_by_feature;
    for (const FullCircleFeatureSource& source :
         sources) {
        source_by_feature.emplace(
            source.feature_id,
            &source);
    }
    FullCirclePairProjectionBundle2 result;
    for (const std::string& encoded_request :
         pair_requests) {
        const FullCirclePairRequest2 request =
            FullCirclePairRequest2::decode(
                encoded_request);
        const auto first_found =
            source_by_feature.find(
                request.first_feature_id());
        const auto second_found =
            source_by_feature.find(
                request.second_feature_id());
        if (first_found == source_by_feature.end()
            || second_found
                == source_by_feature.end()) {
            throw InvalidFullCirclePairRequestError(
                "full-circle pair request names an unknown feature");
        }
        const FullCircleFeatureSource& first =
            *first_found->second;
        const FullCircleFeatureSource& second =
            *second_found->second;
        if (first.support_id
                != request.first_support_id()
            || first.support_kind
                != request.first_support_kind()
            || second.support_id
                != request.second_support_id()
            || second.support_kind
                != request.second_support_kind()) {
            throw InvalidFullCirclePairRequestError(
                "full-circle pair request has inconsistent support identity");
        }
        const auto center_found =
            std::find(
                CENTER_CHART_IDS.begin(),
                CENTER_CHART_IDS.end(),
                request.center_chart_id());
        const std::size_t center_index =
            static_cast<std::size_t>(
                std::distance(
                    CENTER_CHART_IDS.begin(),
                    center_found));
        const Polynomial3Q first_pullback =
            pullback(
                pullback_for(
                    pullbacks,
                    first,
                    request.center_chart_id(),
                    request.first_rim_chart_id()),
                1,
                first.support_kind);
        const Polynomial3Q second_pullback =
            pullback(
                pullback_for(
                    pullbacks,
                    second,
                    request.center_chart_id(),
                    request.second_rim_chart_id()),
                2,
                second.support_kind);
        const RimCoordinates first_rim =
            rim_coordinates(
                1,
                request.first_rim_chart_id());
        const RimCoordinates second_rim =
            rim_coordinates(
                2,
                request.second_rim_chart_id());
        append_pair_projection(
            result,
            "full-circle-pair-orientation-v1",
            first,
            second,
            request.center_chart_id(),
            center_index,
            request.first_rim_chart_id(),
            request.second_rim_chart_id(),
            eliminate(
                first_pullback,
                second_pullback,
                orientation_predicate(
                    first_rim,
                    second_rim)));
        if (cap_ratio < Rational(4)) {
            append_pair_projection(
                result,
                "full-circle-pair-cap-v2",
                first,
                second,
                request.center_chart_id(),
                center_index,
                request.first_rim_chart_id(),
                request.second_rim_chart_id(),
                eliminate(
                    first_pullback,
                    second_pullback,
                    cap_predicate(
                        first_rim,
                        second_rim,
                        cap_ratio)));
        }
    }
    return result;
}
