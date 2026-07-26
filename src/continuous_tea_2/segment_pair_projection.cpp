#include "segment_pair_projection.h"

#include "event_certificate.h"

#include <algorithm>
#include <string_view>

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

Polynomial3Q variable(int index)
{
    return CGAL::shift(
        Polynomial3Q(Rational(1)),
        1,
        index);
}

Polynomial3Q pullback(
    const ProjectionRecord2& projection,
    int rim_variable)
{
    const Polynomial3Q t = variable(0);
    const Polynomial3Q rim = variable(rim_variable);
    Polynomial3Q result(Rational(0));
    Polynomial3Q t_power(Rational(1));
    for (const auto& row : projection.coefficient_rows) {
        Polynomial3Q rim_power(Rational(1));
        for (const std::string& coefficient : row) {
            result += Rational(Integer(coefficient))
                * t_power * rim_power;
            rim_power *= rim;
        }
        t_power *= t;
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
    const Rational sign =
        chart == "rim-half-0-v1"
        ? Rational(1)
        : Rational(-1);
    const Polynomial3Q parameter =
        variable(variable_index);
    const Polynomial3Q square =
        parameter * parameter;
    return {
        sign * (Polynomial3Q(Rational(1)) - square),
        Rational(2) * sign * parameter,
        Polynomial3Q(Rational(1)) + square,
    };
}

Polynomial3Q determinant(
    const RimCoordinates& first,
    const RimCoordinates& second)
{
    return first.x * second.y
        - first.y * second.x;
}

Polynomial3Q cap_predicate(
    const RimCoordinates& first,
    const RimCoordinates& second,
    const Rational& cap_ratio)
{
    const Polynomial3Q chord_x =
        first.x * second.denominator
        - second.x * first.denominator;
    const Polynomial3Q chord_y =
        first.y * second.denominator
        - second.y * first.denominator;
    return chord_x * chord_x
        + chord_y * chord_y
        - cap_ratio
            * first.denominator
            * first.denominator
            * second.denominator
            * second.denominator;
}

Polynomial1Z eliminate(
    const Polynomial3Q& first_pullback,
    const Polynomial3Q& second_pullback,
    const Polynomial3Q& predicate)
{
    const Polynomial2Q inner =
        typename CGAL::Polynomial_traits_d<
            Polynomial3Q>::Resultant()(
            second_pullback,
            predicate);
    const Polynomial1Q result =
        typename CGAL::Polynomial_traits_d<
            Polynomial2Q>::Resultant()(
            first_pullback[0],
            inner);
    using FractionTraits =
        CGAL::Fraction_traits<Polynomial1Q>;
    typename FractionTraits::Numerator_type numerator;
    typename FractionTraits::Denominator_type denominator;
    typename FractionTraits::Decompose()(
        result,
        numerator,
        denominator);
    if (denominator <= 0) {
        throw InvalidAlgebraicPolynomialError(
            "pair resultant denominator is not positive");
    }
    if (CGAL::is_zero(numerator)) {
        return numerator;
    }
    return typename CGAL::Polynomial_traits_d<
        Polynomial1Z>::Canonicalize()(numerator);
}

std::vector<std::string> coefficients(
    const Polynomial1Z& polynomial)
{
    if (CGAL::is_zero(polynomial)) {
        return {"0"};
    }
    std::vector<std::string> result;
    result.reserve(
        static_cast<std::size_t>(
            CGAL::degree(polynomial) + 1));
    for (int index = 0;
         index <= CGAL::degree(polynomial);
         ++index) {
        result.push_back(
            polynomial[index]
                .convert_to<std::string>());
    }
    return result;
}

PartitionEvent2 pair_event(
    const BoundaryFeatureRecord2& first,
    const BoundaryFeatureRecord2& second,
    const std::string& first_chart,
    const std::string& second_chart,
    const std::string& kind,
    const std::string& disposition)
{
    PartitionEvent2 event(
        kind,
        first.feature_id,
        first.support_id,
        first.trim_predicate,
        {},
        {},
        disposition);
    event.first_feature_id = first.feature_id;
    event.second_feature_id = second.feature_id;
    event.first_chart_id = first_chart;
    event.second_chart_id = second_chart;
    return event;
}

ProjectionRecord2 constant_record(
    const std::string& projection_id,
    const std::vector<std::string>& values,
    int bound)
{
    return {
        projection_id,
        {values},
        {},
        0,
        0,
        bound,
        0,
        "segment-pair-resultant-("
            + std::to_string(bound)
            + ",0)-v1",
        encode_string_sequence(values),
    };
}

const ProjectionRecord2& pullback_for(
    const std::vector<ProjectionRecord2>& pullbacks,
    const std::string& feature_id,
    const std::string& chart)
{
    const std::string id =
        encode_string_sequence(
            {
                "segment-pullback-v1",
                feature_id,
                chart,
            });
    const auto found = std::find_if(
        pullbacks.begin(),
        pullbacks.end(),
        [&id](const ProjectionRecord2& projection) {
            return projection.projection_id == id;
        });
    if (found == pullbacks.end()) {
        throw InvalidAlgebraicPolynomialError(
            "pair projection is missing its pullback");
    }
    return *found;
}

void append_pair_projection(
    SegmentPairProjectionBundle2& result,
    const BoundaryFeatureRecord2& first,
    const BoundaryFeatureRecord2& second,
    const std::string& first_chart,
    const std::string& second_chart,
    const std::string& kind,
    const std::string& disposition,
    const Polynomial1Z& polynomial,
    int bound)
{
    const std::string projection_id =
        encode_string_sequence(
            {
                kind == "cap-crossing"
                    ? "segment-pair-cap-v1"
                    : "segment-pair-endpoint-v1",
                first.feature_id,
                second.feature_id,
                first.support_id,
                second.support_id,
                first_chart,
                second_chart,
            });
    std::vector<PartitionEvent2> events{
        pair_event(
            first,
            second,
            first_chart,
            second_chart,
            kind,
            disposition),
        pair_event(
            second,
            first,
            second_chart,
            first_chart,
            kind,
            disposition),
    };
    const std::vector<std::string> values =
        coefficients(polynomial);
    if (CGAL::is_zero(polynomial)) {
        if (kind == "cap-crossing") {
            result.overlaps.push_back(
                {
                    "identically-equal-cap-interval",
                    "0",
                    "1",
                    "1",
                    "2",
                    "ccw",
                });
        }
        return;
    }
    if (CGAL::degree(polynomial) == 0) {
        result.constants.push_back(
            constant_record(
                projection_id,
                values,
                bound));
        return;
    }
    result.projections.emplace_back(
        projection_id,
        values,
        std::move(events));
}

} // namespace

SegmentPairProjectionBundle2
derive_segment_pair_projections(
    const std::vector<BoundaryFeatureRecord2>& records,
    const SegmentEventSource2& source,
    const std::vector<ProjectionRecord2>& pullbacks)
{
    SegmentPairProjectionBundle2 result;
    const Rational cap_ratio =
        parse_rational(
            source.cap_chord_ratio().text(),
            "cap chord ratio");
    const std::vector<std::string> charts{
        "rim-half-0-v1",
        "rim-half-1-v1",
    };
    for (std::size_t first_index = 0;
         first_index < records.size();
         ++first_index) {
        for (std::size_t second_index = first_index + 1;
             second_index < records.size();
             ++second_index) {
            const BoundaryFeatureRecord2& first =
                records[first_index];
            const BoundaryFeatureRecord2& second =
                records[second_index];
            for (const std::string& first_chart : charts) {
                for (const std::string& second_chart :
                     charts) {
                    const Polynomial3Q first_pullback =
                        pullback(
                            pullback_for(
                                pullbacks,
                                first.feature_id,
                                first_chart),
                            1);
                    const Polynomial3Q second_pullback =
                        pullback(
                            pullback_for(
                                pullbacks,
                                second.feature_id,
                                second_chart),
                            2);
                    const RimCoordinates first_rim =
                        rim_coordinates(1, first_chart);
                    const RimCoordinates second_rim =
                        rim_coordinates(2, second_chart);
                    append_pair_projection(
                        result,
                        first,
                        second,
                        first_chart,
                        second_chart,
                        "cap-crossing",
                        "equal-cap",
                        eliminate(
                            first_pullback,
                            second_pullback,
                            cap_predicate(
                                first_rim,
                                second_rim,
                                cap_ratio)),
                        64);
                    append_pair_projection(
                        result,
                        first,
                        second,
                        first_chart,
                        second_chart,
                        "endpoint-order",
                        "merge",
                        eliminate(
                            first_pullback,
                            second_pullback,
                            determinant(
                                first_rim,
                                second_rim)),
                        32);
                }
            }
        }
    }
    return result;
}
