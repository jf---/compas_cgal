#include "segment_projection.h"

#include "event_certificate.h"
#include "parameter_charts.h"
#include "segment_pair_projection.h"
#include "segment_partition.h"

#include <algorithm>
#include <iterator>
#include <map>
#include <optional>
#include <sstream>
#include <string_view>
#include <utility>

#include <CGAL/CORE/BigRat.h>
#include <CGAL/Fraction_traits.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/number_utils.h>

namespace {

using Integer = CORE::BigInt;
using Rational = CORE::BigRat;
using Polynomial = std::vector<Rational>;
using RationalPolynomial = CGAL::Polynomial<Rational>;

constexpr std::size_t
    PAIR_RESULTANTS_PER_FEATURE_PAIR = 8;

struct CoordinateParts {
    Rational base;
    Rational radical_coefficient;
    Rational radicand;
};

struct VertexInput {
    GpsPoint point;
    std::vector<const BoundaryFeatureRecord2*> incidents;
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

Rational exact_rational(const Epeck::FT& value)
{
    return CGAL::exact(value);
}

CoordinateParts coordinate_parts(
    const GpsPoint::CoordNT& coordinate)
{
    return {
        exact_rational(coordinate.a0()),
        exact_rational(coordinate.a1()),
        exact_rational(coordinate.root()),
    };
}

Polynomial trim(Polynomial polynomial)
{
    while (polynomial.size() > 1
           && polynomial.back() == 0) {
        polynomial.pop_back();
    }
    return polynomial;
}

Polynomial add(
    const Polynomial& left,
    const Polynomial& right)
{
    Polynomial result(
        std::max(left.size(), right.size()),
        Rational(0));
    for (std::size_t index = 0;
         index < result.size();
         ++index) {
        if (index < left.size()) {
            result[index] += left[index];
        }
        if (index < right.size()) {
            result[index] += right[index];
        }
    }
    return trim(std::move(result));
}

Polynomial scale(
    const Polynomial& polynomial,
    const Rational& factor)
{
    Polynomial result = polynomial;
    for (Rational& coefficient : result) {
        coefficient *= factor;
    }
    return trim(std::move(result));
}

Polynomial subtract(
    const Polynomial& left,
    const Polynomial& right)
{
    return add(left, scale(right, Rational(-1)));
}

Polynomial multiply(
    const Polynomial& left,
    const Polynomial& right)
{
    Polynomial result(
        left.size() + right.size() - 1,
        Rational(0));
    for (std::size_t first = 0;
         first < left.size();
         ++first) {
        for (std::size_t second = 0;
             second < right.size();
             ++second) {
            result[first + second] +=
                left[first] * right[second];
        }
    }
    return trim(std::move(result));
}

RationalPolynomial as_cgal_polynomial(
    const Polynomial& coefficients)
{
    return RationalPolynomial(
        coefficients.begin(),
        coefficients.end());
}

Polynomial from_cgal_polynomial(
    const RationalPolynomial& polynomial)
{
    const int degree = CGAL::degree(polynomial);
    if (degree < 0) {
        return {Rational(0)};
    }
    Polynomial result;
    result.reserve(
        static_cast<std::size_t>(degree + 1));
    for (int index = 0; index <= degree; ++index) {
        result.push_back(polynomial[index]);
    }
    return trim(std::move(result));
}

Polynomial polynomial_gcd(
    const Polynomial& left,
    const Polynomial& right)
{
    const RationalPolynomial first =
        as_cgal_polynomial(left);
    const RationalPolynomial second =
        as_cgal_polynomial(right);
    if (CGAL::is_zero(first)) {
        return from_cgal_polynomial(second);
    }
    if (CGAL::is_zero(second)) {
        return from_cgal_polynomial(first);
    }
    return from_cgal_polynomial(
        typename CGAL::Polynomial_traits_d<
            RationalPolynomial>::
            Gcd_up_to_constant_factor()(
                first,
                second));
}

Integer greatest_common_divisor(
    Integer left,
    Integer right)
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

std::vector<std::string> primitive_coefficients(
    const Polynomial& raw)
{
    const Polynomial polynomial = trim(raw);
    Integer denominator_lcm = 1;
    for (const Rational& coefficient : polynomial) {
        denominator_lcm = least_common_multiple(
            denominator_lcm,
            CORE::denominator(coefficient));
    }
    std::vector<Integer> coefficients;
    coefficients.reserve(polynomial.size());
    Integer divisor = 0;
    for (const Rational& coefficient : polynomial) {
        const Integer integer = CORE::numerator(
            coefficient * Rational(denominator_lcm));
        coefficients.push_back(integer);
        divisor = greatest_common_divisor(
            divisor,
            integer);
    }
    if (divisor == 0) {
        return {"0"};
    }
    for (Integer& coefficient : coefficients) {
        coefficient /= divisor;
    }
    if (coefficients.back() < 0) {
        for (Integer& coefficient : coefficients) {
            coefficient = -coefficient;
        }
    }
    std::vector<std::string> result;
    result.reserve(coefficients.size());
    for (const Integer& coefficient : coefficients) {
        result.push_back(
            coefficient.convert_to<std::string>());
    }
    return result;
}

Polynomial coordinate_motion(
    const ExactBinary64Rational2& start,
    const ExactBinary64Rational2& end)
{
    const Rational first =
        parse_rational(start.text(), "segment coordinate");
    const Rational second =
        parse_rational(end.text(), "segment coordinate");
    return {first, second - first};
}

std::optional<Rational> rational_square_root(
    const Rational& value)
{
    if (value < 0) {
        return std::nullopt;
    }
    using FractionTraits =
        CGAL::Fraction_traits<Rational>;
    Integer numerator;
    Integer denominator;
    typename FractionTraits::Decompose()(
        value,
        numerator,
        denominator);
    Integer numerator_root;
    Integer denominator_root;
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

PartitionEvent2 support_event(
    const BoundaryFeatureRecord2& record,
    std::string kind,
    std::string vertex_id,
    std::string disposition)
{
    return {
        std::move(kind),
        record.feature_id,
        record.support_id,
        record.trim_predicate,
        std::move(vertex_id),
        {},
        std::move(disposition),
    };
}

Polynomial line_tangency(
    const BoundaryFeatureRecord2& record,
    const Polynomial& center_x,
    const Polynomial& center_y,
    const Rational& radius_squared)
{
    const Rational a =
        parse_rational(
            record.primitive_coefficients[0],
            "line support");
    const Rational b =
        parse_rational(
            record.primitive_coefficients[1],
            "line support");
    const Rational c =
        parse_rational(
            record.primitive_coefficients[2],
            "line support");
    const Polynomial signed_distance =
        add(
            add(
                scale(center_x, a),
                scale(center_y, b)),
            {c});
    return subtract(
        multiply(signed_distance, signed_distance),
        {radius_squared * (a * a + b * b)});
}

std::vector<Polynomial> oriented_line_tangencies(
    const BoundaryFeatureRecord2& record,
    const Polynomial& center_x,
    const Polynomial& center_y,
    const Rational& radius)
{
    const Rational a =
        parse_rational(
            record.primitive_coefficients[0],
            "line support");
    const Rational b =
        parse_rational(
            record.primitive_coefficients[1],
            "line support");
    const Rational c =
        parse_rational(
            record.primitive_coefficients[2],
            "line support");
    const Polynomial signed_distance =
        add(
            add(
                scale(center_x, a),
                scale(center_y, b)),
            {c});
    const std::optional<Rational> normal =
        rational_square_root(a * a + b * b);
    if (!normal.has_value()) {
        return {
            subtract(
                multiply(
                    signed_distance,
                    signed_distance),
                {
                    radius * radius
                        * (a * a + b * b),
                }),
        };
    }
    return {
        add(
            signed_distance,
            {radius * *normal}),
        subtract(
            signed_distance,
            {radius * *normal}),
    };
}

Polynomial line_cap_crossing(
    const BoundaryFeatureRecord2& record,
    const Polynomial& center_x,
    const Polynomial& center_y,
    const Rational& radius_squared,
    const Rational& cap_chord_ratio)
{
    const Rational a =
        parse_rational(
            record.primitive_coefficients[0],
            "line support");
    const Rational b =
        parse_rational(
            record.primitive_coefficients[1],
            "line support");
    const Rational c =
        parse_rational(
            record.primitive_coefficients[2],
            "line support");
    const Rational normal_squared = a * a + b * b;
    const Polynomial signed_distance =
        add(
            add(
                scale(center_x, a),
                scale(center_y, b)),
            {c});
    return subtract(
        {
            (Rational(4) - cap_chord_ratio)
                * radius_squared
                * normal_squared,
        },
        scale(
            multiply(
                signed_distance,
                signed_distance),
            Rational(4)));
}

struct CircleSupport {
    Rational center_x;
    Rational center_y;
    Rational radius_squared;
};

CircleSupport circle_support(
    const BoundaryFeatureRecord2& record)
{
    const Rational quadratic =
        parse_rational(
            record.primitive_coefficients[0],
            "circle support");
    const Rational x_linear =
        parse_rational(
            record.primitive_coefficients[1],
            "circle support");
    const Rational y_linear =
        parse_rational(
            record.primitive_coefficients[2],
            "circle support");
    const Rational constant =
        parse_rational(
            record.primitive_coefficients[3],
            "circle support");
    const Rational center_x =
        -x_linear / (Rational(2) * quadratic);
    const Rational center_y =
        -y_linear / (Rational(2) * quadratic);
    return {
        center_x,
        center_y,
        center_x * center_x
            + center_y * center_y
            - constant / quadratic,
    };
}

Rational minimum_segment_distance_squared(
    const Rational& point_x,
    const Rational& point_y,
    const SegmentEventSource2& source)
{
    const Rational start_x =
        parse_rational(source.x0().text(), "segment x0");
    const Rational start_y =
        parse_rational(source.y0().text(), "segment y0");
    const Rational direction_x =
        parse_rational(source.x1().text(), "segment x1")
        - start_x;
    const Rational direction_y =
        parse_rational(source.y1().text(), "segment y1")
        - start_y;
    const Rational length_squared =
        direction_x * direction_x
        + direction_y * direction_y;
    if (length_squared == 0) {
        throw IncompleteSegmentPartitionError(
            "support reachability requires nondegenerate segment motion");
    }
    const Rational offset_x = point_x - start_x;
    const Rational offset_y = point_y - start_y;
    const Rational projection =
        offset_x * direction_x
        + offset_y * direction_y;
    if (projection <= 0) {
        return offset_x * offset_x
            + offset_y * offset_y;
    }
    if (projection >= length_squared) {
        const Rational end_offset_x =
            point_x - start_x - direction_x;
        const Rational end_offset_y =
            point_y - start_y - direction_y;
        return end_offset_x * end_offset_x
            + end_offset_y * end_offset_y;
    }
    return offset_x * offset_x
        + offset_y * offset_y
        - projection * projection / length_squared;
}

bool support_can_reach_tool(
    const BoundaryFeatureRecord2& record,
    const SegmentEventSource2& source)
{
    if (record.support_kind != "circle") {
        return true;
    }
    const CircleSupport circle =
        circle_support(record);
    const Rational tool_radius =
        parse_rational(
            source.tool_radius().text(),
            "tool radius");
    const Rational distance_squared =
        minimum_segment_distance_squared(
            circle.center_x,
            circle.center_y,
            source);
    const Rational comparison =
        distance_squared
        - circle.radius_squared
        - tool_radius * tool_radius;
    if (comparison <= 0) {
        return true;
    }
    return comparison * comparison
        <= Rational(4)
            * circle.radius_squared
            * tool_radius * tool_radius;
}

Polynomial circle_tangency(
    const BoundaryFeatureRecord2& record,
    const Polynomial& center_x,
    const Polynomial& center_y,
    const Rational& radius_squared)
{
    const CircleSupport circle =
        circle_support(record);
    const Polynomial dx =
        subtract(center_x, {circle.center_x});
    const Polynomial dy =
        subtract(center_y, {circle.center_y});
    const Polynomial distance_squared =
        add(multiply(dx, dx), multiply(dy, dy));
    const Polynomial radial_difference =
        subtract(
            distance_squared,
            {
                circle.radius_squared
                    + radius_squared,
            });
    return subtract(
        multiply(
            radial_difference,
            radial_difference),
        {
            Rational(4)
                * circle.radius_squared
                * radius_squared,
        });
}

Polynomial circle_support_overlap(
    const BoundaryFeatureRecord2& record,
    const Polynomial& center_x,
    const Polynomial& center_y,
    const Rational& radius_squared)
{
    const CircleSupport circle =
        circle_support(record);
    if (circle.radius_squared != radius_squared) {
        return {Rational(1)};
    }
    return polynomial_gcd(
        subtract(center_x, {circle.center_x}),
        subtract(center_y, {circle.center_y}));
}

Polynomial circle_cap_crossing(
    const BoundaryFeatureRecord2& record,
    const Polynomial& center_x,
    const Polynomial& center_y,
    const Rational& radius_squared,
    const Rational& cap_chord_ratio)
{
    const CircleSupport circle =
        circle_support(record);
    const Polynomial dx =
        subtract(center_x, {circle.center_x});
    const Polynomial dy =
        subtract(center_y, {circle.center_y});
    const Polynomial distance_squared =
        add(multiply(dx, dx), multiply(dy, dy));
    const Polynomial radical_axis =
        add(
            distance_squared,
            {
                radius_squared
                    - circle.radius_squared,
            });
    return subtract(
        scale(
            distance_squared,
            (Rational(4) - cap_chord_ratio)
                * radius_squared),
        multiply(radical_axis, radical_axis));
}

Polynomial vertex_passage(
    const GpsPoint& point,
    const Polynomial& center_x,
    const Polynomial& center_y,
    const Rational& radius_squared)
{
    const CoordinateParts x =
        coordinate_parts(point.x());
    const CoordinateParts y =
        coordinate_parts(point.y());
    Rational radicand = 0;
    if (x.radical_coefficient != 0) {
        radicand = x.radicand;
    }
    if (y.radical_coefficient != 0) {
        if (radicand != 0
            && radicand != y.radicand) {
            throw UnsupportedAlgebraicVertexProjectionError(
                "vertex coordinates use distinct radicals");
        }
        radicand = y.radicand;
    }

    const Polynomial dx =
        subtract(center_x, {x.base});
    const Polynomial dy =
        subtract(center_y, {y.base});
    Polynomial rational_part = subtract(
        add(multiply(dx, dx), multiply(dy, dy)),
        {radius_squared});
    Polynomial radical_part = {Rational(0)};
    if (radicand != 0) {
        rational_part = add(
            rational_part,
            {
                radicand
                    * (
                        x.radical_coefficient
                            * x.radical_coefficient
                        + y.radical_coefficient
                            * y.radical_coefficient),
            });
        radical_part = scale(
            add(
                scale(
                    dx,
                    x.radical_coefficient),
                scale(
                    dy,
                    y.radical_coefficient)),
            Rational(-2));
    }
    if (radicand == 0) {
        return rational_part;
    }
    return subtract(
        multiply(rational_part, rational_part),
        scale(
            multiply(
                radical_part,
                radical_part),
            radicand));
}

void append_projection(
    std::vector<ProjectionInput2>& projections,
    std::string projection_id,
    const Polynomial& polynomial,
    std::vector<PartitionEvent2> events)
{
    const std::vector<std::string> coefficients =
        primitive_coefficients(polynomial);
    if (coefficients.size() < 2) {
        return;
    }
    projections.emplace_back(
        std::move(projection_id),
        coefficients,
        std::move(events));
}

std::vector<std::string> pullback_support(
    const BoundaryFeatureRecord2& record)
{
    if (record.support_kind == "line") {
        return record.primitive_coefficients;
    }
    const CircleSupport circle = circle_support(record);
    std::ostringstream radius;
    radius << circle.radius_squared;
    std::ostringstream center_x;
    center_x << circle.center_x;
    std::ostringstream center_y;
    center_y << circle.center_y;
    return {
        center_x.str(),
        center_y.str(),
        radius.str(),
    };
}

} // namespace

SegmentProjectionBundle2 derive_segment_projections(
    const std::vector<BoundaryFeatureRecord2>& records,
    const SegmentEventSource2& source)
{
    SegmentProjectionBundle2 result;
    std::vector<BoundaryFeatureRecord2> event_records;
    event_records.reserve(records.size());
    std::copy_if(
        records.begin(),
        records.end(),
        std::back_inserter(event_records),
        [&source](const BoundaryFeatureRecord2& record) {
            return support_can_reach_tool(
                record,
                source);
        });
    if (event_records.empty()) {
        throw IncompleteSegmentPartitionError(
            "segment projection reachability removed every boundary support");
    }
    const std::size_t unfiltered_pair_resultants =
        PAIR_RESULTANTS_PER_FEATURE_PAIR
        * records.size() * (records.size() - 1)
        / 2;
    const std::size_t scheduled_pair_resultants =
        PAIR_RESULTANTS_PER_FEATURE_PAIR
        * event_records.size()
        * (event_records.size() - 1)
        / 2;
    if (event_records.size() < records.size()
        && scheduled_pair_resultants
            >= unfiltered_pair_resultants) {
        throw IncompleteSegmentPartitionError(
            "exact support filtering did not reduce pair-resultant work");
    }
    const Polynomial center_x =
        coordinate_motion(source.x0(), source.x1());
    const Polynomial center_y =
        coordinate_motion(source.y0(), source.y1());
    const Rational radius =
        parse_rational(
            source.tool_radius().text(),
            "tool radius");
    const Rational radius_squared = radius * radius;

    std::map<std::string, VertexInput> vertices;
    for (const BoundaryFeatureRecord2& record :
         event_records) {
        const std::vector<std::string> support =
            pullback_support(record);
        for (const std::string& rim_chart :
             {"rim-half-0-v1", "rim-half-1-v1"}) {
            ProjectionRecord2 pullback =
                construct_pullback(
                    "segment",
                    source.motion_data(),
                    record.support_kind,
                    support,
                    source.tool_radius().text(),
                    {},
                    rim_chart);
            pullback.projection_id =
                encode_string_sequence(
                    {
                        "segment-pullback-v1",
                        record.feature_id,
                        rim_chart,
                    });
            result.pullbacks.push_back(
                std::move(pullback));
        }

        const Polynomial tangency =
            record.support_kind == "line"
            ? line_tangency(
                  record,
                  center_x,
                  center_y,
                  radius_squared)
            : circle_tangency(
                  record,
                  center_x,
                  center_y,
                  radius_squared);
        append_projection(
            result.event_projections,
            encode_string_sequence(
                {
                    "segment-tangency-v1",
                    record.feature_id,
                }),
            tangency,
            {
                support_event(
                    record,
                    "tangent",
                    {},
                    "contact"),
            });

        if (record.support_kind == "line") {
            const std::vector<Polynomial>
                orientation_boundaries =
                    oriented_line_tangencies(
                        record,
                        center_x,
                        center_y,
                        radius);
            for (std::size_t orientation_index = 0;
                 orientation_index
                     < orientation_boundaries.size();
                 ++orientation_index) {
                append_projection(
                    result.event_projections,
                    encode_string_sequence(
                        {
                            "segment-pair-orientation-v1",
                            record.feature_id,
                            std::to_string(
                                orientation_index),
                        }),
                    orientation_boundaries[
                        orientation_index],
                    {
                        support_event(
                            record,
                            "tangent",
                            {},
                            "contact"),
                    });
            }
            append_projection(
                result.event_projections,
                encode_string_sequence(
                    {
                        "segment-cap-crossing-v1",
                        record.feature_id,
                    }),
                line_cap_crossing(
                    record,
                    center_x,
                    center_y,
                    radius_squared,
                    parse_rational(
                        source.cap_chord_ratio().text(),
                        "cap chord ratio")),
                {
                    support_event(
                        record,
                        "cap-crossing",
                        {},
                        "equal-cap"),
                });
        } else {
            append_projection(
                result.event_projections,
                encode_string_sequence(
                    {
                        "segment-support-overlap-v1",
                        record.feature_id,
                    }),
                circle_support_overlap(
                    record,
                    center_x,
                    center_y,
                    radius_squared),
                {
                    support_event(
                        record,
                        "support-overlap",
                        {},
                        "coincident-support"),
                });
            append_projection(
                result.event_projections,
                encode_string_sequence(
                    {
                        "segment-cap-crossing-v1",
                        record.feature_id,
                    }),
                circle_cap_crossing(
                    record,
                    center_x,
                    center_y,
                    radius_squared,
                    parse_rational(
                        source.cap_chord_ratio().text(),
                        "cap chord ratio")),
                {
                    support_event(
                        record,
                        "cap-crossing",
                        {},
                        "equal-cap"),
                });
        }

        auto [source_vertex, source_inserted] =
            vertices.try_emplace(
                record.source_vertex_id,
                VertexInput{
                    record.curve.source(),
                    {},
                });
        source_vertex->second.incidents.push_back(
            &record);
        auto [target_vertex, target_inserted] =
            vertices.try_emplace(
                record.target_vertex_id,
                VertexInput{
                    record.curve.target(),
                    {},
                });
        target_vertex->second.incidents.push_back(
            &record);
        static_cast<void>(source_inserted);
        static_cast<void>(target_inserted);
    }

    for (const auto& [vertex_id, vertex] : vertices) {
        const Polynomial passage =
            vertex_passage(
                vertex.point,
                center_x,
                center_y,
                radius_squared);
        std::vector<PartitionEvent2> events;
        events.reserve(vertex.incidents.size());
        for (const BoundaryFeatureRecord2* incident :
             vertex.incidents) {
            events.push_back(
                support_event(
                    *incident,
                    "endpoint-order",
                    vertex_id,
                    "merge"));
        }
        append_projection(
            result.event_projections,
            encode_string_sequence(
                {
                    "segment-vertex-passage-v1",
                    vertex_id,
                }),
            passage,
            events);
        for (const BoundaryFeatureRecord2* incident :
             vertex.incidents) {
            if (incident->support_kind != "line") {
                continue;
            }
            const std::vector<Polynomial>
                endpoint_boundaries =
                    oriented_line_tangencies(
                        *incident,
                        center_x,
                        center_y,
                        radius);
            for (std::size_t endpoint_index = 0;
                 endpoint_index
                     < endpoint_boundaries.size();
                 ++endpoint_index) {
                append_projection(
                    result.event_projections,
                    encode_string_sequence(
                        {
                            "segment-endpoint-order-v1",
                            vertex_id,
                            incident->feature_id,
                            std::to_string(
                                endpoint_index),
                        }),
                    polynomial_gcd(
                        passage,
                        endpoint_boundaries[
                            endpoint_index]),
                    {
                        support_event(
                            *incident,
                            "endpoint-order",
                            vertex_id,
                            "merge"),
                    });
            }
        }
    }

    SegmentPairProjectionBundle2 pair_projections =
        derive_segment_pair_projections(
            event_records,
            source,
            result.pullbacks);
    result.event_projections.insert(
        result.event_projections.end(),
        std::make_move_iterator(
            pair_projections.projections.begin()),
        std::make_move_iterator(
            pair_projections.projections.end()));
    result.constant_event_projections =
        std::move(pair_projections.constants);
    result.overlaps =
        std::move(pair_projections.overlaps);

    std::sort(
        result.pullbacks.begin(),
        result.pullbacks.end(),
        [](const ProjectionRecord2& left,
           const ProjectionRecord2& right) {
            return left.projection_id
                < right.projection_id;
        });
    std::sort(
        result.event_projections.begin(),
        result.event_projections.end(),
        [](const ProjectionInput2& left,
           const ProjectionInput2& right) {
            return left.projection_id
                < right.projection_id;
        });
    return result;
}
