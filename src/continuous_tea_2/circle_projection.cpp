#include "circle_projection.h"

#include "circle_fibre_state.h"
#include "event_certificate.h"
#include "event_partition.h"
#include "parameter_charts.h"

#include <algorithm>
#include <array>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <CGAL/CORE/BigRat.h>
#include <CGAL/number_utils.h>

using Integer = CORE::BigInt;
using Rational = CORE::BigRat;
using IntegerPolynomial = std::vector<Integer>;

struct CircleVertexSource {
    std::string x;
    std::string y;
    std::vector<std::vector<std::string>> incidents;
};

constexpr std::array<const char*, 4>
    CENTER_CHART_IDS{
        "center-quarter-0-v1",
        "center-quarter-1-v1",
        "center-quarter-2-v1",
        "center-quarter-3-v1",
    };

std::string exact_rational_text(const Epeck::FT& value)
{
    std::ostringstream stream;
    stream << CGAL::exact(value);
    return stream.str();
}

Rational parse_rational_text(const std::string& text)
{
    try {
        const std::size_t separator = text.find('/');
        if (separator == std::string::npos) {
            return Rational(Integer(text));
        }
        if (text.find('/', separator + 1)
                != std::string::npos
            || Integer(text.substr(separator + 1))
                == 0) {
            throw EventPartitionVerificationError(
                "full-circle source contains a malformed rational");
        }
        return Rational(
            Integer(text.substr(0, separator)),
            Integer(text.substr(separator + 1)));
    } catch (const EventSubstrateError&) {
        throw;
    } catch (const std::exception&) {
        throw EventPartitionVerificationError(
            "full-circle source contains a malformed rational");
    }
}

std::string rational_text(const Rational& value)
{
    std::ostringstream stream;
    stream << value;
    return stream.str();
}

std::optional<Rational> rational_square_root(
    const Rational& value)
{
    if (value < 0) {
        return std::nullopt;
    }
    Integer numerator_root;
    Integer denominator_root;
    if (!CGAL::is_square(
            CORE::numerator(value),
            numerator_root)
        || !CGAL::is_square(
            CORE::denominator(value),
            denominator_root)) {
        return std::nullopt;
    }
    return Rational(
        numerator_root,
        denominator_root);
}

std::optional<std::string> rational_coordinate_text(
    const GpsPoint::CoordNT& coordinate)
{
    const Rational base =
        parse_rational_text(
            exact_rational_text(
                coordinate.a0()));
    if (CGAL::is_zero(coordinate.a1())) {
        return rational_text(base);
    }
    const std::optional<Rational> radical =
        rational_square_root(
            parse_rational_text(
                exact_rational_text(
                    coordinate.root())));
    if (!radical.has_value()) {
        return std::nullopt;
    }
    return rational_text(
        base
        + parse_rational_text(
              exact_rational_text(
                  coordinate.a1()))
            * *radical);
}

struct RationalCircleSupport {
    Rational center_x;
    Rational center_y;
    Rational radius_squared;
    std::optional<Rational> radius;
};

RationalCircleSupport rational_circle_support(
    const std::vector<std::string>& coefficients)
{
    if (coefficients.size() != 4) {
        throw EventPartitionVerificationError(
            "full-circle circle support requires four coefficients");
    }
    const Rational quadratic =
        parse_rational_text(coefficients[0]);
    if (quadratic == 0) {
        throw EventPartitionVerificationError(
            "full-circle circle support is degenerate");
    }
    const Rational center_x =
        -parse_rational_text(coefficients[1])
        / (Rational(2) * quadratic);
    const Rational center_y =
        -parse_rational_text(coefficients[2])
        / (Rational(2) * quadratic);
    const Rational radius_squared =
        center_x * center_x
        + center_y * center_y
        - parse_rational_text(coefficients[3])
            / quadratic;
    return {
        center_x,
        center_y,
        radius_squared,
        rational_square_root(radius_squared),
    };
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
    return left / greatest_common_divisor(
                      left,
                      right)
        * right;
}

std::vector<std::string> primitive_factor(
    IntegerPolynomial polynomial,
    bool normalize_sign = true)
{
    while (polynomial.size() > 1
           && polynomial.back() == 0) {
        polynomial.pop_back();
    }
    Integer divisor = 0;
    for (const Integer& coefficient :
         polynomial) {
        divisor = greatest_common_divisor(
            divisor,
            coefficient);
    }
    if (divisor != 0) {
        for (Integer& coefficient :
             polynomial) {
            coefficient /= divisor;
        }
    }
    if (normalize_sign
        && polynomial.back() < 0) {
        for (Integer& coefficient :
             polynomial) {
            coefficient = -coefficient;
        }
    }
    std::vector<std::string> result;
    result.reserve(polynomial.size());
    for (const Integer& coefficient :
         polynomial) {
        result.push_back(
            coefficient.convert_to<std::string>());
    }
    return result;
}

struct RadialPassageFactors {
    std::vector<std::string> zero_set;
    std::vector<std::string> signed_predicate;
};

IntegerPolynomial polynomial_column(
    const ProjectionRecord2& projection,
    std::size_t column)
{
    IntegerPolynomial result;
    result.reserve(
        projection.coefficient_rows.size());
    for (const auto& row :
         projection.coefficient_rows) {
        result.emplace_back(
            column < row.size()
            ? row[column]
            : "0");
    }
    return result;
}

IntegerPolynomial multiply_polynomials(
    const IntegerPolynomial& left,
    const IntegerPolynomial& right)
{
    IntegerPolynomial result(
        left.size() + right.size() - 1,
        Integer(0));
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
    return result;
}

IntegerPolynomial add_polynomials(
    IntegerPolynomial left,
    const IntegerPolynomial& right,
    const Integer& scale)
{
    left.resize(
        std::max(left.size(), right.size()),
        Integer(0));
    for (std::size_t index = 0;
         index < right.size();
         ++index) {
        left[index] += scale * right[index];
    }
    return left;
}

IntegerPolynomial compose_global_chart(
    const IntegerPolynomial& local,
    std::size_t chart)
{
    const IntegerPolynomial affine{
        -Integer(chart),
        Integer(4),
    };
    IntegerPolynomial result{Integer(0)};
    IntegerPolynomial power{Integer(1)};
    for (const Integer& coefficient : local) {
        result = add_polynomials(
            std::move(result),
            power,
            coefficient);
        power = multiply_polynomials(
            power,
            affine);
    }
    while (result.size() > 1
           && result.back() == 0) {
        result.pop_back();
    }
    return result;
}

std::pair<IntegerPolynomial, IntegerPolynomial>
quarter_unit_numerators(std::size_t chart)
{
    const IntegerPolynomial cosine{
        Integer(1),
        Integer(0),
        Integer(-1),
    };
    const IntegerPolynomial sine{
        Integer(0),
        Integer(2),
        Integer(0),
    };
    if (chart == 0) {
        return {cosine, sine};
    }
    if (chart == 1) {
        return {
            add_polynomials(
                {Integer(0)},
                sine,
                Integer(-1)),
            cosine,
        };
    }
    if (chart == 2) {
        return {
            add_polynomials(
                {Integer(0)},
                cosine,
                Integer(-1)),
            add_polynomials(
                {Integer(0)},
                sine,
                Integer(-1)),
        };
    }
    if (chart == 3) {
        return {
            sine,
            add_polynomials(
                {Integer(0)},
                cosine,
                Integer(-1)),
        };
    }
    throw EventPartitionVerificationError(
        "full-circle vertex passage has an unknown chart");
}

std::pair<IntegerPolynomial, IntegerPolynomial>
scaled_motion_center(
    const std::vector<Integer>& scaled,
    std::size_t chart)
{
    const IntegerPolynomial denominator{
        Integer(1),
        Integer(0),
        Integer(1),
    };
    const auto [unit_x, unit_y] =
        quarter_unit_numerators(chart);
    const IntegerPolynomial center_x =
        add_polynomials(
            add_polynomials(
                add_polynomials(
                    {Integer(0)},
                    denominator,
                    scaled[0]),
                unit_x,
                scaled[2]),
            unit_y,
            -scaled[3]);
    const IntegerPolynomial center_y =
        add_polynomials(
            add_polynomials(
                add_polynomials(
                    {Integer(0)},
                    denominator,
                    scaled[1]),
                unit_x,
                scaled[3]),
            unit_y,
            scaled[2]);
    return {center_x, center_y};
}

RadialPassageFactors radial_passage_factor(
    const std::vector<std::string>& motion_data,
    const std::string& radial_distance,
    const std::string& vertex_x,
    const std::string& vertex_y,
    std::size_t chart)
{
    std::vector<Rational> values{
        parse_rational_text(motion_data[0]),
        parse_rational_text(motion_data[1]),
        parse_rational_text(motion_data[2]),
        parse_rational_text(motion_data[3]),
        parse_rational_text(radial_distance),
        parse_rational_text(vertex_x),
        parse_rational_text(vertex_y),
    };
    Integer denominator = 1;
    for (const Rational& value : values) {
        denominator = least_common_multiple(
            denominator,
            CORE::denominator(value));
    }
    std::vector<Integer> scaled;
    scaled.reserve(values.size());
    for (const Rational& value : values) {
        scaled.push_back(
            CORE::numerator(
                value * Rational(denominator)));
    }
    const IntegerPolynomial unit_denominator{
        Integer(1),
        Integer(0),
        Integer(1),
    };
    const auto [center_x, center_y] =
        scaled_motion_center(scaled, chart);
    const IntegerPolynomial dx =
        add_polynomials(
            center_x,
            unit_denominator,
            -scaled[5]);
    const IntegerPolynomial dy =
        add_polynomials(
            center_y,
            unit_denominator,
            -scaled[6]);
    IntegerPolynomial passage =
        add_polynomials(
            multiply_polynomials(dx, dx),
            multiply_polynomials(dy, dy),
            Integer(1));
    passage = add_polynomials(
        std::move(passage),
        multiply_polynomials(
            unit_denominator,
            unit_denominator),
        -scaled[4] * scaled[4]);
    IntegerPolynomial global =
        compose_global_chart(
            passage,
            chart);
    return {
        primitive_factor(global),
        primitive_factor(
            std::move(global),
            false),
    };
}

std::vector<std::string> circle_cap_factor(
    const std::vector<std::string>& motion_data,
    const std::string& cutter_radius,
    const RationalCircleSupport& circle,
    const std::string& cap_chord_ratio,
    std::size_t chart)
{
    std::vector<Rational> values{
        parse_rational_text(motion_data[0]),
        parse_rational_text(motion_data[1]),
        parse_rational_text(motion_data[2]),
        parse_rational_text(motion_data[3]),
        parse_rational_text(cutter_radius),
        circle.center_x,
        circle.center_y,
        circle.radius_squared,
    };
    Integer denominator = 1;
    for (const Rational& value : values) {
        denominator = least_common_multiple(
            denominator,
            CORE::denominator(value));
    }
    std::vector<Integer> scaled;
    scaled.reserve(7);
    for (std::size_t index = 0;
         index < 7;
         ++index) {
        scaled.push_back(
            CORE::numerator(
                values[index]
                * Rational(denominator)));
    }
    const Integer support_radius_squared =
        CORE::numerator(
            values[7]
            * Rational(
                denominator * denominator));
    const IntegerPolynomial unit_denominator{
        Integer(1),
        Integer(0),
        Integer(1),
    };
    const auto [center_x, center_y] =
        scaled_motion_center(scaled, chart);
    const IntegerPolynomial dx =
        add_polynomials(
            center_x,
            unit_denominator,
            -scaled[5]);
    const IntegerPolynomial dy =
        add_polynomials(
            center_y,
            unit_denominator,
            -scaled[6]);
    const IntegerPolynomial distance_squared =
        add_polynomials(
            multiply_polynomials(dx, dx),
            multiply_polynomials(dy, dy),
            Integer(1));
    const IntegerPolynomial denominator_squared =
        multiply_polynomials(
            unit_denominator,
            unit_denominator);
    const Integer cutter_radius_squared =
        scaled[4] * scaled[4];
    const IntegerPolynomial radical_axis =
        add_polynomials(
            distance_squared,
            denominator_squared,
            cutter_radius_squared
                - support_radius_squared);
    const Rational cap =
        parse_rational_text(cap_chord_ratio);
    const Integer cap_numerator =
        CORE::numerator(cap);
    const Integer cap_denominator =
        CORE::denominator(cap);
    const IntegerPolynomial radical_squared =
        multiply_polynomials(
            radical_axis,
            radical_axis);
    IntegerPolynomial first_term =
        multiply_polynomials(
            distance_squared,
            denominator_squared);
    for (Integer& coefficient :
         first_term) {
        coefficient *=
            (Integer(4) * cap_denominator
                - cap_numerator)
            * cutter_radius_squared;
    }
    IntegerPolynomial equality = add_polynomials(
        std::move(first_term),
        radical_squared,
        -cap_denominator);
    return primitive_factor(
        compose_global_chart(
            equality,
            chart));
}

std::vector<std::string> line_tangency_factor(
    const ProjectionRecord2& pullback,
    std::size_t chart)
{
    const IntegerPolynomial constant =
        polynomial_column(pullback, 0);
    const IntegerPolynomial linear =
        polynomial_column(pullback, 1);
    const IntegerPolynomial quadratic =
        polynomial_column(pullback, 2);
    IntegerPolynomial discriminant =
        add_polynomials(
            multiply_polynomials(
                linear,
                linear),
            multiply_polynomials(
                constant,
                quadratic),
            Integer(-4));
    discriminant =
        compose_global_chart(
            discriminant,
            chart);
    return primitive_factor(
        std::move(discriminant));
}

std::vector<std::string> line_cap_factor(
    const std::vector<std::string>& motion_data,
    const std::vector<std::string>& support,
    const std::string& cutter_radius,
    const std::string& cap_chord_ratio,
    std::size_t chart)
{
    std::vector<Rational> values{
        parse_rational_text(motion_data[0]),
        parse_rational_text(motion_data[1]),
        parse_rational_text(motion_data[2]),
        parse_rational_text(motion_data[3]),
        parse_rational_text(cutter_radius),
        parse_rational_text(support[0]),
        parse_rational_text(support[1]),
        parse_rational_text(support[2]),
    };
    Integer denominator = 1;
    for (const Rational& value : values) {
        denominator = least_common_multiple(
            denominator,
            CORE::denominator(value));
    }
    std::vector<Integer> scaled;
    scaled.reserve(values.size());
    for (const Rational& value : values) {
        scaled.push_back(
            CORE::numerator(
                value * Rational(denominator)));
    }
    const IntegerPolynomial unit_denominator{
        Integer(1),
        Integer(0),
        Integer(1),
    };
    const auto [center_x, center_y] =
        scaled_motion_center(scaled, chart);
    IntegerPolynomial signed_distance =
        add_polynomials(
            multiply_polynomials(
                {scaled[5]},
                center_x),
            multiply_polynomials(
                {scaled[6]},
                center_y),
            Integer(1));
    signed_distance = add_polynomials(
        std::move(signed_distance),
        unit_denominator,
        scaled[7] * denominator);
    const Rational cap =
        parse_rational_text(cap_chord_ratio);
    const Integer cap_numerator =
        CORE::numerator(cap);
    const Integer cap_denominator =
        CORE::denominator(cap);
    IntegerPolynomial radial_term =
        multiply_polynomials(
            unit_denominator,
            unit_denominator);
    for (Integer& coefficient : radial_term) {
        coefficient *=
            (Integer(4) * cap_denominator
                - cap_numerator)
            * scaled[4] * scaled[4]
            * (scaled[5] * scaled[5]
               + scaled[6] * scaled[6]);
    }
    IntegerPolynomial equality =
        add_polynomials(
            std::move(radial_term),
            multiply_polynomials(
                signed_distance,
                signed_distance),
            -Integer(4) * cap_denominator);
    return primitive_factor(
        compose_global_chart(
            equality,
            chart));
}

EventPartitionCertificate2
partition_full_circle_boundary_geometry(
    const std::vector<std::string>& motion_data,
    const std::string& cutter_radius,
    const std::string& cap_chord_ratio,
    const std::vector<std::string>& line_sources,
    const std::vector<std::string>& circle_sources)
{
    if (motion_data.size() != 4
        || cutter_radius.empty()
        || cap_chord_ratio.empty()
        || (line_sources.empty()
            && circle_sources.empty())) {
        throw EventPartitionVerificationError(
            "full-circle boundary source is malformed");
    }
    std::vector<ProjectionRecord2> pullbacks;
    std::vector<ProjectionInput2> tangencies;
    std::vector<OverlapInterval2> cap_overlaps;
    std::map<std::string, CircleVertexSource>
        vertices;
    const auto append_vertices =
        [&vertices](
            const std::vector<std::string>& source) {
        for (const std::size_t offset : {4U, 7U}) {
            auto [vertex, inserted] =
                vertices.try_emplace(
                    source[offset],
                    CircleVertexSource{
                        source[offset + 1],
                        source[offset + 2],
                        {},
                    });
            if (!inserted
                && (vertex->second.x
                        != source[offset + 1]
                    || vertex->second.y
                        != source[offset + 2])) {
                throw EventPartitionVerificationError(
                    "full-circle vertex source has inconsistent coordinates");
            }
            vertex->second.incidents.push_back(
                {
                    source[0],
                    source[1],
                    source[2],
                });
        }
    };
    for (const std::string& encoded_source :
         line_sources) {
        const std::vector<std::string> source =
            decode_string_sequence(encoded_source);
        if (source.size() != 10) {
            throw EventPartitionVerificationError(
                "full-circle line support source is malformed");
        }
        const std::vector<std::string> support =
            decode_string_sequence(source[3]);
        if (support.size() != 3) {
            throw EventPartitionVerificationError(
                "full-circle line support requires three coefficients");
        }
        append_vertices(source);
        for (std::size_t center = 0;
             center < CENTER_CHART_IDS.size();
             ++center) {
            for (std::size_t rim = 0;
                 rim < 2;
                 ++rim) {
                const std::string rim_chart =
                    "rim-half-"
                    + std::to_string(rim)
                    + "-v1";
                ProjectionRecord2 pullback =
                    construct_pullback(
                        "full-circle",
                        motion_data,
                        "line",
                        support,
                        cutter_radius,
                        CENTER_CHART_IDS[center],
                        rim_chart);
                pullback.projection_id =
                    encode_string_sequence(
                        {
                            "full-circle-line-pullback-v1",
                            source[0],
                            CENTER_CHART_IDS[center],
                            rim_chart,
                        });
                if (rim == 0) {
                    PartitionEvent2 event{
                        "tangent",
                        source[0],
                        source[1],
                        source[2],
                        {},
                        encode_string_sequence(
                            {
                                "full-circle-line-tangency-v1",
                                source[0],
                                CENTER_CHART_IDS[center],
                            }),
                        "contact",
                    };
                    tangencies.push_back(
                        {
                            event.branch_id,
                            line_tangency_factor(
                                pullback,
                                center),
                            {std::move(event)},
                        });
                    const std::string cap_branch_id =
                        encode_string_sequence(
                            {
                                "full-circle-line-cap-v1",
                                source[0],
                                CENTER_CHART_IDS[center],
                            });
                    std::vector<std::string> cap_factor =
                        line_cap_factor(
                            motion_data,
                            support,
                            cutter_radius,
                            cap_chord_ratio,
                            center);
                    if (cap_factor.size() > 1) {
                        tangencies.push_back(
                            {
                                cap_branch_id,
                                std::move(cap_factor),
                                {
                                    {
                                        "cap-crossing",
                                        source[0],
                                        source[1],
                                        source[2],
                                        {},
                                        cap_branch_id,
                                        "equal-cap",
                                    },
                                },
                            });
                    } else if (
                        cap_factor.front() == "0") {
                        const Rational low(
                            Integer(center),
                            Integer(4));
                        const Rational high(
                            Integer(center + 1),
                            Integer(4));
                        const Rational witness(
                            Integer(2 * center + 1),
                            Integer(8));
                        cap_overlaps.push_back(
                            {
                                "identically-equal-cap-interval",
                                rational_text(low),
                                rational_text(high),
                                CORE::numerator(witness)
                                    .convert_to<std::string>(),
                                CORE::denominator(witness)
                                    .convert_to<std::string>(),
                                "equal-cap",
                                source[0],
                                source[1],
                                source[2],
                                cap_branch_id,
                            });
                    }
                }
                pullbacks.push_back(
                    std::move(pullback));
            }
        }
    }
    const Rational exact_cutter_radius =
        parse_rational_text(cutter_radius);
    const Rational exact_cap_chord_ratio =
        parse_rational_text(cap_chord_ratio);
    if (exact_cutter_radius <= 0
        || exact_cap_chord_ratio <= 0
        || exact_cap_chord_ratio > 4) {
        throw EventPartitionVerificationError(
            "full-circle boundary source has invalid cutter or cap data");
    }
    for (const std::string& encoded_source :
         circle_sources) {
        const std::vector<std::string> source =
            decode_string_sequence(encoded_source);
        if (source.size() != 10) {
            throw EventPartitionVerificationError(
                "full-circle circle support source is malformed");
        }
        const std::vector<std::string> support =
            decode_string_sequence(source[3]);
        append_vertices(source);
        const RationalCircleSupport circle =
            rational_circle_support(support);
        const std::vector<std::string>
            pullback_support{
                rational_text(circle.center_x),
                rational_text(circle.center_y),
                rational_text(
                    circle.radius_squared),
            };
        for (std::size_t center = 0;
             center < CENTER_CHART_IDS.size();
             ++center) {
            for (std::size_t rim = 0;
                 rim < 2;
                 ++rim) {
                const std::string rim_chart =
                    "rim-half-"
                    + std::to_string(rim)
                    + "-v1";
                ProjectionRecord2 pullback =
                    construct_pullback(
                        "full-circle",
                        motion_data,
                        "circle",
                        pullback_support,
                        cutter_radius,
                        CENTER_CHART_IDS[center],
                        rim_chart);
                pullback.projection_id =
                    encode_string_sequence(
                        {
                            "full-circle-circle-pullback-v1",
                            source[0],
                            CENTER_CHART_IDS[center],
                            rim_chart,
                        });
                pullbacks.push_back(
                    std::move(pullback));
            }
            const std::string cap_branch_id =
                encode_string_sequence(
                    {
                        "full-circle-circle-cap-v1",
                        source[0],
                        CENTER_CHART_IDS[center],
                    });
            tangencies.push_back(
                {
                    cap_branch_id,
                    circle_cap_factor(
                        motion_data,
                        cutter_radius,
                        circle,
                        cap_chord_ratio,
                        center),
                    {
                        {
                            "cap-crossing",
                            source[0],
                            source[1],
                            source[2],
                            {},
                            cap_branch_id,
                            "equal-cap",
                        },
                    },
                });
            if (!circle.radius.has_value()) {
                continue;
            }
            const Rational external_radius =
                *circle.radius
                + exact_cutter_radius;
            const Rational internal_radius =
                *circle.radius
                    >= exact_cutter_radius
                ? *circle.radius
                    - exact_cutter_radius
                : exact_cutter_radius
                    - *circle.radius;
            for (const auto& [kind, radius] :
                 std::array<
                     std::pair<const char*, Rational>,
                     2>{{
                     {
                         "external",
                         external_radius,
                     },
                     {
                         "internal",
                         internal_radius,
                     },
                 }}) {
                const std::string branch_id =
                    encode_string_sequence(
                        {
                            "full-circle-circle-tangency-v1",
                            kind,
                            source[0],
                            CENTER_CHART_IDS[center],
                        });
                tangencies.push_back(
                    {
                        branch_id,
                        radial_passage_factor(
                            motion_data,
                            rational_text(radius),
                            rational_text(
                                circle.center_x),
                            rational_text(
                                circle.center_y),
                            center)
                            .zero_set,
                        {
                            {
                                "tangent",
                                source[0],
                                source[1],
                                source[2],
                                {},
                                branch_id,
                                std::string(kind)
                                    + "-contact",
                            },
                        },
                    });
            }
            if (circle.radius_squared
                == exact_cutter_radius
                    * exact_cutter_radius) {
                const std::string branch_id =
                    encode_string_sequence(
                        {
                            "full-circle-circle-coincidence-v1",
                            source[0],
                            CENTER_CHART_IDS[center],
                        });
                tangencies.push_back(
                    {
                        branch_id,
                        radial_passage_factor(
                            motion_data,
                            "0",
                            rational_text(
                                circle.center_x),
                            rational_text(
                                circle.center_y),
                            center)
                            .zero_set,
                        {
                            {
                                "support-overlap",
                                source[0],
                                source[1],
                                source[2],
                                {},
                                branch_id,
                                "coincident-support",
                            },
                        },
                    });
            }
        }
    }
    for (const auto& [vertex_id, vertex] :
         vertices) {
        std::vector<std::string> support_ids;
        support_ids.reserve(
            vertex.incidents.size());
        for (const auto& incident :
             vertex.incidents) {
            support_ids.push_back(incident[1]);
        }
        std::sort(
            support_ids.begin(),
            support_ids.end());
        support_ids.erase(
            std::unique(
                support_ids.begin(),
                support_ids.end()),
            support_ids.end());
        const std::string disposition =
            support_ids.size() > 1
            ? "merge-split"
            : "seam";
        for (std::size_t center = 0;
             center < CENTER_CHART_IDS.size();
             ++center) {
            std::vector<PartitionEvent2> events;
            events.reserve(
                vertex.incidents.size());
            for (const auto& incident :
                 vertex.incidents) {
                events.push_back(
                    {
                        "endpoint-order",
                        incident[0],
                        incident[1],
                        incident[2],
                        vertex_id,
                        encode_string_sequence(
                            {
                                "full-circle-vertex-passage-v1",
                                vertex_id,
                                CENTER_CHART_IDS[center],
                                incident[0],
                            }),
                        disposition,
                    });
            }
            RadialPassageFactors passage =
                radial_passage_factor(
                    motion_data,
                    cutter_radius,
                    vertex.x,
                    vertex.y,
                    center);
            tangencies.push_back(
                {
                    encode_string_sequence(
                        {
                            "full-circle-vertex-passage-v1",
                            vertex_id,
                            CENTER_CHART_IDS[center],
                        }),
                    std::move(passage.zero_set),
                    std::move(events),
                    std::move(
                        passage.signed_predicate),
                });
        }
    }
    EventPartitionCertificate2 certificate =
        partition_projections(tangencies);
    classify_full_circle_endpoint_fibres(
        certificate);
    certificate.overlaps.insert(
        certificate.overlaps.end(),
        cap_overlaps.begin(),
        cap_overlaps.end());
    certificate.projections.insert(
        certificate.projections.end(),
        pullbacks.begin(),
        pullbacks.end());
    for (ParameterCell2& cell :
         certificate.cells) {
        cell.disposition = "sign-invariant";
    }
    return certificate;
}
