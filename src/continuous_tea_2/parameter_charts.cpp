#include "parameter_charts.h"

#include "../exact_algebraic_1.h"

#include <algorithm>
#include <cstdint>
#include <string_view>
#include <utility>

#include <CGAL/CORE/BigRat.h>

namespace {

using Integer = CORE::BigInt;
using Rational = CORE::BigRat;
using Polynomial2 = std::vector<std::vector<Rational>>;

struct DegreeBound {
    int motion_degree;
    int rim_degree;
    std::string identifier;
};

struct RationalMotion {
    Polynomial2 x_numerator;
    Polynomial2 y_numerator;
    Polynomial2 denominator;
};

struct RationalCircle {
    Polynomial2 x_numerator;
    Polynomial2 y_numerator;
    Polynomial2 denominator;
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

std::string ccan_sequence(
    const std::vector<std::string>& values)
{
    std::string payload = u64_record(values.size());
    for (const std::string& value : values) {
        payload += ccan_bytes(value);
    }
    return ccan_node('S', payload);
}

std::string ccan_tagged(
    std::string_view tag,
    const std::string& payload)
{
    return ccan_node(
        'T',
        ccan_bytes(std::string(tag))
            + ccan_bytes(payload));
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

Integer least_common_multiple(
    const Integer& left,
    const Integer& right)
{
    if (left == 0 || right == 0) {
        return 0;
    }
    return left / greatest_common_divisor(left, right) * right;
}

Rational parse_rational(
    const std::string& value,
    std::string_view role)
{
    const std::size_t separator = value.find('/');
    try {
        if (separator == std::string::npos) {
            return Rational(Integer(value));
        }
        if (value.find('/', separator + 1)
            != std::string::npos) {
            throw InvalidAlgebraicPolynomialError(
                std::string(role)
                + " is not an exact rational");
        }
        const Integer numerator(value.substr(0, separator));
        const Integer denominator(
            value.substr(separator + 1));
        if (denominator == 0) {
            throw InvalidAlgebraicPolynomialError(
                std::string(role)
                + " has zero denominator");
        }
        return Rational(numerator, denominator);
    } catch (const InvalidAlgebraicPolynomialError&) {
        throw;
    } catch (const std::exception&) {
        throw InvalidAlgebraicPolynomialError(
            std::string(role)
            + " is not an exact rational");
    }
}

Polynomial2 zero_polynomial()
{
    return {{Rational(0)}};
}

Polynomial2 constant_polynomial(const Rational& value)
{
    return {{value}};
}

Rational coefficient(
    const Polynomial2& polynomial,
    std::size_t motion_degree,
    std::size_t rim_degree)
{
    if (motion_degree >= polynomial.size()
        || rim_degree >= polynomial[motion_degree].size()) {
        return Rational(0);
    }
    return polynomial[motion_degree][rim_degree];
}

Polynomial2 add(
    const Polynomial2& left,
    const Polynomial2& right)
{
    const std::size_t rows =
        std::max(left.size(), right.size());
    std::size_t columns = 0;
    for (const auto& row : left) {
        columns = std::max(columns, row.size());
    }
    for (const auto& row : right) {
        columns = std::max(columns, row.size());
    }
    Polynomial2 result(
        rows,
        std::vector<Rational>(columns, Rational(0)));
    for (std::size_t row = 0; row < rows; ++row) {
        for (std::size_t column = 0;
             column < columns;
             ++column) {
            result[row][column] =
                coefficient(left, row, column)
                + coefficient(right, row, column);
        }
    }
    return result;
}

Polynomial2 scale(
    const Polynomial2& polynomial,
    const Rational& factor)
{
    Polynomial2 result = polynomial;
    for (auto& row : result) {
        for (Rational& value : row) {
            value *= factor;
        }
    }
    return result;
}

Polynomial2 multiply(
    const Polynomial2& left,
    const Polynomial2& right)
{
    std::size_t left_columns = 0;
    std::size_t right_columns = 0;
    for (const auto& row : left) {
        left_columns =
            std::max(left_columns, row.size());
    }
    for (const auto& row : right) {
        right_columns =
            std::max(right_columns, row.size());
    }
    Polynomial2 result(
        left.size() + right.size() - 1,
        std::vector<Rational>(
            left_columns + right_columns - 1,
            Rational(0)));
    for (std::size_t left_row = 0;
         left_row < left.size();
         ++left_row) {
        for (std::size_t left_column = 0;
             left_column < left[left_row].size();
             ++left_column) {
            for (std::size_t right_row = 0;
                 right_row < right.size();
                 ++right_row) {
                for (std::size_t right_column = 0;
                     right_column
                     < right[right_row].size();
                     ++right_column) {
                    result[left_row + right_row]
                          [left_column + right_column] +=
                        left[left_row][left_column]
                        * right[right_row][right_column];
                }
            }
        }
    }
    return result;
}

Polynomial2 motion_polynomial(
    const Rational& constant,
    const Rational& linear)
{
    return {
        {constant},
        {linear},
    };
}

Polynomial2 motion_quadratic(
    const Rational& constant,
    const Rational& linear,
    const Rational& quadratic)
{
    return {
        {constant},
        {linear},
        {quadratic},
    };
}

Polynomial2 rim_quadratic(
    const Rational& constant,
    const Rational& linear,
    const Rational& quadratic)
{
    return {{
        constant,
        linear,
        quadratic,
    }};
}

std::vector<Rational> parse_data(
    const std::vector<std::string>& data,
    std::size_t expected_size,
    std::string_view role)
{
    if (data.size() != expected_size) {
        throw InvalidAlgebraicPolynomialError(
            std::string(role) + " requires "
            + std::to_string(expected_size)
            + " exact coordinates");
    }
    std::vector<Rational> parsed;
    parsed.reserve(data.size());
    for (const std::string& value : data) {
        parsed.push_back(parse_rational(value, role));
    }
    return parsed;
}

RationalMotion segment_motion(
    const std::vector<std::string>& data)
{
    const std::vector<Rational> values =
        parse_data(data, 4, "segment motion");
    return {
        motion_polynomial(
            values[0],
            values[2] - values[0]),
        motion_polynomial(
            values[1],
            values[3] - values[1]),
        constant_polynomial(Rational(1)),
    };
}

std::pair<Polynomial2, Polynomial2>
quarter_circle_numerators(const std::string& chart_id)
{
    const Polynomial2 one_minus_square =
        motion_quadratic(1, 0, -1);
    const Polynomial2 twice_parameter =
        motion_quadratic(0, 2, 0);
    if (chart_id == "center-quarter-0-v1") {
        return {one_minus_square, twice_parameter};
    }
    if (chart_id == "center-quarter-1-v1") {
        return {
            scale(twice_parameter, -1),
            one_minus_square,
        };
    }
    if (chart_id == "center-quarter-2-v1") {
        return {
            scale(one_minus_square, -1),
            scale(twice_parameter, -1),
        };
    }
    if (chart_id == "center-quarter-3-v1") {
        return {
            twice_parameter,
            scale(one_minus_square, -1),
        };
    }
    throw ChartCoverageError(
        "full-circle motion requires a frozen center chart");
}

RationalMotion full_circle_motion(
    const std::vector<std::string>& data,
    const std::string& center_chart)
{
    const std::vector<Rational> values =
        parse_data(data, 3, "full-circle motion");
    const Polynomial2 denominator =
        motion_quadratic(1, 0, 1);
    const auto [unit_x, unit_y] =
        quarter_circle_numerators(center_chart);
    return {
        add(
            scale(denominator, values[0]),
            scale(unit_x, values[2])),
        add(
            scale(denominator, values[1]),
            scale(unit_y, values[2])),
        denominator,
    };
}

RationalCircle rim_circle(
    const std::string& rim_chart)
{
    const Polynomial2 denominator =
        rim_quadratic(1, 0, 1);
    Polynomial2 x_numerator =
        rim_quadratic(1, 0, -1);
    Polynomial2 y_numerator =
        rim_quadratic(0, 2, 0);
    if (rim_chart == "rim-half-1-v1") {
        x_numerator = scale(x_numerator, -1);
        y_numerator = scale(y_numerator, -1);
    } else if (rim_chart != "rim-half-0-v1") {
        throw ChartCoverageError(
            "pullback requires a frozen rim chart");
    }
    return {
        std::move(x_numerator),
        std::move(y_numerator),
        denominator,
    };
}

DegreeBound degree_bound(
    const std::string& motion_kind,
    const std::string& support_kind)
{
    if (motion_kind == "segment"
        && support_kind == "line") {
        return {1, 2, "segment-line-(1,2)-v1"};
    }
    if (motion_kind == "segment"
        && support_kind == "circle") {
        return {2, 4, "segment-circle-(2,4)-v1"};
    }
    if (motion_kind == "full-circle"
        && support_kind == "line") {
        return {2, 2, "full-circle-line-(2,2)-v1"};
    }
    if (motion_kind == "full-circle"
        && support_kind == "circle") {
        return {
            4,
            4,
            "full-circle-circle-(4,4)-v1",
        };
    }
    throw InvalidAlgebraicPolynomialError(
        "unsupported motion/support pullback pair");
}

DegreeBound degree_bound_from_id(
    const std::string& identifier)
{
    for (const std::string& motion :
         {"segment", "full-circle"}) {
        for (const std::string& support :
             {"line", "circle"}) {
            const DegreeBound candidate =
                degree_bound(motion, support);
            if (candidate.identifier == identifier) {
                return candidate;
            }
        }
    }
    throw ProjectionDegreeBoundError(
        "unknown projection degree bound");
}

std::pair<int, int> actual_degree(
    const std::vector<std::vector<Integer>>& rows)
{
    int motion_degree = -1;
    int rim_degree = -1;
    for (std::size_t row = 0; row < rows.size(); ++row) {
        for (std::size_t column = 0;
             column < rows[row].size();
             ++column) {
            if (rows[row][column] != 0) {
                motion_degree =
                    std::max(
                        motion_degree,
                        static_cast<int>(row));
                rim_degree =
                    std::max(
                        rim_degree,
                        static_cast<int>(column));
            }
        }
    }
    return {motion_degree, rim_degree};
}

std::string coefficient_bytes(
    const std::vector<std::vector<std::string>>& rows)
{
    std::vector<std::string> encoded_rows;
    encoded_rows.reserve(rows.size());
    for (const auto& row : rows) {
        encoded_rows.push_back(ccan_sequence(row));
    }
    return ccan_tagged(
        "projection-coefficient-grid-v1",
        ccan_sequence(encoded_rows));
}

ProjectionRecord2 normalized_projection(
    const std::string& projection_id,
    const Polynomial2& polynomial,
    const DegreeBound& bound)
{
    Integer denominator_lcm = 1;
    for (const auto& row : polynomial) {
        for (const Rational& value : row) {
            denominator_lcm = least_common_multiple(
                denominator_lcm,
                CORE::denominator(value));
        }
    }
    std::vector<std::vector<Integer>> integers(
        static_cast<std::size_t>(bound.motion_degree + 1),
        std::vector<Integer>(
            static_cast<std::size_t>(bound.rim_degree + 1),
            Integer(0)));
    for (std::size_t row = 0;
         row < polynomial.size();
         ++row) {
        for (std::size_t column = 0;
             column < polynomial[row].size();
             ++column) {
            if ((static_cast<int>(row)
                    > bound.motion_degree
                 || static_cast<int>(column)
                    > bound.rim_degree)
                && polynomial[row][column] != 0) {
                throw ProjectionDegreeBoundError(
                    "projection exceeds declared degree bound");
            }
            if (static_cast<int>(row)
                    <= bound.motion_degree
                && static_cast<int>(column)
                    <= bound.rim_degree) {
                const Rational scaled =
                    polynomial[row][column]
                    * Rational(denominator_lcm);
                if (CORE::denominator(scaled) != 1) {
                    throw InvalidAlgebraicPolynomialError(
                        "rational coefficient normalization failed");
                }
                integers[row][column] =
                    CORE::numerator(scaled);
            }
        }
    }

    Integer common_divisor = 0;
    for (const auto& row : integers) {
        for (const Integer& value : row) {
            common_divisor =
                greatest_common_divisor(
                    common_divisor,
                    value);
        }
    }
    if (common_divisor != 0) {
        for (auto& row : integers) {
            for (Integer& value : row) {
                value /= common_divisor;
            }
        }
    }

    const auto [motion_degree, rim_degree] =
        actual_degree(integers);

    std::vector<std::vector<std::string>> rows;
    rows.reserve(integers.size());
    for (const auto& integer_row : integers) {
        std::vector<std::string> row;
        row.reserve(integer_row.size());
        for (const Integer& value : integer_row) {
            row.push_back(value.convert_to<std::string>());
        }
        rows.push_back(std::move(row));
    }
    return {
        projection_id,
        rows,
        {},
        motion_degree,
        rim_degree,
        bound.motion_degree,
        bound.rim_degree,
        bound.identifier,
        coefficient_bytes(rows),
    };
}

Polynomial2 parse_grid(
    const std::vector<std::vector<std::string>>& rows)
{
    if (rows.empty()) {
        throw InvalidAlgebraicPolynomialError(
            "projection grid requires at least one row");
    }
    Polynomial2 result;
    result.reserve(rows.size());
    for (const auto& row : rows) {
        if (row.empty()) {
            throw InvalidAlgebraicPolynomialError(
                "projection grid rows cannot be empty");
        }
        std::vector<Rational> parsed;
        parsed.reserve(row.size());
        for (const std::string& value : row) {
            parsed.push_back(
                parse_rational(
                    value,
                    "projection coefficient"));
        }
        result.push_back(std::move(parsed));
    }
    return result;
}

std::string seam_id(
    std::string_view family,
    std::size_t ordinal)
{
    return ccan_tagged(
        "parameter-chart-seam-v1",
        ccan_tagged(
            family,
            std::to_string(ordinal)));
}

bool charts_equal(
    const ParameterChart2& left,
    const ParameterChart2& right)
{
    return left.chart_id == right.chart_id
        && left.family == right.family
        && left.domain_low == right.domain_low
        && left.domain_high == right.domain_high
        && left.orientation == right.orientation
        && left.start_seam_id == right.start_seam_id
        && left.end_seam_id == right.end_seam_id
        && left.owns_start_seam == right.owns_start_seam
        && left.owns_end_seam == right.owns_end_seam;
}

} // namespace

std::vector<ParameterChart2> parameter_charts()
{
    std::vector<ParameterChart2> charts;
    charts.reserve(6);
    for (std::size_t index = 0; index < 4; ++index) {
        charts.push_back(
            {
                "center-quarter-" + std::to_string(index)
                    + "-v1",
                "center-circle",
                "0",
                "1",
                "ccw",
                seam_id("center-circle", index),
                seam_id("center-circle", (index + 1) % 4),
                true,
                false,
            });
    }
    for (std::size_t index = 0; index < 2; ++index) {
        charts.push_back(
            {
                "rim-half-" + std::to_string(index) + "-v1",
                "cutter-rim",
                "-1",
                "1",
                "ccw",
                seam_id("cutter-rim", index),
                seam_id("cutter-rim", (index + 1) % 2),
                true,
                false,
            });
    }
    return charts;
}

VerifiedEventPartition2 verify_chart_coverage(
    const std::vector<ParameterChart2>& charts)
{
    const std::vector<ParameterChart2> expected =
        parameter_charts();
    const bool valid =
        charts.size() == expected.size()
        && std::equal(
            charts.begin(),
            charts.end(),
            expected.begin(),
            charts_equal);

    EventPartitionCertificate2 certificate;
    certificate.build_evidence =
        exact_algebraic_backend_evidence();
    certificate.charts = charts;
    certificate.source_kind = "parameter-chart-coverage-v1";
    if (valid) {
        certificate.seams.reserve(charts.size());
        for (const ParameterChart2& chart : charts) {
            certificate.seams.push_back(
                {
                    chart.start_seam_id,
                    chart.chart_id,
                });
        }
    }
    return {
        valid
            ? ContinuousTeaVerdict::CERTIFIED
            : ContinuousTeaVerdict::UNRESOLVED_DEGENERACY,
        std::move(certificate),
    };
}

ProjectionRecord2 construct_pullback(
    const std::string& motion_kind,
    const std::vector<std::string>& motion_data,
    const std::string& support_kind,
    const std::vector<std::string>& support_data,
    const std::string& cutter_radius,
    const std::string& center_chart,
    const std::string& rim_chart)
{
    const DegreeBound bound =
        degree_bound(motion_kind, support_kind);
    RationalMotion motion;
    if (motion_kind == "segment") {
        if (!center_chart.empty()) {
            throw ChartCoverageError(
                "segment motion does not accept a center chart");
        }
        motion = segment_motion(motion_data);
    } else if (motion_kind == "full-circle") {
        motion =
            full_circle_motion(
                motion_data,
                center_chart);
    } else {
        throw InvalidAlgebraicPolynomialError(
            "unsupported motion kind");
    }

    const RationalCircle rim = rim_circle(rim_chart);
    const Rational cutter =
        parse_rational(cutter_radius, "cutter radius");
    const Polynomial2 denominator =
        multiply(motion.denominator, rim.denominator);
    const Polynomial2 x_numerator = add(
        multiply(
            motion.x_numerator,
            rim.denominator),
        scale(
            multiply(
                motion.denominator,
                rim.x_numerator),
            cutter));
    const Polynomial2 y_numerator = add(
        multiply(
            motion.y_numerator,
            rim.denominator),
        scale(
            multiply(
                motion.denominator,
                rim.y_numerator),
            cutter));

    Polynomial2 pullback;
    if (support_kind == "line") {
        const std::vector<Rational> support =
            parse_data(support_data, 3, "line support");
        pullback = add(
            add(
                scale(x_numerator, support[0]),
                scale(y_numerator, support[1])),
            scale(denominator, support[2]));
    } else if (support_kind == "circle") {
        const std::vector<Rational> support =
            parse_data(
                support_data,
                3,
                "circle support");
        const Polynomial2 relative_x = add(
            x_numerator,
            scale(denominator, -support[0]));
        const Polynomial2 relative_y = add(
            y_numerator,
            scale(denominator, -support[1]));
        pullback = add(
            add(
                multiply(relative_x, relative_x),
                multiply(relative_y, relative_y)),
            scale(
                multiply(denominator, denominator),
                -support[2]));
    } else {
        throw InvalidAlgebraicPolynomialError(
            "unsupported support kind");
    }

    return normalized_projection(
        motion_kind + "-" + support_kind + "-"
            + center_chart + "-" + rim_chart,
        pullback,
        bound);
}

ProjectionRecord2 projection_from_grid(
    const std::string& projection_id,
    const std::vector<std::vector<std::string>>& coefficient_rows,
    const std::string& degree_bound_id)
{
    return normalized_projection(
        projection_id,
        parse_grid(coefficient_rows),
        degree_bound_from_id(degree_bound_id));
}
