#include "event_partition.h"

#include "../exact_algebraic_1.h"
#include "parameter_charts.h"

#include <algorithm>
#include <cstdint>
#include <sstream>
#include <string_view>
#include <utility>

#include <boost/multiprecision/integer.hpp>
#include <CGAL/CORE/BigRat.h>
#include <CGAL/Polynomial_traits_d.h>

namespace {

using Integer = CORE::BigInt;
using Rational = CORE::BigRat;
using Polynomial = ExactAlgebraicKernel1::Polynomial_1;

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
            throw TrimFilterError(
                std::string(role)
                + " is not an exact rational");
        }
        const Integer numerator(
            text.substr(0, separator));
        const Integer denominator(
            text.substr(separator + 1));
        if (denominator == 0) {
            throw TrimFilterError(
                std::string(role)
                + " has zero denominator");
        }
        return Rational(numerator, denominator);
    } catch (const TrimFilterError&) {
        throw;
    } catch (const std::exception&) {
        throw TrimFilterError(
            std::string(role)
            + " is not an exact rational");
    }
}

std::vector<Rational> parse_values(
    const std::vector<std::string>& values,
    std::size_t expected_size,
    std::string_view role)
{
    if (values.size() != expected_size) {
        throw TrimFilterError(
            std::string(role) + " requires "
            + std::to_string(expected_size)
            + " exact values");
    }
    std::vector<Rational> result;
    result.reserve(values.size());
    for (const std::string& value : values) {
        result.push_back(parse_rational(value, role));
    }
    return result;
}

std::string rational_text(const Rational& value)
{
    std::ostringstream stream;
    stream << value;
    return stream.str();
}

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

Polynomial polynomial_from_integers(
    const std::vector<Integer>& coefficients)
{
    return typename CGAL::Polynomial_traits_d<
        Polynomial>::Construct_polynomial()(
        coefficients.begin(),
        coefficients.end());
}

std::vector<Integer> primitive_coefficients(
    const Polynomial& polynomial)
{
    if (CGAL::is_zero(polynomial)) {
        return {0};
    }
    std::vector<Integer> result;
    Integer divisor = 0;
    const int degree = CGAL::degree(polynomial);
    result.reserve(static_cast<std::size_t>(degree + 1));
    for (int index = 0; index <= degree; ++index) {
        result.push_back(polynomial[index]);
        divisor = greatest_common_divisor(
            divisor,
            polynomial[index]);
    }
    for (Integer& value : result) {
        value /= divisor;
    }
    if (result.back() < 0) {
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
    for (const Integer& coefficient : coefficients) {
        result.push_back(
            coefficient.convert_to<std::string>());
    }
    return result;
}

bool projection_is_zero(const ProjectionRecord2& projection)
{
    return std::all_of(
        projection.coefficient_rows.begin(),
        projection.coefficient_rows.end(),
        [](const auto& row) {
            return std::all_of(
                row.begin(),
                row.end(),
                [](const std::string& value) {
                    return Integer(value) == 0;
                });
        });
}

std::vector<std::string> motion_coefficient_gcd(
    const ProjectionRecord2& projection)
{
    std::size_t columns = 0;
    for (const auto& row : projection.coefficient_rows) {
        columns = std::max(columns, row.size());
    }
    Polynomial common(0);
    typename CGAL::Polynomial_traits_d<
        Polynomial>::Gcd_up_to_constant_factor gcd;
    for (std::size_t column = 0;
         column < columns;
         ++column) {
        std::vector<Integer> coefficients;
        coefficients.reserve(
            projection.coefficient_rows.size());
        for (const auto& row : projection.coefficient_rows) {
            coefficients.emplace_back(
                column < row.size()
                    ? row[column]
                    : "0");
        }
        const Polynomial candidate =
            polynomial_from_integers(coefficients);
        if (CGAL::is_zero(candidate)) {
            continue;
        }
        common = CGAL::is_zero(common)
            ? candidate
            : gcd(common, candidate);
    }
    return coefficient_text(
        primitive_coefficients(common));
}

Integer integer_square_root(const Integer& value)
{
    if (value < 0) {
        throw TrimFilterError(
            "negative discriminant has no rational branch");
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

std::vector<Rational> rational_quadratic_roots(
    const std::vector<Rational>& coefficients)
{
    Integer denominator_lcm = 1;
    for (const Rational& coefficient : coefficients) {
        denominator_lcm = least_common_multiple(
            denominator_lcm,
            CORE::denominator(coefficient));
    }
    std::vector<Integer> integers;
    integers.reserve(coefficients.size());
    for (const Rational& coefficient : coefficients) {
        integers.push_back(
            CORE::numerator(
                coefficient
                * Rational(denominator_lcm)));
    }
    Integer divisor = 0;
    for (const Integer& value : integers) {
        divisor = greatest_common_divisor(
            divisor,
            value);
    }
    for (Integer& value : integers) {
        value /= divisor;
    }

    std::vector<Rational> roots;
    if (integers[2] == 0) {
        if (integers[1] == 0) {
            throw TrimFilterError(
                "trim branch equation is constant");
        }
        roots.emplace_back(
            -integers[0],
            integers[1]);
        return roots;
    }
    const Integer discriminant =
        integers[1] * integers[1]
        - Integer(4) * integers[2] * integers[0];
    const Integer square_root =
        integer_square_root(discriminant);
    if (square_root * square_root != discriminant) {
        throw TrimFilterError(
            "trim branch requires a non-rational rim root");
    }
    const Integer denominator = Integer(2) * integers[2];
    roots.emplace_back(
        -integers[1] + square_root,
        denominator);
    if (square_root != 0) {
        roots.emplace_back(
            -integers[1] - square_root,
            denominator);
    }
    std::sort(
        roots.begin(),
        roots.end(),
        std::greater<>());
    return roots;
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

std::vector<std::string> primitive_rational_polynomial(
    const std::vector<Rational>& coefficients)
{
    Integer denominator_lcm = 1;
    for (const Rational& coefficient : coefficients) {
        denominator_lcm = least_common_multiple(
            denominator_lcm,
            CORE::denominator(coefficient));
    }
    std::vector<Integer> integers;
    integers.reserve(coefficients.size());
    for (const Rational& coefficient : coefficients) {
        integers.push_back(
            CORE::numerator(
                coefficient
                * Rational(denominator_lcm)));
    }
    return coefficient_text(
        primitive_coefficients(
            polynomial_from_integers(integers)));
}

} // namespace

EventPartitionCertificate2 partition_projections(
    const std::vector<ProjectionInput2>& projections)
{
    EventPartitionCertificate2 certificate =
        partition_integer_projections(projections);
    certificate.charts = parameter_charts();
    return certificate;
}

EventPartitionCertificate2 partition_pullback_overlap(
    const ProjectionRecord2& projection,
    const std::vector<PartitionEvent2>& events)
{
    if (projection_is_zero(projection)) {
        EventPartitionCertificate2 certificate;
        certificate.build_evidence =
            exact_algebraic_backend_evidence();
        certificate.charts = parameter_charts();
        certificate.projections = {projection};
        certificate.overlaps.push_back(
            {
                "positive-width-motion-overlap",
                "0",
                "1",
                "1",
                "2",
                "support-identical",
            });
        certificate.source_kind =
            "pullback-overlap-v1";
        return certificate;
    }

    const std::vector<std::string> common =
        motion_coefficient_gcd(projection);
    if (common.size() < 2) {
        EventPartitionCertificate2 certificate;
        certificate.build_evidence =
            exact_algebraic_backend_evidence();
        certificate.charts = parameter_charts();
        certificate.projections = {projection};
        certificate.source_kind =
            "pullback-overlap-v1";
        return certificate;
    }
    EventPartitionCertificate2 certificate =
        partition_projections(
            {
                {
                    projection.projection_id
                        + "-coefficient-gcd",
                    common,
                    events,
                },
            });
    certificate.projections = {projection};
    certificate.source_kind =
        "pullback-overlap-v1";
    return certificate;
}

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
    const Rational start_evaluation =
        line[0] * motion[0]
        + line[1] * motion[1] + line[2];
    const Rational end_evaluation =
        line[0] * motion[2]
        + line[1] * motion[3] + line[2];
    if (start_evaluation != end_evaluation) {
        throw TrimFilterError(
            "native line branch solver requires motion parallel to support");
    }
    Rational orientation = 1;
    if (rim_chart == "rim-half-1-v1") {
        orientation = -1;
    } else if (rim_chart != "rim-half-0-v1") {
        throw ChartCoverageError(
            "trim branch requires a frozen rim chart");
    }
    const std::vector<Rational> equation = {
        start_evaluation
            + radius * orientation * line[0],
        radius * orientation * Rational(2) * line[1],
        start_evaluation
            - radius * orientation * line[0],
    };
    const std::vector<Rational> roots =
        rational_quadratic_roots(equation);
    const Rational trim_dx = end[0] - start[0];
    const Rational trim_dy = end[1] - start[1];
    if (trim_dx == 0 && trim_dy == 0) {
        throw TrimFilterError(
            "trim endpoints must be distinct");
    }
    const Rational velocity_x =
        motion[2] - motion[0];
    const Rational velocity_y =
        motion[3] - motion[1];
    if (velocity_x * trim_dy
        != velocity_y * trim_dx) {
        throw TrimFilterError(
            "motion is not parallel to the closed trim");
    }

    std::vector<TrimmedLineBranch2> branches;
    for (const Rational& root : roots) {
        if (root < -1 || root > 1) {
            continue;
        }
        const Rational denominator =
            Rational(1) + root * root;
        const Rational rim_x =
            orientation
            * (Rational(1) - root * root)
            / denominator;
        const Rational rim_y =
            orientation * Rational(2) * root
            / denominator;
        const Rational point_x =
            motion[0] + radius * rim_x;
        const Rational point_y =
            motion[1] + radius * rim_y;
        if ((point_x - start[0]) * trim_dy
            != (point_y - start[1]) * trim_dx) {
            throw TrimFilterError(
                "rim root does not lie on the exact trim support");
        }
        const bool use_x = trim_dx != 0;
        const Rational lambda_zero =
            use_x
            ? (point_x - start[0]) / trim_dx
            : (point_y - start[1]) / trim_dy;
        const Rational lambda_velocity =
            use_x
            ? velocity_x / trim_dx
            : velocity_y / trim_dy;
        const auto [domain_low, domain_high] =
            clipped_unit_interval(
                lambda_zero,
                lambda_velocity);
        const std::string rim_parameter =
            rational_text(root);
        const std::string feature_id = stable_record(
            "trimmed-line-feature-v1",
            line_support);
        const std::string trim_id = stable_record(
            "closed-line-trim-v1",
            {
                rational_text(start[0]),
                rational_text(start[1]),
                rational_text(end[0]),
                rational_text(end[1]),
            });
        branches.push_back(
            {
                rim_parameter,
                rational_text(domain_low),
                rational_text(domain_high),
                "accepted",
                true,
                feature_id,
                trim_id,
                stable_record(
                    "trimmed-line-branch-v1",
                    {
                        feature_id,
                        trim_id,
                        rim_chart,
                        rim_parameter,
                    }),
            });
    }
    return branches;
}

ProjectedRegularizationVertex2 project_regularization_vertex(
    const Stock2& stock,
    std::size_t first_index,
    std::size_t second_index,
    const std::string& vertex_id)
{
    const std::vector<BoundaryFeatureRecord2> records =
        extract_boundary_records(stock);
    if (first_index >= records.size()
        || second_index >= records.size()) {
        throw BoundaryFeatureIndexError(
            "boundary feature index is out of range");
    }
    const BoundaryFeatureRecord2& first =
        records[first_index];
    const BoundaryFeatureRecord2& second =
        records[second_index];
    if (first.support_kind != "circle"
        || second.support_kind != "circle") {
        throw TrimFilterError(
            "regularization projection requires two circle supports");
    }
    const std::vector<Rational> first_circle =
        parse_values(
            first.primitive_coefficients,
            4,
            "first circle");
    const std::vector<Rational> second_circle =
        parse_values(
            second.primitive_coefficients,
            4,
            "second circle");
    const Rational first_x =
        first_circle[1] / first_circle[0];
    const Rational first_y =
        first_circle[2] / first_circle[0];
    const Rational first_constant =
        first_circle[3] / first_circle[0];
    const Rational second_x =
        second_circle[1] / second_circle[0];
    const Rational second_y =
        second_circle[2] / second_circle[0];
    const Rational second_constant =
        second_circle[3] / second_circle[0];
    const Rational radical_x = first_x - second_x;
    const Rational radical_y = first_y - second_y;
    const Rational radical_constant =
        first_constant - second_constant;
    if (radical_x == 0 && radical_y == 0) {
        throw TrimFilterError(
            "coincident circle supports have no isolated vertex");
    }

    Rational substitution_linear;
    Rational substitution_constant;
    Rational quadratic;
    Rational linear;
    Rational constant;
    if (radical_x != 0) {
        substitution_linear = -radical_y / radical_x;
        substitution_constant =
            -radical_constant / radical_x;
        quadratic =
            substitution_linear * substitution_linear + 1;
        linear =
            Rational(2) * substitution_linear
                * substitution_constant
            + first_x * substitution_linear + first_y;
        constant =
            substitution_constant * substitution_constant
            + first_x * substitution_constant
            + first_constant;
    } else {
        substitution_linear = -radical_x / radical_y;
        substitution_constant =
            -radical_constant / radical_y;
        quadratic =
            substitution_linear * substitution_linear + 1;
        linear =
            Rational(2) * substitution_linear
                * substitution_constant
            + first_y * substitution_linear + first_x;
        constant =
            substitution_constant * substitution_constant
            + first_y * substitution_constant
            + first_constant;
    }
    const std::vector<std::string> depressed =
        primitive_rational_polynomial(
            {
                constant
                    - linear * linear
                        / (Rational(4) * quadratic),
                Rational(0),
                quadratic,
            });

    std::vector<BoundaryEvent2> vertices;
    for (const BoundaryEvent2& event :
         classify_boundary_pair(
             stock,
             first_index,
             second_index)) {
        if (event.kind == "vertex") {
            vertices.push_back(event);
        }
    }
    const auto found = std::find_if(
        vertices.begin(),
        vertices.end(),
        [&vertex_id](const BoundaryEvent2& event) {
            return event.vertex_id == vertex_id;
        });
    if (found == vertices.end()) {
        throw TrimFilterError(
            "requested vertex is not on both closed trims");
    }
    if (found != std::prev(vertices.end())) {
        throw TrimFilterError(
            "requested conjugate lies outside the canonical positive chart");
    }

    const EventPartitionCertificate2 partition =
        partition_projections(
            {
                {
                    "regularization-vertex",
                    depressed,
                    {},
                },
            });
    if (partition.roots.size() != 1) {
        throw TrimFilterError(
            "regularization projection did not isolate one positive root");
    }
    AlgebraicRootRecord2 root = partition.roots.front();
    root.interval_low = rational_text(
        parse_rational(
            root.interval_low,
            "root interval")
        / Rational(2));
    root.interval_high = rational_text(
        (parse_rational(
             root.interval_high,
             "root interval")
         + Rational(1))
        / Rational(2));
    return {
        std::move(root),
        vertex_id,
        "accepted",
        "accepted",
        "rejected-vertex-identity",
    };
}
