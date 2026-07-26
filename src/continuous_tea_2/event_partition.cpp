#include "event_partition.h"

#include "../exact_algebraic_1.h"
#include "event_partition_internal.h"
#include "parameter_charts.h"

#include <algorithm>
#include <sstream>
#include <string_view>
#include <utility>

#include <CGAL/CORE/BigRat.h>
#include <CGAL/Fraction_traits.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>

namespace continuous_tea_2::event_partition_internal {

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
    const Polynomial canonical =
        typename CGAL::Polynomial_traits_d<
            Polynomial>::Canonicalize()(
            polynomial);
    if (CGAL::is_zero(canonical)) {
        return {0};
    }
    std::vector<Integer> result;
    const int degree = CGAL::degree(canonical);
    result.reserve(static_cast<std::size_t>(degree + 1));
    for (int index = 0; index <= degree; ++index) {
        result.push_back(canonical[index]);
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

} // namespace continuous_tea_2::event_partition_internal

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

std::vector<std::string> primitive_rational_polynomial(
    const std::vector<Rational>& coefficients)
{
    using RationalPolynomial =
        CGAL::Polynomial<Rational>;
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
    return coefficient_text(
        primitive_coefficients(numerator));
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
