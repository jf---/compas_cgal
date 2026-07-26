#include "circle_oracle.h"

#include "../exact_algebraic_1.h"
#include "../stock_2.h"
#include "boundary_events.h"
#include "cap_partition.h"
#include "event_certificate.h"
#include "event_partition.h"
#include "parameter_charts.h"
#include "sha256.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <CGAL/CORE/BigRat.h>
#include <CGAL/number_utils.h>

namespace {

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

struct FibreBlock {
    std::size_t ccw_order;
    std::string global_fibre_id;
    std::vector<EventTraceEvent2> events;
};

const ParameterChart2& require_center_chart(
    const EventPartitionCertificate2& partition,
    const std::string& chart_id)
{
    const auto chart = std::find_if(
        partition.charts.begin(),
        partition.charts.end(),
        [&chart_id](const ParameterChart2& candidate) {
            return candidate.chart_id == chart_id;
        });
    if (chart == partition.charts.end()) {
        throw EventTraceVerificationError(
            "full-circle event order requires all four center charts");
    }
    if (chart->family != "center-circle"
        || chart->domain_low != "0"
        || chart->domain_high != "1"
        || chart->orientation != "ccw"
        || !chart->owns_start_seam
        || chart->owns_end_seam) {
        throw EventTraceVerificationError(
            "full-circle center chart violates the frozen topology");
    }
    return *chart;
}

void require_owned_seam(
    const EventPartitionCertificate2& partition,
    const ParameterChart2& chart)
{
    const std::size_t owned_count =
        static_cast<std::size_t>(
            std::count_if(
                partition.seams.begin(),
                partition.seams.end(),
                [&chart](const ChartSeam2& seam) {
                    return seam.owner_chart_id
                            == chart.chart_id
                        && seam.seam_id
                            == chart.start_seam_id;
                }));
    if (owned_count != 1) {
        throw EventTraceVerificationError(
            "full-circle center seam must have one exact owner");
    }
}

std::string require_circle_topology(
    const VerifiedEventPartition2& verified_partition)
{
    if (verified_partition.verdict
        != ContinuousTeaVerdict::CERTIFIED) {
        throw EventTraceVerificationError(
            "full-circle event order requires a verified partition");
    }
    std::vector<std::string> seam_ids;
    seam_ids.reserve(CENTER_CHART_IDS.size());
    for (const char* chart_id : CENTER_CHART_IDS) {
        const ParameterChart2& chart =
            require_center_chart(
                verified_partition.partition,
                chart_id);
        require_owned_seam(
            verified_partition.partition,
            chart);
        if (std::find(
                seam_ids.begin(),
                seam_ids.end(),
                chart.start_seam_id)
            != seam_ids.end()) {
            throw EventTraceVerificationError(
                "full-circle center seams must be distinct");
        }
        seam_ids.push_back(chart.start_seam_id);
    }
    return seam_ids.front();
}

FibreBlock& event_fibre(
    std::vector<FibreBlock>& fibres,
    std::unordered_map<std::string, std::size_t>&
        fibre_indices,
    const EventTraceEvent2& event)
{
    const auto existing =
        fibre_indices.find(
            event.global_fibre_id);
    if (existing != fibre_indices.end()) {
        FibreBlock& fibre =
            fibres[existing->second];
        if (fibre.ccw_order
            != event.motion_order) {
            throw EventTraceVerificationError(
                "one full-circle fibre has conflicting CCW orders");
        }
        return fibre;
    }
    fibres.push_back(
        FibreBlock{
            event.motion_order,
            event.global_fibre_id,
            {},
        });
    fibre_indices.emplace(
        event.global_fibre_id,
        fibres.size() - 1);
    return fibres.back();
}

bool event_identity_less(
    const EventTraceEvent2& left,
    const EventTraceEvent2& right)
{
    return left.canonical_id < right.canonical_id;
}

bool has_anchor_event(
    const FibreBlock& fibre,
    const std::string& anchor_seam_id)
{
    return std::any_of(
        fibre.events.begin(),
        fibre.events.end(),
        [&anchor_seam_id](
            const EventTraceEvent2& event) {
            return event.kind == "seam"
                && std::find(
                       event.feature_ids.begin(),
                       event.feature_ids.end(),
                       anchor_seam_id)
                    != event.feature_ids.end();
        });
}

std::string binary64_identity(double value)
{
    return std::to_string(
        std::bit_cast<std::uint64_t>(value));
}

void require_finite(
    double value,
    const char* role)
{
    if (!std::isfinite(value)) {
        throw IncompleteFullCircleOracleError(
            std::string(role)
            + " must be finite binary64");
    }
}

std::string motion_identity(
    double center_x,
    double center_y,
    double phase_dx,
    double phase_dy,
    bool clockwise)
{
    return encode_canonical_record(
        "full-circle-motion-binary64-v1",
        {
            binary64_identity(center_x),
            binary64_identity(center_y),
            binary64_identity(phase_dx),
            binary64_identity(phase_dy),
            clockwise ? "clockwise" : "counterclockwise",
        });
}

const ParameterChart2& center_chart(
    const std::vector<ParameterChart2>& charts,
    std::size_t index)
{
    const std::string chart_id =
        CENTER_CHART_IDS.at(index);
    const auto found = std::find_if(
        charts.begin(),
        charts.end(),
        [&chart_id](
            const ParameterChart2& chart) {
            return chart.chart_id == chart_id;
        });
    if (found == charts.end()) {
        throw IncompleteFullCircleOracleError(
            "full-circle scaffold lacks a center chart");
    }
    return *found;
}

std::string stock_identity(
    const std::vector<BoundaryFeatureRecord2>& records)
{
    std::vector<std::string> feature_ids;
    feature_ids.reserve(records.size());
    for (const BoundaryFeatureRecord2& record : records) {
        feature_ids.push_back(record.feature_id);
    }
    return sha256_bytes(
        encode_string_sequence(feature_ids));
}

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
    IntegerPolynomial polynomial)
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
    if (polynomial.back() < 0) {
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

std::vector<std::string> radial_passage_factor(
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
    const auto [unit_x, unit_y] =
        quarter_unit_numerators(chart);
    const IntegerPolynomial dx =
        add_polynomials(
            add_polynomials(
                {Integer(0)},
                unit_denominator,
                scaled[0] - scaled[4]),
            unit_x,
            scaled[2]);
    const IntegerPolynomial dy =
        add_polynomials(
            add_polynomials(
                {Integer(0)},
                unit_denominator,
                scaled[1] - scaled[5]),
            unit_y,
            scaled[2]);
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
        -scaled[3] * scaled[3]);
    return primitive_factor(
        compose_global_chart(
            passage,
            chart));
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
    scaled.reserve(6);
    for (std::size_t index = 0;
         index < 6;
         ++index) {
        scaled.push_back(
            CORE::numerator(
                values[index]
                * Rational(denominator)));
    }
    const Integer support_radius_squared =
        CORE::numerator(
            values[6]
            * Rational(
                denominator * denominator));
    const IntegerPolynomial unit_denominator{
        Integer(1),
        Integer(0),
        Integer(1),
    };
    const auto [unit_x, unit_y] =
        quarter_unit_numerators(chart);
    const IntegerPolynomial dx =
        add_polynomials(
            add_polynomials(
                {Integer(0)},
                unit_denominator,
                scaled[0] - scaled[4]),
            unit_x,
            scaled[2]);
    const IntegerPolynomial dy =
        add_polynomials(
            add_polynomials(
                {Integer(0)},
                unit_denominator,
                scaled[1] - scaled[5]),
            unit_y,
            scaled[2]);
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
        scaled[3] * scaled[3];
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

GpsPolygon circle_disk(
    const EPoint& center,
    const Epeck::FT& radius)
{
    const ECircle circle(
        center,
        radius * radius,
        CGAL::COUNTERCLOCKWISE);
    const GpsPoint minimum(
        center.x() - radius,
        center.y());
    const GpsPoint maximum(
        center.x() + radius,
        center.y());
    GpsPolygon polygon;
    polygon.push_back(
        GpsXCurve(
            circle,
            minimum,
            maximum,
            CGAL::COUNTERCLOCKWISE));
    polygon.push_back(
        GpsXCurve(
            circle,
            maximum,
            minimum,
            CGAL::COUNTERCLOCKWISE));
    return polygon;
}

std::optional<std::string> uniform_disposition(
    const Stock2& stock,
    const std::vector<BoundaryFeatureRecord2>& records,
    double center_x,
    double center_y,
    double phase_dx,
    double phase_dy,
    double tool_radius)
{
    if (records.empty()) {
        return "clear";
    }
    const Epeck::FT phase_x(phase_dx);
    const Epeck::FT phase_y(phase_dy);
    Epeck::FT guide_radius;
    if (CGAL::is_zero(phase_y)) {
        guide_radius = CGAL::abs(phase_x);
    } else if (CGAL::is_zero(phase_x)) {
        guide_radius = CGAL::abs(phase_y);
    } else {
        return std::nullopt;
    }
    Gps outer_disk;
    outer_disk.insert(
        circle_disk(
            EPoint(center_x, center_y),
            guide_radius + Epeck::FT(tool_radius)));
    outer_disk.difference(stock.set());
    if (outer_disk.is_empty()) {
        return "material";
    }
    return std::nullopt;
}

bool has_material_rational_chart_witness(
    const Stock2& stock,
    const Epeck::FT& center_x,
    const Epeck::FT& center_y,
    const Epeck::FT& guide_radius,
    const Epeck::FT& tool_radius)
{
    // Quarter-chart t=1/2 gives the exact unit direction (3/5, 4/5).
    const Epeck::FT chart_denominator(5);
    const Epeck::FT witness_x =
        center_x
        + guide_radius
            * Epeck::FT(3)
            / chart_denominator;
    const Epeck::FT witness_y =
        center_y
        + guide_radius
            * Epeck::FT(4)
            / chart_denominator;
    Gps cutter_disk;
    cutter_disk.insert(
        circle_disk(
            EPoint(witness_x, witness_y),
            tool_radius));
    cutter_disk.difference(stock.set());
    return cutter_disk.is_empty();
}

} // namespace

EventPartitionCertificate2
construct_full_circle_uniform_partition(
    const std::string& stock_identity_value,
    const std::string& disposition)
{
    if (stock_identity_value.empty()
        || (disposition != "clear"
            && disposition != "material")) {
        throw EventPartitionVerificationError(
            "full-circle uniform source is malformed");
    }
    const std::vector<ParameterChart2> charts =
        parameter_charts();
    const std::array<std::array<const char*, 2>, 4>
        seam_factors{{
            {{"0", "1"}},
            {{"-1", "4"}},
            {{"-1", "2"}},
            {{"-3", "4"}},
        }};
    std::vector<ProjectionInput2> projections;
    projections.reserve(CENTER_CHART_IDS.size());
    for (std::size_t index = 0;
         index < CENTER_CHART_IDS.size();
         ++index) {
        const ParameterChart2& chart =
            center_chart(charts, index);
        PartitionEvent2 event{
            "seam",
            chart.start_seam_id,
            "full-circle-center-support-v1",
            "full-circle-center-domain-v1",
            chart.start_seam_id,
            chart.chart_id,
            "owned",
        };
        event.original_equations_rechecked = true;
        event.orientation_rechecked = true;
        event.trim_disposition = "accepted";
        projections.push_back(
            {
                "full-circle-seam-"
                    + std::to_string(index)
                    + "-v1",
                {
                    seam_factors[index][0],
                    seam_factors[index][1],
                },
                {std::move(event)},
            });
    }
    EventPartitionCertificate2 certificate =
        partition_projections(projections);
    if (certificate.roots.size() != 4
        || certificate.fibres.size() != 4
        || certificate.cells.size() != 4) {
        throw EventPartitionVerificationError(
            "full-circle uniform partition lacks four cyclic strata");
    }
    for (ParameterCell2& cell : certificate.cells) {
        cell.disposition =
            disposition == "clear"
            ? "below-cap"
            : "cap-exceeds";
    }
    certificate.seams.reserve(charts.size());
    for (const ParameterChart2& chart : charts) {
        certificate.seams.push_back(
            {
                chart.start_seam_id,
                chart.chart_id,
            });
    }
    certificate.source_kind =
        "full-circle-uniform-v1";
    certificate.source_payload =
        encode_string_sequence(
            {
                stock_identity_value,
                disposition,
            });
    finalize_event_partition(certificate);
    return certificate;
}

EventPartitionCertificate2
construct_full_circle_boundary_pullback_partition(
    const std::string& stock_identity_value,
    const std::vector<std::string>& motion_data,
    const std::string& cutter_radius,
    const std::string& cap_chord_ratio,
    const std::vector<std::string>& line_sources,
    const std::vector<std::string>& circle_sources)
{
    if (stock_identity_value.empty()
        || motion_data.size() != 3
        || cutter_radius.empty()
        || cap_chord_ratio.empty()
        || (line_sources.empty()
            && circle_sources.empty())) {
        throw EventPartitionVerificationError(
            "full-circle boundary source is malformed");
    }
    std::vector<ProjectionRecord2> pullbacks;
    std::vector<ProjectionInput2> tangencies;
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
                            center),
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
                            center),
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
            tangencies.push_back(
                {
                    encode_string_sequence(
                        {
                            "full-circle-vertex-passage-v1",
                            vertex_id,
                            CENTER_CHART_IDS[center],
                        }),
                    radial_passage_factor(
                        motion_data,
                        cutter_radius,
                        vertex.x,
                        vertex.y,
                        center),
                    std::move(events),
                });
        }
    }
    EventPartitionCertificate2 certificate =
        partition_projections(tangencies);
    certificate.projections.insert(
        certificate.projections.end(),
        pullbacks.begin(),
        pullbacks.end());
    for (ParameterCell2& cell :
         certificate.cells) {
        cell.disposition = "sign-invariant";
    }
    certificate.seams.reserve(
        certificate.charts.size());
    for (const ParameterChart2& chart :
         certificate.charts) {
        certificate.seams.push_back(
            {
                chart.start_seam_id,
                chart.chart_id,
            });
    }
    certificate.source_kind =
        "full-circle-boundary-pullbacks-v1";
    certificate.source_payload =
        encode_string_sequence(
            {
                stock_identity_value,
                encode_string_sequence(motion_data),
                cutter_radius,
                cap_chord_ratio,
                encode_string_sequence(line_sources),
                encode_string_sequence(circle_sources),
            });
    finalize_event_partition(certificate);
    return certificate;
}

std::vector<EventTraceEvent2> order_full_circle_events(
    const VerifiedEventPartition2& verified_partition,
    bool clockwise,
    std::vector<EventTraceEvent2> events)
{
    const std::string anchor_seam_id =
        require_circle_topology(verified_partition);
    if (events.empty()) {
        throw EventTraceVerificationError(
            "full-circle event order requires trace events");
    }

    std::vector<FibreBlock> fibres;
    std::unordered_map<std::string, std::size_t>
        fibre_indices;
    fibre_indices.reserve(events.size());
    std::unordered_set<std::string> root_ids;
    root_ids.reserve(
        verified_partition.partition.roots.size());
    for (const AlgebraicRootRecord2& root :
         verified_partition.partition.roots) {
        root_ids.insert(root.root_id);
    }
    for (EventTraceEvent2& event : events) {
        if (root_ids.find(event.root_id)
            == root_ids.end()) {
            throw EventTraceVerificationError(
                "full-circle trace event root is absent from the verified partition");
        }
        event_fibre(
            fibres,
            fibre_indices,
            event)
            .events.push_back(std::move(event));
    }
    std::sort(
        fibres.begin(),
        fibres.end(),
        [](const FibreBlock& left,
           const FibreBlock& right) {
            if (left.ccw_order
                != right.ccw_order) {
                return left.ccw_order
                    < right.ccw_order;
            }
            return left.global_fibre_id
                < right.global_fibre_id;
        });
    for (std::size_t index = 1;
         index < fibres.size();
         ++index) {
        if (fibres[index - 1].ccw_order
            == fibres[index].ccw_order) {
            throw EventTraceVerificationError(
                "full-circle CCW order assigns two fibres one rank");
        }
    }
    if (fibres.front().ccw_order != 0
        || !has_anchor_event(
            fibres.front(),
            anchor_seam_id)) {
        throw EventTraceVerificationError(
            "full-circle order must start at the owned phase seam");
    }
    for (FibreBlock& fibre : fibres) {
        std::sort(
            fibre.events.begin(),
            fibre.events.end(),
            event_identity_less);
        for (std::size_t index = 1;
             index < fibre.events.size();
             ++index) {
            if (fibre.events[index - 1].canonical_id
                == fibre.events[index].canonical_id) {
                throw EventTraceVerificationError(
                    "full-circle fibre contains a duplicate event identity");
            }
        }
    }
    if (clockwise && fibres.size() > 1) {
        std::reverse(
            fibres.begin() + 1,
            fibres.end());
    }

    std::vector<EventTraceEvent2> ordered;
    ordered.reserve(events.size());
    for (std::size_t order = 0;
         order < fibres.size();
         ++order) {
        for (EventTraceEvent2& event :
             fibres[order].events) {
            event.motion_order = order;
            ordered.push_back(std::move(event));
        }
    }
    return ordered;
}

std::pair<std::string, EventTrace2>
audit_full_circle_tea_event_exact(
    const Stock2& stock,
    double center_x,
    double center_y,
    double phase_dx,
    double phase_dy,
    bool clockwise,
    double tool_radius,
    double cap_chord_ratio)
{
    require_finite(center_x, "circle center x");
    require_finite(center_y, "circle center y");
    require_finite(phase_dx, "circle phase dx");
    require_finite(phase_dy, "circle phase dy");
    require_finite(tool_radius, "tool radius");
    require_finite(cap_chord_ratio, "cap chord ratio");
    if (phase_dx == 0.0 && phase_dy == 0.0) {
        throw IncompleteFullCircleOracleError(
            "full-circle phase vector must be nonzero");
    }
    if (!(tool_radius > 0.0)) {
        throw IncompleteFullCircleOracleError(
            "tool radius must be positive");
    }
    if (!(cap_chord_ratio > 0.0
          && cap_chord_ratio <= 4.0)) {
        throw IncompleteFullCircleOracleError(
            "cap chord ratio must lie in (0, 4]");
    }

    const std::vector<BoundaryFeatureRecord2>
        boundary_records =
            extract_boundary_records(stock);
    const std::optional<std::string> uniform =
        uniform_disposition(
            stock,
            boundary_records,
            center_x,
            center_y,
            phase_dx,
            phase_dy,
            tool_radius);
    if (uniform.has_value()) {
        const EventPartitionCertificate2 partition =
            construct_full_circle_uniform_partition(
                stock_identity(boundary_records),
                *uniform);
        const VerifiedEventPartition2 verified =
            verify_event_partition(partition);
        if (verified.verdict
            != ContinuousTeaVerdict::CERTIFIED) {
            throw IncompleteFullCircleOracleError(
                "full-circle uniform partition failed replay");
        }
        std::vector<EventTraceEvent2> events;
        events.reserve(
            verified.partition.fibres.size());
        for (std::size_t index = 0;
             index < verified.partition.fibres.size();
             ++index) {
            const EventFibre2& fibre =
                verified.partition.fibres[index];
            if (fibre.events.size() != 1) {
                throw IncompleteFullCircleOracleError(
                    "full-circle uniform seam fibre is not singular");
            }
            const auto root = std::find_if(
                verified.partition.roots.begin(),
                verified.partition.roots.end(),
                [&fibre](
                    const AlgebraicRootRecord2& candidate) {
                    return candidate.root_id
                        == fibre.root_id;
                });
            if (root
                == verified.partition.roots.end()) {
                throw IncompleteFullCircleOracleError(
                    "full-circle uniform seam fibre lacks its exact root");
            }
            const PartitionEvent2& event =
                fibre.events.front();
            events.push_back(
                make_event_trace_event(
                    root->root_id,
                    encode_canonical_record(
                        "full-circle-global-fibre-v1",
                        {
                            root->root_id,
                            event.branch_id,
                        }),
                    "seam",
                    {event.feature_id},
                    {event.branch_id},
                    root->multiplicity,
                    event.disposition,
                    index));
        }
        events = order_full_circle_events(
            verified,
            clockwise,
            std::move(events));
        const ContinuousTeaVerdict verdict =
            *uniform == "clear"
            ? ContinuousTeaVerdict::CERTIFIED
            : ContinuousTeaVerdict::CAP_EXCEEDED;
        EventTrace2 trace =
            build_event_trace(
                verified.partition,
                "full-circle-four-chart-v1",
                motion_identity(
                    center_x,
                    center_y,
                    phase_dx,
                    phase_dy,
                    clockwise),
                encode_canonical_record(
                    "cap-chord-ratio-binary64-v1",
                    {
                        binary64_identity(
                            cap_chord_ratio),
                    }),
                verdict,
                *uniform,
                "full-circle-uniform-event-exact-v1",
                std::move(events));
        return {
            *uniform == "clear"
                ? "certified"
                : "cap_exceeded",
            std::move(trace),
        };
    }

    std::vector<std::string> line_sources;
    std::vector<std::string> circle_sources;
    bool rational_boundary_vertices = true;
    for (const BoundaryFeatureRecord2& record :
         boundary_records) {
        if (record.support_kind == "line") {
            const std::optional<std::string> source_x =
                rational_coordinate_text(
                    record.curve.source().x());
            const std::optional<std::string> source_y =
                rational_coordinate_text(
                    record.curve.source().y());
            const std::optional<std::string> target_x =
                rational_coordinate_text(
                    record.curve.target().x());
            const std::optional<std::string> target_y =
                rational_coordinate_text(
                    record.curve.target().y());
            if (!source_x.has_value()
                || !source_y.has_value()
                || !target_x.has_value()
                || !target_y.has_value()) {
                rational_boundary_vertices = false;
                break;
            }
            line_sources.push_back(
                encode_string_sequence(
                    {
                        record.feature_id,
                        record.support_id,
                        record.trim_predicate,
                        encode_string_sequence(
                            record.primitive_coefficients),
                        record.source_vertex_id,
                        *source_x,
                        *source_y,
                        record.target_vertex_id,
                        *target_x,
                        *target_y,
                    }));
        } else if (record.support_kind == "circle") {
            const std::optional<std::string> source_x =
                rational_coordinate_text(
                    record.curve.source().x());
            const std::optional<std::string> source_y =
                rational_coordinate_text(
                    record.curve.source().y());
            const std::optional<std::string> target_x =
                rational_coordinate_text(
                    record.curve.target().x());
            const std::optional<std::string> target_y =
                rational_coordinate_text(
                    record.curve.target().y());
            if (!source_x.has_value()
                || !source_y.has_value()
                || !target_x.has_value()
                || !target_y.has_value()) {
                rational_boundary_vertices = false;
                break;
            }
            circle_sources.push_back(
                encode_string_sequence(
                    {
                        record.feature_id,
                        record.support_id,
                        record.trim_predicate,
                        encode_string_sequence(
                            record.primitive_coefficients),
                        record.source_vertex_id,
                        *source_x,
                        *source_y,
                        record.target_vertex_id,
                        *target_x,
                        *target_y,
                    }));
        }
    }
    const Epeck::FT exact_phase_x(phase_dx);
    const Epeck::FT exact_phase_y(phase_dy);
    if ((!line_sources.empty()
            || !circle_sources.empty())
        && rational_boundary_vertices
        && (CGAL::is_zero(exact_phase_x)
            || CGAL::is_zero(exact_phase_y))) {
        const Epeck::FT guide_radius =
            CGAL::is_zero(exact_phase_y)
            ? CGAL::abs(exact_phase_x)
            : CGAL::abs(exact_phase_y);
        EventPartitionCertificate2 partition =
            construct_full_circle_boundary_pullback_partition(
                stock_identity(boundary_records),
                {
                    exact_rational_text(
                        Epeck::FT(center_x)),
                    exact_rational_text(
                        Epeck::FT(center_y)),
                    exact_rational_text(guide_radius),
                },
                exact_rational_text(
                    Epeck::FT(tool_radius)),
                exact_rational_text(
                    Epeck::FT(cap_chord_ratio)),
                line_sources,
                circle_sources);
        const bool cap_exceeded =
            line_sources.size()
                == boundary_records.size()
            && has_material_rational_chart_witness(
                stock,
                Epeck::FT(center_x),
                Epeck::FT(center_y),
                guide_radius,
                Epeck::FT(tool_radius));
        EventTrace2 trace =
            build_event_trace(
                std::move(partition),
                "full-circle-four-chart-v1",
                motion_identity(
                    center_x,
                    center_y,
                    phase_dx,
                    phase_dy,
                    clockwise),
                encode_canonical_record(
                    "cap-chord-ratio-binary64-v1",
                    {
                        binary64_identity(
                            cap_chord_ratio),
                    }),
                cap_exceeded
                    ? ContinuousTeaVerdict::
                          CAP_EXCEEDED
                    : ContinuousTeaVerdict::
                          UNRESOLVED_DEGENERACY,
                cap_exceeded
                    ? "partial"
                    : "unresolved",
                cap_exceeded
                    ? "full-circle-line-vertex-witness-exact-v1"
                    : "full-circle-line-pullbacks-exact-v1",
                {});
        return {
            cap_exceeded
                ? "cap_exceeded"
                : "unresolved",
            std::move(trace),
        };
    }

    const std::vector<ParameterChart2> charts =
        parameter_charts();
    const ParameterChart2& anchor =
        center_chart(charts, 0);
    const PartitionEvent2 seam_event{
        "seam",
        anchor.start_seam_id,
        "full-circle-center-support-v1",
        "full-circle-center-domain-v1",
        anchor.start_seam_id,
        anchor.chart_id,
        "owned",
    };
    const EventPartitionCertificate2 scaffold =
        partition_cap_crossings(
            {"0", "1"},
            {"1"},
            {"0", "3"},
            {"1"},
            "64",
            "65",
            seam_event);
    const VerifiedEventPartition2 verified =
        verify_event_partition(scaffold);
    if (verified.verdict
        != ContinuousTeaVerdict::CERTIFIED) {
        throw IncompleteFullCircleOracleError(
            "full-circle scaffold failed exact reconstruction");
    }

    std::vector<EventTraceEvent2> events;
    events.reserve(
        verified.partition.roots.size());
    for (std::size_t index = 0;
         index < verified.partition.roots.size();
         ++index) {
        const AlgebraicRootRecord2& root =
            verified.partition.roots[index];
        const ParameterChart2& chart =
            center_chart(
                charts,
                index % CENTER_CHART_IDS.size());
        events.push_back(
            make_event_trace_event(
                root.root_id,
                encode_string_sequence(
                    {
                        "full-circle-scaffold-fibre-v1",
                        root.root_id,
                    }),
                "seam",
                {chart.start_seam_id},
                {chart.chart_id},
                root.multiplicity,
                "owned",
                index));
    }
    events = order_full_circle_events(
        verified,
        clockwise,
        std::move(events));

    const std::string identity =
        motion_identity(
            center_x,
            center_y,
            phase_dx,
            phase_dy,
            clockwise);
    EventTrace2 trace =
        build_event_trace(
            verified.partition,
            "full-circle-four-chart-v1",
            identity,
            encode_canonical_record(
                "cap-chord-ratio-binary64-v1",
                {
                    binary64_identity(
                        cap_chord_ratio),
                }),
            ContinuousTeaVerdict::
                UNRESOLVED_DEGENERACY,
            "unresolved",
            "full-circle-task5-blocked-v1",
            std::move(events));
    return {"unresolved", std::move(trace)};
}

bool full_circle_rational_probe_exceeds_cap_exact(
    const Stock2& stock,
    double center_x,
    double center_y,
    double phase_dx,
    double phase_dy,
    std::size_t chart,
    std::size_t numerator,
    std::size_t denominator,
    double tool_radius,
    double cap_chord_ratio)
{
    require_finite(center_x, "circle center x");
    require_finite(center_y, "circle center y");
    require_finite(phase_dx, "circle phase dx");
    require_finite(phase_dy, "circle phase dy");
    require_finite(tool_radius, "tool radius");
    require_finite(cap_chord_ratio, "cap chord ratio");
    if (chart >= CENTER_CHART_IDS.size()
        || denominator == 0
        || numerator >= denominator) {
        throw IncompleteFullCircleOracleError(
            "invalid full-circle rational probe");
    }
    if (phase_dx == 0.0 && phase_dy == 0.0) {
        throw IncompleteFullCircleOracleError(
            "full-circle phase vector must be nonzero");
    }
    if (!(tool_radius > 0.0)) {
        throw IncompleteFullCircleOracleError(
            "tool radius must be positive");
    }
    if (!(cap_chord_ratio > 0.0
          && cap_chord_ratio <= 4.0)) {
        throw IncompleteFullCircleOracleError(
            "cap chord ratio must lie in (0, 4]");
    }
    if (extract_boundary_records(stock).empty()) {
        return false;
    }
    throw IncompleteFullCircleOracleError(
        "exact rational probe requires the Task 5 full-circle pullback substrate");
}
