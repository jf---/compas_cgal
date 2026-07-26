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
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using IntegerPolynomial = std::vector<Integer>;

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
    std::vector<std::string> result;
    result.reserve(discriminant.size());
    for (const Integer& coefficient :
         discriminant) {
        result.push_back(
            coefficient.convert_to<std::string>());
    }
    return result;
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
construct_full_circle_line_pullback_partition(
    const std::string& stock_identity_value,
    const std::vector<std::string>& motion_data,
    const std::string& cutter_radius,
    const std::vector<std::string>& line_sources)
{
    if (stock_identity_value.empty()
        || motion_data.size() != 3
        || cutter_radius.empty()
        || line_sources.empty()) {
        throw EventPartitionVerificationError(
            "full-circle line source is malformed");
    }
    std::vector<ProjectionRecord2> pullbacks;
    std::vector<ProjectionInput2> tangencies;
    for (const std::string& encoded_source :
         line_sources) {
        const std::vector<std::string> source =
            decode_string_sequence(encoded_source);
        if (source.size() != 4) {
            throw EventPartitionVerificationError(
                "full-circle line support source is malformed");
        }
        const std::vector<std::string> support =
            decode_string_sequence(source[3]);
        if (support.size() != 3) {
            throw EventPartitionVerificationError(
                "full-circle line support requires three coefficients");
        }
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
        "full-circle-line-pullbacks-v1";
    certificate.source_payload =
        encode_string_sequence(
            {
                stock_identity_value,
                encode_string_sequence(motion_data),
                cutter_radius,
                encode_string_sequence(line_sources),
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
    for (const BoundaryFeatureRecord2& record :
         boundary_records) {
        if (record.support_kind == "line") {
            line_sources.push_back(
                encode_string_sequence(
                    {
                        record.feature_id,
                        record.support_id,
                        record.trim_predicate,
                        encode_string_sequence(
                            record.primitive_coefficients),
                    }));
        }
    }
    const Epeck::FT exact_phase_x(phase_dx);
    const Epeck::FT exact_phase_y(phase_dy);
    if (!line_sources.empty()
        && (CGAL::is_zero(exact_phase_x)
            || CGAL::is_zero(exact_phase_y))) {
        const Epeck::FT guide_radius =
            CGAL::is_zero(exact_phase_y)
            ? CGAL::abs(exact_phase_x)
            : CGAL::abs(exact_phase_y);
        EventPartitionCertificate2 partition =
            construct_full_circle_line_pullback_partition(
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
                line_sources);
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
                ContinuousTeaVerdict::
                    UNRESOLVED_DEGENERACY,
                "unresolved",
                "full-circle-line-pullbacks-exact-v1",
                {});
        return {"unresolved", std::move(trace)};
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
