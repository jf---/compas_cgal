#include "circle_source.h"

#include "../exact_stock_region_2.h"
#include "../exact_sweep_2.h"
#include "circle_projection.h"
#include "event_certificate.h"
#include "event_partition.h"
#include "parameter_charts.h"
#include "sha256.h"

#include <algorithm>
#include <array>
#include <string>
#include <vector>

constexpr std::array<const char*, 4>
    CENTER_CHART_IDS{
        "center-quarter-0-v1",
        "center-quarter-1-v1",
        "center-quarter-2-v1",
        "center-quarter-3-v1",
    };

const ParameterChart2& uniform_center_chart(
    const std::vector<ParameterChart2>& charts,
    std::size_t index)
{
    const auto found = std::find_if(
        charts.begin(),
        charts.end(),
        [index](const ParameterChart2& chart) {
            return chart.chart_id
                == CENTER_CHART_IDS.at(index);
        });
    if (found == charts.end()) {
        throw EventPartitionVerificationError(
            "full-circle uniform source lacks a center chart");
    }
    return *found;
}

std::string full_circle_stock_identity(
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

std::optional<std::string> full_circle_uniform_disposition(
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
    ReachSet overlap =
        reach_full_circle_sweep(
            ReachKernelPoint(
                ReachFT(center_x),
                ReachFT(center_y)),
            ReachKernelVector(
                ReachFT(phase_dx),
                ReachFT(phase_dy)),
            ReachFT(tool_radius));
    overlap.intersection(
        lift_exact_stock_region(stock));
    if (overlap.is_empty()) {
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
    const Epeck::FT& phase_x,
    const Epeck::FT& phase_y,
    const Epeck::FT& tool_radius)
{
    // Quarter-chart t=1/2 gives the exact unit direction (3/5, 4/5).
    const Epeck::FT chart_denominator(5);
    const Epeck::FT witness_x =
        center_x
        + (
              phase_x * Epeck::FT(3)
              - phase_y * Epeck::FT(4))
            / chart_denominator;
    const Epeck::FT witness_y =
        center_y
        + (
              phase_y * Epeck::FT(3)
              + phase_x * Epeck::FT(4))
            / chart_denominator;
    Gps cutter_disk;
    cutter_disk.insert(
        circle_disk(
            EPoint(witness_x, witness_y),
            tool_radius));
    cutter_disk.difference(stock.set());
    return cutter_disk.is_empty();
}


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
            uniform_center_chart(charts, index);
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
    if (stock_identity_value.empty()) {
        throw EventPartitionVerificationError(
            "full-circle boundary source lacks stock identity");
    }
    EventPartitionCertificate2 certificate =
        partition_full_circle_boundary_geometry(
            motion_data,
            cutter_radius,
            cap_chord_ratio,
            line_sources,
            circle_sources);
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
        "full-circle-boundary-pullbacks-v2";
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
