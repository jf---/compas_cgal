#include "circle_oracle.h"

#include "../exact_algebraic_1.h"
#include "../stock_2.h"
#include "boundary_events.h"
#include "cap_partition.h"
#include "circle_projection.h"
#include "circle_source.h"
#include "circle_strata.h"
#include "circle_vertex_source.h"
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
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

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

std::vector<EventTraceEvent2> boundary_trace_events(
    const EventPartitionCertificate2& partition,
    bool clockwise)
{
    std::vector<EventTraceEvent2> result;
    for (std::size_t fibre_index = 0;
         fibre_index < partition.fibres.size();
         ++fibre_index) {
        const EventFibre2& fibre =
            partition.fibres[fibre_index];
        const auto root = std::find_if(
            partition.roots.begin(),
            partition.roots.end(),
            [&fibre](
                const AlgebraicRootRecord2& candidate) {
                return candidate.root_id
                    == fibre.root_id;
            });
        if (root == partition.roots.end()) {
            throw IncompleteFullCircleOracleError(
                "full-circle fibre lacks its exact root");
        }
        std::vector<std::string> active_branch_ids;
        for (const ActiveBoundaryBranch2& branch :
             fibre.left_active_branches) {
            active_branch_ids.push_back(
                branch.branch_id);
        }
        for (const ActiveBoundaryBranch2& branch :
             fibre.right_active_branches) {
            active_branch_ids.push_back(
                branch.branch_id);
        }
        const std::size_t motion_order =
            clockwise && fibre_index != 0
            ? partition.fibres.size()
                - fibre_index
            : fibre_index;
        const auto same_physical_incidence =
            [](const PartitionEvent2& first,
               const PartitionEvent2& second) {
                return std::tie(
                           first.kind,
                           first.feature_id,
                           first.support_id,
                           first.trim_id,
                           first.vertex_id,
                           first.endpoint_role)
                    == std::tie(
                           second.kind,
                           second.feature_id,
                           second.support_id,
                           second.trim_id,
                           second.vertex_id,
                           second.endpoint_role);
            };
        for (std::size_t event_index = 0;
             event_index < fibre.events.size();
             ++event_index) {
            const PartitionEvent2& event =
                fibre.events[event_index];
            if (event_index != 0
                && !fibre.seam_id.empty()
                && same_physical_incidence(
                    fibre.events[event_index - 1],
                    event)) {
                continue;
            }
            std::vector<std::string> branch_ids =
                active_branch_ids;
            branch_ids.push_back(event.branch_id);
            unsigned int multiplicity =
                root->multiplicity;
            for (const LocalEventWitness2& witness :
                 fibre.local_event_witnesses) {
                if (witness.kind == event.kind
                    && witness.feature_id
                        == event.feature_id
                    && witness.support_id
                        == event.support_id
                    && witness.trim_id == event.trim_id
                    && witness.vertex_id
                        == event.vertex_id
                    && witness.endpoint_role
                        == event.endpoint_role
                    && (!fibre.seam_id.empty()
                        || witness.local_branch_id
                            == event.branch_id)) {
                    branch_ids.push_back(
                        witness.local_branch_id);
                    multiplicity =
                        witness.multiplicity;
                }
            }
            const std::string disposition =
                event.kind == "endpoint-order"
                    && !fibre.ccw_direction.empty()
                    && fibre.ccw_direction != "mixed"
                ? (
                      clockwise
                          ? fibre.cw_direction
                          : fibre.ccw_direction)
                : event.disposition;
            result.push_back(
                make_event_trace_event(
                    fibre.root_id,
                    encode_canonical_record(
                        fibre.seam_id.empty()
                            ? "full-circle-global-fibre-v2"
                            : "full-circle-phase-seam-fibre-v3",
                        fibre.seam_id.empty()
                            ? std::vector<std::string>{
                                  fibre.root_id,
                              }
                            : std::vector<std::string>{
                                  fibre.seam_id,
                                  encode_string_sequence(
                                      fibre.local_root_ids),
                              }),
                    event.kind,
                    {event.feature_id},
                    std::move(branch_ids),
                    multiplicity,
                    disposition,
                    motion_order));
        }
    }
    return result;
}

} // namespace

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
        full_circle_uniform_disposition(
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
                full_circle_stock_identity(boundary_records),
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
                "full-circle-uniform-event-exact-v2",
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
    std::vector<std::string> rational_line_sources;
    std::vector<std::string> rational_circle_sources;
    bool rational_boundary_vertices = true;
    for (const BoundaryFeatureRecord2& record :
         boundary_records) {
        const std::string source_x =
            FullCircleCoordinate2::from_exact(
                record.curve.source().x())
                .canonical_source();
        const std::string source_y =
            FullCircleCoordinate2::from_exact(
                record.curve.source().y())
                .canonical_source();
        const std::string target_x =
            FullCircleCoordinate2::from_exact(
                record.curve.target().x())
                .canonical_source();
        const std::string target_y =
            FullCircleCoordinate2::from_exact(
                record.curve.target().y())
                .canonical_source();
        const std::optional<std::string>
            rational_source_x =
                rational_coordinate_text(
                    record.curve.source().x());
        const std::optional<std::string>
            rational_source_y =
                rational_coordinate_text(
                    record.curve.source().y());
        const std::optional<std::string>
            rational_target_x =
                rational_coordinate_text(
                    record.curve.target().x());
        const std::optional<std::string>
            rational_target_y =
                rational_coordinate_text(
                    record.curve.target().y());
        const bool record_has_rational_vertices =
            rational_source_x.has_value()
            && rational_source_y.has_value()
            && rational_target_x.has_value()
            && rational_target_y.has_value();
        rational_boundary_vertices =
            rational_boundary_vertices
            && record_has_rational_vertices;
        const std::vector<std::string> source_fields{
            record.feature_id,
            record.support_id,
            record.trim_predicate,
            encode_string_sequence(
                record.primitive_coefficients),
            record.source_vertex_id,
            source_x,
            source_y,
            record.target_vertex_id,
            target_x,
            target_y,
        };
        std::vector<std::string>
            rational_source_fields;
        if (record_has_rational_vertices) {
            rational_source_fields = {
                record.feature_id,
                record.support_id,
                record.trim_predicate,
                encode_string_sequence(
                    record.primitive_coefficients),
                record.source_vertex_id,
                *rational_source_x,
                *rational_source_y,
                record.target_vertex_id,
                *rational_target_x,
                *rational_target_y,
            };
        }
        if (record.support_kind == "line") {
            line_sources.push_back(
                encode_string_sequence(
                    source_fields));
            if (record_has_rational_vertices) {
                rational_line_sources.push_back(
                    encode_string_sequence(
                        rational_source_fields));
            }
        } else if (record.support_kind == "circle") {
            circle_sources.push_back(
                encode_string_sequence(
                    source_fields));
            if (record_has_rational_vertices) {
                rational_circle_sources.push_back(
                    encode_string_sequence(
                        rational_source_fields));
            }
        }
    }
    const Epeck::FT exact_phase_x(phase_dx);
    const Epeck::FT exact_phase_y(phase_dy);
    if (!line_sources.empty()
        || !circle_sources.empty()) {
        const std::vector<std::string> motion_data{
            exact_rational_text(
                Epeck::FT(center_x)),
            exact_rational_text(
                Epeck::FT(center_y)),
            exact_rational_text(
                exact_phase_x),
            exact_rational_text(
                exact_phase_y),
        };
        const std::string cutter_radius =
            exact_rational_text(
                Epeck::FT(tool_radius));
        const std::string cap_ratio =
            exact_rational_text(
                Epeck::FT(cap_chord_ratio));
        EventPartitionCertificate2 partition =
            rational_boundary_vertices
            ? construct_full_circle_boundary_pullback_partition(
                  full_circle_stock_identity(
                      boundary_records),
                  motion_data,
                  cutter_radius,
                  cap_ratio,
                  rational_line_sources,
                  rational_circle_sources)
            : construct_full_circle_pair_closed_partition(
                  stock,
                  full_circle_stock_identity(
                      boundary_records),
                  motion_data,
                  cutter_radius,
                  cap_ratio,
                  line_sources,
                  circle_sources,
                  center_x,
                  center_y,
                  phase_dx,
                  phase_dy,
                  tool_radius,
                  cap_chord_ratio);
        const bool cap_exceeded =
            line_sources.size()
                == boundary_records.size()
            && has_material_rational_chart_witness(
                stock,
                Epeck::FT(center_x),
                Epeck::FT(center_y),
                exact_phase_x,
                exact_phase_y,
                Epeck::FT(tool_radius));
        std::vector<EventTraceEvent2> events =
            boundary_trace_events(
                partition,
                clockwise);
        const std::string exact_motion_identity =
            motion_identity(
                center_x,
                center_y,
                phase_dx,
                phase_dy,
                clockwise);
        const std::string exact_cap_identity =
            encode_canonical_record(
                "cap-chord-ratio-binary64-v1",
                {
                    binary64_identity(
                        cap_chord_ratio),
                });
        if (cap_exceeded) {
            EventTrace2 trace =
                build_event_trace(
                    std::move(partition),
                    "full-circle-four-chart-v1",
                    exact_motion_identity,
                    exact_cap_identity,
                    ContinuousTeaVerdict::
                        CAP_EXCEEDED,
                    "partial",
                    "full-circle-line-vertex-witness-exact-v1",
                    std::move(events));
            return {
                "cap_exceeded",
                std::move(trace),
            };
        }

        VerifiedEventPartition2 verified_partition =
            verify_event_partition(partition);
        if (verified_partition.verdict
            != ContinuousTeaVerdict::CERTIFIED) {
            throw IncompleteFullCircleOracleError(
                "full-circle boundary partition failed exact reconstruction");
        }
        const FullCircleCellAuthority2 authority =
            construct_full_circle_cell_authority(
                stock,
                verified_partition,
                center_x,
                center_y,
                phase_dx,
                phase_dy,
                tool_radius,
                cap_chord_ratio);
        EventTrace2 trace =
            build_authority_event_trace(
                std::move(verified_partition),
                authority.canonical_bytes,
                authority.canonical_digest,
                "full-circle-four-chart-v1",
                exact_motion_identity,
                exact_cap_identity,
                authority.verdict,
                authority.whole_rim_disposition,
                "full-circle-cell-strata-exact-v3",
                std::move(events));
        const std::string verdict =
            authority.verdict
                    == ContinuousTeaVerdict::CERTIFIED
                ? "certified"
                : (
                      authority.verdict
                              == ContinuousTeaVerdict::
                                  CAP_EXCEEDED
                          ? "cap_exceeded"
                          : "unresolved");
        return {
            verdict,
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
