#include "circle_oracle.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>
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
