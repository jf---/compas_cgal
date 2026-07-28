#pragma once

#include "partition_certificate.h"

#include <cstddef>
#include <string>
#include <vector>

struct EventTraceEvent2 {
    std::string root_id;
    std::string global_fibre_id;
    std::string kind;
    std::vector<std::string> feature_ids;
    std::vector<std::string> branch_ids;
    unsigned int multiplicity;
    std::string disposition;
    std::size_t motion_order;
    std::string canonical_bytes;
    std::string canonical_id;
};

struct EventTrace2 {
    ContinuousTeaVerdict verdict;
    EventPartitionCertificate2 partition;
    std::string decision_authority_bytes;
    std::string decision_authority_digest;
    std::vector<EventTraceEvent2> events;
    std::string motion_chart_id;
    std::string motion_identity;
    std::string effective_cap_bytes;
    std::string whole_rim_disposition;
    std::string oracle_strategy_version;
    std::size_t event_cell_count;
    std::string canonical_bytes;
    std::string canonical_digest;
};

class EventTraceVerificationError
    : public EventSubstrateError {
public:
    using EventSubstrateError::EventSubstrateError;
};

const std::string& event_oracle_component_version();

EventTraceEvent2 make_event_trace_event(
    const std::string& root_id,
    const std::string& global_fibre_id,
    const std::string& kind,
    std::vector<std::string> feature_ids,
    std::vector<std::string> branch_ids,
    unsigned int multiplicity,
    const std::string& disposition,
    std::size_t motion_order);

EventTrace2 build_event_trace(
    const EventPartitionCertificate2& partition,
    const std::string& motion_chart_id,
    const std::string& motion_identity,
    const std::string& effective_cap_bytes,
    ContinuousTeaVerdict verdict,
    const std::string& whole_rim_disposition,
    const std::string& oracle_strategy_version,
    std::vector<EventTraceEvent2> events);

EventTrace2 build_authority_event_trace(
    VerifiedEventPartition2 verified_partition,
    const std::string& decision_authority_bytes,
    const std::string& decision_authority_digest,
    const std::string& motion_chart_id,
    const std::string& motion_identity,
    const std::string& effective_cap_bytes,
    ContinuousTeaVerdict verdict,
    const std::string& whole_rim_disposition,
    const std::string& oracle_strategy_version,
    std::vector<EventTraceEvent2> events);
