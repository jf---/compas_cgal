#pragma once

#include "event_trace.h"
#include "partition_certificate.h"
#include "circle_strata.h"

#include <string>
#include <optional>
#include <utility>
#include <vector>

class Stock2;
class FullCircleEventSource2;

struct FullCircleTeaAudit2 {
    std::string verdict;
    EventTrace2 trace;
    std::optional<FullCircleAuthorityParameter2> violating_parameter;
};

std::vector<EventTraceEvent2> order_full_circle_events(
    const VerifiedEventPartition2& verified_partition,
    bool clockwise,
    std::vector<EventTraceEvent2> events);

std::pair<std::string, EventTrace2>
audit_full_circle_tea_event_exact(
    const Stock2& stock,
    double center_x,
    double center_y,
    double phase_dx,
    double phase_dy,
    bool clockwise,
    double tool_radius,
    double cap_chord_ratio);

FullCircleTeaAudit2
audit_full_circle_tea_event_exact(
    const Stock2& stock,
    const FullCircleEventSource2& source);

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
    double cap_chord_ratio);

class IncompleteFullCircleOracleError
    : public EventSubstrateError {
public:
    using EventSubstrateError::EventSubstrateError;
};
