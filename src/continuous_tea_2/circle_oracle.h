#pragma once

#include "event_trace.h"
#include "partition_certificate.h"

#include <string>
#include <utility>
#include <vector>

class Stock2;

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
