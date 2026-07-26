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

EventPartitionCertificate2
construct_full_circle_uniform_partition(
    const std::string& stock_identity,
    const std::string& disposition);

EventPartitionCertificate2
construct_full_circle_line_pullback_partition(
    const std::string& stock_identity,
    const std::vector<std::string>& motion_data,
    const std::string& cutter_radius,
    const std::vector<std::string>& line_sources);

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
