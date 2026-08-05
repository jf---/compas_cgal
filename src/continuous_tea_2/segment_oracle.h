#pragma once

#include "event_trace.h"

#include <cstddef>
#include <string>

class Stock2;
class SegmentEventSource2;

struct SegmentTeaAudit2 {
    ContinuousTeaVerdict verdict;
    EventTrace2 trace;
};

struct SweptPrefixSegmentTeaAudit2 {
    ContinuousTeaVerdict verdict;
    std::string canonical_bytes;
    std::string canonical_digest;
    std::string strategy_version;
    std::string theorem_version;
    std::string source_canonical_bytes;
    std::string stock_boundary_digest;
    std::size_t motion_stratum_count;
    bool start_disk_clear;
    bool exact_pi_cap;
};

class IncompleteSegmentOracleError
    : public EventSubstrateError {
public:
    using EventSubstrateError::EventSubstrateError;
};

// Task 6 consumes, but never reconstructs or approximates, Task 5's verified
// segment partition. The implementation must walk every exact open cell and
// event fibre, bind the resulting decision to the verified partition digest,
// and fail loudly when any stratum lacks a resolved run/cap disposition.
SegmentTeaAudit2 audit_segment_tea_event_exact(
    const Stock2& stock,
    const SegmentEventSource2& source);

// An advancing cutter consumes its own swept prefix. If the start disk has no
// positive-area stock overlap, every strictly trailing rim point at t > 0 was
// inside an earlier cutter disk. The possibly engaged rim is therefore a
// subset of the closed forward semicircle. This independent theorem certifies
// exactly a pi cap; it does not reinterpret an unresolved frozen-stock audit.
SweptPrefixSegmentTeaAudit2 audit_swept_prefix_segment_tea_exact(
    const Stock2& stock,
    const SegmentEventSource2& source);

const std::string& swept_prefix_strategy_version();
const std::string& swept_prefix_theorem_version();
bool swept_prefix_segment_tea_audit_is_self_consistent(
    const SweptPrefixSegmentTeaAudit2& audit);

bool segment_station_cap_exceeded_exact(
    const Stock2& stock,
    const SegmentEventSource2& source,
    const std::string& witness_numerator,
    const std::string& witness_denominator);
