#pragma once

#include "event_trace.h"

class Stock2;
class SegmentEventSource2;

struct SegmentTeaAudit2 {
    ContinuousTeaVerdict verdict;
    EventTrace2 trace;
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

bool segment_station_cap_exceeded_exact(
    const Stock2& stock,
    const SegmentEventSource2& source,
    const std::string& witness_numerator,
    const std::string& witness_denominator);
