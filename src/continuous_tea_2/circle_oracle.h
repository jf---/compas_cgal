#pragma once

#include "event_trace.h"
#include "partition_certificate.h"

#include <vector>

std::vector<EventTraceEvent2> order_full_circle_events(
    const VerifiedEventPartition2& verified_partition,
    bool clockwise,
    std::vector<EventTraceEvent2> events);
