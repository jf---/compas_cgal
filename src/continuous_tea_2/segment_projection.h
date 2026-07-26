#pragma once

#include "boundary_events.h"
#include "partition_certificate.h"
#include "segment_source.h"

#include <vector>

struct SegmentProjectionBundle2 {
    std::vector<ProjectionRecord2> pullbacks;
    std::vector<ProjectionInput2> event_projections;
    std::vector<ProjectionRecord2>
        constant_event_projections;
    std::vector<OverlapInterval2> overlaps;
};

SegmentProjectionBundle2 derive_segment_projections(
    const std::vector<BoundaryFeatureRecord2>& records,
    const SegmentEventSource2& source);
