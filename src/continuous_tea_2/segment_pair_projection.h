#pragma once

#include "boundary_events.h"
#include "partition_certificate.h"
#include "segment_source.h"

#include <vector>

struct SegmentPairProjectionBundle2 {
    std::vector<ProjectionInput2> projections;
    std::vector<ProjectionRecord2> constants;
    std::vector<OverlapInterval2> overlaps;
};

SegmentPairProjectionBundle2
derive_segment_pair_projections(
    const std::vector<BoundaryFeatureRecord2>& records,
    const SegmentEventSource2& source,
    const std::vector<ProjectionRecord2>& pullbacks);
