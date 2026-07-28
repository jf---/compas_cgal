#pragma once

#include "segment_partition.h"
#include "station_source.h"

#include <vector>

class Stock2;

enum class StationCellDecision {
    CLEAR,
    MATERIAL,
    SAFE_PARTIAL,
    CAP_EXCEEDED,
    UNRESOLVED,
};

struct StationCellClassification2 {
    StationCellDecision decision;
    bool reference_material;
    bool reference_resolved;
    std::vector<BranchPairDisposition2>
        material_runs;
};

StationCellClassification2 classify_station_cell(
    const Stock2& stock,
    const StationEventSource2& source,
    const SegmentEventStratum2& stratum);
