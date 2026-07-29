#pragma once

#include "partition_certificate.h"

#include <string>
#include <vector>

struct FullCirclePairProjectionBundle2 {
    std::vector<ProjectionInput2> projections;
    std::vector<ProjectionRecord2> constants;
};

FullCirclePairProjectionBundle2
derive_full_circle_pair_cap_projections(
    const std::vector<std::string>& line_sources,
    const std::vector<std::string>& circle_sources,
    const std::string& cap_chord_ratio,
    const std::vector<ProjectionRecord2>& pullbacks);
