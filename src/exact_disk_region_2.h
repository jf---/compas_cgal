#pragma once

#include "stock_2.h"

#include <stdexcept>
#include <vector>

class ExactDiskRegionEmptyCentersError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class ExactDiskRegionRadiusError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

Gps build_exact_disk_union_region_2(
    const std::vector<EPoint>& centers,
    const Epeck::FT& radius);
