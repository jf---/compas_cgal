#pragma once

#include "exact_region_2.h"
#include "stock_2.h"

#include <stdexcept>

class ExactStockRegionLiftError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

ReachSet lift_exact_stock_region(const Stock2& stock);
