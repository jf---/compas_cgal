#pragma once

#include <array>
#include <stdexcept>
#include <string>

class ExactCircleChartAtlasError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

struct ExactCircleChartAtlasRecord2 {
    const char* identifier;
    std::array<int, 3> x_numerator;
    std::array<int, 3> y_numerator;
    std::array<int, 3> denominator;
    int domain_low;
    int domain_high;
    bool increasing;
    bool parameter_zero_owns_seam;
};

inline constexpr const char* EXACT_CIRCLE_CHART_ATLAS_VERSION =
    "exact-quarter-pythagorean-chart-v1";

inline constexpr std::array<ExactCircleChartAtlasRecord2, 4>
    EXACT_CIRCLE_CHART_ATLAS{{
        {"center-quarter-0-v1", {1, 0, -1}, {0, 2, 0}, {1, 0, 1}, 0, 1, true, true},
        {"center-quarter-1-v1", {0, -2, 0}, {1, 0, -1}, {1, 0, 1}, 0, 1, true, true},
        {"center-quarter-2-v1", {-1, 0, 1}, {0, -2, 0}, {1, 0, 1}, 0, 1, true, true},
        {"center-quarter-3-v1", {0, 2, 0}, {-1, 0, 1}, {1, 0, 1}, 0, 1, true, true},
    }};

inline const ExactCircleChartAtlasRecord2& exact_circle_chart_record(int chart)
{
    if (chart < 0 || chart >= static_cast<int>(EXACT_CIRCLE_CHART_ATLAS.size())) {
        throw ExactCircleChartAtlasError("exact circle chart index is outside [0, 3]");
    }
    return EXACT_CIRCLE_CHART_ATLAS[static_cast<std::size_t>(chart)];
}

inline std::string exact_circle_chart_id(int chart)
{
    return exact_circle_chart_record(chart).identifier;
}
