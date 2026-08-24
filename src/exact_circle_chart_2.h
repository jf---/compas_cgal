#pragma once

#include "exact_circle_chart_atlas_2.h"
#include "exact_motion_2.h"

#include <stdexcept>
#include <string>

class ExactCircleChartError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class ExactCircleChartParameter2 {
public:
    static ExactCircleChartParameter2 build(
        int chart,
        const Epeck::FT& parameter);

    int chart() const noexcept;
    const Epeck::FT& parameter() const noexcept;

private:
    ExactCircleChartParameter2(
        int chart,
        Epeck::FT parameter);

    int chart_;
    Epeck::FT parameter_;
};

EVector exact_circle_chart_vector(
    const EVector& phase,
    const ExactCircleChartParameter2& parameter);

EPoint exact_circle_chart_point(
    const EPoint& center,
    const EVector& phase,
    const ExactCircleChartParameter2& parameter);

const std::string& exact_circle_chart_strategy_version();
