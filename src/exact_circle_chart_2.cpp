#include "exact_circle_chart_2.h"

#include <CGAL/number_utils.h>

#include <utility>

ExactCircleChartParameter2 ExactCircleChartParameter2::build(
    int chart,
    const Epeck::FT& parameter)
{
    if (chart < 0 || chart > 3) {
        throw ExactCircleChartError("exact circle chart index is outside [0, 3]");
    }
    if (CGAL::compare(parameter, Epeck::FT(0)) == CGAL::SMALLER
        || CGAL::compare(parameter, Epeck::FT(1)) == CGAL::LARGER) {
        throw ExactCircleChartError("exact circle chart parameter is outside [0, 1]");
    }
    return ExactCircleChartParameter2(chart, parameter);
}

ExactCircleChartParameter2::ExactCircleChartParameter2(
    int chart,
    Epeck::FT parameter)
    : chart_(chart), parameter_(std::move(parameter))
{
}

int ExactCircleChartParameter2::chart() const noexcept
{
    return chart_;
}

const Epeck::FT& ExactCircleChartParameter2::parameter() const noexcept
{
    return parameter_;
}

EVector exact_circle_chart_vector(
    const EVector& phase,
    const ExactCircleChartParameter2& parameter)
{
    const Epeck::FT& t = parameter.parameter();
    const Epeck::FT t_squared = t * t;
    const ExactCircleChartAtlasRecord2& chart =
        exact_circle_chart_record(parameter.chart());
    const auto evaluate = [&t, &t_squared](const std::array<int, 3>& values) {
        return Epeck::FT(values[0])
            + Epeck::FT(values[1]) * t
            + Epeck::FT(values[2]) * t_squared;
    };
    const Epeck::FT unit_x = evaluate(chart.x_numerator);
    const Epeck::FT unit_y = evaluate(chart.y_numerator);
    const Epeck::FT denominator = evaluate(chart.denominator);
    return EVector(
        (unit_x * phase.x() - unit_y * phase.y()) / denominator,
        (unit_x * phase.y() + unit_y * phase.x()) / denominator);
}

EPoint exact_circle_chart_point(
    const EPoint& center,
    const EVector& phase,
    const ExactCircleChartParameter2& parameter)
{
    return center + exact_circle_chart_vector(phase, parameter);
}

const std::string& exact_circle_chart_strategy_version()
{
    static const std::string version = EXACT_CIRCLE_CHART_ATLAS_VERSION;
    return version;
}
