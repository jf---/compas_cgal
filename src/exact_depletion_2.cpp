#include "exact_depletion_2.h"

#include <CGAL/Kernel/global_functions_2.h>

#include <algorithm>
#include <limits>
#include <utility>

namespace {

bool is_positive(const Epeck::FT& value)
{
    return CGAL::compare(value, Epeck::FT(0)) == CGAL::LARGER;
}

void validate_inputs(
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit)
{
    if (!is_positive(tool_radius)) {
        throw ExactDepletionConstructionError("exact depletion tool radius must be positive.");
    }
    if (!is_positive(max_chord)) {
        throw ExactDepletionConstructionError("exact depletion maximum chord must be positive.");
    }
    if (center_count_limit == 0) {
        throw ExactDepletionCenterLimitError("exact depletion center-count limit must be positive.");
    }
}

void require_segment_count(std::size_t intervals, std::size_t center_count_limit)
{
    if (intervals == std::numeric_limits<std::size_t>::max()) {
        throw ExactDepletionCenterLimitError("exact segment center count overflows size_t.");
    }
    const std::size_t center_count = intervals + 1;
    if (center_count > center_count_limit) {
        throw ExactDepletionCenterLimitError("exact segment center-count limit exceeded.");
    }
}

std::size_t circle_center_count(std::size_t intervals, std::size_t center_count_limit)
{
    if (intervals > std::numeric_limits<std::size_t>::max() / 4) {
        throw ExactDepletionCenterLimitError("exact circle center count overflows size_t.");
    }
    const std::size_t center_count = 4 * intervals;
    if (center_count > center_count_limit) {
        throw ExactDepletionCenterLimitError("exact circle center-count limit exceeded.");
    }
    return center_count;
}

EPoint segment_center(
    const ExactSegmentMotion2& motion,
    std::size_t numerator,
    std::size_t denominator)
{
    const Epeck::FT start_weight(denominator - numerator);
    const Epeck::FT end_weight(numerator);
    return CGAL::barycenter(motion.start, start_weight, motion.end, end_weight);
}

EVector quarter_chart_vector(
    const EVector& phase,
    int chart,
    std::size_t numerator,
    std::size_t denominator)
{
    const Epeck::FT t = Epeck::FT(numerator) / Epeck::FT(denominator);
    const Epeck::FT t_squared = t * t;
    const Epeck::FT chart_denominator = Epeck::FT(1) + t_squared;
    const Epeck::FT cosine = (Epeck::FT(1) - t_squared) / chart_denominator;
    const Epeck::FT sine = (Epeck::FT(2) * t) / chart_denominator;
    const Epeck::FT x = cosine * phase.x() - sine * phase.y();
    const Epeck::FT y = sine * phase.x() + cosine * phase.y();

    switch (chart) {
    case 0:
        return EVector(x, y);
    case 1:
        return EVector(-y, x);
    case 2:
        return EVector(-x, -y);
    case 3:
        return EVector(y, -x);
    default:
        throw ExactDepletionConstructionError("exact circle chart index is outside [0, 3].");
    }
}

EPoint circle_center(
    const ExactCircleMotion2& motion,
    int chart,
    std::size_t numerator,
    std::size_t denominator)
{
    return motion.center + quarter_chart_vector(
        motion.phase_vector,
        chart,
        numerator,
        denominator);
}

ExactCenterParameter2 ordered_circle_parameter(
    std::size_t ordered_index,
    std::size_t intervals,
    std::size_t center_count,
    bool clockwise)
{
    const std::size_t canonical_index =
        clockwise && ordered_index != 0
        ? center_count - ordered_index
        : ordered_index;
    return {
        static_cast<int>(canonical_index / intervals),
        canonical_index % intervals,
        intervals,
    };
}

bool segment_chords_hold(
    const ExactSegmentMotion2& motion,
    std::size_t intervals,
    const Epeck::FT& max_chord_squared)
{
    EPoint previous = segment_center(motion, 0, intervals);
    for (std::size_t index = 1; index <= intervals; ++index) {
        const EPoint current = segment_center(motion, index, intervals);
        if (CGAL::compare(CGAL::squared_distance(previous, current), max_chord_squared) == CGAL::LARGER) {
            return false;
        }
        previous = current;
    }
    return true;
}

bool circle_chords_hold(
    const ExactCircleMotion2& motion,
    std::size_t intervals,
    std::size_t center_count,
    const Epeck::FT& max_chord_squared,
    bool seam_only)
{
    const ExactCenterParameter2 first_parameter =
        ordered_circle_parameter(0, intervals, center_count, motion.clockwise);
    EPoint first = circle_center(
        motion,
        first_parameter.chart,
        first_parameter.numerator,
        first_parameter.denominator);
    EPoint previous = first;
    for (std::size_t index = 1; index < center_count; ++index) {
        const ExactCenterParameter2 parameter =
            ordered_circle_parameter(index, intervals, center_count, motion.clockwise);
        const EPoint current = circle_center(
            motion,
            parameter.chart,
            parameter.numerator,
            parameter.denominator);
        if (!seam_only && CGAL::compare(CGAL::squared_distance(previous, current), max_chord_squared) == CGAL::LARGER) {
            return false;
        }
        previous = current;
    }
    return CGAL::compare(CGAL::squared_distance(previous, first), max_chord_squared) != CGAL::LARGER;
}

bool segment_parameters_in_range(const std::vector<ExactCenterParameter2>& parameters)
{
    return std::all_of(
        parameters.begin(),
        parameters.end(),
        [](const ExactCenterParameter2& parameter) {
            return parameter.chart == -1
                && parameter.denominator != 0
                && parameter.numerator <= parameter.denominator;
        });
}

bool circle_parameters_in_range(const std::vector<ExactCenterParameter2>& parameters)
{
    return std::all_of(
        parameters.begin(),
        parameters.end(),
        [](const ExactCenterParameter2& parameter) {
            return parameter.chart >= 0
                && parameter.chart <= 3
                && parameter.denominator != 0
                && parameter.numerator < parameter.denominator;
        });
}

bool circle_anchors_present(const std::vector<ExactCenterParameter2>& parameters)
{
    for (int chart = 0; chart < 4; ++chart) {
        const bool present = std::any_of(
            parameters.begin(),
            parameters.end(),
            [chart](const ExactCenterParameter2& parameter) {
                return parameter.chart == chart && parameter.numerator == 0;
            });
        if (!present) {
            return false;
        }
    }
    return true;
}

} // namespace

const std::string& exact_depletion_strategy_version()
{
    static const std::string version = "exact-pythagorean-guide-v1";
    return version;
}

bool DepletionTrace::matches_exact_inputs(
    const Epeck::FT& expected_tool_radius,
    const Epeck::FT& expected_max_chord,
    std::size_t expected_center_count_limit) const
{
    return CGAL::compare(removal_radius, expected_tool_radius) == CGAL::EQUAL
        && CGAL::compare(max_chord, expected_max_chord) == CGAL::EQUAL
        && center_count <= expected_center_count_limit;
}

bool exact_segment_point_is_incident(
    const ExactSegmentMotion2& motion,
    const EPoint& point)
{
    return CGAL::orientation(motion.start, motion.end, point) == CGAL::COLLINEAR
        && CGAL::collinear_are_ordered_along_line(motion.start, point, motion.end);
}

bool exact_circle_point_is_incident(
    const ExactCircleMotion2& motion,
    const EPoint& point)
{
    return CGAL::compare(
               CGAL::squared_distance(motion.center, point),
               motion.phase_vector.squared_length())
        == CGAL::EQUAL;
}

ExactDepletionConstruction2 construct_exact_segment_depletion(
    const ExactSegmentMotion2& motion,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit)
{
    validate_inputs(tool_radius, max_chord, center_count_limit);
    if (motion.start == motion.end) {
        throw ExactDepletionConstructionError("exact segment motion requires distinct endpoints.");
    }

    const Epeck::FT max_chord_squared = max_chord * max_chord;
    std::size_t intervals = 1;
    while (true) {
        require_segment_count(intervals, center_count_limit);
        if (segment_chords_hold(motion, intervals, max_chord_squared)) {
            break;
        }
        if (intervals > std::numeric_limits<std::size_t>::max() / 2) {
            throw ExactDepletionCenterLimitError("exact segment refinement overflows size_t.");
        }
        intervals *= 2;
    }

    const std::size_t center_count = intervals + 1;
    std::vector<EPoint> centers;
    std::vector<ExactCenterParameter2> parameters;
    centers.reserve(center_count);
    parameters.reserve(center_count);
    for (std::size_t index = 0; index <= intervals; ++index) {
        centers.push_back(segment_center(motion, index, intervals));
        parameters.push_back({-1, index, intervals});
    }

    const bool exact_incidence = std::all_of(
        centers.begin(),
        centers.end(),
        [&motion](const EPoint& point) {
            return exact_segment_point_is_incident(motion, point);
        });
    DepletionTrace trace{
        std::move(parameters),
        center_count,
        max_chord,
        tool_radius,
        exact_depletion_strategy_version(),
        false,
        exact_incidence,
        false,
        false,
        is_positive(tool_radius),
        segment_chords_hold(motion, intervals, max_chord_squared),
        true,
    };
    trace.exact_parameters_in_range = segment_parameters_in_range(trace.center_parameters);
    trace.exact_anchors_present =
        trace.center_parameters.front().numerator == 0
        && trace.center_parameters.back().numerator == intervals;
    return {std::move(centers), std::move(trace)};
}

ExactDepletionConstruction2 construct_exact_full_circle_depletion(
    const ExactCircleMotion2& motion,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit)
{
    validate_inputs(tool_radius, max_chord, center_count_limit);
    if (motion.phase_vector == EVector(CGAL::NULL_VECTOR)) {
        throw ExactDepletionConstructionError("exact circle motion requires a nonzero phase vector.");
    }

    const Epeck::FT max_chord_squared = max_chord * max_chord;
    std::size_t intervals = 1;
    std::size_t center_count = 0;
    while (true) {
        center_count = circle_center_count(intervals, center_count_limit);
        if (circle_chords_hold(motion, intervals, center_count, max_chord_squared, false)) {
            break;
        }
        if (intervals > std::numeric_limits<std::size_t>::max() / 2) {
            throw ExactDepletionCenterLimitError("exact circle refinement overflows size_t.");
        }
        intervals *= 2;
    }

    std::vector<EPoint> centers;
    std::vector<ExactCenterParameter2> parameters;
    centers.reserve(center_count);
    parameters.reserve(center_count);
    for (std::size_t ordered_index = 0; ordered_index < center_count; ++ordered_index) {
        const ExactCenterParameter2 parameter =
            ordered_circle_parameter(ordered_index, intervals, center_count, motion.clockwise);
        centers.push_back(circle_center(
            motion,
            parameter.chart,
            parameter.numerator,
            parameter.denominator));
        parameters.push_back(parameter);
    }

    const bool exact_incidence = std::all_of(
        centers.begin(),
        centers.end(),
        [&motion](const EPoint& point) {
            return exact_circle_point_is_incident(motion, point);
        });
    DepletionTrace trace{
        std::move(parameters),
        center_count,
        max_chord,
        tool_radius,
        exact_depletion_strategy_version(),
        true,
        exact_incidence,
        false,
        false,
        is_positive(tool_radius),
        circle_chords_hold(motion, intervals, center_count, max_chord_squared, false),
        circle_chords_hold(motion, intervals, center_count, max_chord_squared, true),
    };
    trace.exact_parameters_in_range = circle_parameters_in_range(trace.center_parameters);
    trace.exact_anchors_present = circle_anchors_present(trace.center_parameters);
    return {std::move(centers), std::move(trace)};
}
