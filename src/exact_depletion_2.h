#pragma once

#include "exact_motion_2.h"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

class ExactDepletionConstructionError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class ExactDepletionCenterLimitError : public ExactDepletionConstructionError {
public:
    using ExactDepletionConstructionError::ExactDepletionConstructionError;
};

struct ExactCenterParameter2 {
    int chart;
    std::size_t numerator;
    std::size_t denominator;
};

struct DepletionTrace {
    std::vector<ExactCenterParameter2> center_parameters;
    std::size_t center_count;
    Epeck::FT max_chord;
    Epeck::FT removal_radius;
    std::string strategy_version;
    bool cyclic;
    bool exact_incidence;
    bool exact_parameters_in_range;
    bool exact_anchors_present;
    bool exact_removal_radius_valid;
    bool exact_chord_bound_holds;
    bool exact_seam_chord_bound_holds;

    bool matches_exact_inputs(
        const Epeck::FT& expected_tool_radius,
        const Epeck::FT& expected_max_chord,
        std::size_t expected_center_count_limit) const;
};

struct ExactDepletionConstruction2 {
    std::vector<EPoint> centers;
    DepletionTrace trace;
};

ExactDepletionConstruction2 construct_exact_segment_depletion(
    const ExactSegmentMotion2& motion,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit);

ExactDepletionConstruction2 construct_exact_full_circle_depletion(
    const ExactCircleMotion2& motion,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit);

bool exact_segment_point_is_incident(
    const ExactSegmentMotion2& motion,
    const EPoint& point);

bool exact_circle_point_is_incident(
    const ExactCircleMotion2& motion,
    const EPoint& point);

bool exact_segment_structural_density_holds(
    const ExactSegmentMotion2& motion,
    const Epeck::FT& max_chord,
    const std::vector<ExactCenterParameter2>& parameters);

bool exact_full_circle_structural_density_holds(
    const ExactCircleMotion2& motion,
    const Epeck::FT& max_chord,
    const std::vector<ExactCenterParameter2>& parameters);

const std::string& exact_depletion_strategy_version();
