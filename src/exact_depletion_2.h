#pragma once

#include "audit_arc_motion_2.h"
#include "exact_circle_chart_2.h"
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

class ExactArcDepletionPolicyError : public ExactDepletionConstructionError {
public:
    using ExactDepletionConstructionError::ExactDepletionConstructionError;
};

class ExactArcForgedTraceError : public ExactDepletionConstructionError {
public:
    using ExactDepletionConstructionError::ExactDepletionConstructionError;
};

class NonFiniteExactArcDepletionInputError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class ExactArcDepletionTrace2 {
public:
    static ExactArcDepletionTrace2 build(
        const AuditArcMotion2& motion,
        const Epeck::FT& tool_radius,
        const Epeck::FT& max_chord,
        std::size_t center_count_limit,
        std::vector<ExactCircleChartParameter2> parameters);

    const std::vector<ExactCircleChartParameter2>& parameters() const noexcept;
    const std::string& canonical_bytes() const noexcept;
    const ExactDepletionTraceDigest2& digest() const noexcept;
    const std::string& strategy_version() const noexcept;
    bool cyclic() const noexcept;
    bool matches_exact_inputs(
        const Epeck::FT& expected_tool_radius,
        const Epeck::FT& expected_max_chord,
        std::size_t expected_center_count_limit) const;
    bool matches_motion(const AuditArcMotion2& expected_motion) const noexcept;

private:
    ExactArcDepletionTrace2(
        NativeMotionDigest2 motion_digest,
        Epeck::FT tool_radius,
        Epeck::FT max_chord,
        std::size_t center_count_limit,
        std::vector<ExactCircleChartParameter2> parameters,
        std::string canonical_bytes,
        ExactDepletionTraceDigest2 digest,
        bool cyclic);

    NativeMotionDigest2 motion_digest_;
    Epeck::FT tool_radius_;
    Epeck::FT max_chord_;
    std::size_t center_count_limit_;
    std::vector<ExactCircleChartParameter2> parameters_;
    std::string canonical_bytes_;
    ExactDepletionTraceDigest2 digest_;
    bool cyclic_;
};

struct ExactArcDepletionConstruction2 {
    std::vector<EPoint> centers;
    ExactArcDepletionTrace2 trace;
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

ExactArcDepletionConstruction2 construct_exact_arc_depletion(
    const AuditArcMotion2& motion,
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

bool exact_arc_structural_density_holds(
    const AuditArcMotion2& motion,
    const Epeck::FT& max_chord,
    const std::vector<ExactCircleChartParameter2>& parameters);

const std::string& exact_depletion_strategy_version();
const std::string& exact_arc_depletion_strategy_version();
