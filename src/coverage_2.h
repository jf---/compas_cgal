#pragma once

#include "exact_motion_2.h"
#include "reachable_domain_2.h"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

class CoverageConstructionError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct CoverageSweepRecord2 {
    std::string strategy_version;
    std::string structural_record;
    bool exact_reconstruction;
    bool segment;
    double center_x;
    double center_y;
    double first_x;
    double first_y;
    double tool_radius;

    bool matches_exact_segment(
        double x0,
        double y0,
        double x1,
        double y1,
        double expected_tool_radius) const;

    bool matches_exact_full_circle(
        double cx,
        double cy,
        double phase_x,
        double phase_y,
        double expected_tool_radius) const;
};

class Coverage2 {
public:
    Coverage2(
        const ExactRegion2& reachable_material,
        double precleared_x,
        double precleared_y,
        double precleared_radius);

    Coverage2 clone() const;
    CoverageSweepRecord2 add_segment_sweep(
        double x0,
        double y0,
        double x1,
        double y1,
        double tool_radius);
    CoverageSweepRecord2 add_full_circle_sweep(
        double cx,
        double cy,
        double phase_x,
        double phase_y,
        double tool_radius);
    bool residual_is_empty() const;
    std::size_t residual_component_count() const;
    ExactRegion2 residual() const;
    ExactRegion2 accumulated_sweeps() const;
    std::vector<std::string> residual_component_records() const;
    std::vector<std::string> sweep_records() const;
    bool exact_residual_relation() const;
    const std::string& strategy_version() const;

private:
    void apply_sweep(
        const ReachSet& sweep,
        const std::string& structural_record);

    ExactRegion2 reachable_material_;
    ExactRegion2 accumulated_sweeps_;
    ExactRegion2 residual_;
    std::vector<std::string> sweep_records_;
    bool exact_residual_relation_;
};
