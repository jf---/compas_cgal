#pragma once

#include "exact_build_audit_2.h"
#include "exact_region_2.h"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

class InvalidCoverageGeometryError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class CoverageTransitionError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct CoverageSweepRecord2 {
    std::string strategy_version;
    std::string structural_record;
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
    const CoverageTransitionAudit2&
        last_transition_audit_for_native_gate() const;
    bool shares_pretransition_storage_with_for_native_gate(
        const Coverage2& other) const;

private:
    struct State {
        ExactRegion2 reachable_material;
        ExactRegion2 accumulated_sweeps;
        ExactRegion2 residual;
        std::vector<std::string> sweep_records;
        bool exact_residual_relation;
    };

    explicit Coverage2(State state);
    static State build_initial_state(
        const ExactRegion2& reachable_material,
        double precleared_x,
        double precleared_y,
        double precleared_radius);
    void apply_sweep(
        ReachSet sweep,
        std::string structural_record,
        CoverageTransitionAudit2 audit);

    State state_;
    CoverageTransitionAudit2 last_transition_audit_;
};
