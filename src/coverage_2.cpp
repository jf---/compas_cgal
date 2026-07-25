#include "coverage_2.h"

#include "exact_sweep_2.h"

#include <cmath>
#include <exception>
#include <string_view>
#include <type_traits>
#include <utility>

namespace {

constexpr std::string_view COVERAGE_STRATEGY_VERSION =
    "incremental-exact-coverage-v2";

ReachFT exact_positive(double value, const std::string& name)
{
    if (!std::isfinite(value)) {
        throw InvalidCoverageGeometryError(
            name + " must be finite.");
    }
    const ReachFT exact(value);
    if (CGAL::compare(exact, ReachFT(0)) != CGAL::LARGER) {
        throw InvalidCoverageGeometryError(
            name + " must be positive.");
    }
    return exact;
}

double finite_coordinate(double value, const std::string& name)
{
    if (!std::isfinite(value)) {
        throw InvalidCoverageGeometryError(
            name + " must be finite.");
    }
    return value == 0.0 ? 0.0 : value;
}

bool finite_sweep_inputs(
    double first_x,
    double first_y,
    double second_x,
    double second_y,
    double tool_radius)
{
    return std::isfinite(first_x)
        && std::isfinite(first_y)
        && std::isfinite(second_x)
        && std::isfinite(second_y)
        && std::isfinite(tool_radius);
}

std::string point_record(double x, double y)
{
    return reach_tagged_record(
        "coverage-point-v2",
        {reach_binary64_record(x), reach_binary64_record(y)});
}

std::string segment_record(
    double x0,
    double y0,
    double x1,
    double y1,
    double tool_radius)
{
    return reach_tagged_record(
        "coverage-segment-sweep-v2",
        {
            point_record(x0, y0),
            point_record(x1, y1),
            reach_binary64_record(tool_radius),
        });
}

std::string circle_record(
    double cx,
    double cy,
    double phase_x,
    double phase_y,
    double tool_radius)
{
    return reach_tagged_record(
        "coverage-full-circle-sweep-v2",
        {
            point_record(cx, cy),
            point_record(phase_x, phase_y),
            reach_binary64_record(tool_radius),
        });
}

bool same_binary64(double left, double right)
{
    return std::isfinite(left)
        && std::isfinite(right)
        && reach_binary64_record(left) == reach_binary64_record(right);
}

[[noreturn]] void throw_transition_error(
    std::string_view operation,
    const std::exception& error)
{
    throw CoverageTransitionError(
        std::string(operation) + " failed: " + error.what());
}

} // namespace

bool CoverageSweepRecord2::matches_exact_segment(
    double x0,
    double y0,
    double x1,
    double y1,
    double expected_tool_radius) const
{
    if (!finite_sweep_inputs(
            x0,
            y0,
            x1,
            y1,
            expected_tool_radius)) {
        return false;
    }
    return strategy_version == COVERAGE_STRATEGY_VERSION
        && segment
        && same_binary64(center_x, x0)
        && same_binary64(center_y, y0)
        && same_binary64(first_x, x1)
        && same_binary64(first_y, y1)
        && same_binary64(tool_radius, expected_tool_radius)
        && structural_record
            == segment_record(
                x0,
                y0,
                x1,
                y1,
                expected_tool_radius);
}

bool CoverageSweepRecord2::matches_exact_full_circle(
    double cx,
    double cy,
    double phase_x,
    double phase_y,
    double expected_tool_radius) const
{
    if (!finite_sweep_inputs(
            cx,
            cy,
            phase_x,
            phase_y,
            expected_tool_radius)) {
        return false;
    }
    return strategy_version == COVERAGE_STRATEGY_VERSION
        && !segment
        && same_binary64(center_x, cx)
        && same_binary64(center_y, cy)
        && same_binary64(first_x, phase_x)
        && same_binary64(first_y, phase_y)
        && same_binary64(tool_radius, expected_tool_radius)
        && structural_record
            == circle_record(
                cx,
                cy,
                phase_x,
                phase_y,
                expected_tool_radius);
}

Coverage2::Coverage2(
    const ExactRegion2& reachable_material,
    double precleared_x,
    double precleared_y,
    double precleared_radius)
    : Coverage2(build_initial_state(
        reachable_material,
        precleared_x,
        precleared_y,
        precleared_radius))
{
}

Coverage2::Coverage2(State state)
    : state_(std::move(state))
{
}

Coverage2::State Coverage2::build_initial_state(
    const ExactRegion2& reachable_material,
    double precleared_x,
    double precleared_y,
    double precleared_radius)
{
    if (reachable_material.role()
        != ExactRegionRole2::ReachableMaterial) {
        throw CoverageTransitionError(
            "coverage requires an exact reachable-material region.");
    }
    const double center_x = finite_coordinate(
        precleared_x,
        "precleared center x");
    const double center_y = finite_coordinate(
        precleared_y,
        "precleared center y");
    const ReachFT exact_precleared_radius = exact_positive(
        precleared_radius,
        "precleared radius");

    try {
        const std::string accumulated_recipe = reach_tagged_record(
            "coverage-precleared-v2",
            {
                reachable_material.recipe_record(),
                reach_binary64_record(center_x),
                reach_binary64_record(center_y),
                reach_binary64_record(precleared_radius),
            });
        const std::string residual_recipe = reach_tagged_record(
            "coverage-residual-v2",
            {
                reachable_material.recipe_record(),
                accumulated_recipe,
            });
        ReachSet accumulated;
        accumulated.insert(reach_disk_polygon(
            ReachKernelPoint(center_x, center_y),
            exact_precleared_radius));
        ReachSet residual(reachable_material.set());
        residual.difference(accumulated);
        return State{
            reachable_material.clone(),
            ExactRegion2::build(
                std::move(accumulated),
                ExactRegionRole2::AccumulatedSweeps,
                accumulated_recipe),
            ExactRegion2::build(
                std::move(residual),
                ExactRegionRole2::CoverageResidual,
                residual_recipe),
            {},
            true,
        };
    }
    catch (const std::exception& error) {
        throw_transition_error(
            "initial exact coverage construction",
            error);
    }
}

Coverage2 Coverage2::clone() const
{
    return *this;
}

CoverageSweepRecord2 Coverage2::add_segment_sweep(
    double x0,
    double y0,
    double x1,
    double y1,
    double tool_radius)
{
    x0 = finite_coordinate(x0, "segment start x");
    y0 = finite_coordinate(y0, "segment start y");
    x1 = finite_coordinate(x1, "segment end x");
    y1 = finite_coordinate(y1, "segment end y");
    const ReachFT radius = exact_positive(
        tool_radius,
        "segment tool radius");
    const ReachKernelPoint start(x0, y0);
    const ReachKernelPoint end(x1, y1);
    if (start == end) {
        throw InvalidCoverageGeometryError(
            "exact segment sweep requires distinct endpoints.");
    }
    const std::string record = segment_record(
        x0,
        y0,
        x1,
        y1,
        tool_radius);
    CoverageSweepRecord2 result{
        std::string(COVERAGE_STRATEGY_VERSION),
        record,
        true,
        x0,
        y0,
        x1,
        y1,
        tool_radius,
    };
    CoverageTransitionAudit2 audit;
    try {
        ReachSet sweep = reach_join_parts(
            reach_capsule_parts(start, end, radius),
            {});
        ++audit.sweep_constructions;
        apply_sweep(
            std::move(sweep),
            record,
            audit);
    }
    catch (const CoverageTransitionError&) {
        throw;
    }
    catch (const std::exception& error) {
        throw_transition_error(
            "exact segment coverage transition",
            error);
    }
    static_assert(
        std::is_nothrow_move_constructible_v<
            CoverageSweepRecord2>);
    return result;
}

CoverageSweepRecord2 Coverage2::add_full_circle_sweep(
    double cx,
    double cy,
    double phase_x,
    double phase_y,
    double tool_radius)
{
    cx = finite_coordinate(cx, "circle center x");
    cy = finite_coordinate(cy, "circle center y");
    phase_x = finite_coordinate(phase_x, "circle phase x");
    phase_y = finite_coordinate(phase_y, "circle phase y");
    const ReachFT radius = exact_positive(
        tool_radius,
        "circle tool radius");
    const ReachKernelVector phase(phase_x, phase_y);
    if (phase == ReachKernelVector(CGAL::NULL_VECTOR)) {
        throw InvalidCoverageGeometryError(
            "exact full-circle sweep requires a nonzero phase vector.");
    }
    const std::string record = circle_record(
        cx,
        cy,
        phase_x,
        phase_y,
        tool_radius);
    CoverageSweepRecord2 result{
        std::string(COVERAGE_STRATEGY_VERSION),
        record,
        false,
        cx,
        cy,
        phase_x,
        phase_y,
        tool_radius,
    };
    CoverageTransitionAudit2 audit;
    try {
        ReachSet sweep = reach_full_circle_sweep(
            ReachKernelPoint(cx, cy),
            phase,
            radius);
        ++audit.sweep_constructions;
        apply_sweep(
            std::move(sweep),
            record,
            audit);
    }
    catch (const CoverageTransitionError&) {
        throw;
    }
    catch (const std::exception& error) {
        throw_transition_error(
            "exact full-circle coverage transition",
            error);
    }
    static_assert(
        std::is_nothrow_move_constructible_v<
            CoverageSweepRecord2>);
    return result;
}

void Coverage2::apply_sweep(
    ReachSet sweep,
    std::string structural_record,
    CoverageTransitionAudit2 audit)
{
    if (!state_.exact_residual_relation) {
        throw CoverageTransitionError(
            "coverage transition requires an exact residual induction state.");
    }

    ReachSet next_accumulated(state_.accumulated_sweeps.set());
    next_accumulated.join(sweep);
    ++audit.accumulated_unions;

    ReachSet next_residual(state_.residual.set());
    next_residual.difference(sweep);
    ++audit.residual_differences;

    std::vector<std::string> next_records =
        state_.sweep_records;
    next_records.push_back(structural_record);
    const std::string accumulated_recipe = reach_tagged_record(
        "coverage-accumulated-transition-v2",
        {
            state_.accumulated_sweeps.recipe_record(),
            structural_record,
        });
    const std::string residual_recipe = reach_tagged_record(
        "coverage-residual-transition-v2",
        {
            state_.residual.recipe_record(),
            structural_record,
        });
    State next_state{
        state_.reachable_material.clone(),
        ExactRegion2::build(
            std::move(next_accumulated),
            ExactRegionRole2::AccumulatedSweeps,
            accumulated_recipe),
        ExactRegion2::build(
            std::move(next_residual),
            ExactRegionRole2::CoverageResidual,
            residual_recipe),
        std::move(next_records),
        true,
    };
    static_assert(std::is_nothrow_move_assignable_v<State>);
    state_ = std::move(next_state);
    last_transition_audit_ = audit;
}

bool Coverage2::residual_is_empty() const
{
    return state_.residual.is_empty();
}

std::size_t Coverage2::residual_component_count() const
{
    return state_.residual.component_count();
}

ExactRegion2 Coverage2::residual() const
{
    return state_.residual.clone();
}

ExactRegion2 Coverage2::accumulated_sweeps() const
{
    return state_.accumulated_sweeps.clone();
}

std::vector<std::string> Coverage2::residual_component_records() const
{
    return reach_component_records(
        state_.residual.set(),
        "residual-component-v1");
}

std::vector<std::string> Coverage2::sweep_records() const
{
    return state_.sweep_records;
}

bool Coverage2::exact_residual_relation() const
{
    return state_.exact_residual_relation;
}

const std::string& Coverage2::strategy_version() const
{
    static const std::string version(COVERAGE_STRATEGY_VERSION);
    return version;
}

const CoverageTransitionAudit2&
Coverage2::last_transition_audit_for_native_gate() const
{
    return last_transition_audit_;
}

bool Coverage2::shares_pretransition_storage_with_for_native_gate(
    const Coverage2& other) const
{
    return state_.reachable_material.shares_storage_with_for_audit(
               other.state_.reachable_material)
        && state_.accumulated_sweeps.shares_storage_with_for_audit(
            other.state_.accumulated_sweeps)
        && state_.residual.shares_storage_with_for_audit(
            other.state_.residual);
}
