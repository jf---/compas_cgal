#include "coverage_2.h"

#include "exact_sweep_2.h"

#include <cmath>
#include <string_view>
#include <utility>

namespace {

constexpr std::string_view COVERAGE_STRATEGY_VERSION =
    "exact-full-sweep-coverage-v1";

std::string tagged_record(
    std::string_view tag,
    const std::vector<std::string>& fields)
{
    std::string result(tag);
    result.push_back('\0');
    result += reach_u64_record(fields.size());
    for (const std::string& field : fields) {
        result += reach_length_prefixed(field);
    }
    return result;
}

ReachFT exact_positive(double value, const std::string& name)
{
    if (!std::isfinite(value)) {
        throw CoverageConstructionError(name + " must be finite.");
    }
    const ReachFT exact(value);
    if (CGAL::compare(exact, ReachFT(0)) != CGAL::LARGER) {
        throw CoverageConstructionError(name + " must be positive.");
    }
    return exact;
}

double finite_coordinate(double value, const std::string& name)
{
    if (!std::isfinite(value)) {
        throw CoverageConstructionError(name + " must be finite.");
    }
    return value == 0.0 ? 0.0 : value;
}

std::string point_record(double x, double y)
{
    return tagged_record(
        "coverage-point-v1",
        {reach_binary64_record(x), reach_binary64_record(y)});
}

std::string segment_record(
    double x0,
    double y0,
    double x1,
    double y1,
    double tool_radius)
{
    return tagged_record(
        "segment-sweep-v1",
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
    return tagged_record(
        "full-circle-sweep-v1",
        {
            point_record(cx, cy),
            point_record(phase_x, phase_y),
            reach_binary64_record(tool_radius),
        });
}

bool same_binary64(double left, double right)
{
    return reach_binary64_record(left) == reach_binary64_record(right);
}

} // namespace

bool CoverageSweepRecord2::matches_exact_segment(
    double x0,
    double y0,
    double x1,
    double y1,
    double expected_tool_radius) const
{
    return segment
        && same_binary64(center_x, x0)
        && same_binary64(center_y, y0)
        && same_binary64(first_x, x1)
        && same_binary64(first_y, y1)
        && same_binary64(tool_radius, expected_tool_radius)
        && structural_record
            == segment_record(x0, y0, x1, y1, expected_tool_radius);
}

bool CoverageSweepRecord2::matches_exact_full_circle(
    double cx,
    double cy,
    double phase_x,
    double phase_y,
    double expected_tool_radius) const
{
    return !segment
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
    : reachable_material_(reachable_material)
    , accumulated_sweeps_(
        ExactRegion2::build(
            ReachSet(),
            ExactRegionRole2::AccumulatedSweeps,
            ""))
    , residual_(
        ExactRegion2::build(
            ReachSet(),
            ExactRegionRole2::CoverageResidual,
            ""))
    , exact_residual_relation_(false)
{
    if (reachable_material.role() != ExactRegionRole2::ReachableMaterial) {
        throw CoverageConstructionError(
            "coverage requires an exact reachable-material region.");
    }
    const double cx =
        finite_coordinate(precleared_x, "precleared center x");
    const double cy =
        finite_coordinate(precleared_y, "precleared center y");
    const ReachFT radius =
        exact_positive(precleared_radius, "precleared radius");
    ReachSet precleared;
    precleared.insert(reach_disk_polygon(
        ReachKernelPoint(cx, cy),
        radius));
    accumulated_sweeps_ = ExactRegion2::build(
        std::move(precleared),
        ExactRegionRole2::AccumulatedSweeps,
        tagged_record(
            "coverage-precleared-v1",
            {
                reachable_material.recipe_record(),
                point_record(cx, cy),
                reach_binary64_record(precleared_radius),
            }));
    ReachSet residual(reachable_material_.set());
    residual.difference(accumulated_sweeps_.set());
    residual_ = ExactRegion2::build(
        std::move(residual),
        ExactRegionRole2::CoverageResidual,
        tagged_record(
            "coverage-residual-v1",
            {reachable_material.recipe_record()}));
    exact_residual_relation_ = exact_residual_relation();
    if (!exact_residual_relation_) {
        throw CoverageConstructionError(
            "precleared coverage residual failed exact reconstruction.");
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
    const ReachFT radius =
        exact_positive(tool_radius, "segment tool radius");
    const ReachKernelPoint start(x0, y0);
    const ReachKernelPoint end(x1, y1);
    if (start == end) {
        throw CoverageConstructionError(
            "exact segment sweep requires distinct endpoints.");
    }
    const ReachSet sweep = reach_join_parts(
        reach_capsule_parts(start, end, radius),
        {});
    const ReachSet replay = reach_join_parts(
        reach_capsule_parts(start, end, radius),
        {});
    const bool reconstructed = reach_exact_equal(sweep, replay);
    if (!reconstructed) {
        throw CoverageConstructionError(
            "exact segment sweep reconstruction changed.");
    }
    const std::string record =
        segment_record(x0, y0, x1, y1, tool_radius);
    apply_sweep(sweep, record);
    return {
        std::string(COVERAGE_STRATEGY_VERSION),
        record,
        reconstructed,
        true,
        x0,
        y0,
        x1,
        y1,
        tool_radius,
    };
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
    const ReachFT radius =
        exact_positive(tool_radius, "circle tool radius");
    const ReachKernelVector phase(phase_x, phase_y);
    if (phase == ReachKernelVector(CGAL::NULL_VECTOR)) {
        throw CoverageConstructionError(
            "exact full-circle sweep requires a nonzero phase vector.");
    }
    const ReachKernelPoint center(cx, cy);
    const ReachSet sweep =
        reach_full_circle_sweep(center, phase, radius);
    const ReachSet replay =
        reach_full_circle_sweep(center, phase, radius);
    const bool reconstructed = reach_exact_equal(sweep, replay);
    if (!reconstructed) {
        throw CoverageConstructionError(
            "exact full-circle sweep reconstruction changed.");
    }
    const std::string record =
        circle_record(cx, cy, phase_x, phase_y, tool_radius);
    apply_sweep(sweep, record);
    return {
        std::string(COVERAGE_STRATEGY_VERSION),
        record,
        reconstructed,
        false,
        cx,
        cy,
        phase_x,
        phase_y,
        tool_radius,
    };
}

void Coverage2::apply_sweep(
    const ReachSet& sweep,
    const std::string& structural_record)
{
    ReachSet accumulated(accumulated_sweeps_.set());
    accumulated.join(sweep);
    ReachSet residual(reachable_material_.set());
    residual.difference(accumulated);
    ReachSet replay(reachable_material_.set());
    replay.difference(accumulated);
    if (!reach_exact_equal(residual, replay)) {
        throw CoverageConstructionError(
            "coverage residual failed exact reconstruction.");
    }
    std::vector<std::string> next_records = sweep_records_;
    next_records.push_back(structural_record);
    accumulated_sweeps_ = ExactRegion2::build(
        std::move(accumulated),
        ExactRegionRole2::AccumulatedSweeps,
        tagged_record(
            "coverage-accumulated-sweeps-v1",
            next_records));
    residual_ = ExactRegion2::build(
        std::move(residual),
        ExactRegionRole2::CoverageResidual,
        tagged_record(
            "coverage-residual-v1",
            {
                reachable_material_.recipe_record(),
                tagged_record("coverage-sweep-lineage-v1", next_records),
            }));
    sweep_records_ = std::move(next_records);
    exact_residual_relation_ = exact_residual_relation();
    if (!exact_residual_relation_) {
        throw CoverageConstructionError(
            "coverage mutation broke the exact residual relation.");
    }
}

bool Coverage2::residual_is_empty() const
{
    return residual_.is_empty();
}

std::size_t Coverage2::residual_component_count() const
{
    return residual_.component_count();
}

ExactRegion2 Coverage2::residual() const
{
    return residual_.clone();
}

ExactRegion2 Coverage2::accumulated_sweeps() const
{
    return accumulated_sweeps_.clone();
}

std::vector<std::string> Coverage2::residual_component_records() const
{
    return reach_component_records(
        residual_.set(),
        "residual-component-v1");
}

std::vector<std::string> Coverage2::sweep_records() const
{
    return sweep_records_;
}

bool Coverage2::exact_residual_relation() const
{
    ReachSet expected(reachable_material_.set());
    expected.difference(accumulated_sweeps_.set());
    return reach_exact_equal(expected, residual_.set());
}

const std::string& Coverage2::strategy_version() const
{
    static const std::string version(COVERAGE_STRATEGY_VERSION);
    return version;
}
