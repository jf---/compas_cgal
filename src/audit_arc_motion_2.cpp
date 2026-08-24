#include "audit_arc_motion_2.h"

#include "continuous_tea_2/sha256.h"

#include <CGAL/number_utils.h>
#include <CGAL/Fraction_traits.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <utility>

namespace {

constexpr double QUARTER_TURN = std::numbers::pi / 2.0;
constexpr double FULL_TURN = std::numbers::pi * 2.0;

struct NormalizedChartCoordinate2 {
    int chart;
    Epeck::FT parameter;
};

struct ExactTurnCount2 {
    Epeck::FT exact;
    std::int64_t binary64_exact;
};

ExactTurnCount2 exact_turn_count(const Epeck::FT& value)
{
    using Exact = Epeck::FT::Exact_type;
    using Traits = CGAL::Fraction_traits<Exact>;
    typename Traits::Numerator_type numerator;
    typename Traits::Denominator_type denominator;
    typename Traits::Decompose()(value.exact(), numerator, denominator);
    typename Traits::Numerator_type quotient = numerator / denominator;
    if (numerator < 0 && numerator % denominator != 0) {
        --quotient;
    }
    const Epeck::FT exact_quotient{Exact(quotient)};
    constexpr std::int64_t MAX_BINARY64_EXACT_INTEGER = INT64_C(1) << 53;
    if (CGAL::compare(
            exact_quotient,
            Epeck::FT(MAX_BINARY64_EXACT_INTEGER))
            == CGAL::LARGER
        || CGAL::compare(
               exact_quotient,
               Epeck::FT(-MAX_BINARY64_EXACT_INTEGER))
            == CGAL::SMALLER) {
        throw AuditArcAngleNormalizationError(
            "arc whole-turn index is not exactly representable at the binary64 seam");
    }
    return {
        exact_quotient,
        quotient.template convert_to<std::int64_t>(),
    };
}

NormalizedChartCoordinate2 normalized_chart_coordinate(double angle)
{
    const Epeck::FT exact_tau(FULL_TURN);
    const Epeck::FT exact_angle(angle);
    const ExactTurnCount2 whole_turns = exact_turn_count(exact_angle / exact_tau);
    const Epeck::FT normalized = exact_angle - whole_turns.exact * exact_tau;

    const Epeck::FT exact_quarter(QUARTER_TURN);
    int chart = 0;
    while (chart < 3
           && CGAL::compare(normalized, Epeck::FT(chart + 1) * exact_quarter)
               != CGAL::SMALLER) {
        ++chart;
    }
    const Epeck::FT local = normalized - Epeck::FT(chart) * exact_quarter;
    if (CGAL::sign(local) == CGAL::ZERO) {
        return {chart, Epeck::FT(0)};
    }
    const double local_angle =
        (angle - static_cast<double>(whole_turns.binary64_exact) * FULL_TURN)
        - static_cast<double>(chart) * QUARTER_TURN;
    const double parameter = std::tan(local_angle / 2.0);
    if (!std::isfinite(parameter)) {
        throw AuditArcAngleNormalizationError(
            "arc angle-to-chart seam produced a non-finite parameter");
    }
    const Epeck::FT exact_parameter(parameter);
    if (CGAL::compare(exact_parameter, Epeck::FT(0)) != CGAL::LARGER
        || CGAL::compare(exact_parameter, Epeck::FT(1)) != CGAL::SMALLER) {
        throw AuditArcAngleNormalizationError(
            "arc angle-to-chart seam produced an out-of-chart parameter");
    }
    return {chart, exact_parameter};
}

void append_interval_bytes(
    std::string& bytes,
    const ExactArcChartInterval2& interval)
{
    append_audit_bytes(bytes, std::string(1, static_cast<char>(interval.chart())));
    append_audit_bytes(bytes, canonical_audit_rational_bytes(interval.start_parameter()));
    append_audit_bytes(bytes, canonical_audit_rational_bytes(interval.end_parameter()));
    append_audit_bytes(bytes, std::string(1, interval.increasing() ? '\1' : '\0'));
    append_audit_bytes(bytes, std::string(1, interval.owns_start_seam() ? '\1' : '\0'));
    append_audit_bytes(bytes, std::string(1, interval.owns_end_seam() ? '\1' : '\0'));
}

std::vector<ExactArcChartInterval2> build_intervals(
    const NormalizedChartCoordinate2& start,
    const NormalizedChartCoordinate2& end,
    bool clockwise,
    bool full_turn)
{
    std::vector<ExactArcChartInterval2> intervals;
    intervals.reserve(5);
    int chart = start.chart;
    Epeck::FT parameter = start.parameter;
    bool crossed = false;

    for (std::size_t guard = 0; guard < 6; ++guard) {
        if (!full_turn && chart == end.chart) {
            const CGAL::Comparison_result order = CGAL::compare(parameter, end.parameter);
            if ((!clockwise && order == CGAL::SMALLER)
                || (clockwise && order == CGAL::LARGER)) {
                intervals.push_back(ExactArcChartInterval2::build(
                    chart,
                    parameter,
                    end.parameter,
                    !clockwise,
                    CGAL::sign(parameter) == CGAL::ZERO,
                    CGAL::sign(end.parameter) == CGAL::ZERO));
                return intervals;
            }
        }
        if (crossed && chart == end.chart) {
            if (CGAL::compare(parameter, end.parameter) == CGAL::EQUAL) {
                return intervals;
            }
            intervals.push_back(ExactArcChartInterval2::build(
                chart,
                parameter,
                end.parameter,
                !clockwise,
                CGAL::sign(parameter) == CGAL::ZERO,
                CGAL::sign(end.parameter) == CGAL::ZERO));
            return intervals;
        }

        const Epeck::FT boundary = clockwise ? Epeck::FT(0) : Epeck::FT(1);
        if (CGAL::compare(parameter, boundary) != CGAL::EQUAL) {
            intervals.push_back(ExactArcChartInterval2::build(
                chart,
                parameter,
                boundary,
                !clockwise,
                CGAL::sign(parameter) == CGAL::ZERO,
                CGAL::sign(boundary) == CGAL::ZERO));
        }
        chart = clockwise ? (chart + 3) % 4 : (chart + 1) % 4;
        parameter = clockwise ? Epeck::FT(1) : Epeck::FT(0);
        crossed = true;
    }
    throw AuditArcTraversalError("arc chart traversal did not terminate canonically");
}

} // namespace

ExactArcChartInterval2 ExactArcChartInterval2::build(
    int chart,
    const Epeck::FT& start_parameter,
    const Epeck::FT& end_parameter,
    bool increasing,
    bool owns_start_seam,
    bool owns_end_seam)
{
    static_cast<void>(ExactCircleChartParameter2::build(chart, start_parameter));
    static_cast<void>(ExactCircleChartParameter2::build(chart, end_parameter));
    const CGAL::Comparison_result order = CGAL::compare(start_parameter, end_parameter);
    if (order == CGAL::EQUAL
        || (increasing && order != CGAL::SMALLER)
        || (!increasing && order != CGAL::LARGER)) {
        throw AuditArcIntervalDirectionError(
            "arc chart interval direction is inconsistent");
    }
    if (owns_start_seam != (CGAL::sign(start_parameter) == CGAL::ZERO)
        || owns_end_seam != (CGAL::sign(end_parameter) == CGAL::ZERO)) {
        throw AuditArcSeamOwnershipError(
            "arc chart interval seam ownership is inconsistent");
    }
    return ExactArcChartInterval2(
        chart,
        start_parameter,
        end_parameter,
        increasing,
        owns_start_seam,
        owns_end_seam);
}

ExactArcChartInterval2::ExactArcChartInterval2(
    int chart,
    Epeck::FT start_parameter,
    Epeck::FT end_parameter,
    bool increasing,
    bool owns_start_seam,
    bool owns_end_seam)
    : chart_(chart),
      start_parameter_(std::move(start_parameter)),
      end_parameter_(std::move(end_parameter)),
      increasing_(increasing),
      owns_start_seam_(owns_start_seam),
      owns_end_seam_(owns_end_seam)
{
}

int ExactArcChartInterval2::chart() const noexcept { return chart_; }
const Epeck::FT& ExactArcChartInterval2::start_parameter() const noexcept { return start_parameter_; }
const Epeck::FT& ExactArcChartInterval2::end_parameter() const noexcept { return end_parameter_; }
bool ExactArcChartInterval2::increasing() const noexcept { return increasing_; }
bool ExactArcChartInterval2::owns_start_seam() const noexcept { return owns_start_seam_; }
bool ExactArcChartInterval2::owns_end_seam() const noexcept { return owns_end_seam_; }

AuditArcMotion2 AuditArcMotion2::build(
    const EPoint& center,
    const EVector& zero_phase,
    const Epeck::FT& guide_radius,
    double authored_start_angle,
    double authored_end_angle,
    bool clockwise,
    const Epeck::FT& cut_z)
{
    if (!std::isfinite(authored_start_angle) || !std::isfinite(authored_end_angle)) {
        throw AuditArcNonFiniteInputError(
            "arc authored angles must be finite binary64 values");
    }
    if (CGAL::sign(guide_radius) != CGAL::POSITIVE
        || CGAL::compare(zero_phase.squared_length(), guide_radius * guide_radius)
            != CGAL::EQUAL) {
        throw AuditArcRadiusIncidenceError(
            "arc zero phase must be exactly incident on its positive guide radius");
    }
    const double authored_signed_sweep =
        authored_end_angle - authored_start_angle;
    if (!std::isfinite(authored_signed_sweep)) {
        throw AuditArcNonFiniteInputError(
            "arc authored binary64 sweep observation must be finite");
    }
    const Epeck::FT signed_sweep(authored_signed_sweep);
    const Epeck::FT exact_tau(FULL_TURN);
    const CGAL::Sign sweep_sign = CGAL::sign(signed_sweep);
    if (sweep_sign == CGAL::ZERO
        || CGAL::compare(CGAL::abs(signed_sweep), exact_tau) == CGAL::LARGER) {
        throw AuditArcSweepExtentError(
            "arc sweep must lie in the injected interval [-2*pi, 2*pi] excluding zero");
    }
    if (clockwise != (sweep_sign == CGAL::NEGATIVE)) {
        throw AuditArcOrientationError("arc orientation contradicts exact signed sweep");
    }
    const bool full_turn = CGAL::compare(CGAL::abs(signed_sweep), exact_tau) == CGAL::EQUAL;
    const NormalizedChartCoordinate2 start = normalized_chart_coordinate(authored_start_angle);
    const NormalizedChartCoordinate2 end = full_turn
        ? start
        : normalized_chart_coordinate(authored_end_angle);
    const ExactCircleChartParameter2 canonical_start =
        ExactCircleChartParameter2::build(start.chart, start.parameter);
    const ExactCircleChartParameter2 canonical_end =
        ExactCircleChartParameter2::build(end.chart, end.parameter);
    std::vector<ExactArcChartInterval2> intervals = build_intervals(start, end, clockwise, full_turn);
    const EPoint start_point = exact_circle_chart_point(
        center,
        zero_phase,
        canonical_start);
    const EPoint end_point = full_turn
        ? start_point
        : exact_circle_chart_point(
              center,
              zero_phase,
              canonical_end);
    if (!full_turn && start_point == end_point) {
        throw AuditArcEndpointCollapseError("distinct authored arc endpoints collapse in the chart surrogate");
    }

    std::string canonical("audit-arc-motion-v1");
    append_audit_bytes(canonical, audit_arc_motion_strategy_version());
    append_audit_bytes(canonical, exact_circle_chart_strategy_version());
    append_audit_bytes(canonical, canonical_audit_rational_bytes(center.x()));
    append_audit_bytes(canonical, canonical_audit_rational_bytes(center.y()));
    append_audit_bytes(canonical, canonical_audit_rational_bytes(zero_phase.x()));
    append_audit_bytes(canonical, canonical_audit_rational_bytes(zero_phase.y()));
    append_audit_bytes(canonical, canonical_audit_rational_bytes(guide_radius));
    append_audit_bytes(canonical, canonical_audit_binary64_bytes(authored_start_angle));
    append_audit_bytes(canonical, canonical_audit_binary64_bytes(authored_end_angle));
    append_audit_bytes(canonical, canonical_audit_binary64_bytes(authored_signed_sweep));
    append_audit_bytes(canonical, std::string(1, clockwise ? '\1' : '\0'));
    append_audit_bytes(canonical, canonical_audit_rational_bytes(cut_z));
    for (const ExactArcChartInterval2& interval : intervals) {
        append_interval_bytes(canonical, interval);
    }
    return AuditArcMotion2(
        center,
        zero_phase,
        guide_radius,
        CGAL::abs(signed_sweep),
        clockwise,
        full_turn,
        cut_z,
        std::move(intervals),
        canonical_start,
        canonical_end,
        start_point,
        end_point,
        NativeMotionDigestAuthority2::hash_canonical(canonical));
}

AuditArcMotion2::AuditArcMotion2(
    EPoint center,
    EVector zero_phase,
    Epeck::FT guide_radius,
    Epeck::FT sweep,
    bool clockwise,
    bool full_turn,
    Epeck::FT cut_z,
    std::vector<ExactArcChartInterval2> intervals,
    ExactCircleChartParameter2 start_parameter,
    ExactCircleChartParameter2 end_parameter,
    EPoint start_point,
    EPoint end_point,
    NativeMotionDigest2 digest)
    : center_(std::move(center)),
      zero_phase_(std::move(zero_phase)),
      guide_radius_(std::move(guide_radius)),
      sweep_(std::move(sweep)),
      clockwise_(clockwise),
      full_turn_(full_turn),
      cut_z_(std::move(cut_z)),
      intervals_(std::move(intervals)),
      start_parameter_(std::move(start_parameter)),
      end_parameter_(std::move(end_parameter)),
      start_point_(std::move(start_point)),
      end_point_(std::move(end_point)),
      digest_(std::move(digest))
{
}

const EPoint& AuditArcMotion2::center() const noexcept { return center_; }
const EVector& AuditArcMotion2::zero_phase() const noexcept { return zero_phase_; }
const Epeck::FT& AuditArcMotion2::guide_radius() const noexcept { return guide_radius_; }
const Epeck::FT& AuditArcMotion2::sweep() const noexcept { return sweep_; }
bool AuditArcMotion2::clockwise() const noexcept { return clockwise_; }
bool AuditArcMotion2::full_turn() const noexcept { return full_turn_; }
const Epeck::FT& AuditArcMotion2::cut_z() const noexcept { return cut_z_; }
const std::vector<ExactArcChartInterval2>& AuditArcMotion2::intervals() const noexcept { return intervals_; }
const ExactCircleChartParameter2& AuditArcMotion2::start_parameter() const noexcept { return start_parameter_; }
const ExactCircleChartParameter2& AuditArcMotion2::end_parameter() const noexcept { return end_parameter_; }
const EPoint& AuditArcMotion2::start_point() const noexcept { return start_point_; }
const EPoint& AuditArcMotion2::end_point() const noexcept { return end_point_; }
const NativeMotionDigest2& AuditArcMotion2::digest() const noexcept { return digest_; }

bool exact_arc_point_is_incident(
    const AuditArcMotion2& motion,
    const EPoint& point)
{
    return CGAL::compare(
               CGAL::squared_distance(motion.center(), point),
               motion.guide_radius() * motion.guide_radius())
        == CGAL::EQUAL;
}

const std::string& audit_arc_motion_strategy_version()
{
    static const std::string version = "audit-arc-quarter-chart-binary64-v2";
    return version;
}
