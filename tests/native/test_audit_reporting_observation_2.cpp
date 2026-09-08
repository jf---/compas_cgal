#include "audit_replay_fixtures_2.h"

#include "audit_reporting_observation_2.h"
#include "audit_reporting_schedule_2.h"
#include "audit_stock_state_identity_2.h"
#include "canonical_encoding.h"
#include "engagement_2.h"

#include <CGAL/number_utils.h>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>
#include <type_traits>
#include <vector>

using namespace audit_replay_fixtures;

namespace
{

template <class Authority>
concept PublicDigestRetagger = requires(const std::string& bytes)
{
    Authority::from_external_bytes(bytes);
};

static_assert(!std::is_default_constructible_v<AuditTeaReportingObservation2>);
static_assert(!std::is_default_constructible_v<AuditReportingRadian2>);
static_assert(!std::is_constructible_v<AuditTeaReportingObservation2,
                                       AuditReportingRadian2,
                                       std::size_t>);
static_assert(!PublicDigestRetagger<TeaReportingObservationDigestAuthority2>);

void require_schedule_points(const AuditReportingSchedule2& schedule,
                             const std::vector<EPoint>& expected,
                             const char* message)
{
    require(schedule.stations().size() == expected.size(), message);
    for (std::size_t index = 0; index < expected.size(); ++index)
    {
        require(schedule.stations()[index].point() == expected[index], message);
    }
    for (std::size_t left = 0; left < schedule.stations().size(); ++left)
    {
        for (std::size_t right = left + 1; right < schedule.stations().size(); ++right)
        {
            require(schedule.stations()[left].point() != schedule.stations()[right].point(),
                    "reporting schedule retains a duplicate exact point");
        }
    }
}

bool schedule_contains(
    const AuditReportingSchedule2& schedule,
    const EPoint& point)
{
    return std::any_of(
        schedule.stations().begin(),
        schedule.stations().end(),
        [&point](const AuditReportingStationWorldXY2& station)
        {
            return station.point() == point;
        });
}

void deterministic_exact_schedule_gate()
{
    // Production mutations caught: reporting schedules use binary64 geometry,
    // lose traversal order, omit owned seams, or retain duplicate endpoints.
    const AuditSegmentMotion2 segment_motion = segment();
    const AuditReportingSchedule2 first_segment =
        build_audit_reporting_schedule(segment_motion);
    const AuditReportingSchedule2 second_segment =
        build_audit_reporting_schedule(segment_motion);
    require_schedule_points(first_segment,
                            {EPoint(-1, 0), EPoint(0, 0), EPoint(1, 0)},
                            "segment reporting schedule is not start/midpoint/end");
    require(first_segment.canonical_bytes() == second_segment.canonical_bytes() &&
                first_segment.digest().bytes() == second_segment.digest().bytes(),
            "segment reporting schedule is not deterministic");

    const AuditCircleMotion2 circle_motion = circle(false);
    const AuditReportingSchedule2 circle_schedule =
        build_audit_reporting_schedule(circle_motion);
    require(circle_schedule.stations().size() == 8,
            "full-circle reporting schedule does not own four anchors and four midpoints");
    require(circle_schedule.motion_digest().bytes() == circle_motion.digest().bytes(),
            "circle reporting schedule does not bind its motion");

    const AuditArcMotion2 quarter_motion = arc();
    const AuditReportingSchedule2 quarter_schedule =
        build_audit_reporting_schedule(quarter_motion);
    require(quarter_schedule.stations().size() == 3,
            "quarter-arc reporting schedule is not start/midpoint/terminal");

    const AuditArcMotion2 full_turn = std::get<AuditArcMotion2>(classify_audit_arc(
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 1.0, 0.0,
        2.0 * std::numbers::pi, false, 0.0, 5.0, "cut"));
    const AuditReportingSchedule2 full_turn_schedule =
        build_audit_reporting_schedule(full_turn);
    require(full_turn_schedule.stations().size() == 8,
            "full-turn arc did not deduplicate its terminal/start anchor");
    require(full_turn_schedule.strategy_version() ==
                "audit-motion-reporting-schedule-v1",
            "reporting schedule lost its strategy identity");

    const AuditArcMotion2 seam_crossing = std::get<AuditArcMotion2>(
        classify_audit_arc(
            {0.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            1.0,
            3.0 * std::numbers::pi / 4.0,
            7.0 * std::numbers::pi / 4.0,
            false,
            0.0,
            5.0,
            "cut"));
    const AuditReportingSchedule2 seam_schedule =
        build_audit_reporting_schedule(seam_crossing);
    require(
        schedule_contains(seam_schedule, seam_crossing.start_point()) &&
            schedule_contains(seam_schedule, seam_crossing.end_point()),
        "partial-arc schedule omits an authored terminal");
    for (const ExactArcChartInterval2& interval : seam_crossing.intervals())
    {
        const Epeck::FT midpoint =
            (interval.start_parameter() + interval.end_parameter()) /
            Epeck::FT(2);
        require(
            schedule_contains(
                seam_schedule,
                exact_circle_chart_point(
                    seam_crossing.center(),
                    seam_crossing.zero_phase(),
                    ExactCircleChartParameter2::build(
                        interval.chart(), midpoint))),
            "partial-arc schedule omits a clipped-chart midpoint");
        if (interval.owns_start_seam())
        {
            require(
                schedule_contains(
                    seam_schedule,
                    exact_circle_chart_point(
                        seam_crossing.center(),
                        seam_crossing.zero_phase(),
                        ExactCircleChartParameter2::build(
                            interval.chart(), interval.start_parameter()))),
                "partial-arc schedule omits an owned start seam");
        }
        if (interval.owns_end_seam())
        {
            require(
                schedule_contains(
                    seam_schedule,
                    exact_circle_chart_point(
                        seam_crossing.center(),
                        seam_crossing.zero_phase(),
                        ExactCircleChartParameter2::build(
                            interval.chart(), interval.end_parameter()))),
                "partial-arc schedule omits an owned end seam");
        }
    }
}

void result_owned_observation_and_parity_gate()
{
    // Production mutations caught: a caller splices reporting onto a result,
    // reporting uses post-depletion stock, or the new span projection diverges
    // from the retained reporter on the same binary64 stations.
    const compas::RowMatrixXd boundary = rectangle(-5.0, -5.0, 5.0, -0.4);
    const AuditPolicy2 audit_policy = policy();
    const AuditDecisionLimits2 limits = decision_limits();
    const AuditSegmentMotion2 motion = segment();
    const AuthenticatedOperationDigest2 operation =
        operation_digest("reporting-owned");
    AuditReplay2 replay = one_motion_replay(
        boundary, audit_policy, limits, motion, operation);
    const AuditLateralResult2 result =
        audit_deplete_segment(replay, motion, operation);
    const AuditTeaReportingObservation2& observation =
        result.reporting_observation();

    require(std::isfinite(observation.max_tea().value()) &&
                observation.max_tea().value() >= 0.0 &&
                observation.max_tea().value() <= 2.0 * std::numbers::pi &&
                observation.station_count() == 3,
            "owned reporting observation violates finite range or schedule count");
    require(observation.request_digest().bytes() == result.request_digest().bytes() &&
                observation.cursor() == result.cursor() &&
                observation.authenticated_operation_digest().bytes() ==
                    result.authenticated_operation_digest().bytes() &&
                observation.motion_digest().bytes() == result.motion_digest().bytes() &&
                observation.policy_digest().bytes() == audit_policy.digest().bytes() &&
                observation.pre_stock_state_digest().bytes() ==
                    AuditStockStateIdentity2::build(Stock2(boundary, {})).digest().bytes() &&
                observation.pre_lineage().bytes() == result.pre_lineage().bytes() &&
                observation.native_decision_digest().bytes() ==
                    result.decision_digest().bytes() &&
                observation.native_result_digest().bytes() == result.digest().bytes(),
            "owned reporting observation lost native transaction context");

    Stock2 stock(boundary, {});
    double expected_max = 0.0;
    for (const AuditReportingStationWorldXY2& station : observation.schedule().stations())
    {
        const EPoint& point = station.point();
        expected_max = std::max(
            expected_max,
            engagement_at(stock,
                          CGAL::to_double(point.x()),
                          CGAL::to_double(point.y()),
                          0.5,
                          audit_cap_chord_ratio(std::numbers::pi / 2.0))
                .max_run_tea);
    }
    require(observation.max_tea().value() == expected_max,
            "reporting projection lost parity on shared binary64 stations");

    AuditReplay2 repeated = one_motion_replay(
        boundary, audit_policy, limits, motion, operation);
    const AuditLateralResult2 repeated_result =
        audit_deplete_segment(repeated, motion, operation);
    require(repeated_result.reporting_observation().canonical_bytes() ==
                observation.canonical_bytes() &&
                repeated_result.reporting_observation().digest().bytes() ==
                    observation.digest().bytes(),
            "reporting observation is not repeatable from identical native context");
}

void observation_component_mutation_gate()
{
    // Production mutation caught: observation authentication omits one native
    // context, schedule, reporting value, count, maximum, or strategy field.
    const compas::RowMatrixXd boundary = rectangle(-5.0, -5.0, 5.0, -0.4);
    const AuditPolicy2 audit_policy = policy();
    const AuditPolicy2 changed_policy = policy(2048);
    const AuditDecisionLimits2 limits = decision_limits();
    const AuditSegmentMotion2 motion = segment();
    const AuditSegmentMotion2 changed_motion = segment(1.0);
    const AuthenticatedOperationDigest2 operation = operation_digest("mutation");
    const AuditNativeRequestIdentity2 request =
        one_motion_request(boundary, audit_policy, limits, motion);
    const AuditNativeRequestIdentity2 changed_request =
        one_motion_request(boundary, audit_policy, limits, changed_motion);
    AuditReplay2 replay = one_motion_replay(
        boundary, audit_policy, limits, motion, operation);
    const AuditLateralResult2 result =
        audit_deplete_segment(replay, motion, operation);
    const AuditTeaReportingObservation2& observation =
        result.reporting_observation();

    const AuthenticatedOperationDigest2 changed_context_operation =
        operation_digest("changed-reporting-context");
    AuditReplay2 changed_context_replay = one_motion_replay(
        boundary,
        audit_policy,
        limits,
        changed_motion,
        changed_context_operation);
    const AuditLateralResult2 changed_context_result = audit_deplete_segment(
        changed_context_replay,
        changed_motion,
        changed_context_operation);

    Stock2 stock(boundary, {});
    Stock2 changed_stock(rectangle(-6.0, -5.0, 5.0, -0.4), {});
    const AuditStockStateDigest2& stock_digest =
        AuditStockStateIdentity2::build(stock).digest();
    const AuditStockStateDigest2& changed_stock_digest =
        AuditStockStateIdentity2::build(changed_stock).digest();
    const AuditLineage2 foreign_lineage = AuditReplayTestAuthority2::seed_lineage(
        input_digest("foreign-reporting-lineage"), request);

    std::vector<AuditReportingRadian2> changed_station_values =
        observation.station_maxima();
    changed_station_values[0] = AuditReportingRadian2::build(
        changed_station_values[0].value() + 0.125);
    const TeaReportingObservationDigest2 baseline =
        AuditReplayTestAuthority2::reporting_observation_digest(
            request.digest(), 0, operation, motion.digest(), audit_policy.digest(),
            stock_digest, result.pre_lineage(), result.decision_digest(),
            result.digest(), observation.schedule().digest(),
            observation.station_maxima(), observation.max_tea(),
            observation.station_count(), observation.strategy_version());
    const std::vector<TeaReportingObservationDigest2> mutations{
        AuditReplayTestAuthority2::reporting_observation_digest(
            changed_request.digest(), 0, operation, motion.digest(),
            audit_policy.digest(), stock_digest, result.pre_lineage(),
            result.decision_digest(), result.digest(), observation.schedule().digest(),
            observation.station_maxima(), observation.max_tea(),
            observation.station_count(), observation.strategy_version()),
        AuditReplayTestAuthority2::reporting_observation_digest(
            request.digest(), 1, operation, motion.digest(), audit_policy.digest(),
            stock_digest, result.pre_lineage(), result.decision_digest(),
            result.digest(), observation.schedule().digest(),
            observation.station_maxima(), observation.max_tea(),
            observation.station_count(), observation.strategy_version()),
        AuditReplayTestAuthority2::reporting_observation_digest(
            request.digest(), 0, operation_digest("changed-reporting-operation"),
            motion.digest(), audit_policy.digest(), stock_digest,
            result.pre_lineage(), result.decision_digest(), result.digest(),
            observation.schedule().digest(), observation.station_maxima(),
            observation.max_tea(), observation.station_count(),
            observation.strategy_version()),
        AuditReplayTestAuthority2::reporting_observation_digest(
            request.digest(), 0, operation, changed_motion.digest(),
            audit_policy.digest(), stock_digest, result.pre_lineage(),
            result.decision_digest(), result.digest(), observation.schedule().digest(),
            observation.station_maxima(), observation.max_tea(),
            observation.station_count(), observation.strategy_version()),
        AuditReplayTestAuthority2::reporting_observation_digest(
            request.digest(), 0, operation, motion.digest(), changed_policy.digest(),
            stock_digest, result.pre_lineage(), result.decision_digest(),
            result.digest(), observation.schedule().digest(),
            observation.station_maxima(), observation.max_tea(),
            observation.station_count(), observation.strategy_version()),
        AuditReplayTestAuthority2::reporting_observation_digest(
            request.digest(), 0, operation, motion.digest(), audit_policy.digest(),
            changed_stock_digest, result.pre_lineage(), result.decision_digest(),
            result.digest(), observation.schedule().digest(),
            observation.station_maxima(), observation.max_tea(),
            observation.station_count(), observation.strategy_version()),
        AuditReplayTestAuthority2::reporting_observation_digest(
            request.digest(), 0, operation, motion.digest(), audit_policy.digest(),
            stock_digest, foreign_lineage.digest(), result.decision_digest(),
            result.digest(), observation.schedule().digest(),
            observation.station_maxima(), observation.max_tea(),
            observation.station_count(), observation.strategy_version()),
        AuditReplayTestAuthority2::reporting_observation_digest(
            request.digest(), 0, operation, motion.digest(), audit_policy.digest(),
            stock_digest, result.pre_lineage(),
            changed_context_result.decision_digest(), result.digest(),
            observation.schedule().digest(), observation.station_maxima(),
            observation.max_tea(), observation.station_count(),
            observation.strategy_version()),
        AuditReplayTestAuthority2::reporting_observation_digest(
            request.digest(), 0, operation, motion.digest(), audit_policy.digest(),
            stock_digest, result.pre_lineage(), result.decision_digest(),
            changed_context_result.digest(), observation.schedule().digest(),
            observation.station_maxima(), observation.max_tea(),
            observation.station_count(), observation.strategy_version()),
        AuditReplayTestAuthority2::reporting_observation_digest(
            request.digest(), 0, operation, motion.digest(), audit_policy.digest(),
            stock_digest, result.pre_lineage(), result.decision_digest(),
            result.digest(), build_audit_reporting_schedule(changed_motion).digest(),
            observation.station_maxima(), observation.max_tea(),
            observation.station_count(), observation.strategy_version()),
        AuditReplayTestAuthority2::reporting_observation_digest(
            request.digest(), 0, operation, motion.digest(), audit_policy.digest(),
            stock_digest, result.pre_lineage(), result.decision_digest(),
            result.digest(), observation.schedule().digest(), changed_station_values,
            observation.max_tea(), observation.station_count(),
            observation.strategy_version()),
        AuditReplayTestAuthority2::reporting_observation_digest(
            request.digest(), 0, operation, motion.digest(), audit_policy.digest(),
            stock_digest, result.pre_lineage(), result.decision_digest(),
            result.digest(), observation.schedule().digest(),
            observation.station_maxima(), AuditReportingRadian2::build(
                observation.max_tea().value() + 0.125),
            observation.station_count(), observation.strategy_version()),
        AuditReplayTestAuthority2::reporting_observation_digest(
            request.digest(), 0, operation, motion.digest(), audit_policy.digest(),
            stock_digest, result.pre_lineage(), result.decision_digest(),
            result.digest(), observation.schedule().digest(),
            observation.station_maxima(), observation.max_tea(),
            observation.station_count() + 1, observation.strategy_version()),
        AuditReplayTestAuthority2::reporting_observation_digest(
            request.digest(), 0, operation, motion.digest(), audit_policy.digest(),
            stock_digest, result.pre_lineage(), result.decision_digest(),
            result.digest(), observation.schedule().digest(),
            observation.station_maxima(), observation.max_tea(),
            observation.station_count(), "audit-motion-reporting-schedule-v2"),
    };
    require(baseline.bytes() == observation.digest().bytes(),
            "observation differs from typed digest authority");
    for (const TeaReportingObservationDigest2& mutation : mutations)
    {
        require(mutation.bytes() != baseline.bytes(),
                "reporting observation digest ignores one typed component");
    }
}

} // namespace

void audit_reporting_observation_gate()
{
    deterministic_exact_schedule_gate();
    result_owned_observation_and_parity_gate();
    observation_component_mutation_gate();
}
