#include "audit_replay_fixtures_2.h"

#include "audit_motion_result_2.h"

#include <cstddef>
#include <string>
#include <vector>

using namespace audit_replay_fixtures;

namespace
{

struct SixMotionValues2
{
    AuditSegmentMotion2 segment_motion;
    AuditArcMotion2 arc_motion;
    AuditCircleMotion2 circle_motion;
    AuditVerticalPlunge2 plunge_motion;
    AuditVerticalRetract2 retract_motion;
    AuditClearanceTransport2 clearance_motion;
    std::vector<AuthenticatedOperationDigest2> operation_digests;
};

SixMotionValues2 six_motion_values()
{
    std::vector<AuthenticatedOperationDigest2> digests;
    for (const std::string& name : {
             "finish-segment",
             "finish-arc",
             "finish-circle",
             "finish-plunge",
             "finish-retract",
             "finish-clearance",
         })
    {
        digests.push_back(operation_digest(name));
    }
    return {
        segment(), arc(), circle(), plunge(), retract(), clearance(), std::move(digests),
    };
}

AuditNativeRequestIdentity2 six_motion_request(const compas::RowMatrixXd& boundary,
                                               const AuditPolicy2& audit_policy,
                                               const AuditDecisionLimits2& limits,
                                               const SixMotionValues2& values)
{
    return AuditNativeRequestIdentity2::build(AuditNativeStockIdentity2::build(boundary, {}),
                                              audit_policy, limits,
                                              {
                                                  values.segment_motion.digest(),
                                                  values.arc_motion.digest(),
                                                  values.circle_motion.digest(),
                                                  values.plunge_motion.digest(),
                                                  values.retract_motion.digest(),
                                                  values.clearance_motion.digest(),
                                              });
}

AuditReplay2 six_motion_replay(const compas::RowMatrixXd& boundary,
                               const AuditPolicy2& audit_policy, const AuditDecisionLimits2& limits,
                               const SixMotionValues2& values)
{
    return AuditReplay2::build(boundary, {}, input_digest(),
                               six_motion_request(boundary, audit_policy, limits, values),
                               audit_policy, limits, values.operation_digests);
}

void consume_prefix(AuditReplay2& replay, const SixMotionValues2& values, std::size_t prefix)
{
    if (prefix >= 1)
    {
        static_cast<void>(
            audit_deplete_segment(replay, values.segment_motion, values.operation_digests[0]));
    }
    if (prefix >= 2)
    {
        static_cast<void>(
            audit_deplete_arc(replay, values.arc_motion, values.operation_digests[1]));
    }
    if (prefix >= 3)
    {
        static_cast<void>(
            audit_deplete_circle(replay, values.circle_motion, values.operation_digests[2]));
    }
    if (prefix >= 4)
    {
        static_cast<void>(
            deplete_audit_plunge(replay, values.plunge_motion, values.operation_digests[3]));
    }
    if (prefix >= 5)
    {
        static_cast<void>(
            record_audit_retract(replay, values.retract_motion, values.operation_digests[4]));
    }
    if (prefix >= 6)
    {
        static_cast<void>(
            record_audit_clearance(replay, values.clearance_motion, values.operation_digests[5]));
    }
}

void every_incomplete_prefix_gate()
{
    // Kills a finalizer that accepts any proper prefix of the authenticated
    // operation sequence or changes state while rejecting it.
    const compas::RowMatrixXd boundary = rectangle(-5.0, -5.0, 5.0, 5.0);
    const AuditPolicy2 audit_policy = policy();
    const AuditDecisionLimits2 limits = decision_limits();
    const SixMotionValues2 values = six_motion_values();

    for (std::size_t prefix = 0; prefix < 6; ++prefix)
    {
        AuditReplay2 replay = six_motion_replay(boundary, audit_policy, limits, values);
        consume_prefix(replay, values, prefix);
        const AuditReplayInspection2 before = inspect_audit_replay_for_test(replay);
        require_throws<AuditReplayIncompleteError>(
            [&] { static_cast<void>(finish_audit_replay(replay)); },
            "proper-prefix finalization was accepted");
        require(inspect_audit_replay_for_test(replay) == before,
                "proper-prefix finalization changed replay state");
        const AuditReplayInstrumentation2 metrics = audit_replay_instrumentation_for_test(replay);
        require(metrics.mutating_swap_count == (prefix < 4 ? prefix : 4) &&
                    metrics.committed_operation_count == prefix &&
                    metrics.committed_finalization_count == 0 &&
                    metrics.failed_finalization_count == 0,
                "incomplete-prefix rejection reported false counters");
    }
}

void post_finalization_routes_fail_closed_gate(AuditReplay2& replay, const SixMotionValues2& values)
{
    // Kills route-specific dispatch that checks cursor, digest, or motion before
    // the terminal finalized state.
    const AuditReplayInspection2 finalized = inspect_audit_replay_for_test(replay);
    require_throws<AuditReplayFinalizedError>(
        [&]
        {
            static_cast<void>(
                audit_deplete_segment(replay, values.segment_motion, values.operation_digests[0]));
        },
        "post-finalization segment call was accepted");
    require_throws<AuditReplayFinalizedError>(
        [&]
        {
            static_cast<void>(
                audit_deplete_arc(replay, values.arc_motion, values.operation_digests[1]));
        },
        "post-finalization arc call was accepted");
    require_throws<AuditReplayFinalizedError>(
        [&]
        {
            static_cast<void>(
                audit_deplete_circle(replay, values.circle_motion, values.operation_digests[2]));
        },
        "post-finalization circle call was accepted");
    require_throws<AuditReplayFinalizedError>(
        [&]
        {
            static_cast<void>(
                deplete_audit_plunge(replay, values.plunge_motion, values.operation_digests[3]));
        },
        "post-finalization plunge call was accepted");
    require_throws<AuditReplayFinalizedError>(
        [&]
        {
            static_cast<void>(
                record_audit_retract(replay, values.retract_motion, values.operation_digests[4]));
        },
        "post-finalization retract call was accepted");
    require_throws<AuditReplayFinalizedError>(
        [&]
        {
            static_cast<void>(record_audit_clearance(replay, values.clearance_motion,
                                                     values.operation_digests[5]));
        },
        "post-finalization clearance call was accepted");
    require(inspect_audit_replay_for_test(replay) == finalized,
            "post-finalization route changed replay state");
}

void finalization_failure_retry_and_terminal_gate()
{
    // Kills finalization that commits before completion construction, cannot be
    // retried, or exposes a completion inconsistent with a clean replay.
    const compas::RowMatrixXd boundary = rectangle(-5.0, -5.0, 5.0, 5.0);
    const AuditPolicy2 audit_policy = policy();
    const AuditDecisionLimits2 limits = decision_limits();
    const SixMotionValues2 values = six_motion_values();
    const AuditNativeRequestIdentity2 request =
        six_motion_request(boundary, audit_policy, limits, values);

    AuditReplay2 clean = six_motion_replay(boundary, audit_policy, limits, values);
    consume_prefix(clean, values, 6);
    const AuditReplayCompletion2 clean_completion = finish_audit_replay(clean);

    AuditReplay2 replay = six_motion_replay(boundary, audit_policy, limits, values);
    consume_prefix(replay, values, 6);
    inject_audit_replay_failure_for_test(replay, AuditReplayFailurePoint2::FINALIZATION);
    const AuditReplayInspection2 before = inspect_audit_replay_for_test(replay);
    require_throws<AuditReplayInjectedFailure>([&]
                                               { static_cast<void>(finish_audit_replay(replay)); },
                                               "controlled finalization failure did not escape");
    require(inspect_audit_replay_for_test(replay) == before,
            "failed finalization changed replay state");
    const AuditReplayInstrumentation2 failed = audit_replay_instrumentation_for_test(replay);
    require(failed.mutating_clone_count == 4 && failed.mutating_swap_count == 4 &&
                failed.non_engaging_clone_count == 0 &&
                failed.committed_operation_count == 6 && failed.failed_transaction_count == 1 &&
                failed.failed_finalization_count == 1 && failed.committed_finalization_count == 0,
            "failed finalization reported false transaction counters");

    clear_audit_replay_failure_for_test(replay);
    const AuditReplayCompletion2 completion = finish_audit_replay(replay);
    const AuditLineage2 expected_seed =
        AuditReplayTestAuthority2::seed_lineage(input_digest(), request);
    require(completion.digest().bytes() == clean_completion.digest().bytes(),
            "retry after finalization failure diverged from clean completion");
    require(completion.input_digest().bytes() == input_digest().bytes() &&
                completion.seed_lineage().bytes() == expected_seed.digest().bytes() &&
                completion.operation_count() == 6 &&
                completion.request_digest().bytes() == request.digest().bytes() &&
                completion.terminal_lineage().bytes() ==
                    clean_completion.terminal_lineage().bytes() &&
                completion.digest().bytes().size() == 32,
            "completion lost input, seed, count, request, terminal lineage, or digest");

    const AuditReplayInstrumentation2 committed = audit_replay_instrumentation_for_test(replay);
    require(committed.mutating_swap_count == 4 &&
                committed.committed_operation_count == 6 &&
                committed.failed_transaction_count == 1 &&
                committed.failed_finalization_count == 1 &&
                committed.committed_finalization_count == 1,
            "finalization retry did not commit exactly once");
    const AuditReplayInspection2 finalized = inspect_audit_replay_for_test(replay);
    require_throws<AuditReplayFinalizedError>([&]
                                              { static_cast<void>(finish_audit_replay(replay)); },
                                              "duplicate finalization was accepted");
    require(inspect_audit_replay_for_test(replay) == finalized,
            "duplicate finalization changed replay state");
    post_finalization_routes_fail_closed_gate(replay, values);

    const AuditReplayInstrumentation2 terminal = audit_replay_instrumentation_for_test(replay);
    require(terminal.mutating_clone_count == 4 && terminal.mutating_swap_count == 4 &&
                terminal.non_engaging_clone_count == 0 &&
                terminal.committed_operation_count == 6 && terminal.failed_transaction_count == 1 &&
                terminal.failed_finalization_count == 1 &&
                terminal.committed_finalization_count == 1,
            "terminal rejections changed clone, commit, or failure counters");
}

} // namespace

void audit_replay_finalization_gate()
{
    every_incomplete_prefix_gate();
    finalization_failure_retry_and_terminal_gate();
}
