#include "audit_replay_fixtures_2.h"

#include "audit_motion_result_2.h"

#include <initializer_list>
#include <string>
#include <utility>
#include <vector>

using namespace audit_replay_fixtures;

namespace
{

using FailureExpectation2 = std::pair<AuditReplayFailurePoint2, std::size_t>;

template <class Motion, class Invoke>
void mutating_failure_route_gate(const char* label, const Motion& motion, Invoke invoke,
                                 std::initializer_list<FailureExpectation2> failures)
{
    // Kills a route-specific throwing stage that commits authority or clones at
    // the wrong point in the transaction chronology.
    const compas::RowMatrixXd boundary = rectangle(-5.0, -5.0, 5.0, 5.0);
    const AuditPolicy2 audit_policy = policy();
    const AuditDecisionLimits2 limits = decision_limits();
    const AuthenticatedOperationDigest2 digest = operation_digest(label);

    AuditReplay2 clean = one_motion_replay(boundary, audit_policy, limits, motion, digest);
    const auto clean_result = invoke(clean, motion, digest);

    for (const auto& [failure, expected_failed_clones] : failures)
    {
        AuditReplay2 replay = one_motion_replay(boundary, audit_policy, limits, motion, digest);
        const AuditReplayInspection2 before = inspect_audit_replay_for_test(replay);
        inject_audit_replay_failure_for_test(replay, failure);
        require_throws<AuditReplayInjectedFailure>(
            [&] { static_cast<void>(invoke(replay, motion, digest)); },
            "controlled mutating-route failure did not escape");
        require(inspect_audit_replay_for_test(replay) == before,
                "failed mutating route changed authority state");

        const AuditReplayInstrumentation2 failed = audit_replay_instrumentation_for_test(replay);
        require(failed.mutating_clone_count == expected_failed_clones &&
                    failed.mutating_swap_count == 0 &&
                    failed.non_engaging_clone_count == 0 && failed.committed_operation_count == 0 &&
                    failed.failed_transaction_count == 1,
                "failed mutating route reported false clone/commit counts");

        clear_audit_replay_failure_for_test(replay);
        const auto retried = invoke(replay, motion, digest);
        require(retried.digest().bytes() == clean_result.digest().bytes(),
                "retry after mutating-route failure diverged from clean replay");
        if constexpr (requires { retried.reporting_observation(); })
        {
            require(retried.reporting_observation().digest().bytes() ==
                        clean_result.reporting_observation().digest().bytes(),
                    "retry after REPORTING failure changed observation identity");
        }
        const AuditReplayInstrumentation2 retried_metrics =
            audit_replay_instrumentation_for_test(replay);
        require(retried_metrics.mutating_clone_count == expected_failed_clones + 1 &&
                    retried_metrics.mutating_swap_count == 1 &&
                    retried_metrics.non_engaging_clone_count == 0 &&
                    retried_metrics.committed_operation_count == 1 &&
                    retried_metrics.failed_transaction_count == 1,
                "mutating-route retry did not perform exactly one clean commit");
    }
}

template <class Motion, class Invoke>
void non_engaging_failure_route_gate(const char* label, const Motion& motion, Invoke invoke)
{
    // Kills stock cloning or progress mutation before a non-engaging result is
    // completely constructed.
    const compas::RowMatrixXd boundary = rectangle(-5.0, -5.0, 5.0, 5.0);
    const AuditPolicy2 audit_policy = policy();
    const AuditDecisionLimits2 limits = decision_limits();
    const AuthenticatedOperationDigest2 digest = operation_digest(label);

    AuditReplay2 clean = one_motion_replay(boundary, audit_policy, limits, motion, digest);
    const auto clean_result = invoke(clean, motion, digest);

    for (const AuditReplayFailurePoint2 failure : {
             AuditReplayFailurePoint2::RESULT,
             AuditReplayFailurePoint2::PROGRESS,
         })
    {
        AuditReplay2 replay = one_motion_replay(boundary, audit_policy, limits, motion, digest);
        const AuditReplayInspection2 before = inspect_audit_replay_for_test(replay);
        inject_audit_replay_failure_for_test(replay, failure);
        require_throws<AuditReplayInjectedFailure>(
            [&] { static_cast<void>(invoke(replay, motion, digest)); },
            "controlled non-engaging failure did not escape");
        require(inspect_audit_replay_for_test(replay) == before,
                "failed non-engaging route changed authority state");

        const AuditReplayInstrumentation2 failed = audit_replay_instrumentation_for_test(replay);
        require(failed.mutating_clone_count == 0 && failed.non_engaging_clone_count == 0 &&
                    failed.mutating_swap_count == 0 &&
                    failed.committed_operation_count == 0 && failed.failed_transaction_count == 1,
                "failed non-engaging route cloned stock or committed progress");

        clear_audit_replay_failure_for_test(replay);
        const auto retried = invoke(replay, motion, digest);
        require(retried.digest().bytes() == clean_result.digest().bytes(),
                "retry after non-engaging failure diverged from clean replay");
        const AuditReplayInstrumentation2 retried_metrics =
            audit_replay_instrumentation_for_test(replay);
        require(retried_metrics.mutating_clone_count == 0 &&
                    retried_metrics.mutating_swap_count == 0 &&
                    retried_metrics.non_engaging_clone_count == 0 &&
                    retried_metrics.committed_operation_count == 1 &&
                    retried_metrics.failed_transaction_count == 1,
                "non-engaging retry did not perform one clone-free commit");
    }
}

void wrong_identity_preserves_state_gate()
{
    // Kills cursor advancement, cloning, or stock mutation on foreign identity.
    const compas::RowMatrixXd boundary = rectangle(-5.0, -5.0, 5.0, 5.0);
    const AuditPolicy2 audit_policy = policy();
    const AuditDecisionLimits2 limits = decision_limits();
    const AuditSegmentMotion2 motion = segment();
    const AuthenticatedOperationDigest2 expected = operation_digest("segment");
    AuditReplay2 replay = one_motion_replay(boundary, audit_policy, limits, motion, expected);
    const AuditReplayInspection2 initial = inspect_audit_replay_for_test(replay);

    require_throws<AuditReplayOperationIdentityError>(
        [&]
        {
            static_cast<void>(
                audit_deplete_segment(replay, segment(1.0), operation_digest("foreign")));
        },
        "motion validation ran before authenticated-operation validation");
    require(inspect_audit_replay_for_test(replay) == initial,
            "combined foreign identities mutated replay state");

    require_throws<AuditReplayMotionIdentityError>(
        [&] { static_cast<void>(audit_deplete_segment(replay, segment(1.0), expected)); },
        "wrong actual motion identity was accepted");
    require(inspect_audit_replay_for_test(replay) == initial,
            "wrong motion identity mutated replay state");

    const AuditReplayInstrumentation2 metrics = audit_replay_instrumentation_for_test(replay);
    require(metrics.mutating_clone_count == 0 && metrics.mutating_swap_count == 0 &&
                metrics.committed_operation_count == 0,
            "identity preflight cloned stock or committed an operation");
}

void mutating_failure_matrix_gate()
{
    const std::initializer_list<FailureExpectation2> lateral_failures{
        {AuditReplayFailurePoint2::DECISION, 0},
        {AuditReplayFailurePoint2::TRIAL_DEPLETION, 1},
        {AuditReplayFailurePoint2::WITNESS_VALIDATION, 1},
        {AuditReplayFailurePoint2::LINEAGE, 1},
        {AuditReplayFailurePoint2::RESULT, 1},
        {AuditReplayFailurePoint2::REPORTING, 1},
        {AuditReplayFailurePoint2::PROGRESS, 1},
    };
    mutating_failure_route_gate(
        "segment-failure-matrix", segment(),
        [](AuditReplay2& replay, const AuditSegmentMotion2& motion,
           const AuthenticatedOperationDigest2& digest)
        { return audit_deplete_segment(replay, motion, digest); },
        lateral_failures);
    mutating_failure_route_gate(
        "circle-failure-matrix", circle(),
        [](AuditReplay2& replay, const AuditCircleMotion2& motion,
           const AuthenticatedOperationDigest2& digest)
        { return audit_deplete_circle(replay, motion, digest); },
        lateral_failures);
    mutating_failure_route_gate(
        "arc-failure-matrix", arc(),
        [](AuditReplay2& replay, const AuditArcMotion2& motion,
           const AuthenticatedOperationDigest2& digest)
        { return audit_deplete_arc(replay, motion, digest); },
        lateral_failures);

    const std::initializer_list<FailureExpectation2> plunge_failures{
        {AuditReplayFailurePoint2::TRIAL_DEPLETION, 1},
        {AuditReplayFailurePoint2::WITNESS_VALIDATION, 1},
        {AuditReplayFailurePoint2::LINEAGE, 1},
        {AuditReplayFailurePoint2::RESULT, 1},
        {AuditReplayFailurePoint2::PROGRESS, 1},
    };
    mutating_failure_route_gate(
        "plunge-failure-matrix", plunge(),
        [](AuditReplay2& replay, const AuditVerticalPlunge2& motion,
           const AuthenticatedOperationDigest2& digest)
        { return deplete_audit_plunge(replay, motion, digest); },
        plunge_failures);
}

void non_engaging_failure_matrix_gate()
{
    non_engaging_failure_route_gate("retract-failure-matrix", retract(),
                                    [](AuditReplay2& replay, const AuditVerticalRetract2& motion,
                                       const AuthenticatedOperationDigest2& digest)
                                    { return record_audit_retract(replay, motion, digest); });
    non_engaging_failure_route_gate("clearance-failure-matrix", clearance(),
                                    [](AuditReplay2& replay, const AuditClearanceTransport2& motion,
                                       const AuthenticatedOperationDigest2& digest)
                                    { return record_audit_clearance(replay, motion, digest); });
}

void out_of_order_duplicate_and_exhaustion_gate()
{
    // Kills cursor-independent dispatch, duplicate consumption, or identity
    // checks that run before the exhausted-state check.
    const compas::RowMatrixXd boundary = rectangle(-5.0, -5.0, 5.0, 5.0);
    const AuditPolicy2 audit_policy = policy();
    const AuditDecisionLimits2 limits = decision_limits();
    const AuditSegmentMotion2 first = segment(0.0);
    const AuditSegmentMotion2 second = segment(1.0);
    const AuthenticatedOperationDigest2 first_digest = operation_digest("first");
    const AuthenticatedOperationDigest2 second_digest = operation_digest("second");
    const AuditNativeRequestIdentity2 request =
        AuditNativeRequestIdentity2::build(AuditNativeStockIdentity2::build(boundary, {}),
                                           audit_policy, limits, {first.digest(), second.digest()});
    AuditReplay2 replay = AuditReplay2::build(boundary, {}, input_digest(), request, audit_policy,
                                              limits, {first_digest, second_digest});
    const AuditReplayInspection2 initial = inspect_audit_replay_for_test(replay);

    require_throws<AuditReplayOperationIdentityError>(
        [&] { static_cast<void>(audit_deplete_segment(replay, second, second_digest)); },
        "out-of-order operation was accepted");
    require_throws<AuditReplayMotionIdentityError>(
        [&] { static_cast<void>(audit_deplete_segment(replay, second, first_digest)); },
        "out-of-order motion was accepted with the current operation digest");
    require(inspect_audit_replay_for_test(replay) == initial,
            "out-of-order preflight changed replay state");

    static_cast<void>(audit_deplete_segment(replay, first, first_digest));
    const AuditReplayInspection2 after_first = inspect_audit_replay_for_test(replay);
    require_throws<AuditReplayOperationIdentityError>(
        [&] { static_cast<void>(audit_deplete_segment(replay, first, first_digest)); },
        "duplicate operation consumption was accepted");
    require(inspect_audit_replay_for_test(replay) == after_first,
            "duplicate operation attempt changed replay state");

    static_cast<void>(audit_deplete_segment(replay, second, second_digest));
    const AuditReplayInspection2 exhausted = inspect_audit_replay_for_test(replay);
    require_throws<AuditReplayOperationExhaustedError>(
        [&]
        {
            static_cast<void>(
                audit_deplete_segment(replay, first, operation_digest("foreign-after-end")));
        },
        "identity validation ran before exhausted-state validation");
    require(inspect_audit_replay_for_test(replay) == exhausted,
            "exhausted replay attempt changed replay state");

    const AuditReplayInstrumentation2 metrics = audit_replay_instrumentation_for_test(replay);
    require(metrics.mutating_clone_count == 2 && metrics.mutating_swap_count == 2 &&
                metrics.non_engaging_clone_count == 0 &&
                metrics.committed_operation_count == 2,
            "preflight rejection changed clone or commit counts");
}

void closed_motion_domain_gate()
{
    // Kills a missing typed route, a no-op mutating route, or non-engaging clone.
    const compas::RowMatrixXd boundary = rectangle(-5.0, -5.0, 5.0, 5.0);
    const AuditPolicy2 audit_policy = policy();
    const AuditDecisionLimits2 limits = decision_limits();
    const AuditSegmentMotion2 segment_motion = segment();
    const AuditArcMotion2 arc_motion = arc();
    const AuditCircleMotion2 circle_motion = circle();
    const AuditVerticalPlunge2 plunge_motion = plunge();
    const AuditVerticalRetract2 retract_motion = retract();
    const AuditClearanceTransport2 clearance_motion = clearance();
    const std::vector<NativeMotionDigest2> motion_digests{
        segment_motion.digest(), arc_motion.digest(),     circle_motion.digest(),
        plunge_motion.digest(),  retract_motion.digest(), clearance_motion.digest(),
    };
    std::vector<AuthenticatedOperationDigest2> operation_digests;
    for (const std::string& name : {"segment", "arc", "circle", "plunge", "retract", "clearance"})
    {
        operation_digests.push_back(operation_digest(name));
    }
    const AuditNativeRequestIdentity2 request = AuditNativeRequestIdentity2::build(
        AuditNativeStockIdentity2::build(boundary, {}), audit_policy, limits, motion_digests);
    AuditReplay2 replay = AuditReplay2::build(boundary, {}, input_digest(), request, audit_policy,
                                              limits, operation_digests);

    AuditReplayInspection2 before = inspect_audit_replay_for_test(replay);
    const AuditLateralResult2 segment_result =
        audit_deplete_segment(replay, segment_motion, operation_digests[0]);
    AuditReplayInspection2 after = inspect_audit_replay_for_test(replay);
    require(after.stock_digest != before.stock_digest,
            "intersecting segment did not change exact stock state");
    before = after;
    const AuditLateralResult2 arc_result =
        audit_deplete_arc(replay, arc_motion, operation_digests[1]);
    after = inspect_audit_replay_for_test(replay);
    require(after.stock_digest != before.stock_digest,
            "intersecting arc did not change exact stock state");
    before = after;
    const AuditLateralResult2 circle_result =
        audit_deplete_circle(replay, circle_motion, operation_digests[2]);
    after = inspect_audit_replay_for_test(replay);
    require(after.stock_digest != before.stock_digest,
            "intersecting circle did not change exact stock state");
    before = after;
    const AuditPlungeResult2 plunge_result =
        deplete_audit_plunge(replay, plunge_motion, operation_digests[3]);
    after = inspect_audit_replay_for_test(replay);
    require(after.stock_digest != before.stock_digest,
            "intersecting plunge did not change exact stock state");
    before = after;
    const AuditNonEngagingResult2 retract_result =
        record_audit_retract(replay, retract_motion, operation_digests[4]);
    after = inspect_audit_replay_for_test(replay);
    require(after.stock_digest == before.stock_digest, "retract changed exact stock state");
    before = after;
    const AuditNonEngagingResult2 clearance_result =
        record_audit_clearance(replay, clearance_motion, operation_digests[5]);
    after = inspect_audit_replay_for_test(replay);
    require(after.stock_digest == before.stock_digest, "clearance changed exact stock state");

    require(segment_result.depletion_witness_digest().bytes().size() == 32 &&
                arc_result.depletion_witness_digest().bytes().size() == 32 &&
                circle_result.depletion_witness_digest().bytes().size() == 32 &&
                plunge_result.depletion_witness_digest().bytes().size() == 32,
            "mutating route omitted full depletion evidence");
    require(retract_result.pre_lineage().bytes() == retract_result.post_lineage().bytes() &&
                clearance_result.pre_lineage().bytes() == clearance_result.post_lineage().bytes(),
            "non-engaging route advanced lineage");
    const AuditReplayInstrumentation2 metrics = audit_replay_instrumentation_for_test(replay);
    require(metrics.mutating_clone_count == 4 && metrics.mutating_swap_count == 4 &&
                metrics.non_engaging_clone_count == 0 &&
                metrics.committed_operation_count == 6,
            "closed-domain replay clone/commit counts are false");
}

} // namespace

void audit_replay_transaction_state_gate()
{
    wrong_identity_preserves_state_gate();
    mutating_failure_matrix_gate();
    non_engaging_failure_matrix_gate();
    out_of_order_duplicate_and_exhaustion_gate();
    closed_motion_domain_gate();
}
