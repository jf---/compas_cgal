#include "audit_replay_2.h"

#include "audit_replay_internal_2.h"
#include "stock_exact_depletion_2.h"

#include <utility>

template <class Motion, class Apply, class Validate>
AuditLateralResult2 AuditReplay2::deplete_lateral(
    const Motion& motion,
    const AuthenticatedOperationDigest2& operation_digest,
    Apply&& apply,
    Validate&& validate)
{
    AuditReplayStorage2& storage = *storage_;
    AuditReplayPrimitiveScope2 primitive_scope(
        storage.diagnostics, AuditReplayRouteKind2::MUTATING);
    require_audit_replay_route(storage, operation_digest, motion.digest());
    try {
        const AuditDecisionWitness2 decision = certify_audit_tea_exact(
            storage.stock,
            motion,
            storage.request.policy,
            storage.request.decision_limits);
        throw_if_audit_replay_failure(
            storage, AuditReplayFailurePoint2::DECISION);

        Stock2 trial = storage.stock.clone();
        const AuditDepletionWitness2 depletion = std::forward<Apply>(apply)(
            storage.stock, trial, motion, storage.request.policy);
        throw_if_audit_replay_failure(
            storage, AuditReplayFailurePoint2::TRIAL_DEPLETION);

        if (!audit_decision_witness_is_self_consistent(
                storage.stock,
                motion,
                storage.request.policy,
                storage.request.decision_limits,
                decision)) {
            throw AuditReplayDecisionEvidenceError(
                "native decision witness failed exact replay");
        }
        if (!std::forward<Validate>(validate)(
                storage.stock,
                trial,
                motion,
                storage.request.policy,
                depletion)) {
            throw AuditReplayDepletionEvidenceError(
                "native depletion witness failed exact replay");
        }
        throw_if_audit_replay_failure(
            storage, AuditReplayFailurePoint2::WITNESS_VALIDATION);

        const AuditLineage2 post_lineage =
            AuditLineage2::transition_lateral(
                storage.request.native_request,
                storage.progress.cursor,
                storage.progress.lineage,
                operation_digest,
                motion.digest(),
                decision,
                depletion);
        throw_if_audit_replay_failure(
            storage, AuditReplayFailurePoint2::LINEAGE);

        AuditLateralResult2 result = AuditMotionResult2::lateral(
            storage.request.native_request,
            storage.progress.cursor,
            operation_digest,
            motion.digest(),
            decision,
            depletion,
            storage.progress.lineage,
            post_lineage);
        throw_if_audit_replay_failure(
            storage, AuditReplayFailurePoint2::RESULT);

        ReplayProgress2 next_progress{
            storage.progress.cursor + 1,
            post_lineage,
        };
        throw_if_audit_replay_failure(
            storage, AuditReplayFailurePoint2::PROGRESS);

        ++storage.diagnostics.instrumentation.committed_operation_count;
        storage.stock.swap(trial);
        storage.progress.swap(next_progress);
        return result;
    } catch (...) {
        ++storage.diagnostics.instrumentation.failed_transaction_count;
        throw;
    }
}

AuditLateralResult2 AuditReplay2::deplete_segment(
    const AuditSegmentMotion2& motion,
    const AuthenticatedOperationDigest2& operation_digest)
{
    return deplete_lateral(
        motion,
        operation_digest,
        apply_audit_segment_depletion_to_trial,
        validate_audit_segment_depletion_witness);
}

AuditLateralResult2 AuditReplay2::deplete_circle(
    const AuditCircleMotion2& motion,
    const AuthenticatedOperationDigest2& operation_digest)
{
    return deplete_lateral(
        motion,
        operation_digest,
        apply_audit_circle_depletion_to_trial,
        validate_audit_circle_depletion_witness);
}

AuditLateralResult2 AuditReplay2::deplete_arc(
    const AuditArcMotion2& motion,
    const AuthenticatedOperationDigest2& operation_digest)
{
    return deplete_lateral(
        motion,
        operation_digest,
        apply_audit_arc_depletion_to_trial,
        validate_audit_arc_depletion_witness);
}

AuditPlungeResult2 AuditReplay2::deplete_plunge(
    const AuditVerticalPlunge2& motion,
    const AuthenticatedOperationDigest2& operation_digest)
{
    AuditReplayStorage2& storage = *storage_;
    AuditReplayPrimitiveScope2 primitive_scope(
        storage.diagnostics, AuditReplayRouteKind2::MUTATING);
    require_audit_replay_route(storage, operation_digest, motion.digest());
    try {
        Stock2 trial = storage.stock.clone();
        const AuditDepletionWitness2 depletion =
            apply_audit_plunge_depletion_to_trial(
                storage.stock, trial, motion, storage.request.policy);
        throw_if_audit_replay_failure(
            storage, AuditReplayFailurePoint2::TRIAL_DEPLETION);

        if (!validate_audit_plunge_depletion_witness(
                storage.stock,
                trial,
                motion,
                storage.request.policy,
                depletion)) {
            throw AuditReplayDepletionEvidenceError(
                "native plunge depletion witness failed exact replay");
        }
        throw_if_audit_replay_failure(
            storage, AuditReplayFailurePoint2::WITNESS_VALIDATION);

        const AuditLineage2 post_lineage =
            AuditLineage2::transition_plunge(
                storage.request.native_request,
                storage.progress.cursor,
                storage.progress.lineage,
                operation_digest,
                motion.digest(),
                depletion);
        throw_if_audit_replay_failure(
            storage, AuditReplayFailurePoint2::LINEAGE);

        AuditPlungeResult2 result = AuditMotionResult2::plunge(
            storage.request.native_request,
            storage.progress.cursor,
            operation_digest,
            motion.digest(),
            depletion,
            storage.progress.lineage,
            post_lineage);
        throw_if_audit_replay_failure(
            storage, AuditReplayFailurePoint2::RESULT);

        ReplayProgress2 next_progress{
            storage.progress.cursor + 1,
            post_lineage,
        };
        throw_if_audit_replay_failure(
            storage, AuditReplayFailurePoint2::PROGRESS);

        ++storage.diagnostics.instrumentation.committed_operation_count;
        storage.stock.swap(trial);
        storage.progress.swap(next_progress);
        return result;
    } catch (...) {
        ++storage.diagnostics.instrumentation.failed_transaction_count;
        throw;
    }
}

AuditLateralResult2 audit_deplete_segment(
    AuditReplay2& replay,
    const AuditSegmentMotion2& motion,
    const AuthenticatedOperationDigest2& operation_digest)
{
    return replay.deplete_segment(motion, operation_digest);
}

AuditLateralResult2 audit_deplete_circle(
    AuditReplay2& replay,
    const AuditCircleMotion2& motion,
    const AuthenticatedOperationDigest2& operation_digest)
{
    return replay.deplete_circle(motion, operation_digest);
}

AuditLateralResult2 audit_deplete_arc(
    AuditReplay2& replay,
    const AuditArcMotion2& motion,
    const AuthenticatedOperationDigest2& operation_digest)
{
    return replay.deplete_arc(motion, operation_digest);
}

AuditPlungeResult2 deplete_audit_plunge(
    AuditReplay2& replay,
    const AuditVerticalPlunge2& motion,
    const AuthenticatedOperationDigest2& operation_digest)
{
    return replay.deplete_plunge(motion, operation_digest);
}
