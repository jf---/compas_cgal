#include "audit_replay_2.h"

#include "audit_replay_internal_2.h"

template <class Motion>
AuditNonEngagingResult2 AuditReplay2::record_non_engaging(
    const Motion& motion,
    const AuthenticatedOperationDigest2& operation_digest,
    AuditNonEngagingReason2 reason)
{
    AuditReplayStorage2& storage = *storage_;
    AuditReplayPrimitiveScope2 primitive_scope(
        storage.diagnostics, AuditReplayRouteKind2::NON_ENGAGING);
    require_audit_replay_route(storage, operation_digest, motion.digest());
    try {
        AuditNonEngagingResult2 result = AuditMotionResult2::non_engaging(
            storage.request.native_request,
            storage.progress.cursor,
            operation_digest,
            motion.digest(),
            reason,
            storage.progress.lineage);
        throw_if_audit_replay_failure(
            storage, AuditReplayFailurePoint2::RESULT);

        ReplayProgress2 next_progress{
            storage.progress.cursor + 1,
            storage.progress.lineage,
        };
        throw_if_audit_replay_failure(
            storage, AuditReplayFailurePoint2::PROGRESS);

        ++storage.diagnostics.instrumentation.committed_operation_count;
        storage.progress.swap(next_progress);
        return result;
    } catch (...) {
        ++storage.diagnostics.instrumentation.failed_transaction_count;
        throw;
    }
}

AuditNonEngagingResult2 AuditReplay2::record_retract(
    const AuditVerticalRetract2& motion,
    const AuthenticatedOperationDigest2& operation_digest)
{
    return record_non_engaging(
        motion,
        operation_digest,
        AuditNonEngagingReason2::VERTICAL_RETRACT);
}

AuditNonEngagingResult2 AuditReplay2::record_clearance(
    const AuditClearanceTransport2& motion,
    const AuthenticatedOperationDigest2& operation_digest)
{
    return record_non_engaging(
        motion,
        operation_digest,
        AuditNonEngagingReason2::CLEARANCE_TRANSPORT);
}

AuditNonEngagingResult2 record_audit_retract(
    AuditReplay2& replay,
    const AuditVerticalRetract2& motion,
    const AuthenticatedOperationDigest2& operation_digest)
{
    return replay.record_retract(motion, operation_digest);
}

AuditNonEngagingResult2 record_audit_clearance(
    AuditReplay2& replay,
    const AuditClearanceTransport2& motion,
    const AuthenticatedOperationDigest2& operation_digest)
{
    return replay.record_clearance(motion, operation_digest);
}
