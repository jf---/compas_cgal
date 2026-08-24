#include "audit_replay_2.h"

#include "audit_replay_internal_2.h"

AuditReplayCompletion2 AuditReplay2::finish()
{
    AuditReplayStorage2& storage = *storage_;
    if (storage.lifecycle.finalized) {
        throw AuditReplayFinalizedError("audit replay is already finalized");
    }
    if (storage.progress.cursor != storage.request.operation_digests.size()) {
        throw AuditReplayIncompleteError(
            "audit replay cannot finalize before every operation is consumed");
    }

    try {
        AuditReplayCompletion2 completion = AuditMotionResult2::completion(
            storage.request.native_request,
            storage.progress.cursor,
            storage.progress.lineage);
        throw_if_audit_replay_failure(
            storage, AuditReplayFailurePoint2::FINALIZATION);

        ++storage.diagnostics.instrumentation.committed_finalization_count;
        storage.lifecycle.finalized = true;
        return completion;
    } catch (...) {
        ++storage.diagnostics.instrumentation.failed_transaction_count;
        ++storage.diagnostics.instrumentation.failed_finalization_count;
        throw;
    }
}

AuditReplayCompletion2 finish_audit_replay(AuditReplay2& replay)
{
    return replay.finish();
}
