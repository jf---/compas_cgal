#include "audit_replay_2.h"

#include "audit_replay_internal_2.h"
#include "audit_stock_identity_2.h"
#include "audit_stock_state_identity_2.h"

#include <memory>
#include <utility>

namespace {

thread_local AuditReplayDiagnostics2* active_replay_diagnostics = nullptr;
thread_local AuditReplayRouteKind2 active_replay_route_kind =
    AuditReplayRouteKind2::MUTATING;
thread_local AuditReplayConstructionInstrumentation2
    replay_construction_instrumentation;

bool same_bytes(const std::string& left, const std::string& right) noexcept
{
    return left == right;
}

} // namespace

AuditReplayPrimitiveScope2::AuditReplayPrimitiveScope2(
    AuditReplayDiagnostics2& diagnostics,
    AuditReplayRouteKind2 route_kind) noexcept
    : previous_(active_replay_diagnostics),
      previous_route_kind_(active_replay_route_kind)
{
    active_replay_diagnostics = &diagnostics;
    active_replay_route_kind = route_kind;
}

AuditReplayPrimitiveScope2::~AuditReplayPrimitiveScope2()
{
    active_replay_diagnostics = previous_;
    active_replay_route_kind = previous_route_kind_;
}

void note_audit_replay_stock_clone_for_test() noexcept
{
    ++replay_construction_instrumentation.stock_clone_count;
    if (active_replay_diagnostics != nullptr) {
        if (active_replay_route_kind == AuditReplayRouteKind2::MUTATING) {
            ++active_replay_diagnostics->instrumentation.mutating_clone_count;
        } else {
            ++active_replay_diagnostics->instrumentation
                  .non_engaging_clone_count;
        }
    }
}

void note_audit_replay_stock_swap_for_test() noexcept
{
    if (active_replay_diagnostics != nullptr) {
        ++active_replay_diagnostics->instrumentation.mutating_swap_count;
    }
}

AuditReplay2::AuditReplay2(
    AuditReplayRequestState2 request,
    Stock2 stock,
    ReplayProgress2 progress)
    : storage_(std::make_unique<AuditReplayStorage2>(
          AuditReplayStorage2{
              std::move(request),
              std::move(stock),
              std::move(progress),
              AuditReplayLifecycle2{false},
              AuditReplayDiagnostics2{},
          }))
{
}

AuditReplay2::AuditReplay2(AuditReplay2&&) noexcept = default;
AuditReplay2& AuditReplay2::operator=(AuditReplay2&&) noexcept = default;
AuditReplay2::~AuditReplay2() = default;

AuditReplay2 AuditReplay2::build(
    Eigen::Ref<const compas::RowMatrixXd> boundary,
    const std::vector<compas::RowMatrixXd>& holes,
    const AuditInputDigest2& input_digest,
    const AuditNativeRequestIdentity2& native_request,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& decision_limits,
    std::vector<AuthenticatedOperationDigest2>
        authenticated_operation_digests)
{
    if (native_request.motion_digests_.size()
        != authenticated_operation_digests.size()) {
        throw AuditReplayCardinalityError(
            "authenticated operations and native motions must have equal cardinality");
    }

    const AuditNativeStockIdentity2 actual_stock =
        AuditNativeStockIdentity2::build(boundary, holes);
    const AuditNativeRequestIdentity2 actual_request =
        AuditNativeRequestIdentity2::build(
            actual_stock,
            policy,
            decision_limits,
            native_request.motion_digests_);
    if (!same_bytes(
            actual_request.canonical_bytes(),
            native_request.canonical_bytes())
        || !same_bytes(
            actual_request.digest().bytes(), native_request.digest().bytes())) {
        throw AuditReplayRequestIdentityError(
            "actual stock, policy, limits, or motion sequence does not match native request");
    }

    Stock2 stock(boundary, holes);
    AuditLineage2 seed = AuditLineage2::seed(input_digest, native_request);
    return AuditReplay2(
        AuditReplayRequestState2{
            native_request,
            policy,
            decision_limits,
            native_request.motion_digests_,
            std::move(authenticated_operation_digests),
        },
        std::move(stock),
        ReplayProgress2{0, std::move(seed)});
}

AuditReplay2 begin_audit_replay(
    Eigen::Ref<const compas::RowMatrixXd> boundary,
    const std::vector<compas::RowMatrixXd>& holes,
    std::string input_digest_bytes,
    const AuditNativeRequestIdentity2& native_request,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& decision_limits,
    std::vector<std::string> authenticated_operation_digest_bytes)
{
    std::vector<AuthenticatedOperationDigest2> operation_digests;
    operation_digests.reserve(authenticated_operation_digest_bytes.size());
    for (std::string& bytes : authenticated_operation_digest_bytes) {
        operation_digests.push_back(
            AuthenticatedOperationDigestAuthority2::from_external_bytes(
                std::move(bytes)));
    }
    return AuditReplay2::build(
        boundary,
        holes,
        AuditInputDigestAuthority2::from_external_bytes(
            std::move(input_digest_bytes)),
        native_request,
        policy,
        decision_limits,
        std::move(operation_digests));
}

void require_audit_replay_route(
    const AuditReplayStorage2& storage,
    const AuthenticatedOperationDigest2& operation_digest,
    const NativeMotionDigest2& motion_digest)
{
    if (storage.lifecycle.finalized) {
        throw AuditReplayFinalizedError("audit replay is already finalized");
    }
    if (storage.progress.cursor >= storage.request.operation_digests.size()) {
        throw AuditReplayOperationExhaustedError(
            "audit replay has consumed every authenticated operation");
    }
    if (!same_bytes(
            storage.request.operation_digests[storage.progress.cursor].bytes(),
            operation_digest.bytes())) {
        throw AuditReplayOperationIdentityError(
            "authenticated operation does not match replay cursor");
    }
    if (storage.progress.cursor >= storage.request.motion_digests.size()
        || !same_bytes(
            storage.request.motion_digests[storage.progress.cursor].bytes(),
            motion_digest.bytes())) {
        throw AuditReplayMotionIdentityError(
            "native motion does not match replay cursor");
    }
}

void throw_if_audit_replay_failure(
    const AuditReplayStorage2& storage,
    AuditReplayFailurePoint2 point)
{
    if (storage.diagnostics.injected_failure == point) {
        throw AuditReplayInjectedFailure(
            "controlled audit replay failure was injected");
    }
}

void inject_audit_replay_failure_for_test(
    AuditReplay2& replay,
    AuditReplayFailurePoint2 point) noexcept
{
    replay.storage_->diagnostics.injected_failure = point;
}

void clear_audit_replay_failure_for_test(AuditReplay2& replay) noexcept
{
    replay.storage_->diagnostics.injected_failure.reset();
}

AuditReplayInspection2 inspect_audit_replay_for_test(const AuditReplay2& replay)
{
    const AuditReplayStorage2& storage = *replay.storage_;
    return {
        storage.progress.cursor,
        AuditStockStateIdentity2::build(storage.stock).digest().bytes(),
        storage.progress.lineage.digest().bytes(),
        storage.lifecycle.finalized,
    };
}

AuditReplayInstrumentation2 audit_replay_instrumentation_for_test(
    const AuditReplay2& replay) noexcept
{
    return replay.storage_->diagnostics.instrumentation;
}

void reset_audit_replay_construction_instrumentation_for_test() noexcept
{
    replay_construction_instrumentation = {};
}

AuditReplayConstructionInstrumentation2
audit_replay_construction_instrumentation_for_test() noexcept
{
    return replay_construction_instrumentation;
}

AuditLineage2 AuditReplayTestAuthority2::seed_lineage(
    const AuditInputDigest2& input_digest,
    const AuditNativeRequestIdentity2& request)
{
    return AuditLineage2::seed(input_digest, request);
}

AuditLineage2 AuditReplayTestAuthority2::transition_lineage(
    const AuditNativeRequestIdentity2& request,
    std::size_t cursor,
    const AuditLineage2& pre_lineage,
    const AuthenticatedOperationDigest2& operation_digest,
    const NativeMotionDigest2& motion_digest,
    const AuditDecisionWitness2& decision,
    const AuditDepletionWitness2& depletion)
{
    return AuditLineage2::transition_lateral(
        request,
        cursor,
        pre_lineage,
        operation_digest,
        motion_digest,
        decision,
        depletion);
}

AuditLineage2 AuditReplayTestAuthority2::transition_plunge_lineage(
    const AuditNativeRequestIdentity2& request,
    std::size_t cursor,
    const AuditLineage2& pre_lineage,
    const AuthenticatedOperationDigest2& operation_digest,
    const NativeMotionDigest2& motion_digest,
    const AuditDepletionWitness2& depletion)
{
    return AuditLineage2::transition_plunge(
        request,
        cursor,
        pre_lineage,
        operation_digest,
        motion_digest,
        depletion);
}

AuditResultDigest2 AuditReplayTestAuthority2::lateral_result_digest(
    const AuditNativeRequestIdentity2& request,
    std::size_t cursor,
    const AuthenticatedOperationDigest2& operation_digest,
    const NativeMotionDigest2& motion_digest,
    const AuditDecisionWitness2& decision,
    const AuditDepletionWitness2& depletion,
    const StockLineageDigest2& pre_lineage,
    const StockLineageDigest2& post_lineage)
{
    return AuditMotionResult2::lateral_digest(
        request,
        cursor,
        operation_digest,
        motion_digest,
        decision,
        depletion,
        pre_lineage,
        post_lineage);
}

AuditResultDigest2 AuditReplayTestAuthority2::plunge_result_digest(
    const AuditNativeRequestIdentity2& request,
    std::size_t cursor,
    const AuthenticatedOperationDigest2& operation_digest,
    const NativeMotionDigest2& motion_digest,
    const AuditDepletionWitness2& depletion,
    const StockLineageDigest2& pre_lineage,
    const StockLineageDigest2& post_lineage)
{
    return AuditMotionResult2::plunge_digest(
        request,
        cursor,
        operation_digest,
        motion_digest,
        depletion,
        pre_lineage,
        post_lineage);
}

AuditResultDigest2 AuditReplayTestAuthority2::non_engaging_result_digest(
    const AuditNativeRequestIdentity2& request,
    std::size_t cursor,
    const AuthenticatedOperationDigest2& operation_digest,
    const NativeMotionDigest2& motion_digest,
    AuditNonEngagingReason2 reason,
    const StockLineageDigest2& unchanged_lineage)
{
    return AuditMotionResult2::non_engaging_digest(
        request,
        cursor,
        operation_digest,
        motion_digest,
        reason,
        unchanged_lineage);
}

AuditResultDigest2 AuditReplayTestAuthority2::completion_digest(
    const AuditNativeRequestIdentity2& request,
    std::size_t operation_count,
    const StockLineageDigest2& terminal_lineage)
{
    return AuditMotionResult2::completion_digest(
        request, operation_count, terminal_lineage);
}
