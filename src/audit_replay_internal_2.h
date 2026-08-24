#pragma once

#include "audit_lineage_2.h"
#include "audit_motion_result_2.h"

#include <cstddef>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

class AuditReplay2;

struct ReplayProgress2 {
    std::size_t cursor;
    AuditLineage2 lineage;

    void swap(ReplayProgress2& other) noexcept
    {
        using std::swap;
        swap(cursor, other.cursor);
        swap(lineage, other.lineage);
    }
};

inline void swap(ReplayProgress2& left, ReplayProgress2& right) noexcept
{
    left.swap(right);
}

static_assert(std::is_nothrow_swappable_v<ReplayProgress2>);

enum class AuditReplayFailurePoint2 {
    DECISION,
    TRIAL_DEPLETION,
    WITNESS_VALIDATION,
    LINEAGE,
    RESULT,
    PROGRESS,
    FINALIZATION,
};

enum class AuditReplayRouteKind2 {
    MUTATING,
    NON_ENGAGING,
};

struct AuditReplayInstrumentation2 {
    std::size_t mutating_clone_count = 0;
    std::size_t mutating_swap_count = 0;
    std::size_t non_engaging_clone_count = 0;
    std::size_t committed_operation_count = 0;
    std::size_t failed_transaction_count = 0;
    std::size_t committed_finalization_count = 0;
    std::size_t failed_finalization_count = 0;
};

struct AuditReplayConstructionInstrumentation2 {
    std::size_t stock_clone_count = 0;
};

struct AuditReplayInspection2 {
    std::size_t cursor;
    std::string stock_digest;
    std::string lineage_digest;
    bool finalized;

    bool operator==(const AuditReplayInspection2&) const = default;
};

struct AuditReplayRequestState2 {
    AuditNativeRequestIdentity2 native_request;
    AuditPolicy2 policy;
    AuditDecisionLimits2 decision_limits;
    std::vector<NativeMotionDigest2> motion_digests;
    std::vector<AuthenticatedOperationDigest2> operation_digests;
};

struct AuditReplayLifecycle2 {
    bool finalized;
};

struct AuditReplayDiagnostics2 {
    AuditReplayInstrumentation2 instrumentation;
    std::optional<AuditReplayFailurePoint2> injected_failure;
};

struct AuditReplayStorage2 {
    AuditReplayRequestState2 request;
    Stock2 stock;
    ReplayProgress2 progress;
    AuditReplayLifecycle2 lifecycle;
    AuditReplayDiagnostics2 diagnostics;
};

class AuditReplayPrimitiveScope2 {
public:
    explicit AuditReplayPrimitiveScope2(
        AuditReplayDiagnostics2& diagnostics,
        AuditReplayRouteKind2 route_kind) noexcept;
    ~AuditReplayPrimitiveScope2();
    AuditReplayPrimitiveScope2(const AuditReplayPrimitiveScope2&) = delete;
    AuditReplayPrimitiveScope2& operator=(
        const AuditReplayPrimitiveScope2&) = delete;

private:
    AuditReplayDiagnostics2* previous_;
    AuditReplayRouteKind2 previous_route_kind_;
};

void note_audit_replay_stock_clone_for_test() noexcept;
void note_audit_replay_stock_swap_for_test() noexcept;
void require_audit_replay_route(
    const AuditReplayStorage2& storage,
    const AuthenticatedOperationDigest2& operation_digest,
    const NativeMotionDigest2& motion_digest);
void throw_if_audit_replay_failure(
    const AuditReplayStorage2& storage,
    AuditReplayFailurePoint2 point);

class AuditReplayTestAuthority2 {
public:
    static AuditLineage2 seed_lineage(
        const AuditInputDigest2& input_digest,
        const AuditNativeRequestIdentity2& request);
    static AuditLineage2 transition_lineage(
        const AuditNativeRequestIdentity2& request,
        std::size_t cursor,
        const AuditLineage2& pre_lineage,
        const AuthenticatedOperationDigest2& operation_digest,
        const NativeMotionDigest2& motion_digest,
        const AuditDecisionWitness2& decision,
        const AuditDepletionWitness2& depletion);
    static AuditLineage2 transition_plunge_lineage(
        const AuditNativeRequestIdentity2& request,
        std::size_t cursor,
        const AuditLineage2& pre_lineage,
        const AuthenticatedOperationDigest2& operation_digest,
        const NativeMotionDigest2& motion_digest,
        const AuditDepletionWitness2& depletion);
    static AuditResultDigest2 lateral_result_digest(
        const AuditNativeRequestIdentity2& request,
        std::size_t cursor,
        const AuthenticatedOperationDigest2& operation_digest,
        const NativeMotionDigest2& motion_digest,
        const AuditDecisionWitness2& decision,
        const AuditDepletionWitness2& depletion,
        const StockLineageDigest2& pre_lineage,
        const StockLineageDigest2& post_lineage);
    static AuditResultDigest2 plunge_result_digest(
        const AuditNativeRequestIdentity2& request,
        std::size_t cursor,
        const AuthenticatedOperationDigest2& operation_digest,
        const NativeMotionDigest2& motion_digest,
        const AuditDepletionWitness2& depletion,
        const StockLineageDigest2& pre_lineage,
        const StockLineageDigest2& post_lineage);
    static AuditResultDigest2 non_engaging_result_digest(
        const AuditNativeRequestIdentity2& request,
        std::size_t cursor,
        const AuthenticatedOperationDigest2& operation_digest,
        const NativeMotionDigest2& motion_digest,
        AuditNonEngagingReason2 reason,
        const StockLineageDigest2& unchanged_lineage);
    static AuditResultDigest2 completion_digest(
        const AuditNativeRequestIdentity2& request,
        std::size_t operation_count,
        const StockLineageDigest2& terminal_lineage);
};

void inject_audit_replay_failure_for_test(
    AuditReplay2& replay,
    AuditReplayFailurePoint2 point) noexcept;
void clear_audit_replay_failure_for_test(AuditReplay2& replay) noexcept;
AuditReplayInspection2 inspect_audit_replay_for_test(const AuditReplay2& replay);
AuditReplayInstrumentation2 audit_replay_instrumentation_for_test(
    const AuditReplay2& replay) noexcept;
void reset_audit_replay_construction_instrumentation_for_test() noexcept;
AuditReplayConstructionInstrumentation2
audit_replay_construction_instrumentation_for_test() noexcept;

struct AuditTrialDepletionInstrumentation2 {
    std::size_t internal_clone_count;
    std::size_t internal_swap_count;
};

void reset_audit_trial_depletion_instrumentation_for_test() noexcept;
AuditTrialDepletionInstrumentation2
audit_trial_depletion_instrumentation_for_test() noexcept;
void note_audit_trial_stock_clone_for_test() noexcept;
void note_audit_trial_stock_swap_for_test() noexcept;
void reset_legacy_arc_sweep_reachability_for_test() noexcept;
std::size_t legacy_arc_sweep_reachability_for_test() noexcept;
void note_legacy_arc_sweep_reached_for_test() noexcept;
