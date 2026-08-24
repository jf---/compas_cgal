#pragma once

#include "audit_lineage_2.h"

#include <cstddef>
#include <stdexcept>

enum class AuditNonEngagingReason2 {
    VERTICAL_RETRACT,
    CLEARANCE_TRANSPORT,
};

class AuditResultEvidenceError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class AuditLateralResult2 {
public:
    AuditLateralResult2(AuditLateralResult2&&) noexcept = default;
    AuditLateralResult2& operator=(AuditLateralResult2&&) noexcept = default;
    AuditLateralResult2(const AuditLateralResult2&) = delete;
    AuditLateralResult2& operator=(const AuditLateralResult2&) = delete;

    AuditTeaVerdict2 verdict() const noexcept;
    std::size_t evidence_count() const noexcept;
    const AuthenticatedOperationDigest2& authenticated_operation_digest()
        const noexcept;
    const NativeDecisionDigest2& decision_digest() const noexcept;
    const DepletionWitnessDigest2& depletion_witness_digest() const noexcept;
    const StockLineageDigest2& pre_lineage() const noexcept;
    const StockLineageDigest2& post_lineage() const noexcept;
    const AuditResultDigest2& digest() const noexcept;

private:
    AuditLateralResult2(
        AuditTeaVerdict2 verdict,
        std::size_t evidence_count,
        AuthenticatedOperationDigest2 operation_digest,
        NativeDecisionDigest2 decision_digest,
        DepletionWitnessDigest2 depletion_digest,
        StockLineageDigest2 pre_lineage,
        StockLineageDigest2 post_lineage,
        AuditResultDigest2 digest);

    AuditTeaVerdict2 verdict_;
    std::size_t evidence_count_;
    AuthenticatedOperationDigest2 operation_digest_;
    NativeDecisionDigest2 decision_digest_;
    DepletionWitnessDigest2 depletion_digest_;
    StockLineageDigest2 pre_lineage_;
    StockLineageDigest2 post_lineage_;
    AuditResultDigest2 digest_;
    friend class AuditMotionResult2;
};

class AuditPlungeResult2 {
public:
    AuditPlungeResult2(AuditPlungeResult2&&) noexcept = default;
    AuditPlungeResult2& operator=(AuditPlungeResult2&&) noexcept = default;
    AuditPlungeResult2(const AuditPlungeResult2&) = delete;
    AuditPlungeResult2& operator=(const AuditPlungeResult2&) = delete;

    const DepletionWitnessDigest2& depletion_witness_digest() const noexcept;
    const AuthenticatedOperationDigest2& authenticated_operation_digest()
        const noexcept;
    const StockLineageDigest2& pre_lineage() const noexcept;
    const StockLineageDigest2& post_lineage() const noexcept;
    const AuditResultDigest2& digest() const noexcept;

private:
    AuditPlungeResult2(
        AuthenticatedOperationDigest2 operation_digest,
        DepletionWitnessDigest2 depletion_digest,
        StockLineageDigest2 pre_lineage,
        StockLineageDigest2 post_lineage,
        AuditResultDigest2 digest);
    AuthenticatedOperationDigest2 operation_digest_;
    DepletionWitnessDigest2 depletion_digest_;
    StockLineageDigest2 pre_lineage_;
    StockLineageDigest2 post_lineage_;
    AuditResultDigest2 digest_;
    friend class AuditMotionResult2;
};

class AuditNonEngagingResult2 {
public:
    AuditNonEngagingResult2(AuditNonEngagingResult2&&) noexcept = default;
    AuditNonEngagingResult2& operator=(AuditNonEngagingResult2&&) noexcept = default;
    AuditNonEngagingResult2(const AuditNonEngagingResult2&) = delete;
    AuditNonEngagingResult2& operator=(const AuditNonEngagingResult2&) = delete;

    AuditNonEngagingReason2 reason() const noexcept;
    const AuthenticatedOperationDigest2& authenticated_operation_digest()
        const noexcept;
    const StockLineageDigest2& pre_lineage() const noexcept;
    const StockLineageDigest2& post_lineage() const noexcept;
    const AuditResultDigest2& digest() const noexcept;

private:
    AuditNonEngagingResult2(
        AuditNonEngagingReason2 reason,
        AuthenticatedOperationDigest2 operation_digest,
        StockLineageDigest2 unchanged_lineage,
        AuditResultDigest2 digest);
    AuditNonEngagingReason2 reason_;
    AuthenticatedOperationDigest2 operation_digest_;
    StockLineageDigest2 unchanged_lineage_;
    AuditResultDigest2 digest_;
    friend class AuditMotionResult2;
};

class AuditReplayCompletion2 {
public:
    AuditReplayCompletion2(AuditReplayCompletion2&&) noexcept = default;
    AuditReplayCompletion2& operator=(AuditReplayCompletion2&&) noexcept = default;
    AuditReplayCompletion2(const AuditReplayCompletion2&) = delete;
    AuditReplayCompletion2& operator=(const AuditReplayCompletion2&) = delete;

    const AuditNativeRequestDigest2& request_digest() const noexcept;
    std::size_t operation_count() const noexcept;
    const StockLineageDigest2& terminal_lineage() const noexcept;
    const AuditResultDigest2& digest() const noexcept;

private:
    AuditReplayCompletion2(
        AuditNativeRequestDigest2 request_digest,
        std::size_t operation_count,
        StockLineageDigest2 terminal_lineage,
        AuditResultDigest2 digest);
    AuditNativeRequestDigest2 request_digest_;
    std::size_t operation_count_;
    StockLineageDigest2 terminal_lineage_;
    AuditResultDigest2 digest_;
    friend class AuditMotionResult2;
};

class AuditMotionResult2 {
private:
    static AuditLateralResult2 lateral(
        const AuditNativeRequestIdentity2& request,
        std::size_t cursor,
        const AuthenticatedOperationDigest2& operation_digest,
        const NativeMotionDigest2& motion_digest,
        const AuditDecisionWitness2& decision,
        const AuditDepletionWitness2& depletion,
        const AuditLineage2& pre_lineage,
        const AuditLineage2& post_lineage);
    static AuditPlungeResult2 plunge(
        const AuditNativeRequestIdentity2& request,
        std::size_t cursor,
        const AuthenticatedOperationDigest2& operation_digest,
        const NativeMotionDigest2& motion_digest,
        const AuditDepletionWitness2& depletion,
        const AuditLineage2& pre_lineage,
        const AuditLineage2& post_lineage);
    static AuditNonEngagingResult2 non_engaging(
        const AuditNativeRequestIdentity2& request,
        std::size_t cursor,
        const AuthenticatedOperationDigest2& operation_digest,
        const NativeMotionDigest2& motion_digest,
        AuditNonEngagingReason2 reason,
        const AuditLineage2& unchanged_lineage);
    static AuditReplayCompletion2 completion(
        const AuditNativeRequestIdentity2& request,
        std::size_t operation_count,
        const AuditLineage2& terminal_lineage);
    static AuditResultDigest2 lateral_digest(
        const AuditNativeRequestIdentity2& request,
        std::size_t cursor,
        const AuthenticatedOperationDigest2& operation_digest,
        const NativeMotionDigest2& motion_digest,
        const AuditDecisionWitness2& decision,
        const AuditDepletionWitness2& depletion,
        const StockLineageDigest2& pre_lineage,
        const StockLineageDigest2& post_lineage);
    static AuditResultDigest2 plunge_digest(
        const AuditNativeRequestIdentity2& request,
        std::size_t cursor,
        const AuthenticatedOperationDigest2& operation_digest,
        const NativeMotionDigest2& motion_digest,
        const AuditDepletionWitness2& depletion,
        const StockLineageDigest2& pre_lineage,
        const StockLineageDigest2& post_lineage);
    static AuditResultDigest2 non_engaging_digest(
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

    friend class AuditReplay2;
    friend class AuditReplayTestAuthority2;
};
