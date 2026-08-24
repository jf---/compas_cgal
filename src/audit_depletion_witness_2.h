#pragma once

#include "audit_motion_identity_2.h"
#include "audit_policy_2.h"
#include "audit_stock_state_identity_2.h"

#include <stdexcept>
#include <string>

enum class AuditDepletionKind2 {
    SEGMENT,
    FULL_CIRCLE,
    ARC,
    PLUNGE,
};

class AuditTrialStockMismatchError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class AuditDepletionWitnessError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class AuditDepletionWitness2 {
public:
    AuditDepletionKind2 kind() const noexcept;
    const AuditStockStateDigest2& pre_stock_digest() const noexcept;
    const AuditStockStateDigest2& post_stock_digest() const noexcept;
    const DepletionWitnessDigest2& digest() const noexcept;

private:
    AuditDepletionWitness2(
        AuditDepletionKind2 kind,
        AuditStockStateDigest2 pre_stock_digest,
        AuditStockStateDigest2 post_stock_digest,
        NativeMotionDigest2 motion_digest,
        AuditPolicyDigest2 policy_digest,
        ExactDepletionTraceDigest2 construction_digest,
        std::string strategy_version,
        std::string canonical_bytes,
        DepletionWitnessDigest2 digest);

    AuditDepletionKind2 kind_;
    AuditStockStateDigest2 pre_stock_digest_;
    AuditStockStateDigest2 post_stock_digest_;
    NativeMotionDigest2 motion_digest_;
    AuditPolicyDigest2 policy_digest_;
    ExactDepletionTraceDigest2 construction_digest_;
    std::string strategy_version_;
    std::string canonical_bytes_;
    DepletionWitnessDigest2 digest_;

    friend class AuditDepletionWitnessFactory2;
};
