#include "audit_depletion_witness_2.h"

#include <utility>

AuditDepletionWitness2::AuditDepletionWitness2(
    AuditDepletionKind2 kind,
    AuditStockStateDigest2 pre_stock_digest,
    AuditStockStateDigest2 post_stock_digest,
    NativeMotionDigest2 motion_digest,
    AuditPolicyDigest2 policy_digest,
    ExactDepletionTraceDigest2 construction_digest,
    std::string strategy_version,
    std::string canonical_bytes,
    DepletionWitnessDigest2 digest)
    : kind_(kind),
      pre_stock_digest_(std::move(pre_stock_digest)),
      post_stock_digest_(std::move(post_stock_digest)),
      motion_digest_(std::move(motion_digest)),
      policy_digest_(std::move(policy_digest)),
      construction_digest_(std::move(construction_digest)),
      strategy_version_(std::move(strategy_version)),
      canonical_bytes_(std::move(canonical_bytes)),
      digest_(std::move(digest))
{
}

AuditDepletionKind2 AuditDepletionWitness2::kind() const noexcept
{
    return kind_;
}

const AuditStockStateDigest2&
AuditDepletionWitness2::pre_stock_digest() const noexcept
{
    return pre_stock_digest_;
}

const AuditStockStateDigest2&
AuditDepletionWitness2::post_stock_digest() const noexcept
{
    return post_stock_digest_;
}

const DepletionWitnessDigest2&
AuditDepletionWitness2::digest() const noexcept
{
    return digest_;
}
