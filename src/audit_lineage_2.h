#pragma once

#include "audit_certification_2.h"
#include "audit_depletion_witness_2.h"
#include "audit_request_identity_2.h"

#include <cstddef>
#include <string>

class AuditReplay2;
class AuditReplayTestAuthority2;

class AuditLineage2 {
public:
    const StockLineageDigest2& digest() const noexcept;

private:
    AuditLineage2(std::string canonical_bytes, StockLineageDigest2 digest);
    static AuditLineage2 seed(
        const AuditInputDigest2& input_digest,
        const AuditNativeRequestIdentity2& request);
    static AuditLineage2 transition_lateral(
        const AuditNativeRequestIdentity2& request,
        std::size_t cursor,
        const AuditLineage2& pre_lineage,
        const AuthenticatedOperationDigest2& operation_digest,
        const NativeMotionDigest2& motion_digest,
        const AuditDecisionWitness2& decision,
        const AuditDepletionWitness2& depletion);
    static AuditLineage2 transition_plunge(
        const AuditNativeRequestIdentity2& request,
        std::size_t cursor,
        const AuditLineage2& pre_lineage,
        const AuthenticatedOperationDigest2& operation_digest,
        const NativeMotionDigest2& motion_digest,
        const AuditDepletionWitness2& depletion);

    std::string canonical_bytes_;
    StockLineageDigest2 digest_;

    friend class AuditReplay2;
    friend class AuditReplayTestAuthority2;
};
