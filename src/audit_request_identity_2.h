#pragma once

#include "audit_policy_2.h"
#include "audit_stock_identity_2.h"

#include <stdexcept>
#include <string>
#include <vector>

class AuditNativeRequestMotionError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class AuditNativeRequestIdentity2 {
public:
    static AuditNativeRequestIdentity2 build(
        const AuditNativeStockIdentity2& stock,
        const AuditPolicy2& policy,
        std::vector<NativeMotionDigest2> motion_digests);

    const std::string& canonical_bytes() const noexcept;
    const AuditNativeRequestDigest2& digest() const noexcept;

private:
    AuditNativeRequestIdentity2(
        AuditNativeStockIdentity2 stock,
        AuditPolicyDigest2 policy_digest,
        std::vector<NativeMotionDigest2> motion_digests,
        std::string canonical_bytes,
        AuditNativeRequestDigest2 digest);

    AuditNativeStockIdentity2 stock_;
    AuditPolicyDigest2 policy_digest_;
    std::vector<NativeMotionDigest2> motion_digests_;
    std::string canonical_bytes_;
    AuditNativeRequestDigest2 digest_;
    friend class AuditReplay2;
};
