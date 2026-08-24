#include "audit_request_identity_2.h"

#include "canonical_encoding.h"

#include <string>
#include <utility>
#include <vector>

AuditNativeRequestIdentity2 AuditNativeRequestIdentity2::build(
    const AuditNativeStockIdentity2& stock,
    const AuditPolicy2& policy,
    std::vector<NativeMotionDigest2> motion_digests)
{
    if (motion_digests.empty()) {
        throw AuditNativeRequestMotionError(
            "native audit request requires at least one opaque motion");
    }
    std::vector<std::string> digest_bytes;
    digest_bytes.reserve(motion_digests.size());
    for (const NativeMotionDigest2& digest : motion_digests) {
        digest_bytes.push_back(digest.bytes());
    }
    std::string canonical = canonical_encode_tagged_union(
        "audit-native-request-v1",
        canonical_encode_component_map({
            {"native-motion-digests", canonical_encode_sequence(digest_bytes)},
            {"policy-digest", policy.digest().bytes()},
            {"stock-digest", stock.digest().bytes()},
        }));
    AuditNativeRequestDigest2 digest =
        AuditNativeRequestDigestAuthority2::hash_canonical(canonical);
    return AuditNativeRequestIdentity2(
        stock,
        policy.digest(),
        std::move(motion_digests),
        std::move(canonical),
        std::move(digest));
}

AuditNativeRequestIdentity2::AuditNativeRequestIdentity2(
    AuditNativeStockIdentity2 stock,
    AuditPolicyDigest2 policy_digest,
    std::vector<NativeMotionDigest2> motion_digests,
    std::string canonical_bytes,
    AuditNativeRequestDigest2 digest)
    : stock_(std::move(stock)),
      policy_digest_(std::move(policy_digest)),
      motion_digests_(std::move(motion_digests)),
      canonical_bytes_(std::move(canonical_bytes)),
      digest_(std::move(digest))
{
}

const std::string& AuditNativeRequestIdentity2::canonical_bytes() const noexcept
{
    return canonical_bytes_;
}

const AuditNativeRequestDigest2& AuditNativeRequestIdentity2::digest() const noexcept
{
    return digest_;
}
