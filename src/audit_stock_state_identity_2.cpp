#include "audit_certification_2.h"

#include "canonical_encoding.h"
#include "continuous_tea_2/boundary_events.h"

#include <algorithm>
#include <utility>

AuditStockStateIdentity2 AuditStockStateIdentity2::build(
    const Stock2& stock)
{
    const std::vector<BoundaryFeatureRecord2> records =
        extract_boundary_records(stock);
    std::vector<std::string> exact_features;
    exact_features.reserve(records.size());
    for (const BoundaryFeatureRecord2& record : records) {
        exact_features.push_back(record.feature_id);
    }
    std::sort(exact_features.begin(), exact_features.end());
    std::string canonical = canonical_encode_tagged_union(
        "audit-current-exact-stock-state-v1",
        canonical_encode_component_map({
            {"boundary-feature-identities",
             canonical_encode_sequence(exact_features)},
            {"empty", canonical_encode_boolean(stock.is_empty())},
        }));
    return AuditStockStateIdentity2(
        canonical,
        AuditStockStateDigestAuthority2::hash_canonical(canonical));
}

AuditStockStateIdentity2::AuditStockStateIdentity2(
    std::string canonical_bytes,
    AuditStockStateDigest2 digest)
    : canonical_bytes_(std::move(canonical_bytes)),
      digest_(std::move(digest))
{
}

const std::string&
AuditStockStateIdentity2::canonical_bytes() const noexcept
{
    return canonical_bytes_;
}

const AuditStockStateDigest2&
AuditStockStateIdentity2::digest() const noexcept
{
    return digest_;
}
