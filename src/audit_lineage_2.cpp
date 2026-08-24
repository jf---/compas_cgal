#include "audit_lineage_2.h"

#include "canonical_encoding.h"

#include <utility>

namespace {

std::string exact_cursor_bytes(std::size_t cursor)
{
    return canonical_audit_rational_bytes(Epeck::FT(cursor));
}

} // namespace

AuditLineage2::AuditLineage2(
    std::string canonical_bytes,
    StockLineageDigest2 digest)
    : canonical_bytes_(std::move(canonical_bytes)),
      digest_(std::move(digest))
{
}

AuditLineage2 AuditLineage2::seed(
    const AuditInputDigest2& input_digest,
    const AuditNativeRequestIdentity2& request)
{
    std::string canonical = canonical_encode_tagged_union(
        "stock-lineage-seed-v1",
        canonical_encode_component_map({
            {"input-digest", input_digest.bytes()},
            {"native-request-digest", request.digest().bytes()},
        }));
    return AuditLineage2(
        canonical,
        StockLineageDigestAuthority2::hash_canonical(canonical));
}

AuditLineage2 AuditLineage2::transition_lateral(
    const AuditNativeRequestIdentity2& request,
    std::size_t cursor,
    const AuditLineage2& pre_lineage,
    const AuthenticatedOperationDigest2& operation_digest,
    const NativeMotionDigest2& motion_digest,
    const AuditDecisionWitness2& decision,
    const AuditDepletionWitness2& depletion)
{
    std::string canonical = canonical_encode_tagged_union(
        "stock-lineage-lateral-transition-v1",
        canonical_encode_component_map({
            {"authenticated-operation-digest", operation_digest.bytes()},
            {"cursor", exact_cursor_bytes(cursor)},
            {"decision-witness-digest", decision.digest().bytes()},
            {"depletion-witness-digest", depletion.digest().bytes()},
            {"motion-digest", motion_digest.bytes()},
            {"native-request-digest", request.digest().bytes()},
            {"pre-lineage-digest", pre_lineage.digest().bytes()},
        }));
    return AuditLineage2(
        canonical,
        StockLineageDigestAuthority2::hash_canonical(canonical));
}

AuditLineage2 AuditLineage2::transition_plunge(
    const AuditNativeRequestIdentity2& request,
    std::size_t cursor,
    const AuditLineage2& pre_lineage,
    const AuthenticatedOperationDigest2& operation_digest,
    const NativeMotionDigest2& motion_digest,
    const AuditDepletionWitness2& depletion)
{
    std::string canonical = canonical_encode_tagged_union(
        "stock-lineage-plunge-transition-v1",
        canonical_encode_component_map({
            {"authenticated-operation-digest", operation_digest.bytes()},
            {"cursor", exact_cursor_bytes(cursor)},
            {"depletion-witness-digest", depletion.digest().bytes()},
            {"motion-digest", motion_digest.bytes()},
            {"native-request-digest", request.digest().bytes()},
            {"pre-lineage-digest", pre_lineage.digest().bytes()},
        }));
    return AuditLineage2(
        canonical,
        StockLineageDigestAuthority2::hash_canonical(canonical));
}

const StockLineageDigest2& AuditLineage2::digest() const noexcept
{
    return digest_;
}
