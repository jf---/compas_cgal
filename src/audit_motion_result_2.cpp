#include "audit_motion_result_2.h"

#include "canonical_encoding.h"

#include <limits>
#include <utility>

namespace {

std::string exact_count_bytes(std::size_t count)
{
    return canonical_audit_rational_bytes(Epeck::FT(count));
}

std::string verdict_bytes(AuditTeaVerdict2 verdict)
{
    switch (verdict) {
    case AuditTeaVerdict2::CAP_EXCEEDED:
        return canonical_encode_bytes("cap-exceeded");
    case AuditTeaVerdict2::CERTIFIED:
        return canonical_encode_bytes("certified");
    case AuditTeaVerdict2::UNRESOLVED:
        return canonical_encode_bytes("unresolved");
    }
    throw AuditResultEvidenceError("unknown audit decision verdict");
}

std::string non_engaging_reason_bytes(AuditNonEngagingReason2 reason)
{
    switch (reason) {
    case AuditNonEngagingReason2::VERTICAL_RETRACT:
        return canonical_encode_bytes("vertical-retract");
    case AuditNonEngagingReason2::CLEARANCE_TRANSPORT:
        return canonical_encode_bytes("clearance-transport");
    }
    throw AuditResultEvidenceError("unknown non-engaging reason");
}

std::size_t evidence_count(const AuditDecisionWitness2& decision)
{
    const AuditDecisionCounters2& counters = decision.counters();
    const std::size_t station = counters.exact_station_replays();
    const std::size_t coverage = counters.exact_coverage_replays();
    if (station > std::numeric_limits<std::size_t>::max() - coverage) {
        throw AuditResultEvidenceError("audit evidence count overflow");
    }
    const std::size_t count = station + coverage;
    if (count == 0) {
        throw AuditResultEvidenceError(
            "audit result requires independently replayed evidence");
    }
    return count;
}

std::string lateral_canonical_bytes(
    const AuditNativeRequestIdentity2& request,
    std::size_t cursor,
    const AuthenticatedOperationDigest2& operation_digest,
    const NativeMotionDigest2& motion_digest,
    const AuditDecisionWitness2& decision,
    const AuditDepletionWitness2& depletion,
    const StockLineageDigest2& pre_lineage,
    const StockLineageDigest2& post_lineage)
{
    return canonical_encode_tagged_union(
        "audit-lateral-result-v1",
        canonical_encode_component_map({
            {"authenticated-operation-digest", operation_digest.bytes()},
            {"cursor", exact_count_bytes(cursor)},
            {"decision-witness-digest", decision.digest().bytes()},
            {"depletion-witness-digest", depletion.digest().bytes()},
            {"evidence-count", exact_count_bytes(evidence_count(decision))},
            {"motion-digest", motion_digest.bytes()},
            {"native-request-digest", request.digest().bytes()},
            {"post-lineage-digest", post_lineage.bytes()},
            {"pre-lineage-digest", pre_lineage.bytes()},
            {"verdict", verdict_bytes(decision.verdict())},
        }));
}

std::string completion_canonical_bytes(
    const AuditNativeRequestIdentity2& request,
    std::size_t operation_count,
    const StockLineageDigest2& terminal_lineage)
{
    return canonical_encode_tagged_union(
        "audit-replay-completion-v1",
        canonical_encode_component_map({
            {"native-request-digest", request.digest().bytes()},
            {"operation-count", exact_count_bytes(operation_count)},
            {"terminal-lineage-digest", terminal_lineage.bytes()},
        }));
}

} // namespace

AuditLateralResult2::AuditLateralResult2(
    AuditTeaVerdict2 verdict,
    std::size_t evidence_count,
    AuthenticatedOperationDigest2 operation_digest,
    NativeDecisionDigest2 decision_digest,
    DepletionWitnessDigest2 depletion_digest,
    StockLineageDigest2 pre_lineage,
    StockLineageDigest2 post_lineage,
    AuditResultDigest2 digest)
    : verdict_(verdict),
      evidence_count_(evidence_count),
      operation_digest_(std::move(operation_digest)),
      decision_digest_(std::move(decision_digest)),
      depletion_digest_(std::move(depletion_digest)),
      pre_lineage_(std::move(pre_lineage)),
      post_lineage_(std::move(post_lineage)),
      digest_(std::move(digest))
{
}

AuditTeaVerdict2 AuditLateralResult2::verdict() const noexcept
{
    return verdict_;
}

std::size_t AuditLateralResult2::evidence_count() const noexcept
{
    return evidence_count_;
}

const AuthenticatedOperationDigest2&
AuditLateralResult2::authenticated_operation_digest() const noexcept
{
    return operation_digest_;
}

const NativeDecisionDigest2& AuditLateralResult2::decision_digest() const noexcept
{
    return decision_digest_;
}

const DepletionWitnessDigest2&
AuditLateralResult2::depletion_witness_digest() const noexcept
{
    return depletion_digest_;
}

const StockLineageDigest2& AuditLateralResult2::pre_lineage() const noexcept
{
    return pre_lineage_;
}

const StockLineageDigest2& AuditLateralResult2::post_lineage() const noexcept
{
    return post_lineage_;
}

const AuditResultDigest2& AuditLateralResult2::digest() const noexcept
{
    return digest_;
}

AuditPlungeResult2::AuditPlungeResult2(
    AuthenticatedOperationDigest2 operation_digest,
    DepletionWitnessDigest2 depletion_digest,
    StockLineageDigest2 pre_lineage,
    StockLineageDigest2 post_lineage,
    AuditResultDigest2 digest)
    : operation_digest_(std::move(operation_digest)),
      depletion_digest_(std::move(depletion_digest)),
      pre_lineage_(std::move(pre_lineage)),
      post_lineage_(std::move(post_lineage)),
      digest_(std::move(digest))
{
}

const AuthenticatedOperationDigest2&
AuditPlungeResult2::authenticated_operation_digest() const noexcept
{
    return operation_digest_;
}

const DepletionWitnessDigest2&
AuditPlungeResult2::depletion_witness_digest() const noexcept
{
    return depletion_digest_;
}

const StockLineageDigest2& AuditPlungeResult2::pre_lineage() const noexcept
{
    return pre_lineage_;
}

const StockLineageDigest2& AuditPlungeResult2::post_lineage() const noexcept
{
    return post_lineage_;
}

const AuditResultDigest2& AuditPlungeResult2::digest() const noexcept
{
    return digest_;
}

AuditNonEngagingResult2::AuditNonEngagingResult2(
    AuditNonEngagingReason2 reason,
    AuthenticatedOperationDigest2 operation_digest,
    StockLineageDigest2 unchanged_lineage,
    AuditResultDigest2 digest)
    : reason_(reason),
      operation_digest_(std::move(operation_digest)),
      unchanged_lineage_(std::move(unchanged_lineage)),
      digest_(std::move(digest))
{
}

AuditNonEngagingReason2 AuditNonEngagingResult2::reason() const noexcept
{
    return reason_;
}

const AuthenticatedOperationDigest2&
AuditNonEngagingResult2::authenticated_operation_digest() const noexcept
{
    return operation_digest_;
}

const StockLineageDigest2& AuditNonEngagingResult2::pre_lineage() const noexcept
{
    return unchanged_lineage_;
}

const StockLineageDigest2& AuditNonEngagingResult2::post_lineage() const noexcept
{
    return unchanged_lineage_;
}

const AuditResultDigest2& AuditNonEngagingResult2::digest() const noexcept
{
    return digest_;
}

AuditReplayCompletion2::AuditReplayCompletion2(
    AuditNativeRequestDigest2 request_digest,
    std::size_t operation_count,
    StockLineageDigest2 terminal_lineage,
    AuditResultDigest2 digest)
    : request_digest_(std::move(request_digest)),
      operation_count_(operation_count),
      terminal_lineage_(std::move(terminal_lineage)),
      digest_(std::move(digest))
{
}

const AuditNativeRequestDigest2&
AuditReplayCompletion2::request_digest() const noexcept
{
    return request_digest_;
}

std::size_t AuditReplayCompletion2::operation_count() const noexcept
{
    return operation_count_;
}

const StockLineageDigest2&
AuditReplayCompletion2::terminal_lineage() const noexcept
{
    return terminal_lineage_;
}

const AuditResultDigest2& AuditReplayCompletion2::digest() const noexcept
{
    return digest_;
}

AuditResultDigest2 AuditMotionResult2::lateral_digest(
    const AuditNativeRequestIdentity2& request,
    std::size_t cursor,
    const AuthenticatedOperationDigest2& operation_digest,
    const NativeMotionDigest2& motion_digest,
    const AuditDecisionWitness2& decision,
    const AuditDepletionWitness2& depletion,
    const StockLineageDigest2& pre_lineage,
    const StockLineageDigest2& post_lineage)
{
    return AuditResultDigestAuthority2::hash_canonical(
        lateral_canonical_bytes(
            request,
            cursor,
            operation_digest,
            motion_digest,
            decision,
            depletion,
            pre_lineage,
            post_lineage));
}

AuditLateralResult2 AuditMotionResult2::lateral(
    const AuditNativeRequestIdentity2& request,
    std::size_t cursor,
    const AuthenticatedOperationDigest2& operation_digest,
    const NativeMotionDigest2& motion_digest,
    const AuditDecisionWitness2& decision,
    const AuditDepletionWitness2& depletion,
    const AuditLineage2& pre_lineage,
    const AuditLineage2& post_lineage)
{
    const std::size_t count = evidence_count(decision);
    return AuditLateralResult2(
        decision.verdict(),
        count,
        operation_digest,
        decision.digest(),
        depletion.digest(),
        pre_lineage.digest(),
        post_lineage.digest(),
        lateral_digest(
            request,
            cursor,
            operation_digest,
            motion_digest,
            decision,
            depletion,
            pre_lineage.digest(),
            post_lineage.digest()));
}

AuditPlungeResult2 AuditMotionResult2::plunge(
    const AuditNativeRequestIdentity2& request,
    std::size_t cursor,
    const AuthenticatedOperationDigest2& operation_digest,
    const NativeMotionDigest2& motion_digest,
    const AuditDepletionWitness2& depletion,
    const AuditLineage2& pre_lineage,
    const AuditLineage2& post_lineage)
{
    return AuditPlungeResult2(
        operation_digest,
        depletion.digest(),
        pre_lineage.digest(),
        post_lineage.digest(),
        plunge_digest(
            request,
            cursor,
            operation_digest,
            motion_digest,
            depletion,
            pre_lineage.digest(),
            post_lineage.digest()));
}

AuditResultDigest2 AuditMotionResult2::plunge_digest(
    const AuditNativeRequestIdentity2& request,
    std::size_t cursor,
    const AuthenticatedOperationDigest2& operation_digest,
    const NativeMotionDigest2& motion_digest,
    const AuditDepletionWitness2& depletion,
    const StockLineageDigest2& pre_lineage,
    const StockLineageDigest2& post_lineage)
{
    const std::string canonical = canonical_encode_tagged_union(
        "audit-plunge-result-v1",
        canonical_encode_component_map({
            {"authenticated-operation-digest", operation_digest.bytes()},
            {"cursor", exact_count_bytes(cursor)},
            {"depletion-witness-digest", depletion.digest().bytes()},
            {"motion-digest", motion_digest.bytes()},
            {"native-request-digest", request.digest().bytes()},
            {"post-lineage-digest", post_lineage.bytes()},
            {"pre-lineage-digest", pre_lineage.bytes()},
        }));
    return AuditResultDigestAuthority2::hash_canonical(canonical);
}

AuditNonEngagingResult2 AuditMotionResult2::non_engaging(
    const AuditNativeRequestIdentity2& request,
    std::size_t cursor,
    const AuthenticatedOperationDigest2& operation_digest,
    const NativeMotionDigest2& motion_digest,
    AuditNonEngagingReason2 reason,
    const AuditLineage2& unchanged_lineage)
{
    return AuditNonEngagingResult2(
        reason,
        operation_digest,
        unchanged_lineage.digest(),
        non_engaging_digest(
            request,
            cursor,
            operation_digest,
            motion_digest,
            reason,
            unchanged_lineage.digest()));
}

AuditResultDigest2 AuditMotionResult2::non_engaging_digest(
    const AuditNativeRequestIdentity2& request,
    std::size_t cursor,
    const AuthenticatedOperationDigest2& operation_digest,
    const NativeMotionDigest2& motion_digest,
    AuditNonEngagingReason2 reason,
    const StockLineageDigest2& unchanged_lineage)
{
    const std::string canonical = canonical_encode_tagged_union(
        "audit-non-engaging-result-v1",
        canonical_encode_component_map({
            {"authenticated-operation-digest", operation_digest.bytes()},
            {"cursor", exact_count_bytes(cursor)},
            {"motion-digest", motion_digest.bytes()},
            {"native-request-digest", request.digest().bytes()},
            {"reason", non_engaging_reason_bytes(reason)},
            {"unchanged-lineage-digest", unchanged_lineage.bytes()},
        }));
    return AuditResultDigestAuthority2::hash_canonical(canonical);
}

AuditResultDigest2 AuditMotionResult2::completion_digest(
    const AuditNativeRequestIdentity2& request,
    std::size_t operation_count,
    const StockLineageDigest2& terminal_lineage)
{
    return AuditResultDigestAuthority2::hash_canonical(
        completion_canonical_bytes(request, operation_count, terminal_lineage));
}

AuditReplayCompletion2 AuditMotionResult2::completion(
    const AuditNativeRequestIdentity2& request,
    std::size_t operation_count,
    const AuditLineage2& terminal_lineage)
{
    return AuditReplayCompletion2(
        request.digest(),
        operation_count,
        terminal_lineage.digest(),
        completion_digest(
            request, operation_count, terminal_lineage.digest()));
}
