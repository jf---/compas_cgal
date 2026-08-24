#include "audit_certification_internal_2.h"

#include "audit_strategy_versions_2.h"
#include "canonical_encoding.h"
#include "continuous_tea_2/sha256.h"

#include <utility>

SegmentAuthoritySnapshot validate_audit_segment_authority(
    const SegmentTeaAudit2& authority)
{
    if (authority.trace.canonical_digest
        != sha256_bytes(authority.trace.canonical_bytes)) {
        throw AuditAuthorityWitnessError(
            "segment authority trace digest does not authenticate its canonical bytes");
    }
    const bool cap = authority.verdict == ContinuousTeaVerdict::CAP_EXCEEDED;
    if (authority.trace.verdict != authority.verdict
        || cap != authority.violating_parameter.has_value()) {
        throw AuditAuthorityWitnessError(
            "segment authority violates its closed verdict/parameter shape");
    }
    return {
        authority.verdict,
        authority.trace.canonical_digest,
        authority.violating_parameter,
    };
}

FullCircleAuthoritySnapshot validate_audit_full_circle_authority(
    const FullCircleTeaAudit2& authority)
{
    if (authority.trace.canonical_digest
        != sha256_bytes(authority.trace.canonical_bytes)) {
        throw AuditAuthorityWitnessError(
            "full-circle authority trace digest does not authenticate its canonical bytes");
    }
    ContinuousTeaVerdict trace_verdict;
    if (authority.verdict == "certified") {
        trace_verdict = ContinuousTeaVerdict::CERTIFIED;
    } else if (authority.verdict == "cap_exceeded") {
        trace_verdict = ContinuousTeaVerdict::CAP_EXCEEDED;
    } else if (authority.verdict == "unresolved") {
        trace_verdict = ContinuousTeaVerdict::UNRESOLVED_DEGENERACY;
    } else {
        throw AuditAuthorityWitnessError(
            "full-circle authority has a foreign verdict");
    }
    const bool cap = authority.verdict == "cap_exceeded";
    if (authority.trace.verdict != trace_verdict
        || cap != authority.violating_parameter.has_value()) {
        throw AuditAuthorityWitnessError(
            "full-circle authority violates its closed verdict/parameter shape");
    }
    return {
        authority.verdict,
        authority.trace.canonical_digest,
        authority.violating_parameter,
    };
}

namespace {
std::string coverage_text(AuditCertifiedCoverageKind2 kind)
{
    switch (kind) {
        case AuditCertifiedCoverageKind2::SEGMENT_EVENT_PARTITION:
            return "segment-event-partition";
        case AuditCertifiedCoverageKind2::FULL_CIRCLE_EVENT_PARTITION:
            return "full-circle-event-partition";
        case AuditCertifiedCoverageKind2::PARTIAL_ARC_FULL_CIRCLE_SAFE_SUPERSET:
            return "partial-arc-full-circle-safe-superset";
    }
    throw AuditDecisionEvidenceReplayError(
        "coverage kind is outside the closed evidence union");
}
} // namespace

AuditCertifiedCoverage2::AuditCertifiedCoverage2(
    AuditCertifiedCoverageKind2 kind,
    std::string stock_digest,
    std::string motion_digest,
    std::string policy_digest,
    std::string limits_bytes,
    std::string authority_bytes)
    : kind_(kind),
      stock_digest_(std::move(stock_digest)),
      motion_digest_(std::move(motion_digest)),
      policy_digest_(std::move(policy_digest)),
      limits_bytes_(std::move(limits_bytes)),
      authority_bytes_(std::move(authority_bytes))
{
}

AuditCertifiedCoverageKind2 AuditCertifiedCoverage2::kind() const noexcept
{
    return kind_;
}

AuditDecisionWitness2 AuditDecisionEvidenceFactory2::make_certified_decision(
    const AuditStockStateIdentity2& stock,
    const NativeMotionDigest2& motion_digest,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    AuditCertifiedCoverage2 coverage,
    AuditDecisionCounters2 counters)
{
    validate_counters(limits, counters);
    const std::string evidence = canonical_encode_tagged_union(
        "audit-certified-coverage-v1",
        canonical_encode_component_map({
            {"authority", coverage.authority_bytes_},
            {"kind", canonical_encode_bytes(coverage_text(coverage.kind_))},
        }));
    const std::string canonical = canonical_encode_tagged_union(
        "audit-native-decision-witness-v1",
        canonical_encode_component_map({
            {"counters", canonical_audit_decision_counters_bytes(counters)},
            {"evidence", evidence},
            {"limits", limits.canonical_bytes()},
            {"motion-digest", motion_digest.bytes()},
            {"policy-digest", policy.digest().bytes()},
            {"stock-state-digest", stock.digest().bytes()},
            {"strategy", canonical_encode_bytes(
                 audit_native_decision_contract_version())},
            {"verdict", canonical_encode_bytes("certified")},
        }));
    return AuditDecisionWitness2(
        AuditTeaVerdict2::CERTIFIED,
        std::nullopt,
        std::move(coverage),
        std::nullopt,
        std::move(counters),
        stock.digest().bytes(),
        motion_digest.bytes(),
        policy.digest().bytes(),
        limits.canonical_bytes(),
        canonical);
}

AuditDecisionWitness2 AuditDecisionEvidenceFactory2::certified_decision(
    const Stock2& stock,
    const AuditSegmentMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const SegmentTeaAudit2& authority,
    AuditDecisionAccounting2 accounting)
{
    static_cast<void>(validate_audit_segment_authority(authority));
    if (authority.verdict != ContinuousTeaVerdict::CERTIFIED) {
        throw AuditDecisionEvidenceReplayError(
            "segment certified evidence requires a certified authority");
    }
    const AuditStockStateIdentity2 identity =
        AuditStockStateIdentity2::build(stock);
    AuditCertifiedCoverage2 coverage(
        AuditCertifiedCoverageKind2::SEGMENT_EVENT_PARTITION,
        identity.digest().bytes(),
        motion.digest().bytes(),
        policy.digest().bytes(),
        limits.canonical_bytes(),
        authority.trace.canonical_bytes);
    if (!replay_audit_certified_coverage(
            stock, motion, policy, limits, coverage)) {
        throw AuditDecisionEvidenceReplayError(
            "segment certified coverage failed exact authority replay");
    }
    accounting.record_coverage_replay();
    return make_certified_decision(
        identity, motion.digest(), policy, limits, std::move(coverage),
        accounting.snapshot());
}

AuditDecisionWitness2 AuditDecisionEvidenceFactory2::certified_decision(
    const Stock2& stock,
    const AuditCircleMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const FullCircleTeaAudit2& authority,
    AuditDecisionAccounting2 accounting)
{
    static_cast<void>(validate_audit_full_circle_authority(authority));
    if (authority.verdict != "certified") {
        throw AuditDecisionEvidenceReplayError(
            "circle certified evidence requires a certified authority");
    }
    const AuditStockStateIdentity2 identity =
        AuditStockStateIdentity2::build(stock);
    AuditCertifiedCoverage2 coverage(
        AuditCertifiedCoverageKind2::FULL_CIRCLE_EVENT_PARTITION,
        identity.digest().bytes(),
        motion.digest().bytes(),
        policy.digest().bytes(),
        limits.canonical_bytes(),
        authority.trace.canonical_bytes);
    if (!replay_audit_certified_coverage(
            stock, motion, policy, limits, coverage)) {
        throw AuditDecisionEvidenceReplayError(
            "circle certified coverage failed exact authority replay");
    }
    accounting.record_coverage_replay();
    return make_certified_decision(
        identity, motion.digest(), policy, limits, std::move(coverage),
        accounting.snapshot());
}

AuditDecisionWitness2 AuditDecisionEvidenceFactory2::certified_decision(
    const Stock2& stock,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const FullCircleTeaAudit2& authority,
    AuditDecisionAccounting2 accounting)
{
    static_cast<void>(validate_audit_full_circle_authority(authority));
    if (authority.verdict != "certified") {
        throw AuditDecisionEvidenceReplayError(
            "arc certified evidence requires a certified authority");
    }
    const AuditStockStateIdentity2 identity =
        AuditStockStateIdentity2::build(stock);
    AuditCertifiedCoverage2 coverage(
        motion.full_turn()
            ? AuditCertifiedCoverageKind2::FULL_CIRCLE_EVENT_PARTITION
            : AuditCertifiedCoverageKind2::
                  PARTIAL_ARC_FULL_CIRCLE_SAFE_SUPERSET,
        identity.digest().bytes(),
        motion.digest().bytes(),
        policy.digest().bytes(),
        limits.canonical_bytes(),
        authority.trace.canonical_bytes);
    if (!replay_audit_certified_coverage(
            stock, motion, policy, limits, coverage)) {
        throw AuditDecisionEvidenceReplayError(
            "arc certified coverage failed exact authority replay");
    }
    accounting.record_coverage_replay();
    return make_certified_decision(
        identity, motion.digest(), policy, limits, std::move(coverage),
        accounting.snapshot());
}
