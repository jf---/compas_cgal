#pragma once

#include "audit_certification_2.h"
#include "continuous_tea_2/circle_oracle.h"
#include "continuous_tea_2/segment_oracle.h"

#include <memory>
#include <variant>
#include <vector>

class AuditDecisionAccounting2 {
private:
    AuditDecisionAccounting2() = default;
    void record_station_replay() noexcept;
    void record_coverage_replay() noexcept;
    void record_refinement_node(std::size_t depth) noexcept;
    std::size_t visited_nodes() const noexcept;
    std::size_t deepest_level() const noexcept;

private:
    AuditDecisionCounters2 snapshot() const;
    friend class AuditDecisionEvidenceFactory2;
    friend class AuditCertificationAdapter2;
    std::size_t visited_nodes_ = 0;
    std::size_t deepest_level_ = 0;
    std::size_t station_replays_ = 0;
    std::size_t coverage_replays_ = 0;
};

class AuditArcRefinementCell2 {
public:
    static AuditArcRefinementCell2 root(
        const AuditArcMotion2& motion,
        std::size_t interval_ordinal);
    static AuditArcRefinementCell2 first_child(
        const AuditArcMotion2& motion,
        const AuditArcRefinementCell2& parent);
    static AuditArcRefinementCell2 second_child(
        const AuditArcMotion2& motion,
        const AuditArcRefinementCell2& parent);
    const NativeMotionDigest2& motion_digest() const noexcept;
    std::size_t interval_ordinal() const noexcept;
    const Epeck::FT& start_parameter() const noexcept;
    const Epeck::FT& end_parameter() const noexcept;
    std::size_t depth() const noexcept;
    Epeck::FT squared_chord_length(const AuditArcMotion2& motion) const;
    const std::string& canonical_bytes() const noexcept;

private:
    AuditArcRefinementCell2(
        NativeMotionDigest2 motion_digest,
        std::size_t interval_ordinal,
        Epeck::FT start_parameter,
        Epeck::FT end_parameter,
        std::size_t depth,
        std::string canonical_bytes);
    NativeMotionDigest2 motion_digest_;
    std::size_t interval_ordinal_;
    Epeck::FT start_parameter_;
    Epeck::FT end_parameter_;
    std::size_t depth_;
    std::string canonical_bytes_;
};

struct SegmentAuthoritySnapshot {
    ContinuousTeaVerdict verdict;
    std::string trace_digest;
    std::optional<SegmentAuthorityParameter2> violating_parameter;
};

struct FullCircleAuthoritySnapshot {
    std::string verdict;
    std::string trace_digest;
    std::optional<FullCircleAuthorityParameter2> violating_parameter;
};

struct AuditArcRefinementCause2 {
    AuditArcRefinementCell2 terminal_cell;
    FullCircleAuthoritySnapshot partial_closure_authority;
};

struct AuditSegmentEventUnresolvedCause2 {
    SegmentAuthoritySnapshot authority;
};

struct AuditCircleEventUnresolvedCause2 {
    FullCircleAuthoritySnapshot authority;
};

struct AuditFullTurnArcEventUnresolvedCause2 {
    FullCircleAuthoritySnapshot authority;
};

SegmentAuthoritySnapshot validate_audit_segment_authority(
    const SegmentTeaAudit2& authority);
FullCircleAuthoritySnapshot validate_audit_full_circle_authority(
    const FullCircleTeaAudit2& authority);

std::string canonical_audit_decision_counters_bytes(
    const AuditDecisionCounters2& counters);

using AuditTypedUnresolvedCause2 = std::variant<
    AuditArcRefinementCause2,
    AuditSegmentEventUnresolvedCause2,
    AuditCircleEventUnresolvedCause2,
    AuditFullTurnArcEventUnresolvedCause2>;

class AuditUnresolvedEvidence2::Cause {
public:
    AuditUnresolvedReason2 reason() const noexcept;
    const AuditTypedUnresolvedCause2& payload() const noexcept;
    const std::string& canonical_bytes() const noexcept;

private:
    Cause(
        AuditUnresolvedReason2 reason,
        AuditTypedUnresolvedCause2 payload,
        std::string canonical_bytes);
    friend class AuditDecisionEvidenceFactory2;
    AuditUnresolvedReason2 reason_;
    AuditTypedUnresolvedCause2 payload_;
    std::string canonical_bytes_;
};

// Internal typed evidence seam. Production certification and native adversarial
// tests share this replaying factory; no raw bytes or forged public constructor
// can mint generated decision evidence.
class AuditDecisionEvidenceFactory2 {
public:
    static AuditDecisionWitness2 exact_station_decision(
        const Stock2& stock,
        const AuditSegmentMotion2& motion,
        const AuditPolicy2& policy,
        const AuditDecisionLimits2& limits,
        const AuditSegmentStationParameter2& parameter,
        AuditExactStationDisposition2 claimed_disposition,
        AuditDecisionAccounting2 accounting);
    static AuditDecisionWitness2 exact_station_decision(
        const Stock2& stock,
        const AuditCircleMotion2& motion,
        const AuditPolicy2& policy,
        const AuditDecisionLimits2& limits,
        const AuditCircleStationParameter2& parameter,
        AuditExactStationDisposition2 claimed_disposition,
        AuditDecisionAccounting2 accounting);
    static AuditDecisionWitness2 exact_station_decision(
        const Stock2& stock,
        const AuditArcMotion2& motion,
        const AuditPolicy2& policy,
        const AuditDecisionLimits2& limits,
        const AuditArcStationParameter2& parameter,
        AuditExactStationDisposition2 claimed_disposition);
    static AuditDecisionWitness2 exact_station_decision(
        const Stock2& stock,
        const AuditArcMotion2& motion,
        const AuditPolicy2& policy,
        const AuditDecisionLimits2& limits,
        const AuditArcStationParameter2& parameter,
        AuditExactStationDisposition2 claimed_disposition,
        AuditDecisionAccounting2 accounting);

    static AuditDecisionWitness2 certified_decision(
        const Stock2& stock,
        const AuditSegmentMotion2& motion,
        const AuditPolicy2& policy,
        const AuditDecisionLimits2& limits,
        const SegmentTeaAudit2& authority,
        AuditDecisionAccounting2 accounting);
    static AuditDecisionWitness2 certified_decision(
        const Stock2& stock,
        const AuditCircleMotion2& motion,
        const AuditPolicy2& policy,
        const AuditDecisionLimits2& limits,
        const FullCircleTeaAudit2& authority,
        AuditDecisionAccounting2 accounting);
    static AuditDecisionWitness2 certified_decision(
        const Stock2& stock,
        const AuditArcMotion2& motion,
        const AuditPolicy2& policy,
        const AuditDecisionLimits2& limits,
        const FullCircleTeaAudit2& authority,
        AuditDecisionAccounting2 accounting);

    static AuditDecisionWitness2 node_limit_exhausted(
        const Stock2& stock,
        const AuditArcMotion2& motion,
        const AuditPolicy2& policy,
        const AuditDecisionLimits2& limits,
        const FullCircleTeaAudit2& partial_closure_authority,
        const AuditArcRefinementCell2& rejected_cell,
        AuditDecisionAccounting2 accounting);
    static AuditDecisionWitness2 depth_limit_exhausted(
        const Stock2& stock,
        const AuditArcMotion2& motion,
        const AuditPolicy2& policy,
        const AuditDecisionLimits2& limits,
        const FullCircleTeaAudit2& partial_closure_authority,
        const AuditArcRefinementCell2& terminal_cell,
        AuditDecisionAccounting2 accounting);
    static AuditDecisionWitness2 spatial_floor_reached(
        const Stock2& stock,
        const AuditArcMotion2& motion,
        const AuditPolicy2& policy,
        const AuditDecisionLimits2& limits,
        const FullCircleTeaAudit2& partial_closure_authority,
        const AuditArcRefinementCell2& terminal_cell,
        AuditDecisionAccounting2 accounting);
    static AuditDecisionWitness2 event_authority_unresolved(
        const Stock2& stock,
        const AuditSegmentMotion2& motion,
        const AuditPolicy2& policy,
        const AuditDecisionLimits2& limits,
        const SegmentTeaAudit2& authority,
        AuditDecisionAccounting2 accounting);
    static AuditDecisionWitness2 event_authority_unresolved(
        const Stock2& stock,
        const AuditCircleMotion2& motion,
        const AuditPolicy2& policy,
        const AuditDecisionLimits2& limits,
        const FullCircleTeaAudit2& authority,
        AuditDecisionAccounting2 accounting);
    static AuditDecisionWitness2 event_authority_unresolved(
        const Stock2& stock,
        const AuditArcMotion2& motion,
        const AuditPolicy2& policy,
        const AuditDecisionLimits2& limits,
        const FullCircleTeaAudit2& authority,
        AuditDecisionAccounting2 accounting);

private:
    static AuditExactStationWitness2 replay_station(
        const AuditStockStateIdentity2& stock_identity,
        const NativeMotionDigest2& motion_digest,
        const AuditPolicy2& policy,
        AuditMotionStationParameter2 parameter,
        AuditExactStationDisposition2 claimed_disposition,
        const Stock2& stock);
    static AuditDecisionWitness2 make_station_decision(
        const AuditStockStateIdentity2& stock,
        const NativeMotionDigest2& motion_digest,
        const AuditPolicy2& policy,
        const AuditDecisionLimits2& limits,
        AuditExactStationWitness2 witness,
        AuditDecisionCounters2 counters);
    static AuditDecisionWitness2 make_certified_decision(
        const AuditStockStateIdentity2& stock,
        const NativeMotionDigest2& motion_digest,
        const AuditPolicy2& policy,
        const AuditDecisionLimits2& limits,
        AuditCertifiedCoverage2 coverage,
        AuditDecisionCounters2 counters);
    static AuditDecisionWitness2 make_unresolved_decision(
        const AuditStockStateIdentity2& stock,
        const NativeMotionDigest2& motion_digest,
        const AuditPolicy2& policy,
        const AuditDecisionLimits2& limits,
        std::shared_ptr<const AuditUnresolvedEvidence2::Cause> cause,
        AuditDecisionCounters2 counters);
    static void validate_counters(
        const AuditDecisionLimits2& limits,
        const AuditDecisionCounters2& counters);
};

bool replay_audit_unresolved_evidence(
    const Stock2& stock,
    const AuditSegmentMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const AuditUnresolvedEvidence2& evidence,
    const AuditDecisionCounters2& counters);
bool replay_audit_unresolved_evidence(
    const Stock2& stock,
    const AuditCircleMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const AuditUnresolvedEvidence2& evidence,
    const AuditDecisionCounters2& counters);
bool replay_audit_unresolved_evidence(
    const Stock2& stock,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const AuditUnresolvedEvidence2& evidence,
    const AuditDecisionCounters2& counters);
