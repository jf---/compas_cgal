#include "audit_certification_fixtures_2.h"

#include "audit_certification_internal_2.h"
#include "continuous_tea_2/circle_source.h"
#include "continuous_tea_2/segment_oracle.h"
#include "continuous_tea_2/segment_source.h"
#include "continuous_tea_2/sha256.h"

#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <type_traits>

namespace {

using namespace audit_certification_fixtures;

template <class Factory>
concept RawSegmentStationEvidenceMinting = requires(
    const Stock2& stock,
    const AuditSegmentMotion2& motion,
    const AuditPolicy2& audit_policy,
    const AuditDecisionLimits2& decision_limits,
    const EPoint& point,
    const std::string& bytes,
    const AuditDecisionCounters2& counters) {
    Factory::exact_station_decision(
        stock,
        motion,
        audit_policy,
        decision_limits,
        point,
        bytes,
        AuditExactStationDisposition2::CAP_EXCEEDED,
        counters);
};

static_assert(!RawSegmentStationEvidenceMinting<AuditDecisionEvidenceFactory2>);

template <class Factory>
concept GenericUnresolvedMinting = requires(
    const Stock2& stock,
    const AuditArcMotion2& motion,
    const AuditPolicy2& audit_policy,
    const AuditDecisionLimits2& decision_limits,
    const AuditDecisionCounters2& counters) {
    Factory::unresolved_decision(
        stock,
        motion,
        audit_policy,
        decision_limits,
        AuditUnresolvedReason2::NODE_LIMIT_EXHAUSTED,
        counters);
};

static_assert(!GenericUnresolvedMinting<AuditDecisionEvidenceFactory2>);
static_assert(!std::is_default_constructible_v<AuditDecisionAccounting2>);

template <class Counters>
concept RawCounterConstruction = requires {
    Counters::build(0, 0, 0, 0);
};

static_assert(!RawCounterConstruction<AuditDecisionCounters2>);

template <class Cell>
concept ArbitraryRefinementCellConstruction = requires(
    const AuditArcMotion2& motion,
    const Epeck::FT& parameter) {
    Cell::build(motion, 0, parameter, parameter, 999);
};

static_assert(!ArbitraryRefinementCellConstruction<AuditArcRefinementCell2>);
static_assert(sizeof(SegmentAuthoritySnapshot) < sizeof(SegmentTeaAudit2));
static_assert(
    sizeof(FullCircleAuthoritySnapshot) < sizeof(FullCircleTeaAudit2));
static_assert(
    sizeof(AuditTypedUnresolvedCause2)
    < sizeof(FullCircleTeaAudit2) + sizeof(AuditArcRefinementCell2));

using ArcUnresolvedReplaySignature = bool (*)(
    const Stock2&,
    const AuditArcMotion2&,
    const AuditPolicy2&,
    const AuditDecisionLimits2&,
    const AuditUnresolvedEvidence2&,
    const AuditDecisionCounters2&);

static_assert(std::is_same_v<
    decltype(static_cast<ArcUnresolvedReplaySignature>(
        &replay_audit_unresolved_evidence)),
    ArcUnresolvedReplaySignature>);

void stock_and_motion_identity_gate()
{
    // Production mutation caught: evidence binds initial stock or guide
    // geometry without binding post-depletion state and exact motion interval.
    const AuditPolicy2 audit_policy = policy();
    const AuditDecisionLimits2 decision_limits = limits();
    Stock2 stock = annular_rib_stock();
    const AuditSegmentMotion2 motion = short_segment();
    const AuditDecisionWitness2 decision = certify_audit_tea_exact(
        stock,
        motion,
        audit_policy,
        decision_limits);
    require_cap_exceeded(
        stock,
        motion,
        audit_policy,
        decision_limits,
        decision);

    Stock2 depleted = stock.clone();
    depleted.subtract_exact_segment(
        motion.xy(),
        audit_policy.tool_radius_mm(),
        audit_policy.depletion_chord_bound_mm(),
        audit_policy.center_count_limit());
    require(
        AuditStockStateIdentity2::build(depleted).digest().bytes()
                != AuditStockStateIdentity2::build(stock).digest().bytes()
            && !audit_decision_witness_is_self_consistent(
                depleted,
                motion,
                audit_policy,
                decision_limits,
                decision),
        "pre-depletion witness replayed against post-depletion stock state");

    Stock2 endpoint_stock(rectangle(-0.6, 1.4, 0.6, 2.6), {});
    const AuditArcMotion2 first = partial_arc(
        0.0,
        std::numbers::pi / 2.0,
        false);
    const AuditArcMotion2 same_guide_foreign = partial_arc(
        std::numbers::pi,
        std::numbers::pi / 2.0,
        true);
    const AuditDecisionLimits2 exhausted_limits = forced_exhaustion_limits();
    const AuditDecisionWitness2 arc_decision = certify_audit_tea_exact(
        endpoint_stock,
        first,
        audit_policy,
        exhausted_limits);
    require_cap_exceeded(
        endpoint_stock,
        first,
        audit_policy,
        exhausted_limits,
        arc_decision);
    require(
        !audit_decision_witness_is_self_consistent(
            endpoint_stock,
            same_guide_foreign,
            audit_policy,
            exhausted_limits,
            arc_decision),
        "arc witness replayed on a foreign same-guide interval");
}

AuditArcStationParameter2 owned_seam(
    const AuditArcMotion2& arc,
    int chart)
{
    if (arc.end_parameter().chart() == chart
        && arc.end_parameter().parameter() == Epeck::FT(0)) {
        return AuditArcStationParameter2::from_terminal_anchor(arc);
    }
    for (std::size_t ordinal = 0;
         ordinal < arc.intervals().size();
         ++ordinal) {
        const ExactArcChartInterval2& interval = arc.intervals()[ordinal];
        if (interval.chart() != chart) {
            continue;
        }
        const bool owns_zero =
            (interval.start_parameter() == Epeck::FT(0)
                && interval.owns_start_seam())
            || (interval.end_parameter() == Epeck::FT(0)
                && interval.owns_end_seam());
        if (owns_zero) {
            return AuditArcStationParameter2::from_interval(
                arc,
                ordinal,
                chart,
                Epeck::FT(0));
        }
    }
    throw std::runtime_error("arc fixture lacks requested owning seam");
}

void typed_decision_evidence_factory_gate()
{
    // Production mutation caught: typed evidence accepts foreign parameters,
    // false dispositions, or mutated execution context without exact replay.
    const AuditPolicy2 audit_policy = policy();
    const AuditDecisionLimits2 decision_limits = limits();
    Stock2 stock(rectangle(-0.6, 1.4, 0.6, 2.6), {});
    const AuditArcMotion2 motion = partial_arc(
        0.0,
        std::numbers::pi / 2.0,
        false);
    const AuditArcStationParameter2 safe_start = owned_seam(motion, 0);
    const AuditArcStationParameter2 live_terminal = owned_seam(motion, 1);
    require(
        replay_audit_unguarded_station_exact(
            stock,
            motion.start_point(),
            audit_policy.tool_radius_mm(),
            audit_policy.engagement_cap().chord_ratio())
                == AuditExactStationDisposition2::WITHIN_CAP
            && replay_audit_unguarded_station_exact(
                   stock,
                   motion.end_point(),
                   audit_policy.tool_radius_mm(),
                   audit_policy.engagement_cap().chord_ratio())
                == AuditExactStationDisposition2::CAP_EXCEEDED,
        "typed evidence fixture lacks safe start and live terminal");

    const AuditDecisionWitness2 decision =
        AuditDecisionEvidenceFactory2::exact_station_decision(
            stock,
            motion,
            audit_policy,
            decision_limits,
            live_terminal,
            AuditExactStationDisposition2::CAP_EXCEEDED);
    require_cap_exceeded(
        stock,
        motion,
        audit_policy,
        decision_limits,
        decision);
    require(
        decision.digest().bytes() == sha256_bytes(decision.canonical_bytes()),
        "native decision digest is not SHA-256 of canonical evidence");

    bool false_cap_rejected = false;
    try {
        static_cast<void>(
            AuditDecisionEvidenceFactory2::exact_station_decision(
                stock,
                motion,
                audit_policy,
                decision_limits,
                safe_start,
                AuditExactStationDisposition2::CAP_EXCEEDED));
    } catch (const AuditDecisionEvidenceReplayError&) {
        false_cap_rejected = true;
    }
    require(false_cap_rejected, "safe typed station forged CAP_EXCEEDED evidence");

    bool safe_clear_decision_rejected = false;
    try {
        static_cast<void>(
            AuditDecisionEvidenceFactory2::exact_station_decision(
                stock,
                motion,
                audit_policy,
                decision_limits,
                safe_start,
                AuditExactStationDisposition2::WITHIN_CAP));
    } catch (const AuditDecisionEvidenceReplayError&) {
        safe_clear_decision_rejected = true;
    }
    require(
        safe_clear_decision_rejected,
        "safe WITHIN_CAP replay minted a CAP_EXCEEDED decision");

    bool false_clear_rejected = false;
    try {
        static_cast<void>(
            AuditDecisionEvidenceFactory2::exact_station_decision(
                stock,
                motion,
                audit_policy,
                decision_limits,
                live_terminal,
                AuditExactStationDisposition2::WITHIN_CAP));
    } catch (const AuditDecisionEvidenceReplayError&) {
        false_clear_rejected = true;
    }
    require(false_clear_rejected, "live typed station forged WITHIN_CAP evidence");

    const AuditArcMotion2 foreign_same_guide = partial_arc(
        std::numbers::pi,
        std::numbers::pi / 2.0,
        true);
    const AuditArcStationParameter2 foreign_parameter =
        owned_seam(foreign_same_guide, 1);
    bool foreign_parameter_rejected = false;
    try {
        static_cast<void>(
            AuditDecisionEvidenceFactory2::exact_station_decision(
                stock,
                motion,
                audit_policy,
                decision_limits,
                foreign_parameter,
                AuditExactStationDisposition2::CAP_EXCEEDED));
    } catch (const AuditDecisionEvidenceReplayError&) {
        foreign_parameter_rejected = true;
    }
    require(
        foreign_parameter_rejected,
        "foreign typed arc parameter forged exact station evidence");

    Stock2 foreign_stock(rectangle(20.0, 20.0, 30.0, 30.0), {});
    const AuditPolicy2 foreign_policy = policy(std::numbers::pi / 3.0);
    const AuditDecisionLimits2 foreign_limits = limits(
        Epeck::FT(1) / Epeck::FT(8192),
        17,
        8193);
    require(
        !audit_decision_witness_is_self_consistent(
            foreign_stock,
            motion,
            audit_policy,
            decision_limits,
            decision)
            && !audit_decision_witness_is_self_consistent(
                stock,
                foreign_same_guide,
                audit_policy,
                decision_limits,
                decision)
            && !audit_decision_witness_is_self_consistent(
                stock,
                motion,
                foreign_policy,
                decision_limits,
                decision)
            && !audit_decision_witness_is_self_consistent(
                stock,
                motion,
                audit_policy,
                foreign_limits,
                decision),
        "typed decision evidence accepted mutated execution context");
}

void authority_snapshot_shape_gate()
{
    // Production mutations caught: certified/unresolved authorities carry a
    // forged CAP parameter, or CAP authority omits its required parameter.
    const AuditPolicy2 audit_policy = policy();
    Stock2 stock(rectangle(-10.0, -10.0, 10.0, -2.4), {});
    const AuditSegmentMotion2 segment_motion = short_segment();
    const SegmentTeaAudit2 segment_authority = audit_segment_tea_event_exact(
        stock,
        SegmentEventSource2::from_exact(
            segment_motion.xy(),
            audit_policy.tool_radius_mm(),
            audit_policy.engagement_cap().chord_ratio()));
    require(
        validate_audit_segment_authority(segment_authority).trace_digest
            == segment_authority.trace.canonical_digest,
        "segment authority snapshot lost its authenticated trace digest");

    SegmentTeaAudit2 forged_segment = segment_authority;
    forged_segment.violating_parameter =
        SegmentAuthorityParameter2{Epeck::FT(1) / Epeck::FT(2)};
    bool certified_segment_rejected = false;
    try {
        static_cast<void>(validate_audit_segment_authority(forged_segment));
    } catch (const AuditAuthorityWitnessError&) {
        certified_segment_rejected = true;
    }
    require(
        certified_segment_rejected,
        "certified segment authority accepted a forged CAP parameter");

    forged_segment.verdict = ContinuousTeaVerdict::UNRESOLVED_DEGENERACY;
    forged_segment.trace.verdict =
        ContinuousTeaVerdict::UNRESOLVED_DEGENERACY;
    bool unresolved_segment_rejected = false;
    try {
        static_cast<void>(validate_audit_segment_authority(forged_segment));
    } catch (const AuditAuthorityWitnessError&) {
        unresolved_segment_rejected = true;
    }
    require(
        unresolved_segment_rejected,
        "unresolved segment authority accepted a forged CAP parameter");
    forged_segment.verdict = ContinuousTeaVerdict::CAP_EXCEEDED;
    forged_segment.trace.verdict = ContinuousTeaVerdict::CAP_EXCEEDED;
    forged_segment.violating_parameter.reset();
    bool parameterless_segment_cap_rejected = false;
    try {
        static_cast<void>(validate_audit_segment_authority(forged_segment));
    } catch (const AuditAuthorityWitnessError&) {
        parameterless_segment_cap_rejected = true;
    }
    require(
        parameterless_segment_cap_rejected,
        "segment CAP authority omitted its required typed parameter");

    const AuditCircleMotion2 circle_motion = full_circle(false);
    const FullCircleTeaAudit2 circle_authority =
        audit_full_circle_tea_event_exact(
            stock,
            FullCircleEventSource2::from_exact(
                circle_motion.xy(),
                audit_policy.tool_radius_mm(),
                audit_policy.engagement_cap().chord_ratio()));
    require(
        validate_audit_full_circle_authority(circle_authority).trace_digest
            == circle_authority.trace.canonical_digest,
        "circle authority snapshot lost its authenticated trace digest");

    FullCircleTeaAudit2 forged_circle = circle_authority;
    forged_circle.violating_parameter = FullCircleAuthorityParameter2{
        0, Epeck::FT(1) / Epeck::FT(2)};
    bool certified_circle_rejected = false;
    try {
        static_cast<void>(validate_audit_full_circle_authority(forged_circle));
    } catch (const AuditAuthorityWitnessError&) {
        certified_circle_rejected = true;
    }
    require(
        certified_circle_rejected,
        "certified circle authority accepted a forged CAP parameter");

    forged_circle.verdict = "unresolved";
    forged_circle.trace.verdict =
        ContinuousTeaVerdict::UNRESOLVED_DEGENERACY;
    bool unresolved_circle_rejected = false;
    try {
        static_cast<void>(validate_audit_full_circle_authority(forged_circle));
    } catch (const AuditAuthorityWitnessError&) {
        unresolved_circle_rejected = true;
    }
    require(
        unresolved_circle_rejected,
        "unresolved circle authority accepted a forged CAP parameter");
    forged_circle.verdict = "cap_exceeded";
    forged_circle.trace.verdict = ContinuousTeaVerdict::CAP_EXCEEDED;
    forged_circle.violating_parameter.reset();
    bool parameterless_circle_cap_rejected = false;
    try {
        static_cast<void>(validate_audit_full_circle_authority(forged_circle));
    } catch (const AuditAuthorityWitnessError&) {
        parameterless_circle_cap_rejected = true;
    }
    require(
        parameterless_circle_cap_rejected,
        "circle CAP authority omitted its required typed parameter");
}

} // namespace

void audit_certification_identity_gate()
{
    stock_and_motion_identity_gate();
    typed_decision_evidence_factory_gate();
    authority_snapshot_shape_gate();
}
