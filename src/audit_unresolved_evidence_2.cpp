#include "audit_certification_internal_2.h"

#include "audit_strategy_versions_2.h"
#include "canonical_encoding.h"
#include "continuous_tea_2/circle_source.h"
#include "continuous_tea_2/segment_source.h"
#include "exact_circle_chart_2.h"

#include <algorithm>
#include <utility>

namespace {
std::string unresolved_text(AuditUnresolvedReason2 reason)
{
    switch (reason) {
        case AuditUnresolvedReason2::NODE_LIMIT_EXHAUSTED:
            return "node-limit-exhausted";
        case AuditUnresolvedReason2::DEPTH_LIMIT_EXHAUSTED:
            return "depth-limit-exhausted";
        case AuditUnresolvedReason2::SPATIAL_FLOOR_REACHED:
            return "spatial-floor-reached";
        case AuditUnresolvedReason2::EVENT_AUTHORITY_UNRESOLVED:
            return "event-authority-unresolved";
    }
    throw AuditUnresolvedReasonError(
        "unresolved reason is outside the closed evidence union");
}

std::string full_circle_authority_bytes(
    const FullCircleAuthoritySnapshot& authority)
{
    const std::string parameter = authority.violating_parameter
        ? canonical_encode_tagged_union(
              "audit-full-circle-authority-parameter-v1",
              canonical_encode_component_map({
                  {"chart", canonical_audit_rational_bytes(
                       Epeck::FT(authority.violating_parameter->chart))},
                  {"parameter", canonical_audit_rational_bytes(
                       authority.violating_parameter->parameter)},
              }))
        : canonical_encode_bytes("none");
    return canonical_encode_tagged_union(
        "audit-full-circle-authority-result-v1",
        canonical_encode_component_map({
            {"parameter", parameter},
            {"trace-digest", authority.trace_digest},
            {"verdict", canonical_encode_bytes(authority.verdict)},
        }));
}

std::string full_circle_authority_bytes(
    const FullCircleTeaAudit2& authority)
{
    return full_circle_authority_bytes(
        validate_audit_full_circle_authority(authority));
}
void require_safe_station(
    const Stock2& stock,
    const AuditPolicy2& policy,
    const EPoint& point)
{
    const AuditExactStationDisposition2 disposition =
        replay_audit_unguarded_station_exact(
            stock,
            point,
            policy.tool_radius_mm(),
            policy.engagement_cap().chord_ratio());
    if (disposition == AuditExactStationDisposition2::CAP_EXCEEDED) {
        throw AuditDecisionEvidenceReplayError(
            "unresolved evidence is preempted by a mandatory live station");
    }
}

std::size_t require_mandatory_stations_safe(
    const Stock2& stock,
    const AuditSegmentMotion2& motion,
    const AuditPolicy2& policy)
{
    for (const Epeck::FT parameter : {
             Epeck::FT(0),
             Epeck::FT(1),
             Epeck::FT(1) / Epeck::FT(2),
         }) {
        require_safe_station(
            stock,
            policy,
            AuditSegmentStationParameter2::build(motion, parameter).point());
    }
    return 3;
}

std::size_t require_mandatory_stations_safe(
    const Stock2& stock,
    const AuditCircleMotion2& motion,
    const AuditPolicy2& policy)
{
    for (int chart = 0; chart < 4; ++chart) {
        require_safe_station(
            stock,
            policy,
            AuditCircleStationParameter2::build(
                motion, chart, Epeck::FT(0)).point());
    }
    return 4;
}

std::size_t require_mandatory_stations_safe(
    const Stock2& stock,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy)
{
    std::vector<EPoint> points {motion.start_point()};
    for (std::size_t ordinal = 0; ordinal < motion.intervals().size(); ++ordinal) {
        const ExactArcChartInterval2& interval = motion.intervals()[ordinal];
        const bool owns_zero =
            (interval.start_parameter() == Epeck::FT(0)
                && interval.owns_start_seam())
            || (interval.end_parameter() == Epeck::FT(0)
                && interval.owns_end_seam());
        if (!owns_zero) {
            continue;
        }
        const EPoint point = AuditArcStationParameter2::from_interval(
            motion, ordinal, interval.chart(), Epeck::FT(0)).point();
        if (std::find(points.begin(), points.end(), point) == points.end()) {
            points.push_back(point);
        }
    }
    if (std::find(points.begin(), points.end(), motion.end_point())
        == points.end()) {
        points.push_back(motion.end_point());
    }
    for (const EPoint& point : points) {
        require_safe_station(stock, policy, point);
    }
    return points.size();
}

bool interval_contains_parameter(
    const ExactArcChartInterval2& interval,
    const Epeck::FT& parameter)
{
    return interval.increasing()
        ? parameter >= interval.start_parameter()
            && parameter <= interval.end_parameter()
        : parameter <= interval.start_parameter()
            && parameter >= interval.end_parameter();
}

bool partial_arc_owns_authority_parameter(
    const AuditArcMotion2& motion,
    const FullCircleAuthorityParameter2& parameter)
{
    const EPoint point = exact_circle_chart_point(
        motion.center(),
        motion.zero_phase(),
        ExactCircleChartParameter2::build(
            parameter.chart, parameter.parameter));
    if (point == motion.start_point() || point == motion.end_point()) {
        return true;
    }
    for (std::size_t ordinal = 0; ordinal < motion.intervals().size(); ++ordinal) {
        const ExactArcChartInterval2& interval = motion.intervals()[ordinal];
        if (interval.chart() != parameter.chart
            || !interval_contains_parameter(interval, parameter.parameter)) {
            continue;
        }
        try {
            static_cast<void>(AuditArcStationParameter2::from_interval(
                motion,
                ordinal,
                parameter.chart,
                parameter.parameter));
            return true;
        } catch (const AuditStationOutsideMotionError&) {
            return false;
        }
    }
    return false;
}

void require_partial_closure_authority_replays(
    const Stock2& stock,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy,
    const FullCircleTeaAudit2& authority)
{
    if (motion.full_turn() || authority.verdict == "certified") {
        throw AuditDecisionEvidenceReplayError(
            "partial refinement requires an uncertified full-circle closure authority");
    }
    if (authority.verdict == "cap_exceeded") {
        if (!authority.violating_parameter
            || partial_arc_owns_authority_parameter(
                motion, *authority.violating_parameter)) {
            throw AuditDecisionEvidenceReplayError(
                "partial refinement closure CAP must be outside the owned arc");
        }
    } else if (authority.verdict != "unresolved") {
        throw AuditDecisionEvidenceReplayError(
            "partial refinement closure authority returned a foreign verdict");
    }
    const FullCircleTeaAudit2 replayed = audit_full_circle_tea_event_exact(
        stock,
        FullCircleEventSource2::from_exact(
            ExactCircleMotion2{
                motion.center(), motion.zero_phase(), motion.clockwise()},
            policy.tool_radius_mm(),
            policy.engagement_cap().chord_ratio()));
    if (full_circle_authority_bytes(replayed)
        != full_circle_authority_bytes(authority)) {
        throw AuditDecisionEvidenceReplayError(
            "partial refinement closure authority did not replay");
    }
}
} // namespace

AuditUnresolvedEvidence2::Cause::Cause(
    AuditUnresolvedReason2 reason,
    AuditTypedUnresolvedCause2 payload,
    std::string canonical_bytes)
    : reason_(reason),
      payload_(std::move(payload)),
      canonical_bytes_(std::move(canonical_bytes))
{
}

AuditUnresolvedReason2
AuditUnresolvedEvidence2::Cause::reason() const noexcept
{
    return reason_;
}

const AuditTypedUnresolvedCause2&
AuditUnresolvedEvidence2::Cause::payload() const noexcept
{
    return payload_;
}

const std::string&
AuditUnresolvedEvidence2::Cause::canonical_bytes() const noexcept
{
    return canonical_bytes_;
}

AuditUnresolvedEvidence2::AuditUnresolvedEvidence2(
    std::shared_ptr<const Cause> cause)
    : reason_(cause->reason()),
      cause_(std::move(cause))
{
}

AuditUnresolvedReason2 AuditUnresolvedEvidence2::reason() const noexcept
{
    return reason_;
}

AuditDecisionWitness2 AuditDecisionEvidenceFactory2::make_unresolved_decision(
    const AuditStockStateIdentity2& stock,
    const NativeMotionDigest2& motion_digest,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    std::shared_ptr<const AuditUnresolvedEvidence2::Cause> cause,
    AuditDecisionCounters2 counters)
{
    validate_counters(limits, counters);
    if (!cause) {
        throw AuditDecisionEvidenceReplayError(
            "unresolved evidence requires a typed cause");
    }
    const AuditUnresolvedReason2 reason = cause->reason();
    static_cast<void>(unresolved_text(reason));
    AuditUnresolvedEvidence2 unresolved(cause);
    const std::string evidence = canonical_encode_tagged_union(
        "audit-unresolved-evidence-v1",
        canonical_encode_component_map({
            {"cause", cause->canonical_bytes()},
            {"reason", canonical_encode_bytes(unresolved_text(reason))},
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
            {"verdict", canonical_encode_bytes("unresolved")},
        }));
    return AuditDecisionWitness2(
        AuditTeaVerdict2::UNRESOLVED,
        std::nullopt,
        std::nullopt,
        std::move(unresolved),
        std::move(counters),
        stock.digest().bytes(),
        motion_digest.bytes(),
        policy.digest().bytes(),
        limits.canonical_bytes(),
        canonical);
}

AuditDecisionWitness2 AuditDecisionEvidenceFactory2::node_limit_exhausted(
    const Stock2& stock,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const FullCircleTeaAudit2& partial_closure_authority,
    const AuditArcRefinementCell2& rejected_cell,
    AuditDecisionAccounting2 accounting)
{
    for (std::size_t count = require_mandatory_stations_safe(
             stock, motion, policy); count > 0; --count) {
        accounting.record_station_replay();
    }
    if (rejected_cell.motion_digest().bytes() != motion.digest().bytes()
        || accounting.visited_nodes() != limits.max_nodes()) {
        throw AuditDecisionCounterError(
            "node exhaustion requires a motion-owned rejected cell at the node limit");
    }
    require_partial_closure_authority_replays(
        stock, motion, policy, partial_closure_authority);
    const FullCircleAuthoritySnapshot closure_snapshot =
        validate_audit_full_circle_authority(partial_closure_authority);
    accounting.record_coverage_replay();
    const std::string cause_bytes = canonical_encode_tagged_union(
        "audit-arc-node-limit-cause-v1",
        canonical_encode_component_map({
            {"closure", full_circle_authority_bytes(closure_snapshot)},
            {"rejected-cell", rejected_cell.canonical_bytes()},
        }));
    std::shared_ptr<const AuditUnresolvedEvidence2::Cause> cause(
        new AuditUnresolvedEvidence2::Cause(
            AuditUnresolvedReason2::NODE_LIMIT_EXHAUSTED,
            AuditArcRefinementCause2{
                rejected_cell, closure_snapshot},
            cause_bytes));
    const AuditStockStateIdentity2 identity = AuditStockStateIdentity2::build(stock);
    return make_unresolved_decision(
        identity,
        motion.digest(),
        policy,
        limits,
        std::move(cause),
        accounting.snapshot());
}

AuditDecisionWitness2 AuditDecisionEvidenceFactory2::depth_limit_exhausted(
    const Stock2& stock,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const FullCircleTeaAudit2& partial_closure_authority,
    const AuditArcRefinementCell2& terminal_cell,
    AuditDecisionAccounting2 accounting)
{
    for (std::size_t count = require_mandatory_stations_safe(
             stock, motion, policy); count > 0; --count) {
        accounting.record_station_replay();
    }
    if (terminal_cell.motion_digest().bytes() != motion.digest().bytes()
        || terminal_cell.depth() != limits.max_depth() + 1
        || accounting.deepest_level() != limits.max_depth()) {
        throw AuditDecisionCounterError(
            "depth exhaustion requires a motion-owned rejected child beyond the depth limit");
    }
    require_partial_closure_authority_replays(
        stock, motion, policy, partial_closure_authority);
    const FullCircleAuthoritySnapshot closure_snapshot =
        validate_audit_full_circle_authority(partial_closure_authority);
    accounting.record_coverage_replay();
    const std::string cause_bytes = canonical_encode_tagged_union(
        "audit-arc-depth-limit-cause-v1",
        canonical_encode_component_map({
            {"closure", full_circle_authority_bytes(closure_snapshot)},
            {"rejected-child", terminal_cell.canonical_bytes()},
        }));
    std::shared_ptr<const AuditUnresolvedEvidence2::Cause> cause(
        new AuditUnresolvedEvidence2::Cause(
            AuditUnresolvedReason2::DEPTH_LIMIT_EXHAUSTED,
            AuditArcRefinementCause2{
                terminal_cell, closure_snapshot},
            cause_bytes));
    const AuditStockStateIdentity2 identity = AuditStockStateIdentity2::build(stock);
    return make_unresolved_decision(
        identity,
        motion.digest(),
        policy,
        limits,
        std::move(cause),
        accounting.snapshot());
}

AuditDecisionWitness2 AuditDecisionEvidenceFactory2::spatial_floor_reached(
    const Stock2& stock,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const FullCircleTeaAudit2& partial_closure_authority,
    const AuditArcRefinementCell2& terminal_cell,
    AuditDecisionAccounting2 accounting)
{
    for (std::size_t count = require_mandatory_stations_safe(
             stock, motion, policy); count > 0; --count) {
        accounting.record_station_replay();
    }
    if (terminal_cell.motion_digest().bytes() != motion.digest().bytes()
        || terminal_cell.squared_chord_length(motion)
            > limits.squared_spatial_floor_mm().value()) {
        throw AuditDecisionCounterError(
            "spatial-floor evidence requires a motion-owned cell at or below the floor");
    }
    require_partial_closure_authority_replays(
        stock, motion, policy, partial_closure_authority);
    const FullCircleAuthoritySnapshot closure_snapshot =
        validate_audit_full_circle_authority(partial_closure_authority);
    accounting.record_coverage_replay();
    const std::string cause_bytes = canonical_encode_tagged_union(
        "audit-arc-spatial-floor-cause-v1",
        canonical_encode_component_map({
            {"closure", full_circle_authority_bytes(closure_snapshot)},
            {"terminal-cell", terminal_cell.canonical_bytes()},
        }));
    std::shared_ptr<const AuditUnresolvedEvidence2::Cause> cause(
        new AuditUnresolvedEvidence2::Cause(
            AuditUnresolvedReason2::SPATIAL_FLOOR_REACHED,
            AuditArcRefinementCause2{
                terminal_cell, closure_snapshot},
            cause_bytes));
    const AuditStockStateIdentity2 identity = AuditStockStateIdentity2::build(stock);
    return make_unresolved_decision(
        identity,
        motion.digest(),
        policy,
        limits,
        std::move(cause),
        accounting.snapshot());
}

AuditDecisionWitness2 AuditDecisionEvidenceFactory2::event_authority_unresolved(
    const Stock2& stock,
    const AuditSegmentMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const SegmentTeaAudit2& authority,
    AuditDecisionAccounting2 accounting)
{
    const SegmentAuthoritySnapshot authority_snapshot =
        validate_audit_segment_authority(authority);
    for (std::size_t count = require_mandatory_stations_safe(
             stock, motion, policy); count > 0; --count) {
        accounting.record_station_replay();
    }
    const SegmentTeaAudit2 replayed = audit_segment_tea_event_exact(
        stock,
        SegmentEventSource2::from_exact(
            motion.xy(),
            policy.tool_radius_mm(),
            policy.engagement_cap().chord_ratio()));
    if (authority.verdict != ContinuousTeaVerdict::UNRESOLVED_DEGENERACY
        || replayed.verdict != authority.verdict
        || replayed.trace.canonical_bytes != authority.trace.canonical_bytes) {
        throw AuditDecisionEvidenceReplayError(
            "segment unresolved authority did not replay");
    }
    accounting.record_coverage_replay();
    const std::string cause_bytes = canonical_encode_tagged_union(
        "audit-segment-event-unresolved-cause-v1",
        canonical_encode_component_map({
            {"authority-trace", authority.trace.canonical_bytes},
        }));
    std::shared_ptr<const AuditUnresolvedEvidence2::Cause> cause(
        new AuditUnresolvedEvidence2::Cause(
            AuditUnresolvedReason2::EVENT_AUTHORITY_UNRESOLVED,
            AuditSegmentEventUnresolvedCause2{authority_snapshot},
            cause_bytes));
    const AuditStockStateIdentity2 identity = AuditStockStateIdentity2::build(stock);
    return make_unresolved_decision(
        identity,
        motion.digest(),
        policy,
        limits,
        std::move(cause),
        accounting.snapshot());
}

AuditDecisionWitness2 AuditDecisionEvidenceFactory2::event_authority_unresolved(
    const Stock2& stock,
    const AuditCircleMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const FullCircleTeaAudit2& authority,
    AuditDecisionAccounting2 accounting)
{
    const FullCircleAuthoritySnapshot authority_snapshot =
        validate_audit_full_circle_authority(authority);
    for (std::size_t count = require_mandatory_stations_safe(
             stock, motion, policy); count > 0; --count) {
        accounting.record_station_replay();
    }
    const FullCircleTeaAudit2 replayed = audit_full_circle_tea_event_exact(
        stock,
        FullCircleEventSource2::from_exact(
            motion.xy(),
            policy.tool_radius_mm(),
            policy.engagement_cap().chord_ratio()));
    if (authority.verdict != "unresolved"
        || replayed.verdict != authority.verdict
        || replayed.trace.canonical_bytes != authority.trace.canonical_bytes) {
        throw AuditDecisionEvidenceReplayError(
            "circle unresolved authority did not replay");
    }
    accounting.record_coverage_replay();
    const std::string cause_bytes = canonical_encode_tagged_union(
        "audit-circle-event-unresolved-cause-v1",
        canonical_encode_component_map({
            {"authority", full_circle_authority_bytes(authority_snapshot)},
        }));
    std::shared_ptr<const AuditUnresolvedEvidence2::Cause> cause(
        new AuditUnresolvedEvidence2::Cause(
            AuditUnresolvedReason2::EVENT_AUTHORITY_UNRESOLVED,
            AuditCircleEventUnresolvedCause2{authority_snapshot},
            cause_bytes));
    const AuditStockStateIdentity2 identity = AuditStockStateIdentity2::build(stock);
    return make_unresolved_decision(
        identity,
        motion.digest(),
        policy,
        limits,
        std::move(cause),
        accounting.snapshot());
}

AuditDecisionWitness2 AuditDecisionEvidenceFactory2::event_authority_unresolved(
    const Stock2& stock,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const FullCircleTeaAudit2& authority,
    AuditDecisionAccounting2 accounting)
{
    const FullCircleAuthoritySnapshot authority_snapshot =
        validate_audit_full_circle_authority(authority);
    for (std::size_t count = require_mandatory_stations_safe(
             stock, motion, policy); count > 0; --count) {
        accounting.record_station_replay();
    }
    const FullCircleTeaAudit2 replayed = audit_full_circle_tea_event_exact(
        stock,
        FullCircleEventSource2::from_exact(
            ExactCircleMotion2{
                motion.center(), motion.zero_phase(), motion.clockwise()},
            policy.tool_radius_mm(),
            policy.engagement_cap().chord_ratio()));
    if (!motion.full_turn() || authority.verdict != "unresolved"
        || replayed.verdict != authority.verdict
        || replayed.trace.canonical_bytes != authority.trace.canonical_bytes) {
        throw AuditDecisionEvidenceReplayError(
            "full-turn arc unresolved authority did not replay");
    }
    accounting.record_coverage_replay();
    const std::string cause_bytes = canonical_encode_tagged_union(
        "audit-full-turn-arc-event-unresolved-cause-v1",
        canonical_encode_component_map({
            {"authority", full_circle_authority_bytes(authority_snapshot)},
        }));
    std::shared_ptr<const AuditUnresolvedEvidence2::Cause> cause(
        new AuditUnresolvedEvidence2::Cause(
            AuditUnresolvedReason2::EVENT_AUTHORITY_UNRESOLVED,
            AuditFullTurnArcEventUnresolvedCause2{authority_snapshot},
            cause_bytes));
    const AuditStockStateIdentity2 identity = AuditStockStateIdentity2::build(stock);
    return make_unresolved_decision(
        identity,
        motion.digest(),
        policy,
        limits,
        std::move(cause),
        accounting.snapshot());
}
