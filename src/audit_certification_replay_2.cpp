#include "audit_certification_internal_2.h"

#include "continuous_tea_2/circle_oracle.h"
#include "continuous_tea_2/segment_oracle.h"
#include "continuous_tea_2/segment_source.h"
#include "continuous_tea_2/sha256.h"
#include "exact_circle_chart_2.h"

#include <algorithm>
#include <optional>
#include <vector>

namespace {

bool segment_contains_point(
    const ExactSegmentMotion2& motion,
    const EPoint& point)
{
    return CGAL::collinear(motion.start, point, motion.end)
        && CGAL::collinear_are_ordered_along_line(
            motion.start,
            point,
            motion.end);
}

bool circle_contains_point(
    const ExactCircleMotion2& motion,
    const Epeck::FT& guide_radius,
    const EPoint& point)
{
    return CGAL::squared_distance(motion.center, point)
        == guide_radius * guide_radius;
}

FullCircleTeaAudit2 replay_circle_authority(
    const Stock2& stock,
    const ExactCircleMotion2& motion,
    const AuditPolicy2& policy)
{
    return audit_full_circle_tea_event_exact(
        stock,
        FullCircleEventSource2::from_exact(
            motion,
            policy.tool_radius_mm(),
            policy.engagement_cap().chord_ratio()));
}

bool same_full_circle_authority(
    const FullCircleTeaAudit2& first,
    const FullCircleTeaAudit2& second)
{
    if (first.verdict != second.verdict
        || first.trace.canonical_bytes != second.trace.canonical_bytes
        || first.violating_parameter.has_value()
            != second.violating_parameter.has_value()) {
        return false;
    }
    if (!first.violating_parameter) {
        return true;
    }
    return first.violating_parameter->chart
            == second.violating_parameter->chart
        && first.violating_parameter->parameter
            == second.violating_parameter->parameter;
}

bool same_full_circle_authority(
    const FullCircleAuthoritySnapshot& snapshot,
    const FullCircleTeaAudit2& authority)
{
    if (snapshot.verdict != authority.verdict
        || snapshot.trace_digest != authority.trace.canonical_digest
        || snapshot.violating_parameter.has_value()
            != authority.violating_parameter.has_value()) {
        return false;
    }
    if (!snapshot.violating_parameter) {
        return true;
    }
    return snapshot.violating_parameter->chart
            == authority.violating_parameter->chart
        && snapshot.violating_parameter->parameter
            == authority.violating_parameter->parameter;
}

bool same_segment_authority(
    const SegmentAuthoritySnapshot& snapshot,
    const SegmentTeaAudit2& authority)
{
    if (snapshot.verdict != authority.verdict
        || snapshot.trace_digest != authority.trace.canonical_digest
        || snapshot.violating_parameter.has_value()
            != authority.violating_parameter.has_value()) {
        return false;
    }
    return !snapshot.violating_parameter
        || snapshot.violating_parameter->parameter
            == authority.violating_parameter->parameter;
}

bool station_within_cap(
    const Stock2& stock,
    const AuditPolicy2& policy,
    const EPoint& point)
{
    return replay_audit_unguarded_station_exact(
               stock,
               point,
               policy.tool_radius_mm(),
               policy.engagement_cap().chord_ratio())
        == AuditExactStationDisposition2::WITHIN_CAP;
}

std::vector<EPoint> mandatory_points(const AuditSegmentMotion2& motion)
{
    return {
        AuditSegmentStationParameter2::build(motion, Epeck::FT(0)).point(),
        AuditSegmentStationParameter2::build(motion, Epeck::FT(1)).point(),
        AuditSegmentStationParameter2::build(
            motion, Epeck::FT(1) / Epeck::FT(2)).point(),
    };
}

std::vector<EPoint> mandatory_points(const AuditCircleMotion2& motion)
{
    std::vector<EPoint> points;
    for (int chart = 0; chart < 4; ++chart) {
        points.push_back(AuditCircleStationParameter2::build(
            motion, chart, Epeck::FT(0)).point());
    }
    return points;
}

std::vector<EPoint> mandatory_points(const AuditArcMotion2& motion)
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
    return points;
}

template <class Motion>
bool mandatory_points_are_safe(
    const Stock2& stock,
    const Motion& motion,
    const AuditPolicy2& policy)
{
    const std::vector<EPoint> points = mandatory_points(motion);
    return std::all_of(
        points.begin(), points.end(), [&](const EPoint& point) {
            return station_within_cap(stock, policy, point);
        });
}

bool counters_equal(
    const AuditDecisionCounters2& counters,
    std::size_t visited_nodes,
    std::size_t deepest_level,
    std::size_t exact_station_replays,
    std::size_t exact_coverage_replays)
{
    return counters.visited_nodes() == visited_nodes
        && counters.deepest_level() == deepest_level
        && counters.exact_station_replays() == exact_station_replays
        && counters.exact_coverage_replays() == exact_coverage_replays;
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

bool arc_owns_authority_parameter(
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

bool partial_closure_replays(
    const Stock2& stock,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy,
    const FullCircleAuthoritySnapshot& authority)
{
    if (motion.full_turn() || authority.verdict == "certified") {
        return false;
    }
    if (authority.verdict == "cap_exceeded") {
        if (!authority.violating_parameter
            || arc_owns_authority_parameter(
                motion, *authority.violating_parameter)) {
            return false;
        }
    } else if (authority.verdict != "unresolved") {
        return false;
    }
    try {
        return same_full_circle_authority(
            authority,
            replay_circle_authority(
                stock,
                ExactCircleMotion2{
                    motion.center(), motion.zero_phase(), motion.clockwise()},
                policy));
    } catch (const EventSubstrateError&) {
        return false;
    }
}

template <class Motion>
bool context_matches(
    const Stock2& stock,
    const Motion& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const std::string& stock_digest,
    const std::string& motion_digest,
    const std::string& policy_digest,
    const std::string& limits_bytes)
{
    return AuditStockStateIdentity2::build(stock).digest().bytes()
               == stock_digest
        && motion.digest().bytes() == motion_digest
        && policy.digest().bytes() == policy_digest
        && limits.canonical_bytes() == limits_bytes;
}

} // namespace

class AuditDecisionReplay2 {
public:
    template <class Motion>
    static bool context(
        const Stock2& stock,
        const Motion& motion,
        const AuditPolicy2& policy,
        const AuditDecisionLimits2& limits,
        const AuditDecisionWitness2& decision)
    {
        return context_matches(
                   stock,
                   motion,
                   policy,
                   limits,
                   decision.stock_digest_,
                   decision.motion_digest_,
                   decision.policy_digest_,
                   decision.limits_bytes_)
            && decision.digest_.bytes()
                == sha256_bytes(decision.canonical_bytes_);
    }

    static const AuditUnresolvedEvidence2::Cause* unresolved_cause(
        const AuditUnresolvedEvidence2& evidence)
    {
        return evidence.cause_.get();
    }
};

bool replay_audit_exact_station_witness(
    const Stock2& stock,
    const AuditSegmentMotion2& motion,
    const AuditPolicy2& policy,
    const AuditExactStationWitness2& witness)
{
    const auto* parameter =
        std::get_if<AuditSegmentStationParameter2>(&witness.parameter_);
    if (parameter == nullptr) {
        return false;
    }
    return witness.stock_state_digest().bytes()
               == AuditStockStateIdentity2::build(stock).digest().bytes()
        && witness.motion_digest().bytes() == motion.digest().bytes()
        && witness.policy_digest().bytes() == policy.digest().bytes()
        && parameter->motion_digest().bytes() == motion.digest().bytes()
        && segment_contains_point(motion.xy(), parameter->point())
        && replay_audit_unguarded_station_exact(
               stock,
               parameter->point(),
               policy.tool_radius_mm(),
               policy.engagement_cap().chord_ratio())
            == witness.disposition();
}

bool replay_audit_exact_station_witness(
    const Stock2& stock,
    const AuditCircleMotion2& motion,
    const AuditPolicy2& policy,
    const AuditExactStationWitness2& witness)
{
    const auto* parameter =
        std::get_if<AuditCircleStationParameter2>(&witness.parameter_);
    if (parameter == nullptr) {
        return false;
    }
    return witness.stock_state_digest().bytes()
               == AuditStockStateIdentity2::build(stock).digest().bytes()
        && witness.motion_digest().bytes() == motion.digest().bytes()
        && witness.policy_digest().bytes() == policy.digest().bytes()
        && parameter->motion_digest().bytes() == motion.digest().bytes()
        && circle_contains_point(
            motion.xy(), motion.guide_radius(), parameter->point())
        && replay_audit_unguarded_station_exact(
               stock,
               parameter->point(),
               policy.tool_radius_mm(),
               policy.engagement_cap().chord_ratio())
            == witness.disposition();
}

bool replay_audit_exact_station_witness(
    const Stock2& stock,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy,
    const AuditExactStationWitness2& witness)
{
    const auto* parameter =
        std::get_if<AuditArcStationParameter2>(&witness.parameter_);
    if (parameter == nullptr) {
        return false;
    }
    return witness.stock_state_digest().bytes()
               == AuditStockStateIdentity2::build(stock).digest().bytes()
        && witness.motion_digest().bytes() == motion.digest().bytes()
        && witness.policy_digest().bytes() == policy.digest().bytes()
        && parameter->motion_digest().bytes() == motion.digest().bytes()
        && exact_arc_point_is_incident(motion, parameter->point())
        && replay_audit_unguarded_station_exact(
               stock,
               parameter->point(),
               policy.tool_radius_mm(),
               policy.engagement_cap().chord_ratio())
            == witness.disposition();
}

bool replay_audit_certified_coverage(
    const Stock2& stock,
    const AuditSegmentMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const AuditCertifiedCoverage2& coverage)
{
    if (!context_matches(
            stock,
            motion,
            policy,
            limits,
            coverage.stock_digest_,
            coverage.motion_digest_,
            coverage.policy_digest_,
            coverage.limits_bytes_)
        || coverage.kind()
            != AuditCertifiedCoverageKind2::SEGMENT_EVENT_PARTITION) {
        return false;
    }
    try {
        const SegmentTeaAudit2 authority = audit_segment_tea_event_exact(
            stock,
            SegmentEventSource2::from_exact(
                motion.xy(),
                policy.tool_radius_mm(),
                policy.engagement_cap().chord_ratio()));
        const SegmentAuthoritySnapshot snapshot =
            validate_audit_segment_authority(authority);
        return snapshot.verdict == ContinuousTeaVerdict::CERTIFIED
            && authority.trace.canonical_bytes == coverage.authority_bytes_;
    } catch (const EventSubstrateError&) {
        return false;
    } catch (const AuditAuthorityWitnessError&) {
        return false;
    }
}

bool replay_audit_unresolved_evidence(
    const Stock2& stock,
    const AuditSegmentMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2&,
    const AuditUnresolvedEvidence2& evidence,
    const AuditDecisionCounters2& counters)
{
    const auto* cause =
        AuditDecisionReplay2::unresolved_cause(evidence);
    if (cause == nullptr
        || evidence.reason() != AuditUnresolvedReason2::EVENT_AUTHORITY_UNRESOLVED) {
        return false;
    }
    const auto* typed = std::get_if<AuditSegmentEventUnresolvedCause2>(
        &cause->payload());
    if (typed == nullptr || !mandatory_points_are_safe(stock, motion, policy)) {
        return false;
    }
    try {
        const SegmentTeaAudit2 replayed = audit_segment_tea_event_exact(
            stock,
            SegmentEventSource2::from_exact(
                motion.xy(),
                policy.tool_radius_mm(),
                policy.engagement_cap().chord_ratio()));
        return typed->authority.verdict
                == ContinuousTeaVerdict::UNRESOLVED_DEGENERACY
            && same_segment_authority(typed->authority, replayed)
            && counters_equal(counters, 0, 0, 6, 2);
    } catch (const EventSubstrateError&) {
        return false;
    }
}

bool replay_audit_unresolved_evidence(
    const Stock2& stock,
    const AuditCircleMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2&,
    const AuditUnresolvedEvidence2& evidence,
    const AuditDecisionCounters2& counters)
{
    const auto* cause =
        AuditDecisionReplay2::unresolved_cause(evidence);
    if (cause == nullptr
        || evidence.reason() != AuditUnresolvedReason2::EVENT_AUTHORITY_UNRESOLVED) {
        return false;
    }
    const auto* typed = std::get_if<AuditCircleEventUnresolvedCause2>(
        &cause->payload());
    if (typed == nullptr || !mandatory_points_are_safe(stock, motion, policy)) {
        return false;
    }
    try {
        const FullCircleTeaAudit2 replayed = replay_circle_authority(
            stock, motion.xy(), policy);
        return typed->authority.verdict == "unresolved"
            && same_full_circle_authority(typed->authority, replayed)
            && counters_equal(counters, 0, 0, 8, 2);
    } catch (const EventSubstrateError&) {
        return false;
    }
}

bool replay_audit_unresolved_evidence(
    const Stock2& stock,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const AuditUnresolvedEvidence2& evidence,
    const AuditDecisionCounters2& counters)
{
    const auto* cause =
        AuditDecisionReplay2::unresolved_cause(evidence);
    if (cause == nullptr || !mandatory_points_are_safe(stock, motion, policy)) {
        return false;
    }
    const std::size_t mandatory_count = mandatory_points(motion).size();
    if (const auto* event =
            std::get_if<AuditFullTurnArcEventUnresolvedCause2>(
                &cause->payload())) {
        if (!motion.full_turn()
            || evidence.reason()
                != AuditUnresolvedReason2::EVENT_AUTHORITY_UNRESOLVED) {
            return false;
        }
        try {
            const FullCircleTeaAudit2 replayed = replay_circle_authority(
                stock,
                ExactCircleMotion2{
                    motion.center(), motion.zero_phase(), motion.clockwise()},
                policy);
            return event->authority.verdict == "unresolved"
                && same_full_circle_authority(event->authority, replayed)
                && counters_equal(
                    counters, 0, 0, mandatory_count * 2, 2);
        } catch (const EventSubstrateError&) {
            return false;
        }
    }
    const auto* refinement = std::get_if<AuditArcRefinementCause2>(
        &cause->payload());
    if (refinement == nullptr
        || !partial_closure_replays(
            stock,
            motion,
            policy,
            refinement->partial_closure_authority)
        || refinement->terminal_cell.motion_digest().bytes()
            != motion.digest().bytes()
        || counters.exact_station_replays()
            != mandatory_count * 2 + counters.visited_nodes()
        || counters.exact_coverage_replays() != 2
        || counters.visited_nodes() == 0
        || counters.deepest_level() > limits.max_depth()) {
        return false;
    }
    if (evidence.reason() == AuditUnresolvedReason2::NODE_LIMIT_EXHAUSTED) {
        return counters.visited_nodes() == limits.max_nodes();
    }
    if (evidence.reason() == AuditUnresolvedReason2::DEPTH_LIMIT_EXHAUSTED) {
        return refinement->terminal_cell.depth() == limits.max_depth() + 1
            && counters.deepest_level() == limits.max_depth();
    }
    if (evidence.reason() == AuditUnresolvedReason2::SPATIAL_FLOOR_REACHED) {
        const AuditArcRefinementCell2& cell = refinement->terminal_cell;
        const Epeck::FT midpoint =
            (cell.start_parameter() + cell.end_parameter()) / Epeck::FT(2);
        const AuditArcStationParameter2 parameter =
            AuditArcStationParameter2::from_interval(
                motion,
                cell.interval_ordinal(),
                motion.intervals()[cell.interval_ordinal()].chart(),
                midpoint);
        return cell.depth() <= counters.deepest_level()
            && cell.squared_chord_length(motion)
                <= limits.squared_spatial_floor_mm().value()
            && station_within_cap(stock, policy, parameter.point());
    }
    return false;
}

bool replay_audit_certified_coverage(
    const Stock2& stock,
    const AuditCircleMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const AuditCertifiedCoverage2& coverage)
{
    if (!context_matches(
            stock,
            motion,
            policy,
            limits,
            coverage.stock_digest_,
            coverage.motion_digest_,
            coverage.policy_digest_,
            coverage.limits_bytes_)
        || coverage.kind()
            != AuditCertifiedCoverageKind2::FULL_CIRCLE_EVENT_PARTITION) {
        return false;
    }
    try {
        const FullCircleTeaAudit2 authority = replay_circle_authority(
            stock, motion.xy(), policy);
        const FullCircleAuthoritySnapshot snapshot =
            validate_audit_full_circle_authority(authority);
        return snapshot.verdict == "certified"
            && authority.trace.canonical_bytes == coverage.authority_bytes_;
    } catch (const EventSubstrateError&) {
        return false;
    } catch (const AuditAuthorityWitnessError&) {
        return false;
    }
}

bool replay_audit_certified_coverage(
    const Stock2& stock,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const AuditCertifiedCoverage2& coverage)
{
    if (!context_matches(
            stock,
            motion,
            policy,
            limits,
            coverage.stock_digest_,
            coverage.motion_digest_,
            coverage.policy_digest_,
            coverage.limits_bytes_)) {
        return false;
    }
    const AuditCertifiedCoverageKind2 required_kind = motion.full_turn()
        ? AuditCertifiedCoverageKind2::FULL_CIRCLE_EVENT_PARTITION
        : AuditCertifiedCoverageKind2::PARTIAL_ARC_FULL_CIRCLE_SAFE_SUPERSET;
    if (coverage.kind() != required_kind) {
        return false;
    }
    try {
        const FullCircleTeaAudit2 authority = replay_circle_authority(
            stock,
            ExactCircleMotion2{
                motion.center(),
                motion.zero_phase(),
                motion.clockwise(),
            },
            policy);
        const FullCircleAuthoritySnapshot snapshot =
            validate_audit_full_circle_authority(authority);
        return snapshot.verdict == "certified"
            && authority.trace.canonical_bytes == coverage.authority_bytes_;
    } catch (const EventSubstrateError&) {
        return false;
    } catch (const AuditAuthorityWitnessError&) {
        return false;
    }
}

template <class Motion>
bool decision_consistent(
    const Stock2& stock,
    const Motion& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const AuditDecisionWitness2& decision)
{
    if (!AuditDecisionReplay2::context(
            stock, motion, policy, limits, decision)) {
        return false;
    }
    if (decision.verdict() == AuditTeaVerdict2::CAP_EXCEEDED) {
        return decision.has_exact_station_witness()
            && !decision.has_certified_coverage()
            && !decision.has_unresolved_evidence()
            && replay_audit_exact_station_witness(
                stock,
                motion,
                policy,
                decision.exact_station_witness());
    }
    if (decision.verdict() == AuditTeaVerdict2::CERTIFIED) {
        return !decision.has_exact_station_witness()
            && decision.has_certified_coverage()
            && !decision.has_unresolved_evidence()
            && replay_audit_certified_coverage(
                stock,
                motion,
                policy,
                limits,
                decision.certified_coverage());
    }
    if (decision.has_exact_station_witness()
        || decision.has_certified_coverage()
        || !decision.has_unresolved_evidence()) {
        return false;
    }
    return replay_audit_unresolved_evidence(
               stock,
               motion,
               policy,
               limits,
               decision.unresolved_evidence(),
               decision.counters())
        && !decision.has_exact_station_witness()
        && !decision.has_certified_coverage()
        && decision.has_unresolved_evidence();
}

bool audit_decision_witness_is_self_consistent(
    const Stock2& stock,
    const AuditSegmentMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const AuditDecisionWitness2& decision)
{
    return decision_consistent(stock, motion, policy, limits, decision);
}

bool audit_decision_witness_is_self_consistent(
    const Stock2& stock,
    const AuditCircleMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const AuditDecisionWitness2& decision)
{
    return decision_consistent(stock, motion, policy, limits, decision);
}

bool audit_decision_witness_is_self_consistent(
    const Stock2& stock,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const AuditDecisionWitness2& decision)
{
    return decision_consistent(stock, motion, policy, limits, decision);
}
