#include "audit_certification_internal_2.h"

#include "continuous_tea_2/circle_oracle.h"
#include "continuous_tea_2/segment_oracle.h"
#include "continuous_tea_2/segment_source.h"
#include "exact_circle_chart_2.h"

#include <algorithm>
#include <optional>
#include <utility>
#include <vector>

class AuditCertificationAdapter2 {
public:
    static AuditDecisionAccounting2 begin()
    {
        return AuditDecisionAccounting2();
    }
    static void record_station(AuditDecisionAccounting2& accounting)
    {
        accounting.record_station_replay();
    }
    static void record_coverage(AuditDecisionAccounting2& accounting)
    {
        accounting.record_coverage_replay();
    }
    static void record_node(
        AuditDecisionAccounting2& accounting,
        std::size_t depth)
    {
        accounting.record_refinement_node(depth);
    }
    static std::size_t visited_nodes(
        const AuditDecisionAccounting2& accounting)
    {
        return accounting.visited_nodes();
    }
};

namespace {

AuditExactStationDisposition2 station_disposition(
    const Stock2& stock,
    const EPoint& point,
    const AuditPolicy2& policy,
    AuditDecisionAccounting2& work)
{
    AuditCertificationAdapter2::record_station(work);
    return replay_audit_unguarded_station_exact(
        stock,
        point,
        policy.tool_radius_mm(),
        policy.engagement_cap().chord_ratio());
}

FullCircleTeaAudit2 circle_authority(
    const Stock2& stock,
    const ExactCircleMotion2& motion,
    const AuditPolicy2& policy,
    AuditDecisionAccounting2& work)
{
    AuditCertificationAdapter2::record_coverage(work);
    return audit_full_circle_tea_event_exact(
        stock,
        FullCircleEventSource2::from_exact(
            motion,
            policy.tool_radius_mm(),
            policy.engagement_cap().chord_ratio()));
}

template <class Motion, class Parameter>
std::optional<AuditDecisionWitness2> replay_live_parameter(
    const Stock2& stock,
    const Motion& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const Parameter& parameter,
    AuditDecisionAccounting2& work)
{
    if (station_disposition(stock, parameter.point(), policy, work)
        != AuditExactStationDisposition2::CAP_EXCEEDED) {
        return std::nullopt;
    }
    return AuditDecisionEvidenceFactory2::exact_station_decision(
        stock,
        motion,
        policy,
        limits,
        parameter,
        AuditExactStationDisposition2::CAP_EXCEEDED,
        work);
}

void append_unique_arc_parameter(
    std::vector<AuditArcStationParameter2>& parameters,
    AuditArcStationParameter2 parameter)
{
    for (const AuditArcStationParameter2& existing : parameters) {
        if (existing.point() == parameter.point()) {
            return;
        }
    }
    parameters.push_back(std::move(parameter));
}

std::vector<AuditArcStationParameter2> arc_mandatory_parameters(
    const AuditArcMotion2& motion)
{
    std::vector<AuditArcStationParameter2> parameters;
    append_unique_arc_parameter(
        parameters,
        AuditArcStationParameter2::from_start_anchor(motion));
    for (std::size_t ordinal = 0;
         ordinal < motion.intervals().size();
         ++ordinal) {
        const ExactArcChartInterval2& interval = motion.intervals()[ordinal];
        const bool owns_zero =
            (interval.start_parameter() == Epeck::FT(0)
                && interval.owns_start_seam())
            || (interval.end_parameter() == Epeck::FT(0)
                && interval.owns_end_seam());
        if (owns_zero) {
            append_unique_arc_parameter(
                parameters,
                AuditArcStationParameter2::from_interval(
                    motion, ordinal, interval.chart(), Epeck::FT(0)));
        }
    }
    append_unique_arc_parameter(
        parameters,
        AuditArcStationParameter2::from_terminal_anchor(motion));
    return parameters;
}

bool parameter_in_interval(
    const ExactArcChartInterval2& interval,
    const Epeck::FT& parameter)
{
    if (interval.increasing()) {
        return parameter >= interval.start_parameter()
            && parameter <= interval.end_parameter();
    }
    return parameter <= interval.start_parameter()
        && parameter >= interval.end_parameter();
}

std::optional<AuditArcStationParameter2> owned_arc_parameter(
    const AuditArcMotion2& motion,
    const FullCircleAuthorityParameter2& candidate)
{
    const EPoint candidate_point = exact_circle_chart_point(
        motion.center(),
        motion.zero_phase(),
        ExactCircleChartParameter2::build(
            candidate.chart, candidate.parameter));
    if (candidate_point == motion.start_point()) {
        return AuditArcStationParameter2::from_start_anchor(motion);
    }
    if (candidate_point == motion.end_point()) {
        return AuditArcStationParameter2::from_terminal_anchor(motion);
    }
    for (std::size_t ordinal = 0;
         ordinal < motion.intervals().size();
         ++ordinal) {
        const ExactArcChartInterval2& interval = motion.intervals()[ordinal];
        if (interval.chart() != candidate.chart
            || !parameter_in_interval(interval, candidate.parameter)) {
            continue;
        }
        try {
            return AuditArcStationParameter2::from_interval(
                motion,
                ordinal,
                candidate.chart,
                candidate.parameter);
        } catch (const AuditStationOutsideMotionError&) {
            return std::nullopt;
        }
    }
    return std::nullopt;
}

struct ArcRefinementResult2 {
    std::optional<AuditArcStationParameter2> violating_parameter;
    std::optional<AuditArcRefinementCell2> terminal_cell;
    AuditUnresolvedReason2 unresolved_reason;
};

ArcRefinementResult2 refine_partial_arc(
    const Stock2& stock,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    AuditDecisionAccounting2& work)
{
    if (motion.intervals().empty()) {
        throw AuditCertificationError(
            "partial arc has no exact refinement intervals");
    }
    std::vector<AuditArcRefinementCell2> stack;
    std::optional<AuditArcRefinementCell2> first_terminal_cell;
    std::optional<AuditUnresolvedReason2> first_terminal_reason;
    for (std::size_t ordinal = motion.intervals().size(); ordinal > 0;
         --ordinal) {
        stack.push_back(AuditArcRefinementCell2::root(
            motion, ordinal - 1));
    }
    while (!stack.empty()) {
        if (AuditCertificationAdapter2::visited_nodes(work)
            == limits.max_nodes()) {
            return {
                std::nullopt,
                stack.back(),
                AuditUnresolvedReason2::NODE_LIMIT_EXHAUSTED,
            };
        }
        AuditArcRefinementCell2 node = std::move(stack.back());
        stack.pop_back();
        AuditCertificationAdapter2::record_node(work, node.depth());

        const Epeck::FT midpoint =
            (node.start_parameter() + node.end_parameter()) / Epeck::FT(2);
        const AuditArcStationParameter2 parameter =
            AuditArcStationParameter2::from_interval(
                motion,
                node.interval_ordinal(),
                motion.intervals()[node.interval_ordinal()].chart(),
                midpoint);
        if (station_disposition(stock, parameter.point(), policy, work)
            == AuditExactStationDisposition2::CAP_EXCEEDED) {
            return {
                parameter,
                std::nullopt,
                AuditUnresolvedReason2::DEPTH_LIMIT_EXHAUSTED,
            };
        }

        if (node.squared_chord_length(motion)
            <= limits.squared_spatial_floor_mm().value()) {
            if (!first_terminal_cell) {
                first_terminal_cell = node;
                first_terminal_reason =
                    AuditUnresolvedReason2::SPATIAL_FLOOR_REACHED;
            }
            continue;
        }
        if (node.depth() == limits.max_depth()) {
            if (!first_terminal_cell) {
                first_terminal_cell =
                    AuditArcRefinementCell2::first_child(motion, node);
                first_terminal_reason =
                    AuditUnresolvedReason2::DEPTH_LIMIT_EXHAUSTED;
            }
            continue;
        }

        stack.push_back(AuditArcRefinementCell2::second_child(motion, node));
        stack.push_back(AuditArcRefinementCell2::first_child(motion, node));
    }
    if (!first_terminal_cell || !first_terminal_reason) {
        throw AuditCertificationError(
            "partial arc refinement finished without a bounded terminal cell");
    }
    return {
        std::nullopt,
        first_terminal_cell,
        *first_terminal_reason,
    };
}

} // namespace

AuditDecisionWitness2 certify_audit_tea_exact(
    const Stock2& stock,
    const AuditSegmentMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits)
{
    AuditDecisionAccounting2 work = AuditCertificationAdapter2::begin();
    for (const Epeck::FT parameter_value : {
             Epeck::FT(0),
             Epeck::FT(1),
             Epeck::FT(1) / Epeck::FT(2),
         }) {
        const AuditSegmentStationParameter2 parameter =
            AuditSegmentStationParameter2::build(motion, parameter_value);
        if (const auto decision = replay_live_parameter(
                stock, motion, policy, limits, parameter, work)) {
            return *decision;
        }
    }

    try {
        AuditCertificationAdapter2::record_coverage(work);
        const SegmentTeaAudit2 authority = audit_segment_tea_event_exact(
            stock,
            SegmentEventSource2::from_exact(
                motion.xy(),
                policy.tool_radius_mm(),
                policy.engagement_cap().chord_ratio()));
        if (authority.verdict == ContinuousTeaVerdict::CAP_EXCEEDED) {
            if (!authority.violating_parameter) {
                throw AuditAuthorityWitnessError(
                    "segment CAP authority omitted its typed parameter");
            }
            const AuditSegmentStationParameter2 parameter =
                AuditSegmentStationParameter2::build(
                    motion, authority.violating_parameter->parameter);
            if (const auto decision = replay_live_parameter(
                    stock, motion, policy, limits, parameter, work)) {
                return *decision;
            }
            throw AuditAuthorityWitnessError(
                "segment CAP authority parameter did not replay live");
        }
        if (authority.verdict == ContinuousTeaVerdict::CERTIFIED) {
            return AuditDecisionEvidenceFactory2::certified_decision(
                stock, motion, policy, limits, authority, work);
        }
        return AuditDecisionEvidenceFactory2::event_authority_unresolved(
            stock, motion, policy, limits, authority, work);
    } catch (const AuditCertificationError&) {
        throw;
    }
    throw AuditDecisionEvidenceReplayError(
        "segment authority returned a foreign verdict");
}

AuditDecisionWitness2 certify_audit_tea_exact(
    const Stock2& stock,
    const AuditCircleMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits)
{
    AuditDecisionAccounting2 work = AuditCertificationAdapter2::begin();
    for (int chart = 0; chart < 4; ++chart) {
        const AuditCircleStationParameter2 parameter =
            AuditCircleStationParameter2::build(
                motion, chart, Epeck::FT(0));
        if (const auto decision = replay_live_parameter(
                stock, motion, policy, limits, parameter, work)) {
            return *decision;
        }
    }

    try {
        const FullCircleTeaAudit2 authority = circle_authority(
            stock, motion.xy(), policy, work);
        if (authority.verdict == "cap_exceeded") {
            if (!authority.violating_parameter) {
                throw AuditAuthorityWitnessError(
                    "circle CAP authority omitted its typed parameter");
            }
            const AuditCircleStationParameter2 parameter =
                AuditCircleStationParameter2::build(
                    motion,
                    authority.violating_parameter->chart,
                    authority.violating_parameter->parameter);
            if (const auto decision = replay_live_parameter(
                    stock, motion, policy, limits, parameter, work)) {
                return *decision;
            }
            throw AuditAuthorityWitnessError(
                "circle CAP authority parameter did not replay live");
        }
        if (authority.verdict == "certified") {
            return AuditDecisionEvidenceFactory2::certified_decision(
                stock, motion, policy, limits, authority, work);
        }
        return AuditDecisionEvidenceFactory2::event_authority_unresolved(
            stock, motion, policy, limits, authority, work);
    } catch (const AuditCertificationError&) {
        throw;
    }
    throw AuditDecisionEvidenceReplayError(
        "circle authority returned a foreign verdict");
}

AuditDecisionWitness2 certify_audit_tea_exact(
    const Stock2& stock,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits)
{
    AuditDecisionAccounting2 work = AuditCertificationAdapter2::begin();
    for (const AuditArcStationParameter2& parameter :
         arc_mandatory_parameters(motion)) {
        if (const auto decision = replay_live_parameter(
                stock, motion, policy, limits, parameter, work)) {
            return *decision;
        }
    }

    std::optional<FullCircleTeaAudit2> partial_closure_authority;
    try {
        const FullCircleTeaAudit2 authority = circle_authority(
            stock,
            ExactCircleMotion2{
                motion.center(), motion.zero_phase(), motion.clockwise()},
            policy,
            work);
        if (authority.verdict == "cap_exceeded"
            && authority.violating_parameter) {
            if (const auto parameter = owned_arc_parameter(
                    motion, *authority.violating_parameter)) {
                if (const auto decision = replay_live_parameter(
                        stock,
                        motion,
                        policy,
                        limits,
                        *parameter,
                        work)) {
                    return *decision;
                }
                throw AuditAuthorityWitnessError(
                    "owned arc authority parameter did not replay live");
            }
            if (motion.full_turn()) {
                throw AuditAuthorityWitnessError(
                    "full-turn CAP authority parameter is not motion-owned");
            }
        } else if (authority.verdict == "cap_exceeded") {
            throw AuditAuthorityWitnessError(
                "arc CAP authority omitted its typed parameter");
        }

        if (authority.verdict == "certified") {
            return AuditDecisionEvidenceFactory2::certified_decision(
                stock, motion, policy, limits, authority, work);
        }
        if (motion.full_turn()) {
            return AuditDecisionEvidenceFactory2::
                event_authority_unresolved(
                    stock,
                    motion,
                    policy,
                    limits,
                    authority,
                    work);
        }
        partial_closure_authority = authority;
    } catch (const AuditCertificationError&) {
        throw;
    }

    if (!partial_closure_authority) {
        throw AuditDecisionEvidenceReplayError(
            "partial arc omitted its typed full-circle closure context");
    }

    const ArcRefinementResult2 refinement = refine_partial_arc(
        stock, motion, policy, limits, work);
    if (refinement.violating_parameter) {
        return AuditDecisionEvidenceFactory2::exact_station_decision(
            stock,
            motion,
            policy,
            limits,
            *refinement.violating_parameter,
            AuditExactStationDisposition2::CAP_EXCEEDED,
            work);
    }
    if (!refinement.terminal_cell) {
        throw AuditDecisionEvidenceReplayError(
            "partial arc refinement terminated without a typed cause");
    }
    if (refinement.unresolved_reason
        == AuditUnresolvedReason2::NODE_LIMIT_EXHAUSTED) {
        return AuditDecisionEvidenceFactory2::node_limit_exhausted(
            stock,
            motion,
            policy,
            limits,
            *partial_closure_authority,
            *refinement.terminal_cell,
            work);
    }
    if (refinement.unresolved_reason
        == AuditUnresolvedReason2::DEPTH_LIMIT_EXHAUSTED) {
        return AuditDecisionEvidenceFactory2::depth_limit_exhausted(
            stock,
            motion,
            policy,
            limits,
            *partial_closure_authority,
            *refinement.terminal_cell,
            work);
    }
    if (refinement.unresolved_reason
        == AuditUnresolvedReason2::SPATIAL_FLOOR_REACHED) {
        return AuditDecisionEvidenceFactory2::spatial_floor_reached(
            stock,
            motion,
            policy,
            limits,
            *partial_closure_authority,
            *refinement.terminal_cell,
            work);
    }
    throw AuditDecisionEvidenceReplayError(
        "partial arc refinement returned a foreign unresolved reason");
}
