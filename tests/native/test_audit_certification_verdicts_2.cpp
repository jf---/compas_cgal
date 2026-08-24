#include "audit_certification_fixtures_2.h"

#include <cstddef>
#include <initializer_list>
#include <numbers>

namespace {

using namespace audit_certification_fixtures;


void nonuniform_certified_and_unresolved_gate()
{
    // Production mutation caught: incomplete coverage collapses the closed
    // CERTIFIED/CAP_EXCEEDED/UNRESOLVED verdict union to a Boolean.
    const AuditPolicy2 audit_policy = policy();
    const AuditDecisionLimits2 full_limits = limits();
    const AuditDecisionLimits2 exhausted_limits = forced_exhaustion_limits();

    Stock2 segment_stock(rectangle(-10.0, -10.0, 10.0, -0.4), {});
    const AuditSegmentMotion2 segment_motion = segment(-1.0, 0.0, 1.0, 0.0);
    require(
        segment_stock.contains(0.0, -0.49)
            && !segment_stock.contains(0.0, 0.0),
        "segment positive control is vacuously clear or fully material");
    require(
        replay_audit_unguarded_station_exact(
            segment_stock,
            EPoint(0, 0),
            audit_policy.tool_radius_mm(),
            audit_policy.engagement_cap().chord_ratio())
            == AuditExactStationDisposition2::WITHIN_CAP,
        "nonuniform segment control exceeds at its central station");
    const AuditDecisionWitness2 segment_certified = certify_audit_tea_exact(
        segment_stock,
        segment_motion,
        audit_policy,
        full_limits);
    require_certified(
        segment_stock,
        segment_motion,
        audit_policy,
        full_limits,
        segment_certified,
        AuditCertifiedCoverageKind2::SEGMENT_EVENT_PARTITION);

    const AuditDecisionWitness2 segment_exhausted = certify_audit_tea_exact(
        segment_stock,
        segment_motion,
        audit_policy,
        exhausted_limits);
    require_certified(
        segment_stock,
        segment_motion,
        audit_policy,
        exhausted_limits,
        segment_exhausted,
        AuditCertifiedCoverageKind2::SEGMENT_EVENT_PARTITION);
    require(
        segment_exhausted.counters().visited_nodes() == 0,
        "closed segment authority fabricated a refinement node");

    Stock2 circle_stock(rectangle(-10.0, -10.0, 10.0, -2.4), {});
    require(
        circle_stock.contains(0.0, -2.49)
            && !circle_stock.contains(0.0, -2.0),
        "circle positive control is vacuously clear or fully material");
    require(
        replay_audit_unguarded_station_exact(
            circle_stock,
            EPoint(0, -2),
            audit_policy.tool_radius_mm(),
            audit_policy.engagement_cap().chord_ratio())
            == AuditExactStationDisposition2::WITHIN_CAP,
        "nonuniform circle control exceeds at its material-contact station");
    for (bool clockwise : {false, true}) {
        const AuditCircleMotion2 circle_motion = full_circle(clockwise);
        const AuditDecisionWitness2 circle_certified = certify_audit_tea_exact(
            circle_stock,
            circle_motion,
            audit_policy,
            full_limits);
        require_certified(
            circle_stock,
            circle_motion,
            audit_policy,
            full_limits,
            circle_certified,
            AuditCertifiedCoverageKind2::FULL_CIRCLE_EVENT_PARTITION);

        const AuditDecisionWitness2 circle_exhausted = certify_audit_tea_exact(
            circle_stock,
            circle_motion,
            audit_policy,
            exhausted_limits);
        require_certified(
            circle_stock,
            circle_motion,
            audit_policy,
            exhausted_limits,
            circle_exhausted,
            AuditCertifiedCoverageKind2::FULL_CIRCLE_EVENT_PARTITION);
        require(
            circle_exhausted.counters().visited_nodes() == 0,
            "closed circle authority fabricated a refinement node");
    }
}


void arc_parameter_ownership_gate()
{
    // Production mutation caught: free chart parameters admit complements or
    // geometrically equal non-owning seam aliases.
    const AuditArcMotion2 quarter = partial_arc(
        0.0,
        std::numbers::pi / 2.0,
        false);
    const AuditArcStationParameter2 start =
        AuditArcStationParameter2::from_start_anchor(quarter);
    const AuditArcStationParameter2 terminal =
        AuditArcStationParameter2::from_terminal_anchor(quarter);
    require(
        start.motion_digest().bytes() == quarter.digest().bytes()
            && terminal.motion_digest().bytes() == quarter.digest().bytes()
            && start.canonical_bytes() != terminal.canonical_bytes(),
        "arc start and terminal anchors are not distinct motion-owned values");
    bool complement_rejected = false;
    try {
        static_cast<void>(AuditArcStationParameter2::from_interval(
            quarter,
            0,
            2,
            Epeck::FT(1) / Epeck::FT(2)));
    } catch (const AuditStationOutsideMotionError&) {
        complement_rejected = true;
    }
    require(complement_rejected, "complement parameter was accepted by quarter arc");

    for (bool clockwise : {false, true}) {
        const AuditArcMotion2 motion = full_turn_arc(clockwise);
        for (int owner_chart = 0; owner_chart < 4; ++owner_chart) {
            std::size_t owner_ordinal = motion.intervals().size();
            std::size_t alias_ordinal = motion.intervals().size();
            const int alias_chart = (owner_chart + 3) % 4;
            for (std::size_t ordinal = 0;
                 ordinal < motion.intervals().size();
                 ++ordinal) {
                const ExactArcChartInterval2& interval =
                    motion.intervals()[ordinal];
                if (interval.chart() == owner_chart
                    && ((interval.start_parameter() == Epeck::FT(0)
                            && interval.owns_start_seam())
                        || (interval.end_parameter() == Epeck::FT(0)
                            && interval.owns_end_seam()))) {
                    owner_ordinal = ordinal;
                }
                if (interval.chart() == alias_chart
                    && (interval.start_parameter() == Epeck::FT(1)
                        || interval.end_parameter() == Epeck::FT(1))) {
                    alias_ordinal = ordinal;
                }
            }
            require(
                owner_ordinal < motion.intervals().size()
                    && alias_ordinal < motion.intervals().size(),
                "full-turn fixture lacks owner and alias intervals");
            const AuditArcStationParameter2 owner =
                AuditArcStationParameter2::from_interval(
                    motion,
                    owner_ordinal,
                    owner_chart,
                    Epeck::FT(0));
            require(
                owner.motion_digest().bytes() == motion.digest().bytes(),
                "arc station parameter does not bind its owning motion");

            bool alias_rejected = false;
            try {
                static_cast<void>(AuditArcStationParameter2::from_interval(
                    motion,
                    alias_ordinal,
                    alias_chart,
                    Epeck::FT(1)));
            } catch (const AuditStationSeamOwnershipError&) {
                alias_rejected = true;
            }
            require(alias_rejected, "non-owning cardinal seam alias was accepted");
        }
    }
}

void nonmidpoint_authority_witness_gate()
{
    // Production mutations caught: adapters discard event-cell witnesses and
    // therefore miss live stations outside endpoints/cardinals/root midpoints.
    const AuditPolicy2 audit_policy = policy();
    const AuditDecisionLimits2 decision_limits = limits(
        Epeck::FT(1) / Epeck::FT(1 << 20), 20, 1 << 16);

    Stock2 segment_stock(rectangle(-1.6, -0.6, -0.4, 0.6), {});
    const AuditSegmentMotion2 segment_motion = segment(-2.0, 0.0, 2.0, 0.0);
    const AuditSegmentStationParameter2 segment_quarter =
        AuditSegmentStationParameter2::build(
            segment_motion, Epeck::FT(1) / Epeck::FT(4));
    require(
        replay_audit_unguarded_station_exact(
            segment_stock,
            segment_motion.xy().start,
            audit_policy.tool_radius_mm(),
            audit_policy.engagement_cap().chord_ratio())
                == AuditExactStationDisposition2::WITHIN_CAP
            && replay_audit_unguarded_station_exact(
                   segment_stock,
                   EPoint(0, 0),
                   audit_policy.tool_radius_mm(),
                   audit_policy.engagement_cap().chord_ratio())
                == AuditExactStationDisposition2::WITHIN_CAP
            && replay_audit_unguarded_station_exact(
                   segment_stock,
                   segment_motion.xy().end,
                   audit_policy.tool_radius_mm(),
                   audit_policy.engagement_cap().chord_ratio())
                == AuditExactStationDisposition2::WITHIN_CAP
            && replay_audit_unguarded_station_exact(
                   segment_stock,
                   segment_quarter.point(),
                   audit_policy.tool_radius_mm(),
                   audit_policy.engagement_cap().chord_ratio())
                == AuditExactStationDisposition2::CAP_EXCEEDED,
        "segment nonmidpoint fixture lacks safe probes and live quarter cell");
    require_cap_exceeded(
        segment_stock,
        segment_motion,
        audit_policy,
        decision_limits,
        certify_audit_tea_exact(
            segment_stock,
            segment_motion,
            audit_policy,
            decision_limits));

    constexpr double live_x = 30.0 / 17.0;
    constexpr double live_y = 16.0 / 17.0;
    Stock2 circle_stock(rectangle(
        live_x - 0.52,
        live_y - 0.52,
        live_x + 0.52,
        live_y + 0.52), {});
    const AuditCircleMotion2 circle_motion = full_circle(false);
    const AuditCircleStationParameter2 circle_quarter =
        AuditCircleStationParameter2::build(
            circle_motion, 0, Epeck::FT(1) / Epeck::FT(4));
    for (int chart = 0; chart < 4; ++chart) {
        for (const Epeck::FT& parameter : {
                 Epeck::FT(0),
                 Epeck::FT(1) / Epeck::FT(2),
             }) {
            const AuditCircleStationParameter2 old_probe =
                AuditCircleStationParameter2::build(
                    circle_motion, chart, parameter);
            require(
                replay_audit_unguarded_station_exact(
                    circle_stock,
                    old_probe.point(),
                    audit_policy.tool_radius_mm(),
                    audit_policy.engagement_cap().chord_ratio())
                    == AuditExactStationDisposition2::WITHIN_CAP,
                "circle nonmidpoint fixture is live at an old fixed probe");
        }
    }
    require_known_station_exceeds(
        circle_stock, circle_quarter.point(), audit_policy);
    require_cap_exceeded(
        circle_stock,
        circle_motion,
        audit_policy,
        decision_limits,
        certify_audit_tea_exact(
            circle_stock,
            circle_motion,
            audit_policy,
            decision_limits));

    const AuditArcMotion2 arc_motion = partial_arc(
        0.0, std::numbers::pi / 2.0, false);
    const AuditArcStationParameter2 arc_quarter =
        AuditArcStationParameter2::from_interval(
            arc_motion, 0, 0, Epeck::FT(1) / Epeck::FT(4));
    const AuditArcStationParameter2 arc_midpoint =
        AuditArcStationParameter2::from_interval(
            arc_motion, 0, 0, Epeck::FT(1) / Epeck::FT(2));
    require(
        replay_audit_unguarded_station_exact(
            circle_stock,
            arc_motion.start_point(),
            audit_policy.tool_radius_mm(),
            audit_policy.engagement_cap().chord_ratio())
                == AuditExactStationDisposition2::WITHIN_CAP
            && replay_audit_unguarded_station_exact(
                   circle_stock,
                   arc_midpoint.point(),
                   audit_policy.tool_radius_mm(),
                   audit_policy.engagement_cap().chord_ratio())
                == AuditExactStationDisposition2::WITHIN_CAP
            && replay_audit_unguarded_station_exact(
                   circle_stock,
                   arc_motion.end_point(),
                   audit_policy.tool_radius_mm(),
                   audit_policy.engagement_cap().chord_ratio())
                == AuditExactStationDisposition2::WITHIN_CAP,
        "partial-arc nonmidpoint fixture is live at an old fixed probe");
    require_known_station_exceeds(
        circle_stock, arc_quarter.point(), audit_policy);
    require_cap_exceeded(
        circle_stock,
        arc_motion,
        audit_policy,
        decision_limits,
        certify_audit_tea_exact(
            circle_stock,
            arc_motion,
            audit_policy,
            decision_limits));
}


void circle_parameter_ownership_gate()
{
    // Production mutation caught: full circles accept t=1 aliases instead of
    // the four canonical t=0 owners.
    for (bool clockwise : {false, true}) {
        const AuditCircleMotion2 motion = full_circle(clockwise);
        for (int owner_chart = 0; owner_chart < 4; ++owner_chart) {
            const AuditCircleStationParameter2 owner =
                AuditCircleStationParameter2::build(
                    motion,
                    owner_chart,
                    Epeck::FT(0));
            require(
                owner.motion_digest().bytes() == motion.digest().bytes(),
                "circle station owner does not bind its motion");

            bool alias_rejected = false;
            try {
                static_cast<void>(AuditCircleStationParameter2::build(
                    motion,
                    (owner_chart + 3) % 4,
                    Epeck::FT(1)));
            } catch (const AuditStationSeamOwnershipError&) {
                alias_rejected = true;
            }
            require(alias_rejected, "full-circle t=1 seam alias was accepted");
        }
    }
}

void full_turn_arc_gate()
{
    // Production mutation caught: Task 3 full-turn arcs lose complete circle
    // coverage or mandatory seam witness replay.
    const AuditPolicy2 audit_policy = policy();
    const AuditDecisionLimits2 full_limits = limits();
    const AuditDecisionLimits2 exhausted_limits = forced_exhaustion_limits();
    Stock2 below_cap(rectangle(-10.0, -10.0, 10.0, -2.4), {});
    Stock2 seam_stock(rectangle(1.4, -0.6, 2.6, 0.6), {});

    for (bool clockwise : {false, true}) {
        const AuditArcMotion2 motion = full_turn_arc(clockwise);
        const AuditDecisionWitness2 certified = certify_audit_tea_exact(
            below_cap,
            motion,
            audit_policy,
            full_limits);
        require_certified(
            below_cap,
            motion,
            audit_policy,
            full_limits,
            certified,
            AuditCertifiedCoverageKind2::FULL_CIRCLE_EVENT_PARTITION);

        const AuditDecisionWitness2 exceeded = certify_audit_tea_exact(
            seam_stock,
            motion,
            audit_policy,
            exhausted_limits);
        require_cap_exceeded(
            seam_stock,
            motion,
            audit_policy,
            exhausted_limits,
            exceeded);
    }
}

void partial_arc_safe_superset_gate()
{
    // Production mutation caught: partial arcs inherit complement violations
    // or replay a safe-superset proof on a foreign same-guide interval.
    const AuditPolicy2 audit_policy = policy();
    const AuditDecisionLimits2 full_limits = limits();
    const AuditArcMotion2 quarter = partial_arc(
        0.0,
        std::numbers::pi / 2.0,
        false);

    Stock2 complement_only(rectangle(-2.6, -0.6, -1.4, 0.6), {});
    const AuditDecisionWitness2 unresolved = certify_audit_tea_exact(
        complement_only,
        quarter,
        audit_policy,
        full_limits);
    require_unresolved(
        complement_only,
        quarter,
        audit_policy,
        full_limits,
        unresolved,
        AuditUnresolvedReason2::SPATIAL_FLOOR_REACHED);

    Stock2 below_cap(rectangle(-10.0, -10.0, 10.0, -2.4), {});
    const AuditDecisionWitness2 certified = certify_audit_tea_exact(
        below_cap,
        quarter,
        audit_policy,
        full_limits);
    require_certified(
        below_cap,
        quarter,
        audit_policy,
        full_limits,
        certified,
        AuditCertifiedCoverageKind2::PARTIAL_ARC_FULL_CIRCLE_SAFE_SUPERSET);

    const AuditArcMotion2 foreign_same_guide = partial_arc(
        std::numbers::pi,
        std::numbers::pi / 2.0,
        true);
    require(
        !replay_audit_certified_coverage(
            below_cap,
            foreign_same_guide,
            audit_policy,
            full_limits,
            certified.certified_coverage()),
        "partial-arc safe-superset proof replayed on a foreign same-guide arc");
}

} // namespace

void audit_certification_verdicts_gate()
{
    nonuniform_certified_and_unresolved_gate();
    arc_parameter_ownership_gate();
    circle_parameter_ownership_gate();
    full_turn_arc_gate();
    partial_arc_safe_superset_gate();
    nonmidpoint_authority_witness_gate();
}
