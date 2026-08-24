#include "audit_certification_fixtures_2.h"
#include "continuous_tea_2/circle_oracle.h"
#include "continuous_tea_2/segment_source.h"

#include <array>
#include <initializer_list>
#include <numbers>

namespace {

using namespace audit_certification_fixtures;

Stock2 polygon_stock(
    std::initializer_list<std::pair<double, double>> vertices)
{
    compas::RowMatrixXd boundary(
        static_cast<Eigen::Index>(vertices.size()), 2);
    Eigen::Index row = 0;
    for (const auto& [x, y] : vertices) {
        boundary(row, 0) = x;
        boundary(row, 1) = y;
        ++row;
    }
    return Stock2(boundary, {});
}

void mandatory_witness_precedence_gate()
{
    // Production mutation caught: forced exhaustion suppresses a mandatory
    // exact seam or endpoint witness that is already live.
    const AuditPolicy2 audit_policy = policy();
    const AuditDecisionLimits2 exhausted_limits = forced_exhaustion_limits();

    Stock2 seam_stock(rectangle(1.4, -0.6, 2.6, 0.6), {});
    require_known_station_exceeds(seam_stock, EPoint(2, 0), audit_policy);
    for (bool clockwise : {false, true}) {
        const AuditCircleMotion2 motion = full_circle(clockwise);
        const AuditDecisionWitness2 decision = certify_audit_tea_exact(
            seam_stock,
            motion,
            audit_policy,
            exhausted_limits);
        require_cap_exceeded(
            seam_stock,
            motion,
            audit_policy,
            exhausted_limits,
            decision);
    }

    Stock2 endpoint_stock(rectangle(-0.6, 1.4, 0.6, 2.6), {});
    require_known_station_exceeds(endpoint_stock, EPoint(0, 2), audit_policy);
    const AuditArcMotion2 ccw = partial_arc(
        0.0,
        std::numbers::pi / 2.0,
        false);
    const AuditArcMotion2 cw = partial_arc(
        std::numbers::pi,
        std::numbers::pi / 2.0,
        true);
    for (const AuditArcMotion2& motion : {ccw, cw}) {
        const AuditDecisionWitness2 decision = certify_audit_tea_exact(
            endpoint_stock,
            motion,
            audit_policy,
            exhausted_limits);
        require_cap_exceeded(
            endpoint_stock,
            motion,
            audit_policy,
            exhausted_limits,
            decision);
    }
}

void partial_arc_budget_reason_gate()
{
    // Production mutations caught: partial-arc sampling ignores sealed depth,
    // node, or spatial-floor budgets, or reports counters unrelated to work.
    const AuditPolicy2 audit_policy = policy();
    const AuditArcMotion2 motion = partial_arc(
        0.0, std::numbers::pi / 2.0, false);
    Stock2 complement_only(rectangle(-2.6, -0.6, -1.4, 0.6), {});

    const struct {
        AuditDecisionLimits2 limits;
        AuditUnresolvedReason2 reason;
    } cases[] {
        {limits(Epeck::FT(1) / Epeck::FT(1 << 20), 0, 64),
         AuditUnresolvedReason2::DEPTH_LIMIT_EXHAUSTED},
        {limits(Epeck::FT(1) / Epeck::FT(1 << 20), 16, 1),
         AuditUnresolvedReason2::NODE_LIMIT_EXHAUSTED},
        {limits(Epeck::FT(100), 16, 64),
         AuditUnresolvedReason2::SPATIAL_FLOOR_REACHED},
    };
    for (const auto& item : cases) {
        const AuditDecisionWitness2 decision = certify_audit_tea_exact(
            complement_only,
            motion,
            audit_policy,
            item.limits);
        require_unresolved(
            complement_only,
            motion,
            audit_policy,
            item.limits,
            decision,
            item.reason);
        require(
            decision.counters().exact_station_replays() >= 2,
            "partial arc did not replay mandatory start and terminal anchors");
    }
}

void later_sibling_witness_precedence_gate()
{
    // Production mutation caught: refinement returns UNRESOLVED at its first
    // safe stopping leaf instead of scanning queued siblings for a live exact
    // witness. Each full-circle authority witness is outside the partial arc,
    // the anchors and earlier DFS probes are safe, and only the named later
    // sibling violates the cap.
    const AuditPolicy2 audit_policy = policy();
    const AuditDecisionLimits2 decision_limits = limits(
        Epeck::FT(1) / Epeck::FT(4096), 1, 64);

    struct SiblingFixture {
        Stock2 stock;
        AuditArcMotion2 motion;
        int chart;
        Epeck::FT live_parameter;
        std::array<Epeck::FT, 2> safe_parameters;
    };
    const SiblingFixture fixtures[] {
        {
            polygon_stock({
                {-3.0 / 2.0, 5.0 / 4.0},
                {-3.0 / 2.0, 3.0},
                {7.0 / 4.0, 3.0},
                {7.0 / 4.0, 17.0 / 16.0},
                {11.0 / 16.0, 17.0 / 16.0},
                {11.0 / 16.0, 17.0 / 8.0},
                {27.0 / 16.0, 17.0 / 8.0},
                {27.0 / 16.0, 47.0 / 16.0},
                {-23.0 / 16.0, 47.0 / 16.0},
                {-23.0 / 16.0, 37.0 / 16.0},
                {-3.0 / 8.0, 37.0 / 16.0},
                {-3.0 / 8.0, 5.0 / 4.0},
            }),
            partial_arc(
                std::numbers::pi,
                std::numbers::pi / 2.0,
                true),
            1,
            Epeck::FT(1) / Epeck::FT(4),
            {
                Epeck::FT(1) / Epeck::FT(2),
                Epeck::FT(3) / Epeck::FT(4),
            },
        },
        {
            polygon_stock({
                {-39.0 / 32.0, -3.0},
                {-39.0 / 32.0, -61.0 / 32.0},
                {-35.0 / 64.0, -61.0 / 32.0},
                {-35.0 / 64.0, -41.0 / 16.0},
                {-37.0 / 32.0, -41.0 / 16.0},
                {-37.0 / 32.0, -47.0 / 16.0},
                {47.0 / 16.0, -47.0 / 16.0},
                {47.0 / 16.0, 47.0 / 16.0},
                {7.0 / 4.0, 47.0 / 16.0},
                {7.0 / 4.0, 17.0 / 16.0},
                {11.0 / 16.0, 17.0 / 16.0},
                {11.0 / 16.0, 17.0 / 8.0},
                {27.0 / 16.0, 17.0 / 8.0},
                {27.0 / 16.0, 3.0},
                {3.0, 3.0},
                {3.0, -3.0},
            }),
            partial_arc(
                std::numbers::pi,
                3.0 * std::numbers::pi / 2.0,
                false),
            2,
            Epeck::FT(3) / Epeck::FT(4),
            {
                Epeck::FT(1) / Epeck::FT(4),
                Epeck::FT(1) / Epeck::FT(2),
            },
        },
    };

    for (const SiblingFixture& fixture : fixtures) {
        const FullCircleTeaAudit2 authority =
            audit_full_circle_tea_event_exact(
                fixture.stock,
                FullCircleEventSource2::from_exact(
                    ExactCircleMotion2{
                        EPoint(0, 0),
                        EVector(2, 0),
                        fixture.motion.clockwise(),
                    },
                    audit_policy.tool_radius_mm(),
                    audit_policy.engagement_cap().chord_ratio()));
        require(
            authority.violating_parameter.has_value()
                && authority.violating_parameter->chart == 0
                && authority.violating_parameter->parameter
                    == Epeck::FT(1) / Epeck::FT(2),
            "sibling fixture full-circle authority witness is not outside");

        for (const AuditArcStationParameter2& anchor : {
                 AuditArcStationParameter2::from_start_anchor(
                     fixture.motion),
                 AuditArcStationParameter2::from_terminal_anchor(
                     fixture.motion),
             }) {
            require(
                replay_audit_unguarded_station_exact(
                    fixture.stock,
                    anchor.point(),
                    audit_policy.tool_radius_mm(),
                    audit_policy.engagement_cap().chord_ratio())
                    == AuditExactStationDisposition2::WITHIN_CAP,
                "sibling fixture has a live motion-owned anchor");
        }

        for (const Epeck::FT& parameter : fixture.safe_parameters) {
            const AuditArcStationParameter2 safe =
                AuditArcStationParameter2::from_interval(
                    fixture.motion, 0, fixture.chart, parameter);
            require(
                replay_audit_unguarded_station_exact(
                    fixture.stock,
                    safe.point(),
                    audit_policy.tool_radius_mm(),
                    audit_policy.engagement_cap().chord_ratio())
                    == AuditExactStationDisposition2::WITHIN_CAP,
                "sibling fixture has a live anchor or earlier DFS probe");
        }
        const AuditArcStationParameter2 live =
            AuditArcStationParameter2::from_interval(
                fixture.motion,
                0,
                fixture.chart,
                fixture.live_parameter);
        require_known_station_exceeds(
            fixture.stock, live.point(), audit_policy);

        const AuditDecisionWitness2 decision = certify_audit_tea_exact(
            fixture.stock,
            fixture.motion,
            audit_policy,
            decision_limits);
        require_cap_exceeded(
            fixture.stock,
            fixture.motion,
            audit_policy,
            decision_limits,
            decision);
        require(
            decision.counters().visited_nodes() == 3
                && decision.counters().deepest_level() == 1
                && decision.counters().exact_station_replays() == 6
                && decision.counters().exact_coverage_replays() == 1,
            "later-sibling witness lost truthful DFS/authority accounting");
    }
}

void clockwise_start_anchor_gate()
{
    // Production mutation caught: clockwise traversal checks only its terminal
    // or aliases the canonical start anchor as a complement seam.
    const AuditPolicy2 audit_policy = policy();
    const AuditDecisionLimits2 exhausted = forced_exhaustion_limits();
    const AuditArcMotion2 motion = partial_arc(
        std::numbers::pi / 2.0, 0.0, true);
    Stock2 start_live(rectangle(-0.6, 1.4, 0.6, 2.6), {});
    require_known_station_exceeds(
        start_live, motion.start_point(), audit_policy);
    require_cap_exceeded(
        start_live,
        motion,
        audit_policy,
        exhausted,
        certify_audit_tea_exact(
            start_live,
            motion,
            audit_policy,
            exhausted));
}

void multi_interval_owned_seam_precedence_gate()
{
    // Production mutation caught: refinement stops at an earlier safe leaf or
    // treats a later canonical t=0 seam as a complement alias.
    const AuditPolicy2 audit_policy = policy();
    const AuditDecisionLimits2 exhausted = forced_exhaustion_limits();
    Stock2 seam_live(rectangle(-5.0 / 8.0, 1.4, 5.0 / 8.0, 2.6), {});
    for (const AuditArcMotion2& motion : {
             partial_arc(0.0, std::numbers::pi, false),
             partial_arc(std::numbers::pi, 0.0, true),
         }) {
        require(
            replay_audit_unguarded_station_exact(
                seam_live,
                motion.start_point(),
                audit_policy.tool_radius_mm(),
                audit_policy.engagement_cap().chord_ratio())
                    == AuditExactStationDisposition2::WITHIN_CAP
                && replay_audit_unguarded_station_exact(
                       seam_live,
                       motion.end_point(),
                       audit_policy.tool_radius_mm(),
                       audit_policy.engagement_cap().chord_ratio())
                    == AuditExactStationDisposition2::WITHIN_CAP,
            "multi-interval seam fixture has a live anchor");
        require_known_station_exceeds(
            seam_live, EPoint(0, 2), audit_policy);
        const AuditDecisionWitness2 decision = certify_audit_tea_exact(
            seam_live,
            motion,
            audit_policy,
            exhausted);
        require_cap_exceeded(
            seam_live,
            motion,
            audit_policy,
            exhausted,
            decision);
        require(
            decision.counters().visited_nodes() == 0,
            "owned seam witness was deferred into discretionary refinement");
    }
}

} // namespace

void audit_certification_refinement_gate()
{
    mandatory_witness_precedence_gate();
    partial_arc_budget_reason_gate();
    later_sibling_witness_precedence_gate();
    clockwise_start_anchor_gate();
    multi_interval_owned_seam_precedence_gate();
}
