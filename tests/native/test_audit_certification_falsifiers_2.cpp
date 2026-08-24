#include "audit_certification_fixtures_2.h"

#include <numbers>

namespace {

using namespace audit_certification_fixtures;

Stock2 machined_rib_stock()
{
    Stock2 stock(rectangle(-6.0, -6.0, 6.0, 6.0), {});
    stock.subtract_disk(
        0.0,
        0.0,
        TOOL_RADIUS - 0.5 * RIB_THICKNESS);
    const Epeck::FT guide_radius = Epeck::FT(501) / Epeck::FT(500);
    const AuditArcMotion2 motion = AuditArcMotion2::build(
        EPoint(0, 0),
        EVector(guide_radius, Epeck::FT(0)),
        guide_radius,
        0.0,
        2.0 * std::numbers::pi,
        false,
        Epeck::FT(0));
    stock.subtract_exact_arc(
        motion,
        Epeck::FT(1) / Epeck::FT(2),
        Epeck::FT(49) / Epeck::FT(100),
        32);
    return stock;
}

AuditArcMotion2 quarter_arc_through(const EPoint& live_point)
{
    const Epeck::FT scale = Epeck::FT(1) / Epeck::FT(32);
    return AuditArcMotion2::build(
        EPoint(
            live_point.x()
                + Epeck::FT(7) * scale / Epeck::FT(5),
            live_point.y()
                - Epeck::FT(24) * scale / Epeck::FT(5)),
        EVector(Epeck::FT(3) * scale, Epeck::FT(4) * scale),
        Epeck::FT(5) * scale,
        0.0,
        std::numbers::pi / 2.0,
        false,
        Epeck::FT(0));
}

AuditArcMotion2 verified_arc_falsifier_motion(
    const Stock2& stock,
    const EPoint& live_point,
    const AuditPolicy2& audit_policy)
{
    // Fixture-validity evidence only. Callers must separately invoke the
    // certification adapter when this is intended as adapter coverage.
    const AuditArcMotion2 motion = quarter_arc_through(live_point);
    const AuditArcStationParameter2 midpoint =
        AuditArcStationParameter2::from_interval(
            motion, 0, 0, Epeck::FT(1) / Epeck::FT(2));
    require(
        midpoint.point() == live_point,
        "arc falsifier midpoint is not the independently derived live point");
    const AuditExactStationDisposition2 start_disposition =
        replay_audit_unguarded_station_exact(
            stock,
            motion.start_point(),
            audit_policy.tool_radius_mm(),
            audit_policy.engagement_cap().chord_ratio());
    const AuditExactStationDisposition2 end_disposition =
        replay_audit_unguarded_station_exact(
            stock,
            motion.end_point(),
            audit_policy.tool_radius_mm(),
            audit_policy.engagement_cap().chord_ratio());
    require(
        start_disposition == AuditExactStationDisposition2::WITHIN_CAP
            && end_disposition == AuditExactStationDisposition2::WITHIN_CAP,
        "arc falsifier lacks exact safe anchors");
    require_known_station_exceeds(stock, midpoint.point(), audit_policy);
    return motion;
}

void require_arc_falsifier(
    const Stock2& stock,
    const EPoint& live_point,
    const AuditPolicy2& audit_policy)
{
    const AuditArcMotion2 motion = verified_arc_falsifier_motion(
        stock, live_point, audit_policy);
    const AuditDecisionLimits2 decision_limits = limits(
        Epeck::FT(1) / Epeck::FT(4096), 4, 31);
    const AuditDecisionWitness2 decision = certify_audit_tea_exact(
        stock, motion, audit_policy, decision_limits);
    require_cap_exceeded(
        stock, motion, audit_policy, decision_limits, decision);
}

void rational_rotation_similarity_gate()
{
    // Production mutation caught: axis-specialized closure changes exact
    // station disposition under integer rotation plus scale.
    for (const RationalRotation2& rotation : RATIONAL_ROTATIONS) {
        require(
            rotation.a * rotation.a + rotation.b * rotation.b
                == rotation.scale * rotation.scale,
            "rotation fixture is not an integer Pythagorean similarity");
        const AuditDecisionLimits2 decision_limits = limits(
            Epeck::FT(rotation.scale * rotation.scale)
            / Epeck::FT(4096));
        const AuditPolicy2 scaled_policy = policy(
            std::numbers::pi / 2.0,
            Epeck::FT(rotation.scale) / Epeck::FT(2),
            Epeck::FT(rotation.scale) / Epeck::FT(50));

        const compas::RowMatrixXd safe_boundary = transform_polygon(
            rectangle(-10.0, -10.0, 10.0, -3.0 / 8.0),
            rotation);
        Stock2 safe_stock(safe_boundary, {});
        const auto [safe_x0, safe_y0] = transform_point(-1.0, 0.0, rotation);
        const auto [safe_x1, safe_y1] = transform_point(1.0, 0.0, rotation);
        const AuditSegmentMotion2 safe_motion = segment(
            safe_x0,
            safe_y0,
            safe_x1,
            safe_y1);

        require(
            replay_audit_unguarded_station_exact(
                safe_stock,
                EPoint(0, 0),
                scaled_policy.tool_radius_mm(),
                scaled_policy.engagement_cap().chord_ratio())
                == AuditExactStationDisposition2::WITHIN_CAP,
            "integer-rotated safe station changed exact disposition");
        const AuditDecisionWitness2 certified = certify_audit_tea_exact(
            safe_stock,
            safe_motion,
            scaled_policy,
            decision_limits);
        require_certified(
            safe_stock,
            safe_motion,
            scaled_policy,
            decision_limits,
            certified,
            AuditCertifiedCoverageKind2::SEGMENT_EVENT_PARTITION);

        Stock2 falsifier = dyadic_square_diamond_rib_stock(rotation);
        const auto [half_x, half_y] = transform_point(
            3.0 / 32.0,
            0.0,
            rotation);
        const AuditSegmentMotion2 falsifier_motion = segment(
            -half_x,
            -half_y,
            half_x,
            half_y);
        const EPoint exact_midpoint(
            (falsifier_motion.xy().start.x()
                + falsifier_motion.xy().end.x())
                / Epeck::FT(2),
            (falsifier_motion.xy().start.y()
                + falsifier_motion.xy().end.y())
                / Epeck::FT(2));
        require(
            replay_audit_unguarded_station_exact(
                falsifier,
                falsifier_motion.xy().start,
                scaled_policy.tool_radius_mm(),
                scaled_policy.engagement_cap().chord_ratio())
                == AuditExactStationDisposition2::WITHIN_CAP
                && replay_audit_unguarded_station_exact(
                       falsifier,
                       exact_midpoint,
                       scaled_policy.tool_radius_mm(),
                       scaled_policy.engagement_cap().chord_ratio())
                    == AuditExactStationDisposition2::CAP_EXCEEDED
                && replay_audit_unguarded_station_exact(
                       falsifier,
                       falsifier_motion.xy().end,
                       scaled_policy.tool_radius_mm(),
                       scaled_policy.engagement_cap().chord_ratio())
                    == AuditExactStationDisposition2::WITHIN_CAP,
            "integer-rotated falsifier lost safe endpoints or live interior");
        const AuditDecisionWitness2 exceeded = certify_audit_tea_exact(
            falsifier,
            falsifier_motion,
            scaled_policy,
            decision_limits);
        require_cap_exceeded(
            falsifier,
            falsifier_motion,
            scaled_policy,
            decision_limits,
            exceeded);
    }
}

void adapter_falsifier_gate()
{
    // Production mutation caught: endpoint-only or concentric-only closure
    // misses annular, sector, machined, or spiral interior witnesses.
    const AuditPolicy2 audit_policy = policy();
    const AuditDecisionLimits2 decision_limits = limits();

    const auto require_segment_falsifier = [&audit_policy, &decision_limits](
                                               const Stock2& stock) {
        const AuditSegmentMotion2 motion = short_segment();
        const EPoint midpoint(
            (motion.xy().start.x() + motion.xy().end.x()) / Epeck::FT(2),
            (motion.xy().start.y() + motion.xy().end.y()) / Epeck::FT(2));
        require(
            replay_audit_unguarded_station_exact(
                stock,
                motion.xy().start,
                audit_policy.tool_radius_mm(),
                audit_policy.engagement_cap().chord_ratio())
                    == AuditExactStationDisposition2::WITHIN_CAP
                && replay_audit_unguarded_station_exact(
                       stock,
                       midpoint,
                       audit_policy.tool_radius_mm(),
                       audit_policy.engagement_cap().chord_ratio())
                    == AuditExactStationDisposition2::CAP_EXCEEDED
                && replay_audit_unguarded_station_exact(
                       stock,
                       motion.xy().end,
                       audit_policy.tool_radius_mm(),
                       audit_policy.engagement_cap().chord_ratio())
                    == AuditExactStationDisposition2::WITHIN_CAP,
            "segment falsifier lacks exact safe endpoints and live midpoint");
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
    };
    require_segment_falsifier(annular_rib_stock());
    require_segment_falsifier(Stock2(sector_rib_polygon(), {}));
    const Stock2 machined = machined_rib_stock();
    require_segment_falsifier(machined);

    require_arc_falsifier(
        dyadic_square_diamond_rib_stock(RATIONAL_ROTATIONS[0]),
        EPoint(0, 0),
        audit_policy);
    // The exact machined fixture remains Task 3 depletion + station evidence.
    // Its full-circle closure exceeded the 120 s bounded native gate and is
    // deliberately not claimed as arc-adapter coverage.
    static_cast<void>(verified_arc_falsifier_motion(
        machined, EPoint(0, 0), audit_policy));

    Stock2 spiral(spiral_rib_polygon(), {});
    constexpr double segment_center_y = 0.0034;
    constexpr double half_dx = 0.0075;
    constexpr double half_dy = 0.01;
    const AuditSegmentMotion2 oblique = segment(
        -half_dx,
        segment_center_y - half_dy,
        half_dx,
        segment_center_y + half_dy);
    const EPoint exact_midpoint(
        (oblique.xy().start.x() + oblique.xy().end.x()) / Epeck::FT(2),
        (oblique.xy().start.y() + oblique.xy().end.y()) / Epeck::FT(2));
    require(
        replay_audit_unguarded_station_exact(
            spiral,
            oblique.xy().start,
            audit_policy.tool_radius_mm(),
            audit_policy.engagement_cap().chord_ratio())
                == AuditExactStationDisposition2::WITHIN_CAP
            && replay_audit_unguarded_station_exact(
                   spiral,
                   exact_midpoint,
                   audit_policy.tool_radius_mm(),
                   audit_policy.engagement_cap().chord_ratio())
                == AuditExactStationDisposition2::CAP_EXCEEDED
            && replay_audit_unguarded_station_exact(
                   spiral,
                   oblique.xy().end,
                   audit_policy.tool_radius_mm(),
                   audit_policy.engagement_cap().chord_ratio())
                == AuditExactStationDisposition2::WITHIN_CAP,
        "spiral falsifier lacks exact safe endpoints and live midpoint");
    const AuditDecisionWitness2 spiral_decision = certify_audit_tea_exact(
        spiral,
        oblique,
        audit_policy,
        decision_limits);
    require_cap_exceeded(
        spiral,
        oblique,
        audit_policy,
        decision_limits,
        spiral_decision);
    // The spiral keeps its exact segment-adapter verdict and arc station
    // validity evidence; its full-circle closure likewise exceeded 120 s.
    static_cast<void>(verified_arc_falsifier_motion(
        spiral, exact_midpoint, audit_policy));
}

} // namespace

void audit_certification_falsifiers_gate()
{
    rational_rotation_similarity_gate();
    adapter_falsifier_gate();
}
