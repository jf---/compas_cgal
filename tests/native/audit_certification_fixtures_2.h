#pragma once

#include "audit_certification_2.h"
#include "audit_classification_core_2.h"
#include "audit_policy_2.h"
#include "stock_2.h"

#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <variant>

namespace audit_certification_fixtures {

inline constexpr double TOOL_RADIUS = 0.5;
inline constexpr double RIB_THICKNESS = 0.004;
inline constexpr std::size_t RIB_FACETS = 256;
inline constexpr double SHORT_MOTION = 0.025;
inline constexpr double SPIRAL_GROWTH = 0.005;
inline constexpr double SPIRAL_HALF_TURN = 2.4;
inline constexpr std::size_t SPIRAL_STEPS = 48;
inline constexpr double SECTOR_EXTENT = 0.75 * std::numbers::pi;
inline constexpr std::size_t SECTOR_FACETS = 96;

struct RationalRotation2 {
    int a;
    int b;
    int scale;
};

inline constexpr RationalRotation2 RATIONAL_ROTATIONS[] {
    {1, 0, 1},
    {3, 4, 5},
    {12, 5, 13},
    {0, 1, 1},
};

inline void require(bool condition, const char* message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

inline compas::RowMatrixXd rectangle(
    double x_min,
    double y_min,
    double x_max,
    double y_max)
{
    compas::RowMatrixXd boundary(4, 2);
    boundary << x_min, y_min,
        x_max, y_min,
        x_max, y_max,
        x_min, y_max;
    return boundary;
}

inline std::pair<double, double> transform_point(
    double x,
    double y,
    const RationalRotation2& rotation)
{
    return {
        static_cast<double>(rotation.a) * x
            - static_cast<double>(rotation.b) * y,
        static_cast<double>(rotation.b) * x
            + static_cast<double>(rotation.a) * y,
    };
}

inline compas::RowMatrixXd transform_polygon(
    const compas::RowMatrixXd& polygon,
    const RationalRotation2& rotation)
{
    compas::RowMatrixXd transformed(polygon.rows(), 2);
    for (Eigen::Index row = 0; row < polygon.rows(); ++row) {
        const auto [x, y] = transform_point(
            polygon(row, 0),
            polygon(row, 1),
            rotation);
        transformed(row, 0) = x;
        transformed(row, 1) = y;
    }
    return transformed;
}

inline compas::RowMatrixXd regular_polygon(double radius)
{
    compas::RowMatrixXd boundary(
        static_cast<Eigen::Index>(RIB_FACETS), 2);
    for (std::size_t index = 0; index < RIB_FACETS; ++index) {
        const double angle = 2.0 * std::numbers::pi
            * static_cast<double>(index)
            / static_cast<double>(RIB_FACETS);
        boundary(static_cast<Eigen::Index>(index), 0)
            = radius * std::cos(angle);
        boundary(static_cast<Eigen::Index>(index), 1)
            = radius * std::sin(angle);
    }
    return boundary;
}

inline compas::RowMatrixXd sector_rib_polygon()
{
    const double outer = TOOL_RADIUS + 0.5 * RIB_THICKNESS;
    const double inner = TOOL_RADIUS - 0.5 * RIB_THICKNESS;
    compas::RowMatrixXd boundary(
        static_cast<Eigen::Index>(2 * (SECTOR_FACETS + 1)), 2);
    for (std::size_t index = 0; index <= SECTOR_FACETS; ++index) {
        const double angle = -0.5 * SECTOR_EXTENT
            + SECTOR_EXTENT * static_cast<double>(index)
                / static_cast<double>(SECTOR_FACETS);
        boundary(static_cast<Eigen::Index>(index), 0)
            = outer * std::cos(angle);
        boundary(static_cast<Eigen::Index>(index), 1)
            = outer * std::sin(angle);
        const std::size_t reverse = SECTOR_FACETS - index;
        const double inner_angle = -0.5 * SECTOR_EXTENT
            + SECTOR_EXTENT * static_cast<double>(reverse)
                / static_cast<double>(SECTOR_FACETS);
        boundary(
            static_cast<Eigen::Index>(SECTOR_FACETS + 1 + index),
            0) = inner * std::cos(inner_angle);
        boundary(
            static_cast<Eigen::Index>(SECTOR_FACETS + 1 + index),
            1) = inner * std::sin(inner_angle);
    }
    return boundary;
}

inline compas::RowMatrixXd spiral_rib_polygon()
{
    compas::RowMatrixXd boundary(
        static_cast<Eigen::Index>(2 * (SPIRAL_STEPS + 1)), 2);
    for (std::size_t index = 0; index <= SPIRAL_STEPS; ++index) {
        const double angle = -SPIRAL_HALF_TURN
            + 2.0 * SPIRAL_HALF_TURN * static_cast<double>(index)
                / static_cast<double>(SPIRAL_STEPS);
        const double outer = TOOL_RADIUS
            + SPIRAL_GROWTH * angle
            + 0.5 * RIB_THICKNESS;
        boundary(static_cast<Eigen::Index>(index), 0)
            = outer * std::cos(angle);
        boundary(static_cast<Eigen::Index>(index), 1)
            = outer * std::sin(angle);

        const std::size_t reverse = SPIRAL_STEPS - index;
        const double inner_angle = -SPIRAL_HALF_TURN
            + 2.0 * SPIRAL_HALF_TURN * static_cast<double>(reverse)
                / static_cast<double>(SPIRAL_STEPS);
        const double inner = TOOL_RADIUS
            + SPIRAL_GROWTH * inner_angle
            - 0.5 * RIB_THICKNESS;
        boundary(
            static_cast<Eigen::Index>(SPIRAL_STEPS + 1 + index),
            0) = inner * std::cos(inner_angle);
        boundary(
            static_cast<Eigen::Index>(SPIRAL_STEPS + 1 + index),
            1) = inner * std::sin(inner_angle);
    }
    return boundary;
}

inline Stock2 annular_rib_stock(const RationalRotation2& rotation)
{
    return Stock2(
        transform_polygon(
            regular_polygon(TOOL_RADIUS + 0.5 * RIB_THICKNESS),
            rotation),
        {transform_polygon(
            regular_polygon(TOOL_RADIUS - 0.5 * RIB_THICKNESS),
            rotation)});
}

inline Stock2 annular_rib_stock()
{
    return annular_rib_stock(RATIONAL_ROTATIONS[0]);
}

inline Stock2 dyadic_square_diamond_rib_stock(
    const RationalRotation2& rotation)
{
    constexpr double outer_half_width = 65.0 / 128.0;
    constexpr double inner_axis_radius = 63.0 / 128.0;
    compas::RowMatrixXd outer(4, 2);
    outer << -outer_half_width, -outer_half_width,
        outer_half_width, -outer_half_width,
        outer_half_width, outer_half_width,
        -outer_half_width, outer_half_width;
    compas::RowMatrixXd inner(4, 2);
    inner << inner_axis_radius, 0.0,
        0.0, inner_axis_radius,
        -inner_axis_radius, 0.0,
        0.0, -inner_axis_radius;
    return Stock2(
        transform_polygon(outer, rotation),
        {transform_polygon(inner, rotation)});
}

inline AuditPolicy2 policy(
    double cap_radians = std::numbers::pi / 2.0,
    const Epeck::FT& tool_radius = Epeck::FT(1) / Epeck::FT(2),
    const Epeck::FT& depletion_chord = Epeck::FT(1) / Epeck::FT(50))
{
    return AuditPolicy2::build(
        AuditCapObservation2::build(
            cap_radians,
            audit_cap_chord_ratio(cap_radians)),
        tool_radius,
        depletion_chord,
        4096);
}

inline AuditDecisionLimits2 limits(
    const Epeck::FT& squared_floor = Epeck::FT(1) / Epeck::FT(4096),
    std::size_t max_depth = 16,
    std::size_t max_nodes = 8192)
{
    return AuditDecisionLimits2::build(
        AuditSquaredSpatialFloorMm2::build(squared_floor),
        max_depth,
        max_nodes);
}

inline AuditDecisionLimits2 forced_exhaustion_limits()
{
    return limits(Epeck::FT(1) / Epeck::FT(4096), 0, 1);
}

inline AuditSegmentMotion2 segment(
    double x0,
    double y0,
    double x1,
    double y1)
{
    return std::get<AuditSegmentMotion2>(classify_audit_line(
        {x0, y0, 0.0},
        {x1, y1, 0.0},
        0.0,
        5.0,
        "cut"));
}

inline AuditSegmentMotion2 short_segment()
{
    return segment(
        -0.5 * SHORT_MOTION,
        0.0,
        0.5 * SHORT_MOTION,
        0.0);
}

inline AuditCircleMotion2 full_circle(bool clockwise)
{
    return std::get<AuditCircleMotion2>(classify_audit_circle(
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        2.0,
        clockwise,
        0.0,
        5.0,
        "cut"));
}

inline AuditArcMotion2 partial_arc(
    double start_angle,
    double end_angle,
    bool clockwise)
{
    return AuditArcMotion2::build(
        EPoint(0, 0),
        EVector(2, 0),
        Epeck::FT(2),
        start_angle,
        end_angle,
        clockwise,
        Epeck::FT(0));
}

inline AuditArcMotion2 full_turn_arc(bool clockwise)
{
    return partial_arc(
        0.0,
        clockwise ? -2.0 * std::numbers::pi : 2.0 * std::numbers::pi,
        clockwise);
}

template <class Motion>
void require_certified(
    const Stock2& stock,
    const Motion& motion,
    const AuditPolicy2& audit_policy,
    const AuditDecisionLimits2& decision_limits,
    const AuditDecisionWitness2& decision,
    AuditCertifiedCoverageKind2 expected_kind)
{
    require(
        decision.verdict() == AuditTeaVerdict2::CERTIFIED,
        "complete exact coverage did not produce CERTIFIED");
    require(
        !decision.has_exact_station_witness()
            && decision.has_certified_coverage(),
        "CERTIFIED decision does not carry exactly one coverage proof");
    const AuditCertifiedCoverage2& coverage = decision.certified_coverage();
    require(
        coverage.kind() == expected_kind,
        "CERTIFIED decision carries the wrong closed coverage proof kind");
    require(
        replay_audit_certified_coverage(
            stock,
            motion,
            audit_policy,
            decision_limits,
            coverage),
        "CERTIFIED coverage proof did not replay against exact authority");
    require(
        decision.counters().exact_coverage_replays() > 0,
        "CERTIFIED decision records no exact coverage replay");
    require(
        audit_decision_witness_is_self_consistent(
            stock,
            motion,
            audit_policy,
            decision_limits,
            decision),
        "CERTIFIED decision did not replay as a whole");
}

template <class Motion>
void require_cap_exceeded(
    const Stock2& stock,
    const Motion& motion,
    const AuditPolicy2& audit_policy,
    const AuditDecisionLimits2& decision_limits,
    const AuditDecisionWitness2& decision)
{
    require(
        decision.verdict() == AuditTeaVerdict2::CAP_EXCEEDED,
        "live exact station did not produce CAP_EXCEEDED");
    require(
        decision.has_exact_station_witness()
            && !decision.has_certified_coverage(),
        "CAP_EXCEEDED does not carry exactly one station witness");
    const AuditExactStationWitness2& witness =
        decision.exact_station_witness();
    const AuditStockStateIdentity2 stock_identity =
        AuditStockStateIdentity2::build(stock);
    require(
        witness.stock_state_digest().bytes()
                == stock_identity.digest().bytes()
            && witness.motion_digest().bytes() == motion.digest().bytes()
            && witness.policy_digest().bytes() == audit_policy.digest().bytes()
            && witness.disposition()
                == AuditExactStationDisposition2::CAP_EXCEEDED,
        "station witness omits stock, motion, policy, or disposition identity");
    require(
        replay_audit_exact_station_witness(
            stock,
            motion,
            audit_policy,
            witness),
        "CAP_EXCEEDED witness did not replay at its motion-owned parameter");
    require(
        decision.counters().exact_station_replays() > 0,
        "CAP_EXCEEDED decision records no exact station replay");
    require(
        audit_decision_witness_is_self_consistent(
            stock,
            motion,
            audit_policy,
            decision_limits,
            decision),
        "CAP_EXCEEDED decision did not replay as a whole");
}

template <class Motion>
void require_unresolved(
    const Stock2& stock,
    const Motion& motion,
    const AuditPolicy2& audit_policy,
    const AuditDecisionLimits2& decision_limits,
    const AuditDecisionWitness2& decision,
    AuditUnresolvedReason2 expected_reason)
{
    require(
        decision.verdict() == AuditTeaVerdict2::UNRESOLVED,
        "incomplete authority did not fail closed to UNRESOLVED");
    require(
        !decision.has_exact_station_witness()
            && !decision.has_certified_coverage()
            && decision.has_unresolved_evidence(),
        "UNRESOLVED decision fabricated witness or coverage evidence");
    require(
        decision.unresolved_evidence().reason() == expected_reason,
        "UNRESOLVED decision lost its typed reason");
    require(
        decision.counters().visited_nodes() <= decision_limits.max_nodes()
            && decision.counters().deepest_level()
                <= decision_limits.max_depth(),
        "UNRESOLVED counters exceed the sealed integer budgets");
    require(
        audit_decision_witness_is_self_consistent(
            stock,
            motion,
            audit_policy,
            decision_limits,
            decision),
        "UNRESOLVED evidence did not replay as a whole");
}

inline void require_known_station_exceeds(
    const Stock2& stock,
    const EPoint& center,
    const AuditPolicy2& audit_policy)
{
    require(
        replay_audit_unguarded_station_exact(
            stock,
            center,
            audit_policy.tool_radius_mm(),
            audit_policy.engagement_cap().chord_ratio())
            == AuditExactStationDisposition2::CAP_EXCEEDED,
        "falsifier fixture lacks its independently replayed exact station");
}

} // namespace audit_certification_fixtures
