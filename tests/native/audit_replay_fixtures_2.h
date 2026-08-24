#pragma once

#include "audit_replay_2.h"
#include "audit_replay_internal_2.h"

#include "audit_classification_core_2.h"
#include "audit_policy_2.h"
#include "audit_request_identity_2.h"
#include "audit_stock_identity_2.h"
#include "continuous_tea_2/sha256.h"

#include <numbers>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace audit_replay_fixtures {

inline void require(bool condition, const char* message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

template <class Error, class Function>
void require_throws(Function&& function, const char* message)
{
    try {
        std::forward<Function>(function)();
    } catch (const Error&) {
        return;
    }
    throw std::runtime_error(message);
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

inline compas::RowMatrixXd square_hole(
    double x_min = -3.0,
    double y_min = -3.0,
    double x_max = -2.0,
    double y_max = -2.0)
{
    compas::RowMatrixXd hole(4, 2);
    hole << x_min, y_min,
        x_min, y_max,
        x_max, y_max,
        x_max, y_min;
    return hole;
}

inline compas::RowMatrixXd reversed_rotated_ring(
    const compas::RowMatrixXd& ring)
{
    compas::RowMatrixXd equivalent(ring.rows(), ring.cols());
    for (Eigen::Index row = 0; row < ring.rows(); ++row) {
        const Eigen::Index source =
            (ring.rows() - row) % ring.rows();
        equivalent.row(row) = ring.row(source);
    }
    return equivalent;
}

inline AuditPolicy2 policy(std::size_t center_limit = 4096)
{
    constexpr double CAP = std::numbers::pi / 2.0;
    return AuditPolicy2::build(
        AuditCapObservation2::build(CAP, audit_cap_chord_ratio(CAP)),
        Epeck::FT(1) / Epeck::FT(2),
        Epeck::FT(1) / Epeck::FT(16),
        center_limit);
}

inline AuditDecisionLimits2 decision_limits(
    const Epeck::FT& squared_floor = Epeck::FT(1) / Epeck::FT(4096),
    std::size_t max_depth = 8,
    std::size_t max_nodes = 256)
{
    return AuditDecisionLimits2::build(
        AuditSquaredSpatialFloorMm2::build(squared_floor),
        max_depth,
        max_nodes);
}

inline std::size_t exact_evidence_count(
    const AuditDecisionCounters2& counters) noexcept
{
    // Refinement nodes are work accounting, not additional evidence: every
    // examined node already contributes its exact station replay.
    return counters.exact_station_replays()
        + counters.exact_coverage_replays();
}

inline AuditSegmentMotion2 segment(double y = 0.0)
{
    return std::get<AuditSegmentMotion2>(classify_audit_line(
        {-1.0, y, 0.0}, {1.0, y, 0.0}, 0.0, 5.0, "cut"));
}

inline AuditCircleMotion2 circle(bool clockwise = false)
{
    return std::get<AuditCircleMotion2>(classify_audit_circle(
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        1.0,
        clockwise,
        0.0,
        5.0,
        "cut"));
}

inline AuditArcMotion2 arc()
{
    return std::get<AuditArcMotion2>(classify_audit_arc(
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        1.0,
        0.0,
        std::numbers::pi / 2.0,
        false,
        0.0,
        5.0,
        "cut"));
}

inline AuditVerticalPlunge2 plunge()
{
    return std::get<AuditVerticalPlunge2>(classify_audit_line(
        {3.0, 3.0, 5.0}, {3.0, 3.0, 0.0}, 0.0, 5.0, "plunge"));
}

inline AuditVerticalRetract2 retract()
{
    return std::get<AuditVerticalRetract2>(classify_audit_line(
        {0.0, 0.0, 0.0}, {0.0, 0.0, 5.0}, 0.0, 5.0, "retract"));
}

inline AuditClearanceTransport2 clearance()
{
    return std::get<AuditClearanceTransport2>(classify_audit_line(
        {-1.0, 0.0, 5.0}, {1.0, 0.0, 5.0}, 0.0, 5.0, "link"));
}

inline AuthenticatedOperationDigest2 operation_digest(
    const std::string& label)
{
    return AuthenticatedOperationDigestAuthority2::from_external_bytes(
        sha256_bytes("audit-replay-operation:" + label));
}

inline AuditInputDigest2 input_digest(const std::string& label = "baseline")
{
    return AuditInputDigestAuthority2::from_external_bytes(
        sha256_bytes("audit-replay-input:" + label));
}

template <class Motion>
AuditNativeRequestIdentity2 one_motion_request(
    const compas::RowMatrixXd& boundary,
    const AuditPolicy2& audit_policy,
    const AuditDecisionLimits2& limits,
    const Motion& motion)
{
    return AuditNativeRequestIdentity2::build(
        AuditNativeStockIdentity2::build(boundary, {}),
        audit_policy,
        limits,
        {motion.digest()});
}

template <class Motion>
AuditReplay2 one_motion_replay(
    const compas::RowMatrixXd& boundary,
    const AuditPolicy2& audit_policy,
    const AuditDecisionLimits2& limits,
    const Motion& motion,
    const AuthenticatedOperationDigest2& digest)
{
    return AuditReplay2::build(
        boundary,
        {},
        input_digest(),
        one_motion_request(boundary, audit_policy, limits, motion),
        audit_policy,
        limits,
        {digest});
}

} // namespace audit_replay_fixtures
