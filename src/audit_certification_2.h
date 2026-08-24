#pragma once

#include "audit_digest_2.h"
#include "audit_motion_identity_2.h"
#include "audit_policy_2.h"
#include "audit_stock_state_identity_2.h"
#include "stock_2.h"

#include <cstddef>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

class AuditCertificationError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class ExactOneRootCoordinateMismatchError : public AuditCertificationError {
public:
    using AuditCertificationError::AuditCertificationError;
};

class AuditSquaredSpatialFloorError : public AuditCertificationError {
public:
    using AuditCertificationError::AuditCertificationError;
};

class AuditDecisionLimitsNonFiniteInputError : public AuditCertificationError {
public:
    using AuditCertificationError::AuditCertificationError;
};

class AuditDecisionDepthLimitError : public AuditCertificationError {
public:
    using AuditCertificationError::AuditCertificationError;
};

class AuditDecisionNodeLimitError : public AuditCertificationError {
public:
    using AuditCertificationError::AuditCertificationError;
};

class AuditStationOutsideMotionError : public AuditCertificationError {
public:
    using AuditCertificationError::AuditCertificationError;
};

class AuditStationSeamOwnershipError : public AuditStationOutsideMotionError {
public:
    using AuditStationOutsideMotionError::AuditStationOutsideMotionError;
};

class AuditDecisionEvidenceReplayError : public AuditCertificationError {
public:
    using AuditCertificationError::AuditCertificationError;
};

class AuditAuthorityWitnessError : public AuditCertificationError {
public:
    using AuditCertificationError::AuditCertificationError;
};

class AuditDecisionCounterError : public AuditCertificationError {
public:
    using AuditCertificationError::AuditCertificationError;
};

class AuditUnresolvedReasonError : public AuditCertificationError {
public:
    using AuditCertificationError::AuditCertificationError;
};

class AuditExactStationToolRadiusError : public AuditCertificationError {
public:
    using AuditCertificationError::AuditCertificationError;
};

class AuditExactStationCapRatioError : public AuditCertificationError {
public:
    using AuditCertificationError::AuditCertificationError;
};

class AuditExactStationGapRatioError : public AuditCertificationError {
public:
    using AuditCertificationError::AuditCertificationError;
};

enum class AuditTeaVerdict2 {
    CERTIFIED,
    CAP_EXCEEDED,
    UNRESOLVED,
};

enum class AuditExactStationDisposition2 {
    WITHIN_CAP,
    CAP_EXCEEDED,
};

enum class AuditCertifiedCoverageKind2 {
    SEGMENT_EVENT_PARTITION,
    FULL_CIRCLE_EVENT_PARTITION,
    PARTIAL_ARC_FULL_CIRCLE_SAFE_SUPERSET,
};

enum class AuditUnresolvedReason2 {
    NODE_LIMIT_EXHAUSTED,
    DEPTH_LIMIT_EXHAUSTED,
    SPATIAL_FLOOR_REACHED,
    EVENT_AUTHORITY_UNRESOLVED,
};

class ExactOneRootPoint2 {
public:
    const Epeck::FT& root() const noexcept;

private:
    explicit ExactOneRootPoint2(Epeck::FT root);
    friend ExactOneRootPoint2 decompose_exact_one_root_point(const GpsPoint&);
    Epeck::FT root_;
};

ExactOneRootPoint2 decompose_exact_one_root_point(const GpsPoint& point);

CGAL::Sign audit_sign_mixed_radical_exact(
    const Epeck::FT& a,
    const Epeck::FT& b,
    const Epeck::FT& c,
    const Epeck::FT& d,
    const Epeck::FT& alpha,
    const Epeck::FT& beta);

class AuditExactStationClassification2 {
public:
    AuditExactStationDisposition2 disposition() const noexcept;
    const std::vector<GpsPoint>& boundary_intersections() const noexcept;
    const std::vector<std::vector<std::pair<GpsPoint, GpsPoint>>>&
    true_run_arcs() const noexcept;

private:
    AuditExactStationClassification2(
        AuditExactStationDisposition2 disposition,
        std::vector<GpsPoint> boundary_intersections,
        std::vector<std::vector<std::pair<GpsPoint, GpsPoint>>> true_run_arcs);
    friend AuditExactStationClassification2
    classify_audit_unguarded_station_exact(
        const Stock2&,
        const EPoint&,
        const Epeck::FT&,
        const Epeck::FT&,
        const Epeck::FT&);
    AuditExactStationDisposition2 disposition_;
    std::vector<GpsPoint> boundary_intersections_;
    std::vector<std::vector<std::pair<GpsPoint, GpsPoint>>> true_run_arcs_;
};

AuditExactStationClassification2 classify_audit_unguarded_station_exact(
    const Stock2& stock,
    const EPoint& center,
    const Epeck::FT& tool_radius,
    const Epeck::FT& cap_chord_ratio,
    const Epeck::FT& gap_close_ratio = Epeck::FT(0));

AuditExactStationDisposition2 replay_audit_unguarded_station_exact(
    const Stock2& stock,
    const EPoint& center,
    const Epeck::FT& tool_radius,
    const Epeck::FT& cap_chord_ratio);

class AuditSquaredSpatialFloorMm2 {
public:
    static AuditSquaredSpatialFloorMm2 build(const Epeck::FT& value);
    const Epeck::FT& value() const noexcept;
    const std::string& canonical_bytes() const noexcept;

private:
    AuditSquaredSpatialFloorMm2(Epeck::FT value, std::string canonical_bytes);
    Epeck::FT value_;
    std::string canonical_bytes_;
};

class AuditDecisionLimits2 {
public:
    static AuditDecisionLimits2 build(
        const AuditSquaredSpatialFloorMm2& squared_spatial_floor_mm,
        std::size_t max_depth,
        std::size_t max_nodes);
    const AuditSquaredSpatialFloorMm2& squared_spatial_floor_mm() const noexcept;
    std::size_t max_depth() const noexcept;
    std::size_t max_nodes() const noexcept;
    const std::string& canonical_bytes() const noexcept;

private:
    AuditDecisionLimits2(
        AuditSquaredSpatialFloorMm2 squared_spatial_floor_mm,
        std::size_t max_depth,
        std::size_t max_nodes,
        std::string canonical_bytes);
    AuditSquaredSpatialFloorMm2 squared_spatial_floor_mm_;
    std::size_t max_depth_;
    std::size_t max_nodes_;
    std::string canonical_bytes_;
};

class AuditSegmentStationParameter2 {
public:
    static AuditSegmentStationParameter2 build(
        const AuditSegmentMotion2& motion,
        const Epeck::FT& parameter);
    const NativeMotionDigest2& motion_digest() const noexcept;
    const EPoint& point() const noexcept;
    const Epeck::FT& parameter() const noexcept;
    const std::string& canonical_bytes() const noexcept;

private:
    AuditSegmentStationParameter2(
        NativeMotionDigest2 motion_digest,
        EPoint point,
        Epeck::FT parameter,
        std::string canonical_bytes);
    NativeMotionDigest2 motion_digest_;
    EPoint point_;
    Epeck::FT parameter_;
    std::string canonical_bytes_;
};

class AuditCircleStationParameter2 {
public:
    static AuditCircleStationParameter2 build(
        const AuditCircleMotion2& motion,
        int chart,
        const Epeck::FT& parameter);
    const NativeMotionDigest2& motion_digest() const noexcept;
    const EPoint& point() const noexcept;
    const std::string& canonical_bytes() const noexcept;

private:
    AuditCircleStationParameter2(
        NativeMotionDigest2 motion_digest,
        EPoint point,
        std::string canonical_bytes);
    NativeMotionDigest2 motion_digest_;
    EPoint point_;
    std::string canonical_bytes_;
};

class AuditArcStationParameter2 {
public:
    static AuditArcStationParameter2 from_interval(
        const AuditArcMotion2& motion,
        std::size_t interval_ordinal,
        int chart,
        const Epeck::FT& parameter);
    static AuditArcStationParameter2 from_start_anchor(
        const AuditArcMotion2& motion);
    static AuditArcStationParameter2 from_terminal_anchor(
        const AuditArcMotion2& motion);
    const NativeMotionDigest2& motion_digest() const noexcept;
    const EPoint& point() const noexcept;
    const std::string& canonical_bytes() const noexcept;

private:
    AuditArcStationParameter2(
        NativeMotionDigest2 motion_digest,
        EPoint point,
        std::string canonical_bytes);
    NativeMotionDigest2 motion_digest_;
    EPoint point_;
    std::string canonical_bytes_;
};

using AuditMotionStationParameter2 = std::variant<
    AuditSegmentStationParameter2,
    AuditCircleStationParameter2,
    AuditArcStationParameter2>;

class AuditDecisionCounters2 {
public:
    std::size_t visited_nodes() const noexcept;
    std::size_t deepest_level() const noexcept;
    std::size_t exact_station_replays() const noexcept;
    std::size_t exact_coverage_replays() const noexcept;

private:
    static AuditDecisionCounters2 build(
        std::size_t visited_nodes,
        std::size_t deepest_level,
        std::size_t exact_station_replays,
        std::size_t exact_coverage_replays);
    AuditDecisionCounters2(
        std::size_t visited_nodes,
        std::size_t deepest_level,
        std::size_t exact_station_replays,
        std::size_t exact_coverage_replays);
    friend class AuditDecisionWitness2;
    friend class AuditDecisionEvidenceFactory2;
    friend class AuditDecisionAccounting2;
    std::size_t visited_nodes_;
    std::size_t deepest_level_;
    std::size_t exact_station_replays_;
    std::size_t exact_coverage_replays_;
};

class AuditCertifiedCoverage2 {
public:
    AuditCertifiedCoverageKind2 kind() const noexcept;

private:
    AuditCertifiedCoverage2(
        AuditCertifiedCoverageKind2 kind,
        std::string stock_digest,
        std::string motion_digest,
        std::string policy_digest,
        std::string limits_bytes,
        std::string authority_bytes);
    friend class AuditDecisionWitness2;
    friend class AuditDecisionEvidenceFactory2;
    friend bool replay_audit_certified_coverage(
        const Stock2&,
        const AuditSegmentMotion2&,
        const AuditPolicy2&,
        const AuditDecisionLimits2&,
        const AuditCertifiedCoverage2&);
    friend bool replay_audit_certified_coverage(
        const Stock2&,
        const AuditCircleMotion2&,
        const AuditPolicy2&,
        const AuditDecisionLimits2&,
        const AuditCertifiedCoverage2&);
    friend bool replay_audit_certified_coverage(
        const Stock2&,
        const AuditArcMotion2&,
        const AuditPolicy2&,
        const AuditDecisionLimits2&,
        const AuditCertifiedCoverage2&);
    AuditCertifiedCoverageKind2 kind_;
    std::string stock_digest_;
    std::string motion_digest_;
    std::string policy_digest_;
    std::string limits_bytes_;
    std::string authority_bytes_;
};

class AuditExactStationWitness2 {
public:
    const AuditStockStateDigest2& stock_state_digest() const noexcept;
    const NativeMotionDigest2& motion_digest() const noexcept;
    const AuditPolicyDigest2& policy_digest() const noexcept;
    AuditExactStationDisposition2 disposition() const noexcept;

private:
    AuditExactStationWitness2(
        AuditStockStateDigest2 stock_state_digest,
        NativeMotionDigest2 motion_digest,
        AuditPolicyDigest2 policy_digest,
        AuditMotionStationParameter2 parameter,
        AuditExactStationDisposition2 disposition);
    friend class AuditDecisionWitness2;
    friend class AuditDecisionEvidenceFactory2;
    friend bool replay_audit_exact_station_witness(
        const Stock2&,
        const AuditSegmentMotion2&,
        const AuditPolicy2&,
        const AuditExactStationWitness2&);
    friend bool replay_audit_exact_station_witness(
        const Stock2&,
        const AuditCircleMotion2&,
        const AuditPolicy2&,
        const AuditExactStationWitness2&);
    friend bool replay_audit_exact_station_witness(
        const Stock2&,
        const AuditArcMotion2&,
        const AuditPolicy2&,
        const AuditExactStationWitness2&);
    AuditStockStateDigest2 stock_state_digest_;
    NativeMotionDigest2 motion_digest_;
    AuditPolicyDigest2 policy_digest_;
    AuditMotionStationParameter2 parameter_;
    AuditExactStationDisposition2 disposition_;
};

class AuditUnresolvedEvidence2 {
public:
    AuditUnresolvedReason2 reason() const noexcept;

private:
    class Cause;
    explicit AuditUnresolvedEvidence2(std::shared_ptr<const Cause> cause);
    friend class AuditDecisionWitness2;
    friend class AuditDecisionEvidenceFactory2;
    friend class AuditDecisionReplay2;
    AuditUnresolvedReason2 reason_;
    std::shared_ptr<const Cause> cause_;
};

class AuditDecisionWitness2 {
public:
    AuditTeaVerdict2 verdict() const noexcept;
    bool has_exact_station_witness() const noexcept;
    bool has_certified_coverage() const noexcept;
    bool has_unresolved_evidence() const noexcept;
    const AuditExactStationWitness2& exact_station_witness() const;
    const AuditCertifiedCoverage2& certified_coverage() const;
    const AuditUnresolvedEvidence2& unresolved_evidence() const;
    const AuditDecisionCounters2& counters() const noexcept;
    const std::string& canonical_bytes() const noexcept;
    const NativeDecisionDigest2& digest() const noexcept;

private:
    AuditDecisionWitness2(
        AuditTeaVerdict2 verdict,
        std::optional<AuditExactStationWitness2> station,
        std::optional<AuditCertifiedCoverage2> coverage,
        std::optional<AuditUnresolvedEvidence2> unresolved,
        AuditDecisionCounters2 counters,
        std::string stock_digest,
        std::string motion_digest,
        std::string policy_digest,
        std::string limits_bytes,
        std::string canonical_bytes);
    friend class AuditDecisionEvidenceFactory2;
    friend class AuditDecisionReplay2;
    AuditTeaVerdict2 verdict_;
    std::optional<AuditExactStationWitness2> station_;
    std::optional<AuditCertifiedCoverage2> coverage_;
    std::optional<AuditUnresolvedEvidence2> unresolved_;
    AuditDecisionCounters2 counters_;
    std::string stock_digest_;
    std::string motion_digest_;
    std::string policy_digest_;
    std::string limits_bytes_;
    std::string canonical_bytes_;
    NativeDecisionDigest2 digest_;
};

bool replay_audit_exact_station_witness(
    const Stock2& stock,
    const AuditSegmentMotion2& motion,
    const AuditPolicy2& policy,
    const AuditExactStationWitness2& witness);
bool replay_audit_exact_station_witness(
    const Stock2& stock,
    const AuditCircleMotion2& motion,
    const AuditPolicy2& policy,
    const AuditExactStationWitness2& witness);
bool replay_audit_exact_station_witness(
    const Stock2& stock,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy,
    const AuditExactStationWitness2& witness);

bool replay_audit_certified_coverage(
    const Stock2& stock,
    const AuditSegmentMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const AuditCertifiedCoverage2& coverage);
bool replay_audit_certified_coverage(
    const Stock2& stock,
    const AuditCircleMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const AuditCertifiedCoverage2& coverage);
bool replay_audit_certified_coverage(
    const Stock2& stock,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const AuditCertifiedCoverage2& coverage);

bool audit_decision_witness_is_self_consistent(
    const Stock2& stock,
    const AuditSegmentMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const AuditDecisionWitness2& decision);
bool audit_decision_witness_is_self_consistent(
    const Stock2& stock,
    const AuditCircleMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const AuditDecisionWitness2& decision);
bool audit_decision_witness_is_self_consistent(
    const Stock2& stock,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const AuditDecisionWitness2& decision);

AuditDecisionWitness2 certify_audit_tea_exact(
    const Stock2& stock,
    const AuditSegmentMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits);
AuditDecisionWitness2 certify_audit_tea_exact(
    const Stock2& stock,
    const AuditCircleMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits);
AuditDecisionWitness2 certify_audit_tea_exact(
    const Stock2& stock,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits);
