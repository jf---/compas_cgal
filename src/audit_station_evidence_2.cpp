#include "audit_certification_internal_2.h"

#include "audit_strategy_versions_2.h"
#include "canonical_encoding.h"

#include <utility>

namespace {
std::string disposition_text(AuditExactStationDisposition2 disposition)
{
    return disposition == AuditExactStationDisposition2::CAP_EXCEEDED
        ? "cap-exceeded"
        : "within-cap";
}
std::string point_bytes(const EPoint& point)
{
    return canonical_encode_tagged_union(
        "audit-exact-station-point-v1",
        canonical_encode_component_map({
            {"x", canonical_audit_rational_bytes(point.x())},
            {"y", canonical_audit_rational_bytes(point.y())},
        }));
}

const EPoint& station_point(const AuditMotionStationParameter2& parameter)
{
    return std::visit(
        [](const auto& typed) -> const EPoint& { return typed.point(); },
        parameter);
}

const NativeMotionDigest2& station_motion_digest(
    const AuditMotionStationParameter2& parameter)
{
    return std::visit(
        [](const auto& typed) -> const NativeMotionDigest2& {
            return typed.motion_digest();
        },
        parameter);
}

const std::string& station_parameter_bytes(
    const AuditMotionStationParameter2& parameter)
{
    return std::visit(
        [](const auto& typed) -> const std::string& {
            return typed.canonical_bytes();
        },
        parameter);
}

} // namespace

AuditDecisionWitness2 AuditDecisionEvidenceFactory2::make_station_decision(
    const AuditStockStateIdentity2& stock,
    const NativeMotionDigest2& motion_digest,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    AuditExactStationWitness2 witness,
    AuditDecisionCounters2 counters)
{
    validate_counters(limits, counters);
    const std::string evidence = canonical_encode_tagged_union(
        "audit-exact-station-witness-v1",
        canonical_encode_component_map({
            {"disposition", canonical_encode_bytes(
                 disposition_text(witness.disposition_))},
            {"motion-digest", witness.motion_digest_.bytes()},
            {"parameter", station_parameter_bytes(witness.parameter_)},
            {"point", point_bytes(station_point(witness.parameter_))},
            {"policy-digest", witness.policy_digest_.bytes()},
            {"stock-state-digest", witness.stock_state_digest_.bytes()},
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
            {"verdict", canonical_encode_bytes("cap-exceeded")},
        }));
    return AuditDecisionWitness2(
        AuditTeaVerdict2::CAP_EXCEEDED,
        std::move(witness),
        std::nullopt,
        std::nullopt,
        std::move(counters),
        stock.digest().bytes(),
        motion_digest.bytes(),
        policy.digest().bytes(),
        limits.canonical_bytes(),
        canonical);
}

AuditExactStationWitness2 AuditDecisionEvidenceFactory2::replay_station(
    const AuditStockStateIdentity2& stock_identity,
    const NativeMotionDigest2& motion_digest,
    const AuditPolicy2& policy,
    AuditMotionStationParameter2 parameter,
    AuditExactStationDisposition2 claimed_disposition,
    const Stock2& stock)
{
    if (claimed_disposition
        != AuditExactStationDisposition2::CAP_EXCEEDED) {
        throw AuditDecisionEvidenceReplayError(
            "only an exact CAP_EXCEEDED replay can mint a CAP decision");
    }
    if (station_motion_digest(parameter).bytes()
        != motion_digest.bytes()) {
        throw AuditDecisionEvidenceReplayError(
            "station parameter belongs to a foreign motion");
    }
    const AuditExactStationDisposition2 actual =
        replay_audit_unguarded_station_exact(
            stock,
            station_point(parameter),
            policy.tool_radius_mm(),
            policy.engagement_cap().chord_ratio());
    if (actual != claimed_disposition) {
        throw AuditDecisionEvidenceReplayError(
            "claimed exact station disposition does not replay");
    }
    return AuditExactStationWitness2(
        stock_identity.digest(),
        motion_digest,
        policy.digest(),
        std::move(parameter),
        actual);
}

AuditExactStationWitness2::AuditExactStationWitness2(
    AuditStockStateDigest2 stock_state_digest,
    NativeMotionDigest2 motion_digest,
    AuditPolicyDigest2 policy_digest,
    AuditMotionStationParameter2 parameter,
    AuditExactStationDisposition2 disposition)
    : stock_state_digest_(std::move(stock_state_digest)),
      motion_digest_(std::move(motion_digest)),
      policy_digest_(std::move(policy_digest)),
      parameter_(std::move(parameter)),
      disposition_(disposition)
{
}

const AuditStockStateDigest2&
AuditExactStationWitness2::stock_state_digest() const noexcept
{
    return stock_state_digest_;
}

const NativeMotionDigest2&
AuditExactStationWitness2::motion_digest() const noexcept
{
    return motion_digest_;
}

const AuditPolicyDigest2&
AuditExactStationWitness2::policy_digest() const noexcept
{
    return policy_digest_;
}

AuditExactStationDisposition2
AuditExactStationWitness2::disposition() const noexcept
{
    return disposition_;
}

AuditDecisionWitness2 AuditDecisionEvidenceFactory2::exact_station_decision(
    const Stock2& stock,
    const AuditSegmentMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const AuditSegmentStationParameter2& parameter,
    AuditExactStationDisposition2 claimed_disposition,
    AuditDecisionAccounting2 accounting)
{
    const AuditStockStateIdentity2 stock_identity =
        AuditStockStateIdentity2::build(stock);
    AuditExactStationWitness2 witness = replay_station(
        stock_identity,
        motion.digest(),
        policy,
        parameter,
        claimed_disposition,
        stock);
    accounting.record_station_replay();
    return make_station_decision(
        stock_identity,
        motion.digest(),
        policy,
        limits,
        std::move(witness),
        accounting.snapshot());
}

AuditDecisionWitness2 AuditDecisionEvidenceFactory2::exact_station_decision(
    const Stock2& stock,
    const AuditCircleMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const AuditCircleStationParameter2& parameter,
    AuditExactStationDisposition2 claimed_disposition,
    AuditDecisionAccounting2 accounting)
{
    if (parameter.motion_digest().bytes() != motion.digest().bytes()) {
        throw AuditDecisionEvidenceReplayError(
            "circle station parameter belongs to a foreign motion");
    }
    const AuditStockStateIdentity2 stock_identity =
        AuditStockStateIdentity2::build(stock);
    AuditExactStationWitness2 witness = replay_station(
        stock_identity,
        motion.digest(),
        policy,
        parameter,
        claimed_disposition,
        stock);
    accounting.record_station_replay();
    return make_station_decision(
        stock_identity,
        motion.digest(),
        policy,
        limits,
        std::move(witness),
        accounting.snapshot());
}

AuditDecisionWitness2 AuditDecisionEvidenceFactory2::exact_station_decision(
    const Stock2& stock,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const AuditArcStationParameter2& parameter,
    AuditExactStationDisposition2 claimed_disposition)
{
    return exact_station_decision(
        stock,
        motion,
        policy,
        limits,
        parameter,
        claimed_disposition,
        AuditDecisionAccounting2{});
}

AuditDecisionWitness2 AuditDecisionEvidenceFactory2::exact_station_decision(
    const Stock2& stock,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& limits,
    const AuditArcStationParameter2& parameter,
    AuditExactStationDisposition2 claimed_disposition,
    AuditDecisionAccounting2 accounting)
{
    if (parameter.motion_digest().bytes() != motion.digest().bytes()) {
        throw AuditDecisionEvidenceReplayError(
            "arc station parameter belongs to a foreign motion");
    }
    const AuditStockStateIdentity2 stock_identity =
        AuditStockStateIdentity2::build(stock);
    AuditExactStationWitness2 witness = replay_station(
        stock_identity,
        motion.digest(),
        policy,
        parameter,
        claimed_disposition,
        stock);
    accounting.record_station_replay();
    return make_station_decision(
        stock_identity,
        motion.digest(),
        policy,
        limits,
        std::move(witness),
        accounting.snapshot());
}
