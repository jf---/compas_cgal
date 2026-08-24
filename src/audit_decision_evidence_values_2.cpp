#include "audit_certification_internal_2.h"

#include "canonical_encoding.h"

#include <algorithm>
#include <utility>

std::string canonical_audit_decision_counters_bytes(
    const AuditDecisionCounters2& counters)
{
    return canonical_encode_tagged_union(
        "audit-decision-counters-v1",
        canonical_encode_component_map({
            {"coverage-replays", canonical_audit_rational_bytes(
                 Epeck::FT(counters.exact_coverage_replays()))},
            {"deepest-level", canonical_audit_rational_bytes(
                 Epeck::FT(counters.deepest_level()))},
            {"station-replays", canonical_audit_rational_bytes(
                 Epeck::FT(counters.exact_station_replays()))},
            {"visited-nodes", canonical_audit_rational_bytes(
                 Epeck::FT(counters.visited_nodes()))},
        }));
}

void AuditDecisionAccounting2::record_station_replay() noexcept
{
    ++station_replays_;
}

void AuditDecisionAccounting2::record_coverage_replay() noexcept
{
    ++coverage_replays_;
}

void AuditDecisionAccounting2::record_refinement_node(
    std::size_t depth) noexcept
{
    ++visited_nodes_;
    deepest_level_ = std::max(deepest_level_, depth);
}

std::size_t AuditDecisionAccounting2::visited_nodes() const noexcept
{
    return visited_nodes_;
}

std::size_t AuditDecisionAccounting2::deepest_level() const noexcept
{
    return deepest_level_;
}

AuditDecisionCounters2 AuditDecisionAccounting2::snapshot() const
{
    return AuditDecisionCounters2::build(
        visited_nodes_, deepest_level_, station_replays_, coverage_replays_);
}

AuditDecisionCounters2 AuditDecisionCounters2::build(
    std::size_t visited_nodes,
    std::size_t deepest_level,
    std::size_t exact_station_replays,
    std::size_t exact_coverage_replays)
{
    return AuditDecisionCounters2(
        visited_nodes,
        deepest_level,
        exact_station_replays,
        exact_coverage_replays);
}

AuditDecisionCounters2::AuditDecisionCounters2(
    std::size_t visited_nodes,
    std::size_t deepest_level,
    std::size_t exact_station_replays,
    std::size_t exact_coverage_replays)
    : visited_nodes_(visited_nodes),
      deepest_level_(deepest_level),
      exact_station_replays_(exact_station_replays),
      exact_coverage_replays_(exact_coverage_replays)
{
}

std::size_t AuditDecisionCounters2::visited_nodes() const noexcept
{
    return visited_nodes_;
}

std::size_t AuditDecisionCounters2::deepest_level() const noexcept
{
    return deepest_level_;
}

std::size_t AuditDecisionCounters2::exact_station_replays() const noexcept
{
    return exact_station_replays_;
}

std::size_t AuditDecisionCounters2::exact_coverage_replays() const noexcept
{
    return exact_coverage_replays_;
}

AuditDecisionWitness2::AuditDecisionWitness2(
    AuditTeaVerdict2 verdict,
    std::optional<AuditExactStationWitness2> station,
    std::optional<AuditCertifiedCoverage2> coverage,
    std::optional<AuditUnresolvedEvidence2> unresolved,
    AuditDecisionCounters2 counters,
    std::string stock_digest,
    std::string motion_digest,
    std::string policy_digest,
    std::string limits_bytes,
    std::string canonical_bytes)
    : verdict_(verdict),
      station_(std::move(station)),
      coverage_(std::move(coverage)),
      unresolved_(std::move(unresolved)),
      counters_(std::move(counters)),
      stock_digest_(std::move(stock_digest)),
      motion_digest_(std::move(motion_digest)),
      policy_digest_(std::move(policy_digest)),
      limits_bytes_(std::move(limits_bytes)),
      canonical_bytes_(std::move(canonical_bytes)),
      digest_(NativeDecisionDigestAuthority2::hash_canonical(canonical_bytes_))
{
}

AuditTeaVerdict2 AuditDecisionWitness2::verdict() const noexcept
{
    return verdict_;
}

bool AuditDecisionWitness2::has_exact_station_witness() const noexcept
{
    return station_.has_value();
}

bool AuditDecisionWitness2::has_certified_coverage() const noexcept
{
    return coverage_.has_value();
}

bool AuditDecisionWitness2::has_unresolved_evidence() const noexcept
{
    return unresolved_.has_value();
}

const AuditExactStationWitness2&
AuditDecisionWitness2::exact_station_witness() const
{
    if (!station_) {
        throw AuditDecisionEvidenceReplayError(
            "decision has no exact station witness");
    }
    return *station_;
}

const AuditCertifiedCoverage2&
AuditDecisionWitness2::certified_coverage() const
{
    if (!coverage_) {
        throw AuditDecisionEvidenceReplayError(
            "decision has no certified coverage");
    }
    return *coverage_;
}

const AuditUnresolvedEvidence2&
AuditDecisionWitness2::unresolved_evidence() const
{
    if (!unresolved_) {
        throw AuditDecisionEvidenceReplayError(
            "decision has no unresolved evidence");
    }
    return *unresolved_;
}

const AuditDecisionCounters2&
AuditDecisionWitness2::counters() const noexcept
{
    return counters_;
}

const std::string& AuditDecisionWitness2::canonical_bytes() const noexcept
{
    return canonical_bytes_;
}

const NativeDecisionDigest2& AuditDecisionWitness2::digest() const noexcept
{
    return digest_;
}

void AuditDecisionEvidenceFactory2::validate_counters(
    const AuditDecisionLimits2& limits,
    const AuditDecisionCounters2& counters)
{
    if (counters.visited_nodes() > limits.max_nodes()
        || counters.deepest_level() > limits.max_depth()) {
        throw AuditDecisionCounterError(
            "decision counters exceed the sealed node or depth budget");
    }
}
