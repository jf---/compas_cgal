#include "audit_certification_2.h"

#include "canonical_encoding.h"

#include <utility>

namespace {

constexpr std::size_t MAX_AUDIT_DECISION_DEPTH = 64;
constexpr std::size_t MAX_AUDIT_DECISION_NODES = 1U << 24U;

std::string exact_size_bytes(std::size_t value)
{
    return canonical_audit_rational_bytes(Epeck::FT(value));
}

} // namespace

AuditSquaredSpatialFloorMm2 AuditSquaredSpatialFloorMm2::build(
    const Epeck::FT& value)
{
    if (CGAL::sign(value) != CGAL::POSITIVE) {
        throw AuditSquaredSpatialFloorError(
            "audit squared spatial floor must be positive mm^2");
    }
    std::string canonical = canonical_encode_tagged_union(
        "audit-squared-spatial-floor-mm2-v1",
        canonical_encode_component_map({
            {"value", canonical_audit_rational_bytes(value)},
        }));
    return AuditSquaredSpatialFloorMm2(value, std::move(canonical));
}

AuditSquaredSpatialFloorMm2::AuditSquaredSpatialFloorMm2(
    Epeck::FT value,
    std::string canonical_bytes)
    : value_(std::move(value)),
      canonical_bytes_(std::move(canonical_bytes))
{
}

const Epeck::FT& AuditSquaredSpatialFloorMm2::value() const noexcept
{
    return value_;
}

const std::string&
AuditSquaredSpatialFloorMm2::canonical_bytes() const noexcept
{
    return canonical_bytes_;
}

AuditDecisionLimits2 AuditDecisionLimits2::build(
    const AuditSquaredSpatialFloorMm2& squared_spatial_floor_mm,
    std::size_t max_depth,
    std::size_t max_nodes)
{
    if (max_depth > MAX_AUDIT_DECISION_DEPTH) {
        throw AuditDecisionDepthLimitError(
            "audit decision depth exceeds the sealed finite maximum");
    }
    if (max_nodes == 0 || max_nodes > MAX_AUDIT_DECISION_NODES) {
        throw AuditDecisionNodeLimitError(
            "audit decision nodes must be in the sealed finite range");
    }
    std::string canonical = canonical_encode_tagged_union(
        "audit-decision-limits-v1",
        canonical_encode_component_map({
            {"max-depth", exact_size_bytes(max_depth)},
            {"max-nodes", exact_size_bytes(max_nodes)},
            {"squared-spatial-floor-mm2",
             squared_spatial_floor_mm.canonical_bytes()},
        }));
    return AuditDecisionLimits2(
        squared_spatial_floor_mm,
        max_depth,
        max_nodes,
        std::move(canonical));
}

AuditDecisionLimits2::AuditDecisionLimits2(
    AuditSquaredSpatialFloorMm2 squared_spatial_floor_mm,
    std::size_t max_depth,
    std::size_t max_nodes,
    std::string canonical_bytes)
    : squared_spatial_floor_mm_(std::move(squared_spatial_floor_mm)),
      max_depth_(max_depth),
      max_nodes_(max_nodes),
      canonical_bytes_(std::move(canonical_bytes))
{
}

const AuditSquaredSpatialFloorMm2&
AuditDecisionLimits2::squared_spatial_floor_mm() const noexcept
{
    return squared_spatial_floor_mm_;
}

std::size_t AuditDecisionLimits2::max_depth() const noexcept
{
    return max_depth_;
}

std::size_t AuditDecisionLimits2::max_nodes() const noexcept
{
    return max_nodes_;
}

const std::string& AuditDecisionLimits2::canonical_bytes() const noexcept
{
    return canonical_bytes_;
}
