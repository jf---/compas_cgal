#include "audit_policy_2.h"

#include "canonical_encoding.h"

#include <CGAL/enum.h>
#include <CGAL/number_utils.h>

#include <bit>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <string>
#include <utility>

AuditCapObservation2 AuditCapObservation2::build(
    double authored_radians,
    double supplied_chord_ratio)
{
    if (!std::isfinite(authored_radians)
        || !std::isfinite(supplied_chord_ratio)) {
        throw AuditPolicyNonFiniteInputError(
            "audit policy cap inputs must be finite binary64 values");
    }
    if (!(authored_radians > 0.0
          && authored_radians <= std::numbers::pi)) {
        throw AuditPolicyEngagementCapRangeError(
            "audit policy cap radians must lie in (0, pi]");
    }
    const double recomputed = audit_cap_chord_ratio(authored_radians);
    if (std::bit_cast<std::uint64_t>(recomputed)
        != std::bit_cast<std::uint64_t>(supplied_chord_ratio)) {
        throw AuditPolicyCapSurrogateMismatchError(
            "audit policy chord surrogate is not the bit-identical native cap observation");
    }
    std::string canonical = canonical_encode_tagged_union(
        "audit-cap-observation-v1",
        canonical_encode_component_map({
            {"authored-radians", canonical_encode_binary64(authored_radians)},
            {"chord-ratio", canonical_encode_binary64(supplied_chord_ratio)},
        }));
    return AuditCapObservation2(
        Epeck::FT(authored_radians),
        Epeck::FT(supplied_chord_ratio),
        std::move(canonical));
}

double audit_cap_chord_ratio(double authored_radians)
{
    if (!std::isfinite(authored_radians)) {
        throw AuditPolicyNonFiniteInputError(
            "audit policy cap radians must be finite");
    }
    if (!(authored_radians > 0.0
          && authored_radians <= std::numbers::pi)) {
        throw AuditPolicyEngagementCapRangeError(
            "audit policy cap radians must lie in (0, pi]");
    }
    const double half_cap_sine = std::sin(0.5 * authored_radians);
    const double ratio = 4.0 * half_cap_sine * half_cap_sine;
    if (!(ratio > 0.0 && ratio <= 4.0)) {
        throw AuditPolicyCapSurrogateMismatchError(
            "audit policy cap has no representable chord ratio surrogate in (0, 4]");
    }
    return ratio;
}

AuditCapObservation2::AuditCapObservation2(
    Epeck::FT authored_radians,
    Epeck::FT chord_ratio,
    std::string canonical_bytes)
    : authored_radians_(std::move(authored_radians)),
      chord_ratio_(std::move(chord_ratio)),
      canonical_bytes_(std::move(canonical_bytes))
{
}

const Epeck::FT& AuditCapObservation2::authored_radians() const noexcept
{
    return authored_radians_;
}

const Epeck::FT& AuditCapObservation2::chord_ratio() const noexcept
{
    return chord_ratio_;
}

const std::string& AuditCapObservation2::canonical_bytes() const noexcept
{
    return canonical_bytes_;
}

AuditPolicy2 AuditPolicy2::build(
    const AuditCapObservation2& engagement_cap,
    const Epeck::FT& tool_radius_mm,
    const Epeck::FT& depletion_chord_bound_mm,
    std::size_t center_count_limit)
{
    if (CGAL::sign(tool_radius_mm) != CGAL::POSITIVE) {
        throw AuditPolicyToolRadiusError(
            "audit policy tool radius must be exact positive");
    }
    if (CGAL::sign(depletion_chord_bound_mm) != CGAL::POSITIVE
        || CGAL::compare(depletion_chord_bound_mm, tool_radius_mm)
            != CGAL::SMALLER) {
        throw AuditPolicyDepletionChordBoundError(
            "audit policy depletion chord bound must be exact positive and smaller than tool radius");
    }
    if (center_count_limit == 0) {
        throw AuditPolicyCenterCountLimitError(
            "audit policy center-count limit must be positive");
    }
    std::string canonical = canonical_encode_tagged_union(
        "audit-policy-v1",
        canonical_encode_component_map({
            {"center-count-limit", canonical_encode_integer(
                 ExactAlgebraicInteger1(center_count_limit))},
            {"depletion-chord-bound-mm", canonical_audit_rational_bytes(
                 depletion_chord_bound_mm)},
            {"engagement-cap", engagement_cap.canonical_bytes()},
            {"tool-radius-mm", canonical_audit_rational_bytes(tool_radius_mm)},
        }));
    AuditPolicyDigest2 digest =
        AuditPolicyDigestAuthority2::hash_canonical(canonical);
    return AuditPolicy2(
        engagement_cap,
        tool_radius_mm,
        depletion_chord_bound_mm,
        center_count_limit,
        std::move(canonical),
        std::move(digest));
}

AuditPolicy2::AuditPolicy2(
    AuditCapObservation2 engagement_cap,
    Epeck::FT tool_radius_mm,
    Epeck::FT depletion_chord_bound_mm,
    std::size_t center_count_limit,
    std::string canonical_bytes,
    AuditPolicyDigest2 digest)
    : engagement_cap_(std::move(engagement_cap)),
      tool_radius_mm_(std::move(tool_radius_mm)),
      depletion_chord_bound_mm_(std::move(depletion_chord_bound_mm)),
      center_count_limit_(center_count_limit),
      canonical_bytes_(std::move(canonical_bytes)),
      digest_(std::move(digest))
{
}

const AuditCapObservation2& AuditPolicy2::engagement_cap() const noexcept
{
    return engagement_cap_;
}

const Epeck::FT& AuditPolicy2::tool_radius_mm() const noexcept
{
    return tool_radius_mm_;
}

const Epeck::FT& AuditPolicy2::depletion_chord_bound_mm() const noexcept
{
    return depletion_chord_bound_mm_;
}

std::size_t AuditPolicy2::center_count_limit() const noexcept
{
    return center_count_limit_;
}

const std::string& AuditPolicy2::canonical_bytes() const noexcept
{
    return canonical_bytes_;
}

const AuditPolicyDigest2& AuditPolicy2::digest() const noexcept
{
    return digest_;
}
