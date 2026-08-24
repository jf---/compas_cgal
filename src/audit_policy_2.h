#pragma once

#include "audit_digest_2.h"

#include <cstddef>
#include <stdexcept>
#include <string>

class AuditPolicyNonFiniteInputError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class AuditPolicyEngagementCapRangeError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class AuditPolicyCapSurrogateMismatchError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class AuditPolicyToolRadiusError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class AuditPolicyDepletionChordBoundError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class AuditPolicyCenterCountLimitError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

double audit_cap_chord_ratio(double authored_radians);

class AuditCapObservation2 {
public:
    static AuditCapObservation2 build(
        double authored_radians,
        double supplied_chord_ratio);

    const Epeck::FT& authored_radians() const noexcept;
    const Epeck::FT& chord_ratio() const noexcept;
    const std::string& canonical_bytes() const noexcept;

private:
    AuditCapObservation2(
        Epeck::FT authored_radians,
        Epeck::FT chord_ratio,
        std::string canonical_bytes);

    Epeck::FT authored_radians_;
    Epeck::FT chord_ratio_;
    std::string canonical_bytes_;
};

class AuditPolicy2 {
public:
    static AuditPolicy2 build(
        const AuditCapObservation2& engagement_cap,
        const Epeck::FT& tool_radius_mm,
        const Epeck::FT& depletion_chord_bound_mm,
        std::size_t center_count_limit);

    const AuditCapObservation2& engagement_cap() const noexcept;
    const Epeck::FT& tool_radius_mm() const noexcept;
    const Epeck::FT& depletion_chord_bound_mm() const noexcept;
    std::size_t center_count_limit() const noexcept;
    const std::string& canonical_bytes() const noexcept;
    const AuditPolicyDigest2& digest() const noexcept;

private:
    AuditPolicy2(
        AuditCapObservation2 engagement_cap,
        Epeck::FT tool_radius_mm,
        Epeck::FT depletion_chord_bound_mm,
        std::size_t center_count_limit,
        std::string canonical_bytes,
        AuditPolicyDigest2 digest);

    AuditCapObservation2 engagement_cap_;
    Epeck::FT tool_radius_mm_;
    Epeck::FT depletion_chord_bound_mm_;
    std::size_t center_count_limit_;
    std::string canonical_bytes_;
    AuditPolicyDigest2 digest_;
};
