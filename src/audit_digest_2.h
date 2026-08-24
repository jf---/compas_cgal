#pragma once

#include "canonical_encoding.h"
#include "continuous_tea_2/sha256.h"
#include "exact_motion_2.h"

#include <CGAL/Fraction_traits.h>

#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

class AuditDigestSizeError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class AuditInputDigestAuthority2;
class AuditNativeStockDigestAuthority2;
class AuditNativeRequestDigestAuthority2;
class AuthenticatedOperationDigestAuthority2;
class AuditPolicyDigestAuthority2;
class NativeMotionDigestAuthority2;
class NativeDecisionDigestAuthority2;
class ExactDepletionTraceDigestAuthority2;
class DepletionWitnessDigestAuthority2;
class AuditStockStateDigestAuthority2;
class StockLineageDigestAuthority2;
class AuditResultDigestAuthority2;
class AuditArcMotion2;
class AuditClassificationFactory2;
class AuditPolicy2;
class AuditNativeStockIdentity2;
class AuditNativeRequestIdentity2;
class ExactArcDepletionTrace2;
class AuditDecisionWitness2;
class AuditDepletionWitness2;
class AuditDepletionWitnessFactory2;
class AuditStockStateIdentity2;
class AuditLineage2;
class AuditMotionResult2;

template <class Domain>
class AuditDigest2 {
public:
    const std::string& bytes() const noexcept
    {
        return bytes_;
    }

private:
    explicit AuditDigest2(std::string bytes)
        : bytes_(std::move(bytes))
    {
    }

    friend class AuditInputDigestAuthority2;
    friend class AuditNativeStockDigestAuthority2;
    friend class AuditNativeRequestDigestAuthority2;
    friend class AuthenticatedOperationDigestAuthority2;
    friend class AuditPolicyDigestAuthority2;
    friend class NativeMotionDigestAuthority2;
    friend class NativeDecisionDigestAuthority2;
    friend class ExactDepletionTraceDigestAuthority2;
    friend class DepletionWitnessDigestAuthority2;
    friend class AuditStockStateDigestAuthority2;
    friend class StockLineageDigestAuthority2;
    friend class AuditResultDigestAuthority2;

    std::string bytes_;
};

struct AuditInputDigestDomain;
struct AuditNativeStockDigestDomain;
struct AuditStockStateDigestDomain;
struct AuditNativeRequestDigestDomain;
struct AuthenticatedOperationDigestDomain;
struct AuditPolicyDigestDomain;
struct NativeMotionDigestDomain;
struct NativeDecisionDigestDomain;
struct ExactDepletionTraceDigestDomain;
struct DepletionWitnessDigestDomain;
struct StockLineageDigestDomain;
struct AuditResultDigestDomain;

using AuditInputDigest2 = AuditDigest2<AuditInputDigestDomain>;
using AuditNativeStockDigest2 = AuditDigest2<AuditNativeStockDigestDomain>;
using AuditStockStateDigest2 = AuditDigest2<AuditStockStateDigestDomain>;
using AuditNativeRequestDigest2 = AuditDigest2<AuditNativeRequestDigestDomain>;
using AuthenticatedOperationDigest2 = AuditDigest2<AuthenticatedOperationDigestDomain>;
using AuditPolicyDigest2 = AuditDigest2<AuditPolicyDigestDomain>;
using NativeMotionDigest2 = AuditDigest2<NativeMotionDigestDomain>;
using NativeDecisionDigest2 = AuditDigest2<NativeDecisionDigestDomain>;
using ExactDepletionTraceDigest2 = AuditDigest2<ExactDepletionTraceDigestDomain>;
using DepletionWitnessDigest2 = AuditDigest2<DepletionWitnessDigestDomain>;
using StockLineageDigest2 = AuditDigest2<StockLineageDigestDomain>;
using AuditResultDigest2 = AuditDigest2<AuditResultDigestDomain>;

inline void require_audit_digest_size(const std::string& bytes)
{
    if (bytes.size() != 32) {
        throw AuditDigestSizeError("audit digest must contain exactly 32 bytes");
    }
}

class AuditInputDigestAuthority2 {
public:
    static AuditInputDigest2 from_external_bytes(std::string bytes)
    {
        require_audit_digest_size(bytes);
        return AuditInputDigest2(std::move(bytes));
    }
};

class AuthenticatedOperationDigestAuthority2 {
public:
    static AuthenticatedOperationDigest2 from_external_bytes(std::string bytes)
    {
        require_audit_digest_size(bytes);
        return AuthenticatedOperationDigest2(std::move(bytes));
    }
};

class AuditNativeStockDigestAuthority2 {
    static AuditNativeStockDigest2 hash_canonical(std::string_view canonical)
    {
        return AuditNativeStockDigest2(sha256_bytes(std::string(canonical)));
    }
    friend class AuditNativeStockIdentity2;
};

class AuditNativeRequestDigestAuthority2 {
    static AuditNativeRequestDigest2 hash_canonical(std::string_view canonical)
    {
        return AuditNativeRequestDigest2(sha256_bytes(std::string(canonical)));
    }
    friend class AuditNativeRequestIdentity2;
};

class AuditPolicyDigestAuthority2 {
    static AuditPolicyDigest2 hash_canonical(std::string_view canonical)
    {
        return AuditPolicyDigest2(sha256_bytes(std::string(canonical)));
    }
    friend class AuditPolicy2;
};

class NativeMotionDigestAuthority2 {
    static NativeMotionDigest2 hash_canonical(std::string_view canonical)
    {
        return NativeMotionDigest2(sha256_bytes(std::string(canonical)));
    }
    friend class AuditArcMotion2;
    friend class AuditClassificationFactory2;
};

class NativeDecisionDigestAuthority2 {
    static NativeDecisionDigest2 hash_canonical(std::string_view canonical)
    {
        return NativeDecisionDigest2(sha256_bytes(std::string(canonical)));
    }
    friend class AuditDecisionWitness2;
};

class DepletionWitnessDigestAuthority2 {
    static DepletionWitnessDigest2 hash_canonical(std::string_view canonical)
    {
        return DepletionWitnessDigest2(sha256_bytes(std::string(canonical)));
    }
    friend class AuditDepletionWitness2;
    friend class AuditDepletionWitnessFactory2;
};

class ExactDepletionTraceDigestAuthority2 {
    static ExactDepletionTraceDigest2 hash_canonical(std::string_view canonical)
    {
        return ExactDepletionTraceDigest2(
            sha256_bytes(std::string(canonical)));
    }
    friend class ExactArcDepletionTrace2;
    friend class AuditDepletionWitness2;
    friend class AuditDepletionWitnessFactory2;
};

class AuditStockStateDigestAuthority2 {
    static AuditStockStateDigest2 hash_canonical(std::string_view canonical)
    {
        return AuditStockStateDigest2(sha256_bytes(std::string(canonical)));
    }
    friend class AuditStockStateIdentity2;
};

class StockLineageDigestAuthority2 {
    static StockLineageDigest2 hash_canonical(std::string_view canonical)
    {
        return StockLineageDigest2(sha256_bytes(std::string(canonical)));
    }
    friend class AuditLineage2;
};

class AuditResultDigestAuthority2 {
    static AuditResultDigest2 hash_canonical(std::string_view canonical)
    {
        return AuditResultDigest2(sha256_bytes(std::string(canonical)));
    }
    friend class AuditMotionResult2;
};

inline void append_audit_size(std::string& target, std::size_t size)
{
    static_assert(sizeof(std::size_t) <= sizeof(std::uint64_t));
    for (int shift = 56; shift >= 0; shift -= 8) {
        target.push_back(static_cast<char>(
            (static_cast<std::uint64_t>(size) >> shift) & 0xffU));
    }
}

inline void append_audit_bytes(std::string& target, std::string_view value)
{
    append_audit_size(target, value.size());
    target.append(value);
}

inline std::string canonical_audit_rational_bytes(const Epeck::FT& value)
{
    using Traits = CGAL::Fraction_traits<Epeck::FT>;
    typename Traits::Numerator_type numerator;
    typename Traits::Denominator_type denominator;
    typename Traits::Decompose()(value, numerator, denominator);
    const auto exact_numerator = numerator.exact();
    const auto exact_denominator = denominator.exact();
    return canonical_encode_rational(
        CORE::BigRat(
            CORE::BigInt(exact_numerator),
            CORE::BigInt(exact_denominator)));
}

inline std::string canonical_audit_binary64_bytes(double value)
{
    return canonical_encode_binary64(value);
}
