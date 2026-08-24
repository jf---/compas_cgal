#pragma once

#include "canonical_encoding.h"
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

template <class Domain>
class AuditDigest2 {
public:
    static AuditDigest2 from_bytes(std::string bytes)
    {
        if (bytes.size() != 32) {
            throw AuditDigestSizeError("audit digest must contain exactly 32 bytes");
        }
        return AuditDigest2(std::move(bytes));
    }

    const std::string& bytes() const noexcept
    {
        return bytes_;
    }

private:
    explicit AuditDigest2(std::string bytes)
        : bytes_(std::move(bytes))
    {
    }

    std::string bytes_;
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
