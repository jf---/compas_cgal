#include "segment_source.h"

#include "event_certificate.h"
#include "sha256.h"

#include <bit>
#include <cmath>
#include <cstdint>
#include <utility>

#include <CGAL/CORE/BigRat.h>

namespace {

using Integer = CORE::BigInt;
using Rational = CORE::BigRat;

Rational exact_binary64(double value)
{
    const std::uint64_t bits =
        std::bit_cast<std::uint64_t>(value);
    const bool negative = (bits >> 63U) != 0;
    const std::uint64_t exponent_bits =
        (bits >> 52U) & 0x7ffU;
    const std::uint64_t fraction_bits =
        bits & ((std::uint64_t(1) << 52U) - 1U);
    if (exponent_bits == 0 && fraction_bits == 0) {
        return Rational(0);
    }

    Integer significand = exponent_bits == 0
        ? Integer(fraction_bits)
        : Integer(
              (std::uint64_t(1) << 52U)
              | fraction_bits);
    if (negative) {
        significand = -significand;
    }
    const int exponent = exponent_bits == 0
        ? -1074
        : static_cast<int>(exponent_bits) - 1023 - 52;
    if (exponent >= 0) {
        return Rational(significand << exponent);
    }
    return Rational(
        significand,
        Integer(1) << -exponent);
}

} // namespace

ExactBinary64Rational2::ExactBinary64Rational2(
    std::string numerator,
    std::string denominator)
    : numerator_(std::move(numerator)),
      denominator_(std::move(denominator))
{
}

const std::string&
ExactBinary64Rational2::numerator() const noexcept
{
    return numerator_;
}

const std::string&
ExactBinary64Rational2::denominator() const noexcept
{
    return denominator_;
}

std::string ExactBinary64Rational2::text() const
{
    return denominator_ == "1"
        ? numerator_
        : numerator_ + "/" + denominator_;
}

std::string ExactBinary64Rational2::canonical_bytes() const
{
    return encode_string_sequence(
        {
            "exact-binary64-rational-v1",
            numerator_,
            denominator_,
        });
}

SegmentEventSource2 SegmentEventSource2::from_binary64(
    double x0,
    double y0,
    double x1,
    double y1,
    double tool_radius,
    double cap_chord_ratio)
{
    if (!std::isfinite(x0) || !std::isfinite(y0)
        || !std::isfinite(x1) || !std::isfinite(y1)
        || !std::isfinite(tool_radius)
        || !std::isfinite(cap_chord_ratio)) {
        throw NonFiniteSegmentInputError(
            "segment source values must be finite");
    }
    if (x0 == x1 && y0 == y1) {
        throw ZeroLengthSegmentMotionError(
            "segment motion endpoints must differ");
    }
    if (tool_radius <= 0.0) {
        throw NonPositiveToolRadiusError(
            "tool radius must be positive");
    }
    if (cap_chord_ratio <= 0.0
        || cap_chord_ratio > 4.0) {
        throw InvalidCapChordRatioError(
            "cap chord ratio must be in (0, 4]");
    }
    return SegmentEventSource2(
        lift_binary64(x0),
        lift_binary64(y0),
        lift_binary64(x1),
        lift_binary64(y1),
        lift_binary64(tool_radius),
        lift_binary64(cap_chord_ratio));
}

ExactBinary64Rational2
SegmentEventSource2::lift_binary64(double value)
{
    const Rational exact = exact_binary64(value);
    return ExactBinary64Rational2(
        CORE::numerator(exact).convert_to<std::string>(),
        CORE::denominator(exact).convert_to<std::string>());
}

SegmentEventSource2::SegmentEventSource2(
    ExactBinary64Rational2 x0,
    ExactBinary64Rational2 y0,
    ExactBinary64Rational2 x1,
    ExactBinary64Rational2 y1,
    ExactBinary64Rational2 tool_radius,
    ExactBinary64Rational2 cap_chord_ratio)
    : x0_(std::move(x0)),
      y0_(std::move(y0)),
      x1_(std::move(x1)),
      y1_(std::move(y1)),
      tool_radius_(std::move(tool_radius)),
      cap_chord_ratio_(std::move(cap_chord_ratio))
{
    canonical_bytes_ = encode_string_sequence(
        {
            "segment-event-source-v1",
            x0_.canonical_bytes(),
            y0_.canonical_bytes(),
            x1_.canonical_bytes(),
            y1_.canonical_bytes(),
            tool_radius_.canonical_bytes(),
            cap_chord_ratio_.canonical_bytes(),
        });
    canonical_digest_ = sha256_bytes(canonical_bytes_);
}

const ExactBinary64Rational2&
SegmentEventSource2::x0() const noexcept
{
    return x0_;
}

const ExactBinary64Rational2&
SegmentEventSource2::y0() const noexcept
{
    return y0_;
}

const ExactBinary64Rational2&
SegmentEventSource2::x1() const noexcept
{
    return x1_;
}

const ExactBinary64Rational2&
SegmentEventSource2::y1() const noexcept
{
    return y1_;
}

const ExactBinary64Rational2&
SegmentEventSource2::tool_radius() const noexcept
{
    return tool_radius_;
}

const ExactBinary64Rational2&
SegmentEventSource2::cap_chord_ratio() const noexcept
{
    return cap_chord_ratio_;
}

std::vector<std::string> SegmentEventSource2::motion_data() const
{
    return {
        x0_.text(),
        y0_.text(),
        x1_.text(),
        y1_.text(),
    };
}

const std::string&
SegmentEventSource2::canonical_bytes() const noexcept
{
    return canonical_bytes_;
}

const std::string&
SegmentEventSource2::canonical_digest() const noexcept
{
    return canonical_digest_;
}
