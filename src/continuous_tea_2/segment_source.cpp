#include "segment_source.h"

#include "../canonical_encoding.h"
#include "event_certificate.h"
#include "sha256.h"

#include <bit>
#include <cmath>
#include <cstdint>
#include <utility>

#include <CGAL/CORE/BigRat.h>
#include <CGAL/Fraction_traits.h>
#include <CGAL/number_utils.h>

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

std::string binary64_bits_text(double value)
{
    return std::to_string(std::bit_cast<std::uint64_t>(value));
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
    return from_exact(
        ExactSegmentMotion2{
            EPoint(Epeck::FT(x0), Epeck::FT(y0)),
            EPoint(Epeck::FT(x1), Epeck::FT(y1)),
        },
        Epeck::FT(tool_radius),
        Epeck::FT(cap_chord_ratio));
}

SegmentEventSource2 SegmentEventSource2::from_exact(
    const ExactSegmentMotion2& motion,
    const Epeck::FT& tool_radius,
    const Epeck::FT& cap_chord_ratio)
{
    if (motion.start == motion.end) {
        throw ZeroLengthSegmentMotionError(
            "segment motion endpoints must differ");
    }
    if (CGAL::sign(tool_radius) != CGAL::POSITIVE) {
        throw NonPositiveToolRadiusError(
            "tool radius must be exact positive");
    }
    if (CGAL::sign(cap_chord_ratio) != CGAL::POSITIVE
        || CGAL::compare(cap_chord_ratio, Epeck::FT(4)) == CGAL::LARGER) {
        throw InvalidCapChordRatioError(
            "cap chord ratio must be exact in (0, 4]");
    }
    return SegmentEventSource2(
        lift_exact(motion.start.x()),
        lift_exact(motion.start.y()),
        lift_exact(motion.end.x()),
        lift_exact(motion.end.y()),
        lift_exact(tool_radius),
        lift_exact(cap_chord_ratio));
}

ExactBinary64Rational2
SegmentEventSource2::lift_binary64(double value)
{
    const Rational exact = exact_binary64(value);
    return ExactBinary64Rational2(
        CORE::numerator(exact).convert_to<std::string>(),
        CORE::denominator(exact).convert_to<std::string>());
}

ExactBinary64Rational2 SegmentEventSource2::lift_exact(
    const Epeck::FT& value)
{
    using Traits = CGAL::Fraction_traits<Epeck::FT>;
    typename Traits::Numerator_type numerator;
    typename Traits::Denominator_type denominator;
    typename Traits::Decompose()(value, numerator, denominator);
    return ExactBinary64Rational2(
        CORE::BigInt(numerator.exact()).convert_to<std::string>(),
        CORE::BigInt(denominator.exact()).convert_to<std::string>());
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

FullCircleEventSource2 FullCircleEventSource2::from_binary64(
    double center_x,
    double center_y,
    double phase_dx,
    double phase_dy,
    bool clockwise,
    double tool_radius,
    double cap_chord_ratio)
{
    if (!std::isfinite(center_x) || !std::isfinite(center_y)
        || !std::isfinite(phase_dx) || !std::isfinite(phase_dy)
        || !std::isfinite(tool_radius)
        || !std::isfinite(cap_chord_ratio)) {
        throw NonFiniteFullCircleInputError(
            "full-circle binary64 inputs must be finite");
    }
    if (phase_dx == 0.0 && phase_dy == 0.0) {
        throw ZeroFullCirclePhaseError(
            "full-circle phase vector must be nonzero");
    }
    if (!(tool_radius > 0.0)) {
        throw NonPositiveToolRadiusError(
            "full-circle tool radius must be positive");
    }
    if (!(cap_chord_ratio > 0.0 && cap_chord_ratio <= 4.0)) {
        throw InvalidCapChordRatioError(
            "full-circle cap chord ratio must be in (0, 4]");
    }
    std::string canonical = encode_string_sequence({
        "full-circle-event-source-binary64-v1",
        canonical_encode_binary64(center_x),
        canonical_encode_binary64(center_y),
        canonical_encode_binary64(phase_dx),
        canonical_encode_binary64(phase_dy),
        clockwise ? "clockwise" : "counterclockwise",
        canonical_encode_binary64(tool_radius),
        canonical_encode_binary64(cap_chord_ratio),
    });
    return FullCircleEventSource2(
        ExactCircleMotion2{
            EPoint(Epeck::FT(center_x), Epeck::FT(center_y)),
            EVector(Epeck::FT(phase_dx), Epeck::FT(phase_dy)),
            clockwise,
        },
        Epeck::FT(tool_radius),
        Epeck::FT(cap_chord_ratio),
        std::move(canonical),
        encode_canonical_record(
            "full-circle-motion-binary64-v1",
            {
                binary64_bits_text(center_x),
                binary64_bits_text(center_y),
                binary64_bits_text(phase_dx),
                binary64_bits_text(phase_dy),
                clockwise ? "clockwise" : "counterclockwise",
            }),
        encode_canonical_record(
            "cap-chord-ratio-binary64-v1",
            {binary64_bits_text(cap_chord_ratio)}));
}

FullCircleEventSource2 FullCircleEventSource2::from_exact(
    const ExactCircleMotion2& motion,
    const Epeck::FT& tool_radius,
    const Epeck::FT& cap_chord_ratio)
{
    if (CGAL::sign(motion.phase_vector.squared_length()) != CGAL::POSITIVE) {
        throw ZeroFullCirclePhaseError(
            "full-circle phase vector must be exact nonzero");
    }
    if (CGAL::sign(tool_radius) != CGAL::POSITIVE) {
        throw NonPositiveToolRadiusError(
            "tool radius must be exact positive");
    }
    if (CGAL::sign(cap_chord_ratio) != CGAL::POSITIVE
        || CGAL::compare(cap_chord_ratio, Epeck::FT(4)) == CGAL::LARGER) {
        throw InvalidCapChordRatioError(
            "cap chord ratio must be exact in (0, 4]");
    }
    const auto rational_text = [](const Epeck::FT& value) {
        using Traits = CGAL::Fraction_traits<Epeck::FT>;
        typename Traits::Numerator_type numerator;
        typename Traits::Denominator_type denominator;
        typename Traits::Decompose()(value, numerator, denominator);
        const std::string n =
            CORE::BigInt(numerator.exact()).convert_to<std::string>();
        const std::string d =
            CORE::BigInt(denominator.exact()).convert_to<std::string>();
        return d == "1" ? n : n + "/" + d;
    };
    std::string canonical = encode_string_sequence({
        "full-circle-event-source-exact-v1",
        rational_text(motion.center.x()),
        rational_text(motion.center.y()),
        rational_text(motion.phase_vector.x()),
        rational_text(motion.phase_vector.y()),
        motion.clockwise ? "clockwise" : "counterclockwise",
        rational_text(tool_radius),
        rational_text(cap_chord_ratio),
    });
    return FullCircleEventSource2(
        motion,
        tool_radius,
        cap_chord_ratio,
        std::move(canonical),
        encode_canonical_record(
            "full-circle-motion-exact-v1",
            {
                rational_text(motion.center.x()),
                rational_text(motion.center.y()),
                rational_text(motion.phase_vector.x()),
                rational_text(motion.phase_vector.y()),
                motion.clockwise ? "clockwise" : "counterclockwise",
            }),
        encode_canonical_record(
            "cap-chord-ratio-exact-v1",
            {rational_text(cap_chord_ratio)}));
}

FullCircleEventSource2::FullCircleEventSource2(
    ExactCircleMotion2 motion,
    Epeck::FT tool_radius,
    Epeck::FT cap_chord_ratio,
    std::string canonical_bytes,
    std::string motion_identity_bytes,
    std::string cap_identity_bytes)
    : motion_(std::move(motion)),
      tool_radius_(std::move(tool_radius)),
      cap_chord_ratio_(std::move(cap_chord_ratio)),
      canonical_bytes_(std::move(canonical_bytes)),
      motion_identity_bytes_(std::move(motion_identity_bytes)),
      cap_identity_bytes_(std::move(cap_identity_bytes))
{
}

const ExactCircleMotion2& FullCircleEventSource2::motion() const noexcept
{
    return motion_;
}

const Epeck::FT& FullCircleEventSource2::tool_radius() const noexcept
{
    return tool_radius_;
}

const Epeck::FT& FullCircleEventSource2::cap_chord_ratio() const noexcept
{
    return cap_chord_ratio_;
}

const std::string& FullCircleEventSource2::canonical_bytes() const noexcept
{
    return canonical_bytes_;
}

const std::string&
FullCircleEventSource2::motion_identity_bytes() const noexcept
{
    return motion_identity_bytes_;
}

const std::string&
FullCircleEventSource2::cap_identity_bytes() const noexcept
{
    return cap_identity_bytes_;
}
