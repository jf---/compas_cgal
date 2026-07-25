#pragma once

#include "partition_certificate.h"

#include <string>
#include <vector>

class ExactBinary64Rational2 {
public:
    const std::string& numerator() const noexcept;
    const std::string& denominator() const noexcept;
    std::string text() const;

private:
    friend class SegmentEventSource2;

    ExactBinary64Rational2(
        std::string numerator,
        std::string denominator);

    std::string numerator_;
    std::string denominator_;
};

class SegmentEventSource2 {
public:
    static SegmentEventSource2 from_binary64(
        double x0,
        double y0,
        double x1,
        double y1,
        double tool_radius,
        double cap_chord_ratio);

    const ExactBinary64Rational2& x0() const noexcept;
    const ExactBinary64Rational2& y0() const noexcept;
    const ExactBinary64Rational2& x1() const noexcept;
    const ExactBinary64Rational2& y1() const noexcept;
    const ExactBinary64Rational2& tool_radius() const noexcept;
    const ExactBinary64Rational2& cap_chord_ratio() const noexcept;
    std::vector<std::string> motion_data() const;
    const std::string& canonical_bytes() const noexcept;
    const std::string& canonical_digest() const noexcept;

private:
    static ExactBinary64Rational2 lift_binary64(double value);

    SegmentEventSource2(
        ExactBinary64Rational2 x0,
        ExactBinary64Rational2 y0,
        ExactBinary64Rational2 x1,
        ExactBinary64Rational2 y1,
        ExactBinary64Rational2 tool_radius,
        ExactBinary64Rational2 cap_chord_ratio);

    ExactBinary64Rational2 x0_;
    ExactBinary64Rational2 y0_;
    ExactBinary64Rational2 x1_;
    ExactBinary64Rational2 y1_;
    ExactBinary64Rational2 tool_radius_;
    ExactBinary64Rational2 cap_chord_ratio_;
    std::string canonical_bytes_;
    std::string canonical_digest_;
};

class NonFiniteSegmentInputError : public EventSubstrateError {
public:
    using EventSubstrateError::EventSubstrateError;
};

class ZeroLengthSegmentMotionError : public EventSubstrateError {
public:
    using EventSubstrateError::EventSubstrateError;
};

class NonPositiveToolRadiusError : public EventSubstrateError {
public:
    using EventSubstrateError::EventSubstrateError;
};

class InvalidCapChordRatioError : public EventSubstrateError {
public:
    using EventSubstrateError::EventSubstrateError;
};
