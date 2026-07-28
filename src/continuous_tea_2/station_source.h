#pragma once

#include "partition_certificate.h"

#include <string>

class ExactRational2 {
public:
    static ExactRational2 build(
        const std::string& text);

    const std::string& numerator() const noexcept;
    const std::string& denominator() const noexcept;
    std::string text() const;
    std::string canonical_bytes() const;

private:
    ExactRational2(
        std::string numerator,
        std::string denominator);

    std::string numerator_;
    std::string denominator_;
};

class StationEventSource2 {
public:
    static StationEventSource2 build(
        const std::string& center_x,
        const std::string& center_y,
        const std::string& tool_radius,
        const std::string& cap_chord_ratio);

    const ExactRational2& center_x() const noexcept;
    const ExactRational2& center_y() const noexcept;
    const ExactRational2& tool_radius() const noexcept;
    const ExactRational2& cap_chord_ratio() const noexcept;
    const std::string& canonical_bytes() const noexcept;

private:
    StationEventSource2(
        ExactRational2 center_x,
        ExactRational2 center_y,
        ExactRational2 tool_radius,
        ExactRational2 cap_chord_ratio);

    ExactRational2 center_x_;
    ExactRational2 center_y_;
    ExactRational2 tool_radius_;
    ExactRational2 cap_chord_ratio_;
    std::string canonical_bytes_;
};

class InvalidStationSourceError
    : public EventSubstrateError {
public:
    using EventSubstrateError::EventSubstrateError;
};
