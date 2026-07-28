#include "station_source.h"

#include "event_certificate.h"

#include <string_view>
#include <utility>

#include <CGAL/CORE/BigRat.h>

namespace {

using Integer = CORE::BigInt;
using Rational = CORE::BigRat;

Rational parse_rational(
    const std::string& text,
    std::string_view role)
{
    const std::size_t separator = text.find('/');
    try {
        if (separator == std::string::npos) {
            return Rational(Integer(text));
        }
        if (text.find('/', separator + 1)
            != std::string::npos) {
            throw InvalidStationSourceError(
                std::string(role)
                + " is not one exact rational");
        }
        const Integer numerator(
            text.substr(0, separator));
        const Integer denominator(
            text.substr(separator + 1));
        if (denominator == 0) {
            throw InvalidStationSourceError(
                std::string(role)
                + " has zero denominator");
        }
        return Rational(numerator, denominator);
    } catch (const EventSubstrateError&) {
        throw;
    } catch (const std::exception&) {
        throw InvalidStationSourceError(
            std::string(role)
            + " is not one exact rational");
    }
}

} // namespace

ExactRational2 ExactRational2::build(
    const std::string& text)
{
    const Rational value =
        parse_rational(text, "station value");
    return ExactRational2(
        CORE::numerator(value)
            .convert_to<std::string>(),
        CORE::denominator(value)
            .convert_to<std::string>());
}

ExactRational2::ExactRational2(
    std::string numerator,
    std::string denominator)
    : numerator_(std::move(numerator)),
      denominator_(std::move(denominator))
{
}

const std::string&
ExactRational2::numerator() const noexcept
{
    return numerator_;
}

const std::string&
ExactRational2::denominator() const noexcept
{
    return denominator_;
}

std::string ExactRational2::text() const
{
    return denominator_ == "1"
        ? numerator_
        : numerator_ + "/" + denominator_;
}

std::string ExactRational2::canonical_bytes() const
{
    return encode_string_sequence(
        {
            "exact-rational-v1",
            numerator_,
            denominator_,
        });
}

StationEventSource2 StationEventSource2::build(
    const std::string& center_x,
    const std::string& center_y,
    const std::string& tool_radius,
    const std::string& cap_chord_ratio)
{
    const ExactRational2 exact_center_x =
        ExactRational2::build(center_x);
    const ExactRational2 exact_center_y =
        ExactRational2::build(center_y);
    const ExactRational2 exact_tool_radius =
        ExactRational2::build(tool_radius);
    const ExactRational2 exact_cap_chord_ratio =
        ExactRational2::build(cap_chord_ratio);
    const Rational radius = parse_rational(
        exact_tool_radius.text(),
        "station tool radius");
    const Rational cap = parse_rational(
        exact_cap_chord_ratio.text(),
        "station cap chord ratio");
    if (radius <= 0) {
        throw InvalidStationSourceError(
            "station tool radius must be positive");
    }
    if (cap <= 0 || cap > 4) {
        throw InvalidStationSourceError(
            "station cap chord ratio must lie in (0, 4]");
    }
    return StationEventSource2(
        exact_center_x,
        exact_center_y,
        exact_tool_radius,
        exact_cap_chord_ratio);
}

StationEventSource2::StationEventSource2(
    ExactRational2 center_x,
    ExactRational2 center_y,
    ExactRational2 tool_radius,
    ExactRational2 cap_chord_ratio)
    : center_x_(std::move(center_x)),
      center_y_(std::move(center_y)),
      tool_radius_(std::move(tool_radius)),
      cap_chord_ratio_(std::move(cap_chord_ratio))
{
    canonical_bytes_ = encode_string_sequence(
        {
            "station-event-source-v1",
            center_x_.canonical_bytes(),
            center_y_.canonical_bytes(),
            tool_radius_.canonical_bytes(),
            cap_chord_ratio_.canonical_bytes(),
        });
}

const ExactRational2&
StationEventSource2::center_x() const noexcept
{
    return center_x_;
}

const ExactRational2&
StationEventSource2::center_y() const noexcept
{
    return center_y_;
}

const ExactRational2&
StationEventSource2::tool_radius() const noexcept
{
    return tool_radius_;
}

const ExactRational2&
StationEventSource2::cap_chord_ratio() const noexcept
{
    return cap_chord_ratio_;
}

const std::string&
StationEventSource2::canonical_bytes() const noexcept
{
    return canonical_bytes_;
}
