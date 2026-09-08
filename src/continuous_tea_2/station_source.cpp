#include "station_source.h"

#include "event_certificate.h"
#include "exact/canonical.h"

#include <string_view>
#include <utility>

#include <CGAL/CORE/BigRat.h>
#include <CGAL/number_utils.h>

namespace exact = compas_cgal::exact;

namespace {

using Integer = CORE::BigInt;
using Rational = CORE::BigRat;

// LEGACY DECODER, reachable only through ExactRational2::build, which nothing in
// the pipeline calls since stage 2. Retained so the attestation gate can keep
// comparing the projection against it. Removing both is stage 6.
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

StationAttestation2 project_attestation(
    const exact::Rational& center_x,
    const exact::Rational& center_y,
    const exact::Rational& tool_radius,
    const exact::Rational& cap_chord_ratio)
{
    const ExactRational2 attested_center_x =
        ExactRational2::project(center_x);
    const ExactRational2 attested_center_y =
        ExactRational2::project(center_y);
    const ExactRational2 attested_tool_radius =
        ExactRational2::project(tool_radius);
    const ExactRational2 attested_cap_chord_ratio =
        ExactRational2::project(cap_chord_ratio);
    return {
        attested_center_x,
        attested_center_y,
        attested_tool_radius,
        attested_cap_chord_ratio,
        encode_string_sequence(
            {
                "station-event-source-v1",
                attested_center_x.canonical_bytes(),
                attested_center_y.canonical_bytes(),
                attested_tool_radius.canonical_bytes(),
                attested_cap_chord_ratio.canonical_bytes(),
            }),
    };
}

} // namespace

ExactRational2 ExactRational2::project(
    const exact::Rational& value)
{
    const exact::CanonicalRational canonical =
        exact::to_canonical(value);
    return ExactRational2(
        canonical.numerator(),
        canonical.denominator());
}

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
    const exact::Rational& center_x,
    const exact::Rational& center_y,
    const exact::Rational& tool_radius,
    const exact::Rational& cap_chord_ratio)
{
    // Exact, filtered decisions on the carrier itself. Epeck::FT is
    // RealEmbeddable, so CGAL::sign and CGAL::compare decide from the lazy
    // interval whenever it separates and materialise the exact rational only
    // when it does not -- where the text carrier forced a decode and a full
    // bignum comparison every time.
    if (CGAL::sign(tool_radius) != CGAL::POSITIVE) {
        throw InvalidStationSourceError(
            "station tool radius must be positive");
    }
    if (CGAL::sign(cap_chord_ratio) != CGAL::POSITIVE
        || CGAL::compare(cap_chord_ratio, exact::Rational(4))
            == CGAL::LARGER) {
        throw InvalidStationSourceError(
            "station cap chord ratio must lie in (0, 4]");
    }
    return StationEventSource2(
        center_x,
        center_y,
        tool_radius,
        cap_chord_ratio);
}

StationEventSource2::StationEventSource2(
    exact::Rational center_x,
    exact::Rational center_y,
    exact::Rational tool_radius,
    exact::Rational cap_chord_ratio)
    : center_x_(std::move(center_x)),
      center_y_(std::move(center_y)),
      tool_radius_(std::move(tool_radius)),
      cap_chord_ratio_(std::move(cap_chord_ratio)),
      // Declared last, so initialised last: the four carriers above are already
      // live when the projection reads them.
      attestation_(
          project_attestation(
              center_x_,
              center_y_,
              tool_radius_,
              cap_chord_ratio_))
{
}

const exact::Rational&
StationEventSource2::center_x() const noexcept
{
    return center_x_;
}

const exact::Rational&
StationEventSource2::center_y() const noexcept
{
    return center_y_;
}

const exact::Rational&
StationEventSource2::tool_radius() const noexcept
{
    return tool_radius_;
}

const exact::Rational&
StationEventSource2::cap_chord_ratio() const noexcept
{
    return cap_chord_ratio_;
}

const StationAttestation2&
StationEventSource2::attestation() const noexcept
{
    return attestation_;
}

const std::string&
StationEventSource2::canonical_bytes() const noexcept
{
    return attestation_.canonical_bytes;
}
