#include "circle_vertex_source.h"

#include "event_certificate.h"
#include "partition_certificate.h"

#include <optional>
#include <sstream>
#include <string_view>
#include <utility>
#include <vector>

#include <CGAL/CORE/BigInt.h>
#include <CGAL/number_utils.h>

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
                != std::string::npos
            || Integer(text.substr(separator + 1))
                == 0) {
            throw EventPartitionVerificationError(
                std::string(role)
                + " is not an exact rational");
        }
        return Rational(
            Integer(text.substr(0, separator)),
            Integer(text.substr(separator + 1)));
    } catch (const EventSubstrateError&) {
        throw;
    } catch (const std::exception&) {
        throw EventPartitionVerificationError(
            std::string(role)
            + " is not an exact rational");
    }
}

std::string rational_text(const Rational& value)
{
    std::ostringstream stream;
    stream << value;
    return stream.str();
}

Rational exact_rational(const Epeck::FT& value)
{
    return parse_rational(
        [&value]() {
            std::ostringstream stream;
            stream << CGAL::exact(value);
            return stream.str();
        }(),
        "full-circle coordinate");
}

std::optional<Rational> rational_square_root(
    const Rational& value)
{
    if (value < 0) {
        return std::nullopt;
    }
    Integer numerator_root;
    Integer denominator_root;
    if (!CGAL::is_square(
            CORE::numerator(value),
            numerator_root)
        || !CGAL::is_square(
            CORE::denominator(value),
            denominator_root)) {
        return std::nullopt;
    }
    return Rational(
        numerator_root,
        denominator_root);
}

} // namespace

FullCircleCoordinate2::FullCircleCoordinate2(
    Rational base,
    Rational radical_coefficient,
    Rational radicand)
    : base_(std::move(base)),
      radical_coefficient_(
          std::move(radical_coefficient)),
      radicand_(std::move(radicand))
{
    if (radicand_ < 0) {
        throw EventPartitionVerificationError(
            "full-circle coordinate has a negative radicand");
    }
    if (radical_coefficient_ == 0) {
        radicand_ = 0;
        return;
    }
    if (radicand_ == 0) {
        throw EventPartitionVerificationError(
            "full-circle coordinate has a zero non-rational radicand");
    }
    const std::optional<Rational> root =
        rational_square_root(radicand_);
    if (root.has_value()) {
        base_ += radical_coefficient_ * *root;
        radical_coefficient_ = 0;
        radicand_ = 0;
    }
}

FullCircleCoordinate2
FullCircleCoordinate2::from_exact(
    const GpsPoint::CoordNT& coordinate)
{
    return FullCircleCoordinate2(
        exact_rational(coordinate.a0()),
        exact_rational(coordinate.a1()),
        exact_rational(coordinate.root()));
}

FullCircleCoordinate2 FullCircleCoordinate2::decode(
    const std::string& source)
{
    const std::vector<std::string> fields =
        decode_string_sequence(source);
    if (fields.size() != 4
        || fields[0]
            != "full-circle-one-root-coordinate-v1") {
        throw EventPartitionVerificationError(
            "full-circle coordinate source is malformed");
    }
    FullCircleCoordinate2 coordinate(
        parse_rational(fields[1],
                       "full-circle coordinate base"),
        parse_rational(
            fields[2],
            "full-circle coordinate radical coefficient"),
        parse_rational(
            fields[3],
            "full-circle coordinate radicand"));
    if (coordinate.canonical_source() != source) {
        throw EventPartitionVerificationError(
            "full-circle coordinate source is not canonical");
    }
    return coordinate;
}

std::string FullCircleCoordinate2::canonical_source() const
{
    return encode_string_sequence(
        {
            "full-circle-one-root-coordinate-v1",
            rational_text(base_),
            rational_text(
                radical_coefficient_),
            rational_text(radicand_),
        });
}
