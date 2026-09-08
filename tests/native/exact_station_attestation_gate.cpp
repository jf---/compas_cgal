#include "exact/canonical.h"
#include "exact/errors.h"
#include "exact/rational.h"

#include "continuous_tea_2/station_source.h"

#include <CGAL/CORE/BigRat.h>

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace exact = compas_cgal::exact;

namespace {

using Integer = CORE::BigInt;
using Rational = CORE::BigRat;

/// Raised when a gate check fails. Named so the gate obeys the same
/// named-exceptions-only rule as the module it exercises.
class GateCheckFailedError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// The project builds with CMAKE_BUILD_TYPE Release, which defines NDEBUG and
// erases <cassert>. A gate whose checks vanish reports a vacuous pass, so every
// decision here goes through a helper that is never compiled out.
void require(bool condition, const char* message)
{
    if (!condition) {
        throw GateCheckFailedError(message);
    }
}

/// Byte comparisons raise the module's own drift error rather than the generic
/// gate error. `AttestationByteDriftError` was declared at stage 0 with no
/// raiser, against the day a projection site existed; stage 2 builds that site,
/// and this is the contract test the design says must raise it.
void require_bytes(
    const std::string& produced,
    const std::string& frozen,
    const char* message)
{
    if (produced != frozen) {
        throw exact::AttestationByteDriftError(message);
    }
}

/// The frozen `exact-rational-v1` record for the rational -355/113.
///
/// FROZEN COMPATIBILITY CONTRACT. These 56 bytes are the on-disk shape every
/// stored full-circle replay digest was computed over: an 8-byte big-endian
/// field count, then per field an 8-byte big-endian length followed by its
/// payload. -355/113 is chosen because it is neither dyadic nor an integer, so
/// it exercises a real numerator and denominator, and because its sign lives in
/// the numerator, which is what canonicalisation must guarantee.
///
/// This is the station lane's THREE-field `encode_string_sequence` framing. A
/// second, unrelated framing shares the same tag at boundary_events.cpp:156
/// (one field, tag-plus-NUL prefix, via `tagged_record`); this literal does not
/// speak for that one.
///
/// The literal exists because the differential checks below CANNOT see a change
/// here: legacy and projected paths share the one `encode_string_sequence`, so a
/// single edit to that framing moves both sides identically and leaves the
/// differential green while every stored digest silently stops verifying.
///
/// If a change makes this check fail, that change invalidates every stored
/// full-circle replay digest. It must be a deliberate, versioned decision: bump
/// the "exact-rational-v1" tag and migrate the stored digests. Editing this
/// literal to restore green is the one repair that is always wrong.
constexpr char kFrozenStationValueBytes[] =
    "\x00\x00\x00\x00\x00\x00\x00\x03"                          // 3 fields
    "\x00\x00\x00\x00\x00\x00\x00\x11" "exact-rational-v1"      // 17 bytes
    "\x00\x00\x00\x00\x00\x00\x00\x04" "-355"                   // 4 bytes
    "\x00\x00\x00\x00\x00\x00\x00\x03" "113";                   // 3 bytes

static_assert(
    sizeof(kFrozenStationValueBytes) - 1 == 56,
    "the frozen value literal must be 56 bytes: 8 + (8+17) + (8+4) + (8+3)");

/// An independent re-implementation of the canonical length framing.
///
/// Deliberately NOT `encode_string_sequence`: the outer
/// `station-event-source-v1` record nests four inner records and would be
/// unreadable as a raw literal, so it is anchored structurally instead. The
/// equality check in `frozen_framing_agrees_with_the_literal` ties this helper
/// to the raw literal above, so together they pin the framing bytes, the field
/// order and the field count without either standing alone.
std::string frozen_u64(std::size_t value)
{
    std::string result(8, '\0');
    for (std::size_t index = 0; index < 8; ++index) {
        const int shift = static_cast<int>(56 - 8 * index);
        result[index] = static_cast<char>(
            (static_cast<std::uint64_t>(value) >> shift) & 0xffU);
    }
    return result;
}

std::string frozen_sequence(const std::vector<std::string>& fields)
{
    std::string result = frozen_u64(fields.size());
    for (const std::string& field : fields) {
        result += frozen_u64(field.size());
        result += field;
    }
    return result;
}

std::string frozen_value(
    const std::string& numerator,
    const std::string& denominator)
{
    return frozen_sequence({"exact-rational-v1", numerator, denominator});
}

/// The probe station: centre (-355/113, 22/7), tool radius 1/2, cap 7/2.
///
/// Every component is a reduced non-dyadic fraction where one exists, the
/// radius is positive and the cap lies in (0, 4], so the source's own
/// validation admits it. TASK 2 REWRITES ONLY THIS FUNCTION BODY, to build the
/// same station through the exact API. The frozen expectations below never move.
StationEventSource2 probe_station()
{
    return StationEventSource2::build("-355/113", "22/7", "1/2", "7/2");
}

std::string frozen_station_bytes()
{
    return frozen_sequence({
        "station-event-source-v1",
        frozen_value("-355", "113"),
        frozen_value("22", "7"),
        frozen_value("1", "2"),
        frozen_value("7", "2"),
    });
}

/// Render a CORE::BigRat the way the legacy decoder expects to read it back.
std::string legacy_text(const Rational& value)
{
    const Integer numerator = CORE::numerator(value);
    const Integer denominator = CORE::denominator(value);
    return denominator == 1
        ? numerator.convert_to<std::string>()
        : numerator.convert_to<std::string>()
            + "/"
            + denominator.convert_to<std::string>();
}

/// Tie the structural helper to the raw literal, so neither anchor stands alone.
void frozen_framing_agrees_with_the_literal()
{
    const std::string literal(
        kFrozenStationValueBytes, sizeof(kFrozenStationValueBytes) - 1);
    require_bytes(
        frozen_value("-355", "113"),
        literal,
        "the gate's independent framing helper disagrees with the frozen "
        "literal; one of the two has been edited");
}

/// The `exact-rational-v1` record of a single attested value must not move.
void value_record_matches_the_frozen_literal()
{
    const std::string literal(
        kFrozenStationValueBytes, sizeof(kFrozenStationValueBytes) - 1);
    const ExactRational2 attested = ExactRational2::build("-355/113");
    require_bytes(
        attested.canonical_bytes(),
        literal,
        "exact-rational-v1 bytes differ from the frozen literal: every stored "
        "full-circle replay digest is invalidated");
    require(
        attested.numerator() == "-355",
        "attested numerator is not -355");
    require(
        attested.denominator() == "113",
        "attested denominator is not 113");
    require(
        attested.text() == "-355/113",
        "attested text is not -355/113");
}

/// The outer `station-event-source-v1` record must not move either: its tag,
/// its field count and the order of its four values are all frozen.
void station_record_matches_the_frozen_framing()
{
    const std::string frozen = frozen_station_bytes();
    require(
        frozen.size() == 281,
        "the frozen station record must be 281 bytes: "
        "8 + (8+23) + (8+56) + (8+52) + (8+51) + (8+51)");
    require_bytes(
        probe_station().canonical_bytes(),
        frozen,
        "station-event-source-v1 bytes differ from the frozen framing: every "
        "stored full-circle replay digest is invalidated");
}

/// The projection door must agree with the legacy text decoder on values the
/// station lane actually carries: arbitrary reduced rationals, not just the
/// dyadic ones a binary64 produces. `exact_canonical_gate` covers the dyadic
/// case; this covers the case that gate cannot reach.
void projection_matches_the_legacy_decoder_on_generic_rationals()
{
    const std::vector<std::pair<const char*, const char*>> fractions = {
        {"1", "3"},
        {"-355", "113"},
        {"1000000000000000000000000000001", "3"},
        {"-987654321098765432109876543210987",
         "1000000000000000000000000000000007"},
        {"246", "369"},   // a common factor that must be divided out: 2/3
        {"5", "-8"},      // a sign that must migrate to the numerator
        {"0", "17"},      // zero, whose canonical denominator is 1
        {"7", "1"},       // an integer
    };
    for (const auto& [numerator_text, denominator_text] : fractions) {
        // Named locals, not `Rational(Integer(a), Integer(b))`: the inline form
        // is a function declaration, not a value (most vexing parse).
        const Integer numerator(numerator_text);
        const Integer denominator(denominator_text);
        const Rational value(numerator, denominator);
        const exact::CanonicalRational projected =
            exact::to_canonical(exact::Rational(value));
        const ExactRational2 legacy = ExactRational2::build(legacy_text(value));
        require(
            projected.numerator() == legacy.numerator(),
            "projected numerator differs from the legacy decoder");
        require(
            projected.denominator() == legacy.denominator(),
            "projected denominator differs from the legacy decoder");
        require(
            projected.text() == legacy.text(),
            "projected text differs from the legacy decoder");
    }
}

/// The producers hand the station a value built by a CHAIN of exact arithmetic,
/// not a leaf. The claim stage 2 rests on is that the lazy carrier's exact
/// representative is bit-identical to the eager one, so the same chain computed
/// both ways must canonicalise to the same bytes. Depth 16 is chosen because
/// that is where the measured lazy/eager divergence risk would be largest --
/// see docs/number_types.md, "What the filter is actually worth".
void lazy_and_eager_chains_canonicalise_identically()
{
    exact::Rational lazy(1);
    Rational eager(1);
    for (int step = 0; step < 16; ++step) {
        lazy = (lazy * exact::Rational(3) + exact::Rational(7))
            / (lazy + exact::Rational(5));
        eager = (eager * Rational(3) + Rational(7))
            / (eager + Rational(5));
    }
    const exact::CanonicalRational projected = exact::to_canonical(lazy);
    const ExactRational2 legacy = ExactRational2::build(legacy_text(eager));
    require(
        projected.numerator() == legacy.numerator(),
        "a depth-16 lazy chain canonicalises to a different numerator than the "
        "same chain computed eagerly");
    require(
        projected.denominator() == legacy.denominator(),
        "a depth-16 lazy chain canonicalises to a different denominator than "
        "the same chain computed eagerly");
}

/// The source's own range checks are part of its contract and must survive the
/// signature change unchanged.
void invalid_stations_are_rejected()
{
    bool raised = false;
    try {
        (void)StationEventSource2::build("0", "0", "0", "1");
    } catch (const InvalidStationSourceError&) {
        raised = true;
    }
    require(raised, "a zero tool radius was accepted");

    raised = false;
    try {
        (void)StationEventSource2::build("0", "0", "1", "0");
    } catch (const InvalidStationSourceError&) {
        raised = true;
    }
    require(raised, "a zero cap chord ratio was accepted");

    raised = false;
    try {
        (void)StationEventSource2::build("0", "0", "1", "17/4");
    } catch (const InvalidStationSourceError&) {
        raised = true;
    }
    require(raised, "a cap chord ratio above 4 was accepted");

    // The closed upper end is admitted, and so is a negative coordinate.
    (void)StationEventSource2::build("-1/3", "-7/2", "1", "4");
}

}  // namespace

int main()
{
    try {
        frozen_framing_agrees_with_the_literal();
        value_record_matches_the_frozen_literal();
        station_record_matches_the_frozen_framing();
        projection_matches_the_legacy_decoder_on_generic_rationals();
        lazy_and_eager_chains_canonicalise_identically();
        invalid_stations_are_rejected();
    } catch (const std::exception& error) {
        std::printf("exact_station_attestation_gate FAILED: %s\n", error.what());
        return 1;
    }
    std::printf("exact_station_attestation_gate OK\n");
    return 0;
}
