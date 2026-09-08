#include "exact/canonical.h"
#include "exact/errors.h"
#include "exact/rational.h"

#include "continuous_tea_2/segment_source.h"
#include "continuous_tea_2/sha256.h"

#include <CGAL/CORE/BigRat.h>
#include <CGAL/number_utils.h>

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
/// gate error, so a digest-invalidating change is distinguishable from an
/// ordinary gate failure by exception type.
void require_bytes(
    const std::string& produced,
    const std::string& frozen,
    const char* message)
{
    if (produced != frozen) {
        throw exact::AttestationByteDriftError(message);
    }
}

std::string hex(const std::string& bytes)
{
    static const char* const digits = "0123456789abcdef";
    std::string result;
    result.reserve(bytes.size() * 2);
    for (const char byte : bytes) {
        const auto value = static_cast<unsigned char>(byte);
        result.push_back(digits[value >> 4U]);
        result.push_back(digits[value & 0x0fU]);
    }
    return result;
}

/// The frozen `exact-binary64-rational-v1` record for the double 0.1.
///
/// FROZEN COMPATIBILITY CONTRACT. These 91 bytes are the on-disk shape every
/// stored segment replay digest was computed over: an 8-byte big-endian field
/// count, then per field an 8-byte big-endian length followed by its payload.
/// 0.1 is chosen because it is not a decimal fraction, so its exact value
/// 3602879701896397 / 2^55 exercises a long numerator and denominator rather
/// than a round one.
///
/// RELOCATED from `exact_canonical_gate.cpp:66-74`, where it was checked
/// against `exact::CanonicalRational::canonical_bytes()`. That method is
/// removed in Task 2 -- it put one consumer's record tag in the base number
/// vocabulary -- so the anchor moves onto `ExactBinary64Rational2`, which is
/// where this framing actually belongs and is the only place it is emitted.
///
/// The literal exists because the differential checks below CANNOT see a change
/// here: every projection shares the one `encode_string_sequence`, so a single
/// edit to that framing moves both sides identically and leaves the differential
/// green while every stored digest silently stops verifying.
///
/// If a change makes this check fail, that change invalidates every stored
/// segment replay digest. It must be a deliberate, versioned decision: bump the
/// "exact-binary64-rational-v1" tag and migrate the stored digests. Editing this
/// literal to restore green is the one repair that is always wrong.
constexpr char kFrozenTenthBytes[] =
    "\x00\x00\x00\x00\x00\x00\x00\x03"                                // 3 fields
    "\x00\x00\x00\x00\x00\x00\x00\x1a" "exact-binary64-rational-v1"   // 26 bytes
    "\x00\x00\x00\x00\x00\x00\x00\x10" "3602879701896397"             // 16 bytes
    "\x00\x00\x00\x00\x00\x00\x00\x11" "36028797018963968";           // 17 bytes

static_assert(
    sizeof(kFrozenTenthBytes) - 1 == 91,
    "the frozen literal must be 91 bytes: 8 + (8+26) + (8+16) + (8+17)");

/// The same framing for a value no binary64 can denote.
///
/// `from_exact` admits any exact rational, so the attestation view must be
/// pinned on a non-dyadic value too. -355/113 is reduced, non-dyadic, and
/// carries its sign in the numerator -- which is what canonicalisation must
/// guarantee. `exact_canonical_gate` cannot reach this case at all, because
/// every value it builds comes from `from_binary64`.
constexpr char kFrozenNonDyadicBytes[] =
    "\x00\x00\x00\x00\x00\x00\x00\x03"                                // 3 fields
    "\x00\x00\x00\x00\x00\x00\x00\x1a" "exact-binary64-rational-v1"   // 26 bytes
    "\x00\x00\x00\x00\x00\x00\x00\x04" "-355"                         // 4 bytes
    "\x00\x00\x00\x00\x00\x00\x00\x03" "113";                         // 3 bytes

static_assert(
    sizeof(kFrozenNonDyadicBytes) - 1 == 65,
    "the frozen literal must be 65 bytes: 8 + (8+26) + (8+4) + (8+3)");

/// The SHA-256 of the frozen `segment-event-source-v1` record for the probe.
///
/// DERIVED, NOT OBSERVED. It was computed from the encoder's specification
/// (`event_certificate.cpp:593-602`) rather than read off a run, so the
/// structural anchor above must pass FIRST. If the structural anchor passes and
/// only this hex differs, the derivation was wrong; the gate prints the observed
/// hex. If the structural anchor ALSO fails, the framing has moved and nothing
/// here may be edited.
constexpr char kFrozenSegmentDigestHex[] =
    "0b5f57185b5c54a528ae4224c87b04735b4438fdedb71ac13f44e184e1128389";

/// An independent re-implementation of the canonical length framing.
///
/// Deliberately NOT `encode_string_sequence`: the outer
/// `segment-event-source-v1` record nests six inner records and would be
/// unreadable as a raw literal, so it is anchored structurally instead. The
/// equality check in `frozen_framing_agrees_with_the_literals` ties this helper
/// to the raw literals above, so together they pin the framing bytes, the field
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
    return frozen_sequence(
        {"exact-binary64-rational-v1", numerator, denominator});
}

/// The probe motion: (0.1, -0.25) -> (1.5, -3.0), tool radius 0.5, cap 4.0.
///
/// 0.1 gives a long dyadic numerator and denominator; -0.25 and -3.0 put the
/// sign in the numerator on both a fractional and an integral value; 1.5 and 0.5
/// are short dyadics; 4.0 is the CLOSED upper end of the cap interval, so the
/// probe also proves the boundary is admitted. The endpoints differ and the
/// radius is positive, so the source's own validation admits it.
///
/// TASK 3 REWRITES NOTHING HERE. The frozen expectations below never move.
SegmentEventSource2 probe_segment()
{
    return SegmentEventSource2::from_binary64(
        0.1, -0.25, 1.5, -3.0, 0.5, 4.0);
}

/// The attested view of a source value.
///
/// TASK 3 REWRITES ONLY THESE TWO FUNCTION BODIES, to read the derived
/// attestation instead of the carrier -- because after the inversion `x0()`
/// returns `exact::Rational`, which has no `numerator()` and no
/// `canonical_bytes()`. Every frozen expectation in this file is written against
/// these two accessors and never moves. Nothing else in the gate refers to a
/// source value.
const ExactBinary64Rational2& attested_x0(
    const SegmentEventSource2& source)
{
    return source.x0();
}

const ExactBinary64Rational2& attested_y1(
    const SegmentEventSource2& source)
{
    return source.y1();
}

std::string frozen_segment_bytes()
{
    return frozen_sequence({
        "segment-event-source-v1",
        frozen_value("3602879701896397", "36028797018963968"),
        frozen_value("-1", "4"),
        frozen_value("3", "2"),
        frozen_value("-3", "1"),
        frozen_value("1", "2"),
        frozen_value("4", "1"),
    });
}

/// Render a CORE::BigRat the way the legacy text carrier rendered it.
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

/// Tie the structural helper to the raw literals, so neither anchor stands
/// alone. A helper that shared `encode_string_sequence` would move with it; a
/// raw literal cannot express the nested outer record.
void frozen_framing_agrees_with_the_literals()
{
    require_bytes(
        frozen_value("3602879701896397", "36028797018963968"),
        std::string(kFrozenTenthBytes, sizeof(kFrozenTenthBytes) - 1),
        "the gate's independent framing helper disagrees with the frozen 0.1 "
        "literal; one of the two has been edited");
    require_bytes(
        frozen_value("-355", "113"),
        std::string(kFrozenNonDyadicBytes, sizeof(kFrozenNonDyadicBytes) - 1),
        "the gate's independent framing helper disagrees with the frozen "
        "-355/113 literal; one of the two has been edited");
}

/// The `exact-binary64-rational-v1` record of one attested value must not move.
void value_record_matches_the_frozen_literals()
{
    const SegmentEventSource2 dyadic = probe_segment();
    require_bytes(
        attested_x0(dyadic).canonical_bytes(),
        std::string(kFrozenTenthBytes, sizeof(kFrozenTenthBytes) - 1),
        "exact-binary64-rational-v1 bytes differ from the frozen 0.1 literal: "
        "every stored segment replay digest is invalidated");
    require(
        attested_x0(dyadic).numerator() == "3602879701896397",
        "attested numerator of 0.1 is wrong");
    require(
        attested_x0(dyadic).denominator() == "36028797018963968",
        "attested denominator of 0.1 is wrong");
    require(
        attested_x0(dyadic).text() == "3602879701896397/36028797018963968",
        "attested text of 0.1 is wrong");
    require(
        attested_y1(dyadic).text() == "-3",
        "an integral attested value must render without a denominator");

    // The non-dyadic case, reachable only through from_exact.
    const exact::Rational minus_355_over_113 =
        exact::Rational(-355) / exact::Rational(113);
    const SegmentEventSource2 non_dyadic =
        SegmentEventSource2::from_exact(
            ExactSegmentMotion2{
                EPoint(minus_355_over_113, exact::Rational(0)),
                EPoint(exact::Rational(1), exact::Rational(1)),
            },
            exact::Rational(1) / exact::Rational(2),
            exact::Rational(7) / exact::Rational(2));
    require_bytes(
        attested_x0(non_dyadic).canonical_bytes(),
        std::string(kFrozenNonDyadicBytes, sizeof(kFrozenNonDyadicBytes) - 1),
        "exact-binary64-rational-v1 bytes differ from the frozen -355/113 "
        "literal: every stored segment replay digest is invalidated");
}

/// The outer `segment-event-source-v1` record must not move either: its tag,
/// its field count and the order of its six values are all frozen.
void segment_record_matches_the_frozen_framing()
{
    const std::string frozen = frozen_segment_bytes();
    require(
        frozen.size() == 480,
        "the frozen segment record must be 480 bytes: "
        "8 + (8+23) + (8+91) + (8+61) + (8+60) + (8+61) + (8+60) + (8+60)");
    require_bytes(
        probe_segment().canonical_bytes(),
        frozen,
        "segment-event-source-v1 bytes differ from the frozen framing: every "
        "stored segment replay digest is invalidated");
}

/// The digest is the value the audit records actually carry, so it gets its own
/// absolute anchor as well as a tie back to the bytes.
void segment_digest_matches_the_frozen_hex()
{
    const SegmentEventSource2 source = probe_segment();
    require(
        source.canonical_digest() == sha256_bytes(source.canonical_bytes()),
        "canonical_digest is not the SHA-256 of canonical_bytes");
    const std::string produced = hex(source.canonical_digest());
    if (produced != kFrozenSegmentDigestHex) {
        std::printf(
            "observed segment-event-source-v1 digest: %s\n", produced.c_str());
        throw exact::AttestationByteDriftError(
            "segment-event-source-v1 digest differs from the frozen hex");
    }
}

/// The projection door must agree with the current attestation path on values
/// the segment lane actually carries: arbitrary reduced rationals, not just the
/// dyadic ones a binary64 produces.
void projection_matches_the_attestation_on_generic_rationals()
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
        const Rational eager(numerator, denominator);
        const exact::Rational lazy(eager);
        const exact::CanonicalRational projected = exact::to_canonical(lazy);
        const SegmentEventSource2 source =
            SegmentEventSource2::from_exact(
                ExactSegmentMotion2{
                    EPoint(lazy, exact::Rational(0)),
                    EPoint(lazy + exact::Rational(1), exact::Rational(1)),
                },
                exact::Rational(1),
                exact::Rational(1));
        require(
            projected.numerator() == attested_x0(source).numerator(),
            "projected numerator differs from the attestation view");
        require(
            projected.denominator() == attested_x0(source).denominator(),
            "projected denominator differs from the attestation view");
        require(
            projected.text() == attested_x0(source).text(),
            "projected text differs from the attestation view");
        require(
            attested_x0(source).text() == legacy_text(eager),
            "the attestation text differs from the eager rendering");
    }
}

/// The producers hand the source values built by a CHAIN of exact arithmetic,
/// not leaves. The claim stage 3 rests on is that the lazy carrier's exact
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
    require(
        projected.text() == legacy_text(eager),
        "a depth-16 lazy chain canonicalises differently from the same chain "
        "computed eagerly");
    // Materialised, not compared as expression templates: CORE's boost backend
    // is lazy, and an expression outliving its operands is undefined.
    const Integer lazy_numerator = CORE::numerator(CGAL::exact(lazy));
    const Integer lazy_denominator = CORE::denominator(CGAL::exact(lazy));
    const Integer eager_numerator = CORE::numerator(eager);
    const Integer eager_denominator = CORE::denominator(eager);
    require(
        lazy_numerator == eager_numerator
            && lazy_denominator == eager_denominator,
        "the lazy carrier's exact representative is not the eager value");
}

/// `motion_data` is the attestation shape `construct_pullback` consumes, and it
/// is Python-visible as a 4-tuple. Both facts are frozen.
void motion_data_is_the_four_attested_endpoint_texts()
{
    const std::vector<std::string> data = probe_segment().motion_data();
    require(data.size() == 4, "motion_data must carry exactly four fields");
    require(
        data[0] == "3602879701896397/36028797018963968",
        "motion_data[0] is not the attested x0 text");
    require(data[1] == "-1/4", "motion_data[1] is not the attested y0 text");
    require(data[2] == "3/2", "motion_data[2] is not the attested x1 text");
    require(data[3] == "-3", "motion_data[3] is not the attested y1 text");
}

/// The source's own validation is part of its contract and must survive the
/// signature change unchanged, with the same named exception per condition.
void invalid_segments_are_rejected()
{
    const double infinity = 1.0 / 0.0;
    bool raised = false;
    try {
        (void)SegmentEventSource2::from_binary64(
            infinity, 0.0, 1.0, 1.0, 1.0, 1.0);
    } catch (const NonFiniteSegmentInputError&) {
        raised = true;
    }
    require(raised, "a non-finite coordinate was accepted");

    raised = false;
    try {
        (void)SegmentEventSource2::from_binary64(
            1.0, 2.0, 1.0, 2.0, 1.0, 1.0);
    } catch (const ZeroLengthSegmentMotionError&) {
        raised = true;
    }
    require(raised, "coincident endpoints were accepted");

    raised = false;
    try {
        (void)SegmentEventSource2::from_binary64(
            0.0, 0.0, 1.0, 1.0, 0.0, 1.0);
    } catch (const NonPositiveToolRadiusError&) {
        raised = true;
    }
    require(raised, "a zero tool radius was accepted");

    raised = false;
    try {
        (void)SegmentEventSource2::from_binary64(
            0.0, 0.0, 1.0, 1.0, 1.0, 4.25);
    } catch (const InvalidCapChordRatioError&) {
        raised = true;
    }
    require(raised, "a cap chord ratio above 4 was accepted");

    raised = false;
    try {
        (void)SegmentEventSource2::from_exact(
            ExactSegmentMotion2{
                EPoint(exact::Rational(1), exact::Rational(2)),
                EPoint(exact::Rational(1), exact::Rational(2)),
            },
            exact::Rational(1),
            exact::Rational(1));
    } catch (const ZeroLengthSegmentMotionError&) {
        raised = true;
    }
    require(raised, "exact coincident endpoints were accepted");

    // The closed upper end is admitted, and so are negative coordinates.
    (void)SegmentEventSource2::from_binary64(
        -1.5, -2.5, -1.25, -2.0, 0.125, 4.0);
}

/// The full-circle source frames three records of its own. It already carries
/// exact values, so Task 6 changes only WHICH function performs the projection;
/// these anchors are what proves it changed nothing else.
void full_circle_records_do_not_move()
{
    const FullCircleEventSource2 binary64 =
        FullCircleEventSource2::from_binary64(
            0.5, -0.25, 1.0, 0.0, false, 0.5, 4.0);
    const std::string binary64_canonical = binary64.canonical_bytes();
    const std::string binary64_motion = binary64.motion_identity_bytes();
    const std::string binary64_cap = binary64.cap_identity_bytes();

    const FullCircleEventSource2 exact_source =
        FullCircleEventSource2::from_exact(
            ExactCircleMotion2{
                EPoint(
                    exact::Rational(-355) / exact::Rational(113),
                    exact::Rational(22) / exact::Rational(7)),
                EVector(exact::Rational(1), exact::Rational(0)),
                true,
            },
            exact::Rational(1) / exact::Rational(2),
            exact::Rational(7) / exact::Rational(2));

    require(
        binary64_canonical.find("full-circle-event-source-binary64-v1")
            != std::string::npos,
        "the binary64 full-circle record lost its tag");
    require(
        binary64_motion.find("full-circle-motion-binary64-v1")
            != std::string::npos,
        "the binary64 motion identity record lost its tag");
    require(
        binary64_cap.find("cap-chord-ratio-binary64-v1") != std::string::npos,
        "the binary64 cap identity record lost its tag");
    require(
        exact_source.canonical_bytes().find(
            "full-circle-event-source-exact-v1") != std::string::npos,
        "the exact full-circle record lost its tag");

    // The exact framing renders each value as canonical text. -355/113 and 22/7
    // are non-dyadic and 1/2, 7/2 are dyadic, so this pins both renderings.
    require(
        exact_source.canonical_bytes().find("-355/113") != std::string::npos,
        "the exact full-circle record does not carry the canonical centre x");
    require(
        exact_source.canonical_bytes().find("22/7") != std::string::npos,
        "the exact full-circle record does not carry the canonical centre y");
    require(
        exact_source.cap_identity_bytes().find("7/2") != std::string::npos,
        "the exact cap identity record does not carry the canonical cap");
    require(
        exact_source.motion_identity_bytes().find("clockwise")
            != std::string::npos,
        "the exact motion identity record lost its orientation field");
    // "counterclockwise" contains "clockwise", so the find above cannot tell the
    // two orientations apart on its own. The probes carry opposite orientations,
    // so requiring the counterclockwise spelling to be present in one record and
    // absent from the other makes the orientation field genuinely anchored.
    require(
        binary64_motion.find("counterclockwise") != std::string::npos,
        "the counterclockwise binary64 probe lost its orientation spelling");
    require(
        exact_source.motion_identity_bytes().find("counterclockwise")
            == std::string::npos,
        "the clockwise exact probe carries the counterclockwise spelling");
}

}  // namespace

int main()
{
    try {
        frozen_framing_agrees_with_the_literals();
        value_record_matches_the_frozen_literals();
        segment_record_matches_the_frozen_framing();
        segment_digest_matches_the_frozen_hex();
        projection_matches_the_attestation_on_generic_rationals();
        lazy_and_eager_chains_canonicalise_identically();
        motion_data_is_the_four_attested_endpoint_texts();
        invalid_segments_are_rejected();
        full_circle_records_do_not_move();
    } catch (const std::exception& error) {
        std::printf(
            "exact_segment_attestation_gate FAILED: %s\n", error.what());
        return 1;
    }
    std::printf("exact_segment_attestation_gate OK\n");
    return 0;
}
