#include "exact/canonical.h"
#include "exact/rational.h"

#include "continuous_tea_2/segment_source.h"

#include <cstdio>
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

namespace exact = compas_cgal::exact;

namespace {

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

const std::vector<double>& edge_doubles()
{
    static const std::vector<double> values = {
        5e-324, 2.2250738585072014e-308, 1.0, 0.5,
        4503599627370496.0,      // 2^52
        9007199254740992.0,      // 2^53
        9007199254740994.0,      // 2^53 + 2
        0.1, 3.0, 1e308,
    };
    return values;
}

/// The frozen attestation bytes of `to_canonical(from_binary64(0.1))`.
///
/// FROZEN COMPATIBILITY CONTRACT. These 91 bytes are the on-disk shape every
/// stored replay digest was computed over: an 8-byte big-endian field count,
/// then per field an 8-byte big-endian length followed by its payload. 0.1 is
/// chosen because it is not a decimal fraction, so its exact value
/// 3602879701896397 / 2^55 exercises a long numerator and denominator rather
/// than a round one.
///
/// This literal exists because the differential checks below CANNOT see a
/// change here: both projections share the one `encode_string_sequence`, so a
/// single edit to that framing moves both sides identically and leaves the
/// differential comparison green while every stored digest silently stops
/// verifying. An absolute constant is the only reference that does not move
/// with the code it is checking.
///
/// If a change makes this check fail, that change invalidates every stored
/// replay digest. It must be a deliberate, versioned decision and never a
/// silent one: bump the "exact-binary64-rational-v1" tag and migrate the stored
/// digests. Editing this literal to restore green is the one repair that is
/// always wrong.
constexpr char kFrozenTenthBytes[] =
    "\x00\x00\x00\x00\x00\x00\x00\x03"                                // 3 fields
    "\x00\x00\x00\x00\x00\x00\x00\x1a" "exact-binary64-rational-v1"   // 26 bytes
    "\x00\x00\x00\x00\x00\x00\x00\x10" "3602879701896397"             // 16 bytes
    "\x00\x00\x00\x00\x00\x00\x00\x11" "36028797018963968";           // 17 bytes

static_assert(
    sizeof(kFrozenTenthBytes) - 1 == 91,
    "the frozen literal must be 91 bytes: 8 + (8+26) + (8+16) + (8+17)");

/// Anchor the frozen format to a constant rather than to another function.
void canonical_bytes_match_the_frozen_literal()
{
    const std::string frozen(kFrozenTenthBytes, sizeof(kFrozenTenthBytes) - 1);
    const std::string produced =
        exact::to_canonical(exact::from_binary64(0.1)).canonical_bytes();
    require(
        produced.size() == frozen.size(),
        "canonical bytes changed length: the frozen attestation framing moved, "
        "which invalidates every stored replay digest");
    require(
        produced == frozen,
        "canonical bytes differ from the frozen literal: every stored replay "
        "digest is invalidated");
}

/// The promoted projection must agree with the existing one, byte for byte.
void to_canonical_matches_existing_projection()
{
    for (const double value : edge_doubles()) {
        for (const double signed_value : {value, -value}) {
            const exact::CanonicalRational promoted =
                exact::to_canonical(exact::from_binary64(signed_value));

            // Existing path, reached through the public source API.
            const SegmentEventSource2 source =
                SegmentEventSource2::from_binary64(
                    signed_value, 0.0, signed_value + 1.0, 1.0, 1.0, 1.0);

            require(
                promoted.numerator() == source.x0().numerator(),
                "promoted numerator differs from the existing projection");
            require(
                promoted.denominator() == source.x0().denominator(),
                "promoted denominator differs from the existing projection");
            require(
                promoted.text() == source.x0().text(),
                "promoted text differs from the existing projection");
            require(
                promoted.canonical_bytes() == source.x0().canonical_bytes(),
                "promoted canonical bytes differ from the existing projection");
        }
    }
}

/// Check one reduction probe. The failure text names the whole expectation, so
/// a mismatch in either component reads as the reduction not having happened.
void require_reduces_to(
    const exact::Rational& value,
    const char* numerator,
    const char* denominator,
    const char* failure)
{
    const exact::CanonicalRational canonical = exact::to_canonical(value);
    require(canonical.numerator() == numerator, failure);
    require(canonical.denominator() == denominator, failure);
}

/// Denominators are always positive and fractions always reduced, which is what
/// makes the bytes a function of the value rather than of the carrier.
///
/// CGAL's Fraction_traits concept guarantees only that `value == num / den`. It
/// does NOT guarantee a reduced pair or a positive denominator: both come from
/// the backend rational's auto-normalisation, which is a property of the
/// arithmetic backend rather than of the contract. A backend swap that
/// normalised differently would move canonical bytes for every attested value,
/// so the probes below pin the normalisation instead of trusting it.
void canonical_form_is_reduced_and_positive()
{
    for (const double value : edge_doubles()) {
        const exact::CanonicalRational canonical =
            exact::to_canonical(exact::from_binary64(-value));
        require(
            !canonical.denominator().empty(),
            "canonical denominator is empty");
        require(
            canonical.denominator().front() != '-',
            "canonical denominator is negative");
    }

    // A common factor must be divided out, not carried.
    require_reduces_to(
        exact::Rational(2) / exact::Rational(4), "1", "2",
        "2/4 did not canonicalise to 1/2");
    // A negative denominator must migrate to the numerator.
    require_reduces_to(
        exact::Rational(1) / exact::Rational(-2), "-1", "2",
        "1/-2 did not canonicalise to -1/2");
    // Both at once, and the doubly-negative pair must come back positive.
    require_reduces_to(
        exact::Rational(-4) / exact::Rational(-6), "2", "3",
        "-4/-6 did not canonicalise to 2/3");
}

}  // namespace

int main()
{
    try {
        canonical_bytes_match_the_frozen_literal();
        to_canonical_matches_existing_projection();
        canonical_form_is_reduced_and_positive();
    } catch (const std::exception& error) {
        std::printf("exact_canonical_gate FAILED: %s\n", error.what());
        return 1;
    }
    std::printf("exact_canonical_gate OK\n");
    return 0;
}
