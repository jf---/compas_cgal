#include "exact/canonical.h"
#include "exact/rational.h"

#include "continuous_tea_2/segment_source.h"

#include <cstdio>
#include <exception>
#include <stdexcept>
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

/// Denominators are always positive and reduced, which is what makes the bytes
/// a function of the value rather than of the carrier.
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
}

}  // namespace

int main()
{
    try {
        to_canonical_matches_existing_projection();
        canonical_form_is_reduced_and_positive();
    } catch (const std::exception& error) {
        std::printf("exact_canonical_gate FAILED: %s\n", error.what());
        return 1;
    }
    std::printf("exact_canonical_gate OK\n");
    return 0;
}
