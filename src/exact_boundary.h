#pragma once

#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

// Boundary guards for the exact-kernel seam (docs/exactness.md, "The boundary
// doctrine"). Doubles cross into exact-land ONCE, at a declared API boundary, by
// exact injection -- which PRESUPPOSES they are rationals. NaN and +/-Inf are not
// rationals at all, and a non-positive radius is not a cutter.
//
// These are emphatically NOT epsilons smuggled into an exact pipeline: they run
// strictly BEFORE any Epeck::FT is constructed, and they ask a representability
// and physicality question about a raw double -- never a geometric one. Every
// geometric decision downstream remains an exact predicate on exact quantities.
//
// Shared by engagement_2.cpp and stock_2.cpp so the two halves of the seam
// cannot drift: before this header existed, `Stock2::subtract_disk` accepted an
// infinite radius and silently emptied the entire stock, while its two siblings
// raised on the same input.
//
// What these guards actually change. The exact kernel ALREADY refuses a
// non-rational: `Epeck::FT(nan)` throws "Cannot convert a non-finite number to an
// integer." That message is the exact number type's, not nanobind's -- verified on
// the `_sign_mixed_radical` binding, whose signature is a pure (double x 6) -> int
// with no integer conversion at the boundary, and which raises it anyway. So the
// refusal exists; it just fires late, from deep inside a construction, naming an
// integer the caller never mentioned. These guards move that refusal forward to
// the seam and give it the parameter's name. They add no rule the kernel did not
// already enforce -- which is precisely why they cannot weaken exactness.

namespace exact_boundary {

// Significant decimal digits that ROUND-TRIP a binary64: printing at this
// precision and reading back reproduces the identical double. Named rather than
// inlined as `17` because the number is a round-trip guarantee, not a formatting
// preference -- and taken from <limits> so it is derived rather than asserted.
constexpr int BINARY64_ROUND_TRIP_DIGITS = std::numeric_limits<double>::max_digits10;

// The offending value, formatted for the exception message. The caller is shown
// the EXACT double they passed (and "nan"/"inf" verbatim), never a rounded
// paraphrase of it.
inline std::string format_double(double value)
{
    std::ostringstream os;
    os.precision(BINARY64_ROUND_TRIP_DIGITS);
    os << value;
    return os.str();
}

inline void require_finite(double value, const char* name)
{
    if (!std::isfinite(value))
        throw std::invalid_argument(std::string(name) + " must be finite (got " + format_double(value) + ").");
}

// Spelled `!(r > 0.0)` rather than `r <= 0.0` so the rejection stays NaN-safe
// independently of the finiteness check, matching the ratio guards in
// engagement_2.cpp.
inline void require_positive_radius(double radius, const char* name)
{
    require_finite(radius, name);
    if (!(radius > 0.0))
        throw std::invalid_argument(std::string(name) + " must be strictly positive (got " + format_double(radius) + ").");
}

}  // namespace exact_boundary
