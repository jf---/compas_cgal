#pragma once

// Number types and the single binary64 door of the coverage certificate.
//
// The certificate is CORE-free by design (docs/held_coverage_certificate.md,
// "Where it lives"): every quantity is built from exact rational inputs plus,
// where a piece needs one, an explicit square root. The circle-circle core runs
// on CGAL's 2D Circular Kernel over the project's exact rational
// (`CGAL::Exact_rational`, boost cpp_rational under CGAL_DISABLE_GMP), behind a
// certified interval front (lazy_crossings_2.h).

#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Sqrt_extension.h>

#include <cstdint>
#include <string>

namespace compas_cgal::coverage_certificate {

using FT = CGAL::Exact_rational;
using CK = CGAL::Circular_kernel_2<CGAL::Cartesian<FT>, CGAL::Algebraic_kernel_for_circles_2_2<FT>>;
using Point_2 = CK::Point_2;
using Circle_2 = CK::Circle_2;
using CircularPoint_2 = CK::Circular_arc_point_2;

/// Certified enclosure arithmetic. `Interval_nt<true>` rounds outward on every
/// operation, so an interval that excludes zero decides a sign; it is a filter,
/// never a tolerance.
using Interval = CGAL::Interval_nt<true>;

/// A one-root number a0 + a1*sqrt(root) with rational parts, used for the
/// VERTEX_DISK centre (a convex vertex of the exact centre domain Omega*).
using OneRoot = CGAL::Sqrt_extension<FT, FT, CGAL::Tag_true, CGAL::Tag_true>;

/// How every filtered decision is taken.
///
/// LAZY: interval filter first, exact predicate when the enclosure cannot
/// decide. EXACT_ONLY: every decision goes straight to the exact predicate.
/// The two modes must produce byte-identical certificates; the contract test
/// `test_forced_exact_mode_is_byte_identical_*` holds them to it.
enum class FilterMode : std::uint8_t { Lazy = 0, ExactOnly = 1 };

/// Convert a binary64 to its exact rational value; the certificate's only
/// double-to-exact door. A finite binary64 is a dyadic rational, so this is
/// exact and total: no parsing, no snapping, no tolerance.
///
/// Raises:
///     compas_cgal::exact::NonFiniteBinary64Error: NaN or infinite input.
[[nodiscard]] FT inject_binary64(double value);

/// Parse an exact rational from signed decimal numerator and positive decimal
/// denominator text (the Python `Fraction` wire format).
///
/// Raises:
///     InvalidCertificateInputError: malformed text or a non-positive denominator.
[[nodiscard]] FT rational_from_text(const std::string& numerator, const std::string& denominator);

/// Signed decimal numerator and positive decimal denominator of `value`.
struct RationalText {
    std::string numerator;
    std::string denominator;
};
[[nodiscard]] RationalText rational_to_text(const FT& value);

[[nodiscard]] inline Interval interval_of(const FT& value) { return Interval(CGAL::to_interval(value)); }

}  // namespace compas_cgal::coverage_certificate
