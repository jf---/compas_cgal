#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

namespace compas_cgal::exact {

/// The project's exact rational carrier.
///
/// This is deliberately the field type of Epeck rather than a fresh type: it is
/// already the coefficient type of Gps_circle_segment_traits_2<Epeck>'s CoordNT,
/// so one-root numbers built on it interoperate with the arrangement package
/// with no conversion, and every value carries Lazy_exact_nt's interval filter.
using Rational = CGAL::Exact_predicates_exact_constructions_kernel::FT;

/// Convert a binary64 to its exact rational value.
///
/// This is the ONLY double-to-exact entry point in the codebase. A binary64 is a
/// dyadic rational, so the conversion is exact and total on finite input: there
/// is no parsing, no tolerance and no snapping.
///
/// Args:
///     value: a finite binary64.
///
/// Returns:
///     The exact rational denoted by `value`. Negative zero maps to zero.
///
/// Raises:
///     NonFiniteBinary64Error: if `value` is NaN or infinite.
[[nodiscard]] Rational from_binary64(double value);

}  // namespace compas_cgal::exact
