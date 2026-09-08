#pragma once

#include "exact/rational.h"

#include <CGAL/Sqrt_extension.h>
#include <CGAL/tags.h>

namespace compas_cgal::exact {

/// One-root algebraic numbers, a0 + a1 * sqrt(root).
///
/// Deliberately a bare alias rather than a checking wrapper: this type IS
/// Gps_circle_segment_traits_2<Epeck>::Point_2::CoordNT, and wrapping it would
/// reintroduce a conversion boundary with the arrangement package. The
/// precondition lives in the free functions below instead.
///
/// The two Sqrt_extension tags are what make that identity hold and are not
/// decorative: CGAL spells the arrangement coordinate
/// `Sqrt_extension<NT, NT, Tag_true, Boolean_tag<Filter_>>`
/// (Arr_geometry_traits/Circle_segment_2.h:46) with `Filter_` defaulting to
/// true. Leaving the tags at their defaults would name a DIFFERENT,
/// non-interconvertible type, so an arrangement coordinate could not even be
/// passed here. `exact_one_root_gate` static_asserts the identity.
using OneRoot =
    CGAL::Sqrt_extension<Rational, Rational, CGAL::Tag_true, CGAL::Tag_true>;

/// Add two one-root numbers that share a root.
///
/// CGAL's precondition is that the operands share an extension. It is checked
/// via `is_extended()` first, because `root()` itself is only defined on an
/// extended value; violating either is undefined behaviour rather than a wrong
/// answer, so neither is assumed.
///
/// Raises:
///     CrossRootExtensionError: if the operands carry distinct non-zero roots.
[[nodiscard]] OneRoot same_root_add(const OneRoot& a, const OneRoot& b);

/// Multiply two one-root numbers that share a root.
///
/// Raises:
///     CrossRootExtensionError: if the operands carry distinct non-zero roots.
[[nodiscard]] OneRoot same_root_multiply(const OneRoot& a, const OneRoot& b);

}  // namespace compas_cgal::exact
