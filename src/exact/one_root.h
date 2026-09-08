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

// Stage-0 status: nothing in the pipeline calls the two functions below yet, so
// today's only caller is exact_one_root_gate. The existing one-root site,
// engagement_2.cpp::as_radpoint (engagement_2.cpp:84-97), decomposes the
// arrangement's coordinates into rational parts and states the shared-root
// precondition in a comment rather than combining them through here. Moving
// sign_mixed_radical into this module is stage 1.

/// Add two one-root numbers that share a root.
///
/// CGAL's precondition is that the operands share an extension, and it is
/// checked through `is_extended()` first because `a1()` and `root()` are
/// documented only for an extended value -- the guard the repo already applies
/// at stock_2.cpp:912-915. The two halves fail differently and neither is
/// assumed: every Sqrt_extension constructor initialises `root_`, so reading
/// `root()` off an unextended operand is a silently wrong answer rather than a
/// fault, while performing the arithmetic across distinct roots is the genuine
/// undefined behaviour.
///
/// Raises:
///     CrossRootExtensionError: if both operands are extended and their roots
///         differ. An operand that is not extended is compatible with any root
///         and is always permitted. Note that `Tag_true` admits an EXTENDED
///         value whose root is zero: it denotes a rational, but it is an
///         extended operand, so it is NOT exempt.
[[nodiscard]] OneRoot same_root_add(const OneRoot& a, const OneRoot& b);

/// Multiply two one-root numbers that share a root.
///
/// Raises:
///     CrossRootExtensionError: if both operands are extended and their roots
///         differ.
[[nodiscard]] OneRoot same_root_multiply(const OneRoot& a, const OneRoot& b);

}  // namespace compas_cgal::exact
