#pragma once

#include "exact/rational.h"

#include <CGAL/Sqrt_extension.h>
#include <CGAL/enum.h>
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

// Stage-0 status: nothing in the pipeline calls same_root_add or
// same_root_multiply yet, so their only caller today is exact_one_root_gate.
// (sign_mixed_radical, below, is on the deciding path.) The existing one-root
// site, engagement_2.cpp::as_radpoint, decomposes the arrangement's coordinates
// into rational parts and states the shared-root precondition in a comment
// rather than combining them through here.

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

/// Exact sign of the mixed two-radical form
/// `a + b*sqrt(alpha) + c*sqrt(beta) + d*sqrt(alpha*beta)`, with RATIONAL
/// a, b, c, d and RATIONAL alpha, beta >= 0.
///
/// The exact cap predicate reduces to this: orientation and squared-chord tests
/// of two cutter-circle points p, q whose coordinates live in Q(sqrt alpha) and
/// Q(sqrt beta) respectively expand to exactly this shape.
///
/// Idiom (docs/exactness.md "Numeric comparison is the exact-kernel idiom"):
/// compare at the NUMBER-TYPE level. OneRoot is RealEmbeddable, so CGAL::sign
/// and same-root CGAL::compare (and same-root +/-/*) are exact -- the derived
/// quantities are built INSIDE one extension Q(sqrt alpha) and CGAL decides,
/// rather than hand-rolling a bignum squaring routine. Cross-root OneRoot
/// arithmetic (documented UB) is never formed: sqrt(beta) only ever appears as
/// a squared factor `beta` (a rational), and the degenerate roots are folded
/// away with exact Rational compares before any extension is built.
///
/// Args:
///     a, b, c, d: the rational coefficients of the form.
///     alpha, beta: the rational radicands, required to be nonnegative. This is
///         a precondition, not a check: every caller in the deciding pipeline
///         takes them from `Sqrt_extension::root()`, which is nonnegative by
///         construction, and the binary64 ingress seam
///         (`audit_sign_mixed_radical_exact`) rejects a negative radicand with
///         `InvalidMixedRadicalRootError` before reaching here.
///
/// Returns:
///     CGAL::NEGATIVE, CGAL::ZERO or CGAL::POSITIVE, decided exactly.
[[nodiscard]] CGAL::Sign sign_mixed_radical(const Rational& a, const Rational& b,
                                            const Rational& c, const Rational& d,
                                            const Rational& alpha,
                                            const Rational& beta);

}  // namespace compas_cgal::exact
