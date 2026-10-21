#pragma once

// One-root numbers and points of the clearance oracle: the circle-segment
// traits' own types (CoordNT = CGAL::Sqrt_extension<FT, FT, Tag_true, ...>,
// Point_2 = _One_root_point_2) over Epeck's lazy, filtered FT, as the exact
// centre domain uses them (held_centre_domain_2_kernel.h).
//
// CGAL supplies the arithmetic and the comparisons: CGAL::sign and
// CGAL::compare on Sqrt_extension are exact and filtered, and compare is
// supported across different extensions. This file adds only what the
// circle-segment contract leaves to the caller:
//   * the same-extension precondition of Sqrt_extension arithmetic
//     (a.root()==0 or b.root()==0 or a.root()==b.root()) is ENFORCED in release
//     builds (CrossRootArithmeticError), never assumed;
//   * a vanished coefficient is folded back to a rational, and a constructed
//     radicand that is a rational square is folded at construction
//     (docs/exactness.md, checklist 17), so a value's extension is the one the
//     Python oracle (benchmarks/held_one_root.py) assigns it;
//   * the sign of a two-radical expression u + w*sqrt(beta), u, w in
//     Q(sqrt(alpha)), from sign(u), sign(w) and compare(u^2, beta*w^2): three
//     supported exact calls (CONVENTIONS.md, "Numeric comparison is the exact-kernel
//     idiom").

#include "held_centre_domain_2_kernel.h"

namespace compas_cgal::held_declared_clearance {

using Kernel = held_centre_domain::Kernel;
using FT = held_centre_domain::FT;
using RationalPoint = Kernel::Point_2;
using Vector = Kernel::Vector_2;
using OneRootPoint = held_centre_domain::OneRootPoint;
using Root = held_centre_domain::Traits::CoordNT;

[[nodiscard]] int sign(const FT& value);
[[nodiscard]] int sign(const Root& value);
/// Exact three-way comparison, also across extensions.
[[nodiscard]] int compare(const Root& left, const Root& right);

/// The rational `value` as a one-root number.
[[nodiscard]] Root rational(const FT& value);
/// a + b*sqrt(c), validated (c >= 0) and folded (b == 0, c == 0, or c a rational square).
///
/// Raises:
///     NegativeRadicandError: c < 0.
[[nodiscard]] Root make_root(const FT& a, const FT& b, const FT& c);
/// a + b*sqrt(c) folded only for b == 0 or c == 0: for values that only ever
/// decide (never reported), where the rational-square test would force an
/// exact evaluation for nothing. A square c is a valid, if redundant, extension.
///
/// Raises:
///     NegativeRadicandError: c < 0.
[[nodiscard]] Root make_deciding_root(const FT& a, const FT& b, const FT& c);
/// Whether the value is rational (folded: no extension).
[[nodiscard]] bool is_rational(const Root& value);

/// Same-extension arithmetic, the precondition enforced.
///
/// Raises:
///     CrossRootArithmeticError: the operands lie in different extensions.
[[nodiscard]] Root plus(const Root& left, const Root& right);
[[nodiscard]] Root minus(const Root& left, const Root& right);
[[nodiscard]] Root times(const Root& left, const Root& right);
[[nodiscard]] Root scaled(const FT& scale, const Root& value);
[[nodiscard]] Root shifted(const Root& value, const FT& shift);

/// Exact sign of u + w*sqrt(beta): u, w in one extension, rational beta >= 0.
[[nodiscard]] int sign_of_sum(const Root& u, const Root& w, const FT& beta);

/// A point from two coordinates of one extension.
///
/// Raises:
///     MixedRootPointError: the coordinates lie in different extensions.
[[nodiscard]] OneRootPoint make_point(const Root& x, const Root& y);
[[nodiscard]] OneRootPoint embed(const RationalPoint& point);
/// Whether both coordinates are rational.
[[nodiscard]] bool is_rational(const OneRootPoint& point);
/// The rational point of a point with rational coordinates (precondition: is_rational(point)).
[[nodiscard]] RationalPoint rational_point(const OneRootPoint& point);

/// Exact lexicographic (x, then y) comparison, across extensions.
[[nodiscard]] int compare_xy(const OneRootPoint& first, const OneRootPoint& second);

/// Exact sign of cross(first - centre, second - centre), across extensions.
[[nodiscard]] int orientation(const RationalPoint& centre, const OneRootPoint& first, const OneRootPoint& second);

/// Exact order of the directions of `first` and `second` from `centre`, by
/// swept angle from the direction of `start` in `sense` (+1 CCW, -1 CW): -1, 0
/// or 1. Only directions count: the points need not lie on one circle.
[[nodiscard]] int compare_on_arc(const RationalPoint& centre, int sense, const OneRootPoint& start,
                                 const OneRootPoint& first, const OneRootPoint& second);

/// compare_on_arc for rational directions (vectors from the centre), by kernel predicates.
[[nodiscard]] int compare_on_arc(int sense, const Vector& start, const Vector& first, const Vector& second);

}  // namespace compas_cgal::held_declared_clearance
