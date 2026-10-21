#pragma once

// The depth-2 quadratic tower of the capsule pieces (docs/held_coverage_certificate.md,
// "Predicate list", P1/P2).
//
// A capsule flank lies on a line n.x = k0 + k1 sqrt(L) with rational n, k0, k1, L
// (L = |q - p|^2 of the Omega edge it offsets). Its crossings with a rational
// circle have chart parameters m +/- n sqrt(d) with m, n, d in Q(sqrt(L)): a
// depth-2 tower. Its crossings with another flank line add one rational radical
// sqrt(L'). Every comparison the certificate needs on one piece reduces to
//
//     sign(t + u sqrt(e)),   t = m + n sqrt(d),   m, n, d, u, e in Q(sqrt(L)),
//
// which the documented mixed-radical rule decides with same-root `Sqrt_extension`
// calls only: sign(m), sign(n), compare(m^2, n^2 d), then the same three steps
// one level up. No cross-root arithmetic is ever formed (every operand of + - *
// is checked to share the base root, in release builds too), and no extension
// is built over a zero radicand: `make_one_root` folds a1 == 0 or root == 0 to a
// non-extended value first (exactness checklist question 17).

#include "coverage_certificate_2/kernel_2.h"

#include <optional>

namespace compas_cgal::coverage_certificate {

/// a0 + a1 sqrt(root), folded: non-extended when a1 == 0 or root == 0.
[[nodiscard]] OneRoot make_one_root(const FT& a0, const FT& a1, const FT& root);

/// a0 + a1 sqrt(root), folded further: also rational when root is the square of
/// a rational (decided exactly over the integers).
[[nodiscard]] OneRoot make_reduced_one_root(const FT& a0, const FT& a1, const FT& root);

/// The rational square root of q when q >= 0 is the square of a rational.
[[nodiscard]] std::optional<FT> rational_sqrt(const FT& q);

/// The exact sign `exact` that rule `rule` computed, checked against a certified
/// enclosure of the same value: when the enclosure excludes 0 its sign must be
/// the exact sign (review I1a). A wrong sign rule therefore fails at its first
/// well-separated comparison instead of steering an exact search into a hang.
///
/// Raises:
///     SignRuleContractError: the enclosure excludes 0 and disagrees.
CGAL::Sign require_sign_agrees(CGAL::Sign exact, const Interval& enclosure, const char* rule);

/// Exact sign of u + v sqrt(P) + w sqrt(Q) for rationals, P, Q >= 0: one positive
/// common scaling, then the mixed-radical rule over the integers.
[[nodiscard]] CGAL::Sign sign_rational_two(const FT& u, const FT& v, const FT& P, const FT& w, const FT& Q);

/// The root a value is extended over, 0 when it is rational.
[[nodiscard]] FT base_of(const OneRoot& value);

/// Raises CoverageCertificateNativeError when a and b are both extended over
/// different roots: the precondition of every same-root operation, enforced in
/// release builds.
void require_same_base(const OneRoot& a, const OneRoot& b);

/// Exact sign of a one-root value. A zero root or zero a1 is folded to the
/// rational sign before `Sqrt_extension::sign` is consulted.
[[nodiscard]] CGAL::Sign sign_one_root(const OneRoot& value);

/// v / q for rational q != 0, staying in v's extension.
[[nodiscard]] OneRoot divide(const OneRoot& value, const FT& q);

/// Interval enclosure of a one-root value from its exact parts (certified:
/// `Sqrt_extension::to_interval` rounds outward).
[[nodiscard]] Interval interval_of(const OneRoot& value);

/// A depth-2 tower value m + n sqrt(d), with m, n, d in one field Q(sqrt(L)) and
/// d >= 0 (checked where it is formed).
struct Tower {
    OneRoot m;
    OneRoot n;
    OneRoot d;

    /// The rational value q as a tower.
    [[nodiscard]] static Tower rational(const FT& q);
    /// The one-root value v as a tower (n = 0).
    [[nodiscard]] static Tower of(const OneRoot& v);
};

/// The common base root of t's parts (0 if all are rational).
///
/// Raises:
///     CoverageCertificateNativeError: two parts are extended over different roots.
[[nodiscard]] FT base_of(const Tower& t);

/// Exact sign of m + n sqrt(d).
///
/// Raises:
///     CoverageCertificateNativeError: d < 0, or mixed bases.
[[nodiscard]] CGAL::Sign sign_tower(const Tower& t);

/// Exact sign of t + u sqrt(e) with u, e in t's base field and e >= 0.
///
/// Raises:
///     CoverageCertificateNativeError: e < 0, or mixed bases.
[[nodiscard]] CGAL::Sign sign_tower_plus(const Tower& t, const OneRoot& u, const OneRoot& e);

/// Exact sign of t1 - t2 when the bases are compatible (equal, or one rational).
///
/// Raises:
///     CoverageCertificateNativeError: both bases are irrational and different
///         (the caller must use a provenance-based comparison instead).
[[nodiscard]] CGAL::Sign sign_difference(const Tower& t1, const Tower& t2);

/// True when sign_difference(t1, t2) is defined: equal bases, or one rational.
[[nodiscard]] bool compatible_bases(const Tower& t1, const Tower& t2);

/// A Circular_kernel_2 coordinate (Root_of_2, a rational-base one-root value)
/// as a tower, with a zero root folded away first.
[[nodiscard]] Tower tower_of_root(const CK::Root_of_2& value);

/// Certified interval enclosure of a tower value from its exact parts.
[[nodiscard]] Interval interval_of(const Tower& t);

/// Certified enclosure of a0 + a1 sqrt(r) without cancellation: when the terms
/// have opposite signs, (a0^2 - a1^2 r) / (a0 - a1 sqrt(r)) (exact identity; the
/// numerator is formed exactly, the denominator is a same-sign sum). Costs exact
/// products: for near-cancelling values only (near tangencies).
[[nodiscard]] Interval tight_interval(const OneRoot& value);
/// m + n sqrt(d) with tight enclosures of m, n and d.
[[nodiscard]] Interval tight_interval(const Tower& t);

/// Same-field tower arithmetic: both operands share the base root and the
/// radicand d (a tower with n == 0 joins any d).
///
/// Raises:
///     CoverageCertificateNativeError: different radicands or bases.
[[nodiscard]] Tower tower_add(const Tower& a, const Tower& b);
[[nodiscard]] Tower tower_sub(const Tower& a, const Tower& b);
[[nodiscard]] Tower tower_mul(const Tower& a, const Tower& b);
/// a * s for s in a's base field.
[[nodiscard]] Tower tower_scale(const Tower& a, const OneRoot& s);

/// Exact sign of P + Q sqrt(e) for towers P, Q of one field and rational e >= 0
/// (the mixed-radical rule one level above a tower: sign(P), sign(Q), then
/// compare(P^2, Q^2 e) as a tower sign). sqrt(e) may be algebraically related to
/// the field; the rule only uses that sqrt(e) is a non-negative real.
[[nodiscard]] CGAL::Sign sign_tower_root(const Tower& P, const Tower& Q, const FT& e);

/// The product (a0 + a1 sqrt(da)) (b0 + b1 sqrt(db)) of one-roots over two
/// possibly different radicands, as the tower (a0 b) + (a1 b) sqrt(da) over b's
/// field. The sum of such products with the same (da, db) is tower_add.
[[nodiscard]] Tower cross_product(const OneRoot& a, const OneRoot& b);

}  // namespace compas_cgal::coverage_certificate
