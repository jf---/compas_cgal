#include "exact/one_root.h"

#include "exact/errors.h"

#include <CGAL/number_utils.h>

namespace compas_cgal::exact {

namespace {

void require_same_root(const OneRoot& a, const OneRoot& b)
{
    // a1() and root() are ONLY defined on an EXTENDED Sqrt_extension. A rational
    // coordinate carries a0() alone, and reading the extension parts
    // unconditionally is undefined behaviour -- the same guard as
    // stock_2.cpp::add_coordinate and engagement_2.cpp::as_radpoint. A
    // non-extended operand is compatible with any root, so it short-circuits.
    //
    // This is the same condition CGAL states in Sqrt_extension::check_roots,
    // promoted from a CGAL_precondition (erased under NDEBUG, which this
    // project builds with) to a raise that survives a release build.
    if (!a.is_extended() || !b.is_extended()) {
        return;
    }
    if (a.root() != b.root()) {
        throw CrossRootExtensionError(
            "one-root arithmetic requires operands in the same extension");
    }
}

}  // namespace

OneRoot same_root_add(const OneRoot& a, const OneRoot& b)
{
    require_same_root(a, b);
    return a + b;
}

OneRoot same_root_multiply(const OneRoot& a, const OneRoot& b)
{
    require_same_root(a, b);
    return a * b;
}

CGAL::Sign sign_mixed_radical(
    const Rational& a,
    const Rational& b,
    const Rational& c,
    const Rational& d,
    const Rational& alpha,
    const Rational& beta)
{
    const bool alpha_ext = !CGAL::is_zero(alpha);
    const bool beta_ext = !CGAL::is_zero(beta);

    // Root degeneracies fold to a single same-root value -> one CGAL::sign.
    //
    // The folding is LOAD-BEARING, not stylistic: CGAL::Sqrt_extension::sign_()
    // (Sqrt_extension/Sqrt_extension_type.h:303) returns sign(a1_) whenever
    // sign(a0_) is ZERO, which is the sign of a1*sqrt(root) only for root > 0.
    // Tag_true makes an EXTENDED value with root 0 representable, and for such a
    // value the number IS a0, so the correct answer is ZERO. Constructing a
    // OneRoot over a zero root here would therefore return a nonzero sign for an
    // expression equal to zero. Fold first, with exact Rational compares.
    if (!alpha_ext && !beta_ext) {
        return CGAL::sign(a);                              // fully rational
    }
    if (!beta_ext) {
        return CGAL::sign(OneRoot(a, b, alpha));           // a + b*sqrt(alpha)
    }
    if (!alpha_ext) {
        return CGAL::sign(OneRoot(a, c, beta));            // a + c*sqrt(beta)
    }
    if (alpha == beta) {                                   // sqrt(alpha*beta) = alpha
        // (a + d*alpha) + (b + c)*sqrt(alpha)
        return CGAL::sign(OneRoot(a + d * alpha, b + c, alpha));
    }

    // General case: group over the shared root alpha into u, w in Q(sqrt alpha),
    // so the form is u + sqrt(beta)*w. Its sign follows from sign(u), sign(w),
    // and -- for opposite non-zero signs -- which magnitude dominates, decided
    // exactly by compare(u^2, beta*w^2) (all same-root, so beta enters as the
    // rational OneRoot(beta), never as sqrt(beta)).
    const OneRoot u(a, b, alpha);
    const OneRoot w(c, d, alpha);
    const CGAL::Sign u_sign = CGAL::sign(u);
    const CGAL::Sign w_sign = CGAL::sign(w);
    if (w_sign == CGAL::ZERO) {
        return u_sign;                    // sqrt(beta)*w = 0 -> u
    }
    if (u_sign == CGAL::ZERO) {
        return w_sign;                    // u = 0 -> sqrt(beta)*w, beta > 0
    }
    if (u_sign == w_sign) {
        return u_sign;                    // like signs add
    }
    switch (CGAL::compare(u * u, w * w * OneRoot(beta))) {
        case CGAL::LARGER:
            return u_sign;                // |u| dominates
        case CGAL::SMALLER:
            return w_sign;                // sqrt(beta)*|w| dominates
        default:
            return CGAL::ZERO;            // equal magnitude, opposite sign
    }
}

}  // namespace compas_cgal::exact
