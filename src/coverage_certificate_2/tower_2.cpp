#include "coverage_certificate_2/tower_2.h"

#include "coverage_certificate_2/errors_2.h"

#include <boost/multiprecision/cpp_int.hpp>

#include <string>

namespace compas_cgal::coverage_certificate {

OneRoot make_one_root(const FT& a0, const FT& a1, const FT& root)
{
    if (CGAL::sign(root) == CGAL::NEGATIVE) throw CoverageCertificateNativeError("one-root radicand must be >= 0");
    if (CGAL::sign(a1) == CGAL::ZERO || CGAL::sign(root) == CGAL::ZERO) return OneRoot(a0);
    return OneRoot(a0, a1, root);
}

std::optional<FT> rational_sqrt(const FT& q)
{
    using Traits = CGAL::Fraction_traits<FT>;
    using Integer = Traits::Numerator_type;
    if (CGAL::sign(q) == CGAL::NEGATIVE) return std::nullopt;
    Integer num;
    Integer den;
    Traits::Decompose()(q, num, den);
    const Integer sn = boost::multiprecision::sqrt(num);
    const Integer sd = boost::multiprecision::sqrt(den);
    if (sn * sn != num || sd * sd != den) return std::nullopt;
    return Traits::Compose()(sn, sd);
}

OneRoot make_reduced_one_root(const FT& a0, const FT& a1, const FT& root)
{
    if (CGAL::sign(root) == CGAL::NEGATIVE) throw CoverageCertificateNativeError("one-root radicand must be >= 0");
    if (CGAL::sign(a1) == CGAL::ZERO || CGAL::sign(root) == CGAL::ZERO) return OneRoot(a0);
    if (const auto s = rational_sqrt(root)) return OneRoot(a0 + a1 * *s);
    return OneRoot(a0, a1, root);
}

CGAL::Sign require_sign_agrees(const CGAL::Sign exact, const Interval& enclosure, const char* rule)
{
    const bool positive = enclosure.inf() > 0;
    const bool negative = enclosure.sup() < 0;
    if ((positive && exact != CGAL::POSITIVE) || (negative && exact != CGAL::NEGATIVE)) {
        throw SignRuleContractError(std::string("exact sign rule '") + rule + "' returned " +
                                    std::to_string(static_cast<int>(exact)) + " but its certified enclosure [" +
                                    std::to_string(enclosure.inf()) + ", " + std::to_string(enclosure.sup()) +
                                    "] excludes 0");
    }
    return exact;
}

namespace {
Interval sqrt_enclosure(const Interval& radicand)
{
    // A radicand known exactly to be >= 0 may have an enclosure dipping below 0.
    return CGAL::sqrt(Interval(radicand.inf() < 0 ? 0.0 : radicand.inf(), radicand.sup() < 0 ? 0.0 : radicand.sup()));
}
}  // namespace

FT base_of(const OneRoot& value) { return value.is_extended() ? value.root() : FT(0); }

void require_same_base(const OneRoot& a, const OneRoot& b)
{
    if (a.is_extended() && b.is_extended() && a.root() != b.root()) {
        throw CoverageCertificateNativeError("cross-root one-root arithmetic requested (different extensions)");
    }
}

CGAL::Sign sign_one_root(const OneRoot& value)
{
    if (!value.is_extended() || CGAL::sign(value.a1()) == CGAL::ZERO || CGAL::sign(value.root()) == CGAL::ZERO) {
        return CGAL::sign(value.a0());
    }
    return require_sign_agrees(CGAL::sign(value), interval_of(value), "one-root");
}

OneRoot divide(const OneRoot& value, const FT& q)
{
    if (CGAL::sign(q) == CGAL::ZERO) throw CoverageCertificateNativeError("one-root division by zero");
    if (!value.is_extended()) return OneRoot(value.a0() / q);
    return make_one_root(value.a0() / q, value.a1() / q, value.root());
}

Interval interval_of(const OneRoot& value)
{
    if (!value.is_extended()) return Interval(CGAL::to_interval(value.a0()));
    return Interval(CGAL::to_interval(value));
}

Tower Tower::rational(const FT& q) { return {OneRoot(q), OneRoot(FT(0)), OneRoot(FT(0))}; }

Tower tower_of_root(const CK::Root_of_2& value)
{
    if (!value.is_extended() || CGAL::sign(value.root()) == CGAL::ZERO) return Tower::rational(value.a0());
    return {OneRoot(value.a0()), OneRoot(value.a1()), OneRoot(value.root())};
}

Tower Tower::of(const OneRoot& v) { return {v, OneRoot(FT(0)), OneRoot(FT(0))}; }

FT base_of(const Tower& t)
{
    require_same_base(t.m, t.n);
    require_same_base(t.m, t.d);
    require_same_base(t.n, t.d);
    if (t.m.is_extended()) return t.m.root();
    if (t.n.is_extended()) return t.n.root();
    if (t.d.is_extended()) return t.d.root();
    return FT(0);
}

namespace {
CGAL::Sign sign_tower_unchecked(const Tower& t);
}  // namespace

CGAL::Sign sign_tower(const Tower& t)
{
    const CGAL::Sign exact = sign_tower_unchecked(t);
    return require_sign_agrees(exact, interval_of(t.m) + interval_of(t.n) * sqrt_enclosure(interval_of(t.d)), "tower");
}

namespace {
CGAL::Sign sign_tower_unchecked(const Tower& t)
{
    (void)base_of(t);
    const CGAL::Sign sd = sign_one_root(t.d);
    if (sd == CGAL::NEGATIVE) throw CoverageCertificateNativeError("tower radicand is negative");
    const CGAL::Sign sm = sign_one_root(t.m);
    if (sd == CGAL::ZERO) return sm;
    const CGAL::Sign sn = sign_one_root(t.n);
    if (sn == CGAL::ZERO) return sm;
    if (sm == CGAL::ZERO) return sn;
    if (sm == sn) return sm;
    // Opposite signs: the larger magnitude wins; m^2 - n^2 d is a same-root value.
    const CGAL::Sign c = sign_one_root(t.m * t.m - t.n * t.n * t.d);
    if (c == CGAL::POSITIVE) return sm;
    if (c == CGAL::NEGATIVE) return sn;
    return CGAL::ZERO;
}
}  // namespace

CGAL::Sign sign_tower_plus(const Tower& t, const OneRoot& u, const OneRoot& e)
{
    // One common base for t, u and e (a rational part joins any field).
    FT common = base_of(t);
    for (const OneRoot* v : {&u, &e}) {
        if (!v->is_extended() || CGAL::sign(v->root()) == CGAL::ZERO) continue;
        if (CGAL::sign(common) == CGAL::ZERO) {
            common = v->root();
        } else if (common != v->root()) {
            throw CoverageCertificateNativeError("cross-root one-root arithmetic requested (different extensions)");
        }
    }
    const CGAL::Sign se = sign_one_root(e);
    if (se == CGAL::NEGATIVE) throw CoverageCertificateNativeError("tower extra radicand is negative");
    const CGAL::Sign st = sign_tower(t);
    if (se == CGAL::ZERO) return st;
    const CGAL::Sign su = sign_one_root(u);
    if (su == CGAL::ZERO) return st;
    if (st == CGAL::ZERO) return su;
    if (st == su) return st;
    // Opposite signs: compare t^2 = (m^2 + n^2 d) + 2 m n sqrt(d) with u^2 e.
    const Tower square{t.m * t.m + t.n * t.n * t.d - u * u * e, FT(2) * t.m * t.n, t.d};
    const CGAL::Sign c = sign_tower(square);
    const CGAL::Sign exact = c == CGAL::POSITIVE ? st : (c == CGAL::NEGATIVE ? su : CGAL::ZERO);
    const Interval enclosure = interval_of(t.m) + interval_of(t.n) * sqrt_enclosure(interval_of(t.d)) +
                               interval_of(u) * sqrt_enclosure(interval_of(e));
    return require_sign_agrees(exact, enclosure, "tower-plus");
}

namespace {
using Traits = CGAL::Fraction_traits<FT>;
using Integer = Traits::Numerator_type;

CGAL::Sign sign_of(const Integer& v) { return v > 0 ? CGAL::POSITIVE : (v < 0 ? CGAL::NEGATIVE : CGAL::ZERO); }

/// sign(a + b sqrt(p)) over the integers, p >= 0.
CGAL::Sign sign_int_one(const Integer& a, const Integer& b, const Integer& p)
{
    const CGAL::Sign sa = sign_of(a);
    const CGAL::Sign sb = p > 0 ? sign_of(b) : CGAL::ZERO;
    if (sb == CGAL::ZERO) return sa;
    if (sa == CGAL::ZERO || sa == sb) return sb;
    const Integer lhs = a * a;
    const Integer rhs = b * b * p;
    if (lhs > rhs) return sa;
    if (lhs < rhs) return sb;
    return CGAL::ZERO;
}

/// sign(a + b sqrt(p) + c sqrt(q)) over the integers, p, q >= 0: the same
/// mixed-radical rule as sign_tower_plus, without rational normalisation.
CGAL::Sign sign_int_two(const Integer& a, const Integer& b, const Integer& p, const Integer& c, const Integer& q)
{
    const CGAL::Sign st = sign_int_one(a, b, p);
    const CGAL::Sign sc = q > 0 ? sign_of(c) : CGAL::ZERO;
    if (sc == CGAL::ZERO) return st;
    if (st == CGAL::ZERO || st == sc) return st == CGAL::ZERO ? sc : st;
    // (a + b sqrt(p))^2 - c^2 q = (a^2 + b^2 p - c^2 q) + 2 a b sqrt(p).
    const CGAL::Sign c2 = sign_int_one(a * a + b * b * p - c * c * q, Integer(2) * a * b, p);
    if (c2 == CGAL::POSITIVE) return st;
    if (c2 == CGAL::NEGATIVE) return sc;
    return CGAL::ZERO;
}

bool all_rational(const Tower& t) { return !t.m.is_extended() && !t.n.is_extended() && !t.d.is_extended(); }

/// sign(u + v sqrt(P) + w sqrt(Q)) for rationals: multiply through by the
/// positive common scale ud * vd * Pd * wd * Qd and decide over the integers.
CGAL::Sign sign_rational_two_unchecked(const FT& u, const FT& v, const FT& P, const FT& w, const FT& Q)
{
    Integer un, ud, vn, vd, pn, pd, wn, wd, qn, qd;
    Traits::Decompose()(u, un, ud);
    Traits::Decompose()(v, vn, vd);
    Traits::Decompose()(P, pn, pd);
    Traits::Decompose()(w, wn, wd);
    Traits::Decompose()(Q, qn, qd);
    // v sqrt(pn / pd) = v sqrt(pn pd) / pd; likewise for w.
    const Integer sv = vd * pd;
    const Integer sw = wd * qd;
    const Integer a = un * sv * sw;
    const Integer b = vn * ud * sw;
    const Integer c = wn * ud * sv;
    return sign_int_two(a, b, pn * pd, c, qn * qd);
}
}  // namespace

CGAL::Sign sign_rational_two(const FT& u, const FT& v, const FT& P, const FT& w, const FT& Q)
{
    if (CGAL::sign(P) == CGAL::NEGATIVE || CGAL::sign(Q) == CGAL::NEGATIVE) {
        throw CoverageCertificateNativeError("two-radical sign with a negative radicand");
    }
    const CGAL::Sign exact = sign_rational_two_unchecked(u, v, P, w, Q);
    const Interval enclosure = interval_of(u) + interval_of(v) * sqrt_enclosure(interval_of(P)) +
                               interval_of(w) * sqrt_enclosure(interval_of(Q));
    return require_sign_agrees(exact, enclosure, "rational two-radical");
}

bool compatible_bases(const Tower& t1, const Tower& t2)
{
    const FT b1 = base_of(t1);
    const FT b2 = base_of(t2);
    return CGAL::sign(b1) == CGAL::ZERO || CGAL::sign(b2) == CGAL::ZERO || b1 == b2;
}

CGAL::Sign sign_difference(const Tower& t1, const Tower& t2)
{
    if (all_rational(t1) && all_rational(t2)) {
        // (m1 - m2) + n1 sqrt(d1) - n2 sqrt(d2), all rational: integer fast path.
        if (CGAL::sign(t1.d.a0()) == CGAL::NEGATIVE || CGAL::sign(t2.d.a0()) == CGAL::NEGATIVE) {
            throw CoverageCertificateNativeError("tower radicand is negative");
        }
        return sign_rational_two(t1.m.a0() - t2.m.a0(), t1.n.a0(), t1.d.a0(), -t2.n.a0(), t2.d.a0());
    }
    const FT b1 = base_of(t1);
    const FT b2 = base_of(t2);
    if (CGAL::sign(b1) != CGAL::ZERO && CGAL::sign(b2) != CGAL::ZERO && b1 != b2) {
        throw CoverageCertificateNativeError("tower difference across two different irrational bases");
    }
    if (CGAL::sign(b1) == CGAL::ZERO && CGAL::sign(b2) != CGAL::ZERO) {
        // Evaluate in t2's field: t1 - t2 = -(t2 - t1).
        return CGAL::opposite(sign_tower_plus({t2.m - t1.m, t2.n, t2.d}, -t1.n, t1.d));
    }
    return sign_tower_plus({t1.m - t2.m, t1.n, t1.d}, -t2.n, t2.d);
}

namespace {
/// The common radicand of two towers (one with n == 0 joins the other).
OneRoot common_d(const Tower& a, const Tower& b)
{
    const bool a_has = sign_one_root(a.n) != CGAL::ZERO && sign_one_root(a.d) != CGAL::ZERO;
    const bool b_has = sign_one_root(b.n) != CGAL::ZERO && sign_one_root(b.d) != CGAL::ZERO;
    if (a_has && b_has) {
        require_same_base(a.d, b.d);
        if (!(a.d == b.d)) throw CoverageCertificateNativeError("tower arithmetic across two different radicands");
    }
    return a_has ? a.d : (b_has ? b.d : OneRoot(FT(0)));
}
}  // namespace

Tower tower_add(const Tower& a, const Tower& b)
{
    const OneRoot d = common_d(a, b);
    require_same_base(a.m, b.m);
    return {a.m + b.m, a.n + b.n, d};
}

Tower tower_sub(const Tower& a, const Tower& b)
{
    const OneRoot d = common_d(a, b);
    require_same_base(a.m, b.m);
    return {a.m - b.m, a.n - b.n, d};
}

Tower tower_mul(const Tower& a, const Tower& b)
{
    const OneRoot d = common_d(a, b);
    require_same_base(a.m, b.m);
    return {a.m * b.m + a.n * b.n * d, a.m * b.n + a.n * b.m, d};
}

Tower tower_scale(const Tower& a, const OneRoot& s)
{
    require_same_base(a.m, s);
    return {a.m * s, a.n * s, a.d};
}

CGAL::Sign sign_tower_root(const Tower& P, const Tower& Q, const FT& e)
{
    if (CGAL::sign(e) == CGAL::NEGATIVE) throw CoverageCertificateNativeError("tower-root radicand is negative");
    const CGAL::Sign sp = sign_tower(P);
    if (CGAL::sign(e) == CGAL::ZERO) return sp;
    const CGAL::Sign sq = sign_tower(Q);
    CGAL::Sign exact = CGAL::ZERO;
    if (sq == CGAL::ZERO) {
        exact = sp;
    } else if (sp == CGAL::ZERO || sp == sq) {
        exact = sq;
    } else {
        const CGAL::Sign c = sign_tower(tower_sub(tower_mul(P, P), tower_scale(tower_mul(Q, Q), OneRoot(e))));
        exact = c == CGAL::POSITIVE ? sp : (c == CGAL::NEGATIVE ? sq : CGAL::ZERO);
    }
    return require_sign_agrees(exact, interval_of(P) + interval_of(Q) * sqrt_enclosure(interval_of(e)), "tower-root");
}

Tower cross_product(const OneRoot& a, const OneRoot& b)
{
    if (!a.is_extended()) return Tower::of(b * a.a0());
    return {b * a.a0(), b * a.a1(), OneRoot(a.root())};
}

Interval tight_interval(const OneRoot& v)
{
    if (!v.is_extended()) return interval_of(v.a0());
    const CGAL::Sign sa = CGAL::sign(v.a0());
    const CGAL::Sign sb = CGAL::sign(v.a1());
    const Interval root = sqrt_enclosure(interval_of(v.root()));
    if (sa == CGAL::ZERO || sb == CGAL::ZERO || sa == sb) return interval_of(v.a0()) + interval_of(v.a1()) * root;
    const FT norm = v.a0() * v.a0() - v.a1() * v.a1() * v.root();
    return interval_of(norm) / (interval_of(v.a0()) - interval_of(v.a1()) * root);
}

Interval tight_interval(const Tower& t)
{
    const Interval m = tight_interval(t.m);
    const Interval n = tight_interval(t.n);
    const Interval root = sqrt_enclosure(tight_interval(t.d));
    const bool cancel = (m.inf() > 0 && n.sup() < 0) || (m.sup() < 0 && n.inf() > 0);
    if (!cancel || sign_one_root(t.d) == CGAL::ZERO) return m + n * root;
    // (m^2 - n^2 d) / (m - n sqrt(d)): the numerator is a same-field value.
    return tight_interval(t.m * t.m - t.n * t.n * t.d) / (m - n * root);
}

Interval interval_of(const Tower& t)
{
    const Interval d = interval_of(t.d);
    if (d.sup() < 0) throw CoverageCertificateNativeError("tower radicand enclosure is negative");
    const Interval root = CGAL::sqrt(Interval(0.0 > d.inf() ? 0.0 : d.inf(), d.sup()));
    return interval_of(t.m) + interval_of(t.n) * root;
}

}  // namespace compas_cgal::coverage_certificate
