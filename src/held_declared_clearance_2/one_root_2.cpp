#include "held_declared_clearance_2/one_root_2.h"

#include "held_declared_clearance_2/errors_2.h"

#include <CGAL/Fraction_traits.h>


#include <utility>

namespace compas_cgal::held_declared_clearance {

namespace {

using ExactRational = FT::ET;
using RationalTraits = CGAL::Fraction_traits<ExactRational>;
using Integer = RationalTraits::Numerator_type;

/// The rational square root of `value` when it has one (value >= 0).
bool rational_sqrt(const FT& value, FT& root)
{
    Integer numerator;
    Integer denominator;
    RationalTraits::Decompose()(CGAL::exact(value), numerator, denominator);
    if (denominator < 0) {
        numerator = -numerator;
        denominator = -denominator;
    }
    Integer numerator_root;
    Integer denominator_root;
    if (!CGAL::is_square(numerator, numerator_root) || !CGAL::is_square(denominator, denominator_root)) return false;
    root = FT(RationalTraits::Compose()(numerator_root, denominator_root));
    return true;
}

/// A value of the extension sqrt(c): folded to a rational when b vanished.
Root in_extension(const FT& a, const FT& b, const FT& c)
{
    if (CGAL::is_zero(b)) return Root(a);
    return Root(a, b, c);
}

void require_one_extension(const Root& left, const Root& right)
{
    if (left.is_extended() && right.is_extended() && CGAL::compare(left.root(), right.root()) != CGAL::EQUAL) {
        throw CrossRootArithmeticError("arithmetic across different quadratic extensions");
    }
}

const FT& root_of(const Root& left, const Root& right) { return left.is_extended() ? left.root() : right.root(); }

/// A + B sqrt(alpha) + C sqrt(beta) + D sqrt(alpha) sqrt(beta): a sum of
/// products of numbers of two extensions, kept for its sign only.
struct TwoRootLinear {
    FT a, b, c, d;
    Root alpha_holder;  // a number of the first extension (or rational)
    Root beta_holder;   // a number of the second extension (or rational)

    static TwoRootLinear product(const Root& left, const Root& right)
    {
        return {left.a0() * right.a0(), left.a1() * right.a0(), left.a0() * right.a1(), left.a1() * right.a1(),
                left, right};
    }

    TwoRootLinear plus(const TwoRootLinear& other, int scale) const
    {
        require_one_extension(alpha_holder, other.alpha_holder);
        require_one_extension(beta_holder, other.beta_holder);
        return {a + scale * other.a, b + scale * other.b, c + scale * other.c, d + scale * other.d,
                alpha_holder.is_extended() ? alpha_holder : other.alpha_holder,
                beta_holder.is_extended() ? beta_holder : other.beta_holder};
    }

    /// Exact sign by grouping (A + B sqrt(alpha)) + sqrt(beta) (C + D sqrt(alpha)).
    int sign() const
    {
        const FT alpha = alpha_holder.is_extended() ? alpha_holder.root() : FT(0);
        const FT beta = beta_holder.is_extended() ? beta_holder.root() : FT(0);
        return sign_of_sum(in_extension(a, b, alpha), in_extension(c, d, alpha), beta);
    }
};

std::pair<Root, Root> from_centre(const RationalPoint& centre, const OneRootPoint& point)
{
    return {shifted(point.x(), -centre.x()), shifted(point.y(), -centre.y())};
}

int dot_sign(const RationalPoint& centre, const OneRootPoint& first, const OneRootPoint& second)
{
    const auto [fx, fy] = from_centre(centre, first);
    const auto [sx, sy] = from_centre(centre, second);
    return TwoRootLinear::product(fx, sx).plus(TwoRootLinear::product(fy, sy), 1).sign();
}

/// 0 for a swept angle in [0, pi) from `start` in `sense`, 1 for [pi, 2 pi).
int half_turn(const RationalPoint& centre, int sense, const OneRootPoint& start, const OneRootPoint& point)
{
    const int side = sense * orientation(centre, start, point);
    if (side > 0) return 0;
    if (side < 0) return 1;
    return dot_sign(centre, start, point) > 0 ? 0 : 1;
}

int half_turn(int sense, const Vector& start, const Vector& direction)
{
    const int side = sense * static_cast<int>(CGAL::orientation(start, direction));
    if (side > 0) return 0;
    if (side < 0) return 1;
    return CGAL::sign(start * direction) == CGAL::POSITIVE ? 0 : 1;
}

}  // namespace

int sign(const FT& value) { return static_cast<int>(CGAL::sign(value)); }

int sign(const Root& value) { return static_cast<int>(CGAL::sign(value)); }

int compare(const Root& left, const Root& right) { return static_cast<int>(CGAL::compare(left, right)); }

Root rational(const FT& value) { return Root(value); }

Root make_root(const FT& a, const FT& b, const FT& c)
{
    if (CGAL::is_negative(c)) throw NegativeRadicandError("a one-root number was given a negative radicand");
    if (CGAL::is_zero(b) || CGAL::is_zero(c)) return Root(a);
    FT root;
    if (rational_sqrt(c, root)) return Root(a + b * root);
    return Root(a, b, c);
}

Root make_deciding_root(const FT& a, const FT& b, const FT& c)
{
    if (CGAL::is_negative(c)) throw NegativeRadicandError("a one-root number was given a negative radicand");
    if (CGAL::is_zero(b) || CGAL::is_zero(c)) return Root(a);
    return Root(a, b, c);
}

bool is_rational(const Root& value) { return !value.is_extended(); }

Root plus(const Root& left, const Root& right)
{
    require_one_extension(left, right);
    return in_extension(left.a0() + right.a0(), left.a1() + right.a1(), root_of(left, right));
}

Root minus(const Root& left, const Root& right)
{
    require_one_extension(left, right);
    return in_extension(left.a0() - right.a0(), left.a1() - right.a1(), root_of(left, right));
}

Root times(const Root& left, const Root& right)
{
    require_one_extension(left, right);
    const FT& root = root_of(left, right);
    return in_extension(left.a0() * right.a0() + left.a1() * right.a1() * root,
                        left.a0() * right.a1() + left.a1() * right.a0(), root);
}

Root scaled(const FT& scale, const Root& value)
{
    if (!value.is_extended()) return Root(scale * value.a0());
    return in_extension(scale * value.a0(), scale * value.a1(), value.root());
}

Root shifted(const Root& value, const FT& shift)
{
    if (!value.is_extended()) return Root(value.a0() + shift);
    return Root(value.a0() + shift, value.a1(), value.root());
}

int sign_of_sum(const Root& u, const Root& w, const FT& beta)
{
    const int sign_u = sign(u);
    const int sign_w = sign(w);
    if (CGAL::is_zero(beta) || sign_w == 0) return sign_u;
    if (sign_u == 0 || sign_u == sign_w) return sign_u == 0 ? sign_w : sign_u;
    // Opposite signs: |u| against |w| sqrt(beta), i.e. u^2 against beta w^2 in one extension.
    const int magnitude = compare(times(u, u), scaled(beta, times(w, w)));
    return sign_u * magnitude;
}

OneRootPoint make_point(const Root& x, const Root& y)
{
    if (x.is_extended() && y.is_extended() && CGAL::compare(x.root(), y.root()) != CGAL::EQUAL) {
        throw MixedRootPointError("a one-root point's coordinates lie in different extensions");
    }
    return OneRootPoint(x, y);
}

OneRootPoint embed(const RationalPoint& point) { return OneRootPoint(Root(point.x()), Root(point.y())); }

bool is_rational(const OneRootPoint& point) { return !point.x().is_extended() && !point.y().is_extended(); }

RationalPoint rational_point(const OneRootPoint& point) { return RationalPoint(point.x().a0(), point.y().a0()); }

int compare_xy(const OneRootPoint& first, const OneRootPoint& second)
{
    const int by_x = compare(first.x(), second.x());
    return by_x != 0 ? by_x : compare(first.y(), second.y());
}

int orientation(const RationalPoint& centre, const OneRootPoint& first, const OneRootPoint& second)
{
    if (is_rational(first) && is_rational(second)) {
        return static_cast<int>(CGAL::orientation(centre, rational_point(first), rational_point(second)));
    }
    const auto [fx, fy] = from_centre(centre, first);
    const auto [sx, sy] = from_centre(centre, second);
    return TwoRootLinear::product(fx, sy).plus(TwoRootLinear::product(fy, sx), -1).sign();
}

int compare_on_arc(const RationalPoint& centre, int sense, const OneRootPoint& start, const OneRootPoint& first,
                   const OneRootPoint& second)
{
    const int first_half = half_turn(centre, sense, start, first);
    const int second_half = half_turn(centre, sense, start, second);
    if (first_half != second_half) return first_half < second_half ? -1 : 1;
    return -sense * orientation(centre, first, second);
}

int compare_on_arc(int sense, const Vector& start, const Vector& first, const Vector& second)
{
    const int first_half = half_turn(sense, start, first);
    const int second_half = half_turn(sense, start, second);
    if (first_half != second_half) return first_half < second_half ? -1 : 1;
    return -sense * static_cast<int>(CGAL::orientation(first, second));
}

}  // namespace compas_cgal::held_declared_clearance
