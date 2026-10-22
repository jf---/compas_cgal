#include "coverage_certificate_2/exact_geometry_2.h"

#include "coverage_certificate_2/errors_2.h"

namespace compas_cgal::coverage_certificate {

CGAL::Sign sign_difference_one_roots(const OneRoot& a, const OneRoot& b)
{
    const FT ba = base_of(a);
    const FT bb = base_of(b);
    if (CGAL::sign(ba) == CGAL::ZERO || CGAL::sign(bb) == CGAL::ZERO || ba == bb) return sign_one_root(a - b);
    return sign_rational_two(a.a0() - b.a0(), a.a1(), ba, -b.a1(), bb);
}

Tower cross_about(const Point_2& c, const OneRoot& wx, const OneRoot& wy, const ExactPoint& p)
{
    return tower_sub(cross_product(wx - c.x(), p.y - c.y()), cross_product(wy - c.y(), p.x - c.x()));
}

Interval cross_about_i(const Interval& cx, const Interval& cy, const Interval& wx, const Interval& wy, const Interval& px,
                       const Interval& py)
{
    return (wx - cx) * (py - cy) - (wy - cy) * (px - cx);
}

Tower squared_distance_minus(const ExactPoint& p, const OneRoot& wx, const OneRoot& wy, const FT& r2)
{
    const FT delta = wx.is_extended() ? wx.root() : (wy.is_extended() ? wy.root() : FT(0));
    const FT x1 = wx.is_extended() ? wx.a1() : FT(0);
    const FT y1 = wy.is_extended() ? wy.a1() : FT(0);
    const OneRoot ex = p.x - wx.a0();
    const OneRoot ey = p.y - wy.a0();
    // (ex - x1 s)^2 + (ey - y1 s)^2 - r2 with s = sqrt(delta).
    const OneRoot m = ex * ex + ey * ey + OneRoot((x1 * x1 + y1 * y1) * delta - r2);
    const OneRoot n = FT(-2) * (ex * x1 + ey * y1);
    return {m, n, OneRoot(delta)};
}

bool in_wedge(const ArcSpan span, const CGAL::Sign o0, const CGAL::Sign o1)
{
    switch (span) {
        case ArcSpan::Minor: return o0 != CGAL::NEGATIVE && o1 != CGAL::NEGATIVE;
        case ArcSpan::Half: return o0 != CGAL::NEGATIVE;
        case ArcSpan::Major: return o0 != CGAL::NEGATIVE || o1 != CGAL::NEGATIVE;
    }
    throw CoverageCertificateNativeError("unknown arc span");
}

CGAL::Sign times(const int sense, const CGAL::Sign s) { return sense > 0 ? s : CGAL::opposite(s); }

ArcSpan span_of(const Point_2& c, const OneRoot& ax, const OneRoot& ay, const OneRoot& bx, const OneRoot& by,
                const int sense)
{
    const ExactPoint b{bx, by};
    const CGAL::Sign s = times(sense, sign_tower(cross_about(c, ax, ay, b)));
    if (s == CGAL::POSITIVE) return ArcSpan::Minor;
    if (s == CGAL::NEGATIVE) return ArcSpan::Major;
    if (sign_difference_one_roots(ax, bx) == CGAL::ZERO && sign_difference_one_roots(ay, by) == CGAL::ZERO) {
        throw InvalidCertificateInputError("an arc of the exact target starts and ends at one point");
    }
    return ArcSpan::Half;
}

bool in_open_wedge(const ArcSpan span, const CGAL::Sign o0, const CGAL::Sign o1)
{
    switch (span) {
        case ArcSpan::Minor: return o0 == CGAL::POSITIVE && o1 == CGAL::POSITIVE;
        case ArcSpan::Half: return o0 == CGAL::POSITIVE;
        case ArcSpan::Major: return o0 == CGAL::POSITIVE || o1 == CGAL::POSITIVE;
    }
    throw CoverageCertificateNativeError("unknown arc span");
}

}  // namespace compas_cgal::coverage_certificate
