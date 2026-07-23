#include "engagement_2.h"
#include "stock_2.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <CGAL/enum.h>
#include <CGAL/number_utils.h>

namespace {

using FT = Epeck::FT;
using CoordNT = GpsPoint::CoordNT;   // Sqrt_extension<FT, FT>: a0 + a1*sqrt(root)

// ----------------------------------------------------------------------------
// Exact sign of a mixed two-radical form   A + B*sqrt(alpha) + C*sqrt(beta)
//                                            + D*sqrt(alpha*beta)
// with RATIONAL A, B, C, D and RATIONAL alpha, beta >= 0. The exact cap
// predicate reduces to this: orientation and squared-chord tests of two
// cutter-circle points p, q whose coordinates live in Q(sqrt alpha) and
// Q(sqrt beta) respectively expand to exactly this shape.
//
// Idiom (CLAUDE.md "Numeric comparison is the exact-kernel idiom"): compare at
// the NUMBER-TYPE level. Sqrt_extension is RealEmbeddable, so CGAL::sign and
// same-root CGAL::compare (and same-root +/-/*) are exact -- we build the
// derived quantities INSIDE one extension Q(sqrt alpha) and let CGAL decide,
// rather than hand-rolling a bignum squaring routine. Cross-root Sqrt_extension
// arithmetic (documented UB) is never formed: sqrt(beta) only ever appears as a
// squared factor `beta` (a rational), and the degenerate roots are folded away
// with exact FT compares before any extension is built.
CGAL::Sign sign_mixed_radical(const FT& A, const FT& B, const FT& C, const FT& D,
                              const FT& alpha, const FT& beta)
{
    const bool alpha_ext = !CGAL::is_zero(alpha);
    const bool beta_ext = !CGAL::is_zero(beta);

    // Root degeneracies fold to a single same-root value -> one CGAL::sign.
    if (!alpha_ext && !beta_ext) return CGAL::sign(A);              // fully rational
    if (!beta_ext) return CGAL::sign(CoordNT(A, B, alpha));         // A + B*sqrt(alpha)
    if (!alpha_ext) return CGAL::sign(CoordNT(A, C, beta));         // A + C*sqrt(beta)
    if (alpha == beta)                                             // sqrt(alpha*beta) = alpha
        return CGAL::sign(CoordNT(A + D * alpha, B + C, alpha));    // (A+D*alpha) + (B+C)*sqrt(alpha)

    // General case: group over the shared root alpha into u, w in Q(sqrt alpha),
    // so the form is u + sqrt(beta)*w. Its sign follows from sign(u), sign(w),
    // and -- for opposite non-zero signs -- which magnitude dominates, decided
    // exactly by compare(u^2, beta*w^2) (all same-root, so beta enters as the
    // rational CoordNT(beta), never as sqrt(beta)).
    const CoordNT u(A, B, alpha);
    const CoordNT w(C, D, alpha);
    const CGAL::Sign su = CGAL::sign(u);
    const CGAL::Sign sw = CGAL::sign(w);
    if (sw == CGAL::ZERO) return su;                 // sqrt(beta)*w = 0 -> u
    if (su == CGAL::ZERO) return sw;                 // u = 0 -> sqrt(beta)*w, beta > 0
    if (su == sw) return su;                         // like signs add
    switch (CGAL::compare(u * u, w * w * CoordNT(beta))) {
        case CGAL::LARGER:  return su;               // |u| dominates
        case CGAL::SMALLER: return sw;               // sqrt(beta)*|w| dominates
        default:            return CGAL::ZERO;        // equal magnitude, opposite sign
    }
}

// A cutter-circle point in rational coordinates over its single root:
// (x0 + x1*sqrt(root), y0 + y1*sqrt(root)). CGAL builds every circle/line and
// circle/circle intersection point with a single shared discriminant, so both
// coordinates of one point share one root; a rational coordinate has root 0 and
// a1 == 0. Distinct points may carry distinct roots (alpha for p, beta for q).
struct RadPoint {
    FT x0, x1, y0, y1, root;
};

RadPoint as_radpoint(const GpsPoint& p)
{
    const CoordNT& X = p.x();
    const CoordNT& Y = p.y();
    FT root(0), x1(0), y1(0);
    if (X.is_extended()) { root = X.root(); x1 = X.a1(); }
    if (Y.is_extended()) {
        // Both coordinates extended -> CGAL guarantees the shared root.
        CGAL_assertion(!X.is_extended() || Y.root() == root);
        root = Y.root();
        y1 = Y.a1();
    }
    return { X.a0(), x1, Y.a0(), y1, root };
}

// Exact cap certificate for one maximal engaged run with CCW endpoints p, q on
// the cutter circle centred at (cx, cy). T = cap_chord_ratio * tool_radius^2 is
// the exact squared-chord threshold. Returns true iff the run's angular extent
// exceeds the cap.
bool run_exceeds_cap(const GpsPoint& p, const GpsPoint& q, const FT& cx,
                     const FT& cy, const FT& cap_chord_ratio, const FT& T)
{
    if (p == q) return true;   // closed loop: full 2*pi run exceeds any cap <= pi

    const RadPoint P = as_radpoint(p);
    const RadPoint Q = as_radpoint(q);

    // orientation(center, p, q) = sign of the determinant of (p - center) and
    // (q - center). Writing p - center = (ax + bx*sqrt(alpha), ay + by*sqrt(alpha))
    // and q - center = (c + d*sqrt(beta), e + f*sqrt(beta)), the determinant is
    // A + B*sqrt(alpha) + C*sqrt(beta) + D*sqrt(alpha*beta).
    const FT ax = P.x0 - cx, bx = P.x1;
    const FT ay = P.y0 - cy, by = P.y1;
    const FT c = Q.x0 - cx, d = Q.x1;
    const FT e = Q.y0 - cy, f = Q.y1;
    const FT A = ax * e - ay * c;
    const FT B = bx * e - by * c;
    const FT C = ax * f - ay * d;
    const FT D = bx * f - by * d;
    const CGAL::Sign orient = sign_mixed_radical(A, B, C, D, P.root, Q.root);

    // CCW angle theta from p to q: orient > 0 <=> theta < pi, orient == 0 <=>
    // theta == pi, orient < 0 <=> theta > pi.
    if (orient == CGAL::NEGATIVE) return true;                 // theta > pi >= cap
    if (orient == CGAL::ZERO) return cap_chord_ratio < FT(4);  // theta == pi

    // theta < pi: chord grows monotonically, so exceeded <=> |pq|^2 > T. With
    // p - q coordinates carrying both sqrt(alpha) and sqrt(beta), |pq|^2 expands
    // to A2 + B2*sqrt(alpha) + C2*sqrt(beta) + D2*sqrt(alpha*beta).
    const FT ex = Q.x0 - P.x0;
    const FT ey = Q.y0 - P.y0;
    const FT A2 = ex * ex + ey * ey
                  + (P.x1 * P.x1 + P.y1 * P.y1) * P.root
                  + (Q.x1 * Q.x1 + Q.y1 * Q.y1) * Q.root;
    const FT B2 = -2 * (ex * P.x1 + ey * P.y1);
    const FT C2 = 2 * (ex * Q.x1 + ey * Q.y1);
    const FT D2 = -2 * (P.x1 * Q.x1 + P.y1 * Q.y1);
    return sign_mixed_radical(A2 - T, B2, C2, D2, P.root, Q.root) == CGAL::POSITIVE;
}

// One rim arc, CCW-normalized: the CCW sweep runs ccw_start -> ccw_end on the
// cutter circle. `span` is the arc's angular extent in radians, a REPORTING
// double (atan2) that never feeds a decision.
struct Arc {
    GpsPoint ccw_start;
    GpsPoint ccw_end;
    double span;
};

} // namespace

EngagementSample engagement_at(const Stock2& stock, double cx, double cy,
                               double tool_radius, double cap_chord_ratio)
{
    // API-boundary contract: cap_chord_ratio = 4*sin^2(cap/2) with 0 < cap <= pi
    // lies in (0, 4]. Validate the raw double before exact injection (NaN fails).
    if (!(cap_chord_ratio > 0.0 && cap_chord_ratio <= 4.0))
        throw std::invalid_argument("cap_chord_ratio must be in (0, 4] (= 4*sin^2(cap/2), 0 < cap <= pi).");

    EngagementSample out{0.0, 0.0, false};

    // 1. region = stock intersect tool disk (regularized, exact).
    Gps disk_set;
    disk_set.insert(disk_polygon(EPoint(cx, cy), FT(tool_radius)));
    Gps region = stock.set();
    region.intersection(disk_set);
    if (region.is_empty()) return out;

    // 2. Harvest boundary arcs supported by the cutter circle (rim-in-material).
    //    The cutter disk's circle data is rational and GPS boolean ops preserve
    //    an arc's supporting circle verbatim, so this identity test is exact.
    const EPoint center(cx, cy);
    const FT r_sq = FT(tool_radius) * FT(tool_radius);

    std::vector<Arc> arcs;
    auto harvest = [&](const GpsPolygon& poly) {
        for (auto cit = poly.curves_begin(); cit != poly.curves_end(); ++cit) {
            const GpsXCurve& xc = *cit;
            if (!xc.is_circular()) continue;
            if (xc.supporting_circle().center() != center) continue;
            if (xc.supporting_circle().squared_radius() != r_sq) continue;

            // CCW-normalize: orientation() is the exact winding of source->target
            // about the centre. CLOCKWISE means the CCW sweep runs target->source.
            GpsPoint s = xc.source();
            GpsPoint t = xc.target();
            if (xc.orientation() == CGAL::CLOCKWISE) std::swap(s, t);

            const double sx = CGAL::to_double(s.x()), sy = CGAL::to_double(s.y());
            const double tx = CGAL::to_double(t.x()), ty = CGAL::to_double(t.y());
            double span = std::atan2(ty - cy, tx - cx) - std::atan2(sy - cy, sx - cx);
            if (span <= 0.0) span += 2.0 * std::numbers::pi;
            arcs.push_back({s, t, span});
        }
    };

    std::vector<GpsPolygonWithHoles> pwhs;
    region.polygons_with_holes(std::back_inserter(pwhs));
    for (const auto& pwh : pwhs) {
        harvest(pwh.outer_boundary());
        for (auto hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit) harvest(*hit);
    }
    if (arcs.empty()) return out;   // region present but rim nowhere in material

    // 3. Assemble maximal engaged runs. Abutting rim arcs share their split
    //    vertex EXACTLY (same arrangement point), so adjacency is exact point
    //    equality -- no angular gap tolerance. Sort CCW by start using exact
    //    point comparisons only (half-plane about the rational horizontal line
    //    through the centre, then x within each half), so atan2 stays confined
    //    to reporting.
    const CoordNT cx_c{FT(cx)};
    const CoordNT cy_c{FT(cy)};
    auto upper = [&](const GpsPoint& pt) -> bool {
        // Angle in [0, pi): y > cy, or (y == cy and x > cx) (the +x seam point).
        const CGAL::Comparison_result cyc = CGAL::compare(pt.y(), cy_c);
        if (cyc == CGAL::LARGER) return true;
        if (cyc == CGAL::SMALLER) return false;
        return CGAL::compare(pt.x(), cx_c) == CGAL::LARGER;
    };
    auto ccw_less = [&](const Arc& u, const Arc& v) -> bool {
        const bool uu = upper(u.ccw_start);
        const bool uv = upper(v.ccw_start);
        if (uu != uv) return uu;   // upper half [0, pi) precedes lower half [pi, 2pi)
        const CGAL::Comparison_result cx_cmp =
            CGAL::compare(u.ccw_start.x(), v.ccw_start.x());
        // Upper half: angle grows as x decreases; lower half: as x increases.
        return uu ? (cx_cmp == CGAL::LARGER) : (cx_cmp == CGAL::SMALLER);
    };
    std::sort(arcs.begin(), arcs.end(), ccw_less);

    std::vector<Arc> runs;
    for (const Arc& a : arcs) {
        if (!runs.empty() && runs.back().ccw_end == a.ccw_start) {
            runs.back().ccw_end = a.ccw_end;
            runs.back().span += a.span;
        } else {
            runs.push_back(a);
        }
    }
    // Wrap-around: at most one run crosses the +x seam, joining the last run's
    // end to the first run's start (again exact point equality).
    if (runs.size() > 1 && runs.back().ccw_end == runs.front().ccw_start) {
        runs.front().ccw_start = runs.back().ccw_start;
        runs.front().span += runs.back().span;
        runs.pop_back();
    }

    // 4. Report totals (doubles) and certify the cap (exact) per run.
    const FT ratio_ft(cap_chord_ratio);
    const FT T = ratio_ft * r_sq;
    for (const Arc& run : runs) {
        out.total_tea += run.span;
        out.max_run_tea = std::max(out.max_run_tea, run.span);
        if (run_exceeds_cap(run.ccw_start, run.ccw_end, FT(cx), FT(cy), ratio_ft, T))
            out.cap_exceeded = true;
    }
    return out;
}

void register_engagement(nanobind::module_& m)
{
    m.def("engagement_at",
          [](const Stock2& stock, double cx, double cy, double tool_radius, double cap_chord_ratio) {
              EngagementSample s = engagement_at(stock, cx, cy, tool_radius, cap_chord_ratio);
              return std::make_tuple(s.total_tea, s.max_run_tea, s.cap_exceeded);
          },
          "stock"_a, "cx"_a, "cy"_a, "tool_radius"_a, "cap_chord_ratio"_a);

    // Test-only: exact sign of A + B*sqrt(alpha) + C*sqrt(beta) + D*sqrt(alpha*beta)
    // (returns -1/0/+1) so the cap predicate's core primitive is unit-tested
    // directly against high-precision references.
    m.def("_sign_mixed_radical",
          [](double a, double b, double c, double d, double alpha, double beta) {
              return static_cast<int>(
                  sign_mixed_radical(FT(a), FT(b), FT(c), FT(d), FT(alpha), FT(beta)));
          },
          "a"_a, "b"_a, "c"_a, "d"_a, "alpha"_a, "beta"_a);
}
