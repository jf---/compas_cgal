#include "engagement_2.h"
#include "stock_2.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <tuple>
#include <utility>
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
// Idiom (docs/exactness.md "Numeric comparison is the exact-kernel idiom"):
// compare at the NUMBER-TYPE level. Sqrt_extension is RealEmbeddable, so CGAL::sign and
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

// Gap-closure pessimism: absorb every VOID gap between consecutive engaged runs
// whose angular span does NOT exceed the gap-closure angle gamma (surrogate
// gap_close_ratio, exact threshold T_gamma = gap_close_ratio * r^2), returning
// the resulting PESSIMISTIC runs as (ccw_start, ccw_end) endpoint pairs for the
// cap DECISION. `runs` are the maximal true runs, CCW-sorted, none sharing an
// endpoint (abutting runs were already merged), so between run[i] and
// run[(i+1)%n] lies exactly one void gap -- the CCW arc from run[i].ccw_end to
// run[(i+1)%n].ccw_start (for n == 1 the single gap is the run's complement).
//
// EXACTNESS. A gap is itself an arc with exact one-root endpoints; it is absorbed
// iff its span does not exceed gamma, decided by run_exceeds_cap VERBATIM (the
// identical exact orientation + squared-chord machinery the run cap uses) -- no
// angle is ever summed. Absorbing gap[i] merges run[i] and run[(i+1)%n]; a chain
// of runs joined by absorbed gaps is a SINGLE CCW arc from the chain's first
// ccw_start to its last ccw_end (this is what keeps the whole construction exact:
// the pessimistic run is one arc, so the same endpoint predicates apply directly
// with no per-gap accumulation). Closures chain transitively -- A-gap-B-gap-C with
// both gaps absorbed becomes one run A.start..C.end -- and because each gap's span
// is fixed (its endpoints are frozen run endpoints, unmoved by a neighbour's
// closure) deciding every gap once and then merging connected runs is exactly the
// iterate-to-fixpoint the brief describes, with no iteration needed.
//
// gamma <= pi by contract, so a gap wider than pi (orientation NEGATIVE in
// run_exceeds_cap) never closes. When EVERY gap is absorbed the rim is
// pessimistically full: a degenerate start == end pair is returned, the
// closed-loop case run_exceeds_cap reports as exceeding any cap <= pi.
std::vector<std::pair<GpsPoint, GpsPoint>>
pessimistic_runs(const std::vector<Arc>& runs, const FT& cx, const FT& cy,
                 const FT& gap_close_ratio, const FT& T_gamma)
{
    const std::size_t n = runs.size();
    // A lone run already closing the loop (2*pi, ccw_start == ccw_end) has no gap.
    if (n == 1 && runs[0].ccw_start == runs[0].ccw_end)
        return {{runs[0].ccw_start, runs[0].ccw_end}};

    // Decide every gap once. gap_closed[i] absorbs the void from run[i].ccw_end to
    // run[(i+1)%n].ccw_start. Post-merge these endpoints are distinct, so the gap
    // never degenerates to run_exceeds_cap's p == q full-loop branch.
    std::vector<bool> gap_closed(n);
    std::size_t open_count = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const GpsPoint& g_start = runs[i].ccw_end;
        const GpsPoint& g_end = runs[(i + 1) % n].ccw_start;
        gap_closed[i] = !run_exceeds_cap(g_start, g_end, cx, cy, gap_close_ratio, T_gamma);
        if (!gap_closed[i]) ++open_count;
    }

    // Every gap absorbed => the whole rim is one pessimistic full-circle run.
    if (open_count == 0)
        return {{runs[0].ccw_start, runs[0].ccw_start}};

    // At least one open gap breaks the circle. Start a fresh pessimistic run at the
    // run whose PRECEDING gap is open, then walk all n runs once: an absorbed
    // preceding gap extends the current chain's end, an open one closes the chain
    // and starts the next. The open gap preceding the start run is never rejoined,
    // so no chain spuriously wraps the seam.
    std::size_t s = 0;
    while (gap_closed[(s + n - 1) % n]) ++s;   // guaranteed to halt: open_count > 0

    std::vector<std::pair<GpsPoint, GpsPoint>> pess;
    GpsPoint cur_start = runs[s].ccw_start;
    GpsPoint cur_end = runs[s].ccw_end;
    for (std::size_t step = 1; step < n; ++step) {
        const std::size_t i = (s + step) % n;
        if (gap_closed[(i + n - 1) % n]) {
            cur_end = runs[i].ccw_end;               // absorbed gap: extend the chain
        } else {
            pess.emplace_back(cur_start, cur_end);   // open gap: close off, restart
            cur_start = runs[i].ccw_start;
            cur_end = runs[i].ccw_end;
        }
    }
    pess.emplace_back(cur_start, cur_end);
    return pess;
}

// ----------------------------------------------------------------------------
// Task 5: exact-station TEA cap certificate along a linear cutter motion.
// ----------------------------------------------------------------------------

// Conservative bound on how far a single run's TEA can grow over a center travel
// `d` (tool radius r, stock frozen): the factor-1 analytic lemma. Two mechanisms
// move an existing run's angular extent between two nearby stations.
//
//   (a) ENDPOINT DRIFT. Each end of an existing engaged run sits where the rim
//       crosses a material boundary. Translating the center by d slides such a
//       crossing along the rim by at most the arc subtended by a chord of length
//       d on the radius-r circle, 2*asin(min(1, d/(2r))). A run has two ends, so
//       its extent can grow by up to 4*asin(min(1, d/(2r))).
//   (b) NEWBORN CONTACT. A contact absent at a station can appear in between when
//       a feature first bites the rim. Over travel d the deepest first bite
//       reaches radial depth d into the disk; the rim chord it cuts spans a full
//       angle 2*acos(max(-1, 1 - d/r)).
//
// GROWTH(d, r) = (a) + (b). Both terms are monotone non-decreasing in d and
// saturate at d = 2r, so GROWTH is monotone -- required so that growth over the
// (variable) nearest-station distance is bounded by growth over the half-spacing.
// Evaluated in doubles: an analytic REFINEMENT bound, never a geometric decision
// (docs/exactness.md, "Analytic bounds are not precision handling").
//
//   (c) RUN MERGE -- and why (a)+(b) suffice with NO merge term. A third event
//       changes the LARGEST run: two runs separated by a thin void gap on the rim
//       fuse into one when that gap closes as the cutter advances. This is an O(1)
//       jump in max_run (by min(|A|, |B|) of the two fused runs) reachable within
//       an arbitrarily small step -- unbounded by GROWTH, which is O(sqrt d). A
//       merge term is therefore impossible; the certificate instead removes the
//       event at the source with GAP-CLOSURE PESSIMISM (see pessimistic_runs):
//       each station is measured with every void gap of span <= gamma_guard
//       pre-absorbed, gamma_guard = 2*GROWTH(hs) at half-spacing hs.
//
//       CLAIM: with that pre-closure at both stations, no merge completing within
//       the step is invisible, so max_run at any interior center P is bounded by
//       the PESSIMISTIC max-run at the nearer station S plus the ordinary (a)+(b)
//       growth -- the same shape the single-run lemma already certifies.
//       ARGUMENT (contradiction): suppose the true run through P is the fusion of
//       runs that were SEPARATE at S, across a gap G still open at S. |P - S| <= hs,
//       so G closes from its span sigma(S) > 0 to 0 over travel <= hs. G's two ends
//       are run endpoints; each drifts along the rim by <= 2*asin(min(1, hs/(2r)))
//       (mechanism a, per end) and a newborn bridging G spans <= 2*acos(max(-1,
//       1 - hs/r)) (mechanism b), so the most G can shrink over hs is
//       4*asin(min(1, hs/(2r))) + 2*acos(max(-1, 1 - hs/r)) = GROWTH(hs). Hence
//       sigma(S) <= GROWTH(hs) <= 2*GROWTH(hs) = gamma_guard -- but a gap of span
//       <= gamma_guard is ABSORBED in S's pessimistic measurement, i.e. those runs
//       were ALREADY counted as one at S. That contradicts "separate at S". So the
//       pessimistic run at S already spans the fused arc, and P differs from it by
//       (a)+(b) only. The guarded cap subtracts 2*GROWTH(hs) while this accounting
//       needs only GROWTH(hs): span(run at P) <= pess_max_run(S) + GROWTH(hs) <=
//       (cap - 2*GROWTH(hs)) + GROWTH(hs) = cap - GROWTH(hs) <= cap. QED.
//
//       SAFE FAILURE DIRECTION. gamma_guard = 2*GROWTH(hs) over-closes (only
//       GROWTH(hs) is strictly required), and closing MORE gaps only enlarges the
//       pessimistic runs, making the exact station test STRICTER -- forcing extra
//       refinement or a conservative "uncertified" verdict, never a false pass.
//       Reported total_tea/max_run_tea stay the TRUE (unclosed) measures, so the
//       pessimism inflates only the decision, never the numbers shown to humans.
double tea_growth_bound(double d, double r)
{
    const double a = 4.0 * std::asin(std::min(1.0, d / (2.0 * r)));
    const double b = 2.0 * std::acos(std::max(-1.0, 1.0 - d / r));
    return a + b;
}

// Explicit integer safety factor applied to GROWTH to form the certificate's
// guard. The guard need only bound the TRUE growth; multiplying by 2 buries both
// any looseness in the lemma and the ~1e-15 relative error of asin/acos under
// proof-level slack (CLAUDE.md analytic-bounds clause). SAFE FAILURE DIRECTION:
// too LARGE a guard only forces extra refinement or a conservative "uncertified"
// verdict -- it can NEVER certify a motion that violates the cap. (Too small a
// guard could hide a violation; the factor exists precisely to forbid that.)
constexpr int TEA_GUARD_SAFETY_FACTOR = 2;

double tea_guard(double d, double r)
{
    return TEA_GUARD_SAFETY_FACTOR * tea_growth_bound(d, r);
}

// Refinement stops when a segment is shorter than this fraction of the tool
// radius. At r-scale that leaves a residual guard on the order of 1e-1 rad;
// below it the newborn-contact term (which falls off like sqrt(d)) shrinks so
// slowly that halving again barely moves the guard, so recursion cannot close a
// still-open margin. Reaching the floor with the margin open reports
// uncertified -- the safe direction.
constexpr double STATION_FLOOR_FRACTION = 1e-3;

// Belt-and-braces recursion bound, redundant with the spacing floor for any
// non-degenerate segment but guaranteeing termination regardless of scale.
constexpr int CERTIFY_MAX_DEPTH = 24;

// Certify TEA <= cap for every center on segment [(x0,y0),(x1,y1)], accumulating
// into `acc`. INVARIANT: acc.cap_certified stays true until counter-evidence is
// found, after which callers must stop -- one uncertified sub-span condemns the
// whole motion. cap_chord_ratio is the precomputed FULL-cap surrogate, reused
// for max_tea reporting on segments too coarse to guard.
void certify_recursive(const Stock2& stock, double x0, double y0, double x1,
                       double y1, double r, double cap, double cap_chord_ratio,
                       CertifiedTea& acc, int depth)
{
    const double seg_len = std::hypot(x1 - x0, y1 - y0);
    const double half_spacing = 0.5 * seg_len;
    const double guard = tea_guard(half_spacing, r);
    const double cap_guarded = cap - guard;
    acc.stations += 1;

    if (cap_guarded > 0.0) {
        // BOUNDARY: the guarded transcendental cap becomes its exact rational
        // chord surrogate here, computed as a double and injected exactly by
        // engagement_at. cap_guarded in (0, cap] subset (0, pi] => this ratio is
        // in (0, 4], the exact predicate's valid range. The station VERDICT
        // (cap_exceeded) is exact; only the threshold SELECTION is analytic.
        const double sg = std::sin(0.5 * cap_guarded);
        const double guarded_ratio = 4.0 * sg * sg;
        // GAP-CLOSURE GUARD (merge repair). gamma_guard = 2*GROWTH(half_spacing)
        // == guard (the existing factor-2 tea_guard, reused). Any void gap that
        // could close within one half-spacing of travel has span <= GROWTH(hs) <=
        // gamma_guard (endpoints of the gap are run endpoints moving under the same
        // drift/newborn moduli the growth lemma bounds), so pre-closing every gap
        // <= gamma_guard makes each station's PESSIMISTIC max-run already account
        // for any merge completing before the next station -- the growth lemma then
        // bridges stations with NO merge term (derivation at tea_growth_bound). The
        // cap min(gamma_guard, pi) keeps the surrogate in [0, 4]; gamma_guard >= pi
        // closes all sub-pi gaps (maximally pessimistic, the safe direction).
        const double gg = std::sin(0.5 * std::min(guard, std::numbers::pi));
        const double gap_close_ratio = 4.0 * gg * gg;
        const EngagementSample e0 = engagement_at(stock, x0, y0, r, guarded_ratio, gap_close_ratio);
        const EngagementSample e1 = engagement_at(stock, x1, y1, r, guarded_ratio, gap_close_ratio);
        acc.max_tea = std::max({acc.max_tea, e0.max_run_tea, e1.max_run_tea});
        if (!e0.cap_exceeded && !e1.cap_exceeded)
            return;   // both stations under the guarded cap => whole span certified
    } else {
        // guard >= cap: no positive guarded cap exists at this spacing, so no
        // certification is possible here. Measure only to report max_tea (a
        // reported double; cap_exceeded against the full cap is deliberately
        // ignored -- the verdict comes only from the guarded test), then refine.
        const EngagementSample e0 = engagement_at(stock, x0, y0, r, cap_chord_ratio);
        const EngagementSample e1 = engagement_at(stock, x1, y1, r, cap_chord_ratio);
        acc.max_tea = std::max({acc.max_tea, e0.max_run_tea, e1.max_run_tea});
    }

    // Not certifiable at this spacing: bisect, unless the floor/depth bound is
    // hit -- then the margin is unclosable and the motion is uncertified.
    if (seg_len < STATION_FLOOR_FRACTION * r || depth >= CERTIFY_MAX_DEPTH) {
        acc.cap_certified = false;
        return;
    }
    const double mx = 0.5 * (x0 + x1), my = 0.5 * (y0 + y1);
    certify_recursive(stock, x0, y0, mx, my, r, cap, cap_chord_ratio, acc, depth + 1);
    if (!acc.cap_certified) return;
    certify_recursive(stock, mx, my, x1, y1, r, cap, cap_chord_ratio, acc, depth + 1);
}

} // namespace

EngagementSample engagement_at(const Stock2& stock, double cx, double cy,
                               double tool_radius, double cap_chord_ratio,
                               double gap_close_ratio)
{
    // API-boundary contract: cap_chord_ratio = 4*sin^2(cap/2) with 0 < cap <= pi
    // lies in (0, 4]. Validate the raw double before exact injection (NaN fails).
    if (!(cap_chord_ratio > 0.0 && cap_chord_ratio <= 4.0))
        throw std::invalid_argument("cap_chord_ratio must be in (0, 4] (= 4*sin^2(cap/2), 0 < cap <= pi).");
    // gap_close_ratio = 4*sin^2(gamma/2) with 0 <= gamma <= pi lies in [0, 4];
    // 0 (no gap closed) is the default and the pre-pessimism semantics.
    if (!(gap_close_ratio >= 0.0 && gap_close_ratio <= 4.0))
        throw std::invalid_argument("gap_close_ratio must be in [0, 4] (= 4*sin^2(gamma/2), 0 <= gamma <= pi).");

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

            // Tangent-touch degeneracy (docs/exactness.md "Degeneracy is an
            // ordinary input case"): a zero-measure contact. CGAL never emits a
            // degenerate (source == target) x-monotone arc, but guard exactly so
            // one could never be mistaken for a full 2*pi run by span wrap-around.
            if (s == t) continue;

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

    // 4a. REPORTING (doubles) over the TRUE runs -- never gap-closed. total_tea and
    //     max_run_tea describe the material actually engaged at this station.
    for (const Arc& run : runs) {
        out.total_tea += run.span;
        out.max_run_tea = std::max(out.max_run_tea, run.span);
    }

    // 4b. DECISION (exact) over the PESSIMISTIC runs. Void gaps <= gamma are
    //     absorbed (default gamma == 0 closes none => pessimistic == true runs, so
    //     the verdict is bit-for-bit the pre-pessimism result). The cap predicate
    //     runs against the ACTUAL cap threshold T; the gap-closure predicate uses
    //     the separate gamma threshold T_gamma. A pessimistic run contains its true
    //     runs, and run_exceeds_cap is monotone in the arc, so testing the
    //     pessimistic runs alone is conservative -- any true run over the cap forces
    //     its pessimistic superset over it too.
    const FT ratio_ft(cap_chord_ratio);
    const FT T = ratio_ft * r_sq;
    const FT gap_ratio_ft(gap_close_ratio);
    const FT T_gamma = gap_ratio_ft * r_sq;
    const FT cxf(cx), cyf(cy);
    for (const auto& [start, end] :
         pessimistic_runs(runs, cxf, cyf, gap_ratio_ft, T_gamma)) {
        if (run_exceeds_cap(start, end, cxf, cyf, ratio_ft, T)) {
            out.cap_exceeded = true;
            break;
        }
    }
    return out;
}

CertifiedTea certify_segment_tea(const Stock2& stock, double x0, double y0,
                                 double x1, double y1, double tool_radius,
                                 double cap_radians)
{
    // BOUNDARY (docs/exactness.md, boundary doctrine): validate the ergonomic
    // angular cap and convert it to its exact chord surrogate here, at the one
    // declared seam. Contract 0 < cap <= pi: a single engaged run subtends at
    // most a half turn before the >pi case is an exact orientation verdict, so a
    // cap above pi is meaningless. The full-cap ratio is carried down for max_tea
    // reporting; per-level guarded ratios are derived inside the recursion.
    if (!(cap_radians > 0.0 && cap_radians <= std::numbers::pi))
        throw std::invalid_argument("cap_radians must be in (0, pi].");

    const double sc = std::sin(0.5 * cap_radians);
    const double cap_chord_ratio = 4.0 * sc * sc;

    CertifiedTea acc{0.0, true, 0};
    certify_recursive(stock, x0, y0, x1, y1, tool_radius, cap_radians,
                      cap_chord_ratio, acc, 0);
    return acc;
}

void register_engagement(nanobind::module_& m)
{
    // Nanobind boundary. cx, cy, tool_radius are measured/computed station data:
    // each double IS a rational and enters exact-land by exact injection
    // (Epeck::FT) -- no snapping, no tolerance at the seam. cap_chord_ratio is
    // the caller's exact rational surrogate 4*sin^2(cap/2) for the transcendental
    // cap (docs/exactness.md "Input semantics" and "boundary doctrine").
    m.def("engagement_at",
          [](const Stock2& stock, double cx, double cy, double tool_radius,
             double cap_chord_ratio, double gap_close_ratio) {
              EngagementSample s = engagement_at(stock, cx, cy, tool_radius,
                                                 cap_chord_ratio, gap_close_ratio);
              return std::make_tuple(s.total_tea, s.max_run_tea, s.cap_exceeded);
          },
          "stock"_a, "cx"_a, "cy"_a, "tool_radius"_a, "cap_chord_ratio"_a,
          "gap_close_ratio"_a = 0.0);

    // Nanobind boundary for the motion certificate: x0..y1, tool_radius and the
    // angular cap enter exact-land inside certify_segment_tea (validation + exact
    // chord-surrogate injection). Returns (max_tea, cap_certified, stations),
    // the CertifiedTea fields, matching the tuple style of engagement_at.
    m.def("certify_segment_tea",
          [](const Stock2& stock, double x0, double y0, double x1, double y1,
             double tool_radius, double cap_radians) {
              CertifiedTea c = certify_segment_tea(stock, x0, y0, x1, y1,
                                                   tool_radius, cap_radians);
              return std::make_tuple(c.max_tea, c.cap_certified, c.stations);
          },
          "stock"_a, "x0"_a, "y0"_a, "x1"_a, "y1"_a, "tool_radius"_a, "cap_radians"_a);

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
