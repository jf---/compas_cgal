#include "engagement_2.h"
#include "exact_boundary.h"
#include "stock_2.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <CGAL/enum.h>
#include <CGAL/number_utils.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arrangement_zone_2.h>

#include <boost/variant.hpp>

namespace {

using FT = Epeck::FT;
using CoordNT = GpsPoint::CoordNT;   // Sqrt_extension<FT, FT>: a0 + a1*sqrt(root)

// ----------------------------------------------------------------------------
// Boundary guards for the geometry parameters live in exact_boundary.h, shared
// verbatim with stock_2.cpp so the two halves of the seam cannot drift (the
// doctrine, and what the guards do and do not add, are documented there).
// Refusing both classes at this seam keeps every downstream FT construction
// total, and replaced two defects with named domain errors. An impossible radius
// used to return a confident answer: tool_radius = -1 reported full 2*pi
// immersion, and tool_radius = +Inf reported zero engagement with the cap NOT
// exceeded -- a false pass. A non-finite coordinate used to leak the exact number
// type's internal "Cannot convert a non-finite number to an integer" RuntimeError.
// ----------------------------------------------------------------------------

using exact_boundary::format_double;
using exact_boundary::require_finite;
using exact_boundary::require_positive_radius;

// The mixed-radical primitive's roots are RADICANDS: sqrt(alpha) is real only for
// alpha >= 0, and sign_mixed_radical's contract presumes it. Zero is legal -- it
// is the degenerate "not extended" branch. Spelled `!(v >= 0.0)` to stay NaN-safe
// independently of the finiteness check, matching the ratio guards below.
void require_radicand(double value, const char* name)
{
    require_finite(value, name);
    if (!(value >= 0.0))
        throw std::invalid_argument(std::string(name) + " must be a non-negative radicand (got " + format_double(value) + ").");
}

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
        // LOAD-BEARING PRECONDITION, checked at RUNTIME rather than asserted. It is
        // CGAL's own contract for this number type -- Sqrt_extension::check_roots
        // states it as `CGAL_precondition(a.root() == b.root())` -- and it holds here
        // because CGAL builds every circle/line and circle/circle intersection point
        // from a single shared discriminant. `CGAL_assertion` would state it for a
        // debug build only: NDEBUG makes it `static_cast<void>(0)` (CGAL/assertions.h),
        // and this project compiles Release (pyproject.toml `cmake.build-type`), so an
        // asserted form does not exist in any shipped wheel.
        //
        // Were it ever false, the collapse below would keep the FIRST coordinate's a1
        // while adopting the SECOND coordinate's root -- silently re-interpreting
        // x1*sqrt(alpha) as x1*sqrt(beta). Every orientation and squared-chord sign in
        // run_exceeds_cap is computed from these five rationals, so the cap certificate
        // would come out confidently wrong rather than absent. An exception is the only
        // outcome a caller can see; the exact comparison is on FT, not a tolerance.
        if (X.is_extended() && Y.root() != root)
            throw std::logic_error(
                "as_radpoint: the two coordinates of one arrangement point carry different roots (radicands report as "
                + format_double(CGAL::to_double(root)) + " and " + format_double(CGAL::to_double(Y.root()))
                + "), so CGAL did not build this intersection point from a single shared discriminant. The exact cap "
                  "predicate collapses both coordinates onto one root and cannot do so here; no engagement certificate "
                  "computed from this point would be sound.");
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

// A full turn of engaged rim (radians). Both the reporting CEILING -- the rim is a
// closed curve of that angular measure, so no engaged run can measure more -- and
// the swept bound's SATURATED value: with the cap contractually in (0, pi],
// returning this always forces refinement or a conservative refusal, so it is the
// safe answer whenever the construction cannot bound the swept region below a
// full turn.
constexpr double FULL_TURN = 2.0 * std::numbers::pi;

// Slack on that reporting ceiling, in radians. NOT a geometric tolerance: no
// decision reads it, no geometry is sized by it, and nothing is clamped to it. It
// exists only so that a floating-point SUM of individually-accurate angles is not
// mistaken for the logic error the ceiling check hunts.
//
// FLOOR -- what it must stay above. One arc span is atan2(|u x v|, u . v) over
// station-relative doubles; the products carry a handful of ulps, so a span is
// accurate to a few ulps of pi (ulp(pi) = 4.4e-16 rad), call it 1e-15 rad. A run
// is the sum of the spans of the sub-arcs it was assembled from, and those
// sub-arcs are disjoint pieces of one rim, so a run built from N of them lands
// within N * 1e-15 rad of its true measure. 1e-6 rad absorbs N up to 1e9 sub-arcs
// on a single cutter rim -- nine decades past anything the zone can hand back,
// since every sub-arc costs one rim/boundary crossing in the arrangement.
//
// CEILING -- what it must stay below. The cheapest impossibility this can catch
// is worth a whole turn: a double-counted run, or the wrap normalisation the span
// formula replaced firing on a sub-ulp arc, each add exactly 2*pi = 6.28 rad.
// 1e-6 rad is 6.3 million times smaller, so no real double-count can hide beneath
// it -- measured on the pre-fix build, the three witness stations in
// tests/test_growth_bound.py overshoot this check by 6.28 rad, not by ulps.
constexpr double FULL_TURN_REPORT_SLACK = 1e-6;

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

// Shared engagement tail: given the rim arcs supported by the cutter circle
// (from EITHER harvest -- the overlay boundary or the local zone query), assemble
// the maximal engaged runs by EXACT endpoint equality, report the true TEA
// (doubles), and decide the cap EXACTLY over the pessimistic runs. The
// certificate logic is identical regardless of how `arcs` was produced -- the
// harvest only changes HOW rim arcs are found, never the decision. `arcs` is
// consumed (sorted in place). `r_sq` is FT(tool_radius)^2, the exact cutter
// squared radius already computed by the caller.
EngagementSample finish_engagement(std::vector<Arc>& arcs, double cx, double cy,
                                   const FT& r_sq, double cap_chord_ratio,
                                   double gap_close_ratio)
{
    EngagementSample out{0.0, 0.0, false};
    if (arcs.empty()) return out;   // rim nowhere in material

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
    //
    //     INVARIANT, checked because the reporting path had none. The rim is a
    //     closed curve of angular measure 2*pi and the runs are DISJOINT arcs of
    //     it, so neither a single run nor their sum can measure more than a full
    //     turn -- whatever the stock looks like. A larger reading is a logic error
    //     in the harvest or the assembly, never a tolerance question, so it is
    //     raised rather than clamped: clamping would ship the wrong number quietly,
    //     and this number is the audit's headline engagement metric.
    for (const Arc& run : runs) {
        if (run.span > FULL_TURN + FULL_TURN_REPORT_SLACK)
            throw std::logic_error("engagement run span " + format_double(run.span)
                                   + " rad exceeds a full turn at station ("
                                   + format_double(cx) + ", " + format_double(cy) + ").");
        out.total_tea += run.span;
        out.max_run_tea = std::max(out.max_run_tea, run.span);
    }
    if (out.total_tea > FULL_TURN + FULL_TURN_REPORT_SLACK)
        throw std::logic_error("engagement total_tea " + format_double(out.total_tea)
                               + " rad over " + std::to_string(runs.size())
                               + " disjoint runs exceeds a full turn at station ("
                               + format_double(cx) + ", " + format_double(cy) + ").");

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

// ----------------------------------------------------------------------------
// Local zone-query harvest: read the engaged rim arcs LOCALLY off the Gps's own
// arrangement (faces carry contained() = material), instead of overlaying the
// whole stock with the cutter disk. Validated end to end by the compiled
// prototypes in docs/dev/arrangement-redesign/ (see findings.md).
// ----------------------------------------------------------------------------

using Arr = Gps::Arrangement_2;
using PL = CGAL::Arr_walk_along_line_point_location<Arr>;

// Arrangement_zone_2 visitor (models CGAL/Arrangement_2/Arr_compute_zone_visitor.h,
// verified against the vendored CGAL 6.0.1 header): collect every cutter-circle
// sub-arc whose containing face is contained() -- those are exactly the rim
// pieces in material. Non-inserting (returns an invalid halfedge + "continue"),
// so the zone never mutates the arrangement.
struct EngagementVisitor {
    using X_monotone_curve_2 = Arr::X_monotone_curve_2;
    using Vertex_handle = Arr::Vertex_handle;
    using Halfedge_handle = Arr::Halfedge_handle;
    using Face_handle = Arr::Face_handle;
    using Result = std::pair<Halfedge_handle, bool>;

    std::vector<X_monotone_curve_2> engaged;

    void init(Arr*) {}
    Result found_subcurve(const X_monotone_curve_2& cv, Face_handle face,
                          Vertex_handle, Halfedge_handle, Vertex_handle, Halfedge_handle) {
        if (face->contained()) engaged.push_back(cv);   // cutter rim in material
        return Result(Halfedge_handle(), false);
    }
    Result found_overlap(const X_monotone_curve_2&, Halfedge_handle, Vertex_handle, Vertex_handle) {
        return Result(Halfedge_handle(), false);         // grazing boundary: measure-zero
    }
};

// Harvest the rim arcs of the cutter of radius `tool_radius` centred at (cx, cy)
// by zoning the cutter circle in the stock's own arrangement. Fills the SAME
// Arc{ccw_start, ccw_end, span} vector the earlier whole-stock overlay produced,
// off the SAME split points -- disk_polygon and make_x_monotone_2 both split at
// the x-extreme rational points (cx +/- r, cy), and the remaining splits are the
// same exact cutter/stock crossings -- so finish_engagement DECIDES identically.
// (The per-arc span normalisation is no longer the overlay's; see the span formula
// below. It is a reporting double, so it moves no verdict.)
//
// const_cast: engagement_at holds a const Stock2&, but Arrangement_zone_2 and
// the point-location structure want a non-const Arrangement_2 handle. The zone
// here is strictly READ-ONLY -- a non-inserting visitor, no arrangement mutation
// -- and the kernel is single-threaded, so treating the logically-const stock as
// non-const for this local read is justified and encapsulated here.
void engaged_arcs_zone(const Stock2& stock, double cx, double cy,
                       double tool_radius, std::vector<Arc>& arcs)
{
    Arr& arr = const_cast<Gps&>(stock.set()).arrangement();
    PL pl(arr);

    const EPoint center(cx, cy);
    const FT r_sq = FT(tool_radius) * FT(tool_radius);
    GpsTraits traits;
    // Same exact cutter circle disk_polygon builds: ECircle(center, FT(r)^2).
    const GpsTraits::Curve_2 cutter(ECircle(center, r_sq));

    std::vector<boost::variant<GpsPoint, GpsXCurve>> xmono;
    traits.make_x_monotone_2_object()(cutter, std::back_inserter(xmono));

    EngagementVisitor vis;
    for (const auto& piece : xmono) {
        if (const GpsXCurve* xc = boost::get<GpsXCurve>(&piece)) {
            CGAL::Arrangement_zone_2<Arr, EngagementVisitor> zone(arr, &vis);
            zone.init(*xc, pl);
            zone.compute_zone();
        }
    }

    // Extract Arc{ccw_start, ccw_end, span} from each engaged rim sub-arc:
    // CCW-normalize by orientation, drop tangent-touch degeneracies, report the
    // span as a REPORTING double that never feeds a decision.
    //
    // SPAN FORMULA. The span is the UNSIGNED angle between the two station-relative
    // radius vectors, atan2(|u x v|, u . v) -- never a difference of two atan2
    // headings. Three properties earn it that shape:
    //
    //   * ITS CODOMAIN IS [0, pi] BY CONSTRUCTION, because the first argument is
    //     non-negative. So it needs no wrap normalisation, and none can misfire.
    //     A heading difference does need one: every harvested arc is a sub-arc of
    //     an x-monotone piece, make_x_monotone_2 splits the cutter circle at its
    //     x-extremes (cx +/- r, cy), and the CCW-start of any lower-half arc sits
    //     at heading +pi while its CCW-end sits at a NEGATIVE heading -- so the
    //     honest difference is negative and has to be lifted by 2*pi.
    //     `if (span <= 0) span += 2*pi` did that lift, and could not tell a
    //     genuine wrap from a difference that merely ROUNDED to zero.
    //
    //   * THAT ROUNDING IS REACHABLE, and was the defect. When the stock boundary
    //     crosses the rim within an ulp of an x-extreme split point, the zone hands
    //     back a real sub-arc of ~1e-16 rad whose endpoints are exactly DISTINCT
    //     (so the s == t skip above does not fire) but whose to_double headings are
    //     the identical double. The difference was then exactly 0.0, the lift fired,
    //     and a sliver of material was reported as a FULL TURN of engagement --
    //     folded into the abutting run by the exact endpoint merge, so
    //     total_tea/max_run_tea came out one whole turn too large. The three
    //     witness stations in tests/test_growth_bound.py are that configuration.
    //
    //   * IT IS UNIFORMLY ACCURATE where a heading difference is not. Cancellation
    //     in the cross product costs at most a few ulps of |u||v|, i.e. a few ulps
    //     of angle ABSOLUTE, at any span -- including the two ends. Near 0 the
    //     cross product is the small quantity and atan2 resolves it directly; near
    //     pi the dot product carries the answer. Taking |cross| rather than trusting
    //     its sign is what makes the pi end safe: the exact semicircle gives
    //     ux*vy - uy*vx == -0.0, and atan2(-0.0, negative) is -pi where
    //     atan2(+0.0, negative) is +pi.
    //
    // The [0, pi] codomain is the CORRECT range, not a truncation of it: an
    // x-monotone circular arc never turns back in x, so it lies wholly in one
    // half-disk about the centre and its angular extent cannot exceed pi. Should a
    // future traits change break that, the reporting invariant in
    // finish_engagement is the loud failure.
    //
    // The radius vectors are divided by the tool radius before the products are
    // formed. atan2 is invariant under a positive common scaling of both arguments,
    // so this changes no answer, but it keeps every component O(1) and both products
    // in range for ANY radius the seam admits (require_positive_radius takes any
    // finite positive double, and tool_radius is an independent argument at the
    // binding -- an ordinary stock really can be queried with r = 1e-200). Formed
    // raw, the products of such a cutter underflow to zero and those of an
    // r = 1e+200 cutter overflow to infinity -- a scale regression against the
    // heading difference this replaces, which was scale-free.
    //
    // What NO normalisation rescues is the read-out itself: once r falls below an
    // ulp of the station coordinate, to_double lands both endpoints exactly on the
    // centre, u and v are both zero, and no formula can recover an angle from
    // that. Measured at (5, 5) buried in material, r <= 1e-100: this reports
    // total_tea = 0 where the truth is 2*pi. The heading difference reported 2*pi
    // PER ARC there, i.e. 4*pi -- the same defect in its purest form, since a
    // collapsed difference is exactly the zero its wrap lifted. Both readings are
    // wrong and neither moves a verdict: cap_exceeded is decided on the exact
    // one-root endpoints and stays correct (True buried, False clear) to r = 1e-300.
    // That is the deciding/reporting split doing precisely its job.
    for (const GpsXCurve& xc : vis.engaged) {
        GpsPoint s = xc.source();
        GpsPoint t = xc.target();
        if (s == t) continue;   // tangent-touch degeneracy: zero-measure contact
        if (xc.orientation() == CGAL::CLOCKWISE) std::swap(s, t);
        const double ux = (CGAL::to_double(s.x()) - cx) / tool_radius;
        const double uy = (CGAL::to_double(s.y()) - cy) / tool_radius;
        const double vx = (CGAL::to_double(t.x()) - cx) / tool_radius;
        const double vy = (CGAL::to_double(t.y()) - cy) / tool_radius;
        const double span = std::atan2(std::fabs(ux * vy - uy * vx), ux * vx + uy * vy);
        arcs.push_back({s, t, span});
    }
}

// ----------------------------------------------------------------------------
// Task 7b: swept-annulus geometry -- the certificate's interior guarantee.
//
// The helpers below turn "what can the material inside an annulus subtend at its
// centre" into a double. They are the only place in this file where a geometric
// QUANTITY (never a geometric TRUTH) is read out in doubles for the certificate;
// every one of them is written to over-estimate. The derivation lives at
// swept_run_bound.
// ----------------------------------------------------------------------------

// Relative slack applied wherever a double is compared, divided, or used to size
// exact geometry rather than merely reported: the half-spacing the annulus is
// built from (widened UP), the transfer term's numerator (UP) and denominator
// (DOWN), and -- scaled by r + travel -- both the tangency-membership margin and
// the margin at which two components count as touching. Every quantity it covers
// comes from a to_double or a hypot (<= 1 ulp, 1.1e-16 relative), so 1e-12 clears
// the round-off by four decades while costing a part in 1e12. SAFE FAILURE
// DIRECTION: each use enlarges the bound (a wider annulus, a bigger transfer
// angle, a more eager group, a tangency admitted rather than dropped), which can
// only force extra refinement or a conservative refusal -- never a false
// certificate.
constexpr double SWEPT_BOUND_REL_SLACK = 1e-12;

// Absolute angular inflation (radians), applied both to each direction span
// before the gaps between them are measured and to the assembled bound. Every
// direction is an atan2 of station-relative coordinates of magnitude at most
// r + travel, so it carries a few ulps of that -- under 1e-15 rad -- and the
// transfer term's asin adds as little again. 1e-9 rad = 5.7e-8 deg clears the
// total by six decades and costs a part in 1e9 of the cap. SAFE FAILURE DIRECTION
// as above: widening the spans can only shrink a gap, and a shrunken gap can only
// enlarge the extent.
constexpr double SWEPT_BOUND_ANGULAR_SLACK = 1e-9;

// Axis-aligned box over a component, in coordinates RELATIVE TO THE STATION. Used
// ONLY to decide which components must be bounded TOGETHER (see group_touching) --
// never to measure an angle. An earlier revision read the angular extent off this
// box and was measurably frame-dependent: a box straddles the station in x iff the
// component's direction range covers 90 or 270 deg and in y iff it covers 0 or
// 180 deg, so above 90 deg of span the answer turned on where the feature happened
// to sit relative to the world X axis. Measured: the same wall at h = 0.25 was
// bounded tightly at the four axis-aligned orientations and saturated to a full
// turn at 12 of 24, and the certifier refused a 120 deg pass against a 150 deg cap
// for no reason but the wall's angle. The physics is rotation-invariant; the bound
// now is too (component_direction_spans).
struct Aabb {
    double xmin = std::numeric_limits<double>::infinity();
    double ymin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();

    void add(double x, double y)
    {
        xmin = std::min(xmin, x);
        xmax = std::max(xmax, x);
        ymin = std::min(ymin, y);
        ymax = std::max(ymax, y);
    }

    // The two boxes come within `margin` of each other in BOTH axes.
    bool near(const Aabb& other, double margin) const
    {
        return xmin - margin <= other.xmax && other.xmin - margin <= xmax
               && ymin - margin <= other.ymax && other.ymin - margin <= ymax;
    }
};

// A closed arc of DIRECTIONS seen from the station: `extent` radians counter-
// clockwise from `start`. One is produced per boundary sub-curve; their union over
// a component is that component's direction range.
struct DirSpan {
    double start;
    double extent;
};

// Into (-pi, pi].
double wrap_to_pi(double a)
{
    a = std::fmod(a + std::numbers::pi, FULL_TURN);
    if (a <= 0.0) a += FULL_TURN;
    return a - std::numbers::pi;
}

// Counter-clockwise sweep from `from` to `to`, in [0, 2*pi).
double ccw_delta(double from, double to)
{
    const double d = std::fmod(to - from, FULL_TURN);
    return d < 0.0 ? d + FULL_TURN : d;
}

// Direction of an exact boundary point seen from the station, in (-pi, pi]. The
// station is subtracted EXACTLY first -- it is rational (root 0), so the
// Sqrt_extension same-root precondition holds trivially -- which keeps both
// operands O(r) and the direction accurate to a few ulps of the tool radius
// however far from the origin the stock sits.
double direction_from_station(const GpsPoint& p, const CoordNT& sx, const CoordNT& sy)
{
    return std::atan2(CGAL::to_double(p.y() - sy), CGAL::to_double(p.x() - sx));
}

// The direction range of ONE boundary sub-curve, appended to `spans`.
//
// Every point of a component's boundary lies in the closed annulus, so it is at
// least (r - travel) > 0 from the station: no sub-curve contains the station and
// no direction below is undefined.
//
// THREE CASES, selected by an EXACT predicate, because the direction genuinely
// behaves differently in each:
//
//   * SEGMENT. The direction is monotone along it, and the angle it subtends is
//     strictly below pi (the whole segment stays outside the disk of radius
//     r - travel, which caps the subtended angle at 2*acos((r-travel)/(r+travel))).
//     So the span is the SHORT arc between the endpoint directions, unambiguously.
//   * ARC whose supporting circle CONTAINS the station (|c - S| <= rho). The
//     direction then winds MONOTONICALLY with travel around that circle, by an
//     amount that can exceed pi -- the annulus's own two circles are exactly this
//     case, with c == S. What disambiguates the span is therefore the sub-curve's
//     ORIENTATION, not a shorter-arc rule.
//   * ARC whose supporting circle EXCLUDES the station. Every point of that circle
//     is seen within alpha = asin(rho / |c - S|) of the direction to its centre, so
//     the span lies in a cone of extent 2*alpha < pi. Along the arc the direction
//     is monotone except at the two TANGENCY points (where the line of sight
//     touches the circle), so the span is the smallest arc containing the endpoint
//     directions plus whichever tangency directions lie ON the arc.
//
// `|c - S| <= rho` is decided EXACTLY on FT: squared distance against squared
// radius, no tolerance. Tangency MEMBERSHIP is decided in doubles with an
// INCLUSIVE margin -- admitting a tangency that is not on the arc only widens the
// span, while dropping one that is would narrow it, so the bias is the safe one.
void append_curve_span(std::vector<DirSpan>& spans, const GpsXCurve& cv,
                       const CoordNT& sx, const CoordNT& sy, double margin_scale)
{
    if (cv.is_linear()) {
        const double a = direction_from_station(cv.left(), sx, sy);
        const double b = direction_from_station(cv.right(), sx, sy);
        const double delta = wrap_to_pi(b - a);
        spans.push_back(delta >= 0.0 ? DirSpan{a, delta} : DirSpan{b, -delta});
        return;
    }

    const ECircle circle = cv.supporting_circle();
    const FT dx_ft = circle.center().x() - sx.a0();
    const FT dy_ft = circle.center().y() - sy.a0();
    const double src = direction_from_station(cv.source(), sx, sy);
    const double dst = direction_from_station(cv.target(), sx, sy);

    if (CGAL::compare(dx_ft * dx_ft + dy_ft * dy_ft, circle.squared_radius()) != CGAL::LARGER) {
        const bool ccw = cv.orientation() == CGAL::COUNTERCLOCKWISE;
        const double from = ccw ? src : dst;
        const double to = ccw ? dst : src;
        spans.push_back({from, ccw_delta(from, to)});
        return;
    }

    const double cx = CGAL::to_double(dx_ft);
    const double cy = CGAL::to_double(dy_ft);
    const double rho = std::sqrt(CGAL::to_double(circle.squared_radius()));
    const double dist = std::hypot(cx, cy);
    const double axis = std::atan2(cy, cx);
    const double sin_alpha = std::min(1.0, rho / dist);

    double lo = wrap_to_pi(src - axis);
    double hi = lo;
    const double other = wrap_to_pi(dst - axis);
    lo = std::min(lo, other);
    hi = std::max(hi, other);

    // An x-monotone circle sub-arc lies wholly in the upper or lower half of its
    // supporting circle (CGAL's own _is_upper(), spelled from the two public
    // accessors), and spans the x-range of its two endpoints. A tangency point is
    // on the arc iff it satisfies both.
    const bool upper = (cv.orientation() == CGAL::COUNTERCLOCKWISE) != cv.is_directed_right();
    const double x_lo = CGAL::to_double(cv.left().x() - sx);
    const double x_hi = CGAL::to_double(cv.right().x() - sx);
    const double margin = SWEPT_BOUND_REL_SLACK * margin_scale;
    // Tangency points sit at circle-parameter axis +/- (pi - acos(rho/|c-S|)); the
    // matching line of sight is at axis +/- asin(rho/|c-S|).
    const double psi = std::numbers::pi - std::acos(sin_alpha);
    const double alpha = std::asin(sin_alpha);
    for (const double sign : {-1.0, 1.0}) {
        const double tx = cx + rho * std::cos(axis + sign * psi);
        const double ty = cy + rho * std::sin(axis + sign * psi);
        const bool on_half = upper ? (ty >= cy - margin) : (ty <= cy + margin);
        if (!on_half || tx < x_lo - margin || tx > x_hi + margin) continue;
        lo = std::min(lo, sign * alpha);
        hi = std::max(hi, sign * alpha);
    }
    spans.push_back({axis + lo, hi - lo});
}

// Angular extent of the UNION of a component group's direction spans: a full turn
// minus its largest uncovered gap.
//
// This is where the wrap case is settled numerically as well as topologically --
// spans that leave no gap mean the group reaches every direction from the station,
// and a full turn is the honest answer. It is also where the bound becomes
// rotation-invariant: nothing here refers to a coordinate axis.
//
// ROUNDING. Each span is widened by `slack` on both sides before the gaps are
// measured, so every gap is UNDER-estimated and the returned extent OVER-estimated
// -- the safe direction. A gap narrower than 2*slack is closed, which can only
// merge two genuinely-adjacent runs into one larger claim.
double union_extent(const std::vector<DirSpan>& spans, double slack)
{
    if (spans.empty()) return 0.0;

    std::vector<std::pair<double, double>> arcs;
    arcs.reserve(2 * spans.size());
    for (const DirSpan& span : spans) {
        const double extent = span.extent + 2.0 * slack;
        if (extent >= FULL_TURN) return FULL_TURN;
        double lo = span.start - slack;
        lo -= FULL_TURN * std::floor(lo / FULL_TURN);   // into [0, 2*pi)
        const double hi = lo + extent;
        if (hi <= FULL_TURN) {
            arcs.emplace_back(lo, hi);
        } else {
            arcs.emplace_back(lo, FULL_TURN);
            arcs.emplace_back(0.0, hi - FULL_TURN);
        }
    }

    std::sort(arcs.begin(), arcs.end());
    double covered_to = arcs.front().second;
    double largest_gap = 0.0;
    for (std::size_t i = 1; i < arcs.size(); ++i) {
        largest_gap = std::max(largest_gap, arcs[i].first - covered_to);
        covered_to = std::max(covered_to, arcs[i].second);
    }
    // The seam: from the far end of the covered set back round to the first arc.
    largest_gap = std::max(largest_gap, arcs.front().first + FULL_TURN - covered_to);
    return largest_gap <= 0.0 ? FULL_TURN : FULL_TURN - largest_gap;
}

// One connected component of the swept material, reduced to what the bound needs.
struct SweptComponent {
    Aabb box;                    // adjacency only -- never an angle
    std::vector<DirSpan> spans;  // its direction range, as a union of arcs
    bool encloses_station;       // EXACT: the station lies inside its outer boundary
};

// Read one component: its station-frame box, its boundary's direction spans, and
// the EXACT answer to whether it encloses the station.
//
// THE WRAP DECISION IS A POINT-IN-REGION QUERY, not a numeric proxy. A component
// lies inside the annulus, so it never contains the station; the station can only
// be inside its OUTER boundary by sitting in one of its holes. Only a component
// that HAS a hole can therefore enclose it, and for those the question is settled
// by the same exact `Gps::oriented_side` that `Stock2::contains` uses -- the
// annular rib, whose component is a ring around the station, lands here.
//
// Reading only the outer boundary is sufficient for the extent as well: with the
// station outside that boundary every ray from it that reaches the component
// crosses the outer boundary first, so the component's direction range is
// contained in the boundary's; and when the station is inside it, the wrap test
// above has already returned.
SweptComponent read_component(const GpsPolygonWithHoles& component, const CoordNT& sx,
                              const CoordNT& sy, double margin_scale)
{
    SweptComponent out;
    out.encloses_station = false;
    if (component.holes_begin() != component.holes_end()) {
        Gps outline;
        outline.insert(component.outer_boundary());
        out.encloses_station =
            outline.oriented_side(GpsPoint(sx.a0(), sy.a0())) == CGAL::ON_POSITIVE_SIDE;
        if (out.encloses_station) return out;
    }

    const GpsPolygon& outer = component.outer_boundary();
    for (auto it = outer.curves_begin(); it != outer.curves_end(); ++it) {
        append_curve_span(out.spans, *it, sx, sy, margin_scale);
        out.box.add(CGAL::to_double(it->left().x() - sx), CGAL::to_double(it->left().y() - sy));
        out.box.add(CGAL::to_double(it->right().x() - sx), CGAL::to_double(it->right().y() - sy));
    }
    return out;
}

// Label components so that any two whose boxes come within `margin` share a label,
// transitively (union-find).
//
// WHY THIS IS NOT OPTIONAL. Step (2) of the derivation needs the CONNECTED
// components of the swept material as a POINT SET, but polygons_with_holes
// decomposes by EDGE adjacency: two material lobes meeting at a single point --
// two exactly tangent subtraction disks leave exactly that -- come back as two
// polygons, while a cutter rim can pass straight through the pinch and hold ONE
// engaged run spanning both. Bounding the two lobes separately would then
// UNDER-estimate, the one direction this bound may never fail in. The same
// argument covers the regularization of the Gps boolean, which drops the
// lower-dimensional contact where a run leaves the annulus. Lobes that touch
// necessarily have overlapping boxes, so grouping on box proximity cannot miss
// such a pair.
//
// The converse -- grouping two genuinely separate components whose boxes happen to
// overlap -- only enlarges the bound, and it does not arise for the shape that
// matters: the two banks of a slot, or the two crossings of a rib, sit on opposite
// sides of the station and their boxes are disjoint. Measured: grouping changed no
// verdict anywhere in the suite.
std::vector<std::size_t> group_touching(const std::vector<SweptComponent>& parts, double margin)
{
    std::vector<std::size_t> parent(parts.size());
    for (std::size_t i = 0; i < parent.size(); ++i) parent[i] = i;
    const auto root = [&parent](std::size_t i) {
        while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i]; }
        return i;
    };
    for (std::size_t i = 0; i < parts.size(); ++i)
        for (std::size_t j = i + 1; j < parts.size(); ++j)
            if (parts[i].box.near(parts[j].box, margin)) parent[root(j)] = root(i);
    std::vector<std::size_t> label(parts.size());
    for (std::size_t i = 0; i < parts.size(); ++i) label[i] = root(i);
    return label;
}

// ----------------------------------------------------------------------------
// Task 5: exact-station TEA cap certificate along a linear cutter motion.
// ----------------------------------------------------------------------------

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

// THE INTERIOR CERTIFICATE. True only when no cutter centre on the span
// [(x0,y0),(x1,y1)] can hold an engaged run over `cap`.
//
// Every centre on the span lies within `half_spacing` of the NEARER endpoint, so
// swept_run_bound at both endpoints bounds the largest run at every centre. This
// is the whole soundness argument, and it is a statement about the region the
// cutter actually sweeps -- not about how a run measured at a station grows, the
// premise a rib wrapping the rim breaks (tests/test_false_certificate.py).
//
// half_spacing == 0 is not a degenerate case to guard against but a complete
// certificate on its own: the span IS the single centre (x0,y0) == (x1,y1), which
// the caller's EXACT station predicate decides outright, with no interior to
// bound. (swept_run_bound saturates to a full turn there -- its annulus has empty
// interior -- so it must not be consulted.)
bool interior_run_within_cap(const Stock2& stock, double x0, double y0, double x1,
                             double y1, double r, double half_spacing, double cap)
{
    if (!(half_spacing > 0.0)) return true;
    return swept_run_bound(stock, x0, y0, r, half_spacing) <= cap
           && swept_run_bound(stock, x1, y1, r, half_spacing) <= cap;
}

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
        // The stations are the CHEAP test (a local zone query each) and the swept
        // bound the expensive one (a boolean against the stock), so the stations
        // are asked first: the interior bound is only ever computed for a span the
        // pre-repair certifier would have passed. Both must hold -- (i) the
        // interior bound IS the proof, (ii) the guarded station test is the
        // retained conservative filter (see certify_segment_tea in the header).
        if (!e0.cap_exceeded && !e1.cap_exceeded
            && interior_run_within_cap(stock, x0, y0, x1, y1, r, half_spacing, cap))
            return;   // whole span certified
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

// Exported at file scope (declared in engagement_2.h) so the Python audit layer
// CALLS this guard rather than mirroring it: it has exactly one definition.
//
// RETIRED AS A SOUNDNESS ARGUMENT (Task 7b). This lemma is measurably NOT an
// upper bound on TEA growth: tests/test_growth_bound.py exhibits 18
// certificate-critical violations, by ratios up to 24x, on concave features
// (rho > r) and on the tool sitting in its own hole (rho == r). The false
// certificates that follow are committed as tests/test_false_certificate.py --
// an annular RIB of the tool's own radius, the same rib opened to a 135 deg
// sector (both stations then report NO contact at all, max_tea == 0.0, over stock
// the cutter is 135 deg engaged in), and a SPIRAL rib whose radius merely sweeps
// THROUGH the tool radius rather than matching it, falsely certified in 3 of 3
// directions at stations == 1.
//
// The failure is not looseness. Clauses (a) and (b) below bound how an EXISTING
// engaged run moves; a rib wrapping the rim creates engagement that grows from
// nothing either station can see, so the premise itself does not hold and no
// choice of safety factor repairs it. The certificate's interior guarantee is now
// swept_run_bound, which bounds what the swept REGION can contain and needs no
// such premise. This function and tea_guard survive as that record and as the
// certifier's retained station filter, which can only refuse.
//
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

// ----------------------------------------------------------------------------
// EXACT SWEPT-ANNULUS RUN BOUND -- the certificate's interior guarantee.
//
// Bounds the largest contiguous engaged run that ANY cutter centre C' within
// `travel` (hs below) of the station C = (cx, cy) can see. Where the retired
// growth lemma asked how an EXISTING run moves between two measurements, this
// asks what the region the cutter actually sweeps is even CAPABLE of holding --
// so it has nothing to be blind to.
//
//   (1) CONTAINMENT. For |C' - C| <= hs, every point x of the rim d B(C', r)
//       obeys r - hs <= |x - C| <= r + hs by the triangle inequality. So the rim
//       of EVERY reachable centre lies in the annulus A(C, r - hs, r + hs).
//   (2) CONNECTIVITY. A maximal engaged run of C' is a CONNECTED arc of
//       d B(C', r) lying in material, hence (with 1) a connected subset of
//       material  A -- so it lies inside a single CONNECTED COMPONENT K of it.
//       This is what keeps ordinary cutting certifiable: two banks of a slot are
//       separate components, and no run can span both.
//   (3) ANGULAR TRANSFER. For x on d B(C', r), the directions to x from C' and
//       from C differ by the angle at x in triangle (C, C', x). The side opposite
//       it is |C - C'| <= hs, the other two are r and >= r - hs (by 1), so with
//       hs < r/2 it is the strictly shortest side and its angle the strictly
//       smallest -- below pi/3, so certainly in [0, pi/2] where the law of sines
//       inverts: the deviation is at most asin(hs / (r - hs)). Applying it at
//       both ends of the run,
//
//         max_run(C')  <=  max_K ang(K)  +  2 * asin(hs / (r - hs))
//
//       where ang(K) is K's angular extent seen from C. (Formally: dir_C'(run)
//       lies inside the asin-neighbourhood of dir_C(run)  dir_C(K), and an arc
//       contained in an arc has no greater extent.)
//   (4) ang(K) FROM ITS BOUNDARY. With the station outside K's outer boundary,
//       every ray from the station that reaches K crosses that boundary first, so
//       dir(K) is contained in the boundary's own direction range. That range is
//       assembled sub-curve by sub-curve (append_curve_span, three cases selected
//       by an exact predicate) and unioned; ang(K) is a full turn minus the
//       union's largest gap (union_extent). Nothing in it refers to a coordinate
//       axis, so the bound is ROTATION-INVARIANT, as the geometry it measures is.
//       The WRAP case -- the station inside K's outer boundary, i.e. in one of its
//       holes, which is the annular rib -- is settled first and EXACTLY by
//       Gps::oriented_side (read_component).
//
// EXACTNESS (docs/exactness.md). Every geometric TRUTH consumed here is exact:
// the annulus is built from two exact circles about the exactly-injected centre,
// `material  annulus` is an exact boolean on the stock's own Gps, the connected
// components are that set's own polygons_with_holes decomposition, the station is
// subtracted from every boundary coordinate EXACTLY before it is read out, WHETHER
// A COMPONENT ENCLOSES THE STATION is an exact Gps::oriented_side query (the one
// topological fact in the construction), and the case split inside
// append_curve_span is an exact FT predicate. What is read out in doubles is a
// QUANTITY, not a truth -- directions and two angles -- and every read-out is
// inflated (SWEPT_BOUND_REL_SLACK, SWEPT_BOUND_ANGULAR_SLACK), by four to six
// decades over the accumulated round-off. This is the "analytic bounds are not
// precision handling" clause, and unlike tea_growth_bound this bound IS
// load-bearing, so the slack is not optional decoration: it is what makes "the
// computed number is an upper bound on the exact one" true rather than
// approximately true.
//
// TWO REPRESENTATIONS OF THE SAME CAP now coexist, deliberately. The station
// predicate enforces the cap through its exact rational chord surrogate
// 4*sin^2(cap/2); this bound is compared against cap in RADIANS
// (interior_run_within_cap). The two thresholds differ by ~1e-16 rad and the bound
// carries 1e-9 rad of inflation, so the interior gate is the stricter of the two --
// but do NOT "harmonise" them by relaxing the radian comparison toward the
// surrogate. The safe direction is to tighten the bound, never the threshold.
//
// SAFE FAILURE DIRECTION, stated once for all of it: every approximation here
// enlarges the returned bound. A bound too large forces extra refinement or a
// conservative "uncertified" verdict. A bound too small would be a false
// certificate -- which is precisely what the inflations forbid.
//
// COST. One exact boolean against the stock per call, versus the O(cutter-
// crossings) zone query a station costs. certify_recursive therefore asks the
// stations FIRST and reaches this only for a span that would otherwise certify.
// ----------------------------------------------------------------------------
double swept_run_bound(const Stock2& stock, double cx, double cy,
                       double tool_radius, double travel)
{
    // Station geometry, in signature order: the same seam contract engagement_at
    // enforces, plus `travel` -- a displacement magnitude, so negative is not a
    // shorter motion but a nonsensical one, and it is refused rather than folded.
    require_finite(cx, "cx");
    require_finite(cy, "cy");
    require_positive_radius(tool_radius, "tool_radius");
    require_finite(travel, "travel");
    if (!(travel >= 0.0))
        throw std::invalid_argument("travel must be non-negative (got " + format_double(travel) + ").");

    // travel == 0: the annulus has empty interior, so the construction cannot see
    // the rim at all and 0.0 would be a lie. travel >= r/2: step (3)'s transfer
    // term needs hs < r - hs, and step (4)'s R0 = r - hs must stay positive. Both
    // saturate -- the honest answer for a step this construction cannot bound.
    if (!(travel > 0.0)) return FULL_TURN;
    if (!(travel < 0.5 * tool_radius)) return FULL_TURN;

    // `travel` is itself a rounded quantity where it matters most -- certify_recursive
    // derives it from a hypot of two coordinate differences -- so it can land a
    // couple of ulps BELOW the half-spacing it stands for. Widen it ONCE, here, and
    // use the widened value for both the annulus and the transfer term. That keeps
    // the function's invariant ("every approximation enlarges the returned bound")
    // true of the CONTAINMENT step too, where a too-narrow annulus would hide
    // reachable material instead of over-reporting it. Everything downstream is
    // then exact in FT: the annulus radii are the exact rationals step (1) names,
    // not doubles re-rounded at the seam.
    const FT r_ft(tool_radius);
    const FT hs_ft = FT(travel) * FT(1.0 + SWEPT_BOUND_REL_SLACK);
    const EPoint centre(cx, cy);
    Gps annulus;
    annulus.insert(disk_polygon(centre, r_ft + hs_ft));
    Gps inner;
    inner.insert(disk_polygon(centre, r_ft - hs_ft));
    annulus.difference(inner);

    // The swept region's material, and its CONNECTED COMPONENTS: exactly what
    // polygons_with_holes decomposes an exact Gps into.
    Gps swept;
    swept.intersection(stock.set(), annulus);
    std::vector<GpsPolygonWithHoles> components;
    components.reserve(swept.number_of_polygons_with_holes());
    swept.polygons_with_holes(std::back_inserter(components));
    // No material in the annulus: by (1) no reachable centre's rim touches
    // anything, so the bound is exactly zero -- not the transfer term, which
    // measures the spread of a run that does not exist.
    if (components.empty()) return 0.0;

    // Step (4): the widest connected group, seen from the station. Components that
    // touch are grouped first -- polygons_with_holes splits at a point pinch that a
    // single engaged run can cross (group_touching).
    const CoordNT station_x{FT(cx)};
    const CoordNT station_y{FT(cy)};
    const double scale = tool_radius + travel;   // every station-frame coordinate is at most this
    std::vector<SweptComponent> parts;
    parts.reserve(components.size());
    for (const GpsPolygonWithHoles& component : components) {
        // A bounded stock intersected with a bounded annulus cannot produce an
        // unbounded component: this is an invariant breach, not a case, so it fails
        // loudly rather than returning the most conservative answer and hiding a
        // real bug behind a conservative verdict.
        if (component.is_unbounded())
            throw std::logic_error("swept_run_bound: material ^ annulus produced an unbounded component.");
        parts.push_back(read_component(component, station_x, station_y, scale));
        if (parts.back().encloses_station) return FULL_TURN;   // exact wrap
    }

    const std::vector<std::size_t> label = group_touching(parts, SWEPT_BOUND_REL_SLACK * scale);
    double widest = 0.0;
    for (std::size_t g = 0; g < parts.size(); ++g) {
        if (label[g] != g) continue;   // not a group representative
        std::vector<DirSpan> group;
        for (std::size_t i = 0; i < parts.size(); ++i)
            if (label[i] == g) group.insert(group.end(), parts[i].spans.begin(), parts[i].spans.end());
        widest = std::max(widest, union_extent(group, SWEPT_BOUND_ANGULAR_SLACK));
        if (widest >= FULL_TURN) return FULL_TURN;   // saturated: nothing tighter to learn
    }

    // Step (3)'s transfer term, 2*asin(hs / (r - hs)), on the SAME widened half-
    // spacing the annulus was built from. travel < r/2 was established above, so
    // the ratio is below 1 up to that widening; the numerator is inflated and the
    // denominator deflated, so the quotient is a certain upper bound on the exact
    // ratio, and min against 1.0 keeps asin in domain regardless.
    const double r_inner = CGAL::to_double(r_ft - hs_ft) * (1.0 - SWEPT_BOUND_REL_SLACK);
    const double hs_up = CGAL::to_double(hs_ft) * (1.0 + SWEPT_BOUND_REL_SLACK);
    const double transfer = 2.0 * std::asin(std::min(1.0, hs_up / r_inner));
    return std::min(FULL_TURN, widest + transfer + SWEPT_BOUND_ANGULAR_SLACK);
}

EngagementSample engagement_at(const Stock2& stock, double cx, double cy,
                               double tool_radius, double cap_chord_ratio,
                               double gap_close_ratio)
{
    // Station geometry, in signature order. The centre must be a real point and
    // the cutter must exist before either enters exact-land (see the boundary
    // guards above): FT(cx), FT(cy) and FT(tool_radius)^2 are all constructed
    // downstream of here.
    require_finite(cx, "cx");
    require_finite(cy, "cy");
    require_positive_radius(tool_radius, "tool_radius");

    // API-boundary contract: cap_chord_ratio = 4*sin^2(cap/2) with 0 < cap <= pi
    // lies in (0, 4]. Validate the raw double before exact injection (NaN fails).
    if (!(cap_chord_ratio > 0.0 && cap_chord_ratio <= 4.0))
        throw std::invalid_argument("cap_chord_ratio must be in (0, 4] (= 4*sin^2(cap/2), 0 < cap <= pi).");
    // gap_close_ratio = 4*sin^2(gamma/2) with 0 <= gamma <= pi lies in [0, 4];
    // 0 (no gap closed) is the default and the pre-pessimism semantics.
    if (!(gap_close_ratio >= 0.0 && gap_close_ratio <= 4.0))
        throw std::invalid_argument("gap_close_ratio must be in [0, 4] (= 4*sin^2(gamma/2), 0 <= gamma <= pi).");

    // Harvest the engaged rim arcs LOCALLY off the stock's OWN arrangement (its Gps
    // faces carry contained() = material): zone the cutter circle there and collect
    // the rim sub-arcs in material, then run the exact run-assembly + cap decision.
    // This replaces the earlier whole-stock overlay (region =
    // disk.intersection(stock.set())): the SAME exact certificate, read in
    // O(cutter-crossings) instead of O(stock) -- measured ~4-5x/query, more on a
    // heavily depleted pocket. The overlay was validated to reproduce this exactly:
    // 0 certificate mismatches over 7200 comparisons; the only residual was a
    // <= 1e-15 rad reporting-double representation artifact (algebraically-equal
    // one-root crossings whose representation-sensitive to_double rounds by ulps),
    // never the decision (docs/superpowers/state/engagement-zone-divergence.md).
    const FT r_sq = FT(tool_radius) * FT(tool_radius);
    std::vector<Arc> arcs;
    engaged_arcs_zone(stock, cx, cy, tool_radius, arcs);
    return finish_engagement(arcs, cx, cy, r_sq, cap_chord_ratio, gap_close_ratio);
}

CertifiedTea certify_segment_tea(const Stock2& stock, double x0, double y0,
                                 double x1, double y1, double tool_radius,
                                 double cap_radians)
{
    // Motion geometry, in signature order. Beyond the same station contract
    // engagement_at enforces, a non-physical radius here also broke the
    // REFINEMENT: the recursion's spacing floor is STATION_FLOOR_FRACTION * r,
    // and for r <= 0 or NaN no segment length can ever fall below it, so the
    // bisection descended to CERTIFY_MAX_DEPTH only to report a merely
    // "uncertified" verdict for a tool that cannot exist.
    require_finite(x0, "x0");
    require_finite(y0, "y0");
    require_finite(x1, "x1");
    require_finite(y1, "y1");
    require_positive_radius(tool_radius, "tool_radius");
    // Diagnostic refinement only: the range test below already rejects non-finite
    // cap_radians (NaN fails `> 0.0`, +Inf fails `<= pi`). Naming the actual
    // defect beats reporting NaN as "outside a range", and keeps every double
    // parameter's finiteness stated explicitly rather than implied by a
    // comparison a later edit could relax.
    require_finite(cap_radians, "cap_radians");

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

    // Refinement-bound accessors: the Python audit calls these rather than mirroring
    // them, so the guard has exactly one definition (docs/exactness.md). No exact
    // seam here -- these are pure double arithmetic that selects how far to shrink
    // the cap the exact station predicate then tests, never a geometric decision.
    m.def("tea_growth_bound", &tea_growth_bound, "d"_a, "r"_a);
    m.def("tea_guard", &tea_guard, "d"_a, "r"_a);

    // The certificate's interior guarantee, bound through a VALIDATING lambda
    // rather than handed to nanobind raw: cx, cy, tool_radius and travel all reach
    // exact injection (FT) or a domain assumption (travel >= 0, checked before the
    // annulus is built) inside swept_run_bound, and the sibling _sign_mixed_radical
    // binding records what a binding that trusts its caller returns -- a confident
    // answer for sqrt(-1). Exposed so the falsification harness measures the
    // SHIPPED bound rather than a Python re-derivation of it, exactly as
    // tea_growth_bound is.
    m.def("swept_run_bound",
          [](const Stock2& stock, double cx, double cy, double tool_radius, double travel) {
              require_finite(cx, "cx");
              require_finite(cy, "cy");
              require_positive_radius(tool_radius, "tool_radius");
              require_finite(travel, "travel");
              return swept_run_bound(stock, cx, cy, tool_radius, travel);
          },
          "stock"_a, "cx"_a, "cy"_a, "tool_radius"_a, "travel"_a);

    // Test-only: exact sign of A + B*sqrt(alpha) + C*sqrt(beta) + D*sqrt(alpha*beta)
    // (returns -1/0/+1) so the cap predicate's core primitive is unit-tested
    // directly against high-precision references.
    m.def("_sign_mixed_radical",
          [](double a, double b, double c, double d, double alpha, double beta) {
              // Same seam, same order as the signature: six doubles, each injected
              // as FT(x). The roots additionally carry the primitive's documented
              // precondition alpha, beta >= 0 -- unenforced, it returned a
              // confident sign for sqrt(-1) (_sign_mixed_radical(0,1,0,0,-1,0) -> 1).
              require_finite(a, "a");
              require_finite(b, "b");
              require_finite(c, "c");
              require_finite(d, "d");
              require_radicand(alpha, "alpha");
              require_radicand(beta, "beta");
              return static_cast<int>(
                  sign_mixed_radical(FT(a), FT(b), FT(c), FT(d), FT(alpha), FT(beta)));
          },
          "a"_a, "b"_a, "c"_a, "d"_a, "alpha"_a, "beta"_a);
}
