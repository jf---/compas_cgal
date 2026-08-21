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
// Arc{ccw_start, ccw_end, span} vector the overlay harvest produces, with the
// IDENTICAL per-arc normalization, so finish_engagement yields the identical
// result (the split points coincide: disk_polygon and make_x_monotone_2 both
// split at the x-extreme rational points (cx +/- r, cy), the remaining splits
// are the same exact cutter/stock crossings).
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

    // Extract Arc{ccw_start, ccw_end, span} from each engaged rim sub-arc, with
    // the SAME normalization as the overlay harvest (docs: engagement_at):
    // CCW-normalize by orientation, drop tangent-touch degeneracies, report the
    // span as a REPORTING double (atan2) that never feeds a decision.
    for (const GpsXCurve& xc : vis.engaged) {
        GpsPoint s = xc.source();
        GpsPoint t = xc.target();
        if (s == t) continue;   // tangent-touch degeneracy: zero-measure contact
        if (xc.orientation() == CGAL::CLOCKWISE) std::swap(s, t);
        const double sx = CGAL::to_double(s.x()), sy = CGAL::to_double(s.y());
        const double tx = CGAL::to_double(t.x()), ty = CGAL::to_double(t.y());
        double span = std::atan2(ty - cy, tx - cx) - std::atan2(sy - cy, sx - cx);
        if (span <= 0.0) span += 2.0 * std::numbers::pi;
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

// A full turn of engaged rim (radians). Also the bound's SATURATED value: with
// the cap contractually in (0, pi], returning this always forces refinement or a
// conservative refusal, so it is the safe answer whenever the construction
// cannot bound the swept region below a full turn.
constexpr double FULL_TURN = 2.0 * std::numbers::pi;

// Relative slack applied wherever a double is compared or divided rather than
// merely reported: the transfer term's travel (inflated UP) and annulus inner
// radius (deflated DOWN), and -- scaled by r + travel -- the margin at which two
// component boxes count as touching. Every quantity it covers comes from a
// to_double (<= 1 ulp, 1.1e-16 relative), so 1e-12 clears the round-off by four
// decades while costing a part in 1e12. SAFE FAILURE DIRECTION: each use enlarges
// the bound (a bigger transfer angle, a more eager fuse), which can only force
// extra refinement or a conservative refusal -- never a false certificate.
constexpr double SWEPT_BOUND_REL_SLACK = 1e-12;

// Absolute angular inflation (radians) added to the assembled bound. The
// component extent is read from atan2 of station-relative coordinates whose
// magnitude is at most r + travel, so each corner direction carries at most a few
// ulps of that -- under 1e-15 rad -- and the transfer term's asin adds as little
// again. 1e-9 rad = 5.7e-8 deg clears the total by six decades and costs a part
// in 1e9 of the cap. SAFE FAILURE DIRECTION as above.
constexpr double SWEPT_BOUND_ANGULAR_SLACK = 1e-9;

// Axis-aligned box accumulated over a component's outer boundary, in coordinates
// RELATIVE TO THE STATION. Station-relative is not cosmetic: the box is read out
// in doubles and then turned into angles, so keeping every coordinate O(r) rather
// than O(world position) keeps a corner direction accurate to a few ulps of the
// tool radius no matter where in the model the stock sits.
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

    // Nothing was ever added -- an outer boundary with no sub-curves, which is how
    // an unbounded component presents itself. Never a real box.
    bool empty() const { return xmin > xmax; }

    // The station (the origin of this frame) lies in the closed box. Spelled with
    // <= on both sides so a coordinate that lands exactly on the station -- a
    // cutter riding along a straight wall does exactly that -- saturates rather
    // than being read as clearance.
    bool holds_origin() const
    {
        return xmin <= 0.0 && 0.0 <= xmax && ymin <= 0.0 && 0.0 <= ymax;
    }

    // The two boxes come within `margin` of each other in BOTH axes -- the test
    // that decides whether two components must be bounded together (see
    // merge_touching_boxes).
    bool near(const Aabb& other, double margin) const
    {
        return xmin - margin <= other.xmax && other.xmin - margin <= xmax
               && ymin - margin <= other.ymax && other.ymin - margin <= ymax;
    }

    void absorb(const Aabb& other)
    {
        add(other.xmin, other.ymin);
        add(other.xmax, other.ymax);
    }
};

// Shortest angular distance between two directions, in [0, pi].
double angular_distance(double u, double v)
{
    const double d = std::fabs(u - v);
    return d > std::numbers::pi ? FULL_TURN - d : d;
}

// Grow `box` to contain one boundary sub-curve, in the frame centred on
// (station_x, station_y). The two exact coordinates are subtracted EXACTLY before
// any to_double: the station is rational (root 0), so the CoordNT subtraction
// meets the Sqrt_extension same-root precondition trivially and loses nothing.
//
// A TIGHT arc box, deliberately not CGAL's own X_monotone_curve_2::bbox(): that
// one extends an upper arc's y_max to the whole supporting circle's top whether
// or not the top is on the arc, which for a stock arc of the tool's own radius
// inflates a thin annular sliver's box by up to r -- enough to swamp the bound.
// Here the extremum is admitted only when it is actually on the arc, and THAT is
// decided EXACTLY: an x-monotone circle sub-arc lies wholly in one half of its
// supporting circle, so its only y-extremum beyond the endpoints is that circle's
// top (upper arc) or bottom (lower arc), which lies on the arc iff the arc's
// x-range spans the circle's centre -- one exact CoordNT comparison per side.
void grow_box_with_curve(Aabb& box, const GpsXCurve& cv, const CoordNT& station_x,
                         const CoordNT& station_y)
{
    for (const GpsPoint* pt : {&cv.left(), &cv.right()})
        box.add(CGAL::to_double(pt->x() - station_x), CGAL::to_double(pt->y() - station_y));
    if (!cv.is_circular()) return;   // a segment's extrema ARE its endpoints

    const ECircle circle = cv.supporting_circle();
    const CoordNT ox{circle.center().x()};
    if (CGAL::compare(cv.left().x(), ox) == CGAL::LARGER) return;    // arc left of centre
    if (CGAL::compare(cv.right().x(), ox) == CGAL::SMALLER) return;  // arc right of centre

    // Same predicate CGAL's own _is_upper() uses, spelled from the two public
    // accessors: CCW travelling leftwards, or CW travelling rightwards, is the
    // upper half.
    const bool upper = (cv.orientation() == CGAL::COUNTERCLOCKWISE) != cv.is_directed_right();
    const FT dx = circle.center().x() - station_x.a0();
    const FT dy = circle.center().y() - station_y.a0();
    const double rho = std::sqrt(CGAL::to_double(circle.squared_radius()));
    const double dyd = CGAL::to_double(dy);
    box.add(CGAL::to_double(dx), upper ? dyd + rho : dyd - rho);
}

// Upper bound on the angular extent, seen from the station, of everything inside
// one station-frame box -- step (4) of the derivation at swept_run_bound.
//
// A box is convex, so if the station is outside it the directions from the station
// to the box form an arc spanned by its four CORNERS, and an extent is monotone
// under inclusion. Hence extent(anything inside the box) <= that corner spread.
//
// THE WRAP CASE, decided by an inclusion rather than by topology. A component
// that wraps the station reaches every direction from it, so it holds points on
// opposite sides in both axes and its box MUST contain the station. Saturating
// whenever the box holds the station therefore cannot miss a wrap -- and the
// converse over-estimate it admits (a non-wrapping component whose box straddles
// the station in both axes) costs nothing: such a component spans nearly a half
// turn or more, which no cap in (0, pi] can accept once the transfer term is
// added anyway. An empty box means the outer boundary carried no sub-curves --
// an unbounded component -- and saturates for the same reason: nothing below a
// full turn has been established.
double box_angular_bound(const Aabb& box)
{
    if (box.empty() || box.holds_origin()) return FULL_TURN;

    const double corner[4] = {
        std::atan2(box.ymin, box.xmin), std::atan2(box.ymin, box.xmax),
        std::atan2(box.ymax, box.xmin), std::atan2(box.ymax, box.xmax),
    };
    double widest = 0.0;
    for (int i = 0; i < 4; ++i)
        for (int j = i + 1; j < 4; ++j)
            widest = std::max(widest, angular_distance(corner[i], corner[j]));
    return widest;
}

// Fuse every group of boxes that come within `margin` of one another, to a
// fixpoint, and report how many independent boxes remain (compacted to the front
// of `boxes`).
//
// WHY THIS IS NOT OPTIONAL. Step (2) of the derivation needs the CONNECTED
// components of the swept material as a POINT SET, but polygons_with_holes
// decomposes by EDGE adjacency: two material lobes meeting at a single point --
// two exactly tangent subtraction disks leave exactly that -- come back as two
// polygons, while a cutter rim can pass straight through the pinch and hold ONE
// engaged run spanning both. Taking the max over the two lobes separately would
// then UNDER-estimate, which is the one direction this bound may never fail in.
// Lobes that touch necessarily have overlapping boxes, so fusing on box proximity
// cannot miss such a pair.
//
// The converse -- fusing two genuinely separate components whose boxes happen to
// overlap -- only enlarges the bound, and it does not arise for the shape that
// matters: the two banks of a slot, or the two crossings of a rib, sit on
// opposite sides of the station and their boxes are disjoint. Measured: fusing
// changed no verdict anywhere in the suite.
//
// `margin` absorbs the read-out round-off so an exact touch is never missed by an
// ulp; it is a few decades above the coordinate round-off and its direction is
// safe, since a spurious fuse can only over-estimate.
std::size_t merge_touching_boxes(std::vector<Aabb>& boxes, double margin)
{
    std::size_t count = boxes.size();
    for (std::size_t i = 0; i < count; ++i) {
        for (std::size_t j = i + 1; j < count;) {
            if (boxes[i].near(boxes[j], margin)) {
                boxes[i].absorb(boxes[j]);
                boxes[j] = boxes[--count];
                j = i + 1;   // the grown box may now reach boxes already passed
            } else {
                ++j;
            }
        }
    }
    return count;
}

// Station-frame box of one connected component, read off its OUTER boundary
// alone: the component (holes and all) lies inside that boundary.
Aabb component_box(const GpsPolygon& outer, const CoordNT& station_x,
                   const CoordNT& station_y)
{
    Aabb box;
    for (auto it = outer.curves_begin(); it != outer.curves_end(); ++it)
        grow_box_with_curve(box, *it, station_x, station_y);
    return box;
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
//   (4) ang(K) FROM AN INCLUSION. K lies inside the axis-aligned box of its own
//       outer boundary, taken in the station's frame, and extent is monotone
//       under inclusion -- so the box's corner spread bounds ang(K)
//       (component_angular_bound). The WRAP case falls out of the same inclusion:
//       a component wrapping the station reaches every direction and so straddles
//       it in both axes, which is exactly when the box test saturates. The
//       annular rib (a ring around the station) and the 4.8-rad spiral rib both
//       land there; a component clear of the station never does.
//
// EXACTNESS (docs/exactness.md). Every geometric TRUTH consumed here is exact:
// the annulus is built from two exact circles about the exactly-injected centre,
// `material  annulus` is an exact boolean on the stock's own Gps, the connected
// components are that set's own polygons_with_holes decomposition, the station is
// subtracted from every boundary coordinate EXACTLY before it is read out, and
// the comparisons inside grow_box_with_curve are exact CoordNT predicates. What
// is read out in doubles is a QUANTITY, not a truth -- a box and two angles -- and
// every read-out is inflated (SWEPT_BOUND_REL_SLACK, SWEPT_BOUND_ANGULAR_SLACK),
// by six decades over the accumulated round-off. This is the "analytic bounds are
// not precision handling" clause, and unlike tea_growth_bound this bound IS
// load-bearing, so the slack is not optional decoration: it is what makes "the
// computed number is an upper bound on the exact one" true rather than
// approximately true.
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

    // Exact annulus about the exactly-injected station. r +/- travel are formed in
    // FT, so the annulus radii are the exact rationals the containment argument
    // names -- not doubles re-rounded at the seam.
    const FT r_ft(tool_radius);
    const FT hs_ft(travel);
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
    // touch are fused first -- polygons_with_holes splits at a point pinch that a
    // single engaged run can cross (merge_touching_boxes).
    const CoordNT station_x{FT(cx)};
    const CoordNT station_y{FT(cy)};
    std::vector<Aabb> boxes;
    boxes.reserve(components.size());
    for (const GpsPolygonWithHoles& component : components) {
        // An unbounded component cannot be bounded by a box; it also cannot occur
        // for a bounded stock intersected with a bounded annulus, so it is a
        // contract breach rather than a case -- saturate rather than read a box
        // that does not describe it.
        if (component.is_unbounded()) return FULL_TURN;
        boxes.push_back(component_box(component.outer_boundary(), station_x, station_y));
    }
    const double touch_margin = SWEPT_BOUND_REL_SLACK * (tool_radius + travel);
    const std::size_t groups = merge_touching_boxes(boxes, touch_margin);

    double widest = 0.0;
    for (std::size_t i = 0; i < groups; ++i) {
        widest = std::max(widest, box_angular_bound(boxes[i]));
        if (widest >= FULL_TURN) return FULL_TURN;   // saturated: nothing tighter to learn
    }

    // Step (3)'s transfer term, 2*asin(hs / (r - hs)). travel < r/2 was
    // established above, so the ratio is below 1 and the asin is real. The
    // numerator is inflated and the denominator deflated, so the quotient is a
    // certain upper bound on the exact ratio; min against 1.0 keeps asin in
    // domain regardless.
    const double r_inner = CGAL::to_double(r_ft - hs_ft) * (1.0 - SWEPT_BOUND_REL_SLACK);
    const double transfer = 2.0 * std::asin(std::min(1.0, travel * (1.0 + SWEPT_BOUND_REL_SLACK) / r_inner));
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
