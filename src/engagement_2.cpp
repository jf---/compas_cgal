#include "engagement_2.h"
#include "audit_certification_2.h"
#include "audit_policy_2.h"
#include "stock_2.h"
#include "exact/one_root.h"
#include "exact_cap_certificate_2.h"
#include "audit_arc_motion_2.h"
#include "exact_circle_chart_2.h"
#include "swept_prefix_engagement_2.h"

#include <algorithm>
#include <cmath>
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
using CoordNT = GpsPoint::CoordNT;

// One rim arc, CCW-normalized: the CCW sweep runs ccw_start -> ccw_end on the
// cutter circle. `span` is the arc's angular extent in radians, a REPORTING
// double that never feeds a decision.
struct Arc {
    GpsPoint ccw_start;
    GpsPoint ccw_end;
    double span;
};


// One full turn, the reported span of a rim entirely buried in material.
constexpr double FULL_TURN = 2.0 * std::numbers::pi;

// Reporting-only slack on the "a run cannot exceed a full turn" invariant. The
// shared exact station classifier returns a partition of the cutter rim, so its
// reporting spans sum to at most FULL_TURN exactly; the only error here is the
// rounding of each rim_sub_arc_span evaluation and their summation, a handful of
// ulps of 2*pi (ulp(2*pi) ~ 8.9e-16). 1e-12 rad (~6e-11 deg) is ~1000x that and
// still 12 orders below the whole-turn normalization failure it detects.
//
// NOT A DECISION TOLERANCE (docs/exactness.md, deciding/reporting split). It never
// appears in the cap verdict, which run_exceeds_cap decides exactly on one-root
// endpoints; it bounds the rounding of a human-facing number only.
constexpr double FULL_TURN_REPORTING_SLACK = 1e-12;

// Angular span of ONE x-monotone rim sub-arc, in radians -- a REPORTING double
// that never feeds a decision.
//
// The cutter circle enters make_x_monotone_2 split at its two x-extreme points and
// the zone only subdivides further, so EVERY sub-arc reaching here subtends at
// most a half turn. That bound is what makes the endpoints alone sufficient: with
// a, b the endpoints taken relative to the centre (|a| == |b| == r), the chord
// half-length is r*sin(theta/2) and the chord midpoint's distance to the centre is
// r*cos(theta/2), so
//
//     theta = 2 * atan2(|b - a|, |a + b|)                             (Kahan)
//
// Both arguments are lengths, hence non-negative, so atan2 lands in [0, pi/2] and
// theta lands in [0, pi] BY CONSTRUCTION: there is no branch cut to normalize
// across, and a whole-turn promotion is unrepresentable rather than merely
// unlikely. The formula is also uniformly well conditioned -- theta -> 0 shrinks
// the first argument, theta -> pi shrinks the second, and neither limit is a
// subtractive cancellation.
//
// WHAT THIS REPLACES. The span used to be a difference of two ABSOLUTE atan2
// angles, wrapped by "+= 2*pi when non-positive". That wrap is correct for a
// genuine half-turn sub-arc, whose endpoints straddle the atan2 branch cut at the
// -x extreme and whose difference is exactly -pi. But it could not tell that case
// apart from a sub-arc of true span ~1e-16 whose angle difference merely ROUNDED
// non-positive -- and it promoted the latter by a whole turn, inflating the
// reported max_run_tea/total_tea of every run the degenerate arc merged into.
double rim_sub_arc_span(const GpsPoint& s, const GpsPoint& t, double cx, double cy)
{
    const double ax = CGAL::to_double(s.x()) - cx, ay = CGAL::to_double(s.y()) - cy;
    const double bx = CGAL::to_double(t.x()) - cx, by = CGAL::to_double(t.y()) - cy;
    return 2.0 * std::atan2(std::hypot(bx - ax, by - ay), std::hypot(ax + bx, ay + by));
}
// whose angular span does NOT exceed the gap-closure angle gamma (surrogate
// gap_close_ratio, exact threshold T_gamma = gap_close_ratio * r^2), returning
// the resulting PESSIMISTIC runs as (ccw_start, ccw_end) endpoint pairs for the
// cap DECISION. `runs` are the maximal true runs, CCW-sorted, none sharing an
// endpoint (abutting runs were already merged), so between run[i] and
// run[(i+1)%n] lies exactly one void gap -- the CCW arc from run[i].ccw_end to
// run[(i+1)%n].ccw_start (for n == 1 the single gap is the run's complement).
//
// EXACTNESS. A gap is itself an arc with exact one-root endpoints; it is absorbed
// iff its span does not exceed gamma, decided by the same
// compas_cgal::cap::run_exceeds_cap the run cap uses -- literally the one
// definition, not a second copy of the same machinery -- so no angle is ever
// summed. Absorbing gap[i] merges run[i] and run[(i+1)%n]; a chain
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
        gap_closed[i] = !compas_cgal::cap::run_exceeds_cap(g_start, g_end, cx, cy, gap_close_ratio, T_gamma);
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

    // INVARIANT. The engaged sub-arcs partition the cutter circle, so the runs
    // assembled from them cover it at most once: total_tea <= FULL_TURN, and
    // max_run_tea <= total_tea, so this one comparison covers both. Checked HERE
    // because this is the single choke point every reported engagement number
    // passes through, and because a mis-normalized sub-arc span is invisible until
    // it lands in one of these two fields. Costs one double compare per station.
    if (out.total_tea > FULL_TURN + FULL_TURN_REPORTING_SLACK)
        throw RimSpanNormalizationError(
            "Engaged rim runs sum to " + std::to_string(out.total_tea) +
            " rad, beyond the full turn the cutter circle can offer: a rim sub-arc "
            "span was normalized wrong (" + std::to_string(runs.size()) + " runs).");

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
        if (compas_cgal::cap::run_exceeds_cap(start, end, cxf, cyf, ratio_ft, T)) {
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
    // CCW-normalize by orientation, drop tangent-touch degeneracies, and report
    // the span as a REPORTING double that never feeds a decision.
    for (const GpsXCurve& xc : vis.engaged) {
        GpsPoint s = xc.source();
        GpsPoint t = xc.target();
        if (s == t) continue;   // tangent-touch degeneracy: zero-measure contact
        if (xc.orientation() == CGAL::CLOCKWISE) std::swap(s, t);
        arcs.push_back({s, t, rim_sub_arc_span(s, t, cx, cy)});
    }
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
// proof-level slack (CONVENTIONS.md analytic-bounds clause). SAFE FAILURE DIRECTION:
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
// whole motion. full_cap_chord_ratio is the precomputed FULL-cap surrogate,
// reused for max_tea reporting on segments too coarse to guard.
void certify_recursive(const Stock2& stock, double x0, double y0, double x1,
                       double y1, double r, double cap, double full_cap_chord_ratio,
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
        const double guarded_ratio = cap_chord_ratio(cap_guarded);
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
        const double closed_gap_cap = std::min(guard, std::numbers::pi);
        const double gap_close_ratio =
            closed_gap_cap == 0.0 ? 0.0 : cap_chord_ratio(closed_gap_cap);
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
        const EngagementSample e0 = engagement_at(stock, x0, y0, r, full_cap_chord_ratio);
        const EngagementSample e1 = engagement_at(stock, x1, y1, r, full_cap_chord_ratio);
        acc.max_tea = std::max({acc.max_tea, e0.max_run_tea, e1.max_run_tea});
    }

    // Not certifiable at this spacing: bisect, unless the floor/depth bound is
    // hit -- then the margin is unclosable and the motion is uncertified.
    if (seg_len < STATION_FLOOR_FRACTION * r || depth >= CERTIFY_MAX_DEPTH) {
        acc.cap_certified = false;
        return;
    }
    const double mx = 0.5 * (x0 + x1), my = 0.5 * (y0 + y1);
    certify_recursive(stock, x0, y0, mx, my, r, cap, full_cap_chord_ratio, acc, depth + 1);
    if (!acc.cap_certified) return;
    certify_recursive(stock, mx, my, x1, y1, r, cap, full_cap_chord_ratio, acc, depth + 1);
}

// Reporting view of an exact station classification: the verdict is copied,
// the spans are measured about the station's double centre (REPORTING only).
EngagementSample engagement_sample_from_exact(
    const AuditExactStationClassification2& exact,
    double center_x,
    double center_y)
{
    EngagementSample out{
        0.0,
        0.0,
        exact.disposition()
            == AuditExactStationDisposition2::CAP_EXCEEDED,
    };
    for (const auto& run : exact.true_run_arcs()) {
        double run_span = 0.0;
        for (const auto& [start, end] : run) {
            run_span += rim_sub_arc_span(start, end, center_x, center_y);
        }
        out.total_tea += run_span;
        out.max_run_tea = std::max(out.max_run_tea, run_span);
    }
    if (out.total_tea > FULL_TURN + FULL_TURN_REPORTING_SLACK) {
        throw RimSpanNormalizationError(
            "Engaged rim runs exceed the cutter's full turn");
    }
    return out;
}

// Boundary validation of the prefix parameters: comparisons of the input
// doubles themselves (each IS a rational), not a tolerance.
void validate_prefix_parameters(const std::vector<double>& parameters)
{
    double previous = 0.0;
    for (const double parameter : parameters) {
        if (!std::isfinite(parameter)
            || parameter < 0.0
            || parameter > 1.0
            || parameter < previous) {
            throw InvalidSweptPrefixInputError(
                "swept-prefix parameters must be finite, non-decreasing, and lie in [0, 1]");
        }
        previous = parameter;
    }
}

// An exact rational unit vector along the direction (x, y) != 0, by the
// half-angle parametrization ((1 - t^2), 2t) / (1 + t^2), whose squared length is
// 1 identically. t = y / (|p| + x) is a CONSTRUCTION parameter computed in
// doubles (placement intent, never a decision): any rational t gives an exact
// unit vector, and the double only steers its direction. Reflecting x < 0 into
// x >= 0 keeps |t| <= 1, away from the parametrization's pole.
EVector rational_unit_along(double x, double y)
{
    const bool reflected = x < 0.0;
    const double px = reflected ? -x : x;
    const double py = reflected ? -y : y;
    const FT t(py / (std::hypot(px, py) + px));
    const FT t_squared = t * t;
    const FT denominator = FT(1) + t_squared;
    const EVector unit((FT(1) - t_squared) / denominator, FT(2) * t / denominator);
    return reflected ? -unit : unit;
}

// One circular motion from `start` about `guide`, its stations classified
// against the stock minus the motion's own swept prefix. Parameter 0 (empty
// prefix) is classified at `zero_station`: where the tool physically is.
// `total_sweep` is the travel INTENT in radians that places the interior
// stations; every station is `start` rotated exactly by the quarter-chart
// atlas, so all lie on start's guide circle.
struct ArcPrefixMotion2 {
    EPoint guide;
    EPoint zero_station;
    EPoint start;
    EPoint end;
    EVector phase;  // start - guide, the vector every interior station rotates
    bool full_turn;
    double total_sweep;
    bool clockwise;
};

// One station of a swept-prefix motion: its exact centre and prefix extent,
// placed BEFORE any classification so that a decision-only scan raises on
// exactly the placement errors the full prefix call raises on.
struct ArcPrefixStation2 {
    EPoint centre;
    ArcPrefixExtent2 extent;
};

ArcPrefixStation2 arc_prefix_station(
    const ArcPrefixMotion2& motion, double parameter)
{
    if (parameter == 1.0) {
        return {
            motion.end,
            motion.full_turn
                ? ArcPrefixExtent2::FULL_TURN
                : classify_arc_prefix_extent_exact(
                      motion.guide, motion.start, motion.end, motion.clockwise),
        };
    }
    if (parameter == 0.0) {
        return {motion.zero_station, ArcPrefixExtent2::ZERO};
    }
    const double turned = motion.total_sweep * parameter;
    const EPoint centre = exact_circle_chart_point(
        motion.guide,
        motion.phase,
        audit_arc_chart_parameter_for_angle(
            motion.clockwise ? -turned : turned));
    return {
        centre,
        classify_arc_prefix_extent_exact(
            motion.guide, motion.start, centre, motion.clockwise),
    };
}

AuditExactStationClassification2 classify_arc_prefix_station(
    const Stock2& stock, const ArcPrefixMotion2& motion,
    const ArcPrefixStation2& station, const FT& tool_radius,
    const FT& cap_chord_ratio, const FT& gap_close_ratio)
{
    return classify_arc_swept_prefix_station_exact(
        stock,
        motion.guide,
        motion.start,
        station.centre,
        motion.clockwise,
        station.extent,
        tool_radius,
        cap_chord_ratio,
        gap_close_ratio);
}

std::vector<EngagementSample> arc_prefix_samples(
    const Stock2& stock, const ArcPrefixMotion2& motion,
    const std::vector<double>& parameters, double tool_radius,
    double cap_chord_ratio, double gap_close_ratio)
{
    std::vector<EngagementSample> result;
    result.reserve(parameters.size());
    for (const double parameter : parameters) {
        const ArcPrefixStation2 station = arc_prefix_station(motion, parameter);
        result.push_back(engagement_sample_from_exact(
            classify_arc_prefix_station(
                stock, motion, station, FT(tool_radius),
                FT(cap_chord_ratio), FT(gap_close_ratio)),
            CGAL::to_double(station.centre.x()),
            CGAL::to_double(station.centre.y())));
    }
    return result;
}

// The DECISION-ONLY probing order over n stations: the last station first (a
// refused motion's end is where it bites: the profile of the gate cases found
// a refused bridge's end station exceeded in every case), then the rest in
// order. Any order yields the same verdict: "some station exceeds" is a
// disjunction of the same exact per-station verdicts.
std::vector<std::size_t> decision_probe_order(std::size_t count)
{
    std::vector<std::size_t> order;
    order.reserve(count);
    if (count == 0) return order;
    order.push_back(count - 1);
    for (std::size_t index = 0; index + 1 < count; ++index) {
        order.push_back(index);
    }
    return order;
}

bool arc_prefix_exceeded(
    const Stock2& stock, const ArcPrefixMotion2& motion,
    const std::vector<double>& parameters, double tool_radius,
    double cap_chord_ratio, double gap_close_ratio)
{
    std::vector<ArcPrefixStation2> stations;
    stations.reserve(parameters.size());
    for (const double parameter : parameters) {
        stations.push_back(arc_prefix_station(motion, parameter));
    }
    const FT radius(tool_radius);
    const FT ratio(cap_chord_ratio);
    const FT gap(gap_close_ratio);
    for (const std::size_t index : decision_probe_order(stations.size())) {
        if (arc_swept_prefix_station_exceeds_exact(
                stock, motion.guide, motion.start, stations[index].centre,
                motion.clockwise, stations[index].extent, radius, ratio, gap)) {
            return true;
        }
    }
    return false;
}

struct LinePrefixMotion2 {
    EPoint start;
    EPoint end;
};

LinePrefixMotion2 line_prefix_motion(
    double x0, double y0, double x1, double y1,
    const std::vector<double>& parameters, double tool_radius,
    double cap_chord_ratio, double gap_close_ratio)
{
    validate_engagement_input_binary64(
        x0, y0, tool_radius, cap_chord_ratio, gap_close_ratio);
    validate_engagement_input_binary64(
        x1, y1, tool_radius, cap_chord_ratio, gap_close_ratio);
    validate_prefix_parameters(parameters);
    LinePrefixMotion2 motion{EPoint{FT(x0), FT(y0)}, EPoint{FT(x1), FT(y1)}};
    if (motion.start == motion.end) {
        throw InvalidSweptPrefixInputError(
            "line swept-prefix motion requires distinct endpoints");
    }
    return motion;
}

ArcPrefixMotion2 arc_prefix_motion(
    double guide_cx, double guide_cy,
    double start_x, double start_y, double end_x, double end_y, bool clockwise,
    const std::vector<double>& parameters, double tool_radius,
    double cap_chord_ratio, double gap_close_ratio)
{
    validate_engagement_input_binary64(
        start_x, start_y, tool_radius, cap_chord_ratio, gap_close_ratio);
    validate_engagement_input_binary64(
        end_x, end_y, tool_radius, cap_chord_ratio, gap_close_ratio);
    if (!std::isfinite(guide_cx) || !std::isfinite(guide_cy)) {
        throw InvalidSweptPrefixInputError(
            "arc swept-prefix guide center must be finite");
    }
    validate_prefix_parameters(parameters);
    const EPoint guide{FT(guide_cx), FT(guide_cy)};
    const EPoint start{FT(start_x), FT(start_y)};
    const EPoint end{FT(end_x), FT(end_y)};
    if (guide == start || guide == end) {
        throw InvalidSweptPrefixInputError(
            "arc swept-prefix guide radius must be positive");
    }
    const bool full_turn = start == end;
    // INTENT, not geometry: the travel sweep in radians that places the interior
    // stations. It crosses into exact-land only through the chart seam; no
    // decision reads it.
    double total_sweep = FULL_TURN;
    if (!full_turn) {
        const double start_angle = std::atan2(start_y - guide_cy, start_x - guide_cx);
        const double end_angle = std::atan2(end_y - guide_cy, end_x - guide_cx);
        total_sweep = clockwise ? start_angle - end_angle : end_angle - start_angle;
        if (total_sweep <= 0.0) total_sweep += FULL_TURN;
    }
    return {guide, start, start, end, start - guide, full_turn, total_sweep, clockwise};
}

ArcPrefixMotion2 circle_prefix_motion(
    double guide_cx, double guide_cy, double guide_radius,
    double entry_x, double entry_y, bool clockwise,
    const std::vector<double>& parameters, double tool_radius,
    double cap_chord_ratio, double gap_close_ratio)
{
    validate_engagement_input_binary64(
        guide_cx, guide_cy, tool_radius, cap_chord_ratio, gap_close_ratio);
    validate_prefix_parameters(parameters);
    if (!std::isfinite(guide_radius) || !std::isfinite(entry_x)
        || !std::isfinite(entry_y)) {
        throw InvalidSweptPrefixInputError(
            "circle swept-prefix guide radius and entry must be finite");
    }
    const FT radius(guide_radius);
    if (CGAL::sign(radius) != CGAL::POSITIVE) {
        throw InvalidSweptPrefixInputError(
            "circle swept-prefix guide radius must be positive");
    }
    const EPoint guide{FT(guide_cx), FT(guide_cy)};
    const EPoint entry{FT(entry_x), FT(entry_y)};
    if (entry == guide) {
        throw InvalidSweptPrefixInputError(
            "circle swept-prefix entry must differ from the guide centre");
    }
    // Station 0 is the authored ENTRY itself: the empty prefix sees the tool
    // where it physically is, which is where the previous motion ended (Ruling
    // 0.2-station0). Stations p > 0 lie on the AUTHORED radius: the first is
    // rho times an exact rational unit vector along the entry direction (that
    // direction is placement intent), the rest rotate exactly from it, all on
    // the circle the replay's annulus is built on.
    const EPoint start =
        guide + radius * rational_unit_along(entry_x - guide_cx, entry_y - guide_cy);
    return {guide, entry, start, start, start - guide, true, FULL_TURN, clockwise};
}

} // namespace

double cap_chord_ratio(double cap_radians)
{
    return audit_cap_chord_ratio(cap_radians);
}

bool cap_chord_ratio_le(double lhs, double rhs)
{
    if (!(lhs > 0.0 && lhs <= 4.0))
        throw InvalidEngagementCapRatioError(
            "lhs cap_chord_ratio must be in (0, 4]");
    if (!(rhs > 0.0 && rhs <= 4.0))
        throw InvalidEngagementCapRatioError(
            "rhs cap_chord_ratio must be in (0, 4]");

    return CGAL::compare(FT(lhs), FT(rhs)) != CGAL::LARGER;
}

EngagementSample engagement_at(const Stock2& stock, double cx, double cy,
                               double tool_radius, double cap_chord_ratio,
                               double gap_close_ratio)
{
    validate_engagement_input_binary64(
        cx, cy, tool_radius, cap_chord_ratio, gap_close_ratio);
    // API-boundary contract: cap_chord_ratio = 4*sin^2(cap/2) with 0 < cap <= pi
    // lies in (0, 4]. Validate the raw double before exact injection (NaN fails).
    // gap_close_ratio = 4*sin^2(gamma/2) with 0 <= gamma <= pi lies in [0, 4];
    // 0 (no gap closed) is the default and the pre-pessimism semantics.

    const AuditExactStationClassification2 exact =
        classify_audit_unguarded_station_exact(
            stock,
            EPoint(FT(cx), FT(cy)),
            FT(tool_radius),
            FT(cap_chord_ratio),
            FT(gap_close_ratio));
    return engagement_sample_from_exact(exact, cx, cy);
}

std::vector<EngagementSample> engagement_line_prefixes(
    const Stock2& stock, double x0, double y0, double x1, double y1,
    const std::vector<double>& parameters, double tool_radius,
    double cap_chord_ratio, double gap_close_ratio)
{
    const LinePrefixMotion2 motion = line_prefix_motion(
        x0, y0, x1, y1, parameters, tool_radius, cap_chord_ratio,
        gap_close_ratio);
    std::vector<EngagementSample> result;
    result.reserve(parameters.size());
    for (const double parameter : parameters) {
        const FT exact_parameter(parameter);
        const AuditExactStationClassification2 exact =
            classify_line_swept_prefix_station_exact(
                stock,
                motion.start,
                motion.end,
                exact_parameter,
                FT(tool_radius),
                FT(cap_chord_ratio),
                FT(gap_close_ratio));
        const EPoint station =
            motion.start + exact_parameter * (motion.end - motion.start);
        result.push_back(engagement_sample_from_exact(
            exact,
            CGAL::to_double(station.x()),
            CGAL::to_double(station.y())));
    }
    return result;
}

bool engagement_line_exceeded(
    const Stock2& stock, double x0, double y0, double x1, double y1,
    const std::vector<double>& parameters, double tool_radius,
    double cap_chord_ratio, double gap_close_ratio)
{
    const LinePrefixMotion2 motion = line_prefix_motion(
        x0, y0, x1, y1, parameters, tool_radius, cap_chord_ratio,
        gap_close_ratio);
    const FT radius(tool_radius);
    const FT ratio(cap_chord_ratio);
    const FT gap(gap_close_ratio);
    for (const std::size_t index : decision_probe_order(parameters.size())) {
        if (line_swept_prefix_station_exceeds_exact(
                stock, motion.start, motion.end, FT(parameters[index]),
                radius, ratio, gap)) {
            return true;
        }
    }
    return false;
}

std::vector<EngagementSample> engagement_arc_prefixes(
    const Stock2& stock, double guide_cx, double guide_cy,
    double start_x, double start_y, double end_x, double end_y, bool clockwise,
    const std::vector<double>& parameters, double tool_radius,
    double cap_chord_ratio, double gap_close_ratio)
{
    return arc_prefix_samples(
        stock,
        arc_prefix_motion(
            guide_cx, guide_cy, start_x, start_y, end_x, end_y, clockwise,
            parameters, tool_radius, cap_chord_ratio, gap_close_ratio),
        parameters, tool_radius, cap_chord_ratio, gap_close_ratio);
}

bool engagement_arc_exceeded(
    const Stock2& stock, double guide_cx, double guide_cy,
    double start_x, double start_y, double end_x, double end_y, bool clockwise,
    const std::vector<double>& parameters, double tool_radius,
    double cap_chord_ratio, double gap_close_ratio)
{
    return arc_prefix_exceeded(
        stock,
        arc_prefix_motion(
            guide_cx, guide_cy, start_x, start_y, end_x, end_y, clockwise,
            parameters, tool_radius, cap_chord_ratio, gap_close_ratio),
        parameters, tool_radius, cap_chord_ratio, gap_close_ratio);
}

std::vector<EngagementSample> engagement_circle_prefixes(
    const Stock2& stock, double guide_cx, double guide_cy, double guide_radius,
    double entry_x, double entry_y, bool clockwise,
    const std::vector<double>& parameters, double tool_radius,
    double cap_chord_ratio, double gap_close_ratio)
{
    return arc_prefix_samples(
        stock,
        circle_prefix_motion(
            guide_cx, guide_cy, guide_radius, entry_x, entry_y, clockwise,
            parameters, tool_radius, cap_chord_ratio, gap_close_ratio),
        parameters, tool_radius, cap_chord_ratio, gap_close_ratio);
}

bool engagement_circle_exceeded(
    const Stock2& stock, double guide_cx, double guide_cy, double guide_radius,
    double entry_x, double entry_y, bool clockwise,
    const std::vector<double>& parameters, double tool_radius,
    double cap_chord_ratio, double gap_close_ratio)
{
    return arc_prefix_exceeded(
        stock,
        circle_prefix_motion(
            guide_cx, guide_cy, guide_radius, entry_x, entry_y, clockwise,
            parameters, tool_radius, cap_chord_ratio, gap_close_ratio),
        parameters, tool_radius, cap_chord_ratio, gap_close_ratio);
}

CertifiedTea certify_segment_tea(const Stock2& stock, double x0, double y0,
                                 double x1, double y1, double tool_radius,
                                 double cap_radians)
{
    validate_segment_certification_input_binary64(
        x0, y0, x1, y1, tool_radius, cap_radians);
    // BOUNDARY (docs/exactness.md, boundary doctrine): validate the ergonomic
    // angular cap and convert it to its exact chord surrogate here, at the one
    // declared seam. Contract 0 < cap <= pi: a single engaged run subtends at
    // most a half turn before the >pi case is an exact orientation verdict, so a
    // cap above pi is meaningless. The full-cap ratio is carried down for max_tea
    // reporting; per-level guarded ratios are derived inside the recursion.
    const double full_cap_chord_ratio = cap_chord_ratio(cap_radians);

    CertifiedTea acc{0.0, true, 0};
    certify_recursive(stock, x0, y0, x1, y1, tool_radius, cap_radians,
                      full_cap_chord_ratio, acc, 0);
    return acc;
}

int sign_mixed_radical_for_binding(
    double a,
    double b,
    double c,
    double d,
    double alpha,
    double beta)
{
    validate_mixed_radical_input_binary64(
        a, b, c, d, alpha, beta);
    return static_cast<int>(audit_sign_mixed_radical_exact(
        FT(a), FT(b), FT(c), FT(d), FT(alpha), FT(beta)));
}
