#pragma once

#include "compas.h"

class Stock2;

// Tool-engagement-angle sample at one cutter station. Angles in radians and are
// REPORTING quantities (doubles produced for humans/statistics): total_tea sums
// the angular extent of every rim arc in contact with material, max_run_tea is
// the largest maximal contiguous engaged run. cap_exceeded is the only DECISION
// carried out here and it is decided EXACTLY on the exact arrangement, never
// from these reported doubles (see engagement_2.cpp).
//
// DECIDING / REPORTING SPLIT under gap-closure pessimism (see engagement_at):
// total_tea and max_run_tea always describe the TRUE runs (no gap closed);
// cap_exceeded is decided on the PESSIMISTIC runs (void gaps <= gap_close_ratio
// pre-absorbed). With the default gap_close_ratio == 0 the two coincide and
// every field is bit-for-bit the pre-pessimism result.
struct EngagementSample {
    double total_tea;
    double max_run_tea;
    bool cap_exceeded;
};

// Exact TEA query, read LOCALLY. Zone the cutter circle of radius tool_radius at
// (cx, cy) in the stock's OWN Gps arrangement (faces carry contained() =
// material) and harvest the engaged rim sub-arcs (the rim-in-material pieces),
// assemble maximal engaged runs by EXACT endpoint equality, and certify each run
// against the engagement cap. The query is O(cutter-crossings), not O(stock) --
// no whole-stock overlay (see engagement_2.cpp; the earlier overlay was validated
// to reproduce this exactly before it was removed).
//
// (cx, cy) must be finite and tool_radius finite and strictly positive; all are
// checked BEFORE any exact injection and raise std::invalid_argument otherwise.
// A non-finite double is not a rational, so it has no exact image at this seam,
// and a non-positive radius is not a cutter.
//
// cap_chord_ratio is the dimensionless squared-chord surrogate for the angular
// cap: cap_chord_ratio = 4*sin^2(cap/2), which the CALLER computes as a double
// and which is injected exactly (Epeck::FT(double)). Contractually the cap is
// 0 < cap <= pi, so cap_chord_ratio lies in (0, 4]; values outside that range
// raise std::invalid_argument. The certificate is an exact statement about the
// exact rational threshold T = FT(cap_chord_ratio) * FT(tool_radius)^2; the
// sub-ulp gap between this rational surrogate and the transcendental angle cap
// is API semantics documented here, NOT an in-core correction constant.
//
// gap_close_ratio is the analogous squared-chord surrogate 4*sin^2(gamma/2) for
// a gap-closure angle gamma in [0, pi] (so gap_close_ratio lies in [0, 4]; out
// of range raises std::invalid_argument). GAP-CLOSURE PESSIMISM: before deciding
// the cap, every VOID gap between consecutive engaged runs whose angular span is
// <= gamma is absorbed, merging its two bounding runs (chained to a fixpoint) so
// the cap DECISION sees the PESSIMISTIC runs. This closes the certificate's merge
// hole: two runs a hair of void apart at a station are treated as already merged,
// so a merge that completes within a certifier step cannot slip past unseen. The
// closure test is the SAME exact orientation+chord predicate used for runs (a
// gap absorbed iff its span does not exceed gamma), applied to the gap's exact
// one-root endpoints -- no angle sums. The default 0.0 (gamma = 0) closes no gap:
// pessimistic runs == true runs, and every result is the pre-pessimism value
// bit-for-bit. Reported total_tea/max_run_tea always describe the TRUE runs.
EngagementSample engagement_at(const Stock2& stock, double cx, double cy,
                               double tool_radius, double cap_chord_ratio,
                               double gap_close_ratio = 0.0);

// Result of certifying the engagement cap along one linear cutter motion.
// max_tea is the largest run TEA seen at any station visited -- a REPORTING
// double, best-effort over the sampled stations, never a decision input.
// cap_certified is the DECISION: true iff no cutter center on the motion can
// exceed the cap, established purely from EXACT station verdicts plus the
// analytic guard (see certify_segment_tea); false means a station violated the
// guarded cap and the margin could not be closed by refinement. stations counts
// the sub-intervals examined.
struct CertifiedTea {
    double max_tea;
    bool cap_certified;
    int stations;
};

// Analytic TEA-growth bound and the guard derived from it, exposed so the Python
// audit layer CALLS the certifier's guard instead of mirroring it. A mirrored
// safety constant that drifts turns a conservative certifier unsound in silence.
// REFINEMENT bound only, never a geometric decision (docs/exactness.md).
//
// NOT LOAD-BEARING FOR THE CERTIFICATE ANY MORE. tea_growth_bound is measurably
// NOT an upper bound on TEA growth (tests/test_growth_bound.py falsifies it by up
// to 24x on concave features), and the false certificates that follow are
// committed as tests/test_false_certificate.py. The certificate's interior
// guarantee is now swept_run_bound below; these two survive as the retired
// lemma's record and as the certifier's retained (strictly conservative) station
// filter -- see certify_segment_tea and engagement_2.cpp.
double tea_growth_bound(double d, double r);
double tea_guard(double d, double r);

// EXACT SWEPT-ANNULUS BOUND -- the certificate's interior guarantee.
//
// Upper bound, in radians, on the largest contiguous engaged run ANY cutter
// centre P with |P - (cx, cy)| <= travel can see against the frozen stock. Unlike
// tea_growth_bound this bounds the run ITSELF, not its growth, so it needs no
// premise about what an existing run does between two stations -- the premise the
// rib witnesses break.
//
// The bound is `max_K ang(K) + 2*asin(travel / (tool_radius - travel))`, where K
// ranges over the CONNECTED COMPONENTS of the material inside the annulus
// A((cx,cy), r - travel, r + travel) and ang(K) is K's angular extent seen from
// (cx, cy). Derivation, the exact component/wrap decision, and the safe failure
// direction are in engagement_2.cpp.
//
// Returns 2*pi -- the saturated value, which forces refinement or a refusal
// against any cap in (0, pi] -- whenever the construction cannot see the swept
// rim (travel == 0, a degenerate annulus) or cannot bound a component below a
// full turn (travel >= tool_radius/2, or a component wide enough to wrap the
// centre). Returns 0.0 when no material lies in the annulus at all, which is
// exact: no reachable centre touches anything.
//
// (cx, cy) must be finite, tool_radius finite and strictly positive, and travel
// finite and non-negative; all are checked BEFORE any exact injection and raise
// std::invalid_argument otherwise (same seam contract as engagement_at).
double swept_run_bound(const Stock2& stock, double cx, double cy,
                       double tool_radius, double travel);

// Certify TEA(P) <= cap_radians for EVERY cutter center P on the segment
// (x0,y0)->(x1,y1) with tool radius tool_radius, against the frozen stock.
//
// Method -- adaptive station sampling with an exact station test and an exact
// SWEPT-ANNULUS interior bound. A span is certified iff BOTH hold at
// half-station-spacing d:
//
//   (i)  INTERIOR (the soundness argument). swept_run_bound(station, r, d) <= cap
//        at both stations. Every center on the span lies within d of the nearer
//        station, so this bounds the largest engaged run at EVERY center, not
//        merely at the two measured ones. This alone proves the claim.
//   (ii) STATIONS (retained conservative filter). Each station is measured by the
//        EXACT engagement_at cap predicate against the guarded cap cap - guard(d)
//        with gap-closure pessimism, exactly as before.
//
// (ii) is NOT part of the proof -- guard() rests on tea_growth_bound, which is
// measurably not an upper bound -- and is retained only because it can solely
// REFUSE, never certify: it is the pre-existing path, kept alongside the new one
// until its removal is decided separately. If either test fails the span is
// bisected; on reaching the spacing floor with the margin still open the motion
// is reported uncertified. max_tea stays reporting.
//
// BOUNDARY (docs/exactness.md, boundary doctrine): cap_radians is validated to
// (0, pi] and converted to exact chord surrogates here, at the one declared
// seam; every value crossing into the exact station test is injected exactly.
// Raises std::invalid_argument if any endpoint coordinate is non-finite, if
// tool_radius is not finite and strictly positive, or if cap_radians lies
// outside (0, pi]. tool_radius additionally sets the refinement's spacing floor
// (STATION_FLOOR_FRACTION * tool_radius), so a non-physical radius would leave
// the bisection with no reachable floor.
CertifiedTea certify_segment_tea(const Stock2& stock, double x0, double y0,
                                 double x1, double y1, double tool_radius,
                                 double cap_radians);

// Engagement queries share the _stock_2 nanobind module (NB_STATIC forbids
// cross-module type sharing), so registration is a hook the module macro in
// stock_2.cpp calls. It also exposes a test-only `_sign_mixed_radical` binding
// for direct unit tests of the exact cap predicate's core primitive.
void register_engagement(nanobind::module_& m);
