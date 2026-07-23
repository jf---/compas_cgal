#pragma once

#include "compas.h"

class Stock2;

// Tool-engagement-angle sample at one cutter station. Angles in radians and are
// REPORTING quantities (doubles produced for humans/statistics): total_tea sums
// the angular extent of every rim arc in contact with material, max_run_tea is
// the largest maximal contiguous engaged run. cap_exceeded is the only DECISION
// carried out here and it is decided EXACTLY on the exact arrangement, never
// from these reported doubles (see engagement_2.cpp).
struct EngagementSample {
    double total_tea;
    double max_run_tea;
    bool cap_exceeded;
};

// Exact TEA query. Intersect the stock with the cutter disk of radius
// tool_radius at (cx, cy), harvest the boundary arcs supported by the cutter
// circle (the rim-in-material pieces), assemble maximal engaged runs by EXACT
// endpoint equality, and certify each run against the engagement cap.
//
// cap_chord_ratio is the dimensionless squared-chord surrogate for the angular
// cap: cap_chord_ratio = 4*sin^2(cap/2), which the CALLER computes as a double
// and which is injected exactly (Epeck::FT(double)). Contractually the cap is
// 0 < cap <= pi, so cap_chord_ratio lies in (0, 4]; values outside that range
// raise std::invalid_argument. The certificate is an exact statement about the
// exact rational threshold T = FT(cap_chord_ratio) * FT(tool_radius)^2; the
// sub-ulp gap between this rational surrogate and the transcendental angle cap
// is API semantics documented here, NOT an in-core correction constant.
EngagementSample engagement_at(const Stock2& stock, double cx, double cy,
                               double tool_radius, double cap_chord_ratio);

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

// Certify TEA(P) <= cap_radians for EVERY cutter center P on the segment
// (x0,y0)->(x1,y1) with tool radius tool_radius, against the frozen stock.
//
// Method -- adaptive station sampling with a guarded exact test. Each station is
// measured by the EXACT engagement_at cap predicate, but against a GUARDED cap:
// at half-station-spacing d the station threshold is cap - guard(d), where guard
// is a conservative analytic TEA-growth bound carrying an explicit integer
// safety factor (derivation in engagement_2.cpp). Every center lies within its
// half-spacing of the nearer measured station, so two stations both under the
// guarded cap certify the whole span at TEA <= cap. If the guarded cap is
// exceeded -- or is non-positive, the spacing being too coarse to admit any
// guard -- the span is bisected; on reaching the spacing floor with the margin
// still open the motion is reported uncertified. The verdict is thus EXACT
// station predicates + the analytic guard ONLY; max_tea is reporting.
//
// BOUNDARY (docs/exactness.md, boundary doctrine): cap_radians is validated to
// (0, pi] and converted to exact chord surrogates here, at the one declared
// seam; every value crossing into the exact station test is injected exactly.
// Raises std::invalid_argument if cap_radians lies outside (0, pi].
CertifiedTea certify_segment_tea(const Stock2& stock, double x0, double y0,
                                 double x1, double y1, double tool_radius,
                                 double cap_radians);

// Engagement queries share the _stock_2 nanobind module (NB_STATIC forbids
// cross-module type sharing), so registration is a hook the module macro in
// stock_2.cpp calls. It also exposes a test-only `_sign_mixed_radical` binding
// for direct unit tests of the exact cap predicate's core primitive.
void register_engagement(nanobind::module_& m);
