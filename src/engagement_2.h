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

// Engagement queries share the _stock_2 nanobind module (NB_STATIC forbids
// cross-module type sharing), so registration is a hook the module macro in
// stock_2.cpp calls. It also exposes a test-only `_sign_mixed_radical` binding
// for direct unit tests of the exact cap predicate's core primitive.
void register_engagement(nanobind::module_& m);
