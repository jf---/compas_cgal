#pragma once

// Owner-by-construction TOUCH: the structural decision dist(piece, owner) == r
// for a boundary piece of Omega* and the declared primitive it offsets.
//
// A piece of dOmega* on an offset support carries its owner (benchmarks/
// held_centre_domain_supports.py, BoundaryFeature LINE i / ARC i). The owner
// relation is accepted ONLY where one of the two lemmas below proves the
// distance; their hypotheses are checked exactly, every one of them, and a
// pair that fails any hypothesis is decided like every other pair (filter,
// crossing, closed-form candidates). Every non-owner primitive is always
// decided exactly.
//
// Lemma L (line owner). Let L = {S + lambda u : lambda in [0, Lambda]} with
// |u| = 1 and left normal n = u^perp. Let a segment have endpoints p0, p1
// with n.(p_i - S) == r and u.(p_i - S) in [0, Lambda] exactly (i = 0, 1).
// Then every point x of the segment is a convex combination of p0 and p1, so
// n.(x - S) == r and u.(x - S) in [0, Lambda]: the foot S + (u.(x - S)) u lies
// on L and |x - foot| == r, while every other point of L is farther (Pythagoras
// on the line through the foot). Hence dist(x, L) == r for every x, and
// dist(segment, L) == r: TOUCH.
//
// Lemma A (arc owner). Let A be the arc of the circle (C, rho) from S to E in
// sense sigma (not a full circle). Let an arc element lie on the circle
// (C, R') with the same centre, R' > 0 and |R' - rho| == r, and let its sweep
// from s to t (in its own sense) lie in A's angular span: in A's sense sigma,
// the directions of both ends lie in [S, E], and the element's sweep runs
// forwards from its first to its last direction within that span (the
// element is not a full circle). Then for every point x of the element the
// ray from C through x meets A; the nearest point of the full circle (C, rho)
// to x is on that ray at distance | |x - C| - rho | = |R' - rho| = r, and it
// lies on A. Hence dist(x, A) == r for every x, and dist(element, A) == r:
// TOUCH. Directions are compared exactly (compare_on_arc; only directions
// from C count, so the declared ends S, E and the element ends s, t may lie
// on different circles and in different extensions).
//
// Junction owners (r-circles about a declared vertex) are not covered by a
// lemma here: those pieces are decided exactly against every primitive.

#include "held_declared_clearance_2/elements_2.h"
#include "held_declared_clearance_2/primitives_2.h"

namespace compas_cgal::held_declared_clearance {

/// Lemma L's hypotheses hold exactly for this segment and line owner.
[[nodiscard]] bool line_owner_touches(const SegmentElement& element, const Primitive& owner, const FT& tool_radius);
/// Lemma A's hypotheses hold exactly for this arc element and arc owner.
[[nodiscard]] bool arc_owner_touches(const ArcElement& element, const Primitive& owner, const FT& tool_radius);

}  // namespace compas_cgal::held_declared_clearance
