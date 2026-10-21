#pragma once

// Exact witnesses on an exposed circle edge (design steps P3, P4, P7).
//
// P3: a rational stereographic parameter s strictly inside the edge. A double
//     seed (midpoint angle of the endpoints, from their EXACT values) is only a
//     proposal; it is accepted by exact CCW comparisons against the two exact
//     endpoint vertices. When it is rejected, an exact bisection in s runs: the
//     edge has positive length, so it terminates, and every step is exact.
// P4: the witness w is re-classified against every counted candidate (circle
//     regions and rectangles); it must be in the exterior of all of them. The
//     counters are bookkeeping, this is the decision. Disagreement is a
//     contract breach.
// P7: an uncovered point u = c + (1 +/- h)(w - c), h = 2^-k, on the edge's
//     checked side (outward for an outer/disk/cap boundary, toward the centre
//     for an inner one or a target circle). It is accepted when it lies exactly
//     outside the owner and every candidate region and rectangle (coincident
//     ones included), satisfies the caller's acceptance predicate (u in T for
//     the certificate) and, for an outward offset, lies inside the circle's
//     exact bounding box, which makes the candidate set complete for u: any
//     region containing u has a box overlapping the circle's. The parameters
//     -1, 0, 1 (the extreme points, where an outward offset leaves the box at
//     once) are never used as witnesses.

#include "coverage_certificate_2/circle_events_2.h"
#include "coverage_certificate_2/circle_sweep_2.h"
#include "coverage_certificate_2/counters_2.h"
#include "coverage_certificate_2/lines_2.h"
#include "coverage_certificate_2/regions_2.h"

#include <cstdint>
#include <functional>
#include <vector>

namespace compas_cgal::coverage_certificate {

/// What P4 and P7 test a circle witness against.
struct CircleWitnessContext {
    const CircleRegionSet& set;
    const LineTable& lines;
    const std::vector<RectRegion>& rects;  ///< swept rectangles
    CircleEvents& events;
};

/// P3 alone: a rational stereographic parameter strictly inside edge `edge`.
[[nodiscard]] FT edge_parameter(std::uint32_t k, std::size_t edge, const CircleSweep& sweep, CircleWitnessContext& context,
                                DecisionCounters& counters);

/// Audit oracle for an edge the counters call covered: its P3 witness must lie,
/// by pure exact P4 (no filter), in the interior of some counted candidate.
/// Independent of the counters; the forced check behind `audit=True`.
///
/// Raises:
///     CoverageCertificateNativeError: the witness is in no candidate's
///         interior, or lies on a candidate boundary.
void audit_covered_edge(std::uint32_t k, std::size_t edge, const CircleSweep& sweep, CircleWitnessContext& context,
                        DecisionCounters& counters);

/// P4 alone: the rational witness w on circle k lies in the exterior of every
/// counted candidate.
///
/// Raises:
///     CoverageCertificateNativeError: the counters disagreed with P4.
void confirm_circle_witness(std::uint32_t k, const Point_2& w, const CircleSweep& sweep, CircleWitnessContext& context,
                            DecisionCounters& counters);

/// P3 + P4: the exact on-edge witness of an exposed edge.
///
/// Raises:
///     CoverageCertificateNativeError: the witness lies in, or on, a counted
///         candidate (the counters disagreed with P4).
[[nodiscard]] std::pair<FT, Point_2> confirm_exposed_edge(std::uint32_t k, std::size_t edge, const CircleSweep& sweep,
                                                          CircleWitnessContext& context, DecisionCounters& counters);

/// P7: the uncovered point beside the witness w of an exposed edge.
///
/// Raises:
///     CoverageCertificateNativeError: P7 exceeded its contract bound.
[[nodiscard]] Point_2 uncovered_beside(std::uint32_t k, const Point_2& w, const CircleSweep& sweep,
                                       CircleWitnessContext& context, const std::function<bool(const Point_2&)>& accept,
                                       DecisionCounters& counters);

}  // namespace compas_cgal::coverage_certificate
