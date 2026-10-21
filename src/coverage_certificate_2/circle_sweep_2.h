#pragma once

// The cover-counter sweep around one boundary circle (design steps 2-3).
//
// The crossing events on circle C are put in CCW order from the leftmost point
// (lower half by ascending x, then upper half by descending x; exact
// comparisons behind the interval front, circle_events_2.h). Events at exactly
// equal points form one arrangement vertex. The open sub-arcs between
// consecutive vertices are the arrangement edges on C.
//
// Events come from three kinds of support:
//   * boundary circles of candidate circle regions (cuts, capsule caps): a
//     transversal crossing flips one bit of its region (inside the outer
//     circle / inside the inner circle);
//   * the four lines of every candidate rectangle (lines_2.h): a transversal
//     crossing flips that line's side; a rectangle covers an edge iff all four
//     sides are its inward sides;
//   * target cut supports (target vertex circles and target flank lines, for a
//     swept circle only): they cut C so that T-membership is constant along
//     each edge, and flip nothing.
// Tangencies flip nothing. The start state is decided at a rational witness on
// C that lies on no support (exact P4). The counters are bookkeeping: the
// decision that an edge is exposed is re-confirmed at an exact witness
// (witness_2.h).
//
// Coincident supports (P5/P6): a candidate region whose boundary circle IS C
// never has C in its interior. If it lies on the side of C the certificate
// checks (BoundaryCircle::check_inside) it covers every edge of C, and C is
// covered whole. Otherwise it contributes nothing and is excluded from the
// counters. Target circles are never counted.

#include "coverage_certificate_2/circle_events_2.h"
#include "coverage_certificate_2/counters_2.h"
#include "coverage_certificate_2/lines_2.h"
#include "coverage_certificate_2/regions_2.h"

#include <cstdint>
#include <functional>
#include <vector>

namespace compas_cgal::coverage_certificate {

/// Candidates of one circle beyond its circle regions (CrossingTable::candidates).
struct CircleCandidates {
    std::vector<std::uint32_t> rects;  ///< swept rectangle indices, sorted
    std::vector<std::uint32_t> lines;  ///< every line with events on the circle (rect lines + target cuts), sorted
    std::vector<std::uint32_t> vertices;  ///< one-root vertex circles (target cuts), sorted
};

/// The arrangement of one boundary circle after the sweep.
struct CircleSweep {
    std::vector<CrossingEvent> events;  ///< CCW from the leftmost point
    std::vector<std::uint8_t> halves;   ///< 0 lower, 1 upper, per event
    std::vector<std::uint32_t> group_starts;  ///< vertex v spans events [group_starts[v], group_starts[v + 1])
    /// Edge v runs from vertex v to vertex v + 1; the last edge wraps through
    /// the leftmost point back to vertex 0. With no vertex, one edge is the
    /// whole circle.
    std::vector<std::uint8_t> edge_exposed;
    std::vector<std::uint32_t> counted;   ///< candidate circle regions that take part in the counters
    std::vector<std::uint32_t> excluded;  ///< same-side coincident candidates (C lies on their boundary)
    std::vector<std::uint32_t> rects;     ///< candidate rectangles (all counted: a circle never shares a line)
    std::vector<std::uint32_t> lines;     ///< line supports with events on C
    std::vector<std::uint32_t> vertices;  ///< one-root vertex circles cutting C
    bool covered_by_coincident_support = false;

    [[nodiscard]] std::size_t vertex_count() const { return group_starts.empty() ? 0 : group_starts.size() - 1; }
};

/// Theorem bound on the start-witness attempts of circle C: the witnesses are
/// distinct points (the leftmost point, then s = 1/2, 3/2, ...), and every support
/// not coincident with C meets it in at most two points, so at most
/// 2 (2 counted + target_circles + lines) attempts can fail before one lies on no
/// support (a counted circle region has up to two boundary circles).
[[nodiscard]] std::uint32_t start_witness_attempt_bound(std::size_t counted, std::size_t target_circles, std::size_t lines);

/// A rational start witness on C that lies on no support.
struct StartWitness {
    Point_2 point;
    std::uint8_t half;      ///< 0 lower, 1 upper (stereographic half)
    std::uint32_t attempt;  ///< 0 for the leftmost point
};

/// The start-witness search: the leftmost point, then s = 1/2, 3/2, 5/2, ...,
/// until `on_support` reports a witness on no support.
///
/// Raises:
///     StartWitnessSearchError: more than `bound` attempts lay on a support
///         (an uncut coincident support; see start_witness_attempt_bound).
[[nodiscard]] StartWitness find_start_witness(const BoundaryCircle& C, std::uint32_t bound,
                                              const std::function<bool(const Point_2&)>& on_support,
                                              DecisionCounters& counters);

/// Run the sweep on circle k.
///
/// Raises:
///     CoverageCertificateNativeError: the counters fail parity around C.
///     StartWitnessSearchError: see find_start_witness.
[[nodiscard]] CircleSweep sweep_circle(std::uint32_t k, const CircleRegionSet& set, const LineTable& lines,
                                       const std::vector<RectRegion>& rects, const CircleCandidates& candidates,
                                       CircleEvents& events, DecisionCounters& counters);

/// Position of the rational point q on circle C relative to event e in the CCW
/// order from the leftmost point: SMALLER when e comes first, EQUAL when they
/// are the same point, LARGER when q comes first. `q_half` is q's half.
[[nodiscard]] CGAL::Comparison_result order_event_point(CircleEvents& events, const CrossingEvent& e,
                                                        std::uint8_t e_half, const Point_2& q, const Interval& qx,
                                                        std::uint8_t q_half, DecisionCounters& counters);

}  // namespace compas_cgal::coverage_certificate
