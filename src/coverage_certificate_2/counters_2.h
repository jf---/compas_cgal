#pragma once

#include <cstddef>

namespace compas_cgal::coverage_certificate {

/// Reporting-only tallies of how each decision was taken. They never feed a
/// decision; the contract tests read them to prove that a fixture really
/// exercised the exact path it claims to.
struct DecisionCounters {
    std::size_t candidate_pairs = 0;          ///< bbox-overlapping boundary-circle pairs, distinct owners
    std::size_t pairs_filtered_disjoint = 0;  ///< interval discriminant < 0: no crossing, no exact work
    std::size_t pairs_filtered_transversal = 0;  ///< interval discriminant > 0: two crossings, kept lazy
    std::size_t pairs_exact = 0;              ///< discriminant enclosure straddles 0 (or EXACT_ONLY)
    std::size_t tangencies = 0;               ///< exact double roots
    std::size_t coincident_pairs = 0;         ///< identical support circles (P5)
    std::size_t materialised_pairs = 0;       ///< lazy pairs whose exact points were built on demand
    std::size_t pair_towers = 0;              ///< circle-circle pairs whose exact tower coordinates were built (P2)
    std::size_t exact_branch_matches = 0;     ///< materialised lazy pairs whose two branch boxes overlap:
                                              ///< the exact rational branch matcher decided the labels
    std::size_t exact_half = 0;               ///< exact lower/upper-half decisions
    std::size_t exact_compare_x = 0;          ///< exact x comparisons of crossing points on one circle
    std::size_t exact_membership = 0;         ///< exact squared-distance memberships
    std::size_t exact_target_membership = 0;  ///< exact target-member memberships
    std::size_t start_retries = 0;            ///< start witnesses that fell on a candidate boundary
    std::size_t witness_bisections = 0;       ///< P3 exact bisection steps beyond the double seed
    std::size_t uncovered_halvings = 0;       ///< P7 halvings of the outward offset
    std::size_t boundary_circles = 0;
    std::size_t sub_arcs = 0;
    std::size_t covered_by_coincident_support = 0;  ///< circles covered whole by an opposite-side support
    // Capsule and target pieces (phase 2).
    std::size_t circle_line_pairs = 0;          ///< candidate (circle, line) pairs, on either chart
    std::size_t circle_line_exact = 0;          ///< discriminant enclosure straddles 0 (or EXACT_ONLY)
    std::size_t line_line_exact = 0;            ///< parallelism / crossing decided exactly
    std::size_t line_tangencies = 0;            ///< exact circle-line double roots
    std::size_t materialised_line_events = 0;   ///< line-event towers built on demand
    std::size_t exact_tower_compares = 0;       ///< tower comparisons the enclosures could not decide
    std::size_t cross_base_separations = 0;     ///< comparisons over two irrational bases (provenance + bisection)
    std::size_t exact_line_side = 0;            ///< exact sides of rational points relative to lines
    std::size_t line_pieces = 0;
    std::size_t line_sub_edges = 0;
    std::size_t coincident_line_supports = 0;   ///< (line piece, rect) pairs sharing a support
    std::size_t target_pieces = 0;
    std::size_t sweep_exposed_outside_target = 0;  ///< sweep-exposed sub-edges whose outer side is not in T
    std::size_t target_decisions = 0;              ///< target-side decisions at sweep-exposed witnesses
    // One-root vertex circles of the exact target (vertex_events_2.h).
    std::size_t vertex_circle_pairs = 0;         ///< candidate (rational circle, vertex circle) pairs
    std::size_t vertex_circle_exact = 0;         ///< their discriminant enclosure straddles 0 (or EXACT_ONLY)
    std::size_t vertex_line_pairs = 0;           ///< candidate (line, vertex circle) pairs
    std::size_t vertex_line_exact = 0;           ///< their discriminant enclosure straddles 0 (or EXACT_ONLY)
    std::size_t vertex_circle_tangencies = 0;    ///< exact double roots with a vertex circle
    std::size_t materialised_vertex_events = 0;  ///< vertex-circle events built exactly on demand
    std::size_t exact_target_events = 0;         ///< event comparisons involving a vertex circle, decided exactly

    [[nodiscard]] std::size_t exact_fallbacks() const
    {
        return pairs_exact + materialised_pairs + pair_towers + exact_half + exact_compare_x + exact_membership +
               exact_target_membership + circle_line_exact + line_line_exact + materialised_line_events +
               exact_tower_compares + exact_line_side + vertex_circle_exact + vertex_line_exact + materialised_vertex_events +
               exact_target_events;
    }
};

}  // namespace compas_cgal::coverage_certificate
