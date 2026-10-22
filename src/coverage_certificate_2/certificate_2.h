#pragma once

// Driver of the coverage certificate (docs/held_coverage_certificate.md).
//
// run_certificate (phase 2): the full local theorem on a polygonal centre
// domain. Regions: full-circle cuts, and connector capsules decomposed into
// end disks (circle regions) and rectangles (lines_2.h), after the exact union
// of collinear connector pieces. Target T = Omega u edge rectangles u vertex
// disks. Pieces checked:
//   * every swept boundary circle (cuts and capsule caps) and every connector
//     flank: exposed iff no other sweep covers its outer side AND that side is
//     in T (target_decision_2.h; the piece is cut at the target supports);
//   * every target flank and vertex circle: exposed iff no sweep covers its
//     member side (the target-edge rule; an exposed target piece is a genuine
//     uncovered point whether or not it lies on dT, so no target-target cut is
//     needed, and every dT piece is among them).
// Rectangle end segments are never exposed: away from their corners they lie in
// the interior of the end disk of the same capsule (or of the vertex disk), so
// they are not pieces. Every exposed edge yields an exact uncovered point
// u in T, outside every sweep (P7); none means T is covered.
//
// run_circle_core (phase 1 semantics, kept for its contract tests): circle
// cuts only, nothing of the target cut into the arrangement; an exposed sweep
// edge carries its uncovered point and whether that point is in T.

#include "coverage_certificate_2/counters_2.h"
#include "coverage_certificate_2/kernel_2.h"
#include "coverage_certificate_2/lines_2.h"
#include "coverage_certificate_2/regions_2.h"
#include "coverage_certificate_2/target_2.h"

#include <cstdint>
#include <vector>

namespace compas_cgal::coverage_certificate {

/// One exposed edge on a full-circle boundary, with exact witnesses (phase 1).
struct ExposedEdge {
    std::uint32_t circle;  ///< boundary-circle index
    std::uint32_t region;  ///< full-circle cut index
    PieceKind kind;
    std::uint32_t edge;  ///< edge index on the circle (CCW from the leftmost point)
    bool whole_circle;   ///< the circle has no crossing vertex
    std::vector<std::uint32_t> start_vertex;  ///< boundary circles through the edge's start vertex, sorted
    std::vector<std::uint32_t> end_vertex;    ///< boundary circles through its end vertex, sorted
    FT parameter;       ///< stereographic s of the on-edge witness
    Point_2 on_edge;    ///< exact point of the open edge
    Point_2 uncovered;  ///< exact point beside it, outside every swept region
    bool uncovered_in_target;  ///< exact closed membership of `uncovered` in T
};

struct CircleCoreResult {
    std::vector<ExposedEdge> exposed;
    DecisionCounters counters;
    std::size_t region_count = 0;
    std::size_t piece_count = 0;
    double elapsed_s = 0.0;  ///< reporting only
    // Phase wall times, reporting only.
    double crossings_s = 0.0;  ///< box candidates + filtered P1
    double sweep_s = 0.0;      ///< CCW order + cover counters
    double witness_s = 0.0;    ///< P3 + P4 + P7 on exposed edges
    double target_s = 0.0;     ///< target membership of uncovered points
};

/// Args:
///     audit: also witness every edge the counters call covered and confirm it
///         by pure exact membership (an oracle independent of the counters;
///         slow, for contract tests). A disagreement raises.
[[nodiscard]] CircleCoreResult run_circle_core(const CircleRegionSet& set, const MachinableTarget& target,
                                               FilterMode mode, bool audit);

/// The kind of an exposed piece of the full certificate.
enum class CertificatePieceKind : std::uint8_t {
    AnnulusOuter = 0,
    AnnulusInner = 1,
    Disk = 2,
    CapsuleCap = 3,
    CapsuleFlank = 4,
    TargetBoundary = 5,
    TargetInterior = 6,  ///< the exact target's interior anchor point is uncovered (exact_certificate_2.h)
    ConnectorArc = 7,    ///< a support-circle arc of an exact dOmega* arc connector (exact_connectors_2.h)
};

/// One exposed edge of the full certificate: an exact proof that T is not covered.
struct CertificateWitness {
    CertificatePieceKind kind;
    std::uint32_t region_index;  ///< full circles, then connectors, then target members (1 + edge)
    std::uint32_t piece;         ///< boundary-circle id, or n_circles + flank-piece index
    std::uint32_t edge;          ///< edge index on the piece
    FT parameter;                ///< chart parameter of the on-edge witness (stereographic s, or flank t)
    Point_2 uncovered;           ///< exact point of T outside every swept region
};

struct CertificateResult {
    std::vector<CertificateWitness> exposed;
    std::size_t undecided = 0;  ///< sweep edges whose side could not be decided (VERTEX_DISK members only)
    DecisionCounters counters;
    std::size_t region_count = 0;  ///< swept regions: cuts + capsule disks + capsule rectangles
    std::size_t piece_count = 0;   ///< checked pieces: circles + flanks
    std::size_t edge_count = 0;    ///< checked edges (sub-arcs + flank sub-edges)
    double elapsed_s = 0.0;        ///< reporting only
    double crossings_s = 0.0, circle_sweep_s = 0.0, line_sweep_s = 0.0, witness_s = 0.0;
};

/// Run the full certificate.
///
/// Args:
///     cuts: full-circle cuts (built by CircleRegionSet::build; capsule caps and
///         target circles are appended here).
///     audit: also confirm every edge the counters call covered by pure exact
///         membership (slow; for contract tests).
[[nodiscard]] CertificateResult run_certificate(CircleRegionSet cuts, const MachinableTarget& target,
                                                const std::vector<ConnectorSpec>& connectors, FilterMode mode, bool audit);

}  // namespace compas_cgal::coverage_certificate
