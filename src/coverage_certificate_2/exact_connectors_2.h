#pragma once

// Exact dOmega* connectors as swept regions (docs/held_coverage_certificate.md,
// "Exact connectors along dOmega*").
//
// An exact connector is a CCW sub-chain of dOmega*: sub-pieces of the exact
// target's own pieces (exact_target_2.h), each on an offset line or an offset
// circle, with one-root end points w0, w1. The tool centre runs along the
// sub-piece, so it sweeps exactly
//
//   line sub-piece:  {s0 <= u.x <= s1, |n.x - level| <= r}  u  D_r(w0)  u  D_r(w1)
//   arc sub-piece:   {(R - r)+ <= |x - c| <= R + r, x in the wedge of the arc}
//                    u  D_r(w0)  u  D_r(w1)
//
// (the sub-piece (+) D_r: a point whose nearest sub-piece point is interior lies
// on its normal, a point whose nearest point is an end lies in that end disk).
// The sub-pieces of all connectors must cover dOmega* (verified, below), so their
// union is exactly dOmega* (+) D_r, the union of piece (+) D_r over the WHOLE
// pieces: the certificate sweeps one composite region per dOmega* PIECE, with end
// disks at the dOmega* vertices only (the sites inside a piece add no support).
// The sub-pieces are kept for validation, the cover check and reporting. A
// piece's composite region has its boundary on its SUPPORTS -- the two flank
// lines n.x = level +/- r (rational: n is a rational unit normal), or the
// circles |x - c| = R +/- r (rational) -- and on its two cap circles |x - w| = r
// about the piece's end vertices w: a vertex circle when w is one-root, a
// rational support region when w is rational. The end lines and radial
// lines are never on the boundary: inside the strip or annulus they are
// diameters of a cap. So the region's membership along any piece is constant
// between consecutive crossings with those supports, and every support is cut
// into every piece it meets; the membership itself is decided exactly at the
// edge witness (the connector supports flip nothing in the cover counters).
//
// Checked pieces: the flank segments and the support-circle arcs of each
// piece's composite region, inside the PIECE's range (an out-of-range edge has
// both sides outside the region). The cap circles are NOT checked. That is a
// theorem whose precondition is K = dOmega*, K the union of all sub-pieces. The
// build verifies it relative to the declared pieces, exactly: the pieces form
// closed chains (every end is exactly one start) with no reflex junction
// (verify_boundary_chains), and the sub-pieces cover every piece. That the
// declared pieces ARE dOmega* is the exact target's contract (they come from the
// construction of Omega*, exact_target_2.h) and is not re-derived here.
//
// Proof. Every uncovered point y' of T has dist(y', K) > r (U contains K (+) D_r),
// so y' lies in int Omega* (a point of T outside Omega* is within r of dOmega*).
// Let e be an edge of a cap circle |x - w| = r on dU whose outer side is
// uncovered along a set of positive length. Limits of uncovered points give
// points y on e with dist(y, K) >= r = |y - w| and y in Omega*, so the open disk
// D(y, r) lies in int Omega* and touches dOmega* at w. At a convex vertex w
// (interior angle < pi) the tangent cone of Omega* holds no such disk; at a
// smooth or G1 point w the disk's centre is y = w + r n_in, a single point of e.
// A reflex vertex would admit an arc of such y; Omega* has none, because it is
// an erosion: its complement is P^c (+) D_r, so every boundary point has an
// exterior r-disk, and the junction check above verifies it on the pieces. So
// no cap edge has an uncovered side in T.

#include "coverage_certificate_2/counters_2.h"
#include "coverage_certificate_2/exact_target_2.h"
#include "coverage_certificate_2/kernel_2.h"
#include "coverage_certificate_2/lines_2.h"
#include "coverage_certificate_2/regions_2.h"
#include "coverage_certificate_2/vertex_events_2.h"

#include <CGAL/Bbox_2.h>

#include <array>
#include <cstdint>
#include <vector>

namespace compas_cgal::coverage_certificate {

/// One sub-piece on the wire: part of exact-target piece `piece` from `start` to
/// `end` in travel order; `full` for a whole closed circle piece.
struct ExactConnectorPieceSpec {
    std::uint32_t piece;
    bool full;
    OneRootPoint start, end;
};

/// One exact connector: its sub-pieces in travel order (empty for a fan).
struct ExactConnectorSpec {
    std::vector<ExactConnectorPieceSpec> pieces;
};

/// One sub-piece as a composite swept region.
struct ConnectorRegion {
    std::uint32_t connector;  ///< exact-connector index
    std::uint32_t piece;      ///< exact-target piece index
    bool is_line;
    bool degenerate;  ///< zero length (start == end, not a full turn): the region is its end disk only
    // Line sub-piece: unit interior normal n, level, travel direction u = (n_y, -n_x).
    FT nx, ny, level, ux, uy;
    OneRoot s0, s1;                      ///< u.w0 <= u.w1
    std::array<std::uint32_t, 2> flank;  ///< line ids of n.x = level + r and n.x = level - r
    std::array<FT, 2> t_lo, t_hi;        ///< rational chart enclosure of [t(w0), t(w1)] on each flank
    // Arc sub-piece.
    Point_2 c;
    FT radius;  ///< support radius R
    int sense;
    ArcSpan span;
    bool full;
    std::uint32_t support_region;  ///< support-only CircleRegion of |x - c| = R +/- r
    // End disks.
    OneRootPoint w0, w1;
    std::array<std::int32_t, 2> cap_vertex;  ///< VertexCircle id (irrational centre) or -1
    std::array<std::int32_t, 2> cap_region;  ///< support-only CircleRegion (rational centre) or -1
    CGAL::Bbox_2 bbox;                       ///< certified outward box of the whole region
    // Precomputed exact squares and interval shadows (filters only).
    FT outer_sq, inner_sq;  ///< (R + r)^2, (R - r)^2 (inner only when R > r)
    bool annulus = false;   ///< R > r
    Interval ux_i, uy_i, nx_i, ny_i, level_i, s0_i, s1_i, cx_i, cy_i, outer_sq_i, inner_sq_i;
};

/// The connector regions and the supports they added.
struct ConnectorSet {
    /// One region per dOmega* piece (the swept set; see build), connector = the
    /// lowest connector with a sub-piece on it.
    std::vector<ConnectorRegion> regions;
    /// The connectors' own sub-pieces (validated; reporting and the cover check).
    std::vector<ConnectorRegion> sub_regions;
    std::vector<std::uint32_t> flank_lines;       ///< distinct flank line ids, sorted
    std::vector<std::uint8_t> cap_circle;         ///< per boundary circle of the set: 1 for a rational cap circle
    std::vector<std::uint32_t> support_reports;   ///< per support-only region: lowest connector index using it

    /// Build the regions on the exact target's pieces, register their supports
    /// (flank lines in `lines`, support circles and rational caps as support-only
    /// regions in `set`, irrational caps as vertex circles in `vertices`, shared
    /// with identical ones already there), and verify that the sub-pieces cover
    /// dOmega*. `mode` selects filtered (LAZY) or exact-only (EXACT_ONLY) signs
    /// in the arc ordering of the cover check.
    ///
    /// Raises:
    ///     InvalidCertificateInputError: a sub-piece names no target piece, an end
    ///         point is not exactly on the piece's support or the ends are out of
    ///         travel order.
    ///     ConnectorBoundaryCoverageError: the sub-pieces do not cover dOmega*.
    [[nodiscard]] static ConnectorSet build(const std::vector<ExactConnectorSpec>& connectors, const ExactTarget& target,
                                            CircleRegionSet& set, LineTable& lines, std::vector<VertexCircle>& vertices,
                                            FilterMode mode);
};

/// A vertex circle |x - w| = r for an irrational one-root centre w.
[[nodiscard]] VertexCircle make_vertex_circle(const OneRoot& x, const OneRoot& y, const FT& r, std::uint32_t disk);

/// Exact closed side of p (coordinates in one field Q(sqrt(L))) relative to the
/// composite region (interior of any member: Interior; else on any member:
/// Boundary; else Exterior).
[[nodiscard]] Side connector_membership(const ConnectorRegion& region, const ExactPoint& p, const FT& r, FilterMode mode,
                                        DecisionCounters& counters);

/// Whether p lies strictly inside the sub-piece's range: s0 < u.p < s1 for a line,
/// the open wedge (anything, for a full turn) for an arc.
[[nodiscard]] bool connector_in_range(const ConnectorRegion& region, const ExactPoint& p, FilterMode mode);

}  // namespace compas_cgal::coverage_certificate
