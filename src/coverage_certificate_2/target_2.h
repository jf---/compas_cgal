#pragma once

// The machinable target T and its members (docs/held_coverage_certificate.md).
//
// T = Omega (+) D_r is decomposed into simple closed members:
//   * POLYGON: the centre domain Omega itself (outer loop CCW, holes CW);
//   * EDGE_CAPSULE: dist(x, e_k) <= r for every loop edge e_k (holes included);
//   * VERTEX_DISK: D_r(w) about a convex vertex w of the exact centre domain
//     Omega*, whose coordinates are one-root numbers in Q(sqrt(delta)).
//
// Phase 1 decides MEMBERSHIP of rational points in every member exactly (P4).
// Phase 2 cuts the target into the arrangement as the edge rectangles and
// vertex disks of the capsule decomposition (lines_2.h, regions_2.h) and decides
// target sides with target_decision_2.h. A VERTEX_DISK (one-root centre, the
// exact centre domain Omega*) still enters membership only: its boundary is not
// cut into the arrangement yet, so with vertex disks present the certificate can
// prove non-coverage but not coverage.

#include "coverage_certificate_2/counters_2.h"
#include "coverage_certificate_2/kernel_2.h"

#include <array>
#include <cstdint>
#include <vector>

namespace compas_cgal::coverage_certificate {

/// The kinds of target member, in region-index order after the sweeps.
enum class TargetMemberKind : std::uint8_t { Polygon = 0, EdgeCapsule = 1, VertexDisk = 2 };

/// A validated centre domain: loop 0 is the outer loop (CCW), the rest are
/// holes (CW), every loop simple, holes strictly inside and pairwise disjoint.
struct CentreDomain {
    std::vector<std::vector<Point_2>> loops;

    /// Raises:
    ///     InvalidCentreDomainError: any invariant above fails, a loop has fewer
    ///         than three vertices, or an edge has zero length.
    ///     NonFiniteBinary64Error: a non-finite coordinate.
    [[nodiscard]] static CentreDomain build(const std::vector<std::array<double, 2>>& outer,
                                            const std::vector<std::vector<std::array<double, 2>>>& holes);
};

/// D_r(w) for a one-root centre w = (x0 + x1*sqrt(delta), y0 + y1*sqrt(delta)).
struct VertexDisk {
    OneRoot x;
    OneRoot y;
    Interval x_i, y_i;

    /// Raises:
    ///     InvalidCertificateInputError: delta is negative.
    [[nodiscard]] static VertexDisk build(const FT& x0, const FT& x1, const FT& y0, const FT& y1, const FT& delta);
};

/// The machinable target T for one tool radius.
struct MachinableTarget {
    CentreDomain domain;
    FT tool_radius_squared;
    Interval tool_radius_squared_i;
    std::vector<VertexDisk> vertex_disks;
    std::vector<std::array<Interval, 4>> edge_boxes;  ///< certified (xmin, xmax, ymin, ymax) per edge, loop-major
    std::vector<std::vector<std::array<Interval, 2>>> vertex_intervals;  ///< interval shadows of the loop vertices
    FT tool_radius;                                    ///< r, exact

    [[nodiscard]] static MachinableTarget build(CentreDomain domain, double tool_radius,
                                                std::vector<VertexDisk> vertex_disks);

    /// Exact closed membership of the rational point q in T (P4 over members).
    /// Phase-1 scan over every member; the phase-2 certificate uses the
    /// box-indexed candidates of target_decision_2.h instead.
    [[nodiscard]] bool contains(const Point_2& q, FilterMode mode, DecisionCounters& counters) const;
};

/// sign(|q - w|^2 - r^2) for a one-root vertex-disk centre w and rational q.
[[nodiscard]] CGAL::Sign vertex_disk_sign(const VertexDisk& disk, const Point_2& q, const FT& r_sq);

/// Exact closed membership of a point with one-root coordinates (one shared
/// field) in the polygon with holes (crossing number; interval filter first).
struct ExactPoint;
[[nodiscard]] bool in_closed_domain(const CentreDomain& domain, const std::vector<std::vector<std::array<Interval, 2>>>& shadows,
                                    const ExactPoint& q, FilterMode mode, DecisionCounters& counters);

}  // namespace compas_cgal::coverage_certificate
