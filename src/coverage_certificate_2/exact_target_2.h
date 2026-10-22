#pragma once

// The exact machinable target T = Omega* (+) D_r of the exact centre domain
// (docs/held_exact_centre_domain.md, "Coverage target and the exact path";
// docs/held_coverage_certificate.md, "Exact target Omega*").
//
// Omega* = {x in P : dist(x, dP) >= r} for the declared boundary P (rational
// lines with Pythagorean unit directions, rational arcs). Its boundary pieces lie
// on rational supports (offset lines, offset circles, junction r-circles) and
// end at one-root vertices w. T decomposes exactly into closed members:
//
//   T = Omega*  u  (band of every dOmega* piece)  u  (D_r(w), w a convex vertex),
//
// a band being the normal segments of length r from the piece toward the
// primitive it offsets: a rectangle {level - r <= n.x <= level, u.w0 <= u.x <= u.w1}
// for an offset line, an annular sector between radii R and rho (rho = 0 for a
// junction r-circle) inside the wedge of the piece for a circle. (Nearest-point
// argument: a point of T \ Omega* has its nearest Omega* point in a piece
// interior, then it lies in that piece's band, or at a vertex, whose normal cone
// is non-trivial only at a convex vertex, then it lies in D_r(w).)
//
// dT lies on the declared boundary supports and the vertex circles |x - w| = r:
// a point of dT inside P has its nearest Omega* point at distance exactly r, at a
// convex vertex w. These are the certificate's target cut supports.
//
// Closed membership of a point with coordinates in one field Q(sqrt(L)) (a
// rational witness, or a one-root flank witness) is decided exactly: every
// member test is a sign in Q(sqrt(L))(sqrt(delta)) for the one radicand delta of
// the member's one-root vertex, i.e. a tower over Q(sqrt(L)) (tower_2.h). Each
// test runs an interval filter first.

#include "coverage_certificate_2/counters_2.h"
#include "coverage_certificate_2/kernel_2.h"
#include "coverage_certificate_2/lines_2.h"
#include "coverage_certificate_2/target_2.h"
#include "coverage_certificate_2/tower_2.h"

#include <CGAL/Bbox_2.h>

#include <array>
#include <cstdint>
#include <utility>
#include <vector>

namespace compas_cgal::coverage_certificate {

/// A point whose coordinates share one radicand delta (0: rational), folded.
struct OneRootPoint {
    OneRoot x;
    OneRoot y;
    FT delta;  ///< the common radicand, 0 when both coordinates are rational
    Interval x_i, y_i;

    /// Raises:
    ///     InvalidCertificateInputError: delta is negative.
    [[nodiscard]] static OneRootPoint build(const FT& x0, const FT& x1, const FT& y0, const FT& y1, const FT& delta);
    [[nodiscard]] bool is_rational() const { return !x.is_extended() && !y.is_extended(); }
};

/// The sense of an arc's span: minor (< pi), half (== pi) or major (> pi).
enum class ArcSpan : std::uint8_t { Minor = 0, Half = 1, Major = 2 };

/// A declared boundary line segment a -> b with rational unit direction u; the
/// pocket is on its left (normal (-u_y, u_x)).
struct DeclaredSegment {
    Point_2 a, b;
    FT ux, uy, length;
    Interval ax_i, ay_i, bx_i, by_i;
    CGAL::Bbox_2 bbox;
};

/// A y-monotone sub-arc of a declared arc (rational ends), on the left
/// (x <= c_x) or right side of its circle; for the crossing-number test.
struct MonotoneArc {
    Point_2 p, q;
    int side;  ///< -1 left half, +1 right half
};

/// A declared boundary arc from a to b about rational centre c, radius rho, in
/// sense +1 (CCW) or -1.
struct DeclaredCircleArc {
    Point_2 c, a, b;
    FT rho;
    int sense;
    ArcSpan span;
    std::vector<MonotoneArc> monotone;
    Interval cx_i, cy_i, rho_i;
    CGAL::Bbox_2 bbox;
};

/// The band of an offset-line piece of dOmega*.
struct LineBand {
    FT nx, ny;          ///< interior unit normal of the declared line
    FT level;           ///< offset line n.x = level; declared line n.x = level - r
    FT ux, uy;          ///< travel direction (unit, rational)
    OneRoot s0, s1;     ///< u.w0 <= u.x <= u.w1 (one-root, radicands of w0 and w1)
    CGAL::Bbox_2 bbox;  ///< certified outward box
};

/// The band of a circle piece of dOmega* (offset arc or junction arc).
struct CircleBand {
    Point_2 c;
    FT radius;              ///< the piece's own support radius R (the tool-centre circle)
    FT inner_sq, outer_sq;  ///< squared radii min(R, rho)^2, max(R, rho)^2
    int sense;              ///< travel sense from w0 to w1
    OneRootPoint w0, w1;
    ArcSpan span;
    Interval cx_i, cy_i, inner_sq_i, outer_sq_i;
    CGAL::Bbox_2 bbox;
};

/// One dOmega* piece as it arrives from the wire.
struct ExactPieceSpec {
    bool is_line;
    // Line: interior unit normal (nx, ny), level, direction (ux, uy).
    FT nx, ny, level, ux, uy;
    // Circle: centre, piece radius R, boundary radius rho (0 for a junction), sense.
    Point_2 centre;
    FT radius, boundary_radius;
    int sense;
    OneRootPoint source, target;
};

/// One declared primitive as it arrives from the wire.
struct DeclaredPrimitiveSpec {
    bool is_line;
    Point_2 start, end;
    FT ux, uy, length;  ///< line only
    Point_2 centre;     ///< arc only
    FT radius;          ///< arc only
    int sense;          ///< arc only
};

/// The exact target.
struct ExactTarget {
    FT tool_radius, tool_radius_squared;
    Interval tool_radius_squared_i;
    std::vector<DeclaredSegment> segments;
    std::vector<DeclaredCircleArc> arcs;
    std::vector<LineBand> line_bands;
    std::vector<CircleBand> circle_bands;
    /// Per dOmega* piece in input order: (is_line, index into line_bands or circle_bands).
    std::vector<std::pair<bool, std::uint32_t>> piece_bands;
    /// Per dOmega* piece in input order: its source and target vertices.
    std::vector<OneRootPoint> piece_sources, piece_targets;
    std::vector<VertexDisk> vertex_disks;
    std::vector<CGAL::Bbox_2> vertex_boxes;

    /// Raises:
    ///     InvalidCertificateInputError: a non-positive radius, a degenerate
    ///         primitive, a piece whose band is empty, or no primitive.
    [[nodiscard]] static ExactTarget build(double tool_radius, const std::vector<DeclaredPrimitiveSpec>& primitives,
                                           const std::vector<ExactPieceSpec>& pieces, std::vector<VertexDisk> vertex_disks);

    /// Exact closed membership of p in T.
    [[nodiscard]] bool contains(const ExactPoint& p, FilterMode mode, DecisionCounters& counters) const;

    /// Exact closed membership of p in Omega* (p in P and dist(p, dP) >= r).
    [[nodiscard]] bool in_centre_domain(const ExactPoint& p, FilterMode mode, DecisionCounters& counters) const;

    /// A rational point in the interior of T (strictly inside a line band).
    ///
    /// Raises:
    ///     InvalidCertificateInputError: dOmega* has no offset-line piece.
    [[nodiscard]] Point_2 interior_anchor(DecisionCounters& counters) const;
};

}  // namespace compas_cgal::coverage_certificate
