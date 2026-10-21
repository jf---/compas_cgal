#pragma once

// Exact nearest-point queries on the boundary of the Held centre domain Omega*
// (docs/held_exact_centre_domain.md, "Projection of a hypothesis contact").
//
// Every distance is a one-root number over the rational supports: rational to
// an offset line, (sqrt(D) - R)^2 with D = |p - c|^2 to a circle piece in its
// wedge, and one-root in the vertex's own root to a vertex. Signs of one-root
// and of two-radical sums a + b*sqrt(alpha) + c*sqrt(beta) are decided over
// the rational parts by repeated squaring (docs/exactness.md, "Numeric
// comparison is the exact-kernel idiom"); no root of one extension is ever
// combined arithmetically with another.
//
// Filtering is CGAL's static-filter idiom, not a tolerance: the whole query
// first runs on certified intervals (Interval_nt, outward rounding), where
// every sign that the enclosure cannot decide throws
// Uncertain_conversion_exception; the query is then re-run on the exact
// rationals. A certified double box per piece only prunes: a piece is skipped
// only when a rigorous lower bound on its squared distance exceeds a rigorous
// upper bound on a distance already attained.

#include "held_centre_domain_2_kernel.h"

#include <CGAL/Interval_nt.h>

#include <cstddef>
#include <utility>
#include <vector>

namespace compas_cgal::held_centre_domain {

using Interval = CGAL::Interval_nt<false>;
/// The exact rational type behind FT, used directly by the exact stage of the
/// queries (no lazy DAG per operation).
using ExactRational = FT::ET;

/// a + b * sqrt(r) with b = r = 0 for a rational value.
template <typename NT>
struct OneRootT {
    NT a;
    NT b;
    NT r;
};
using OneRoot = OneRootT<FT>;

/// The numbers of one piece, in one number type.
template <typename NT>
struct PieceNumbers {
    NT a, b, c;         // line a x + b y + c = 0
    NT norm_squared;    // a^2 + b^2
    NT dx, dy;          // line travel direction (a multiple of (-b, a))
    NT cx, cy, radius;  // circle, rational radius
    OneRootT<NT> sx, sy, tx, ty;
};

/// One maximal boundary piece. Its source is vertex `index` (its own global
/// index); its target is vertex `target_vertex` (the next piece's source).
struct QueryPiece {
    bool circular = false;
    bool full = false;  // a full circle: source == target
    int sense = 0;      // circle travel sense, +1 CCW / -1 CW
    int target_half = 0;  // half turn of the target from the source: 0 in (0, pi), 1 in [pi, 2pi)
    std::size_t target_vertex = 0;
    PieceNumbers<FT> exact;
    PieceNumbers<ExactRational> rational;
    PieceNumbers<Interval> approximate;
    double xmin = 0, ymin = 0, xmax = 0, ymax = 0;  // certified outward box
};

/// Fill `rational` and `approximate` from `exact` (exact values, certified enclosures) and decide the
/// query-independent `target_half` exactly.
void finish(QueryPiece& piece);

enum class NearestKind : int { Vertex = 0, Interior = 1 };

/// The decision of a query: which point is nearest, and whether another,
/// distinct point is exactly as near. A vertex is reported on the piece whose
/// source it is.
struct NearestDecision {
    bool ambiguous;
    NearestKind kind;
    std::size_t piece;
    bool beyond_gap;  // only meaningful when a gap was given
};

/// Exact one-root values of a decided nearest point.
struct NearestValues {
    OneRoot squared_distance;
    OneRoot x;
    OneRoot y;
};

[[nodiscard]] int sign(const OneRoot& value);
/// Exact sign of `first - second`, also across different roots.
[[nodiscard]] int compare(const OneRoot& first, const OneRoot& second);

/// Rigorous enclosure [low, high] of a one-root value (filters and boxes only).
[[nodiscard]] std::pair<double, double> enclosure(const OneRoot& value);

/// Scratch storage reused across queries, and the count of queries the
/// interval stage could not decide (reported, never a decision input).
struct QueryWorkspace {
    std::vector<double> lower;
    std::vector<std::size_t> seen_vertices;
    std::size_t exact_fallbacks = 0;
};

/// Decide the nearest point of the union of `pieces` to the rational point
/// `query`; with `gap_squared`, also whether it is farther than the gap.
[[nodiscard]] NearestDecision decide_nearest(const std::vector<QueryPiece>& pieces, const RationalPoint& query,
                                             const FT* gap_squared, QueryWorkspace& workspace);

/// The exact squared distance and coordinates of a decided nearest point.
[[nodiscard]] NearestValues nearest_values(const std::vector<QueryPiece>& pieces, const RationalPoint& query,
                                           const NearestDecision& decision);

}  // namespace compas_cgal::held_centre_domain
