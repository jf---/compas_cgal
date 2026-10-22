#pragma once

// Binary64 machining circles on exact sites of the Held centre domain's
// boundary (docs/held_exact_centre_domain.md, "Native site and placement").
//
// A decided contact (binary64 probe, piece, vertex flag) names an exact site q
// on dOmega*: the piece's source vertex, or the orthogonal projection of the
// probe onto the piece's support. The tangent circle of guide radius t is
// centred at c* = q + t*n (n the inward unit normal): rational on a line,
// C + (rho' - sigma*t) * w / sqrt(|w|^2) on a circle (w = probe - C, one
// root), and in the vertex's own root at a vertex. Every decision is exact:
//
// - the site's placement (source or interior) and the outside-the-piece
//   verdict are exact signs over the declared supports and the one-root
//   vertices, filtered by certified intervals (CGAL::Interval_nt, outward
//   rounding) and decided exactly whenever the enclosure is uncertain;
// - the centre is the correctly rounded (nearest, ties to even) binary64 of
//   each exact coordinate, found by comparing the exact value with the dyadic
//   midpoints of neighbouring doubles: sign((a - m) + b sqrt(r)) decomposed
//   as Sqrt_extension::compare does, cross-multiplied in integers (the exact
//   rational type reduces by gcd on every operation, which dominated);
// - the radius is the largest binary64 <= min(t, bound), where bound is the
//   exact clearance of the rounded centre from the piece's whole support
//   (n.c - level on a line, sqrt(D) -/+ rho' on a circle): the guide radius
//   is kept when it fits, else the bound is rounded down exactly.
//
// Doubles enter once by exact injection (every finite binary64 is a dyadic
// rational); `nextafter` only enumerates the neighbouring candidates whose
// midpoints are compared exactly. No tolerance anywhere.
//
// Every line piece of dOmega* lies on an r-offset of a declared line, so it
// carries that line's exactly unit normal; the frame keeps the support in that
// declared form (n, level), whose rationals are small, and the piece's ends in
// their one-root form for the range decisions.

#include "held_centre_domain_query_2.h"

#include <CGAL/Fraction_traits.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Sqrt_extension.h>

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace compas_cgal::held_centre_domain {

/// A line piece of dOmega* lies on no offset of a declared line, so it has no
/// exactly rational unit normal (a construction contract breach).
class NoRationalUnitNormalError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/// a0 + a1 * sqrt(root) over the exact rationals: CGAL's one-root number.
/// Arithmetic only within one extension (or with a rational); comparison with
/// a rational is exact.
using PlacementRoot = CGAL::Sqrt_extension<ExactRational, ExactRational, CGAL::Tag_true, CGAL::Tag_true>;
/// Certified enclosures of the range filters (protected rounding).
using PlacementInterval = CGAL::Interval_nt<true>;
/// The integer type behind the exact rationals (numerator and denominator).
using PlacementInteger = CGAL::Fraction_traits<ExactRational>::Numerator_type;

/// A one-root value with its certified enclosure.
struct EnclosedRoot {
    PlacementRoot value;
    PlacementInterval box;
};

/// Query-independent placement data of one piece.
///
/// Line: inward unit normal n, `level` = n.x on the support, travel unit
/// u = (n_y, -n_x), and the ends' positions u.s, u.t. Circle: rational centre,
/// radius and travel sense (from the query piece), and the ends' vectors from
/// the centre. Both: the source vertex.
struct PlacementFrame {
    bool circular = false;
    bool has_normal = false;
    ExactRational nx, ny, level;
    PlacementInterval ux_box, uy_box;
    EnclosedRoot source_along, target_along;
    EnclosedRoot source_x, source_y;            // circle: s - C; line: unused
    EnclosedRoot target_x, target_y;            // circle: t - C
    PlacementInterval centre_x_box, centre_y_box;
    double centre_x_view = 0.0, centre_y_view = 0.0, radius_view = 0.0;  // estimates only, never decide
    // Circle: centre and radius as integer numerator / denominator pairs.
    PlacementInteger centre_x_n, centre_x_d, centre_y_n, centre_y_d, radius_n, radius_d;
    PlacementRoot vertex_x, vertex_y;           // the source vertex s
};

/// Placement frames of every piece. A line piece's normal is the unit normal
/// of a parallel declared line (all declared normals are exactly unit, so
/// parallel ones agree up to sign), oriented to the left of its travel.
[[nodiscard]] std::vector<PlacementFrame> placement_frames(const std::vector<QueryPiece>& pieces,
                                                           const std::vector<Kernel::Vector_2>& declared_normals);

/// Verdict of one placement (held_centre_domain_native.py mirrors the codes).
enum class PlacementStatus : int {
    Placed = 0,
    /// The projection of the probe is not strictly before the piece's target.
    OutsidePiece = 1,
    /// The probe is a circle piece's centre and names no direction.
    ProbeAtCircleCentre = 2,
};

enum class SitePlacement : int { Source = 0, Interior = 1 };

/// One placed circle: the exact site's placement (and half turn from the
/// piece's source on a circle), the correctly rounded centre, the radius, and
/// the binary64 view of the inward normal (reporting geometry: the contact
/// point is centre - radius * normal_view).
struct PlacedCircle {
    PlacementStatus status = PlacementStatus::Placed;
    SitePlacement placement = SitePlacement::Interior;
    int half = 0;
    double cx = 0.0;
    double cy = 0.0;
    double radius = 0.0;
    bool limited = false;
    double normal_x = 0.0;
    double normal_y = 0.0;
};

/// Place the tangent circle of guide radius `guide_radius` at the site of the
/// decided contact (`qx`, `qy`) on `piece`.
///
/// Raises:
///     NoRationalUnitNormalError: a line piece without a rational unit normal.
///     InvalidCentreDomainInputError: a non-finite input or a non-positive guide radius.
[[nodiscard]] PlacedCircle place_on_normal(const QueryPiece& piece, const PlacementFrame& frame, double qx, double qy,
                                           bool vertex, double guide_radius);

/// Correctly rounded binary64 nearest an exact rational (ties to even, +0 for 0).
[[nodiscard]] double nearest_binary64(const ExactRational& value);
/// Largest binary64 <= an exact rational.
[[nodiscard]] double round_down_binary64(const ExactRational& value);

}  // namespace compas_cgal::held_centre_domain
