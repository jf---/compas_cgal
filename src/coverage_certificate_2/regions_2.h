#pragma once

// Swept regions of full-circle cuts and their boundary pieces.
//
// A full-circle cut with path radius R about centre c, tool radius r, sweeps
// the closed annulus max(R - r, 0) <= |x - c| <= R + r (a solid disk when
// R <= r). Its boundary is one or two circles with rational centre and rational
// radius, so every piece is a Circular_kernel_2 circle and every stereographic
// witness point on it is rational.
//
// Phase 2: a connector capsule is decomposed into two disks D_r(p), D_r(q) and a
// rectangle (lines_2.h). The disks are ordinary circle regions here
// (PieceKind::CapsuleCap, radius r, rational centre). The target's vertex disks
// D_r(v) contribute boundary circles too (PieceKind::TargetCircle): they are
// PIECES checked by the target-edge rule, never swept regions, so their `owner`
// is a target-member index and `target` is set.

#include "coverage_certificate_2/counters_2.h"
#include "coverage_certificate_2/kernel_2.h"

#include <CGAL/Bbox_2.h>

#include <array>

#include <cstdint>
#include <vector>

namespace compas_cgal::coverage_certificate {

/// Which boundary piece of a full-circle cut a circle is.
enum class PieceKind : std::uint8_t {
    AnnulusOuter = 0,  ///< outer circle of an annulus (R > r); region inside
    AnnulusInner = 1,  ///< inner circle of an annulus; region outside
    Disk = 2,          ///< boundary of a solid disk (R <= r); region inside
    CapsuleCap = 3,    ///< boundary of a connector end disk D_r(p); region inside
    TargetCircle = 4,  ///< boundary of a target vertex disk D_r(v); a target piece, not a sweep
};

/// Exact location of a point relative to a closed region.
enum class Side : std::uint8_t { Exterior = 0, Boundary = 1, Interior = 2 };

/// The swept region of one full-circle cut, exact and with interval shadows.
struct CircleRegion {
    Point_2 centre;
    FT outer_squared;  ///< (R + r)^2
    FT inner_squared;  ///< (R - r)^2 for an annulus, unused for a disk
    bool annulus;      ///< R > r, decided exactly
    Interval centre_x, centre_y, outer_squared_i, inner_squared_i;
    std::uint32_t report_index;  ///< full-circle index, or n_full + connector index for a capsule cap
    PieceKind kind;              ///< kind of the region's outer boundary circle
    /// False for a SUPPORT-ONLY region (exact connectors, exact_connectors_2.h):
    /// its circles cut every piece they meet and may be checked as pieces, but it
    /// never covers anything. Its membership is never consulted for coverage.
    bool covers = true;
};

/// One boundary circle (piece) of a swept region.
struct BoundaryCircle {
    std::uint32_t owner;  ///< index of the CircleRegion, or the target member for a TargetCircle
    PieceKind kind;
    Circle_2 circle;
    FT radius;  ///< exact rational radius (R + r or R - r)
    Interval centre_x, centre_y, radius_squared_i;
    CGAL::Bbox_2 bbox;  ///< certified outward-rounded box of the circle

    [[nodiscard]] bool is_inner() const { return kind == PieceKind::AnnulusInner; }
    [[nodiscard]] bool is_target() const { return kind == PieceKind::TargetCircle; }
    /// The side of the circle whose cover the certificate checks: the outer side
    /// of a swept region's boundary (inside the circle for an inner annulus
    /// circle), the member side of a target piece (inside the vertex disk).
    [[nodiscard]] bool check_inside() const { return is_inner() || is_target(); }
    /// For a swept region's circle: whether the region lies inside it.
    [[nodiscard]] bool region_inside() const { return !is_inner(); }
};

/// Regions and their boundary circles, in input order: region i owns circle
/// indices [first_circle[i], first_circle[i + 1]).
struct CircleRegionSet {
    std::vector<CircleRegion> regions;
    std::vector<BoundaryCircle> circles;
    std::vector<std::uint32_t> first_circle;

    /// Inject the cuts exactly and build their pieces.
    ///
    /// Raises:
    ///     InvalidCertificateInputError: shapes differ, a radius is not positive.
    ///     NonFiniteBinary64Error: a non-finite coordinate or radius.
    [[nodiscard]] static CircleRegionSet build(const std::vector<std::array<double, 2>>& centres,
                                               const std::vector<double>& path_radii, double tool_radius);

    /// Append connector end disks D_r(p) as circle regions (PieceKind::CapsuleCap).
    void add_capsule_caps(const std::vector<Point_2>& centres, const std::vector<std::uint32_t>& report_index,
                          const FT& tool_radius);

    /// Append the target vertex circles |x - v| = r as target pieces (no region).
    void add_target_circles(const std::vector<Point_2>& centres, const std::vector<std::uint32_t>& member,
                            const FT& tool_radius);

    /// Append one target circle |x - c| = radius (a cut support, no region);
    /// returns its circle id.
    std::uint32_t add_target_circle(const Point_2& centre, const FT& radius, std::uint32_t member);

    /// Append a support-only region (covers == false) with the circles
    /// |x - c| = outer and, when inner > 0, |x - c| = inner; returns its index.
    ///
    /// Raises:
    ///     InvalidCertificateInputError: outer <= inner or inner < 0.
    std::uint32_t add_support_region(const Point_2& centre, const FT& outer, const FT& inner, std::uint32_t report_index);

    /// Number of circle regions (cuts + capsule caps); target circles own none.
    [[nodiscard]] std::size_t region_count() const { return regions.size(); }
};

/// Exact side of the rational point q relative to the closed region (P4 for a
/// full-circle cut): compares |q - c|^2 with (R + r)^2 and (R - r)^2.
[[nodiscard]] Side membership(const CircleRegion& region, const Point_2& q, const Interval& qx, const Interval& qy,
                              FilterMode mode, DecisionCounters& counters);

/// The rational stereographic point of circle C at parameter s,
/// c + rho ((1 - s^2), 2 s) / (1 + s^2). Every rational s gives a rational point,
/// and s increasing from -inf to +inf is the CCW order from the leftmost point
/// (s = +/-inf), so s < 0 is the lower half and s >= 0 the upper half.
[[nodiscard]] Point_2 stereographic_point(const BoundaryCircle& circle, const FT& s);

/// 0 (lower half) when s < 0, else 1.
[[nodiscard]] inline std::uint8_t half_of_stereographic(const FT& s)
{
    return CGAL::sign(s) == CGAL::NEGATIVE ? 0 : 1;
}

}  // namespace compas_cgal::coverage_certificate
