#pragma once

// Exact Held centre domain Omega* = P minus (dP + D_r) on the declared
// (Pythagorean) boundary (docs/held_exact_centre_domain.md, "Construction").
//
// Every input is rational: line segments with rational ends, arcs on rational
// circles with rational ends, a rational tool radius. Every member of
// dP + D_r (line band, annular-sector arc band, junction disk) then has
// rational supports, so every vertex of the arrangement is a one-root point
// and every decision stays inside Gps_circle_segment_traits_2<Epeck>'s
// documented one-root contract. CORE-free: this file never defines
// CGAL_USE_CORE and links no exact core library.

#include "held_centre_domain_2_kernel.h"
#include "held_centre_domain_query_2.h"

#include <array>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace compas_cgal::held_centre_domain {

/// Malformed input: a non-positive radius or denominator, an unknown kind, a
/// zero-length primitive, an arc whose ends are not exactly on its circle.
class InvalidCentreDomainInputError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

/// Module tag of `_held_centre_domain_2` for the templated Gps CCB repair
/// errors (gps_ccb_repair_2.h).
struct CentreDomainCcbRepairDomain {};

/// The arrangement returned a boundary cycle that does not chain exactly
/// (target of one halfedge != source of the next). A contract breach of the
/// extraction, never expected.
class BrokenCentreDomainChainError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/// Signed decimal numerator and positive decimal denominator (the Python
/// `Fraction` wire format).
struct RationalText {
    std::string numerator;
    std::string denominator;
};

/// a0 + a1 * sqrt(root); `a1 == root == 0` for a rational value.
struct OneRootText {
    RationalText rational;
    RationalText coefficient;
    RationalText radicand;
};

struct PointText {
    OneRootText x;
    OneRootText y;
};

enum class PrimitiveKind : int { Line = 0, Arc = 1 };

/// One declared primitive, travelled from `start` to `end` with the pocket on
/// its left. A line carries its exactly unit interior (left) normal. An arc
/// carries its declared rational `centre`, `radius` and `sense` (+1 CCW,
/// -1 CW), with |start - centre| = |end - centre| = radius exactly.
struct Primitive {
    PrimitiveKind kind;
    RationalPoint start;
    RationalPoint end;
    Kernel::Vector_2 normal;
    RationalPoint centre;
    FT radius;
    int sense;
};

/// One maximal curve of a boundary cycle of Omega* on a single support.
///
/// Line support: `support = (a, b, c)` with a x + b y + c = 0, as CGAL stores
/// it (not normalised). Circle support: `support = (cx, cy, R^2)`; `sense` is
/// the travel sense along the cycle (+1 CCW, -1 CW), 0 for a line.
struct Piece {
    bool circular;
    std::array<RationalText, 3> support;
    int sense;
    PointText source;
    PointText target;
};

using Cycle = std::vector<Piece>;

/// One face of Omega*: the outer cycle (CCW) first, then its holes (CW). The
/// face lies to the left of every piece.
struct Component {
    std::vector<Cycle> cycles;
};

enum class Location : int { Exterior = -1, Boundary = 0, Interior = 1 };

class CentreDomain2 {
public:
    /// Construct Omega* = P minus join(line bands, arc bands, junction disks).
    ///
    /// Raises:
    ///     InvalidCentreDomainInputError: malformed primitives or radius.
    ///     BrokenCentreDomainChainError: an extracted cycle does not chain.
    static std::unique_ptr<CentreDomain2> build(const std::vector<Primitive>& primitives, const FT& tool_radius);

    [[nodiscard]] const std::vector<Component>& components() const { return components_; }
    [[nodiscard]] std::size_t member_count() const { return member_count_; }
    [[nodiscard]] double construction_seconds() const { return construction_seconds_; }
    /// Outer CCBs the Gps difference misfiled and the build re-filed as holes
    /// (gps_ccb_repair_2.h); reporting only.
    [[nodiscard]] std::size_t refiled_outer_ccbs() const { return refiled_outer_ccbs_; }

    /// Exact location of a rational point relative to the closed set Omega*.
    [[nodiscard]] Location locate(const RationalPoint& point) const;

    /// Exact nearest boundary point of a rational point (pieces numbered in
    /// component, cycle, CCW order, as `components()` lists them); with
    /// `gap_squared`, also whether it is farther than the gap.
    [[nodiscard]] NearestDecision decide_nearest(const RationalPoint& point, const FT* gap_squared,
                                                 QueryWorkspace& workspace) const;
    [[nodiscard]] NearestValues nearest_values(const RationalPoint& point, const NearestDecision& decision) const;

private:
    CentreDomain2() = default;

    // Owned in place and never copied, so the arrangement's traits pointer is
    // the set's own (docs/exactness.md, checklist 18).
    std::unique_ptr<Gps> domain_;
    std::vector<Component> components_;
    std::vector<QueryPiece> query_pieces_;
    std::size_t member_count_ = 0;
    double construction_seconds_ = 0.0;
    std::size_t refiled_outer_ccbs_ = 0;
};

[[nodiscard]] FT rational_from_text(const std::string& numerator, const std::string& denominator);
[[nodiscard]] RationalText rational_to_text(const FT& value);

}  // namespace compas_cgal::held_centre_domain
