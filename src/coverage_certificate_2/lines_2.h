#pragma once

// Line supports and rectangle regions of the capsule decomposition
// (docs/held_coverage_certificate.md, "Capsules as rectangles and disks").
//
// A capsule dist(x, [p, q]) <= r is the union of two disks D_r(p), D_r(q) and the
// rectangle R(p, q) = {0 <= tau <= 1, |n.(x - p)| <= r sqrt(L)}, n = (q - p)^perp,
// L = |q - p|^2. The certificate's theorem holds for any finite family of
// regular-closed regions, so each capsule enters as those three regions. A
// rectangle is the intersection of four closed half-planes whose supports are
// UNRESTRICTED lines, so its membership along a piece is the conjunction of four
// side bits that flip exactly at transversal crossings of the lines: there are
// no junction points and no restricted-segment conventions.
//
// Every connector piece and every target edge capsule lies on an edge e of the
// centre domain. Their rectangles are written in that edge's own frame, so all
// flanks along one edge are the SAME line object, the target flank included:
// coincident supports are recognised by identity, and the exact coincidence test
// below covers the remaining cases (different edges on one line).
//
// A line is n.x = k0 + k1 sqrt(L) with rational n, k0, k1, L; a flank of a
// rational-length edge (L a rational square) is folded to a rational line, so a
// stored k1 != 0 always comes with an irrational sqrt(L). Its chart is
// x(t) = b0 + b1 sqrt(L) + t d with d = (n_y, -n_x), b1 = (k1 / |n|^2) n, so
// d . b1 = 0; on a flank t is the edge parameter.

#include "coverage_certificate_2/counters_2.h"
#include "coverage_certificate_2/kernel_2.h"
#include "coverage_certificate_2/target_2.h"
#include "coverage_certificate_2/tower_2.h"

#include <CGAL/Bbox_2.h>

#include <array>
#include <cstdint>
#include <map>
#include <tuple>
#include <vector>

namespace compas_cgal::coverage_certificate {

/// An exact point whose coordinates share one field Q(sqrt(L)) (rational when
/// both are non-extended): rational witnesses on circles, one-root witnesses on
/// flanks.
struct ExactPoint {
    OneRoot x;
    OneRoot y;

    [[nodiscard]] static ExactPoint rational(const Point_2& p) { return {OneRoot(p.x()), OneRoot(p.y())}; }
};

/// A line n.x = k0 + k1 sqrt(L) with its chart and interval shadows.
struct LineSupport {
    FT nx, ny;
    FT k0, k1, L;  ///< L == 0 and k1 == 0 for a rational line
    FT b0x, b0y, b1x, b1y;
    FT nn;  ///< |n|^2 = |d|^2
    Interval nx_i, ny_i, k0_i, k1_i, sqrt_l_i, offset_i, b0x_i, b0y_i, b1x_i, b1y_i, nn_i;

    [[nodiscard]] FT dx() const { return ny; }
    [[nodiscard]] FT dy() const { return -nx; }
    [[nodiscard]] bool is_rational() const { return CGAL::sign(k1) == CGAL::ZERO; }

    /// The line through the rational point `through` with normal n, shifted by
    /// `shift` * sqrt(radicand) along n / |n|... in the stored form; folds a
    /// rational-square radicand. Used by the flank and end-line factories.
    [[nodiscard]] static LineSupport build(const FT& nx, const FT& ny, const Point_2& through, const FT& k1,
                                           const FT& radicand);
};

/// Exact side of a point relative to a line: sign(n.p - k0 - k1 sqrt(L)).
[[nodiscard]] CGAL::Sign line_side(const LineSupport& line, const ExactPoint& p);

/// Filtered side of a rational point; exact when the enclosure straddles 0.
[[nodiscard]] CGAL::Sign line_side(const LineSupport& line, const Point_2& q, const Interval& qx, const Interval& qy,
                                   FilterMode mode, DecisionCounters& counters);

/// The chart point x(t) of a line (one-root in Q(sqrt(L))).
[[nodiscard]] ExactPoint line_point(const LineSupport& line, const FT& t);

/// Filtered exact comparison of |p - c|^2 with r2 for a one-root point p
/// (`px`, `py` its certified enclosures).
[[nodiscard]] CGAL::Comparison_result compare_to_circle(const ExactPoint& p, const Interval& px, const Interval& py,
                                                        const Point_2& c, const Interval& cx, const Interval& cy,
                                                        const FT& r2, const Interval& r2i, FilterMode mode,
                                                        DecisionCounters& counters);

/// Filtered exact side of a one-root point relative to a line.
[[nodiscard]] CGAL::Sign line_side(const LineSupport& line, const ExactPoint& p, const Interval& px, const Interval& py,
                                   FilterMode mode, DecisionCounters& counters);

/// Exact coincidence of two line supports: 0 when they are different lines,
/// +1 when they are the same line with codirected normals, -1 when the same line
/// with opposite normals.
[[nodiscard]] int line_coincidence(const LineSupport& a, const LineSupport& b);

/// Lines of the certificate, deduplicated by their construction key.
class LineTable {
public:
    /// Flank sigma = +1 (left of the edge, n.x = n.a + r sqrt(L)) or -1 of edge e.
    std::uint32_t flank(std::uint32_t edge, int sigma);
    /// The end line d.x = d.(a + tau d) of edge e at parameter tau.
    std::uint32_t end_line(std::uint32_t edge, const FT& tau);
    /// A rational line n.x = n.through (a target cut support, e.g. a declared
    /// boundary line); not deduplicated.
    std::uint32_t add_rational_line(const FT& nx, const FT& ny, const Point_2& through);

    LineTable(const CentreDomain& domain, const FT& tool_radius);

    [[nodiscard]] const LineSupport& operator[](std::uint32_t id) const { return lines_[id]; }
    [[nodiscard]] std::size_t size() const { return lines_.size(); }

    /// Global edge e: start vertex, direction, and its loop.
    [[nodiscard]] const Point_2& edge_start(std::uint32_t edge) const { return starts_[edge]; }
    [[nodiscard]] const Point_2& edge_end(std::uint32_t edge) const { return ends_[edge]; }
    [[nodiscard]] std::size_t edge_count() const { return starts_.size(); }
    [[nodiscard]] Point_2 edge_point(std::uint32_t edge, const FT& tau) const;

private:
    std::vector<LineSupport> lines_;
    std::vector<Point_2> starts_, ends_;
    FT tool_radius_;
    std::map<std::pair<std::uint32_t, int>, std::uint32_t> flanks_;
    std::map<std::pair<std::uint32_t, FT>, std::uint32_t> ends_by_tau_;
};

/// A closed rectangle region: four half-planes {side(line_i) == inward_i or 0}.
struct RectRegion {
    std::uint32_t edge;  ///< global edge index (loop-major)
    FT tau0, tau1;       ///< 0 <= tau0 < tau1 <= 1 on that edge
    std::array<std::uint32_t, 4> lines;  ///< flank +1, flank -1, end tau0, end tau1
    std::array<CGAL::Sign, 4> inward;    ///< side of each line's normal the rectangle lies on
    CGAL::Bbox_2 bbox;                   ///< certified outward box of the four corners
    bool target;                         ///< target edge capsule (else a merged connector piece)
    std::uint32_t report_index;          ///< connector index, or target member index

    /// Which of the four lines is the flank on side sigma.
    [[nodiscard]] static std::size_t flank_slot(int sigma) { return sigma > 0 ? 0 : 1; }
};

/// Exact side of a point relative to the closed rectangle.
enum class RectSide : std::uint8_t { Exterior = 0, Boundary = 1, Interior = 2 };
[[nodiscard]] RectSide rect_membership(const RectRegion& rect, const LineTable& lines, const ExactPoint& p);

/// Build the rectangle of edge e over [tau0, tau1] (exact, certified box).
[[nodiscard]] RectRegion make_rect(LineTable& lines, std::uint32_t edge, const FT& tau0, const FT& tau1, bool target,
                                   std::uint32_t report_index);

/// One connector: forward along loop `loop` from (segment0, tau0) to (segment1, tau1).
struct ConnectorSpec {
    std::uint32_t loop;
    std::uint32_t segment0;
    FT tau0;
    std::uint32_t segment1;
    FT tau1;
};

/// The connector capsules after the exact union of collinear pieces.
///
/// Every connector is a chain of straight pieces on domain edges. The union of
/// two capsules on one edge whose parameter intervals overlap or abut is the
/// capsule of the union interval, exactly (a point within r of [a, c] is within
/// r of [a, b] or [b, c]). So all connector pieces on one edge merge into the
/// maximal intervals of their union; a disk D_r(p) at an interior point of such
/// an interval lies inside the merged capsule. The swept set U is unchanged.
struct ConnectorCapsules {
    std::vector<RectRegion> rects;
    std::vector<Point_2> disk_centres;            ///< merged-interval endpoints, deduplicated
    std::vector<std::uint32_t> disk_report_index;  ///< lowest connector index with that endpoint
    std::size_t piece_count = 0;                   ///< straight pieces before the union

    /// Raises:
    ///     InvalidCertificateInputError: a site outside its loop, tau outside [0, 1),
    ///         or a connector with no forward progress.
    [[nodiscard]] static ConnectorCapsules build(const CentreDomain& domain, LineTable& lines,
                                                 const std::vector<ConnectorSpec>& connectors);
};

/// Global edge index of (loop, segment).
[[nodiscard]] std::uint32_t global_edge(const CentreDomain& domain, std::uint32_t loop, std::uint32_t segment);

}  // namespace compas_cgal::coverage_certificate
