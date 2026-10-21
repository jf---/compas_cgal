#pragma once

// P1/P2 for circle-circle pieces: crossings of boundary circles behind a
// certified interval front, with CGAL's Circular_kernel_2 as the exact back end.
//
// For a candidate pair (A, B), D = |c_B - c_A|^2, q = D + rho_A^2 - rho_B^2:
//   disc = 4 D rho_A^2 - q^2   (Interval_nt, outward rounded)
//   * disc.sup() < 0: no crossing, certified.
//   * disc.inf() > 0: two transversal crossings, certified. They are carried as
//     enclosures m +/- s perp(d), m = c_A + (q / 2D) d, s = sqrt(disc) / 2D,
//     perp(d) = (-d_y, d_x). Branch 0 is the + branch.
//   * otherwise (or in EXACT_ONLY mode) CK::Intersect_2 runs at once; CGAL
//     reports tangency multiplicity exactly. Concentric pairs are decided by
//     exact centre and radius equality before the kernel is called (a
//     coincident pair is a shared support, P5, not a crossing).
//
// Exact points are built on demand. Their branch labels are decided by an
// exact rational matcher, never by box containment: the two crossings differ
// by 2 s perp(d) with s > 0, so when d_y != 0 the + branch is the point with the
// smaller x iff d_y > 0, and when d_y == 0 it is the point with the larger y iff
// d_x > 0. compare_x/compare_y on the two exact points and the rational signs
// of d_x, d_y decide it. The boxes are then only checked to CONTAIN their exact
// point (a contract check of the filter); overlapping boxes, the near-tangent
// case the spike could only abort on, are decided and counted
// (`exact_branch_matches`).

#include "coverage_certificate_2/counters_2.h"
#include "coverage_certificate_2/kernel_2.h"
#include "coverage_certificate_2/regions_2.h"
#include "coverage_certificate_2/tower_2.h"

#include <array>
#include <cstdint>
#include <memory>
#include <optional>
#include <vector>

namespace compas_cgal::coverage_certificate {

/// One crossing point of a boundary circle with another (an arrangement
/// vertex on the circle), with its certified enclosure.
struct CrossingEvent {
    std::uint32_t pair;         ///< circle-circle pair, or circle-line record (source 1)
    std::uint8_t branch;        ///< 0 = + branch, 1 = - branch; 0 for a tangency
    std::uint8_t multiplicity;  ///< 1 transversal (flips membership), 2 tangent (flips nothing)
    std::uint32_t other;        ///< index of the crossing boundary circle, or line id (source 1)
    Interval x, y;
    std::uint8_t source = 0;    ///< 0: circle-circle (Circular_kernel_2); 1: circle-line (tower)
};

class CrossingTable {
public:
    /// Candidate pairs by bbox (CGAL::box_self_intersection_d over certified
    /// boxes, a superset filter), then the filtered P1 on every pair.
    [[nodiscard]] static CrossingTable build(const CircleRegionSet& set, FilterMode mode, DecisionCounters& counters);

    /// Events on circle k, unordered.
    std::vector<std::vector<CrossingEvent>> events;
    /// Swept-region indices whose boundary-circle boxes overlap circle k's box,
    /// excluding k's own region; sorted, unique. Target circles are not regions:
    /// pairs of two target circles are never formed.
    std::vector<std::vector<std::uint32_t>> candidates;
    /// Boundary circles identical to circle k (same centre, same radius).
    std::vector<std::vector<std::uint32_t>> coincident;
    /// Target circles whose boxes overlap circle k's box (cut supports of a
    /// swept circle; empty for a target circle). Sorted.
    std::vector<std::vector<std::uint32_t>> target_circles;

    /// The exact crossing point of an event (materialised on demand).
    [[nodiscard]] const CircularPoint_2& exact_point(const CrossingEvent& event, DecisionCounters& counters);

    /// The event's exact x and y as rational-base towers m + n sqrt(disc):
    /// p = c_A + (q / 2D) d +/- (sqrt(disc) / 2D) perp(d), the + branch first
    /// (the same labels as the CK points). Built from the two circles with a few
    /// rational operations, without a Circular_kernel_2 construction; this is
    /// the exact back end of every comparison (P2).
    [[nodiscard]] const Tower& exact_x(const CrossingEvent& event, DecisionCounters& counters);
    [[nodiscard]] const Tower& exact_y(const CrossingEvent& event, DecisionCounters& counters);

    /// 0 when the event lies on the lower half of circle C (y < c_y, or y == c_y
    /// and x < c_x: the leftmost point opens the lower half), 1 otherwise.
    [[nodiscard]] int half(const CrossingEvent& event, const BoundaryCircle& circle, DecisionCounters& counters);

    /// Exact comparison of the x coordinates of two events on one circle.
    [[nodiscard]] CGAL::Comparison_result compare_x(const CrossingEvent& a, const CrossingEvent& b,
                                                    DecisionCounters& counters);

    /// Exact comparison of an event's x with the x of a rational point.
    [[nodiscard]] CGAL::Comparison_result compare_x(const CrossingEvent& a, const Point_2& q, const Interval& qx,
                                                    DecisionCounters& counters);

    [[nodiscard]] FilterMode mode() const { return mode_; }

    /// Also check, for every pair tower built, that the certified enclosure
    /// contains the exact value (a contract on the filter; on in audit runs).
    void verify_enclosures(const bool on) { verify_enclosures_ = on; }

private:
    struct Pair {
        std::uint32_t a, b;
        bool lazy;  ///< certified transversal by the filter, exact points not yet built
        std::array<Interval, 2> box_x, box_y;
        std::optional<std::array<CircularPoint_2, 2>> exact;
    };

    void materialise(Pair& pair, DecisionCounters& counters);

    struct PairTowers {
        std::array<Tower, 2> x, y;
    };
    const PairTowers& towers_of(std::uint32_t pair, DecisionCounters& counters);
    std::vector<std::unique_ptr<PairTowers>> towers_;  ///< per pair, on demand

    const CircleRegionSet* set_ = nullptr;
    FilterMode mode_ = FilterMode::Lazy;
    bool verify_enclosures_ = false;
    std::vector<Pair> pairs_;
};

}  // namespace compas_cgal::coverage_certificate
