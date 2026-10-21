#pragma once

// P1 for the one-root vertex circles |x - w| = r of the exact target
// (exact_target_2.h): w = (w0x + w1x s, w0y + w1y s), s = sqrt(delta), delta not a
// rational square (folded, exactness checklist q17). They are target cut
// supports: they cut swept circles and swept flanks, and flip nothing.
//
// On a rational circle C (centre c, radius rho): with e = w - c, D = |e|^2,
// q = D + rho^2 - r^2, disc = 4 D rho^2 - q^2 (all in Q(s)), the crossings are
//     p = c + (q / 2D) e +/- (sqrt(disc) / 2D) perp(e),   perp(e) = (-e_y, e_x),
// towers m + n sqrt(disc) over Q(s) (tower_2.h), branch 0 the + sign: the same
// shape as the circle-circle and circle-line events, over the radicand delta.
//
// On a line support (chart x(t) = b0 + b1 sqrt(L) + t d, lines_2.h):
//     |x(t) - w|^2 = r^2  <=>  Ld t^2 - 2 mu t + gamma = 0,
//     mu = -d.(b0 - w) in Q(s),  gamma, disc = mu^2 - Ld gamma in K = Q(sqrt(L))(s),
//     t = (mu +/- sqrt(disc)) / Ld,
// a depth-3 number. It is never formed as a tower: the sweep needs only
// sign(t - q) for rational q, which is sign(M +/- sqrt(disc)) with M = mu - q Ld in
// Q(s), decided by the mixed-radical rule with one tower sign over Q(sqrt(L))
// for compare(M^2, disc) (elements of K are m + n s with m, n in Q(sqrt(L))).
//
// Two-support coincidences are decided by provenance: a tower point X lies on
// a vertex circle iff |X - w|^2 - r^2 = P + Q s = 0 with towers P, Q
// (sign_tower_root), and on its branch by one more such sign.

#include "coverage_certificate_2/counters_2.h"
#include "coverage_certificate_2/kernel_2.h"
#include "coverage_certificate_2/lines_2.h"
#include "coverage_certificate_2/regions_2.h"
#include "coverage_certificate_2/separation_2.h"
#include "coverage_certificate_2/tower_2.h"

#include <CGAL/Bbox_2.h>

#include <array>
#include <cstdint>
#include <memory>

namespace compas_cgal::coverage_certificate {

/// A vertex circle with an irrational one-root centre.
struct VertexCircle {
    OneRoot x, y;  ///< extended over sqrt(delta) (at least one of them)
    FT delta;
    FT r_sq;
    Interval x_i, y_i, r_i, r_sq_i;
    CGAL::Bbox_2 bbox;
    std::uint32_t disk;  ///< index of the target vertex disk
};

/// P + Q sqrt(delta) for towers P, Q sharing one field and radicand.
struct TowerPlusRoot {
    Tower P, Q;
    FT delta;

    [[nodiscard]] static TowerPlusRoot of(const Tower& t, const FT& delta);
    /// A one-root v over sqrt(delta) (or rational).
    [[nodiscard]] static TowerPlusRoot of(const OneRoot& v, const FT& delta);
    [[nodiscard]] CGAL::Sign sign() const { return sign_tower_root(P, Q, delta); }
};
[[nodiscard]] TowerPlusRoot operator+(const TowerPlusRoot& a, const TowerPlusRoot& b);
[[nodiscard]] TowerPlusRoot operator-(const TowerPlusRoot& a, const TowerPlusRoot& b);
[[nodiscard]] TowerPlusRoot operator*(const TowerPlusRoot& a, const TowerPlusRoot& b);

/// sign(|X - w|^2 - r^2) for the tower point (x, y): zero iff X is on V.
[[nodiscard]] CGAL::Sign vertex_circle_side(const Tower& x, const Tower& y, const VertexCircle& V);

// ---------------------------------------------------------------------------
// On a rational circle
// ---------------------------------------------------------------------------

struct CircleVertexExact {
    std::uint8_t count = 0;   ///< 0, 1 (tangent) or 2
    std::array<Tower, 2> x, y;
    OneRoot mx, my, ex, ey;  ///< chord midpoint m and e = w - c (branch provenance)
};

struct CircleVertexRecord {
    std::uint32_t circle;
    std::uint32_t vertex;
    std::uint8_t count;
    std::array<Interval, 2> x, y;
    std::unique_ptr<CircleVertexExact> exact;
};

[[nodiscard]] CircleVertexExact exact_circle_vertex(const BoundaryCircle& circle, const VertexCircle& V);

/// Filtered record: the interval front decides the count when it can.
[[nodiscard]] CircleVertexRecord make_circle_vertex_record(std::uint32_t circle_id, const BoundaryCircle& circle,
                                                           std::uint32_t vertex_id, const VertexCircle& V, FilterMode mode,
                                                           DecisionCounters& counters);

/// The exact crossings (materialised on demand, checked against the enclosures).
const CircleVertexExact& exact_of(CircleVertexRecord& record, const BoundaryCircle& circle, const VertexCircle& V,
                                  DecisionCounters& counters);

/// Whether the tower point X (on circle C) is the crossing `branch` of C and V.
[[nodiscard]] bool is_circle_vertex_crossing(const Tower& x, const Tower& y, const CircleVertexExact& exact,
                                             std::uint8_t branch, const VertexCircle& V);

// ---------------------------------------------------------------------------
// On a line support
// ---------------------------------------------------------------------------

struct LineVertexExact {
    std::uint8_t count = 0;  ///< 0, 1 (tangent) or 2
    OneRoot mu;              ///< over sqrt(delta)
    OneRoot disc_a, disc_b;  ///< disc = disc_a + disc_b sqrt(delta), parts over sqrt(L)
    FT delta;
    FT ld;  ///< |d|^2 > 0

    /// sign(t_branch - q): branch 0 is the + root.
    [[nodiscard]] CGAL::Sign sign_minus(std::uint8_t branch, const FT& q) const;
    /// sign(t_branch - t) for a tower t over the line's field Q(sqrt(L)): a
    /// constant number of signs in K = Q(sqrt(L))(sqrt(delta)), no bisection.
    [[nodiscard]] CGAL::Sign sign_minus(std::uint8_t branch, const Tower& t) const;
    /// Separation measure of t_branch (separation_2.h) and its radical count.
    [[nodiscard]] SeparationMeasure measure() const;
    static constexpr int RADICALS = 3;  ///< sqrt(L), sqrt(delta), sqrt(disc)
};

/// An element a + b sqrt(delta) of K = Q(sqrt(L))(sqrt(delta)), a, b one-roots
/// over sqrt(L) (one L for all operands; checked by the one-root arithmetic).
struct KNumber {
    OneRoot a, b;
    FT delta;

    [[nodiscard]] CGAL::Sign sign() const { return sign_tower({a, b, OneRoot(delta)}); }
    [[nodiscard]] Interval interval() const;
};
[[nodiscard]] KNumber operator+(const KNumber& x, const KNumber& y);
[[nodiscard]] KNumber operator-(const KNumber& x, const KNumber& y);
[[nodiscard]] KNumber operator*(const KNumber& x, const KNumber& y);

/// Exact sign of A + B sqrt(D) for A, B, D in K, D >= 0 (the mixed-radical rule,
/// compare(A^2, B^2 D) as one K sign).
[[nodiscard]] CGAL::Sign sign_k_root(const KNumber& A, const KNumber& B, const KNumber& D);

struct LineVertexRecord {
    std::uint32_t vertex;
    std::uint8_t count;
    std::array<Interval, 2> t;
    std::unique_ptr<LineVertexExact> exact;
};

[[nodiscard]] LineVertexExact exact_line_vertex(const LineSupport& line, const VertexCircle& V);

[[nodiscard]] LineVertexRecord make_line_vertex_record(const LineSupport& line, std::uint32_t vertex_id,
                                                       const VertexCircle& V, FilterMode mode, DecisionCounters& counters);

const LineVertexExact& exact_of(LineVertexRecord& record, const LineSupport& line, const VertexCircle& V,
                                DecisionCounters& counters);

/// Whether the tower parameter t on `line` is the root `branch` of V on it
/// (the point x(t) lies on V and on that side of the chord midpoint).
[[nodiscard]] bool is_line_vertex_root(const Tower& t, const LineSupport& line, const LineVertexExact& exact,
                                       std::uint8_t branch, const VertexCircle& V);

}  // namespace compas_cgal::coverage_certificate
