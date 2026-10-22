#pragma once

// Exact one-root geometry shared by the exact target (exact_target_2.h) and the
// exact dOmega* connectors (exact_connectors_2.h): signs of differences,
// cross products about a rational centre and squared distances to a one-root
// centre for points with coordinates in one field Q(sqrt(L)), and the wedge of
// an arc. Every sign is the tower/mixed-radical rule of tower_2.h; no
// cross-root arithmetic is formed.

#include "coverage_certificate_2/exact_target_2.h"
#include "coverage_certificate_2/kernel_2.h"
#include "coverage_certificate_2/lines_2.h"
#include "coverage_certificate_2/tower_2.h"

namespace compas_cgal::coverage_certificate {

/// sign(a - b) for one-roots over possibly different radicands.
[[nodiscard]] CGAL::Sign sign_difference_one_roots(const OneRoot& a, const OneRoot& b);

/// (w - c) x (p - c) for w over one radicand and p over another, as a tower over
/// p's field (tower_2.h, cross_product).
[[nodiscard]] Tower cross_about(const Point_2& c, const OneRoot& wx, const OneRoot& wy, const ExactPoint& p);

/// Interval shadow of cross_about.
[[nodiscard]] Interval cross_about_i(const Interval& cx, const Interval& cy, const Interval& wx, const Interval& wy,
                                     const Interval& px, const Interval& py);

/// |p - w|^2 - r2 for a one-root centre w, as a tower over p's field.
[[nodiscard]] Tower squared_distance_minus(const ExactPoint& p, const OneRoot& wx, const OneRoot& wy, const FT& r2);

/// Wedge membership of p for an arc about c from w0 to w1 in `sense` with span
/// `span` (closed; the apex c is inside). o0, o1 are the signs of
/// sense (w0 - c) x (p - c) and sense (p - c) x (w1 - c). A full turn
/// (w0 == w1) is a Major span and contains every point.
[[nodiscard]] bool in_wedge(ArcSpan span, CGAL::Sign o0, CGAL::Sign o1);

/// The open wedge: the same rule with strict signs (the bounding rays excluded).
[[nodiscard]] bool in_open_wedge(ArcSpan span, CGAL::Sign o0, CGAL::Sign o1);

/// s for sense +1, -s for sense -1.
[[nodiscard]] CGAL::Sign times(int sense, CGAL::Sign s);

/// The span of the arc about c from a to b in `sense`.
///
/// Raises:
///     InvalidCertificateInputError: a == b.
[[nodiscard]] ArcSpan span_of(const Point_2& c, const OneRoot& ax, const OneRoot& ay, const OneRoot& bx, const OneRoot& by,
                              int sense);

}  // namespace compas_cgal::coverage_certificate
