#pragma once

// The exact cap certificate: the single deciding predicate for "does this
// engaged run exceed the cap".
//
// It used to exist twice -- engagement_2.cpp (the guarded engagement lane) and
// audit_exact_station_2.cpp (the pure-exact audit lane) -- with identical
// coefficients, identical orientation branching and a line-for-line duplicated
// gap-closure walk, both linked into the same _stock_2 module. A comment at
// engagement_2.cpp conceded the audit lane decided it "VERBATIM". Per the
// Single-Path rule the two LANES stay isolated peers, but the DECIDING half may
// exist only once, and this is it.
//
// The two copies were not in fact equivalent, and the difference decided which
// one survives. engagement_2.cpp read the shared root inline and stated the
// shared-root precondition with CGAL_assertion, which NDEBUG compiles out of the
// Release build this project ships: on a point whose two coordinates carry
// different roots it silently took the y root and fed a wrong radicand to
// sign_mixed_radical -- a wrong cap verdict, with no diagnostic. The audit lane
// routed the same decomposition through decompose_exact_one_root_point, which
// throws ExactOneRootCoordinateMismatchError instead. The loud path is kept.
// exact/one_root.h:30-32 already named as_radpoint as the site that bypassed the
// shared helper.

#include "audit_certification_2.h"  // ExactOneRootPoint2, decompose_exact_one_root_point
#include "exact/one_root.h"         // compas_cgal::exact::sign_mixed_radical
#include "exact_motion_2.h"         // Epeck, EPoint
#include "stock_2.h"                // GpsPoint

#include <CGAL/enum.h>
#include <CGAL/number_utils.h>

namespace compas_cgal::cap {

/// One arrangement point split into rational parts against its single root.
///
/// A Gps_circle_segment_traits_2 point is a circle/circle intersection with one
/// shared discriminant, so BOTH coordinates of one point carry the same root; a
/// rational coordinate has root 0 and a1 == 0. Distinct points may carry
/// distinct roots (alpha for p, beta for q), which is why the predicate below
/// is a mixed two-radical sign and not a same-root compare.
struct RadialPoint2 {
    Epeck::FT x0;
    Epeck::FT x1;
    Epeck::FT y0;
    Epeck::FT y1;
    Epeck::FT root;
};

/// Split `point` into rational parts and its root.
///
/// Raises:
///     ExactOneRootCoordinateMismatchError: if both coordinates are extended
///         with different roots, via decompose_exact_one_root_point. That is a
///         malformed arrangement point, not a tolerance question.
inline RadialPoint2 radial_point(const GpsPoint& point)
{
    const GpsPoint::CoordNT& x = point.x();
    const GpsPoint::CoordNT& y = point.y();
    const ExactOneRootPoint2 decomposition = decompose_exact_one_root_point(point);
    return {
        x.a0(),
        x.is_extended() ? x.a1() : Epeck::FT(0),
        y.a0(),
        y.is_extended() ? y.a1() : Epeck::FT(0),
        decomposition.root(),
    };
}

/// Exact cap certificate for one maximal engaged run with CCW endpoints
/// `start`, `end` on the cutter circle centred at (cx, cy).
///
/// Args:
///     start, end: CCW endpoints of the run on the cutter circle.
///     cx, cy: cutter centre.
///     cap_chord_ratio: the exact rational surrogate for the angular cap,
///         4*sin^2(cap/2), injected by the caller (docs/exactness.md, boundary
///         doctrine for doubles).
///     squared_threshold: cap_chord_ratio * tool_radius^2, the exact squared
///         chord the run is compared against.
///
/// Returns:
///     true iff the run's angular extent exceeds the cap.
inline bool run_exceeds_cap(
    const GpsPoint& start,
    const GpsPoint& end,
    const Epeck::FT& cx,
    const Epeck::FT& cy,
    const Epeck::FT& cap_chord_ratio,
    const Epeck::FT& squared_threshold)
{
    // Closed loop: a full 2*pi run exceeds any cap <= pi.
    if (start == end) {
        return true;
    }

    const RadialPoint2 p = radial_point(start);
    const RadialPoint2 q = radial_point(end);

    // orientation(center, p, q) is the sign of the determinant of (p - center)
    // and (q - center). With p - center = (ax + bx*sqrt(alpha), ay + by*sqrt(alpha))
    // and q - center = (c + d*sqrt(beta), e + f*sqrt(beta)) the determinant is
    // A + B*sqrt(alpha) + C*sqrt(beta) + D*sqrt(alpha*beta).
    const Epeck::FT ax = p.x0 - cx;
    const Epeck::FT ay = p.y0 - cy;
    const Epeck::FT c = q.x0 - cx;
    const Epeck::FT e = q.y0 - cy;
    const CGAL::Sign orientation = compas_cgal::exact::sign_mixed_radical(
        ax * e - ay * c,
        p.x1 * e - p.y1 * c,
        ax * q.y1 - ay * q.x1,
        p.x1 * q.y1 - p.y1 * q.x1,
        p.root,
        q.root);

    // CCW angle theta from start to end: orientation > 0 <=> theta < pi,
    // orientation == 0 <=> theta == pi, orientation < 0 <=> theta > pi.
    if (orientation == CGAL::NEGATIVE) {
        return true;  // theta > pi >= cap
    }
    if (orientation == CGAL::ZERO) {
        return CGAL::compare(cap_chord_ratio, Epeck::FT(4)) == CGAL::SMALLER;
    }

    // theta < pi: the chord grows monotonically in theta, so the run exceeds the
    // cap iff |pq|^2 > squared_threshold. With p - q carrying both sqrt(alpha)
    // and sqrt(beta), |pq|^2 expands to a + b*sqrt(alpha) + c*sqrt(beta) +
    // d*sqrt(alpha*beta).
    const Epeck::FT dx = q.x0 - p.x0;
    const Epeck::FT dy = q.y0 - p.y0;
    const Epeck::FT a2 = dx * dx + dy * dy
        + (p.x1 * p.x1 + p.y1 * p.y1) * p.root
        + (q.x1 * q.x1 + q.y1 * q.y1) * q.root;
    const Epeck::FT b2 = -Epeck::FT(2) * (dx * p.x1 + dy * p.y1);
    const Epeck::FT c2 = Epeck::FT(2) * (dx * q.x1 + dy * q.y1);
    const Epeck::FT d2 = -Epeck::FT(2) * (p.x1 * q.x1 + p.y1 * q.y1);
    return compas_cgal::exact::sign_mixed_radical(
               a2 - squared_threshold, b2, c2, d2, p.root, q.root)
        == CGAL::POSITIVE;
}

/// Cutter centre as an exact point rather than a coordinate pair.
inline bool run_exceeds_cap(
    const GpsPoint& start,
    const GpsPoint& end,
    const EPoint& center,
    const Epeck::FT& cap_chord_ratio,
    const Epeck::FT& squared_threshold)
{
    return run_exceeds_cap(
        start, end, center.x(), center.y(), cap_chord_ratio, squared_threshold);
}

}  // namespace compas_cgal::cap
