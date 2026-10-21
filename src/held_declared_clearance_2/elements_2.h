#pragma once

// Path elements: loci of the tool centre (benchmarks/held_declared_clearance.py,
// SegmentElement / ArcElement / CircleElement).
//
// A segment runs on the rational line normal . x = level between one-root
// endpoints (possibly of different extensions), along the rational direction.
// An arc runs on the rational circle (centre, radius) from start to end in
// sense; start == end is the full circle. A machining circle has a rational
// centre and radius (binary64 values injected exactly).

#include "held_declared_clearance_2/distance_2.h"
#include "held_declared_clearance_2/one_root_2.h"

namespace compas_cgal::held_declared_clearance {

struct SegmentElement {
    OneRootPoint start;
    OneRootPoint end;
    Vector direction;  // along the travel direction
    Vector normal;     // normal . direction == 0
    FT level;          // normal . x == level on the line

    /// Raises:
    ///     DegenerateClearanceElementError: zero direction or normal, or a normal not orthogonal to it.
    [[nodiscard]] static SegmentElement build(const OneRootPoint& start, const OneRootPoint& end,
                                              const Vector& direction, const Vector& normal, const FT& level);
    /// The segment between two distinct rational points (a binary64 polyline piece).
    ///
    /// Raises:
    ///     DegenerateClearanceElementError: the points coincide.
    [[nodiscard]] static SegmentElement between(const RationalPoint& start, const RationalPoint& end);

    /// d . x of a one-root point (exact, in its extension).
    [[nodiscard]] Root position(const OneRootPoint& point) const;
    /// Whether the foot of `point` on the line lies in the closed segment (exact, across extensions).
    [[nodiscard]] bool contains_parameter_of(const RationalPoint& point) const;
    /// Exact distance from a rational point.
    [[nodiscard]] Distance distance_to_point(const RationalPoint& point) const;
};

struct ArcElement {
    RationalPoint centre;
    FT radius;
    int sense = 0;
    OneRootPoint start;
    OneRootPoint end;
    bool full = false;  // start == end exactly

    /// Raises:
    ///     DegenerateClearanceElementError: radius <= 0 or sense not +-1.
    [[nodiscard]] static ArcElement build(const RationalPoint& centre, const FT& radius, int sense,
                                          const OneRootPoint& start, const OneRootPoint& end);

    /// Whether the ray from the centre through `probe` meets the arc (exact).
    [[nodiscard]] bool spans(const OneRootPoint& probe) const;
    /// Exact distance from a rational point.
    [[nodiscard]] Distance distance_to_point(const RationalPoint& point) const;
};

struct CircleElement {
    RationalPoint centre;
    FT radius;
};

}  // namespace compas_cgal::held_declared_clearance
