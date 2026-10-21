#pragma once

// The declared (Pythagorean) boundary primitives and their exact distance
// closed forms (benchmarks/held_declared_clearance.py, _Line / _Arc).
//
// A line is start + lambda*u, lambda in [0, length], with an exactly unit
// rational direction u and left (interior) normal n = u^perp. An arc runs from
// start to end about its rational centre with rational radius, in sense +1
// (CCW) or -1 (CW). Every number is an Epeck FT; every decision is a kernel
// predicate or an exact Sqrt_extension comparison.

#include "held_declared_clearance_2/distance_2.h"
#include "held_declared_clearance_2/one_root_2.h"

#include <cstddef>
#include <vector>

namespace compas_cgal::held_declared_clearance {

enum class PrimitiveKind : int { Line = 0, Arc = 1 };

struct Primitive {
    std::size_t index = 0;
    PrimitiveKind kind = PrimitiveKind::Line;
    RationalPoint start;
    RationalPoint end;
    // Line only.
    Vector unit;    // exactly unit direction
    Vector normal;  // left normal u^perp (counter-clockwise quarter turn)
    FT length;
    // Arc only.
    RationalPoint centre;
    FT radius;
    int sense = 0;

    /// Validate a declared line: |u| == 1, length > 0, end == start + length * u exactly.
    ///
    /// Raises:
    ///     InvalidDeclaredPrimitiveError: any breach.
    [[nodiscard]] static Primitive line(std::size_t index, const RationalPoint& start, const RationalPoint& end,
                                        const Vector& unit, const FT& length);
    /// Validate a declared arc: radius > 0, sense +-1, both ends exactly on the circle.
    ///
    /// Raises:
    ///     InvalidDeclaredPrimitiveError: any breach.
    [[nodiscard]] static Primitive arc(std::size_t index, const RationalPoint& start, const RationalPoint& end,
                                       const RationalPoint& centre, const FT& radius, int sense);

    [[nodiscard]] bool is_arc() const { return kind == PrimitiveKind::Arc; }
};

/// Exact distance from a one-root point to the primitive.
[[nodiscard]] Distance distance_to(const Primitive& primitive, const OneRootPoint& point);

/// Whether the ray from an arc's centre through `probe` meets the arc (exact).
[[nodiscard]] bool spans(const Primitive& arc, const OneRootPoint& probe);
/// spans() for the rational direction from the arc's centre.
[[nodiscard]] bool spans_direction(const Primitive& arc, const Vector& direction);

/// The closed form of the nearest distance from a rational centre (one candidate).
[[nodiscard]] Distance nearest(const Primitive& primitive, const RationalPoint& centre);
/// The closed forms whose maximum is the farthest distance from a rational centre.
[[nodiscard]] std::vector<Distance> farthest(const Primitive& primitive, const RationalPoint& centre);

/// |point - origin|^2 in the point's extension.
[[nodiscard]] Root squared_distance(const OneRootPoint& point, const RationalPoint& origin);
/// Python's min(a, b) on one-root numbers: `first` unless second < first (exact, across extensions).
[[nodiscard]] const Root& smaller(const Root& first, const Root& second);

}  // namespace compas_cgal::held_declared_clearance
