#pragma once

// Crossings (distance zero) of an element with a primitive, decided exactly
// rather than assumed away (benchmarks/held_declared_clearance.py,
// _segment_crosses / _arc_crosses). Intersection points of lines and circles
// with rational supports are one-root points in the discriminant's extension.

#include "held_declared_clearance_2/elements_2.h"
#include "held_declared_clearance_2/primitives_2.h"

namespace compas_cgal::held_declared_clearance {

[[nodiscard]] bool segment_crosses(const SegmentElement& element, const Primitive& primitive);
[[nodiscard]] bool arc_crosses(const ArcElement& element, const Primitive& primitive);

}  // namespace compas_cgal::held_declared_clearance
