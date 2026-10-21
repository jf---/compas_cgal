#pragma once

// Curve-to-primitive distance candidates (benchmarks/held_declared_clearance.py,
// _segment_candidates / _arc_candidates).
//
// The distance between two curves that do not cross is the minimum over the
// closed-form critical pairs: each endpoint against the other curve, and the
// interior-interior pairs with a common normal (the arc point on the normal of
// a line through its centre; the two arc points on the line of centres). Every
// candidate is the distance of a pair of real points of the two curves.
// Crossings (distance zero) are decided separately (crossings_2.h).

#include "held_declared_clearance_2/distance_2.h"
#include "held_declared_clearance_2/elements_2.h"
#include "held_declared_clearance_2/primitives_2.h"

#include <vector>

namespace compas_cgal::held_declared_clearance {

[[nodiscard]] std::vector<Distance> segment_candidates(const SegmentElement& element, const Primitive& primitive);
[[nodiscard]] std::vector<Distance> arc_candidates(const ArcElement& element, const Primitive& primitive);

}  // namespace compas_cgal::held_declared_clearance
