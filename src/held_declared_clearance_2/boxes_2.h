#pragma once

// Outward binary64 boxes of primitives and elements, for pruning only: each
// end is an end of a certified enclosure of the exact bound (CGAL::to_interval
// of the lazy FT or of the Sqrt_extension), so the box contains the exact set.

#include "held_declared_clearance_2/elements_2.h"
#include "held_declared_clearance_2/primitive_index_2.h"
#include "held_declared_clearance_2/primitives_2.h"

namespace compas_cgal::held_declared_clearance {

/// A line's endpoint box; an arc's swept box (its ends, and each axis extreme
/// of its circle that the sweep spans, decided exactly).
[[nodiscard]] Box primitive_box(const Primitive& primitive);
/// The box of a disk (a machining circle).
[[nodiscard]] Box disk_box(const RationalPoint& centre, const FT& radius);
/// An arc element's swept box (its full circle's when it is full).
[[nodiscard]] Box arc_box(const ArcElement& element);
/// A segment element's endpoint box.
[[nodiscard]] Box segment_box(const SegmentElement& element);

}  // namespace compas_cgal::held_declared_clearance
