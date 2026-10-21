#pragma once

// Wire decoding of the clearance oracle's nanobind surface.
//
// Rationals cross as (numerator, denominator) decimal text pairs (the Python
// `Fraction` wire, as held_centre_domain_2_bindings.cpp), one-root numbers as
// (rational, coefficient, radicand) triples of those, points as pairs of
// one-root numbers. Binary64 coordinates cross as doubles and are injected
// exactly (every finite binary64 is a dyadic rational). No tolerance, no
// snapping at the seam.

#include "held_declared_clearance_2/elements_2.h"
#include "held_declared_clearance_2/one_root_2.h"
#include "held_declared_clearance_2/primitives_2.h"

#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace compas_cgal::held_declared_clearance {

using RationalWire = std::pair<std::string, std::string>;
using RootWire = std::tuple<RationalWire, RationalWire, RationalWire>;
using PointWire = std::pair<RootWire, RootWire>;
/// (kind, sense, rationals): a line (0, 0, [sx, sy, ex, ey, ux, uy, length]);
/// an arc (1, sense, [sx, sy, ex, ey, cx, cy, radius]).
using PrimitiveWire = std::tuple<int, int, std::vector<RationalWire>>;
/// (kind, sense, owners, start, end, rationals): a segment (0, 0, owners, s, e,
/// [dx, dy, nx, ny, level]); an arc (1, sense, owners, s, e, [cx, cy, radius]).
using ElementWire = std::tuple<int, int, std::vector<std::size_t>, PointWire, PointWire, std::vector<RationalWire>>;

inline constexpr int LINE_KIND = 0;
inline constexpr int ARC_KIND = 1;

/// Raises: InvalidClearanceWireError on malformed text or a non-positive denominator.
[[nodiscard]] FT rational_from_wire(const RationalWire& wire);
/// Raises: InvalidClearanceWireError, NegativeRadicandError.
[[nodiscard]] Root root_from_wire(const RootWire& wire);
/// Raises: InvalidClearanceWireError, MixedRootPointError.
[[nodiscard]] OneRootPoint point_from_wire(const PointWire& wire);
/// Raises: InvalidClearanceWireError, InvalidDeclaredPrimitiveError.
[[nodiscard]] Primitive primitive_from_wire(std::size_t index, const PrimitiveWire& wire);
/// Raises: InvalidClearanceWireError on a non-finite value.
[[nodiscard]] FT binary64_to_exact(double value);

/// A decoded element: exactly one of `segment` / `arc` is meaningful, by `is_arc`.
struct DecodedElement {
    bool is_arc = false;
    SegmentElement segment;
    ArcElement arc;
    std::vector<std::size_t> owners;
};
/// Raises: InvalidClearanceWireError (also an owner out of range), DegenerateClearanceElementError.
[[nodiscard]] DecodedElement element_from_wire(const ElementWire& wire, std::size_t primitive_count);

}  // namespace compas_cgal::held_declared_clearance
