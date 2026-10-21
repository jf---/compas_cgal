#pragma once

// A distance as one closed form over one-root numbers (benchmarks/
// held_declared_clearance.py, _Root / _Absolute / _ShiftedRoot): the exact
// sign of (distance - threshold) for a rational threshold, and the Python
// oracle's binary64 view of it (reporting only).

#include "held_declared_clearance_2/one_root_2.h"

#include <utility>

namespace compas_cgal::held_declared_clearance {

class Distance {
public:
    /// sqrt(square), square >= 0.
    [[nodiscard]] static Distance root(Root square);
    /// |value|.
    [[nodiscard]] static Distance absolute(Root value);
    /// |sqrt(square) + shift| (a point and a circle, or two circles on their centre line).
    [[nodiscard]] static Distance shifted_root(Root square, FT shift);

    /// Exact sign(distance - threshold).
    [[nodiscard]] int sign_minus(const FT& threshold) const;
    /// The Python oracle's binary64 view (reporting only, never a decision).
    [[nodiscard]] double approx() const;

private:
    enum class Kind { Root, Absolute, ShiftedRoot };
    Distance(Kind kind, Root value, FT shift) : kind_(kind), value_(std::move(value)), shift_(std::move(shift)) {}

    Kind kind_;
    Root value_;  // the square (Root, ShiftedRoot) or the value (Absolute)
    FT shift_;
};

}  // namespace compas_cgal::held_declared_clearance
