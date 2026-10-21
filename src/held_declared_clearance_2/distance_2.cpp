// The binary64 views must evaluate exactly as the Python oracle's.
#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "held_declared_clearance_2/distance_2.h"

#include "held_declared_clearance_2/reporting_2.h"

#include <algorithm>
#include <cmath>

namespace compas_cgal::held_declared_clearance {

namespace {

/// Exact sign(sqrt(square) - threshold) for square >= 0.
int sqrt_sign(const Root& square, const FT& threshold)
{
    if (CGAL::is_negative(threshold)) return 1;
    return sign(shifted(square, -(threshold * threshold)));
}

/// Python's sqrt(max(float(square), 0.0)).
double reported_root(const Root& square) { return std::sqrt(std::max(binary64(square), 0.0)); }

}  // namespace

Distance Distance::root(Root square) { return Distance(Kind::Root, std::move(square), FT(0)); }

Distance Distance::absolute(Root value) { return Distance(Kind::Absolute, std::move(value), FT(0)); }

Distance Distance::shifted_root(Root square, FT shift)
{
    return Distance(Kind::ShiftedRoot, std::move(square), std::move(shift));
}

int Distance::sign_minus(const FT& threshold) const
{
    switch (kind_) {
        case Kind::Root:
            return sqrt_sign(value_, threshold);
        case Kind::Absolute:
            // |v| - m = max(v - m, -v - m).
            return std::max(sign(shifted(value_, -threshold)), sign(shifted(scaled(FT(-1), value_), -threshold)));
        case Kind::ShiftedRoot:
            // |sqrt(D) + k| - m = max(sqrt(D) + k - m, -sqrt(D) - k - m).
            return std::max(sqrt_sign(value_, threshold - shift_), -sqrt_sign(value_, -shift_ - threshold));
    }
    return 1;  // unreachable: every kind returns above
}

double Distance::approx() const
{
    switch (kind_) {
        case Kind::Root:
            return reported_root(value_);
        case Kind::Absolute:
            return std::fabs(binary64(value_));
        case Kind::ShiftedRoot: {
            const double root = reported_root(value_);
            const double shift = binary64(shift_);
            return std::fabs(root + shift);
        }
    }
    return 0.0;  // unreachable: every kind returns above
}

}  // namespace compas_cgal::held_declared_clearance
