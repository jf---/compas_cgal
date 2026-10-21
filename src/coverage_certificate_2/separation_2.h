#pragma once

// Root-separation bounds for the exact search loops (review I1b).
//
// Every exact search loop of the certificate (P3 chart bisection, cross-base
// separation, the circle witness search) terminates because the two algebraic
// numbers it separates are distinct. Its contract bound is the theorem bound on
// how close two distinct numbers of the given bit sizes can be, not a fixed
// count.
//
// The bound (BFMSS-style, proved in docs/held_coverage_certificate.md,
// "Contract bounds"): write a value as xi = U / L, U and L algebraic integers
// built from integer leaves by + - * and square roots, with every conjugate of
// U bounded by u and of L by l (rules below; a conjugate flips the signs of
// radicals, which leaves these modulus bounds intact). U lies in a field of
// degree D <= 2^k over Q, k the number of distinct square roots, so a nonzero
// U has a nonzero integer norm and |U| >= u^(1 - D) (u >= 1); hence
//
//     xi != 0   =>   |xi| >= 2^-(  (D - 1) log2 u + log2 l  ).
//
// Rules (log2 domain, every step rounded up):
//   rational p / q:  u = |p|, l = q;
//   x +- y:          u = u_x l_y + l_x u_y, l = l_x l_y;
//   x * y:           u = u_x u_y, l = l_x l_y;
//   x / y:           u = u_x l_y, l = l_x u_y;
//   sqrt(x):         u = sqrt(u_x l_x), l = l_x   (sqrt(U / L) = sqrt(U L) / L).
//
// The measures are evaluated in doubles from exact integer bit lengths. They
// only bound a loop count (an analytic bound driving refinement, docs/exactness.md):
// the safe failure direction is a larger count, and every bound carries
// SEPARATION_SLACK_STEPS extra steps.

#include "coverage_certificate_2/kernel_2.h"
#include "coverage_certificate_2/tower_2.h"

#include <initializer_list>

namespace compas_cgal::coverage_certificate {

/// Integer slack added to every loop bound, steps. It absorbs the +1 of the
/// halving argument and the ceilings of the log2 sums; it is a fixed integer, so
/// it never makes a bound smaller than the theorem's.
inline constexpr int SEPARATION_SLACK_STEPS = 64;

/// Steps an exact loop runs before its separation bound is computed at all
/// (computing it needs exact parts; correct runs finish far below this).
inline constexpr int SEPARATION_LAZY_STEPS = 256;

/// log2 bounds of xi = U / L (see the header comment).
struct SeparationMeasure {
    double log2_u;
    double log2_l;

    [[nodiscard]] static SeparationMeasure of(const FT& value);
    [[nodiscard]] static SeparationMeasure of(const OneRoot& value);
    [[nodiscard]] static SeparationMeasure of(const Tower& value);
};

[[nodiscard]] SeparationMeasure plus(const SeparationMeasure& a, const SeparationMeasure& b);
[[nodiscard]] SeparationMeasure times(const SeparationMeasure& a, const SeparationMeasure& b);
[[nodiscard]] SeparationMeasure over(const SeparationMeasure& a, const SeparationMeasure& b);
[[nodiscard]] SeparationMeasure square_root(const SeparationMeasure& a);

/// Number of distinct nonzero radicands among `radicands` (an upper bound on k
/// for a value whose radicals are exactly these).
[[nodiscard]] int distinct_radicals(std::initializer_list<const FT*> radicands);

/// Distinct square roots of values in towers `towers`: their base roots and
/// their (possibly extended) tower radicands, each counted once.
[[nodiscard]] int tower_radicals(std::initializer_list<const Tower*> towers);

/// Bits b with |xi| >= 2^-b for any nonzero xi of measure `measure` over at
/// most `radicals` distinct square roots.
[[nodiscard]] double separation_bits(const SeparationMeasure& measure, int radicals);

/// The step bound of a halving loop that starts from a bracket of width
/// `width` and stops once the bracket is shorter than a nonzero gap whose
/// separation is `bits`: ceil(log2 width) + ceil(bits) + slack (and at least the
/// lazy threshold).
[[nodiscard]] long long halving_step_bound(double width, double bits);

/// The step budget of one exact search loop: SEPARATION_LAZY_STEPS free steps,
/// then the root-separation bound (computed once, on demand), then a loud
/// SeparationBoundExceededError naming the loop.
class SeparationBudget {
public:
    explicit SeparationBudget(const char* loop) : loop_(loop) {}

    /// Count one step; `bound` computes the theorem bound when first needed.
    template <class BoundFn>
    void step(BoundFn&& bound)
    {
        if (++steps_ <= limit_) return;
        if (!known_) {
            limit_ = bound();
            known_ = true;
            if (steps_ <= limit_) return;
        }
        fail();
    }

    [[nodiscard]] long long steps() const { return steps_; }

private:
    [[noreturn]] void fail() const;

    const char* loop_;
    long long steps_ = 0;
    long long limit_ = SEPARATION_LAZY_STEPS;
    bool known_ = false;
};

}  // namespace compas_cgal::coverage_certificate
