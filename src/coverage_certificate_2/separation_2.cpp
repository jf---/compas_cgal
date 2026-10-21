#include "coverage_certificate_2/separation_2.h"

#include "coverage_certificate_2/errors_2.h"

#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace compas_cgal::coverage_certificate {

namespace {
using Traits = CGAL::Fraction_traits<FT>;
using Integer = Traits::Numerator_type;

/// log2 |n| <= bit length of n (0 for n == 0, an over-estimate there).
double bit_length(const Integer& n)
{
    if (n == 0) return 0.0;
    return static_cast<double>(boost::multiprecision::msb(boost::multiprecision::abs(n)) + 1);
}

/// log2(2^a + 2^b) <= max(a, b) + 1.
double log2_sum(const double a, const double b) { return std::max(a, b) + 1.0; }
}  // namespace

SeparationMeasure SeparationMeasure::of(const FT& value)
{
    Integer num;
    Integer den;
    Traits::Decompose()(value, num, den);
    return {bit_length(num), bit_length(den)};
}

SeparationMeasure SeparationMeasure::of(const OneRoot& value)
{
    if (!value.is_extended()) return of(value.a0());
    return plus(of(value.a0()), times(of(value.a1()), square_root(of(value.root()))));
}

SeparationMeasure SeparationMeasure::of(const Tower& value)
{
    return plus(of(value.m), times(of(value.n), square_root(of(value.d))));
}

SeparationMeasure plus(const SeparationMeasure& a, const SeparationMeasure& b)
{
    return {log2_sum(a.log2_u + b.log2_l, a.log2_l + b.log2_u), a.log2_l + b.log2_l};
}

SeparationMeasure times(const SeparationMeasure& a, const SeparationMeasure& b)
{
    return {a.log2_u + b.log2_u, a.log2_l + b.log2_l};
}

SeparationMeasure over(const SeparationMeasure& a, const SeparationMeasure& b)
{
    return {a.log2_u + b.log2_l, a.log2_l + b.log2_u};
}

SeparationMeasure square_root(const SeparationMeasure& a) { return {0.5 * (a.log2_u + a.log2_l), a.log2_l}; }

int distinct_radicals(const std::initializer_list<const FT*> radicands)
{
    std::vector<const FT*> seen;
    for (const FT* r : radicands) {
        if (CGAL::sign(*r) == CGAL::ZERO) continue;
        if (std::none_of(seen.begin(), seen.end(), [&](const FT* s) { return *s == *r; })) seen.push_back(r);
    }
    return static_cast<int>(seen.size());
}

int tower_radicals(const std::initializer_list<const Tower*> towers)
{
    std::vector<FT> bases;
    std::vector<const OneRoot*> radicands;
    for (const Tower* t : towers) {
        const FT b = base_of(*t);
        if (CGAL::sign(b) != CGAL::ZERO && std::find(bases.begin(), bases.end(), b) == bases.end()) bases.push_back(b);
        if (sign_one_root(t->d) == CGAL::ZERO) continue;
        if (std::none_of(radicands.begin(), radicands.end(), [&](const OneRoot* r) {
                return compatible_bases(Tower::of(*r), Tower::of(t->d)) &&
                       sign_difference(Tower::of(*r), Tower::of(t->d)) == CGAL::ZERO;
            })) {
            radicands.push_back(&t->d);
        }
    }
    return static_cast<int>(bases.size() + radicands.size());
}

double separation_bits(const SeparationMeasure& measure, const int radicals)
{
    const double degree = std::ldexp(1.0, radicals);
    return (degree - 1.0) * std::max(measure.log2_u, 0.0) + std::max(measure.log2_l, 0.0);
}

long long halving_step_bound(const double width, const double bits)
{
    const double log_width = width > 1.0 ? std::ceil(std::log2(width)) : 0.0;
    const double steps = log_width + std::ceil(bits) + SEPARATION_SLACK_STEPS;
    // Bounds beyond 2^62 steps are not reachable in practice; clamp to keep the
    // integer finite (a larger bound only delays the loud failure).
    const double cap = std::ldexp(1.0, 62);
    return static_cast<long long>(std::max(static_cast<double>(SEPARATION_LAZY_STEPS), std::min(steps, cap)));
}

void SeparationBudget::fail() const
{
    throw SeparationBoundExceededError(std::string(loop_) + " exceeded its root-separation bound (" +
                                       std::to_string(limit_) +
                                       " steps): the two values it separates are equal, so an upstream strict-order "
                                       "decision is wrong");
}

}  // namespace compas_cgal::coverage_certificate
