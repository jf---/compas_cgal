// Reporting views follow Python's IEEE evaluation: no contraction of
// a + b * s into a fused multiply-add.
#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "held_declared_clearance_2/reporting_2.h"

#include "held_declared_clearance_2/errors_2.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>

namespace compas_cgal::held_declared_clearance {

namespace {

/// Nearest binary64 of a non-negative exact value enclosed by 0 <= low <= high,
/// ties to even. CGAL's enclosure is tight (adjacent doubles) for normal
/// values; near the subnormal range it may not be, so the enclosure is first
/// narrowed by bisection over the binary64 encodings (ordered like the values
/// for non-negative doubles), each step one exact comparison.
double nearest_non_negative(const FT::ET& exact, double low, double high)
{
    std::uint64_t lower = std::bit_cast<std::uint64_t>(low);
    std::uint64_t upper = std::bit_cast<std::uint64_t>(high);
    while (upper - lower > 1) {
        const std::uint64_t middle = lower + (upper - lower) / 2;
        if (FT::ET(std::bit_cast<double>(middle)) <= exact) {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    low = std::bit_cast<double>(lower);
    high = std::bit_cast<double>(upper);
    if (FT::ET(low) == exact) return low;
    if (FT::ET(high) == exact) return high;
    const FT::ET midpoint = (FT::ET(low) + FT::ET(high)) / 2;
    switch (CGAL::compare(exact, midpoint)) {
        case CGAL::SMALLER:
            return low;
        case CGAL::LARGER:
            return high;
        case CGAL::EQUAL:
            break;
    }
    return (lower & 1U) == 0 ? low : high;
}

}  // namespace

double binary64(const FT& value)
{
    const FT::ET& exact = CGAL::exact(value);
    const auto [low, high] = CGAL::to_interval(exact);
    if (std::isinf(low) || std::isinf(high)) throw ReportingOverflowError("reported value exceeds binary64 range");
    if (low == high) return low;
    switch (CGAL::sign(exact)) {
        case CGAL::ZERO:
            return 0.0;
        case CGAL::POSITIVE:
            return nearest_non_negative(exact, std::max(low, 0.0), high);
        case CGAL::NEGATIVE:
            break;
    }
    return -nearest_non_negative(-exact, std::max(-high, 0.0), -low);
}

double binary64(const Root& value)
{
    const double a = binary64(value.a0());
    const double b = binary64(value.a1());
    const double root = std::sqrt(binary64(value.root()));
    const double product = b * root;
    return a + product;
}

}  // namespace compas_cgal::held_declared_clearance
