#include "exact/rational.h"

#include "exact/errors.h"

#include <cstdio>
#include <exception>
#include <limits>
#include <stdexcept>

namespace exact = compas_cgal::exact;

namespace {

// The project builds with CMAKE_BUILD_TYPE Release, which defines NDEBUG and
// erases <cassert>. A gate whose checks vanish reports a vacuous pass, so every
// decision here goes through a helper that is never compiled out.
void require(bool condition, const char* message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void from_binary64_is_exact()
{
    // 0.1 is not a decimal fraction; its exact value is 3602879701896397 / 2^55.
    const exact::Rational tenth = exact::from_binary64(0.1);
    const exact::Rational expected =
        exact::Rational(3602879701896397LL) / exact::Rational(1LL << 55);
    require(tenth == expected, "from_binary64(0.1) is not the dyadic value");

    // Smallest positive subnormal must survive exactly.
    const exact::Rational tiny = exact::from_binary64(5e-324);
    require(tiny > exact::Rational(0), "subnormal 5e-324 collapsed to zero");

    // Negative zero and positive zero denote the same rational.
    require(
        exact::from_binary64(-0.0) == exact::from_binary64(0.0),
        "signed zeros denote different rationals");
    require(
        exact::from_binary64(0.0) == exact::Rational(0),
        "from_binary64(0.0) is not zero");
}

void non_finite_input_is_rejected()
{
    bool raised = false;
    try {
        exact::from_binary64(std::numeric_limits<double>::quiet_NaN());
    } catch (const exact::NonFiniteBinary64Error&) {
        raised = true;
    }
    require(raised, "NaN was accepted by from_binary64");

    raised = false;
    try {
        exact::from_binary64(std::numeric_limits<double>::infinity());
    } catch (const exact::NonFiniteBinary64Error&) {
        raised = true;
    }
    require(raised, "infinity was accepted by from_binary64");
}

}  // namespace

int main()
{
    try {
        from_binary64_is_exact();
        non_finite_input_is_rejected();
    } catch (const std::exception& error) {
        std::printf("exact_vocabulary_gate FAILED: %s\n", error.what());
        return 1;
    }
    std::printf("exact_vocabulary_gate OK\n");
    return 0;
}
