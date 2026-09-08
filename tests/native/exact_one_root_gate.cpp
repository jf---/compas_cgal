#include "exact/one_root.h"

#include "exact/errors.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/number_utils.h>

#include <cstdio>
#include <exception>
#include <stdexcept>
#include <type_traits>

namespace exact = compas_cgal::exact;

namespace {

/// Raised when a gate check fails. Named so the gate obeys the same
/// named-exceptions-only rule as the module it exercises.
class GateCheckFailedError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// The project builds with CMAKE_BUILD_TYPE Release, which defines NDEBUG and
// erases <cassert>. A gate whose checks vanish reports a vacuous pass, so every
// decision here goes through a helper that is never compiled out.
void require(bool condition, const char* message)
{
    if (!condition) {
        throw GateCheckFailedError(message);
    }
}

// one_root.h claims OneRoot IS the arrangement package's coordinate type, which
// is the whole reason it is an alias rather than a wrapper. The claim is a type
// identity, so it is checked as one: a mismatch in either Sqrt_extension tag
// would silently reintroduce the conversion boundary the alias exists to avoid,
// and no runtime check could see it.
using GpsTraits = CGAL::Gps_circle_segment_traits_2<
    CGAL::Exact_predicates_exact_constructions_kernel>;
static_assert(
    std::is_same_v<exact::OneRoot, GpsTraits::Point_2::CoordNT>,
    "OneRoot is not Gps_circle_segment_traits_2<Epeck>::Point_2::CoordNT");

void same_root_operations_are_exact()
{
    // (1 + sqrt(2)) + (3 + 2*sqrt(2)) == 4 + 3*sqrt(2)
    const exact::OneRoot a(exact::Rational(1), exact::Rational(1), exact::Rational(2));
    const exact::OneRoot b(exact::Rational(3), exact::Rational(2), exact::Rational(2));
    const exact::OneRoot sum = exact::same_root_add(a, b);
    require(
        sum == exact::OneRoot(exact::Rational(4), exact::Rational(3), exact::Rational(2)),
        "same_root_add is not (1+sqrt2)+(3+2sqrt2) == 4+3sqrt2");
    require(CGAL::sign(sum) == CGAL::POSITIVE, "4+3sqrt2 is not positive");

    // (1 + sqrt(2)) * (1 - sqrt(2)) == -1, so the product of two extended values
    // collapsing back to a rational is exact too.
    const exact::OneRoot conjugate(
        exact::Rational(1), exact::Rational(-1), exact::Rational(2));
    const exact::OneRoot product = exact::same_root_multiply(a, conjugate);
    require(
        product == exact::OneRoot(exact::Rational(-1)),
        "same_root_multiply is not (1+sqrt2)*(1-sqrt2) == -1");
    require(CGAL::sign(product) == CGAL::NEGATIVE, "-1 is not negative");
}

void a_rational_operand_is_always_permitted()
{
    // A non-extended value carries a0() alone and is compatible with any root.
    const exact::OneRoot rational(exact::Rational(5));
    require(!rational.is_extended(), "a rational operand reports itself extended");
    const exact::OneRoot extended(
        exact::Rational(0), exact::Rational(1), exact::Rational(3));
    const exact::OneRoot sum = exact::same_root_add(rational, extended);
    require(CGAL::sign(sum) == CGAL::POSITIVE, "5+sqrt3 is not positive");

    // The same in the other argument position: neither operand may be assumed
    // extended, so both orders must reach the short circuit.
    const exact::OneRoot swapped = exact::same_root_add(extended, rational);
    require(sum == swapped, "same_root_add is not commutative across a rational");
    require(
        CGAL::sign(exact::same_root_multiply(rational, extended)) == CGAL::POSITIVE,
        "5*sqrt3 is not positive");
}

void cross_root_arithmetic_is_rejected()
{
    const exact::OneRoot two(
        exact::Rational(0), exact::Rational(1), exact::Rational(2));
    const exact::OneRoot three(
        exact::Rational(0), exact::Rational(1), exact::Rational(3));
    bool raised = false;
    try {
        static_cast<void>(exact::same_root_add(two, three));
    } catch (const exact::CrossRootExtensionError&) {
        raised = true;
    }
    require(raised, "cross-root addition was accepted");

    raised = false;
    try {
        static_cast<void>(exact::same_root_multiply(two, three));
    } catch (const exact::CrossRootExtensionError&) {
        raised = true;
    }
    require(raised, "cross-root multiplication was accepted");
}

void an_extended_zero_root_is_not_a_rational_operand()
{
    // ACDE_TAG == Tag_true -- the arrangement's tag, and therefore ours -- admits
    // an EXTENDED value whose root is zero. Such a value denotes a rational, yet
    // it is not a rational OPERAND: the exemption is `!is_extended()`, never
    // "the root is zero". Roots 2 and 3 alone cannot tell those two rules apart,
    // so the distinction is pinned here.
    const exact::OneRoot extended_zero(
        exact::Rational(0), exact::Rational(1), exact::Rational(0));
    require(
        extended_zero.is_extended(),
        "an extended value with root 0 reports itself unextended");

    const exact::OneRoot three(
        exact::Rational(0), exact::Rational(1), exact::Rational(3));
    bool raised = false;
    try {
        static_cast<void>(exact::same_root_add(extended_zero, three));
    } catch (const exact::CrossRootExtensionError&) {
        raised = true;
    }
    require(raised, "an extended root-0 operand was added across sqrt(3)");

    raised = false;
    try {
        static_cast<void>(exact::same_root_multiply(three, extended_zero));
    } catch (const exact::CrossRootExtensionError&) {
        raised = true;
    }
    require(raised, "an extended root-0 operand was multiplied across sqrt(3)");

    // The contrast that gives the check its meaning: the genuinely non-extended
    // zero denotes the same number and IS permitted against the same operand.
    const exact::OneRoot unextended_zero(exact::Rational(0));
    require(
        !unextended_zero.is_extended(),
        "an unextended zero reports itself extended");
    require(
        exact::same_root_add(unextended_zero, three) == three,
        "an unextended operand was not permitted against sqrt(3)");
}

}  // namespace

int main()
{
    try {
        same_root_operations_are_exact();
        a_rational_operand_is_always_permitted();
        cross_root_arithmetic_is_rejected();
        an_extended_zero_root_is_not_a_rational_operand();
    } catch (const std::exception& error) {
        std::printf("exact_one_root_gate FAILED: %s\n", error.what());
        return 1;
    }
    std::printf("exact_one_root_gate OK\n");
    return 0;
}
