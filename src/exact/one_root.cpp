#include "exact/one_root.h"

#include "exact/errors.h"

namespace compas_cgal::exact {

namespace {

void require_same_root(const OneRoot& a, const OneRoot& b)
{
    // a1() and root() are ONLY defined on an EXTENDED Sqrt_extension. A rational
    // coordinate carries a0() alone, and reading the extension parts
    // unconditionally is undefined behaviour -- the same guard as
    // stock_2.cpp::add_coordinate and engagement_2.cpp::as_radpoint. A
    // non-extended operand is compatible with any root, so it short-circuits.
    //
    // This is the same condition CGAL states in Sqrt_extension::check_roots,
    // promoted from a CGAL_precondition (erased under NDEBUG, which this
    // project builds with) to a raise that survives a release build.
    if (!a.is_extended() || !b.is_extended()) {
        return;
    }
    if (a.root() != b.root()) {
        throw CrossRootExtensionError(
            "one-root arithmetic requires operands in the same extension");
    }
}

}  // namespace

OneRoot same_root_add(const OneRoot& a, const OneRoot& b)
{
    require_same_root(a, b);
    return a + b;
}

OneRoot same_root_multiply(const OneRoot& a, const OneRoot& b)
{
    require_same_root(a, b);
    return a * b;
}

}  // namespace compas_cgal::exact
