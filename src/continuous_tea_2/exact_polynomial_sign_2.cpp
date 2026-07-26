#include "exact_polynomial_sign_2.h"

#include <iterator>
#include <utility>
#include <vector>

#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/number_utils.h>

CGAL::Sign exact_polynomial_sign_at_2(
    const ExactAlgebraicKernel2::Polynomial_2& polynomial,
    const ExactAlgebraicKernel2::Algebraic_real_2& point)
{
    using Kernel = ExactAlgebraicKernel2;
    using Polynomial = Kernel::Polynomial_2;
    using Traits = CGAL::Polynomial_traits_d<Polynomial>;

    if (CGAL::is_zero(polynomial)) {
        return CGAL::ZERO;
    }

    CGAL::Sign result = CGAL::sign(
        typename Traits::Innermost_leading_coefficient()(
            polynomial));
    std::vector<
        std::pair<
            Polynomial,
            Kernel::Multiplicity_type>>
        factors;
    Kernel kernel;
    kernel.square_free_factorize_2_object()(
        polynomial,
        std::back_inserter(factors));
    for (const auto& [factor, multiplicity] : factors) {
        const CGAL::Sign factor_sign =
            kernel.sign_at_2_object()(
                factor,
                point);
        if (factor_sign == CGAL::ZERO) {
            return CGAL::ZERO;
        }
        if (multiplicity % 2 != 0) {
            result = result * factor_sign;
        }
    }
    return result;
}
