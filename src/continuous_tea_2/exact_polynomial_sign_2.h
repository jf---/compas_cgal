#pragma once

#include "../exact_algebraic_1.h"

#include <CGAL/enum.h>

CGAL::Sign exact_polynomial_sign_at_2(
    const ExactAlgebraicKernel2::Polynomial_2& polynomial,
    const ExactAlgebraicKernel2::Algebraic_real_2& point);
