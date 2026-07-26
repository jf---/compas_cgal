#pragma once

#include "../exact_algebraic_1.h"

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

#include <CGAL/CORE/BigRat.h>

namespace continuous_tea_2::event_partition_internal {

using Integer = CORE::BigInt;
using Rational = CORE::BigRat;
using Polynomial =
    ExactAlgebraicKernel1::Polynomial_1;

Rational parse_rational(
    const std::string& text,
    std::string_view role);

std::vector<Rational> parse_values(
    const std::vector<std::string>& values,
    std::size_t expected_size,
    std::string_view role);

std::string rational_text(
    const Rational& value);

Polynomial polynomial_from_integers(
    const std::vector<Integer>& coefficients);

std::vector<Integer> primitive_coefficients(
    const Polynomial& polynomial);

std::vector<std::string> coefficient_text(
    const std::vector<Integer>& coefficients);

} // namespace continuous_tea_2::event_partition_internal
