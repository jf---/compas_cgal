#include "exact/rational.h"

#include "exact/errors.h"

#include <cmath>

namespace compas_cgal::exact {

Rational from_binary64(const double value)
{
    if (!std::isfinite(value)) {
        throw NonFiniteBinary64Error(
            "exact::from_binary64 requires a finite binary64");
    }
    // Epeck::FT's double constructor is exact: it stores the dyadic rational the
    // binary64 denotes. Constructing through it keeps the lazy interval intact,
    // where a bit-decomposition would build an eager cpp_rational instead.
    return Rational(value);
}

}  // namespace compas_cgal::exact
