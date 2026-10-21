#include "coverage_certificate_2/kernel_2.h"

#include "coverage_certificate_2/errors_2.h"
#include "exact/errors.h"

#include <cmath>

namespace compas_cgal::coverage_certificate {

namespace {
using Traits = CGAL::Fraction_traits<FT>;
using Integer = Traits::Numerator_type;

Integer integer_from_text(const std::string& text)
{
    if (text.empty()) throw InvalidCertificateInputError("empty integer text");
    const std::size_t start = text[0] == '-' ? 1 : 0;
    if (start == text.size()) throw InvalidCertificateInputError("integer text has no digits");
    for (std::size_t i = start; i < text.size(); ++i) {
        if (text[i] < '0' || text[i] > '9') throw InvalidCertificateInputError("integer text must be signed decimal");
    }
    return Integer(text);
}
}  // namespace

FT inject_binary64(const double value)
{
    if (!std::isfinite(value)) {
        throw compas_cgal::exact::NonFiniteBinary64Error("coverage certificate inputs must be finite binary64");
    }
    // cpp_rational's double constructor stores the dyadic value exactly.
    return FT(value);
}

FT rational_from_text(const std::string& numerator, const std::string& denominator)
{
    const Integer num = integer_from_text(numerator);
    const Integer den = integer_from_text(denominator);
    if (CGAL::sign(den) != CGAL::POSITIVE) throw InvalidCertificateInputError("rational denominator must be positive");
    return Traits::Compose()(num, den);
}

RationalText rational_to_text(const FT& value)
{
    Integer num;
    Integer den;
    Traits::Decompose()(value, num, den);
    return {num.str(), den.str()};
}

}  // namespace compas_cgal::coverage_certificate
