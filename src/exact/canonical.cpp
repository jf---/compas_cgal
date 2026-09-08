#include "exact/canonical.h"

#include "exact/errors.h"
#include "continuous_tea_2/event_certificate.h"

#include <CGAL/CORE/BigInt.h>
#include <CGAL/Fraction_traits.h>

#include <utility>
#include <vector>

namespace compas_cgal::exact {

CanonicalRational::CanonicalRational(
    std::string numerator,
    std::string denominator)
    : numerator_(std::move(numerator)),
      denominator_(std::move(denominator))
{
}

const std::string& CanonicalRational::numerator() const noexcept
{
    return numerator_;
}

const std::string& CanonicalRational::denominator() const noexcept
{
    return denominator_;
}

std::string CanonicalRational::text() const
{
    return denominator_ == "1" ? numerator_ : numerator_ + "/" + denominator_;
}

std::string CanonicalRational::canonical_bytes() const
{
    return encode_string_sequence(
        {
            "exact-binary64-rational-v1",
            numerator_,
            denominator_,
        });
}

CanonicalRational to_canonical(const Rational& value)
{
    using Traits = CGAL::Fraction_traits<Rational>;
    typename Traits::Numerator_type numerator;
    typename Traits::Denominator_type denominator;
    typename Traits::Decompose()(value, numerator, denominator);

    const CORE::BigInt exact_denominator(denominator.exact());
    if (exact_denominator <= 0) {
        throw UnreducedCanonicalRationalError(
            "canonical rational requires a positive denominator");
    }
    return CanonicalRational(
        CORE::BigInt(numerator.exact()).convert_to<std::string>(),
        exact_denominator.convert_to<std::string>());
}

}  // namespace compas_cgal::exact
