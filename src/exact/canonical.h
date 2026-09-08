#pragma once

#include "exact/rational.h"

#include <string>

namespace compas_cgal::exact {

/// A rational in canonical attestation form: reduced, positive denominator.
///
/// This is a derived VIEW of an exact value, never a carrier. It exists so that
/// attestation bytes are a function of the mathematical value rather than of
/// whichever number type produced it. Construct it only via `to_canonical`.
class CanonicalRational {
public:
    [[nodiscard]] const std::string& numerator() const noexcept;
    [[nodiscard]] const std::string& denominator() const noexcept;

    /// Decimal text, "n" when the denominator is 1 and "n/d" otherwise.
    [[nodiscard]] std::string text() const;

    /// Frozen attestation bytes. This encoding is a compatibility contract:
    /// changing it invalidates every stored replay digest.
    [[nodiscard]] std::string canonical_bytes() const;

private:
    friend CanonicalRational to_canonical(const Rational& value);

    CanonicalRational(std::string numerator, std::string denominator);

    std::string numerator_;
    std::string denominator_;
};

/// Project an exact rational into canonical attestation form.
///
/// This is the ONLY exact-to-attestation-bytes exit point in the codebase.
///
/// Raises:
///     UnreducedCanonicalRationalError: if the decomposed denominator is not
///         positive, which would make the encoding ambiguous.
[[nodiscard]] CanonicalRational to_canonical(const Rational& value);

}  // namespace compas_cgal::exact
