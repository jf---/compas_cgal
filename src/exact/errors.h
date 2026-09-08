#pragma once

#include <stdexcept>

namespace compas_cgal::exact {

/// Raised when a double that must denote an exact rational is NaN or infinite.
class NonFiniteBinary64Error : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/// Raised when `to_canonical` decomposes a rational into a non-positive
/// denominator, which would leave the attestation bytes ambiguous. That is the
/// only condition raised today. Reducedness is relied on for the same reason
/// but is not re-checked per value: it comes from the backend rational's
/// auto-normalisation, not from CGAL's `Fraction_traits` concept, which
/// guarantees only `value == num / den`. The probes in `exact_canonical_gate`
/// pin that normalisation instead, because a per-value bignum gcd would cost
/// more than it can ever catch.
class UnreducedCanonicalRationalError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/// Raised when one-root arithmetic is attempted across distinct roots. CGAL
/// documents this as a precondition; violating it is undefined behaviour.
class CrossRootExtensionError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/// Raised when a projection produced attestation bytes differing from the
/// frozen contract. Declared here so the module's error model is complete; it
/// has no raiser yet, and gets one when the projections are routed through
/// `to_canonical`.
class AttestationByteDriftError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

}  // namespace compas_cgal::exact
