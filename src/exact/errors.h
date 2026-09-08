#pragma once

#include <stdexcept>

namespace compas_cgal::exact {

/// Raised when a double that must denote an exact rational is NaN or infinite.
class NonFiniteBinary64Error : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/// Raised when a canonical rational violates the reduced / positive-denominator
/// / gcd == 1 invariant that makes attestation bytes value-determined.
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
