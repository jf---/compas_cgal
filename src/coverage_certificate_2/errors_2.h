#pragma once

#include <stdexcept>

namespace compas_cgal::coverage_certificate {

/// Malformed native input that the typed Python facade should have rejected
/// (wrong array shapes, non-positive radii, malformed rational text).
class InvalidCertificateInputError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

/// The centre domain violates a polygon invariant: a loop is not simple, the
/// outer loop is not CCW, a hole is not CW, a hole is not strictly inside the
/// outer loop, loops touch, or an edge has zero length.
class InvalidCentreDomainError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

/// An internal contract breach: the cover counters disagree with the exact
/// witness membership, a certified enclosure does not contain its exact value,
/// the counters fail parity around a circle, or a coincidence appears that the
/// crossing stage did not predict. Always raised, never swallowed.
class CoverageCertificateNativeError : public std::logic_error {
public:
    using std::logic_error::logic_error;
};

/// The circle-sweep start-witness search exceeded its theorem bound
/// 1 + 2 (2 |counted| + |target circles| + |lines|): every rational witness tried
/// lay on a support, so some support coincides with the circle without being
/// recognised as coincident (an uncut coincident support).
class StartWitnessSearchError : public CoverageCertificateNativeError {
public:
    using CoverageCertificateNativeError::CoverageCertificateNativeError;
};

/// An exact sign rule returned a sign that its own certified interval
/// enclosure excludes (the enclosure does not contain 0 and has the other sign):
/// the rule is wrong. Raised at the first such comparison, never deferred.
class SignRuleContractError : public CoverageCertificateNativeError {
public:
    using CoverageCertificateNativeError::CoverageCertificateNativeError;
};

/// An exact search loop (P3 chart bisection, cross-base separation, the circle
/// witness search) ran past its root-separation bound: the two values it
/// separates would be closer than any two distinct algebraic numbers of their
/// bit sizes can be, so they are equal, and a decision upstream that called
/// them strictly ordered is wrong.
class SeparationBoundExceededError : public CoverageCertificateNativeError {
public:
    using CoverageCertificateNativeError::CoverageCertificateNativeError;
};

/// A swept piece lies exactly on a target cut support of the exact target
/// (a declared boundary line or circle, or a vertex circle): its target side
/// changes along the support where no support cuts it, which the exact-target
/// certificate does not model. Exact connectors along dOmega* produce this
/// coincidence by design and need their own target pieces.
class CoincidentTargetSupportError : public CoverageCertificateNativeError {
public:
    using CoverageCertificateNativeError::CoverageCertificateNativeError;
};

}  // namespace compas_cgal::coverage_certificate
