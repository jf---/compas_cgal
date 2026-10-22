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

/// A swept piece lies exactly on a support where no exact side probe exists: an
/// irrational (polygon connector) flank on a support of an exact connector, whose
/// one-root witness cannot seed the rational probe of side_probe_2.h. Rational
/// coincident supports (exact connectors on dOmega*, cuts on target circles) are
/// decided by the probe and never raise this.
class CoincidentTargetSupportError : public CoverageCertificateNativeError {
public:
    using CoverageCertificateNativeError::CoverageCertificateNativeError;
};

/// The exact connectors do not cover dOmega*: some point of dOmega* lies on no
/// sub-piece, so the cap-circle theorem (exact_connectors_2.h) does not apply and
/// the one-root cap circles would have to be checked as pieces, which is not
/// modelled. Raised instead of certifying with an unproven omission. A member of
/// the native family, so a caller catching CoverageCertificateNativeError sees it.
class ConnectorBoundaryCoverageError : public CoverageCertificateNativeError {
public:
    using CoverageCertificateNativeError::CoverageCertificateNativeError;
};

/// The declared dOmega* pieces do not form closed chains with non-reflex
/// junctions: a piece's end is the start of no piece or of more than one, a
/// piece start is reached twice, or the travel turns right (or reverses) at a
/// junction. The cap-circle theorem (exact_connectors_2.h) needs both, so the
/// exact connectors are refused instead of certified on an unproven boundary.
class BoundaryChainError : public CoverageCertificateNativeError {
public:
    using CoverageCertificateNativeError::CoverageCertificateNativeError;
};

}  // namespace compas_cgal::coverage_certificate
