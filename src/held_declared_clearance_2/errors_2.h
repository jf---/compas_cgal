#pragma once

// Named failure modes of the native declared-boundary clearance oracle
// (held_declared_clearance_2). One exception per failure mode; the bindings
// map each to a Python exception of the same name.

#include <stdexcept>

namespace compas_cgal::held_declared_clearance {

/// Base of every failure of the oracle.
class DeclaredClearanceError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/// A declared primitive breaks its contract: a line whose direction is not
/// exactly unit or whose end is not start + length * u, an arc whose ends are
/// not exactly on its circle, a non-positive length or radius, a sense not +-1.
class InvalidDeclaredPrimitiveError : public DeclaredClearanceError {
public:
    using DeclaredClearanceError::DeclaredClearanceError;
};

/// The tool radius is not a positive binary64, or a circle radius is not positive.
class NonPositiveClearanceRadiusError : public DeclaredClearanceError {
public:
    using DeclaredClearanceError::DeclaredClearanceError;
};

/// A segment element has no direction, or an arc element no radius or no sense.
class DegenerateClearanceElementError : public DeclaredClearanceError {
public:
    using DeclaredClearanceError::DeclaredClearanceError;
};

/// A one-root number was given a negative radicand.
class NegativeRadicandError : public DeclaredClearanceError {
public:
    using DeclaredClearanceError::DeclaredClearanceError;
};

/// Arithmetic between numbers of different quadratic extensions was attempted
/// (an internal contract breach: comparisons across roots are supported,
/// arithmetic is not; docs/exactness.md, "One-root numbers").
class CrossRootArithmeticError : public DeclaredClearanceError {
public:
    using DeclaredClearanceError::DeclaredClearanceError;
};

/// A one-root point's coordinates do not share one quadratic extension.
class MixedRootPointError : public DeclaredClearanceError {
public:
    using DeclaredClearanceError::DeclaredClearanceError;
};

/// A wire record is malformed: wrong arity, a non-positive denominator, a
/// non-finite binary64, an owner index out of range.
class InvalidClearanceWireError : public DeclaredClearanceError {
public:
    using DeclaredClearanceError::DeclaredClearanceError;
};

/// A reported binary64 would overflow (reporting only; never a decision).
class ReportingOverflowError : public DeclaredClearanceError {
public:
    using DeclaredClearanceError::DeclaredClearanceError;
};

}  // namespace compas_cgal::held_declared_clearance
