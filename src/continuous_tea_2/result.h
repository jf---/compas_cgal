#pragma once

#include <stdexcept>

enum class ContinuousTeaVerdict {
    CERTIFIED,
    CAP_EXCEEDED,
    UNRESOLVED_DEGENERACY,
};

class BoundaryExtractionError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class DegenerateBoundarySupportError
    : public BoundaryExtractionError {
public:
    using BoundaryExtractionError::BoundaryExtractionError;
};

class MissingBoundaryEndpointError
    : public BoundaryExtractionError {
public:
    using BoundaryExtractionError::BoundaryExtractionError;
};

class MissingBoundaryIntersectionError
    : public BoundaryExtractionError {
public:
    using BoundaryExtractionError::BoundaryExtractionError;
};

class BoundaryFeatureIndexError : public BoundaryExtractionError {
public:
    using BoundaryExtractionError::BoundaryExtractionError;
};
