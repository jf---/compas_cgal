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

class BoundaryFeatureIndexError : public BoundaryExtractionError {
public:
    using BoundaryExtractionError::BoundaryExtractionError;
};
