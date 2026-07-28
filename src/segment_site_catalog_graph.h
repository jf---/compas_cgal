#pragma once

#include "reachable_input_2.h"
#include "segment_site_mat.h"

#include <stdexcept>

class UnsupportedCanonicalMatRectangleGraphError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class IncompleteCanonicalMatRectangleGraphError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class InvalidCanonicalMatRectangleNodeError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class UnknownCanonicalMatGraphSourceError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

MatExactGraph2 canonical_rectangle_mat_graph(
    const CanonicalReachInput2& input,
    const CORE::BigRat& radius_squared);
