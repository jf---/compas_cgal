#pragma once

#include "boundary_normal_circle_2.h"

#include <stdexcept>
#include <utility>

namespace boundary_normal {

class BoundaryCircleToolMismatchError : public std::runtime_error {
public: using std::runtime_error::runtime_error;
};
class BoundaryCircleSpacingError : public std::runtime_error {
public: using std::runtime_error::runtime_error;
};
class UncoveredStationaryCircleError : public std::runtime_error {
public: using std::runtime_error::runtime_error;
};
class InvalidBoundaryEngagementCapError : public std::runtime_error {
public: using std::runtime_error::runtime_error;
};
class BoundaryCircleEngagementGeometryError : public std::runtime_error {
public: using std::runtime_error::runtime_error;
};

// Standard predecessor model, conditional on its entire outer disk being
// cleared. Eq. 4 failures raise independently of the engagement cap. The bool
// compares exact cosines against the injected squared-chord surrogate [0, 4];
// angle radians are reporting only. Geometrically unresolved closest-contact
// branches conservatively report pi, not a proved exact maximum.
std::pair<double, bool> boundary_circle_engagement(
    const BoundaryNormalCircleProposal2& previous,
    const BoundaryNormalCircleProposal2& current,
    double cap_chord_ratio);

} // namespace boundary_normal
