#pragma once

#include "stock_2.h"

#include <array>
#include <stdexcept>
#include <utility>

class InvalidHeldContourInputError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};
class UndefinedPredecessorDirectionError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};
class NoExposedPredecessorArcError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Held-Pfeiffer Section 3.1's union of filled outer disks, in world XY mm.
// This is a placement model, distinct from Stock2's physical annular sweeps.
class HeldDiskContour2 {
public:
    using XY = std::array<double, 2>;
    HeldDiskContour2(const XY& center, double guide_radius, double tool_radius);
    void append(const XY& center, double guide_radius);

    // Exact b on the latest outer circle, moved CW to the exposed contour when
    // covered. The flag records whether b moved. Keep this point native for any
    // subsequent deciding geometry; the Python binding reports coordinates only.
    std::pair<GpsPoint, bool> contact_toward(const XY& candidate_center) const;

private:
    HeldDiskContour2(const EPoint& center, const Epeck::FT& guide_radius,
                     const Epeck::FT& tool_radius);
    Epeck::FT tool_radius_;
    ECircle predecessor_;
    Gps contour_;
};
