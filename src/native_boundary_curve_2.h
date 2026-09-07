#pragma once

#include "reachable_arrangement_2.h"
#include "boundary_normal_circle_2.h"
#include <array>

class InvalidNativeBoundaryCurveError : public std::runtime_error {
public: using std::runtime_error::runtime_error;
};
class InvalidNativeBoundaryChainError : public std::runtime_error {
public: using std::runtime_error::runtime_error;
};

class InvalidNativeBoundaryMedialInputError : public std::runtime_error {
public: using std::runtime_error::runtime_error;
};
class NoPositiveNativeBoundaryCircleError : public std::runtime_error {
public: using std::runtime_error::runtime_error;
};
class NativeBoundaryMedialConstructionError : public std::runtime_error {
public: using std::runtime_error::runtime_error;
};

// World-XY millimetre geometry. Arc fitting preserves endpoints, not tangents
// or a source-error budget. Its displacement property is reporting only.
class NativeBoundaryCurve2 {
public:
    using XY = std::array<double, 2>;
    static NativeBoundaryCurve2 line(const XY& start, const XY& end);
    static NativeBoundaryCurve2 arc(const XY& start, const XY& end, const XY& center, bool counterclockwise);
    const ReachCurve& curve() const { return curve_; }
    ReachPoint start() const { return curve_.source(); }
    ReachPoint end() const { return curve_.target(); }
    double center_adjustment_mm() const;
    XY center_mm() const;
    double radius_mm() const;
    bool is_arc() const { return curve_.is_circular(); }
private:
    NativeBoundaryCurve2(ReachCurve curve, ReachFT adjustment_squared)
        : curve_(std::move(curve)), adjustment_squared_(std::move(adjustment_squared)) {}
    ReachCurve curve_;
    ReachFT adjustment_squared_;
};

class NativeBoundary2 {
public:
    explicit NativeBoundary2(std::vector<NativeBoundaryCurve2> curves);
    const std::vector<NativeBoundaryCurve2>& curves() const { return curves_; }
    const ReachableBoundaryCycle2& cycle() const { return cycle_; }
    ExactRegion2 design_region() const;
    boundary_normal::BoundaryNormalCircleProposal2 circle_on_piece(
        std::int64_t piece_index, double parameter, double tool_radius) const;
private:
    std::vector<NativeBoundaryCurve2> curves_;
    ReachableBoundaryCycle2 cycle_;
    ReachSet design_;
};
