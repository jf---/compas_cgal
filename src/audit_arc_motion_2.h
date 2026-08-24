#pragma once

#include "audit_digest_2.h"
#include "exact_circle_chart_2.h"

#include <stdexcept>
#include <string>
#include <vector>

class AuditArcMotionError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class AuditArcNonFiniteInputError : public AuditArcMotionError {
public:
    using AuditArcMotionError::AuditArcMotionError;
};

class AuditArcRadiusIncidenceError : public AuditArcMotionError {
public:
    using AuditArcMotionError::AuditArcMotionError;
};

class AuditArcSweepExtentError : public AuditArcMotionError {
public:
    using AuditArcMotionError::AuditArcMotionError;
};

class AuditArcTraversalError : public AuditArcMotionError {
public:
    using AuditArcMotionError::AuditArcMotionError;
};

class AuditArcAngleNormalizationError : public AuditArcTraversalError {
public:
    using AuditArcTraversalError::AuditArcTraversalError;
};

class AuditArcIntervalDirectionError : public AuditArcTraversalError {
public:
    using AuditArcTraversalError::AuditArcTraversalError;
};

class AuditArcSeamOwnershipError : public AuditArcTraversalError {
public:
    using AuditArcTraversalError::AuditArcTraversalError;
};

class AuditArcEndpointCollapseError : public AuditArcMotionError {
public:
    using AuditArcMotionError::AuditArcMotionError;
};

class AuditArcOrientationError : public AuditArcMotionError {
public:
    using AuditArcMotionError::AuditArcMotionError;
};

class ExactArcChartInterval2 {
public:
    static ExactArcChartInterval2 build(
        int chart,
        const Epeck::FT& start_parameter,
        const Epeck::FT& end_parameter,
        bool increasing,
        bool owns_start_seam,
        bool owns_end_seam);

    int chart() const noexcept;
    const Epeck::FT& start_parameter() const noexcept;
    const Epeck::FT& end_parameter() const noexcept;
    bool increasing() const noexcept;
    bool owns_start_seam() const noexcept;
    bool owns_end_seam() const noexcept;

private:
    ExactArcChartInterval2(
        int chart,
        Epeck::FT start_parameter,
        Epeck::FT end_parameter,
        bool increasing,
        bool owns_start_seam,
        bool owns_end_seam);

    int chart_;
    Epeck::FT start_parameter_;
    Epeck::FT end_parameter_;
    bool increasing_;
    bool owns_start_seam_;
    bool owns_end_seam_;
};

class AuditArcMotion2 {
public:
    static AuditArcMotion2 build(
        const EPoint& center,
        const EVector& zero_phase,
        const Epeck::FT& guide_radius,
        double authored_start_angle,
        double authored_end_angle,
        bool clockwise,
        const Epeck::FT& cut_z);

    const EPoint& center() const noexcept;
    const EVector& zero_phase() const noexcept;
    const Epeck::FT& guide_radius() const noexcept;
    const Epeck::FT& sweep() const noexcept;
    bool clockwise() const noexcept;
    bool full_turn() const noexcept;
    const Epeck::FT& cut_z() const noexcept;
    const std::vector<ExactArcChartInterval2>& intervals() const noexcept;
    const ExactCircleChartParameter2& start_parameter() const noexcept;
    const ExactCircleChartParameter2& end_parameter() const noexcept;
    const EPoint& start_point() const noexcept;
    const EPoint& end_point() const noexcept;
    const NativeMotionDigest2& digest() const noexcept;

private:
    AuditArcMotion2(
        EPoint center,
        EVector zero_phase,
        Epeck::FT guide_radius,
        Epeck::FT sweep,
        bool clockwise,
        bool full_turn,
        Epeck::FT cut_z,
        std::vector<ExactArcChartInterval2> intervals,
        ExactCircleChartParameter2 start_parameter,
        ExactCircleChartParameter2 end_parameter,
        EPoint start_point,
        EPoint end_point,
        NativeMotionDigest2 digest);

    EPoint center_;
    EVector zero_phase_;
    Epeck::FT guide_radius_;
    Epeck::FT sweep_;
    bool clockwise_;
    bool full_turn_;
    Epeck::FT cut_z_;
    std::vector<ExactArcChartInterval2> intervals_;
    ExactCircleChartParameter2 start_parameter_;
    ExactCircleChartParameter2 end_parameter_;
    EPoint start_point_;
    EPoint end_point_;
    NativeMotionDigest2 digest_;
};

bool exact_arc_point_is_incident(
    const AuditArcMotion2& motion,
    const EPoint& point);

const std::string& audit_arc_motion_strategy_version();
