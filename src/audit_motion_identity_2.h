#pragma once

#include "audit_arc_motion_2.h"
#include "exact_motion_2.h"

#include <string>

class AuditClassificationFactory2;

class AuditSegmentMotion2 {
public:
    const ExactSegmentMotion2& xy() const noexcept;
    const Epeck::FT& cut_z() const noexcept;
    const Epeck::FT& clearance_z() const noexcept;
    const NativeMotionDigest2& digest() const noexcept;

private:
    AuditSegmentMotion2(
        ExactSegmentMotion2 xy,
        Epeck::FT cut_z,
        Epeck::FT clearance_z,
        NativeMotionDigest2 digest);
    friend class AuditClassificationFactory2;

    ExactSegmentMotion2 xy_;
    Epeck::FT cut_z_;
    Epeck::FT clearance_z_;
    NativeMotionDigest2 digest_;
};

class AuditCircleMotion2 {
public:
    const ExactCircleMotion2& xy() const noexcept;
    const Epeck::FT& guide_radius() const noexcept;
    const Epeck::FT& cut_z() const noexcept;
    const Epeck::FT& clearance_z() const noexcept;
    const NativeMotionDigest2& digest() const noexcept;

private:
    AuditCircleMotion2(
        ExactCircleMotion2 xy,
        Epeck::FT guide_radius,
        Epeck::FT cut_z,
        Epeck::FT clearance_z,
        NativeMotionDigest2 digest);
    friend class AuditClassificationFactory2;

    ExactCircleMotion2 xy_;
    Epeck::FT guide_radius_;
    Epeck::FT cut_z_;
    Epeck::FT clearance_z_;
    NativeMotionDigest2 digest_;
};

class AuditVerticalPlunge2 {
public:
    const EPoint& cut_endpoint() const noexcept;
    const Epeck::FT& cut_z() const noexcept;
    const Epeck::FT& clearance_z() const noexcept;
    const NativeMotionDigest2& digest() const noexcept;

private:
    AuditVerticalPlunge2(
        EPoint cut_endpoint,
        Epeck::FT cut_z,
        Epeck::FT clearance_z,
        NativeMotionDigest2 digest);
    friend class AuditClassificationFactory2;

    EPoint cut_endpoint_;
    Epeck::FT cut_z_;
    Epeck::FT clearance_z_;
    NativeMotionDigest2 digest_;
};

class AuditVerticalRetract2 {
public:
    const EPoint& cut_endpoint() const noexcept;
    const Epeck::FT& cut_z() const noexcept;
    const Epeck::FT& clearance_z() const noexcept;
    const NativeMotionDigest2& digest() const noexcept;

private:
    AuditVerticalRetract2(
        EPoint cut_endpoint,
        Epeck::FT cut_z,
        Epeck::FT clearance_z,
        NativeMotionDigest2 digest);
    friend class AuditClassificationFactory2;

    EPoint cut_endpoint_;
    Epeck::FT cut_z_;
    Epeck::FT clearance_z_;
    NativeMotionDigest2 digest_;
};

class AuditClearanceTransport2 {
public:
    const EPoint& start() const noexcept;
    const EPoint& end() const noexcept;
    const Epeck::FT& cut_z() const noexcept;
    const Epeck::FT& clearance_z() const noexcept;
    const NativeMotionDigest2& digest() const noexcept;

private:
    AuditClearanceTransport2(
        EPoint start,
        EPoint end,
        Epeck::FT cut_z,
        Epeck::FT clearance_z,
        NativeMotionDigest2 digest);
    friend class AuditClassificationFactory2;

    EPoint start_;
    EPoint end_;
    Epeck::FT cut_z_;
    Epeck::FT clearance_z_;
    NativeMotionDigest2 digest_;
};

std::string canonical_audit_segment_motion_bytes(
    const ExactSegmentMotion2& motion,
    const Epeck::FT& cut_z,
    const Epeck::FT& clearance_z);

std::string canonical_audit_circle_motion_bytes(
    const ExactCircleMotion2& motion,
    const Epeck::FT& guide_radius,
    const Epeck::FT& cut_z,
    const Epeck::FT& clearance_z);

std::string canonical_audit_vertical_plunge_bytes(
    const EPoint& cut_endpoint,
    const Epeck::FT& cut_z,
    const Epeck::FT& clearance_z);

std::string canonical_audit_vertical_retract_bytes(
    const EPoint& cut_endpoint,
    const Epeck::FT& cut_z,
    const Epeck::FT& clearance_z);

std::string canonical_audit_clearance_line_bytes(
    const EPoint& start,
    const EPoint& end,
    const Epeck::FT& cut_z,
    const Epeck::FT& clearance_z);

std::string canonical_audit_clearance_circle_bytes(
    const ExactCircleMotion2& motion,
    const Epeck::FT& guide_radius,
    const Epeck::FT& cut_z,
    const Epeck::FT& clearance_z);

std::string canonical_audit_clearance_arc_bytes(
    const AuditArcMotion2& motion,
    const Epeck::FT& cut_z,
    const Epeck::FT& clearance_z);
