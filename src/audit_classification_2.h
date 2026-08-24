#pragma once

#include "audit_arc_motion_2.h"
#include "exact_motion_2.h"

#include <array>
#include <stdexcept>
#include <string>
#include <variant>

#include <nanobind/nanobind.h>

class AuditNonFiniteInputError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class AuditInvalidPlaneError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class AuditUnsupportedGeometryError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class AuditOffPlaneError : public AuditUnsupportedGeometryError {
public:
    using AuditUnsupportedGeometryError::AuditUnsupportedGeometryError;
};

class AuditContradictoryRoleError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class AuditContradictoryOrientationError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

struct AuditSegmentMotion2 {
    ExactSegmentMotion2 xy;
    Epeck::FT cut_z;
};

struct AuditCircleMotion2 {
    ExactCircleMotion2 xy;
    Epeck::FT guide_radius;
    Epeck::FT cut_z;
};

struct AuditVerticalPlunge2 {
    EPoint cut_endpoint;
};

struct AuditVerticalRetract2 {
    EPoint cut_endpoint;
};

struct AuditClearanceTransport2 {
    EPoint start;
    EPoint end;
};

using AuditLineClassification2 = std::variant<
    AuditSegmentMotion2,
    AuditVerticalPlunge2,
    AuditVerticalRetract2,
    AuditClearanceTransport2>;
using AuditCircleClassification2 = std::variant<AuditCircleMotion2, AuditClearanceTransport2>;
using AuditArcClassification2 = std::variant<AuditArcMotion2, AuditClearanceTransport2>;

AuditLineClassification2 classify_audit_line(
    const std::array<double, 3>& start,
    const std::array<double, 3>& end,
    double cut_z,
    double clearance_z,
    const std::string& operation_role);

AuditCircleClassification2 classify_audit_circle(
    const std::array<double, 3>& center,
    const std::array<double, 3>& xaxis,
    const std::array<double, 3>& yaxis,
    double radius,
    bool clockwise,
    double cut_z,
    double clearance_z,
    const std::string& operation_role);

AuditArcClassification2 classify_audit_arc(
    const std::array<double, 3>& center,
    const std::array<double, 3>& xaxis,
    const std::array<double, 3>& yaxis,
    double radius,
    double start_angle,
    double end_angle,
    bool clockwise,
    double cut_z,
    double clearance_z,
    const std::string& operation_role);

void register_audit_classification_2(nanobind::module_& module);
