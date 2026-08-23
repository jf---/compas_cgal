#include "audit_classification_2.h"

#include <CGAL/enum.h>
#include <CGAL/number_utils.h>

#include <cmath>
#include <initializer_list>
#include <numbers>
#include <string>
#include <string_view>

#include <nanobind/stl/array.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace {

constexpr double FULL_TURN_SURROGATE = std::numbers::pi * 2.0;
constexpr std::string_view ARC_PHASE_STRATEGY_VERSION
    = "audit-arc-phase-binary64-v1";

enum class CurvePlane {
    Cut,
    Clearance,
};

void require_finite(const std::initializer_list<double> values)
{
    for (const double value : values) {
        if (!std::isfinite(value)) {
            throw AuditNonFiniteInputError(
                "audit classification inputs must be finite binary64 values");
        }
    }
}

void require_finite(const std::array<double, 3>& values)
{
    require_finite({ values[0], values[1], values[2] });
}

void require_valid_plane(double cut_z, double clearance_z)
{
    require_finite({ cut_z, clearance_z });
    if (CGAL::compare(Epeck::FT(clearance_z), Epeck::FT(cut_z)) != CGAL::LARGER) {
        throw AuditInvalidPlaneError(
            "audit clearance Z must be exactly greater than cut Z");
    }
}

void require_known_role(const std::string& role)
{
    if (role != "cut" && role != "lead_in" && role != "lead_out"
        && role != "link" && role != "retract" && role != "plunge") {
        throw AuditContradictoryRoleError(
            "operation role is outside the closed audit domain");
    }
}

void require_role(const std::string& observed, const std::string& expected)
{
    require_known_role(observed);
    if (observed != expected) {
        throw AuditContradictoryRoleError(
            "operation role contradicts exact geometry classification");
    }
}

void require_lateral_role(const std::string& role)
{
    require_known_role(role);
    if (role == "plunge" || role == "retract") {
        throw AuditContradictoryRoleError(
            "vertical operation role cannot label exact lateral motion");
    }
}

void require_canonical_world_xy_frame(
    const std::array<double, 3>& xaxis,
    const std::array<double, 3>& yaxis)
{
    const EVector exact_xaxis {
        Epeck::FT(xaxis[0]), Epeck::FT(xaxis[1])
    };
    const EVector exact_yaxis {
        Epeck::FT(yaxis[0]), Epeck::FT(yaxis[1])
    };
    const Epeck::FT one(1);
    if (CGAL::sign(Epeck::FT(xaxis[2])) != CGAL::ZERO
        || CGAL::sign(Epeck::FT(yaxis[2])) != CGAL::ZERO
        || CGAL::compare(exact_xaxis.squared_length(), one) != CGAL::EQUAL
        || CGAL::compare(exact_yaxis.squared_length(), one) != CGAL::EQUAL
        || CGAL::sign(exact_xaxis * exact_yaxis) != CGAL::ZERO) {
        throw AuditUnsupportedGeometryError(
            "circular motion requires an exact orthonormal world-XY frame");
    }
    const EPoint origin(Epeck::FT(0), Epeck::FT(0));
    const EPoint x_axis_point(exact_xaxis.x(), exact_xaxis.y());
    const EPoint y_axis_point(exact_yaxis.x(), exact_yaxis.y());
    if (CGAL::orientation(origin, x_axis_point, y_axis_point)
        != CGAL::LEFT_TURN) {
        throw AuditUnsupportedGeometryError(
            "circular motion requires a positive-normal world-XY frame");
    }
}

CurvePlane classify_curve_plane(
    const Epeck::FT& center_z,
    const Epeck::FT& cut_z,
    const Epeck::FT& clearance_z,
    const std::string& operation_role)
{
    if (CGAL::compare(center_z, cut_z) == CGAL::EQUAL) {
        require_lateral_role(operation_role);
        return CurvePlane::Cut;
    }
    if (CGAL::compare(center_z, clearance_z) == CGAL::EQUAL) {
        require_role(operation_role, "link");
        return CurvePlane::Clearance;
    }
    throw AuditOffPlaneError(
        "circular motion lies on neither exact declared plane");
}

EPoint translated_point(const EPoint& center, const EVector& phase)
{
    return center + phase;
}

EVector exact_circle_phase(
    const std::array<double, 3>& xaxis,
    double radius)
{
    const Epeck::FT exact_radius(radius);
    return EVector(
        exact_radius * Epeck::FT(xaxis[0]),
        exact_radius * Epeck::FT(xaxis[1]));
}

EVector exact_arc_phase_surrogate(
    const std::array<double, 3>& xaxis,
    const std::array<double, 3>& yaxis,
    double radius,
    double angle)
{
    // ARC_PHASE_STRATEGY_VERSION names this sole transcendental seam. The
    // binary64 phase is injected once and remains opaque Epeck geometry.
    const double cosine = std::cos(angle);
    const double sine = std::sin(angle);
    return EVector(
        Epeck::FT(radius * (cosine * xaxis[0] + sine * yaxis[0])),
        Epeck::FT(radius * (cosine * xaxis[1] + sine * yaxis[1])));
}

} // namespace

AuditLineClassification2 classify_audit_line(
    const std::array<double, 3>& start,
    const std::array<double, 3>& end,
    double cut_z,
    double clearance_z,
    const std::string& operation_role)
{
    require_finite(start);
    require_finite(end);
    require_valid_plane(cut_z, clearance_z);
    require_known_role(operation_role);

    const EPoint start_xy { Epeck::FT(start[0]), Epeck::FT(start[1]) };
    const EPoint end_xy { Epeck::FT(end[0]), Epeck::FT(end[1]) };
    const Epeck::FT start_z(start[2]);
    const Epeck::FT end_z(end[2]);
    const Epeck::FT exact_cut_z(cut_z);
    const Epeck::FT exact_clearance_z(clearance_z);

    if (start_xy == end_xy) {
        if (CGAL::compare(start_z, exact_clearance_z) == CGAL::EQUAL
            && CGAL::compare(end_z, exact_cut_z) == CGAL::EQUAL) {
            require_role(operation_role, "plunge");
            return AuditVerticalPlunge2 { end_xy };
        }
        if (CGAL::compare(start_z, exact_cut_z) == CGAL::EQUAL
            && CGAL::compare(end_z, exact_clearance_z) == CGAL::EQUAL) {
            require_role(operation_role, "retract");
            return AuditVerticalRetract2 { start_xy };
        }
        throw AuditUnsupportedGeometryError(
            "vertical motion must connect the exact declared planes");
    }
    if (CGAL::compare(start_z, end_z) != CGAL::EQUAL) {
        throw AuditUnsupportedGeometryError(
            "mixed-Z XY ramp has no Stage 1 audit semantics");
    }
    if (CGAL::compare(start_z, exact_cut_z) == CGAL::EQUAL) {
        require_lateral_role(operation_role);
        return AuditSegmentMotion2 {
            ExactSegmentMotion2 { start_xy, end_xy },
            exact_cut_z,
        };
    }
    if (CGAL::compare(start_z, exact_clearance_z) == CGAL::EQUAL) {
        require_role(operation_role, "link");
        return AuditClearanceTransport2 { start_xy, end_xy };
    }
    throw AuditOffPlaneError(
        "lateral motion lies on neither exact declared plane");
}

AuditCircleClassification2 classify_audit_circle(
    const std::array<double, 3>& center,
    const std::array<double, 3>& xaxis,
    const std::array<double, 3>& yaxis,
    double radius,
    bool clockwise,
    double cut_z,
    double clearance_z,
    const std::string& operation_role)
{
    require_finite(center);
    require_finite(xaxis);
    require_finite(yaxis);
    require_finite({ radius });
    require_valid_plane(cut_z, clearance_z);
    const Epeck::FT exact_radius(radius);
    if (CGAL::sign(exact_radius) != CGAL::POSITIVE) {
        throw AuditUnsupportedGeometryError(
            "circle radius must be exact positive");
    }
    require_canonical_world_xy_frame(xaxis, yaxis);
    const EPoint exact_center {
        Epeck::FT(center[0]), Epeck::FT(center[1])
    };
    const EVector phase = exact_circle_phase(xaxis, radius);
    const CurvePlane plane = classify_curve_plane(
        Epeck::FT(center[2]), Epeck::FT(cut_z),
        Epeck::FT(clearance_z), operation_role);
    if (plane == CurvePlane::Cut) {
        return AuditCircleMotion2 {
            ExactCircleMotion2 { exact_center, phase, clockwise },
            exact_radius,
            Epeck::FT(cut_z),
        };
    }
    const EPoint seam = translated_point(exact_center, phase);
    return AuditClearanceTransport2 { seam, seam };
}

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
    const std::string& operation_role)
{
    require_finite(center);
    require_finite(xaxis);
    require_finite(yaxis);
    require_finite({ radius, start_angle, end_angle });
    require_valid_plane(cut_z, clearance_z);
    const Epeck::FT exact_radius(radius);
    const Epeck::FT signed_sweep
        = Epeck::FT(end_angle) - Epeck::FT(start_angle);
    if (CGAL::sign(exact_radius) != CGAL::POSITIVE) {
        throw AuditUnsupportedGeometryError(
            "arc radius must be exact positive");
    }
    const CGAL::Sign sweep_sign = CGAL::sign(signed_sweep);
    if (sweep_sign == CGAL::ZERO
        || CGAL::compare(
               CGAL::abs(signed_sweep),
               Epeck::FT(FULL_TURN_SURROGATE))
            == CGAL::LARGER) {
        throw AuditUnsupportedGeometryError(
            "arc sweep must lie in the exact injected interval [-2*pi, 2*pi] excluding zero");
    }
    if (clockwise != (sweep_sign == CGAL::NEGATIVE)) {
        throw AuditContradictoryOrientationError(
            "operation orientation contradicts exact arc sweep sign");
    }
    require_canonical_world_xy_frame(xaxis, yaxis);
    const EPoint exact_center {
        Epeck::FT(center[0]), Epeck::FT(center[1])
    };
    const EVector start_phase = exact_arc_phase_surrogate(
        xaxis, yaxis, radius, start_angle);
    const CurvePlane plane = classify_curve_plane(
        Epeck::FT(center[2]), Epeck::FT(cut_z),
        Epeck::FT(clearance_z), operation_role);
    if (plane == CurvePlane::Cut) {
        return AuditArcMotion2 {
            exact_center,
            start_phase,
            exact_radius,
            CGAL::abs(signed_sweep),
            clockwise,
            Epeck::FT(cut_z),
        };
    }
    const EVector end_phase = exact_arc_phase_surrogate(
        xaxis, yaxis, radius, end_angle);
    return AuditClearanceTransport2 {
        translated_point(exact_center, start_phase),
        translated_point(exact_center, end_phase),
    };
}

void register_audit_classification_2(nb::module_& module)
{
    nb::exception<AuditNonFiniteInputError>(
        module, "AuditNonFiniteInputError", PyExc_ValueError);
    nb::exception<AuditInvalidPlaneError>(
        module, "AuditInvalidPlaneError", PyExc_ValueError);
    nb::exception<AuditUnsupportedGeometryError> unsupported(
        module, "AuditUnsupportedGeometryError", PyExc_ValueError);
    nb::exception<AuditOffPlaneError>(
        module, "AuditOffPlaneError", unsupported.ptr());
    nb::exception<AuditContradictoryRoleError>(
        module, "AuditContradictoryRoleError", PyExc_ValueError);
    nb::exception<AuditContradictoryOrientationError>(
        module, "AuditContradictoryOrientationError", PyExc_ValueError);

    nb::class_<AuditSegmentMotion2>(module, "AuditSegmentMotion2");
    nb::class_<AuditCircleMotion2>(module, "AuditCircleMotion2");
    nb::class_<AuditArcMotion2>(module, "AuditArcMotion2");
    nb::class_<AuditVerticalPlunge2>(module, "AuditVerticalPlunge2");
    nb::class_<AuditVerticalRetract2>(module, "AuditVerticalRetract2");
    nb::class_<AuditClearanceTransport2>(
        module, "AuditClearanceTransport2");

    module.def(
        "audit_arc_phase_strategy_version",
        []() {
            return nb::bytes(
                ARC_PHASE_STRATEGY_VERSION.data(),
                ARC_PHASE_STRATEGY_VERSION.size());
        });

    module.def(
        "classify_audit_line", &classify_audit_line,
        "start"_a, "end"_a, "cut_z"_a, "clearance_z"_a,
        "operation_role"_a);
    module.def(
        "classify_audit_circle", &classify_audit_circle,
        "center"_a, "xaxis"_a, "yaxis"_a, "radius"_a,
        "clockwise"_a, "cut_z"_a, "clearance_z"_a,
        "operation_role"_a);
    module.def(
        "classify_audit_arc", &classify_audit_arc,
        "center"_a, "xaxis"_a, "yaxis"_a, "radius"_a,
        "start_angle"_a, "end_angle"_a, "clockwise"_a,
        "cut_z"_a, "clearance_z"_a, "operation_role"_a);
}
