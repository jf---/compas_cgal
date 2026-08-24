#include "audit_classification_core_2.h"

#include <CGAL/enum.h>
#include <CGAL/number_utils.h>

#include <cmath>
#include <initializer_list>
#include <string>
#include <utility>

class AuditClassificationFactory2 {
public:
    static AuditSegmentMotion2 segment(
        ExactSegmentMotion2 xy,
        const Epeck::FT& cut_z,
        const Epeck::FT& clearance_z)
    {
        NativeMotionDigest2 digest = NativeMotionDigestAuthority2::hash_canonical(
            canonical_audit_segment_motion_bytes(xy, cut_z, clearance_z));
        return AuditSegmentMotion2(
            std::move(xy), cut_z, clearance_z, std::move(digest));
    }

    static AuditCircleMotion2 circle(
        ExactCircleMotion2 xy,
        const Epeck::FT& guide_radius,
        const Epeck::FT& cut_z,
        const Epeck::FT& clearance_z)
    {
        NativeMotionDigest2 digest = NativeMotionDigestAuthority2::hash_canonical(
            canonical_audit_circle_motion_bytes(
                xy, guide_radius, cut_z, clearance_z));
        return AuditCircleMotion2(
            std::move(xy), guide_radius, cut_z, clearance_z,
            std::move(digest));
    }

    static AuditVerticalPlunge2 plunge(
        const EPoint& cut_endpoint,
        const Epeck::FT& cut_z,
        const Epeck::FT& clearance_z)
    {
        return AuditVerticalPlunge2(
            cut_endpoint,
            cut_z,
            clearance_z,
            NativeMotionDigestAuthority2::hash_canonical(
                canonical_audit_vertical_plunge_bytes(
                    cut_endpoint, cut_z, clearance_z)));
    }

    static AuditVerticalRetract2 retract(
        const EPoint& cut_endpoint,
        const Epeck::FT& cut_z,
        const Epeck::FT& clearance_z)
    {
        return AuditVerticalRetract2(
            cut_endpoint,
            cut_z,
            clearance_z,
            NativeMotionDigestAuthority2::hash_canonical(
                canonical_audit_vertical_retract_bytes(
                    cut_endpoint, cut_z, clearance_z)));
    }

    static AuditClearanceTransport2 clearance_line(
        const EPoint& start,
        const EPoint& end,
        const Epeck::FT& cut_z,
        const Epeck::FT& clearance_z)
    {
        return AuditClearanceTransport2(
            start,
            end,
            cut_z,
            clearance_z,
            NativeMotionDigestAuthority2::hash_canonical(
                canonical_audit_clearance_line_bytes(
                    start, end, cut_z, clearance_z)));
    }

    static AuditClearanceTransport2 clearance_circle(
        const ExactCircleMotion2& motion,
        const Epeck::FT& guide_radius,
        const Epeck::FT& cut_z,
        const Epeck::FT& clearance_z)
    {
        const EPoint seam = motion.center + motion.phase_vector;
        return AuditClearanceTransport2(
            seam,
            seam,
            cut_z,
            clearance_z,
            NativeMotionDigestAuthority2::hash_canonical(
                canonical_audit_clearance_circle_bytes(
                    motion, guide_radius, cut_z, clearance_z)));
    }

    static AuditClearanceTransport2 clearance_arc(
        const AuditArcMotion2& motion,
        const Epeck::FT& cut_z,
        const Epeck::FT& clearance_z)
    {
        return AuditClearanceTransport2(
            motion.start_point(),
            motion.end_point(),
            cut_z,
            clearance_z,
            NativeMotionDigestAuthority2::hash_canonical(
                canonical_audit_clearance_arc_bytes(
                    motion, cut_z, clearance_z)));
    }
};

namespace {

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
    require_finite({values[0], values[1], values[2]});
}

void require_valid_plane(double cut_z, double clearance_z)
{
    require_finite({cut_z, clearance_z});
    if (CGAL::compare(Epeck::FT(clearance_z), Epeck::FT(cut_z))
        != CGAL::LARGER) {
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

EVector exact_circle_phase(
    const std::array<double, 3>& xaxis,
    double radius)
{
    const Epeck::FT exact_radius(radius);
    return EVector(
        exact_radius * Epeck::FT(xaxis[0]),
        exact_radius * Epeck::FT(xaxis[1]));
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

    const EPoint start_xy {Epeck::FT(start[0]), Epeck::FT(start[1])};
    const EPoint end_xy {Epeck::FT(end[0]), Epeck::FT(end[1])};
    const Epeck::FT start_z(start[2]);
    const Epeck::FT end_z(end[2]);
    const Epeck::FT exact_cut_z(cut_z);
    const Epeck::FT exact_clearance_z(clearance_z);

    if (start_xy == end_xy) {
        if (CGAL::compare(start_z, exact_clearance_z) == CGAL::EQUAL
            && CGAL::compare(end_z, exact_cut_z) == CGAL::EQUAL) {
            require_role(operation_role, "plunge");
            return AuditClassificationFactory2::plunge(
                end_xy, exact_cut_z, exact_clearance_z);
        }
        if (CGAL::compare(start_z, exact_cut_z) == CGAL::EQUAL
            && CGAL::compare(end_z, exact_clearance_z) == CGAL::EQUAL) {
            require_role(operation_role, "retract");
            return AuditClassificationFactory2::retract(
                start_xy, exact_cut_z, exact_clearance_z);
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
        return AuditClassificationFactory2::segment(
            ExactSegmentMotion2 {start_xy, end_xy},
            exact_cut_z,
            exact_clearance_z);
    }
    if (CGAL::compare(start_z, exact_clearance_z) == CGAL::EQUAL) {
        require_role(operation_role, "link");
        return AuditClassificationFactory2::clearance_line(
            start_xy, end_xy, exact_cut_z, exact_clearance_z);
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
    require_finite({radius});
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
    const Epeck::FT exact_cut_z(cut_z);
    const Epeck::FT exact_clearance_z(clearance_z);
    const CurvePlane plane = classify_curve_plane(
        Epeck::FT(center[2]), exact_cut_z,
        exact_clearance_z, operation_role);
    const ExactCircleMotion2 motion {exact_center, phase, clockwise};
    if (plane == CurvePlane::Cut) {
        return AuditClassificationFactory2::circle(
            motion,
            exact_radius,
            exact_cut_z,
            exact_clearance_z);
    }
    return AuditClassificationFactory2::clearance_circle(
        motion,
        exact_radius,
        exact_cut_z,
        exact_clearance_z);
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
    require_finite({radius, start_angle, end_angle});
    require_valid_plane(cut_z, clearance_z);
    const Epeck::FT exact_radius(radius);
    if (CGAL::sign(exact_radius) != CGAL::POSITIVE) {
        throw AuditUnsupportedGeometryError(
            "arc radius must be exact positive");
    }
    require_canonical_world_xy_frame(xaxis, yaxis);
    const EPoint exact_center {
        Epeck::FT(center[0]), Epeck::FT(center[1])
    };
    const EVector zero_phase = exact_circle_phase(xaxis, radius);
    const Epeck::FT exact_cut_z(cut_z);
    const Epeck::FT exact_clearance_z(clearance_z);
    const CurvePlane plane = classify_curve_plane(
        Epeck::FT(center[2]), exact_cut_z,
        exact_clearance_z, operation_role);
    const Epeck::FT& active_z = plane == CurvePlane::Cut
        ? exact_cut_z
        : exact_clearance_z;
    AuditArcMotion2 motion = [&]() {
        try {
            return AuditArcMotion2::build(
                exact_center,
                zero_phase,
                exact_radius,
                start_angle,
                end_angle,
                clockwise,
                active_z);
        } catch (const AuditArcOrientationError& error) {
            throw AuditContradictoryOrientationError(
                std::string(
                    "operation orientation contradicts exact arc sweep sign: ")
                + error.what());
        } catch (const AuditArcMotionError& error) {
            throw AuditUnsupportedGeometryError(error.what());
        }
    }();
    if (plane == CurvePlane::Cut) {
        return motion;
    }
    return AuditClassificationFactory2::clearance_arc(
        motion,
        exact_cut_z,
        exact_clearance_z);
}
