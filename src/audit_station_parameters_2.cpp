#include "audit_certification_2.h"

#include "canonical_encoding.h"
#include "exact_circle_chart_2.h"

#include <utility>

namespace {

std::string point_bytes(const EPoint& point)
{
    return canonical_encode_tagged_union(
        "audit-exact-point2-v1",
        canonical_encode_component_map({
            {"x", canonical_audit_rational_bytes(point.x())},
            {"y", canonical_audit_rational_bytes(point.y())},
        }));
}

std::string chart_bytes(int chart)
{
    return canonical_audit_rational_bytes(Epeck::FT(chart));
}

std::string ordinal_bytes(std::size_t ordinal)
{
    return canonical_audit_rational_bytes(Epeck::FT(ordinal));
}

bool parameter_in_interval(
    const ExactArcChartInterval2& interval,
    const Epeck::FT& parameter)
{
    if (interval.increasing()) {
        return CGAL::compare(parameter, interval.start_parameter())
                   != CGAL::SMALLER
            && CGAL::compare(parameter, interval.end_parameter())
                   != CGAL::LARGER;
    }
    return CGAL::compare(parameter, interval.end_parameter())
               != CGAL::SMALLER
        && CGAL::compare(parameter, interval.start_parameter())
               != CGAL::LARGER;
}

} // namespace

AuditSegmentStationParameter2 AuditSegmentStationParameter2::build(
    const AuditSegmentMotion2& motion,
    const Epeck::FT& parameter)
{
    if (CGAL::sign(parameter) == CGAL::NEGATIVE
        || CGAL::compare(parameter, Epeck::FT(1)) == CGAL::LARGER) {
        throw AuditStationOutsideMotionError(
            "segment station parameter must lie in [0, 1]");
    }
    const EVector displacement = motion.xy().end - motion.xy().start;
    const EPoint point = motion.xy().start + displacement * parameter;
    std::string canonical = canonical_encode_tagged_union(
        "audit-segment-station-parameter-v1",
        canonical_encode_component_map({
            {"motion-digest", motion.digest().bytes()},
            {"parameter", canonical_audit_rational_bytes(parameter)},
            {"point", point_bytes(point)},
        }));
    return AuditSegmentStationParameter2(
        motion.digest(),
        point,
        parameter,
        std::move(canonical));
}

AuditSegmentStationParameter2::AuditSegmentStationParameter2(
    NativeMotionDigest2 motion_digest,
    EPoint point,
    Epeck::FT parameter,
    std::string canonical_bytes)
    : motion_digest_(std::move(motion_digest)),
      point_(std::move(point)),
      parameter_(std::move(parameter)),
      canonical_bytes_(std::move(canonical_bytes))
{
}

const NativeMotionDigest2&
AuditSegmentStationParameter2::motion_digest() const noexcept
{
    return motion_digest_;
}

const EPoint& AuditSegmentStationParameter2::point() const noexcept
{
    return point_;
}

const Epeck::FT&
AuditSegmentStationParameter2::parameter() const noexcept
{
    return parameter_;
}

const std::string&
AuditSegmentStationParameter2::canonical_bytes() const noexcept
{
    return canonical_bytes_;
}

AuditCircleStationParameter2 AuditCircleStationParameter2::build(
    const AuditCircleMotion2& motion,
    int chart,
    const Epeck::FT& parameter)
{
    if (chart < 0 || chart > 3
        || CGAL::sign(parameter) == CGAL::NEGATIVE
        || CGAL::compare(parameter, Epeck::FT(1)) == CGAL::LARGER) {
        throw AuditStationOutsideMotionError(
            "circle station parameter is outside the frozen atlas");
    }
    if (parameter == Epeck::FT(1)) {
        throw AuditStationSeamOwnershipError(
            "circle cardinal seams are owned only by chart t=0");
    }
    const EPoint point = exact_circle_chart_point(
        motion.xy().center,
        motion.xy().phase_vector,
        ExactCircleChartParameter2::build(chart, parameter));
    std::string canonical = canonical_encode_tagged_union(
        "audit-circle-station-parameter-v1",
        canonical_encode_component_map({
            {"chart", chart_bytes(chart)},
            {"motion-digest", motion.digest().bytes()},
            {"parameter", canonical_audit_rational_bytes(parameter)},
            {"point", point_bytes(point)},
        }));
    return AuditCircleStationParameter2(
        motion.digest(),
        point,
        std::move(canonical));
}

AuditCircleStationParameter2::AuditCircleStationParameter2(
    NativeMotionDigest2 motion_digest,
    EPoint point,
    std::string canonical_bytes)
    : motion_digest_(std::move(motion_digest)),
      point_(std::move(point)),
      canonical_bytes_(std::move(canonical_bytes))
{
}

const NativeMotionDigest2&
AuditCircleStationParameter2::motion_digest() const noexcept
{
    return motion_digest_;
}

const EPoint& AuditCircleStationParameter2::point() const noexcept
{
    return point_;
}

const std::string&
AuditCircleStationParameter2::canonical_bytes() const noexcept
{
    return canonical_bytes_;
}

AuditArcStationParameter2 AuditArcStationParameter2::from_interval(
    const AuditArcMotion2& motion,
    std::size_t interval_ordinal,
    int chart,
    const Epeck::FT& parameter)
{
    if (interval_ordinal >= motion.intervals().size()) {
        throw AuditStationOutsideMotionError(
            "arc station interval ordinal is outside the motion");
    }
    const ExactArcChartInterval2& interval =
        motion.intervals()[interval_ordinal];
    if (interval.chart() != chart
        || !parameter_in_interval(interval, parameter)) {
        throw AuditStationOutsideMotionError(
            "arc station parameter is outside its declared interval");
    }
    if (parameter == Epeck::FT(1)) {
        throw AuditStationSeamOwnershipError(
            "arc cardinal seams are owned only by chart t=0");
    }
    if (parameter == Epeck::FT(0)) {
        const bool owns_zero =
            (interval.start_parameter() == Epeck::FT(0)
                && interval.owns_start_seam())
            || (interval.end_parameter() == Epeck::FT(0)
                && interval.owns_end_seam());
        if (!owns_zero) {
            throw AuditStationSeamOwnershipError(
                "arc t=0 seam is not owned by the selected interval");
        }
    }
    const EPoint point = exact_circle_chart_point(
        motion.center(),
        motion.zero_phase(),
        ExactCircleChartParameter2::build(chart, parameter));
    std::string canonical = canonical_encode_tagged_union(
        "audit-arc-station-parameter-v1",
        canonical_encode_component_map({
            {"chart", chart_bytes(chart)},
            {"interval-ordinal", ordinal_bytes(interval_ordinal)},
            {"motion-digest", motion.digest().bytes()},
            {"parameter", canonical_audit_rational_bytes(parameter)},
            {"point", point_bytes(point)},
        }));
    return AuditArcStationParameter2(
        motion.digest(),
        point,
        std::move(canonical));
}

AuditArcStationParameter2 AuditArcStationParameter2::from_terminal_anchor(
    const AuditArcMotion2& motion)
{
    std::string canonical = canonical_encode_tagged_union(
        "audit-arc-terminal-anchor-v1",
        canonical_encode_component_map({
            {"chart", chart_bytes(motion.end_parameter().chart())},
            {"motion-digest", motion.digest().bytes()},
            {"parameter", canonical_audit_rational_bytes(
                 motion.end_parameter().parameter())},
            {"point", point_bytes(motion.end_point())},
        }));
    return AuditArcStationParameter2(
        motion.digest(),
        motion.end_point(),
        std::move(canonical));
}

AuditArcStationParameter2 AuditArcStationParameter2::from_start_anchor(
    const AuditArcMotion2& motion)
{
    std::string canonical = canonical_encode_tagged_union(
        "audit-arc-start-anchor-v1",
        canonical_encode_component_map({
            {"chart", chart_bytes(motion.start_parameter().chart())},
            {"motion-digest", motion.digest().bytes()},
            {"parameter", canonical_audit_rational_bytes(
                 motion.start_parameter().parameter())},
            {"point", point_bytes(motion.start_point())},
        }));
    return AuditArcStationParameter2(
        motion.digest(),
        motion.start_point(),
        std::move(canonical));
}

AuditArcStationParameter2::AuditArcStationParameter2(
    NativeMotionDigest2 motion_digest,
    EPoint point,
    std::string canonical_bytes)
    : motion_digest_(std::move(motion_digest)),
      point_(std::move(point)),
      canonical_bytes_(std::move(canonical_bytes))
{
}

const NativeMotionDigest2&
AuditArcStationParameter2::motion_digest() const noexcept
{
    return motion_digest_;
}

const EPoint& AuditArcStationParameter2::point() const noexcept
{
    return point_;
}

const std::string&
AuditArcStationParameter2::canonical_bytes() const noexcept
{
    return canonical_bytes_;
}
