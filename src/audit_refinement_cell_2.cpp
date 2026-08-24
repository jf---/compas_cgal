#include "audit_certification_internal_2.h"

#include "canonical_encoding.h"
#include "exact_circle_chart_2.h"

#include <utility>

namespace {
std::string refinement_cell_bytes(
    const AuditArcMotion2& motion,
    std::size_t interval_ordinal,
    const Epeck::FT& start_parameter,
    const Epeck::FT& end_parameter,
    std::size_t depth)
{
    return canonical_encode_tagged_union(
        "audit-arc-refinement-cell-v1",
        canonical_encode_component_map({
            {"depth", canonical_audit_rational_bytes(Epeck::FT(depth))},
            {"end", canonical_audit_rational_bytes(end_parameter)},
            {"interval", canonical_audit_rational_bytes(
                 Epeck::FT(interval_ordinal))},
            {"motion-digest", motion.digest().bytes()},
            {"start", canonical_audit_rational_bytes(start_parameter)},
        }));
}

} // namespace

AuditArcRefinementCell2 AuditArcRefinementCell2::root(
    const AuditArcMotion2& motion,
    std::size_t interval_ordinal)
{
    if (interval_ordinal >= motion.intervals().size()) {
        throw AuditStationOutsideMotionError(
            "arc refinement cell interval is outside the motion");
    }
    const ExactArcChartInterval2& interval =
        motion.intervals()[interval_ordinal];
    if (interval.start_parameter() == interval.end_parameter()) {
        throw AuditStationOutsideMotionError(
            "arc refinement root must have nonzero traversal extent");
    }
    return AuditArcRefinementCell2(
        motion.digest(),
        interval_ordinal,
        interval.start_parameter(),
        interval.end_parameter(),
        0,
        refinement_cell_bytes(
            motion,
            interval_ordinal,
            interval.start_parameter(),
            interval.end_parameter(),
            0));
}

AuditArcRefinementCell2 AuditArcRefinementCell2::first_child(
    const AuditArcMotion2& motion,
    const AuditArcRefinementCell2& parent)
{
    if (parent.motion_digest().bytes() != motion.digest().bytes()) {
        throw AuditDecisionEvidenceReplayError(
            "arc refinement parent belongs to a foreign motion");
    }
    const Epeck::FT midpoint =
        (parent.start_parameter() + parent.end_parameter()) / Epeck::FT(2);
    if (midpoint == parent.start_parameter()
        || midpoint == parent.end_parameter()) {
        throw AuditStationOutsideMotionError(
            "arc refinement child has zero traversal extent");
    }
    const std::size_t depth = parent.depth() + 1;
    return AuditArcRefinementCell2(
        motion.digest(),
        parent.interval_ordinal(),
        parent.start_parameter(),
        midpoint,
        depth,
        refinement_cell_bytes(
            motion,
            parent.interval_ordinal(),
            parent.start_parameter(),
            midpoint,
            depth));
}

AuditArcRefinementCell2 AuditArcRefinementCell2::second_child(
    const AuditArcMotion2& motion,
    const AuditArcRefinementCell2& parent)
{
    if (parent.motion_digest().bytes() != motion.digest().bytes()) {
        throw AuditDecisionEvidenceReplayError(
            "arc refinement parent belongs to a foreign motion");
    }
    const Epeck::FT midpoint =
        (parent.start_parameter() + parent.end_parameter()) / Epeck::FT(2);
    if (midpoint == parent.start_parameter()
        || midpoint == parent.end_parameter()) {
        throw AuditStationOutsideMotionError(
            "arc refinement child has zero traversal extent");
    }
    const std::size_t depth = parent.depth() + 1;
    return AuditArcRefinementCell2(
        motion.digest(),
        parent.interval_ordinal(),
        midpoint,
        parent.end_parameter(),
        depth,
        refinement_cell_bytes(
            motion,
            parent.interval_ordinal(),
            midpoint,
            parent.end_parameter(),
            depth));
}

AuditArcRefinementCell2::AuditArcRefinementCell2(
    NativeMotionDigest2 motion_digest,
    std::size_t interval_ordinal,
    Epeck::FT start_parameter,
    Epeck::FT end_parameter,
    std::size_t depth,
    std::string canonical_bytes)
    : motion_digest_(std::move(motion_digest)),
      interval_ordinal_(interval_ordinal),
      start_parameter_(std::move(start_parameter)),
      end_parameter_(std::move(end_parameter)),
      depth_(depth),
      canonical_bytes_(std::move(canonical_bytes))
{
}

const NativeMotionDigest2&
AuditArcRefinementCell2::motion_digest() const noexcept
{
    return motion_digest_;
}

std::size_t AuditArcRefinementCell2::interval_ordinal() const noexcept
{
    return interval_ordinal_;
}

const Epeck::FT&
AuditArcRefinementCell2::start_parameter() const noexcept
{
    return start_parameter_;
}

const Epeck::FT& AuditArcRefinementCell2::end_parameter() const noexcept
{
    return end_parameter_;
}

std::size_t AuditArcRefinementCell2::depth() const noexcept
{
    return depth_;
}

Epeck::FT AuditArcRefinementCell2::squared_chord_length(
    const AuditArcMotion2& motion) const
{
    if (motion.digest().bytes() != motion_digest_.bytes()) {
        throw AuditDecisionEvidenceReplayError(
            "arc refinement cell belongs to a foreign motion");
    }
    const int chart = motion.intervals()[interval_ordinal_].chart();
    const EPoint start = exact_circle_chart_point(
        motion.center(),
        motion.zero_phase(),
        ExactCircleChartParameter2::build(chart, start_parameter_));
    const EPoint end = exact_circle_chart_point(
        motion.center(),
        motion.zero_phase(),
        ExactCircleChartParameter2::build(chart, end_parameter_));
    return CGAL::squared_distance(start, end);
}

const std::string& AuditArcRefinementCell2::canonical_bytes() const noexcept
{
    return canonical_bytes_;
}
