#include "audit_motion_identity_2.h"

#include "canonical_encoding.h"

#include <string>
#include <utility>
#include <vector>

AuditSegmentMotion2::AuditSegmentMotion2(
    ExactSegmentMotion2 xy,
    Epeck::FT cut_z,
    Epeck::FT clearance_z,
    NativeMotionDigest2 digest)
    : xy_(std::move(xy)),
      cut_z_(std::move(cut_z)),
      clearance_z_(std::move(clearance_z)),
      digest_(std::move(digest))
{
}

const ExactSegmentMotion2& AuditSegmentMotion2::xy() const noexcept { return xy_; }
const Epeck::FT& AuditSegmentMotion2::cut_z() const noexcept { return cut_z_; }
const Epeck::FT& AuditSegmentMotion2::clearance_z() const noexcept { return clearance_z_; }
const NativeMotionDigest2& AuditSegmentMotion2::digest() const noexcept { return digest_; }

AuditCircleMotion2::AuditCircleMotion2(
    ExactCircleMotion2 xy,
    Epeck::FT guide_radius,
    Epeck::FT cut_z,
    Epeck::FT clearance_z,
    NativeMotionDigest2 digest)
    : xy_(std::move(xy)),
      guide_radius_(std::move(guide_radius)),
      cut_z_(std::move(cut_z)),
      clearance_z_(std::move(clearance_z)),
      digest_(std::move(digest))
{
}

const ExactCircleMotion2& AuditCircleMotion2::xy() const noexcept { return xy_; }
const Epeck::FT& AuditCircleMotion2::guide_radius() const noexcept { return guide_radius_; }
const Epeck::FT& AuditCircleMotion2::cut_z() const noexcept { return cut_z_; }
const Epeck::FT& AuditCircleMotion2::clearance_z() const noexcept { return clearance_z_; }
const NativeMotionDigest2& AuditCircleMotion2::digest() const noexcept { return digest_; }

AuditVerticalPlunge2::AuditVerticalPlunge2(
    EPoint cut_endpoint,
    Epeck::FT cut_z,
    Epeck::FT clearance_z,
    NativeMotionDigest2 digest)
    : cut_endpoint_(std::move(cut_endpoint)),
      cut_z_(std::move(cut_z)),
      clearance_z_(std::move(clearance_z)),
      digest_(std::move(digest))
{
}

const EPoint& AuditVerticalPlunge2::cut_endpoint() const noexcept { return cut_endpoint_; }
const Epeck::FT& AuditVerticalPlunge2::cut_z() const noexcept { return cut_z_; }
const Epeck::FT& AuditVerticalPlunge2::clearance_z() const noexcept { return clearance_z_; }
const NativeMotionDigest2& AuditVerticalPlunge2::digest() const noexcept { return digest_; }

AuditVerticalRetract2::AuditVerticalRetract2(
    EPoint cut_endpoint,
    Epeck::FT cut_z,
    Epeck::FT clearance_z,
    NativeMotionDigest2 digest)
    : cut_endpoint_(std::move(cut_endpoint)),
      cut_z_(std::move(cut_z)),
      clearance_z_(std::move(clearance_z)),
      digest_(std::move(digest))
{
}

const EPoint& AuditVerticalRetract2::cut_endpoint() const noexcept { return cut_endpoint_; }
const Epeck::FT& AuditVerticalRetract2::cut_z() const noexcept { return cut_z_; }
const Epeck::FT& AuditVerticalRetract2::clearance_z() const noexcept { return clearance_z_; }
const NativeMotionDigest2& AuditVerticalRetract2::digest() const noexcept { return digest_; }

AuditClearanceTransport2::AuditClearanceTransport2(
    EPoint start,
    EPoint end,
    Epeck::FT cut_z,
    Epeck::FT clearance_z,
    NativeMotionDigest2 digest)
    : start_(std::move(start)),
      end_(std::move(end)),
      cut_z_(std::move(cut_z)),
      clearance_z_(std::move(clearance_z)),
      digest_(std::move(digest))
{
}

const EPoint& AuditClearanceTransport2::start() const noexcept { return start_; }
const EPoint& AuditClearanceTransport2::end() const noexcept { return end_; }
const Epeck::FT& AuditClearanceTransport2::cut_z() const noexcept { return cut_z_; }
const Epeck::FT& AuditClearanceTransport2::clearance_z() const noexcept { return clearance_z_; }
const NativeMotionDigest2& AuditClearanceTransport2::digest() const noexcept { return digest_; }

namespace {

std::string exact_point_bytes(const EPoint& point)
{
    return canonical_encode_tagged_union(
        "audit-exact-point2-v1",
        canonical_encode_sequence({
            canonical_audit_rational_bytes(point.x()),
            canonical_audit_rational_bytes(point.y()),
        }));
}

std::string exact_vector_bytes(const EVector& vector)
{
    return canonical_encode_tagged_union(
        "audit-exact-vector2-v1",
        canonical_encode_sequence({
            canonical_audit_rational_bytes(vector.x()),
            canonical_audit_rational_bytes(vector.y()),
        }));
}

std::string canonical_motion(
    const std::string& tag,
    std::vector<std::pair<std::string, std::string>> fields)
{
    return canonical_encode_tagged_union(
        tag,
        canonical_encode_component_map(fields));
}

std::vector<std::pair<std::string, std::string>> plane_fields(
    const Epeck::FT& cut_z,
    const Epeck::FT& clearance_z)
{
    return {
        {"clearance-z", canonical_audit_rational_bytes(clearance_z)},
        {"cut-z", canonical_audit_rational_bytes(cut_z)},
    };
}

void append_circle_fields(
    std::vector<std::pair<std::string, std::string>>& fields,
    const ExactCircleMotion2& motion,
    const Epeck::FT& guide_radius)
{
    fields.emplace_back("center", exact_point_bytes(motion.center));
    fields.emplace_back("clockwise", canonical_encode_boolean(motion.clockwise));
    fields.emplace_back("guide-radius", canonical_audit_rational_bytes(guide_radius));
    fields.emplace_back("phase-vector", exact_vector_bytes(motion.phase_vector));
}

} // namespace

std::string canonical_audit_segment_motion_bytes(
    const ExactSegmentMotion2& motion,
    const Epeck::FT& cut_z,
    const Epeck::FT& clearance_z)
{
    auto fields = plane_fields(cut_z, clearance_z);
    fields.emplace_back("end", exact_point_bytes(motion.end));
    fields.emplace_back("start", exact_point_bytes(motion.start));
    return canonical_motion("audit-segment-motion-v1", std::move(fields));
}

std::string canonical_audit_circle_motion_bytes(
    const ExactCircleMotion2& motion,
    const Epeck::FT& guide_radius,
    const Epeck::FT& cut_z,
    const Epeck::FT& clearance_z)
{
    auto fields = plane_fields(cut_z, clearance_z);
    append_circle_fields(fields, motion, guide_radius);
    return canonical_motion("audit-circle-motion-v1", std::move(fields));
}

std::string canonical_audit_vertical_plunge_bytes(
    const EPoint& cut_endpoint,
    const Epeck::FT& cut_z,
    const Epeck::FT& clearance_z)
{
    auto fields = plane_fields(cut_z, clearance_z);
    fields.emplace_back("cut-endpoint", exact_point_bytes(cut_endpoint));
    return canonical_motion("audit-vertical-plunge-v1", std::move(fields));
}

std::string canonical_audit_vertical_retract_bytes(
    const EPoint& cut_endpoint,
    const Epeck::FT& cut_z,
    const Epeck::FT& clearance_z)
{
    auto fields = plane_fields(cut_z, clearance_z);
    fields.emplace_back("cut-endpoint", exact_point_bytes(cut_endpoint));
    return canonical_motion("audit-vertical-retract-v1", std::move(fields));
}

std::string canonical_audit_clearance_line_bytes(
    const EPoint& start,
    const EPoint& end,
    const Epeck::FT& cut_z,
    const Epeck::FT& clearance_z)
{
    auto fields = plane_fields(cut_z, clearance_z);
    fields.emplace_back("end", exact_point_bytes(end));
    fields.emplace_back("start", exact_point_bytes(start));
    return canonical_motion("audit-clearance-line-v1", std::move(fields));
}

std::string canonical_audit_clearance_circle_bytes(
    const ExactCircleMotion2& motion,
    const Epeck::FT& guide_radius,
    const Epeck::FT& cut_z,
    const Epeck::FT& clearance_z)
{
    auto fields = plane_fields(cut_z, clearance_z);
    append_circle_fields(fields, motion, guide_radius);
    return canonical_motion("audit-clearance-circle-v1", std::move(fields));
}

std::string canonical_audit_clearance_arc_bytes(
    const AuditArcMotion2& motion,
    const Epeck::FT& cut_z,
    const Epeck::FT& clearance_z)
{
    auto fields = plane_fields(cut_z, clearance_z);
    fields.emplace_back("arc-motion-digest", motion.digest().bytes());
    return canonical_motion("audit-clearance-arc-v1", std::move(fields));
}
