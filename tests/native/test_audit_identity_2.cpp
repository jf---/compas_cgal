#include "audit_classification_2.h"
#include "canonical_encoding.h"
#include "audit_digest_2.h"
#include "audit_motion_identity_2.h"
#include "audit_policy_2.h"
#include "audit_request_identity_2.h"
#include "continuous_tea_2/sha256.h"

#include <array>
#include <numbers>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace {

template <class Digest>
concept HasGenericFromBytes = requires(std::string bytes) {
    Digest::from_bytes(std::move(bytes));
};

template <class Authority>
concept HasPublicCanonicalHash = requires(std::string bytes) {
    Authority::hash_canonical(bytes);
};

static_assert(!HasGenericFromBytes<NativeMotionDigest2>);
static_assert(!HasGenericFromBytes<AuditPolicyDigest2>);
static_assert(!HasGenericFromBytes<AuditNativeRequestDigest2>);
static_assert(!std::is_constructible_v<NativeMotionDigest2, std::string>);
static_assert(!HasPublicCanonicalHash<AuditNativeStockDigestAuthority2>);
static_assert(!HasPublicCanonicalHash<AuditNativeRequestDigestAuthority2>);
static_assert(!HasPublicCanonicalHash<AuditPolicyDigestAuthority2>);
static_assert(!HasPublicCanonicalHash<NativeMotionDigestAuthority2>);
static_assert(!HasPublicCanonicalHash<NativeDecisionDigestAuthority2>);
static_assert(!HasPublicCanonicalHash<DepletionWitnessDigestAuthority2>);
static_assert(!HasPublicCanonicalHash<AuditStockStateDigestAuthority2>);
static_assert(!HasPublicCanonicalHash<StockLineageDigestAuthority2>);
static_assert(!HasPublicCanonicalHash<AuditResultDigestAuthority2>);
static_assert(!std::is_default_constructible_v<AuditSegmentMotion2>);
static_assert(!std::is_default_constructible_v<AuditCircleMotion2>);
static_assert(!std::is_default_constructible_v<AuditVerticalPlunge2>);
static_assert(!std::is_default_constructible_v<AuditVerticalRetract2>);
static_assert(!std::is_default_constructible_v<AuditClearanceTransport2>);

void require(bool condition, const char* message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

std::string expected_exact_point(const EPoint& point)
{
    return canonical_encode_tagged_union(
        "audit-exact-point2-v1",
        canonical_encode_sequence({
            canonical_audit_rational_bytes(point.x()),
            canonical_audit_rational_bytes(point.y()),
        }));
}

std::string expected_exact_vector(const EVector& vector)
{
    return canonical_encode_tagged_union(
        "audit-exact-vector2-v1",
        canonical_encode_sequence({
            canonical_audit_rational_bytes(vector.x()),
            canonical_audit_rational_bytes(vector.y()),
        }));
}

std::vector<std::pair<std::string, std::string>> expected_planes(
    const Epeck::FT& cut_z,
    const Epeck::FT& clearance_z)
{
    return {
        {"clearance-z", canonical_audit_rational_bytes(clearance_z)},
        {"cut-z", canonical_audit_rational_bytes(cut_z)},
    };
}

std::string expected_motion(
    const std::string& tag,
    std::vector<std::pair<std::string, std::string>> fields)
{
    return canonical_encode_tagged_union(
        tag,
        canonical_encode_component_map(fields));
}

std::string expected_segment(
    const ExactSegmentMotion2& motion,
    const Epeck::FT& cut_z,
    const Epeck::FT& clearance_z)
{
    auto fields = expected_planes(cut_z, clearance_z);
    fields.emplace_back("end", expected_exact_point(motion.end));
    fields.emplace_back("start", expected_exact_point(motion.start));
    return expected_motion("audit-segment-motion-v1", std::move(fields));
}

std::string expected_circle(
    const std::string& tag,
    const ExactCircleMotion2& motion,
    const Epeck::FT& guide_radius,
    const Epeck::FT& cut_z,
    const Epeck::FT& clearance_z)
{
    auto fields = expected_planes(cut_z, clearance_z);
    fields.emplace_back("center", expected_exact_point(motion.center));
    fields.emplace_back(
        "clockwise", canonical_encode_boolean(motion.clockwise));
    fields.emplace_back(
        "guide-radius", canonical_audit_rational_bytes(guide_radius));
    fields.emplace_back(
        "phase-vector", expected_exact_vector(motion.phase_vector));
    return expected_motion(tag, std::move(fields));
}

std::string expected_vertical(
    const std::string& tag,
    const EPoint& cut_endpoint,
    const Epeck::FT& cut_z,
    const Epeck::FT& clearance_z)
{
    auto fields = expected_planes(cut_z, clearance_z);
    fields.emplace_back("cut-endpoint", expected_exact_point(cut_endpoint));
    return expected_motion(tag, std::move(fields));
}

std::string expected_clearance_line(
    const EPoint& start,
    const EPoint& end,
    const Epeck::FT& cut_z,
    const Epeck::FT& clearance_z)
{
    auto fields = expected_planes(cut_z, clearance_z);
    fields.emplace_back("end", expected_exact_point(end));
    fields.emplace_back("start", expected_exact_point(start));
    return expected_motion("audit-clearance-line-v1", std::move(fields));
}

std::string expected_quarter_arc_canonical(
    const EPoint& center,
    const EVector& zero_phase,
    const Epeck::FT& guide_radius,
    const Epeck::FT& cut_z)
{
    std::string canonical("audit-arc-motion-v1");
    append_audit_bytes(
        canonical, "audit-arc-quarter-chart-binary64-v2");
    append_audit_bytes(
        canonical, "exact-quarter-pythagorean-chart-v1");
    append_audit_bytes(
        canonical, canonical_audit_rational_bytes(center.x()));
    append_audit_bytes(
        canonical, canonical_audit_rational_bytes(center.y()));
    append_audit_bytes(
        canonical, canonical_audit_rational_bytes(zero_phase.x()));
    append_audit_bytes(
        canonical, canonical_audit_rational_bytes(zero_phase.y()));
    append_audit_bytes(
        canonical, canonical_audit_rational_bytes(guide_radius));
    append_audit_bytes(canonical, canonical_audit_binary64_bytes(0.0));
    append_audit_bytes(
        canonical,
        canonical_audit_binary64_bytes(std::numbers::pi / 2.0));
    append_audit_bytes(
        canonical,
        canonical_audit_binary64_bytes(std::numbers::pi / 2.0));
    append_audit_bytes(canonical, std::string(1, '\0'));
    append_audit_bytes(
        canonical, canonical_audit_rational_bytes(cut_z));
    append_audit_bytes(canonical, std::string(1, '\0'));
    append_audit_bytes(
        canonical, canonical_audit_rational_bytes(Epeck::FT(0)));
    append_audit_bytes(
        canonical, canonical_audit_rational_bytes(Epeck::FT(1)));
    append_audit_bytes(canonical, std::string(1, '\1'));
    append_audit_bytes(canonical, std::string(1, '\1'));
    append_audit_bytes(canonical, std::string(1, '\0'));
    return canonical;
}

std::string expected_clearance_arc(
    const std::string& arc_digest,
    const Epeck::FT& cut_z,
    const Epeck::FT& clearance_z)
{
    auto fields = expected_planes(cut_z, clearance_z);
    fields.emplace_back("arc-motion-digest", arc_digest);
    return expected_motion("audit-clearance-arc-v1", std::move(fields));
}

void require_digest(
    const NativeMotionDigest2& observed,
    const std::string& independently_built_canonical,
    const char* message)
{
    require(
        observed.bytes() == sha256_bytes(independently_built_canonical),
        message);
}

void production_classification_identity_gate()
{
    const Epeck::FT cut_z(0);
    const Epeck::FT clearance_z(5);
    const EPoint line_start(1, 2);
    const EPoint line_end(4, 2);

    const auto segment_classification = classify_audit_line(
        {1.0, 2.0, 0.0}, {4.0, 2.0, 0.0}, 0.0, 5.0, "cut");
    const auto& segment = std::get<AuditSegmentMotion2>(
        segment_classification);
    require(
        segment.xy().start == line_start
            && segment.xy().end == line_end
            && segment.cut_z() == cut_z
            && segment.clearance_z() == clearance_z,
        "classified segment did not retain exact source geometry and planes");
    require_digest(
        segment.digest(),
        expected_segment(
            ExactSegmentMotion2 {line_start, line_end},
            cut_z,
            clearance_z),
        "classified segment digest diverges from independent canonical bytes");

    const auto plunge_classification = classify_audit_line(
        {1.0, 2.0, 5.0}, {1.0, 2.0, 0.0}, 0.0, 5.0, "plunge");
    const auto& plunge = std::get<AuditVerticalPlunge2>(
        plunge_classification);
    require(
        plunge.cut_endpoint() == line_start
            && plunge.cut_z() == cut_z
            && plunge.clearance_z() == clearance_z,
        "classified plunge did not retain exact endpoint and planes");
    require_digest(
        plunge.digest(),
        expected_vertical(
            "audit-vertical-plunge-v1", line_start, cut_z, clearance_z),
        "classified plunge digest diverges from independent canonical bytes");

    const auto retract_classification = classify_audit_line(
        {1.0, 2.0, 0.0}, {1.0, 2.0, 5.0}, 0.0, 5.0, "retract");
    const auto& retract = std::get<AuditVerticalRetract2>(
        retract_classification);
    require(
        retract.cut_endpoint() == line_start
            && retract.cut_z() == cut_z
            && retract.clearance_z() == clearance_z,
        "classified retract did not retain exact endpoint and planes");
    require_digest(
        retract.digest(),
        expected_vertical(
            "audit-vertical-retract-v1", line_start, cut_z, clearance_z),
        "classified retract digest diverges from independent canonical bytes");

    const auto clearance_line_classification = classify_audit_line(
        {1.0, 2.0, 5.0}, {4.0, 2.0, 5.0}, 0.0, 5.0, "link");
    const auto& clearance_line = std::get<AuditClearanceTransport2>(
        clearance_line_classification);
    require(
        clearance_line.start() == line_start
            && clearance_line.end() == line_end
            && clearance_line.cut_z() == cut_z
            && clearance_line.clearance_z() == clearance_z,
        "classified clearance line did not retain exact geometry and planes");
    require_digest(
        clearance_line.digest(),
        expected_clearance_line(
            line_start, line_end, cut_z, clearance_z),
        "classified clearance line digest diverges from independent canonical bytes");
    const auto shifted_clearance_line_classification = classify_audit_line(
        {1.0, 2.0, 7.0}, {4.0, 2.0, 7.0}, -2.0, 7.0, "link");
    const auto& shifted_clearance_line =
        std::get<AuditClearanceTransport2>(
            shifted_clearance_line_classification);
    require(
        clearance_line.digest().bytes()
            != shifted_clearance_line.digest().bytes(),
        "same-XY clearance line digest ignores shifted declared planes");

    const ExactCircleMotion2 circle_xy {
        EPoint(2, 3), EVector(2, 0), false
    };
    const auto circle_classification = classify_audit_circle(
        {2.0, 3.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        2.0,
        false,
        0.0,
        5.0,
        "cut");
    const auto& circle = std::get<AuditCircleMotion2>(
        circle_classification);
    require(
        circle.xy().center == circle_xy.center
            && circle.xy().phase_vector == circle_xy.phase_vector
            && circle.xy().clockwise == circle_xy.clockwise
            && circle.guide_radius() == Epeck::FT(2)
            && circle.cut_z() == cut_z
            && circle.clearance_z() == clearance_z,
        "classified circle did not retain exact source geometry and planes");
    require_digest(
        circle.digest(),
        expected_circle(
            "audit-circle-motion-v1",
            circle_xy,
            Epeck::FT(2),
            cut_z,
            clearance_z),
        "classified circle digest diverges from independent canonical bytes");

    const auto clearance_circle_classification = classify_audit_circle(
        {2.0, 3.0, 5.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        2.0,
        false,
        0.0,
        5.0,
        "link");
    const auto& clearance_circle = std::get<AuditClearanceTransport2>(
        clearance_circle_classification);
    require(
        clearance_circle.start() == EPoint(4, 3)
            && clearance_circle.end() == EPoint(4, 3)
            && clearance_circle.cut_z() == cut_z
            && clearance_circle.clearance_z() == clearance_z,
        "classified clearance circle did not retain exact seam and planes");
    require_digest(
        clearance_circle.digest(),
        expected_circle(
            "audit-clearance-circle-v1",
            circle_xy,
            Epeck::FT(2),
            cut_z,
            clearance_z),
        "classified clearance circle digest omits full source identity");

    const auto arc_classification = classify_audit_arc(
        {2.0, 3.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        2.0,
        0.0,
        std::numbers::pi / 2.0,
        false,
        0.0,
        5.0,
        "cut");
    const auto& arc = std::get<AuditArcMotion2>(arc_classification);
    require(
        arc.center() == EPoint(2, 3)
            && arc.zero_phase() == EVector(2, 0)
            && arc.guide_radius() == Epeck::FT(2)
            && arc.cut_z() == cut_z,
        "classified arc did not retain exact source geometry and plane");
    require_digest(
        arc.digest(),
        expected_quarter_arc_canonical(
            EPoint(2, 3), EVector(2, 0), Epeck::FT(2), cut_z),
        "classified arc digest diverges from independent canonical bytes");

    const auto clearance_arc_classification = classify_audit_arc(
        {2.0, 3.0, 5.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        2.0,
        0.0,
        std::numbers::pi / 2.0,
        false,
        0.0,
        5.0,
        "link");
    const auto& clearance_arc = std::get<AuditClearanceTransport2>(
        clearance_arc_classification);
    require(
        clearance_arc.start() == EPoint(4, 3)
            && clearance_arc.end() == EPoint(2, 5)
            && clearance_arc.cut_z() == cut_z
            && clearance_arc.clearance_z() == clearance_z,
        "classified clearance arc did not retain exact endpoints and planes");
    const std::string expected_clearance_source_digest = sha256_bytes(
        expected_quarter_arc_canonical(
            EPoint(2, 3),
            EVector(2, 0),
            Epeck::FT(2),
            clearance_z));
    require_digest(
        clearance_arc.digest(),
        expected_clearance_arc(
            expected_clearance_source_digest, cut_z, clearance_z),
        "classified clearance arc digest omits full source identity");
}

} // namespace

void audit_identity_gate()
{
    production_classification_identity_gate();
    const std::string segment = canonical_audit_segment_motion_bytes(
        ExactSegmentMotion2 {EPoint(0, 0), EPoint(1, 0)},
        Epeck::FT(0),
        Epeck::FT(5));
    const std::string changed = canonical_audit_segment_motion_bytes(
        ExactSegmentMotion2 {EPoint(0, 0), EPoint(2, 0)},
        Epeck::FT(0),
        Epeck::FT(5));
    require(
        segment != changed,
        "segment digest ignores endpoint geometry");

    const double cap_radians = std::numbers::pi / 2.0;
    const AuditPolicy2 policy = AuditPolicy2::build(
        AuditCapObservation2::build(
            cap_radians,
            audit_cap_chord_ratio(cap_radians)),
        Epeck::FT(2),
        Epeck::FT(1) / Epeck::FT(16),
        4096);
    require(
        policy.digest().bytes().size() == 32,
        "policy digest is not SHA-256 sized");

    compas::RowMatrixXd boundary(4, 2);
    boundary << 0.0, 0.0,
        10.0, 0.0,
        10.0, 10.0,
        0.0, 10.0;
    const AuditNativeStockIdentity2 stock =
        AuditNativeStockIdentity2::build(boundary, {});
    require(
        stock.digest().bytes().size() == 32,
        "stock digest is not SHA-256 sized");

    const AuditArcMotion2 first = AuditArcMotion2::build(
        EPoint(2, 2),
        EVector(1, 0),
        Epeck::FT(1),
        0.0,
        std::numbers::pi / 2.0,
        false,
        Epeck::FT(0));
    const AuditArcMotion2 second = AuditArcMotion2::build(
        EPoint(4, 4),
        EVector(1, 0),
        Epeck::FT(1),
        0.0,
        std::numbers::pi / 2.0,
        false,
        Epeck::FT(0));
    const AuditNativeRequestIdentity2 request =
        AuditNativeRequestIdentity2::build(
            stock,
            policy,
            {first.digest(), second.digest()});
    const AuditNativeRequestIdentity2 reordered =
        AuditNativeRequestIdentity2::build(
            stock,
            policy,
            {second.digest(), first.digest()});
    require(
        request.digest().bytes().size() == 32,
        "request digest is not SHA-256 sized");
    require(
        request.digest().bytes() != reordered.digest().bytes(),
        "request digest ignores native motion order");
}
