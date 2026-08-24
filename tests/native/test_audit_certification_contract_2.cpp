#include "audit_certification_fixtures_2.h"

#include "canonical_encoding.h"
#include "continuous_tea_2/event_certificate.h"
#include "continuous_tea_2/circle_oracle.h"
#include "continuous_tea_2/sha256.h"
#include "continuous_tea_2/segment_source.h"
#include "engagement_2.h"

#include <cstddef>
#include <initializer_list>
#include <limits>
#include <numbers>
#include <string_view>
#include <type_traits>

namespace {

using namespace audit_certification_fixtures;

std::string bytes_from_hex(std::string_view hex)
{
    const auto nibble = [](char value) -> unsigned char {
        if (value >= '0' && value <= '9') {
            return static_cast<unsigned char>(value - '0');
        }
        return static_cast<unsigned char>(value - 'a' + 10);
    };
    std::string bytes;
    bytes.reserve(hex.size() / 2);
    for (std::size_t index = 0; index < hex.size(); index += 2) {
        bytes.push_back(static_cast<char>(
            (nibble(hex[index]) << 4) | nibble(hex[index + 1])));
    }
    return bytes;
}

template <class Motion>
concept AuditCertifiableLateral = requires(
    const Stock2& stock,
    const Motion& motion,
    const AuditPolicy2& audit_policy,
    const AuditDecisionLimits2& decision_limits) {
    certify_audit_tea_exact(stock, motion, audit_policy, decision_limits);
};

template <class Decision>
concept DecisionCarriesReportingDouble = requires(const Decision& decision) {
    decision.reported_max_tea();
};

using ExactStationReplaySignature = AuditExactStationDisposition2 (*)(
    const Stock2&,
    const EPoint&,
    const Epeck::FT&,
    const Epeck::FT&);

static_assert(AuditCertifiableLateral<AuditSegmentMotion2>);
static_assert(AuditCertifiableLateral<AuditCircleMotion2>);
static_assert(AuditCertifiableLateral<AuditArcMotion2>);
static_assert(!AuditCertifiableLateral<AuditVerticalPlunge2>);
static_assert(!AuditCertifiableLateral<AuditVerticalRetract2>);
static_assert(!AuditCertifiableLateral<AuditClearanceTransport2>);
static_assert(!DecisionCarriesReportingDouble<AuditDecisionWitness2>);
static_assert(!std::is_default_constructible_v<AuditDecisionWitness2>);
static_assert(!std::is_default_constructible_v<AuditCertifiedCoverage2>);
static_assert(std::is_same_v<
    decltype(&replay_audit_unguarded_station_exact),
    ExactStationReplaySignature>);
static_assert(requires(
    const ExactSegmentMotion2& motion,
    const Epeck::FT& tool_radius,
    const Epeck::FT& cap_ratio) {
    SegmentEventSource2::from_exact(motion, tool_radius, cap_ratio);
});
static_assert(requires(
    const ExactCircleMotion2& motion,
    const Epeck::FT& tool_radius,
    const Epeck::FT& cap_ratio) {
    FullCircleEventSource2::from_exact(motion, tool_radius, cap_ratio);
});
static_assert(requires(
    const AuditSegmentMotion2& motion,
    const Epeck::FT& parameter) {
    AuditSegmentStationParameter2::build(motion, parameter);
});
static_assert(requires(const AuditArcMotion2& motion) {
    AuditArcStationParameter2::from_start_anchor(motion);
    AuditArcStationParameter2::from_terminal_anchor(motion);
});

void exact_ingress_and_limit_errors_gate()
{
    // Production mutation caught: binary64-only event ingress or unsealed
    // limits bypass exact decision authority and guarded exhaustion.
    const AuditPolicy2 audit_policy = policy();
    const AuditSegmentMotion2 motion = short_segment();
    const SegmentEventSource2 segment_source = SegmentEventSource2::from_exact(
        motion.xy(),
        audit_policy.tool_radius_mm(),
        audit_policy.engagement_cap().chord_ratio());
    require(
        !segment_source.canonical_bytes().empty(),
        "segment exact ingress did not produce a source identity");

    const AuditCircleMotion2 circle = full_circle(false);
    const FullCircleEventSource2 circle_source =
        FullCircleEventSource2::from_exact(
            circle.xy(),
            audit_policy.tool_radius_mm(),
            audit_policy.engagement_cap().chord_ratio());
    require(
        !circle_source.canonical_bytes().empty(),
        "full-circle exact ingress did not produce a source identity");

    const FullCircleEventSource2 binary64_circle_source =
        FullCircleEventSource2::from_binary64(
            0.0,
            0.0,
            2.0,
            0.0,
            false,
            0.5,
            2.0);
    const std::string expected_binary64_circle_bytes =
        encode_string_sequence({
            "full-circle-event-source-binary64-v1",
            canonical_encode_binary64(0.0),
            canonical_encode_binary64(0.0),
            canonical_encode_binary64(2.0),
            canonical_encode_binary64(0.0),
            "counterclockwise",
            canonical_encode_binary64(0.5),
            canonical_encode_binary64(2.0),
        });
    require(
        binary64_circle_source.canonical_bytes()
            == expected_binary64_circle_bytes,
        "legacy full-circle binary64 source identity drifted");

    // Production mutation caught: accepting a nonpositive squared spatial
    // floor erases the unit-bearing termination invariant.
    for (const Epeck::FT& invalid : {Epeck::FT(0), Epeck::FT(-1)}) {
        bool rejected = false;
        try {
            static_cast<void>(AuditSquaredSpatialFloorMm2::build(invalid));
        } catch (const AuditSquaredSpatialFloorError&) {
            rejected = true;
        }
        require(rejected, "nonpositive squared spatial floor was accepted");
    }

    bool depth_rejected = false;
    try {
        static_cast<void>(limits(
            Epeck::FT(1),
            std::numeric_limits<std::size_t>::max(),
            1));
    } catch (const AuditDecisionDepthLimitError&) {
        depth_rejected = true;
    }
    require(depth_rejected, "unbounded decision depth lacked named rejection");

    bool zero_nodes_rejected = false;
    try {
        static_cast<void>(limits(Epeck::FT(1), 0, 0));
    } catch (const AuditDecisionNodeLimitError&) {
        zero_nodes_rejected = true;
    }
    require(zero_nodes_rejected, "zero node budget lacked named rejection");

    bool nodes_rejected = false;
    try {
        static_cast<void>(limits(
            Epeck::FT(1),
            1,
            std::numeric_limits<std::size_t>::max()));
    } catch (const AuditDecisionNodeLimitError&) {
        nodes_rejected = true;
    }
    require(nodes_rejected, "unbounded node budget lacked named rejection");
}

void legacy_full_circle_identity_gate()
{
    // Production mutation caught: exact delegation retags the legacy
    // binary64 motion or effective-cap identity.
    const FullCircleEventSource2 source =
        FullCircleEventSource2::from_binary64(
            5.0, 5.0, 1.0, 0.0, false, 0.5, 4.0);
    require(
        sha256_bytes(source.motion_identity_bytes())
            == bytes_from_hex(
                "8c74544c7b8ee1e95a43cf71a1fb016d"
                "f3bea1bd31153aafae2c068cad932b50"),
        "legacy binary64 full-circle motion identity drifted");
    require(
        sha256_bytes(source.cap_identity_bytes())
            == bytes_from_hex(
                "1b885d0a4746ff7f86899a3915885062"
                "5008cd693bcd4a29b863d8e297791f05"),
        "legacy binary64 full-circle cap identity drifted");

    Stock2 clear(rectangle(0.0, 0.0, 10.0, 10.0), {});
    clear.subtract_disk(5.0, 5.0, 100.0);
    Stock2 virgin(rectangle(0.0, 0.0, 10.0, 10.0), {});
    Stock2 partial(rectangle(0.0, 0.0, 10.0, 10.0), {});
    partial.subtract_disk(5.0, 5.0, 1.375);
    for (const Stock2* stock : {&clear, &virgin, &partial}) {
        const FullCircleTeaAudit2 audit =
            audit_full_circle_tea_event_exact(*stock, source);
        require(
            audit.trace.motion_identity == source.motion_identity_bytes()
                && audit.trace.effective_cap_bytes
                    == source.cap_identity_bytes(),
            "legacy wrapper trace lost binary64 motion or cap identity");
        require(
            audit.trace.canonical_digest
                == sha256_bytes(audit.trace.canonical_bytes),
            "legacy wrapper trace digest does not authenticate canonical bytes");
    }

    const FullCircleEventSource2 signed_zero =
        FullCircleEventSource2::from_binary64(
            -0.0, 5.0, 1.0, 0.0, false, 0.5, 4.0);
    const FullCircleEventSource2 clockwise =
        FullCircleEventSource2::from_binary64(
            5.0, 5.0, 1.0, 0.0, true, 0.5, 4.0);
    const FullCircleEventSource2 next_cap =
        FullCircleEventSource2::from_binary64(
            5.0,
            5.0,
            1.0,
            0.0,
            false,
            0.5,
            std::nextafter(4.0, 0.0));
    require(
        signed_zero.motion_identity_bytes() != source.motion_identity_bytes()
            && clockwise.motion_identity_bytes()
                != source.motion_identity_bytes()
            && next_cap.cap_identity_bytes() != source.cap_identity_bytes(),
        "legacy wrapper collapsed signed-zero, direction, or nextafter identity");
}

void named_ingress_error_gate()
{
    // Production mutations caught: unchecked NaN/Inf reaches exact injection,
    // or a generic invalid_argument erases the boundary failure domain.
    Stock2 stock(rectangle(-2.0, -2.0, 2.0, 2.0), {});
    for (double invalid : {
             std::numeric_limits<double>::quiet_NaN(),
             std::numeric_limits<double>::infinity(),
             -std::numeric_limits<double>::infinity(),
         }) {
        bool engagement_rejected = false;
        try {
            validate_engagement_input_binary64(
                invalid, 0.0, 0.5, 2.0, 0.0);
        } catch (const NonFiniteEngagementInputError&) {
            engagement_rejected = true;
        }
        require(
            engagement_rejected,
            "nonfinite legacy station input lacked named rejection");

        bool segment_rejected = false;
        try {
            validate_segment_certification_input_binary64(
                invalid, 0.0, 1.0, 0.0, 0.5, 1.0);
        } catch (const NonFiniteSegmentCertificationInputError&) {
            segment_rejected = true;
        }
        require(
            segment_rejected,
            "nonfinite legacy segment input lacked named rejection");

        bool radical_rejected = false;
        try {
            validate_mixed_radical_input_binary64(
                invalid, 0.0, 0.0, 0.0, 0.0, 0.0);
        } catch (const NonFiniteMixedRadicalInputError&) {
            radical_rejected = true;
        }
        require(
            radical_rejected,
            "nonfinite mixed-radical input lacked named rejection");
    }

    bool tool_rejected = false;
    try {
        validate_engagement_input_binary64(
            0.0, 0.0, 0.0, 2.0, 0.0);
    } catch (const NonPositiveEngagementToolRadiusError&) {
        tool_rejected = true;
    }
    require(tool_rejected, "nonpositive station tool radius was accepted");

    bool cap_rejected = false;
    try {
        validate_engagement_input_binary64(
            0.0, 0.0, 0.5, 0.0, 0.0);
    } catch (const InvalidEngagementCapRatioError&) {
        cap_rejected = true;
    }
    require(cap_rejected, "out-of-range station cap ratio was accepted");

    bool gap_rejected = false;
    try {
        validate_engagement_input_binary64(
            0.0, 0.0, 0.5, 2.0, 5.0);
    } catch (const InvalidEngagementGapRatioError&) {
        gap_rejected = true;
    }
    require(gap_rejected, "out-of-range station gap ratio was accepted");

    bool radical_root_rejected = false;
    try {
        validate_mixed_radical_input_binary64(
            0.0, 0.0, 0.0, 0.0, -1.0, 0.0);
    } catch (const InvalidMixedRadicalRootError&) {
        radical_root_rejected = true;
    }
    require(
        radical_root_rejected,
        "negative mixed-radical root was accepted");
}

void exact_pi_equality_and_shared_root_gate()
{
    // Production mutation caught: non-strict cap equality or independent
    // coordinate roots corrupt exact station classification.
    const AuditPolicy2 pi_policy = policy(std::numbers::pi);
    require(
        pi_policy.engagement_cap().chord_ratio() == Epeck::FT(4),
        "pi policy did not retain the exact squared-chord surrogate four");
    Stock2 half_rim(rectangle(-2.0, -2.0, 2.0, 0.0), {});

    require(
        replay_audit_unguarded_station_exact(
            half_rim,
            EPoint(0, 0),
            pi_policy.tool_radius_mm(),
            Epeck::FT(4))
            == AuditExactStationDisposition2::WITHIN_CAP,
        "exact pi half-rim equality was reported exceeded");

    Stock2 clear_stock(rectangle(20.0, 20.0, 30.0, 30.0), {});
    const AuditExactStationClassification2 clear_station =
        classify_audit_unguarded_station_exact(
            clear_stock,
            EPoint(0, 0),
            Epeck::FT(1) / Epeck::FT(2),
            Epeck::FT(2));
    require(
        clear_station.disposition()
                == AuditExactStationDisposition2::WITHIN_CAP
            && clear_station.true_run_arcs().empty(),
        "clear exact station fabricated a material run");

    using CoordNT = GpsPoint::CoordNT;
    const GpsPoint legitimate(
        CoordNT(Epeck::FT(0), Epeck::FT(1), Epeck::FT(3)),
        CoordNT(Epeck::FT(1) / Epeck::FT(4)));
    const ExactOneRootPoint2 decomposed =
        decompose_exact_one_root_point(legitimate);
    require(
        decomposed.root() == Epeck::FT(3),
        "shared production one-root decomposition lost its exact radicand");

    const GpsPoint mismatched(
        CoordNT(Epeck::FT(0), Epeck::FT(1), Epeck::FT(2)),
        CoordNT(Epeck::FT(0), Epeck::FT(1), Epeck::FT(3)));
    bool mismatch_rejected = false;
    try {
        static_cast<void>(decompose_exact_one_root_point(mismatched));
    } catch (const ExactOneRootCoordinateMismatchError&) {
        mismatch_rejected = true;
    }
    require(
        mismatch_rejected,
        "cross-root point lacked named release-mode rejection");

    Stock2 one_root_case(rectangle(-2.0, -2.0, 2.0, 0.25), {});
    require(
        replay_audit_unguarded_station_exact(
            one_root_case,
            EPoint(0, 0),
            Epeck::FT(1) / Epeck::FT(2),
            Epeck::FT(2))
            == AuditExactStationDisposition2::CAP_EXCEEDED,
        "real one-root engagement case did not replay through shared decomposition");

    Stock2 generic_crossing(rectangle(-3.0, -3.0, 3.0, 3.0), {});
    generic_crossing.subtract_disk(0.25, 0.375, 0.75);
    const AuditExactStationClassification2 classified =
        classify_audit_unguarded_station_exact(
            generic_crossing,
            EPoint(
                Epeck::FT(-1) / Epeck::FT(8),
                Epeck::FT(1) / Epeck::FT(7)),
            Epeck::FT(5) / Epeck::FT(8),
            Epeck::FT(2));
    bool found_shared_extended_root = false;
    for (const GpsPoint& crossing : classified.boundary_intersections()) {
        if (!crossing.x().is_extended() || !crossing.y().is_extended()) {
            continue;
        }
        const ExactOneRootPoint2 exact_crossing =
            decompose_exact_one_root_point(crossing);
        require(
            exact_crossing.root() == crossing.x().root()
                && exact_crossing.root() == crossing.y().root(),
            "production station classifier split a shared algebraic root");
        found_shared_extended_root = true;
    }
    require(
        found_shared_extended_root,
        "generic exact circle/circle station fixture exposed no two-coordinate root");
}

} // namespace

void audit_certification_contract_gate()
{
    exact_ingress_and_limit_errors_gate();
    legacy_full_circle_identity_gate();
    named_ingress_error_gate();
    exact_pi_equality_and_shared_root_gate();
}
