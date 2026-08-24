#include "audit_arc_motion_2.h"
#include "audit_digest_2.h"
#include "exact_depletion_2.h"

#include <algorithm>
#include <numbers>
#include <limits>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

void require(bool condition, const char* message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

AuditArcMotion2 quarter_arc()
{
    return AuditArcMotion2::build(
        EPoint(Epeck::FT(1), Epeck::FT(2)),
        EVector(Epeck::FT(3), Epeck::FT(4)),
        Epeck::FT(5),
        0.0,
        std::numbers::pi / 2.0,
        false,
        Epeck::FT(0));
}

AuditArcMotion2 arc(double start, double end, bool clockwise)
{
    return AuditArcMotion2::build(
        EPoint(Epeck::FT(2), Epeck::FT(-3)),
        EVector(Epeck::FT(3), Epeck::FT(4)),
        Epeck::FT(5),
        start,
        end,
        clockwise,
        Epeck::FT(0));
}

void require_constructed_arc(
    const AuditArcMotion2& motion,
    bool full_turn)
{
    const Epeck::FT max_chord = Epeck::FT(1) / Epeck::FT(4);
    const ExactArcDepletionConstruction2 construction =
        construct_exact_arc_depletion(
            motion,
            Epeck::FT(1),
            max_chord,
            4096);
    require(
        !construction.trace.cyclic(),
        "partial-arc trace became cyclic");
    require(
        !full_turn
            || (construction.centers.front() == construction.centers.back()
                && construction.trace.parameters().front().chart()
                    == construction.trace.parameters().back().chart()
                && CGAL::compare(
                       construction.trace.parameters().front().parameter(),
                       construction.trace.parameters().back().parameter())
                    == CGAL::EQUAL),
        "full-turn partial arc omitted its explicit terminal anchor");
    require(
        exact_arc_structural_density_holds(
            motion,
            max_chord,
            construction.trace.parameters()),
        "arc variant lost canonical structural density");
    for (const EPoint& center : construction.centers) {
        require(
            exact_arc_point_is_incident(motion, center),
            "arc variant emitted an off-guide center");
    }
}

} // namespace

namespace {

void require_trace_rejected(
    const AuditArcMotion2& motion,
    std::vector<ExactCircleChartParameter2> parameters,
    const char* message)
{
    bool rejected = false;
    try {
        static_cast<void>(ExactArcDepletionTrace2::build(
            motion,
            Epeck::FT(1),
            Epeck::FT(1) / Epeck::FT(4),
            4096,
            std::move(parameters)));
    } catch (const ExactDepletionConstructionError&) {
        rejected = true;
    }
    require(rejected, message);
}

std::string u64(std::size_t value)
{
    std::string result;
    for (int shift = 56; shift >= 0; shift -= 8) {
        result.push_back(static_cast<char>(
            (static_cast<std::uint64_t>(value) >> shift) & 0xffU));
    }
    return result;
}

std::string node(char kind, const std::string& payload)
{
    return std::string("CCAN\0\1", 6) + kind + u64(payload.size()) + payload;
}

std::string integer(unsigned char sign, unsigned char magnitude)
{
    return node('I', std::string(1, static_cast<char>(sign))
        + std::string(1, static_cast<char>(magnitude)));
}

} // namespace

void exact_arc_depletion_gate()
{
    const AuditArcMotion2 motion = quarter_arc();
    const Epeck::FT tool_radius(1);
    const Epeck::FT max_chord = Epeck::FT(1) / Epeck::FT(4);
    const ExactArcDepletionConstruction2 construction =
        construct_exact_arc_depletion(
            motion,
            tool_radius,
            max_chord,
            4096);

    require(
        construction.centers.front() == motion.start_point(),
        "arc depletion omitted exact start anchor");
    require(
        construction.centers.back() == motion.end_point(),
        "arc depletion omitted exact end anchor");
    require(
        exact_arc_structural_density_holds(
            motion,
            max_chord,
            construction.trace.parameters()),
        "arc depletion trace does not prove structural density");
    require(
        construction.trace.matches_exact_inputs(
            tool_radius,
            max_chord,
            4096),
        "arc depletion trace does not bind exact policy inputs");
    require(
        construction.trace.matches_motion(motion),
        "arc depletion trace does not bind native motion identity");

    bool radius_relation_rejected = false;
    try {
        static_cast<void>(construct_exact_arc_depletion(
            motion,
            tool_radius,
            tool_radius,
            4096));
    } catch (const ExactArcDepletionPolicyError&) {
        radius_relation_rejected = true;
    }
    require(
        radius_relation_rejected,
        "arc depletion accepted chord bound at tool radius");

    bool center_limit_rejected = false;
    try {
        static_cast<void>(construct_exact_arc_depletion(
            motion,
            tool_radius,
            max_chord,
            2));
    } catch (const ExactDepletionCenterLimitError&) {
        center_limit_rejected = true;
    }
    require(center_limit_rejected, "arc center limit was not enforced before allocation");

    for (const AuditArcMotion2& variant : {
             arc(0.0, std::numbers::pi / 2.0, false),
             arc(0.37, 4.91, false),
             arc(std::numbers::pi / 2.0, 0.0, true),
             arc(4.91, 0.37, true),
         }) {
        require_constructed_arc(variant, false);
    }
    require_constructed_arc(
        arc(0.37, 0.37 + std::numbers::pi * 2.0, false),
        true);
    require_constructed_arc(
        arc(1.2, 1.2 - std::numbers::pi * 2.0, true),
        true);

    for (const double seam : {
             0.0,
             std::numbers::pi / 2.0,
             std::numbers::pi,
             3.0 * std::numbers::pi / 2.0,
         }) {
        require_constructed_arc(
            arc(seam, seam + std::numbers::pi / 2.0, false),
            false);
        require_constructed_arc(
            arc(seam, seam - std::numbers::pi / 2.0, true),
            false);
    }

    for (const auto& invalid : {
             std::pair<Epeck::FT, Epeck::FT>{Epeck::FT(0), max_chord},
             std::pair<Epeck::FT, Epeck::FT>{Epeck::FT(-1), max_chord},
             std::pair<Epeck::FT, Epeck::FT>{tool_radius, Epeck::FT(0)},
             std::pair<Epeck::FT, Epeck::FT>{tool_radius, Epeck::FT(-1)},
         }) {
        bool rejected = false;
        try {
            static_cast<void>(construct_exact_arc_depletion(
                motion,
                invalid.first,
                invalid.second,
                4096));
        } catch (const ExactArcDepletionPolicyError&) {
            rejected = true;
        }
        require(
            rejected,
            "nonpositive arc depletion policy lacked named rejection");
    }

    const std::vector<ExactCircleChartParameter2> canonical =
        construction.trace.parameters();
    std::vector<ExactCircleChartParameter2> reordered = canonical;
    std::swap(reordered[1], reordered[2]);
    require_trace_rejected(motion, std::move(reordered), "reordered arc parameters were accepted");

    std::vector<ExactCircleChartParameter2> missing = canonical;
    missing.erase(missing.begin() + 1);
    require_trace_rejected(motion, std::move(missing), "missing arc parameter was accepted");

    std::vector<ExactCircleChartParameter2> duplicated = canonical;
    duplicated.insert(duplicated.begin() + 1, duplicated[1]);
    require_trace_rejected(motion, std::move(duplicated), "duplicated arc parameter was accepted");

    std::vector<ExactCircleChartParameter2> complement = canonical;
    complement[1] = ExactCircleChartParameter2::build(
        (complement[1].chart() + 2) % 4,
        complement[1].parameter());
    require_trace_rejected(motion, std::move(complement), "complement arc parameter was accepted");

    const AuditArcMotion2 opposite_direction_complement =
        AuditArcMotion2::build(
            EPoint(Epeck::FT(1), Epeck::FT(2)),
            EVector(Epeck::FT(3), Epeck::FT(4)),
            Epeck::FT(5),
            0.0,
            std::numbers::pi / 2.0 - std::numbers::pi * 2.0,
            true,
            Epeck::FT(0));
    const ExactArcDepletionConstruction2 complement_construction =
        construct_exact_arc_depletion(
            opposite_direction_complement,
            tool_radius,
            max_chord,
            4096);
    require(
        complement_construction.centers.front() == construction.centers.front()
            && complement_construction.centers.back() == construction.centers.back(),
        "test complement does not share the quarter arc endpoints");
    require_trace_rejected(
        motion,
        complement_construction.trace.parameters(),
        "complete opposite-direction complement trace was accepted");

    const auto seam = std::find_if(
        canonical.begin() + 1,
        canonical.end(),
        [](const ExactCircleChartParameter2& parameter) {
            return CGAL::sign(parameter.parameter()) == CGAL::ZERO;
        });
    require(seam != canonical.end(), "test arc did not cross an owned seam");
    std::vector<ExactCircleChartParameter2> non_owner = canonical;
    const std::size_t seam_index = static_cast<std::size_t>(seam - canonical.begin());
    non_owner[seam_index] = ExactCircleChartParameter2::build(
        (seam->chart() + 3) % 4,
        Epeck::FT(1));
    require_trace_rejected(
        motion,
        std::move(non_owner),
        "geometrically equal non-owning seam representation was accepted");

    const AuditArcMotion2 translated_same_path = AuditArcMotion2::build(
        EPoint(Epeck::FT(7), Epeck::FT(-3)),
        EVector(Epeck::FT(3), Epeck::FT(4)),
        Epeck::FT(5),
        0.0,
        std::numbers::pi / 2.0,
        false,
        Epeck::FT(0));
    require(
        !construction.trace.matches_motion(translated_same_path),
        "translated motion with the same chart traversal matched a foreign trace");

    const AuditArcMotion2 full_motion = arc(
        -0.37,
        -0.37 + std::numbers::pi * 2.0,
        false);
    const ExactArcDepletionConstruction2 full_construction =
        construct_exact_arc_depletion(
            full_motion,
            tool_radius,
            max_chord,
            4096);
    require(!full_construction.trace.cyclic(), "full-turn partial arc trace is cyclic");
    std::vector<ExactCircleChartParameter2> missing_terminal =
        full_construction.trace.parameters();
    missing_terminal.pop_back();
    require_trace_rejected(
        full_motion,
        std::move(missing_terminal),
        "full-turn trace without terminal anchor was accepted");
    std::vector<ExactCircleChartParameter2> extra_terminal =
        full_construction.trace.parameters();
    extra_terminal.push_back(extra_terminal.back());
    require_trace_rejected(
        full_motion,
        std::move(extra_terminal),
        "full-turn trace with an extra terminal anchor was accepted");

    require(
        construction.trace.strategy_version()
            == "exact-arc-pythagorean-guide-v1",
        "arc trace strategy is not versioned");
    require(
        construction.trace.digest().bytes().size() == 32,
        "arc trace witness digest is not SHA-256 sized");
    const ExactArcDepletionConstruction2 tighter = construct_exact_arc_depletion(
        motion,
        tool_radius,
        Epeck::FT(1) / Epeck::FT(8),
        4096);
    require(
        tighter.trace.digest().bytes() != construction.trace.digest().bytes(),
        "arc trace digest ignores exact depletion policy");
    const ExactArcDepletionConstruction2 other_tool = construct_exact_arc_depletion(
        motion,
        Epeck::FT(2),
        max_chord,
        4096);
    const ExactArcDepletionConstruction2 other_limit = construct_exact_arc_depletion(
        motion,
        tool_radius,
        max_chord,
        8192);
    const ExactArcDepletionConstruction2 other_motion = construct_exact_arc_depletion(
        arc(0.37, 1.2, false),
        tool_radius,
        max_chord,
        4096);
    require(
        other_tool.trace.digest().bytes() != construction.trace.digest().bytes(),
        "arc trace digest ignores tool radius");
    require(
        other_limit.trace.digest().bytes() != construction.trace.digest().bytes(),
        "arc trace digest ignores center-count limit");
    require(
        other_motion.trace.digest().bytes() != construction.trace.digest().bytes(),
        "arc trace digest ignores native motion identity");
    require(
        construction.trace.canonical_bytes().find(
            construction.trace.strategy_version()) != std::string::npos,
        "arc trace canonical bytes omit strategy version");

    const std::string expected_half = node(
        'R',
        integer(0, 1) + integer(0, 2));
    require(
        canonical_audit_rational_bytes(Epeck::FT(1) / Epeck::FT(2))
            == expected_half,
        "audit exact rational bytes diverge from CCAN rational encoding");
    require(
        canonical_decode_rational(
            canonical_audit_rational_bytes(
                Epeck::FT(-3) / Epeck::FT(5)))
            == CORE::BigRat(-3, 5),
        "negative audit rational does not round-trip through CCAN");
    require(
        canonical_decode_rational(
            canonical_audit_rational_bytes(Epeck::FT(0)))
            == CORE::BigRat(0),
        "zero audit rational does not round-trip through CCAN");
    CORE::BigInt big = 1;
    big <<= 130;
    const Epeck::FT exact_big{Epeck::FT::Exact_type(big)};
    require(
        canonical_decode_rational(canonical_audit_rational_bytes(exact_big))
            == CORE::BigRat(big),
        "large audit integer does not round-trip through CCAN");
    require(
        canonical_audit_binary64_bytes(-3.5)
            == canonical_encode_binary64(-3.5),
        "audit binary64 bytes diverge from canonical CCAN encoding");

    bool input_digest_size_rejected = false;
    try {
        static_cast<void>(
            AuditInputDigestAuthority2::from_external_bytes("short"));
    } catch (const AuditDigestSizeError&) {
        input_digest_size_rejected = true;
    }
    require(
        input_digest_size_rejected,
        "malformed external input digest lacked named rejection");
    bool operation_digest_size_rejected = false;
    try {
        static_cast<void>(
            AuthenticatedOperationDigestAuthority2::from_external_bytes(
                "short"));
    } catch (const AuditDigestSizeError&) {
        operation_digest_size_rejected = true;
    }
    require(
        operation_digest_size_rejected,
        "malformed external operation digest lacked named rejection");

    for (const auto& invalid : {
             std::tuple<Epeck::FT, Epeck::FT, std::size_t>{
                 Epeck::FT(0), max_chord, 4096},
             std::tuple<Epeck::FT, Epeck::FT, std::size_t>{
                 tool_radius, Epeck::FT(0), 4096},
             std::tuple<Epeck::FT, Epeck::FT, std::size_t>{
                 tool_radius, max_chord, 0},
         }) {
        bool rejected = false;
        try {
            static_cast<void>(ExactArcDepletionTrace2::build(
                motion,
                std::get<0>(invalid),
                std::get<1>(invalid),
                std::get<2>(invalid),
                canonical));
        } catch (const ExactDepletionConstructionError&) {
            rejected = true;
        }
        require(rejected, "trace factory bypassed its own policy validation");
    }
}
