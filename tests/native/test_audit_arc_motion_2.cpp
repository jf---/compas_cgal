#include "audit_arc_motion_2.h"
#include "exact_circle_chart_2.h"

#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

void exact_arc_depletion_gate();

namespace {

void require(bool condition, const char* message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

AuditArcMotion2 arc(
    double start,
    double end,
    bool clockwise)
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

void exact_arc_motion_gate()
{
    require(
        audit_arc_motion_strategy_version()
            == "audit-arc-quarter-chart-binary64-v2",
        "arc motion identity does not name the authored-sweep seam version");
    bool nonfinite_rejected = false;
    try {
        static_cast<void>(arc(
            0.0,
            std::numeric_limits<double>::infinity(),
            false));
    } catch (const AuditArcNonFiniteInputError&) {
        nonfinite_rejected = true;
    }
    require(nonfinite_rejected, "non-finite arc input lacked named rejection");

    bool incidence_rejected = false;
    try {
        static_cast<void>(AuditArcMotion2::build(
            EPoint(0, 0),
            EVector(1, 0),
            Epeck::FT(2),
            0.0,
            1.0,
            false,
            Epeck::FT(0)));
    } catch (const AuditArcRadiusIncidenceError&) {
        incidence_rejected = true;
    }
    require(incidence_rejected, "off-radius zero phase lacked named rejection");

    bool extent_rejected = false;
    try {
        static_cast<void>(arc(0.0, 0.0, false));
    } catch (const AuditArcSweepExtentError&) {
        extent_rejected = true;
    }
    require(extent_rejected, "zero sweep lacked named extent rejection");

    bool direction_rejected = false;
    try {
        static_cast<void>(ExactArcChartInterval2::build(
            0,
            Epeck::FT(3) / Epeck::FT(4),
            Epeck::FT(1) / Epeck::FT(4),
            true,
            false,
            false));
    } catch (const AuditArcIntervalDirectionError&) {
        direction_rejected = true;
    }
    require(direction_rejected, "foreign interval direction lacked named rejection");

    bool seam_ownership_rejected = false;
    try {
        static_cast<void>(ExactArcChartInterval2::build(
            0,
            Epeck::FT(0),
            Epeck::FT(1) / Epeck::FT(4),
            true,
            false,
            false));
    } catch (const AuditArcSeamOwnershipError&) {
        seam_ownership_rejected = true;
    }
    require(
        seam_ownership_rejected,
        "foreign seam ownership lacked named rejection");

    const AuditArcMotion2 awkward = arc(0.37, 4.91, false);
    require(
        exact_arc_point_is_incident(awkward, awkward.start_point()),
        "awkward arc start left exact guide");
    require(
        exact_arc_point_is_incident(awkward, awkward.end_point()),
        "awkward arc end left exact guide");
    require(
        awkward.intervals().size() == 4,
        "major CCW arc lost a clipped chart span");
    require(
        awkward.digest().bytes().size() == 32,
        "arc motion digest is not SHA-256 sized");

    const AuditArcMotion2 negative = arc(-std::numbers::pi / 2.0, 0.0, false);
    const AuditArcMotion2 multiturn = arc(3.0 * std::numbers::pi / 2.0, 2.0 * std::numbers::pi, false);
    require(
        negative.start_point() == multiturn.start_point()
            && negative.end_point() == multiturn.end_point(),
        "equivalent normalized arcs disagree geometrically");
    require(
        negative.digest().bytes() != multiturn.digest().bytes(),
        "authored-angle identity was erased by normalization");

    const AuditArcMotion2 ccw = arc(0.0, std::numbers::pi / 2.0, false);
    const AuditArcMotion2 cw = arc(std::numbers::pi / 2.0, 0.0, true);
    require(
        ccw.digest().bytes() != cw.digest().bytes(),
        "direction does not affect native motion identity");
    require(
        ccw.intervals().front().increasing()
            && !cw.intervals().front().increasing(),
        "chart traversal direction contradicts authored orientation");

    for (const double ordinary_start : {0.37, 1.2, -0.37}) {
        const AuditArcMotion2 full = arc(
            ordinary_start,
            ordinary_start + std::numbers::pi * 2.0,
            false);
        const AuditArcMotion2 full_cw = arc(
            ordinary_start,
            ordinary_start - std::numbers::pi * 2.0,
            true);
        require(
            full.full_turn() && full_cw.full_turn(),
            "ordinary non-seam authored full turn was direction-asymmetric");
        require(
            full.start_point() == full.end_point()
                && full_cw.start_point() == full_cw.end_point(),
            "ordinary non-seam full turn did not close exactly");
        require(
            full.intervals().size() == 5
                && full_cw.intervals().size() == 5,
            "ordinary non-seam full turn lost a clipped chart span");
    }

    bool collapsed_rejected = false;
    try {
        constexpr double COLLAPSING_START = 1.56;
        static_cast<void>(arc(
            COLLAPSING_START,
            std::nextafter(
                COLLAPSING_START,
                std::numeric_limits<double>::infinity()),
            false));
    } catch (const AuditArcEndpointCollapseError&) {
        collapsed_rejected = true;
    }
    require(collapsed_rejected, "rounded chart endpoint collapse was accepted");

    const std::vector<double> seams{
        0.0,
        std::numbers::pi / 2.0,
        std::numbers::pi,
        3.0 * std::numbers::pi / 2.0,
    };
    for (const double seam : seams) {
        const AuditArcMotion2 seam_ccw = arc(
            seam,
            seam + std::numbers::pi / 2.0,
            false);
        const AuditArcMotion2 seam_cw = arc(
            seam,
            seam - std::numbers::pi / 2.0,
            true);
        require(
            seam_ccw.intervals().front().owns_start_seam()
                && !seam_ccw.intervals().back().owns_end_seam(),
            "CCW seam ownership is not parameter-zero canonical");
        require(
            !seam_cw.intervals().front().owns_start_seam()
                && seam_cw.intervals().back().owns_end_seam(),
            "CW seam ownership is not parameter-zero canonical");
    }

    const double quarter = std::numbers::pi / 2.0;
    for (const double adjacent : {
             std::nextafter(quarter, -std::numeric_limits<double>::infinity()),
             std::nextafter(quarter, std::numeric_limits<double>::infinity()),
         }) {
        const AuditArcMotion2 near_seam = arc(adjacent, adjacent + 0.25, false);
        require(
            exact_arc_point_is_incident(near_seam, near_seam.start_point())
                && exact_arc_point_is_incident(near_seam, near_seam.end_point()),
            "nextafter seam arc left its exact rational guide");
    }

    const AuditArcMotion2 repeat = arc(0.37, 4.91, false);
    const AuditArcMotion2 translated = AuditArcMotion2::build(
        EPoint(Epeck::FT(5), Epeck::FT(-3)),
        EVector(Epeck::FT(3), Epeck::FT(4)),
        Epeck::FT(5),
        0.37,
        4.91,
        false,
        Epeck::FT(0));
    const AuditArcMotion2 scaled = AuditArcMotion2::build(
        EPoint(Epeck::FT(2), Epeck::FT(-3)),
        EVector(Epeck::FT(6), Epeck::FT(8)),
        Epeck::FT(10),
        0.37,
        4.91,
        false,
        Epeck::FT(0));
    require(
        repeat.digest().bytes() == awkward.digest().bytes(),
        "arc digest is not repeatable");
    require(
        translated.digest().bytes() != awkward.digest().bytes()
            && scaled.digest().bytes() != awkward.digest().bytes(),
        "arc digest ignores rational translation or scale");
    const AuditArcMotion2 changed_start = arc(
        std::nextafter(0.37, std::numeric_limits<double>::infinity()),
        4.91,
        false);
    const AuditArcMotion2 changed_end = arc(
        0.37,
        std::nextafter(4.91, -std::numeric_limits<double>::infinity()),
        false);
    require(
        changed_start.digest().bytes() != awkward.digest().bytes()
            && changed_end.digest().bytes() != awkward.digest().bytes(),
        "arc digest ignores a public interval-parameter mutation");

    for (int chart = 0; chart < 4; ++chart) {
        require(
            exact_circle_chart_id(chart)
                == "center-quarter-" + std::to_string(chart) + "-v1",
            "shared frozen atlas id drifted");
        const ExactCircleChartAtlasRecord2& coefficients =
            exact_circle_chart_record(chart);
        require(
            coefficients.denominator == std::array<int, 3>{1, 0, 1},
            "shared atlas denominator drifted");
    }

    const Epeck::FT fifth = Epeck::FT(1) / Epeck::FT(5);
    const Epeck::FT twenty_fifth = Epeck::FT(1) / Epeck::FT(25);
    const std::array<Epeck::FT, 4> parameters{
        Epeck::FT(0),
        Epeck::FT(1) / Epeck::FT(7),
        Epeck::FT(1) / Epeck::FT(2),
        Epeck::FT(1),
    };
    const std::array<std::array<EPoint, 4>, 4> expected{{
        {EPoint(1, 0), EPoint(Epeck::FT(24) * twenty_fifth, Epeck::FT(7) * twenty_fifth), EPoint(Epeck::FT(3) * fifth, Epeck::FT(4) * fifth), EPoint(0, 1)},
        {EPoint(0, 1), EPoint(Epeck::FT(-7) * twenty_fifth, Epeck::FT(24) * twenty_fifth), EPoint(Epeck::FT(-4) * fifth, Epeck::FT(3) * fifth), EPoint(-1, 0)},
        {EPoint(-1, 0), EPoint(Epeck::FT(-24) * twenty_fifth, Epeck::FT(-7) * twenty_fifth), EPoint(Epeck::FT(-3) * fifth, Epeck::FT(-4) * fifth), EPoint(0, -1)},
        {EPoint(0, -1), EPoint(Epeck::FT(7) * twenty_fifth, Epeck::FT(-24) * twenty_fifth), EPoint(Epeck::FT(4) * fifth, Epeck::FT(-3) * fifth), EPoint(1, 0)},
    }};
    for (int chart = 0; chart < 4; ++chart) {
        for (std::size_t index = 0; index < parameters.size(); ++index) {
            require(
                exact_circle_chart_point(
                    EPoint(0, 0),
                    EVector(1, 0),
                    ExactCircleChartParameter2::build(chart, parameters[index]))
                    == expected[static_cast<std::size_t>(chart)][index],
                "shared atlas polynomial disagrees with exact evaluator");
        }
    }
}

} // namespace

int main()
{
    exact_arc_motion_gate();
    exact_arc_depletion_gate();
    return 0;
}
