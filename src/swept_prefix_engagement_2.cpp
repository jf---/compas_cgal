#include "swept_prefix_engagement_2.h"

#include "engagement_2.h"
#include "exact_rim_order_2.h"

#include <algorithm>
#include <functional>
#include <utility>
#include <vector>

#include <CGAL/enum.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/number_utils.h>

namespace {

using FT = Epeck::FT;
using CoordNT = GpsPoint::CoordNT;
using RimArc = std::pair<GpsPoint, GpsPoint>;

GpsPoint rational_witness_on_ccw_arc(
    const GpsPoint& start,
    const GpsPoint& end,
    const EPoint& center,
    const FT& radius)
{
    CoordNT minimum_x = start.x();
    CoordNT maximum_x = end.x();
    if (CGAL::compare(minimum_x, maximum_x) == CGAL::LARGER) {
        std::swap(minimum_x, maximum_x);
    }
    if (CGAL::compare(minimum_x, maximum_x) == CGAL::EQUAL) {
        throw SweptPrefixRimCellError(
            "distinct swept-prefix rim-cell endpoints share an x coordinate");
    }
    FT lower = center.x() - radius;
    FT upper = center.x() + radius;
    FT witness_x;
    while (true) {
        witness_x = (lower + upper) / FT(2);
        const CoordNT candidate(witness_x);
        if (CGAL::compare(candidate, minimum_x) == CGAL::LARGER
            && CGAL::compare(candidate, maximum_x) == CGAL::SMALLER) {
            break;
        }
        if (CGAL::compare(candidate, minimum_x) != CGAL::LARGER) {
            lower = witness_x;
        } else {
            upper = witness_x;
        }
    }
    const FT offset_x = witness_x - center.x();
    const FT root = radius * radius - offset_x * offset_x;
    if (CGAL::sign(root) != CGAL::POSITIVE) {
        throw SweptPrefixRimCellError(
            "swept-prefix rational witness radicand must be positive");
    }
    const bool upper_branch =
        CGAL::compare(start.x(), end.x()) == CGAL::LARGER;
    const GpsPoint witness(
        CoordNT(witness_x),
        CoordNT(center.y(), upper_branch ? FT(1) : FT(-1), root));
    if (!strictly_inside_ccw_arc(witness, start, end, center)) {
        throw SweptPrefixRimCellError(
            "failed to construct an exact rational witness inside a swept-prefix rim cell");
    }
    return witness;
}

GpsPoint radial_rim_point(
    const EPoint& rim_center,
    const FT& vector_x,
    const FT& vector_y,
    const FT& radius,
    bool outward)
{
    const FT squared_length = vector_x * vector_x + vector_y * vector_y;
    if (CGAL::is_zero(squared_length)) {
        throw SweptPrefixRimCellError("swept-prefix radial direction is degenerate");
    }
    const FT sign = outward ? FT(1) : FT(-1);
    return GpsPoint(
        CoordNT(
            rim_center.x(),
            sign * radius * vector_x / squared_length,
            squared_length),
        CoordNT(
            rim_center.y(),
            sign * radius * vector_y / squared_length,
            squared_length));
}

std::vector<GpsPoint> equal_circle_intersections(
    const EPoint& first,
    const EPoint& second,
    const FT& radius)
{
    const FT dx = second.x() - first.x();
    const FT dy = second.y() - first.y();
    const FT distance_squared = dx * dx + dy * dy;
    if (CGAL::is_zero(distance_squared)) return {};
    const FT height_squared =
        radius * radius - distance_squared / FT(4);
    if (CGAL::sign(height_squared) == CGAL::NEGATIVE) return {};
    const FT middle_x = (first.x() + second.x()) / FT(2);
    const FT middle_y = (first.y() + second.y()) / FT(2);
    if (CGAL::is_zero(height_squared)) {
        return {GpsPoint(middle_x, middle_y)};
    }
    const FT root = height_squared * distance_squared;
    const FT x_coefficient = -dy / distance_squared;
    const FT y_coefficient = dx / distance_squared;
    return {
        GpsPoint(
            CoordNT(middle_x, x_coefficient, root),
            CoordNT(middle_y, y_coefficient, root)),
        GpsPoint(
            CoordNT(middle_x, -x_coefficient, root),
            CoordNT(middle_y, -y_coefficient, root)),
    };
}

std::vector<RimArc> flatten_arcs(
    const AuditExactStationClassification2& classification)
{
    std::vector<RimArc> arcs;
    for (const auto& run : classification.true_run_arcs()) {
        arcs.insert(arcs.end(), run.begin(), run.end());
    }
    return arcs;
}

std::vector<RimArc> clip_rim_arcs(
    const std::vector<RimArc>& arcs,
    const std::vector<GpsPoint>& events,
    const EPoint& center,
    const FT& radius,
    const std::function<bool(const GpsPoint&)>& covered)
{
    std::vector<RimArc> kept;
    for (const auto& [arc_start, arc_end] : arcs) {
        std::vector<GpsPoint> inside;
        for (const GpsPoint& event : events) {
            if (strictly_inside_ccw_arc(event, arc_start, arc_end, center)) {
                inside.push_back(event);
            }
        }
        std::sort(
            inside.begin(),
            inside.end(),
            [&](const GpsPoint& lhs, const GpsPoint& rhs) {
                const bool lhs_wraps = rim_ccw_less(lhs, arc_start, center);
                const bool rhs_wraps = rim_ccw_less(rhs, arc_start, center);
                if (lhs_wraps != rhs_wraps) return !lhs_wraps;
                return rim_ccw_less(lhs, rhs, center);
            });
        inside.erase(std::unique(inside.begin(), inside.end()), inside.end());
        GpsPoint piece_start = arc_start;
        inside.push_back(arc_end);
        for (const GpsPoint& piece_end : inside) {
            if (piece_start == piece_end) continue;
            const GpsPoint witness = rational_witness_on_ccw_arc(
                piece_start, piece_end, center, radius);
            if (!covered(witness)) kept.emplace_back(piece_start, piece_end);
            piece_start = piece_end;
        }
    }
    return kept;
}

} // namespace

ArcPrefixExtent2 classify_arc_prefix_extent_exact(
    const EPoint& guide_center,
    const EPoint& motion_start,
    const EPoint& station_center,
    bool clockwise)
{
    if (guide_center == motion_start || guide_center == station_center) {
        throw InvalidSweptPrefixInputError(
            "arc swept-prefix guide radius must be positive");
    }
    CGAL::Orientation turn =
        CGAL::orientation(guide_center, motion_start, station_center);
    if (clockwise) {
        turn = -turn;
    }
    if (turn == CGAL::LEFT_TURN) return ArcPrefixExtent2::AT_MOST_HALF_TURN;
    if (turn == CGAL::RIGHT_TURN) return ArcPrefixExtent2::BEYOND_HALF_TURN;
    if (CGAL::angle(motion_start, guide_center, station_center) == CGAL::OBTUSE) {
        return ArcPrefixExtent2::AT_MOST_HALF_TURN;
    }
    throw ArcPrefixExtentAmbiguousError(
        "arc swept-prefix station lies on the start ray strictly inside the motion");
}

namespace {

// A line station: its exact centre and the motion direction
// (motion_end - motion_start, not normalised).
struct LinePrefixStation2 {
    EPoint center;
    EVector direction;
};

LinePrefixStation2 line_prefix_station(
    const EPoint& motion_start,
    const EPoint& motion_end,
    const FT& parameter)
{
    if (CGAL::sign(parameter) == CGAL::NEGATIVE
        || CGAL::compare(parameter, FT(1)) == CGAL::LARGER) {
        throw InvalidSweptPrefixInputError(
            "line swept-prefix parameter must lie in [0, 1]");
    }
    const EVector direction = motion_end - motion_start;
    if (CGAL::is_zero(direction.squared_length())) {
        throw InvalidSweptPrefixInputError(
            "line swept-prefix motion must have positive length");
    }
    return {motion_start + parameter * direction, direction};
}

// The static classification clipped to the rim the line prefix has not covered.
AuditExactStationClassification2 clip_line_prefix(
    const AuditExactStationClassification2& base,
    const LinePrefixStation2& station,
    const FT& tool_radius,
    const FT& cap_chord_ratio,
    const FT& gap_close_ratio)
{
    const EPoint& center = station.center;
    const EVector normal = station.direction.perpendicular(CGAL::COUNTERCLOCKWISE);
    const std::vector<GpsPoint> events{
        radial_rim_point(center, normal.x(), normal.y(), tool_radius, true),
        radial_rim_point(center, normal.x(), normal.y(), tool_radius, false),
    };
    const auto covered = [&](const GpsPoint& point) {
        const CoordNT offset_x = point.x() - CoordNT(center.x());
        const CoordNT offset_y = point.y() - CoordNT(center.y());
        return CGAL::sign(
                   CoordNT(station.direction.x()) * offset_x
                   + CoordNT(station.direction.y()) * offset_y)
            != CGAL::POSITIVE;
    };
    return classify_audit_filtered_station_exact(
        center,
        tool_radius,
        cap_chord_ratio,
        gap_close_ratio,
        clip_rim_arcs(
            flatten_arcs(base), events, center, tool_radius, covered));
}

} // namespace

AuditExactStationClassification2 classify_line_swept_prefix_station_exact(
    const Stock2& stock,
    const EPoint& motion_start,
    const EPoint& motion_end,
    const FT& parameter,
    const FT& tool_radius,
    const FT& cap_chord_ratio,
    const FT& gap_close_ratio)
{
    const LinePrefixStation2 station =
        line_prefix_station(motion_start, motion_end, parameter);
    const AuditExactStationClassification2 base =
        classify_audit_unguarded_station_exact(
            stock, station.center, tool_radius, cap_chord_ratio, gap_close_ratio);
    if (CGAL::is_zero(parameter)) return base;
    return clip_line_prefix(
        base, station, tool_radius, cap_chord_ratio, gap_close_ratio);
}

bool line_swept_prefix_station_exceeds_exact(
    const Stock2& stock,
    const EPoint& motion_start,
    const EPoint& motion_end,
    const FT& parameter,
    const FT& tool_radius,
    const FT& cap_chord_ratio,
    const FT& gap_close_ratio)
{
    const LinePrefixStation2 station =
        line_prefix_station(motion_start, motion_end, parameter);
    const AuditExactStationClassification2 base =
        classify_audit_unguarded_station_exact(
            stock, station.center, tool_radius, cap_chord_ratio, gap_close_ratio);
    // CLIP MONOTONICITY (header): a static WITHIN_CAP is the in-process verdict.
    if (base.disposition() != AuditExactStationDisposition2::CAP_EXCEEDED) {
        return false;
    }
    if (CGAL::is_zero(parameter)) return true;
    return clip_line_prefix(
               base, station, tool_radius, cap_chord_ratio, gap_close_ratio)
               .disposition()
        == AuditExactStationDisposition2::CAP_EXCEEDED;
}

namespace {

// The static classification clipped to the rim the arc prefix has not covered.
AuditExactStationClassification2 clip_arc_prefix(
    const AuditExactStationClassification2& base,
    const EPoint& guide_center,
    const EPoint& motion_start,
    const EPoint& station_center,
    bool clockwise,
    ArcPrefixExtent2 extent,
    const FT& tool_radius,
    const FT& cap_chord_ratio,
    const FT& gap_close_ratio)
{
    const FT current_x = station_center.x() - guide_center.x();
    const FT current_y = station_center.y() - guide_center.y();
    const FT guide_radius_squared = current_x * current_x + current_y * current_y;
    std::vector<GpsPoint> events{
        radial_rim_point(
            station_center, current_x, current_y, tool_radius, true),
        radial_rim_point(
            station_center, current_x, current_y, tool_radius, false),
    };
    const std::vector<GpsPoint> start_disk_events =
        equal_circle_intersections(station_center, motion_start, tool_radius);
    events.insert(events.end(), start_disk_events.begin(), start_disk_events.end());
    FT sweep_start_x = motion_start.x() - guide_center.x();
    FT sweep_start_y = motion_start.y() - guide_center.y();
    FT sweep_end_x = current_x;
    FT sweep_end_y = current_y;
    if (clockwise) {
        std::swap(sweep_start_x, sweep_end_x);
        std::swap(sweep_start_y, sweep_end_y);
    }
    const auto covered = [&](const GpsPoint& point) {
        const CoordNT start_dx = point.x() - CoordNT(motion_start.x());
        const CoordNT start_dy = point.y() - CoordNT(motion_start.y());
        if (CGAL::compare(
                start_dx * start_dx + start_dy * start_dy,
                CoordNT(tool_radius * tool_radius))
            != CGAL::LARGER) {
            return true;
        }
        if (extent == ArcPrefixExtent2::FULL_TURN) return true;
        const CoordNT vector_x = point.x() - CoordNT(guide_center.x());
        const CoordNT vector_y = point.y() - CoordNT(guide_center.y());
        if (CGAL::is_zero(vector_x) && CGAL::is_zero(vector_y)) {
            return CGAL::compare(
                       guide_radius_squared,
                       tool_radius * tool_radius)
                != CGAL::LARGER;
        }
        const CGAL::Sign after_start = CGAL::sign(
            CoordNT(sweep_start_x) * vector_y
            - CoordNT(sweep_start_y) * vector_x);
        const CGAL::Sign before_end = CGAL::sign(
            vector_x * CoordNT(sweep_end_y)
            - vector_y * CoordNT(sweep_end_x));
        if (extent == ArcPrefixExtent2::AT_MOST_HALF_TURN) {
            return after_start != CGAL::NEGATIVE
                && before_end != CGAL::NEGATIVE;
        }
        return after_start != CGAL::NEGATIVE
            || before_end != CGAL::NEGATIVE;
    };
    return classify_audit_filtered_station_exact(
        station_center,
        tool_radius,
        cap_chord_ratio,
        gap_close_ratio,
        clip_rim_arcs(
            flatten_arcs(base),
            events,
            station_center,
            tool_radius,
            covered));
}

FT arc_station_guide_radius_squared(
    const EPoint& guide_center,
    const EPoint& station_center)
{
    const FT current_x = station_center.x() - guide_center.x();
    const FT current_y = station_center.y() - guide_center.y();
    const FT guide_radius_squared = current_x * current_x + current_y * current_y;
    if (CGAL::is_zero(guide_radius_squared)) {
        throw InvalidSweptPrefixInputError(
            "arc swept-prefix guide radius must be positive");
    }
    return guide_radius_squared;
}

} // namespace

AuditExactStationClassification2 clip_arc_swept_prefix_station_exact(
    const AuditExactStationClassification2& base,
    const EPoint& guide_center,
    const EPoint& motion_start,
    const EPoint& station_center,
    bool clockwise,
    ArcPrefixExtent2 extent,
    const FT& tool_radius,
    const FT& cap_chord_ratio,
    const FT& gap_close_ratio)
{
    arc_station_guide_radius_squared(guide_center, station_center);
    if (extent == ArcPrefixExtent2::ZERO) return base;
    return clip_arc_prefix(
        base, guide_center, motion_start, station_center, clockwise, extent,
        tool_radius, cap_chord_ratio, gap_close_ratio);
}

AuditExactStationClassification2 classify_arc_swept_prefix_station_exact(
    const Stock2& stock,
    const EPoint& guide_center,
    const EPoint& motion_start,
    const EPoint& station_center,
    bool clockwise,
    ArcPrefixExtent2 extent,
    const FT& tool_radius,
    const FT& cap_chord_ratio,
    const FT& gap_close_ratio)
{
    arc_station_guide_radius_squared(guide_center, station_center);
    const AuditExactStationClassification2 base =
        classify_audit_unguarded_station_exact(
            stock,
            station_center,
            tool_radius,
            cap_chord_ratio,
            gap_close_ratio);
    if (extent == ArcPrefixExtent2::ZERO) return base;
    return clip_arc_prefix(
        base, guide_center, motion_start, station_center, clockwise, extent,
        tool_radius, cap_chord_ratio, gap_close_ratio);
}

bool arc_swept_prefix_station_exceeds_exact(
    const Stock2& stock,
    const EPoint& guide_center,
    const EPoint& motion_start,
    const EPoint& station_center,
    bool clockwise,
    ArcPrefixExtent2 extent,
    const FT& tool_radius,
    const FT& cap_chord_ratio,
    const FT& gap_close_ratio)
{
    arc_station_guide_radius_squared(guide_center, station_center);
    const AuditExactStationClassification2 base =
        classify_audit_unguarded_station_exact(
            stock,
            station_center,
            tool_radius,
            cap_chord_ratio,
            gap_close_ratio);
    // CLIP MONOTONICITY (header): a static WITHIN_CAP is the in-process verdict.
    if (base.disposition() != AuditExactStationDisposition2::CAP_EXCEEDED) {
        return false;
    }
    if (extent == ArcPrefixExtent2::ZERO) return true;
    return clip_arc_prefix(
               base, guide_center, motion_start, station_center, clockwise,
               extent, tool_radius, cap_chord_ratio, gap_close_ratio)
               .disposition()
        == AuditExactStationDisposition2::CAP_EXCEEDED;
}
