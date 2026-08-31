#include "audit_certification_2.h"
#include "engagement_2.h"

#include <CGAL/Arr_trapezoid_ric_point_location.h>
#include <CGAL/Arrangement_zone_2.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>

#include <boost/variant.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <utility>
#include <vector>

namespace {

using FT = Epeck::FT;
using CoordNT = GpsPoint::CoordNT;

CGAL::Sign sign_mixed_radical_impl(
    const FT& a,
    const FT& b,
    const FT& c,
    const FT& d,
    const FT& alpha,
    const FT& beta)
{
    const bool alpha_extended = !CGAL::is_zero(alpha);
    const bool beta_extended = !CGAL::is_zero(beta);
    if (!alpha_extended && !beta_extended) {
        return CGAL::sign(a);
    }
    if (!beta_extended) {
        return CGAL::sign(CoordNT(a, b, alpha));
    }
    if (!alpha_extended) {
        return CGAL::sign(CoordNT(a, c, beta));
    }
    if (alpha == beta) {
        return CGAL::sign(CoordNT(a + d * alpha, b + c, alpha));
    }
    const CoordNT u(a, b, alpha);
    const CoordNT w(c, d, alpha);
    const CGAL::Sign u_sign = CGAL::sign(u);
    const CGAL::Sign w_sign = CGAL::sign(w);
    if (w_sign == CGAL::ZERO) {
        return u_sign;
    }
    if (u_sign == CGAL::ZERO || u_sign == w_sign) {
        return w_sign;
    }
    const CGAL::Comparison_result magnitude =
        CGAL::compare(u * u, w * w * CoordNT(beta));
    if (magnitude == CGAL::LARGER) {
        return u_sign;
    }
    if (magnitude == CGAL::SMALLER) {
        return w_sign;
    }
    return CGAL::ZERO;
}

struct RadPoint2 {
    FT x0;
    FT x1;
    FT y0;
    FT y1;
    FT root;
};

RadPoint2 radial_point(const GpsPoint& point)
{
    const CoordNT& x = point.x();
    const CoordNT& y = point.y();
    const ExactOneRootPoint2 decomposition =
        decompose_exact_one_root_point(point);
    return {
        x.a0(),
        x.is_extended() ? x.a1() : FT(0),
        y.a0(),
        y.is_extended() ? y.a1() : FT(0),
        decomposition.root(),
    };
}

bool run_exceeds_cap(
    const GpsPoint& start,
    const GpsPoint& end,
    const EPoint& center,
    const FT& cap_chord_ratio,
    const FT& squared_threshold)
{
    if (start == end) {
        return true;
    }
    const RadPoint2 p = radial_point(start);
    const RadPoint2 q = radial_point(end);
    const FT px = p.x0 - center.x();
    const FT py = p.y0 - center.y();
    const FT qx = q.x0 - center.x();
    const FT qy = q.y0 - center.y();
    const CGAL::Sign orientation = sign_mixed_radical_impl(
        px * qy - py * qx,
        p.x1 * qy - p.y1 * qx,
        px * q.y1 - py * q.x1,
        p.x1 * q.y1 - p.y1 * q.x1,
        p.root,
        q.root);
    if (orientation == CGAL::NEGATIVE) {
        return true;
    }
    if (orientation == CGAL::ZERO) {
        return CGAL::compare(cap_chord_ratio, FT(4)) == CGAL::SMALLER;
    }

    const FT dx = q.x0 - p.x0;
    const FT dy = q.y0 - p.y0;
    const FT a = dx * dx + dy * dy
        + (p.x1 * p.x1 + p.y1 * p.y1) * p.root
        + (q.x1 * q.x1 + q.y1 * q.y1) * q.root;
    const FT b = -FT(2) * (dx * p.x1 + dy * p.y1);
    const FT c = FT(2) * (dx * q.x1 + dy * q.y1);
    const FT d = -FT(2) * (p.x1 * q.x1 + p.y1 * q.y1);
    return sign_mixed_radical_impl(
               a - squared_threshold,
               b,
               c,
               d,
               p.root,
               q.root)
        == CGAL::POSITIVE;
}

struct ExactRimArc2 {
    GpsPoint start;
    GpsPoint end;
    std::vector<std::pair<GpsPoint, GpsPoint>> reporting_arcs;
};

using Arrangement = Gps::Arrangement_2;
// Regularized stock can place a cutter endpoint in a face with multiple outer
// CCBs. The walk and naive locators assume one outer CCB in paths used here;
// trapezoidal point location supports the valid topology and seeds the same
// exact arrangement-zone traversal without a floating-point decision.
using PointLocation = CGAL::Arr_trapezoid_ric_point_location<Arrangement>;

struct ExactEngagementVisitor2 {
    using X_monotone_curve_2 = Arrangement::X_monotone_curve_2;
    using Vertex_handle = Arrangement::Vertex_handle;
    using Halfedge_handle = Arrangement::Halfedge_handle;
    using Face_handle = Arrangement::Face_handle;
    using Result = std::pair<Halfedge_handle, bool>;

    std::vector<X_monotone_curve_2> engaged;

    void init(Arrangement*) {}
    Result found_subcurve(
        const X_monotone_curve_2& curve,
        Face_handle face,
        Vertex_handle,
        Halfedge_handle,
        Vertex_handle,
        Halfedge_handle)
    {
        if (face->contained()) {
            engaged.push_back(curve);
        }
        return Result(Halfedge_handle(), false);
    }
    Result found_overlap(
        const X_monotone_curve_2&,
        Halfedge_handle,
        Vertex_handle,
        Vertex_handle)
    {
        return Result(Halfedge_handle(), false);
    }
};

std::vector<ExactRimArc2> exact_engaged_arcs(
    const Stock2& stock,
    const EPoint& center,
    const FT& tool_radius,
    std::vector<GpsPoint>& intersections)
{
    Arrangement& arrangement =
        const_cast<Gps&>(stock.set()).arrangement();
    PointLocation point_location(arrangement);
    const GpsTraits::Curve_2 cutter(
        ECircle(center, tool_radius * tool_radius));
    GpsTraits traits;
    std::vector<boost::variant<GpsPoint, GpsXCurve>> pieces;
    traits.make_x_monotone_2_object()(cutter, std::back_inserter(pieces));

    ExactEngagementVisitor2 visitor;
    for (const auto& piece : pieces) {
        if (const GpsXCurve* curve = boost::get<GpsXCurve>(&piece)) {
            CGAL::Arrangement_zone_2<
                Arrangement,
                ExactEngagementVisitor2> zone(arrangement, &visitor);
            zone.init(*curve, point_location);
            zone.compute_zone();
        }
    }

    std::vector<ExactRimArc2> arcs;
    for (const GpsXCurve& curve : visitor.engaged) {
        GpsPoint start = curve.source();
        GpsPoint end = curve.target();
        if (start == end) {
            continue;
        }
        if (curve.orientation() == CGAL::CLOCKWISE) {
            std::swap(start, end);
        }
        intersections.push_back(start);
        intersections.push_back(end);
        arcs.push_back({
            start,
            end,
            {{std::move(start), std::move(end)}},
        });
    }
    std::sort(
        intersections.begin(),
        intersections.end(),
        [](const GpsPoint& lhs, const GpsPoint& rhs) {
            return GpsTraits().compare_xy_2_object()(lhs, rhs)
                == CGAL::SMALLER;
        });
    intersections.erase(
        std::unique(intersections.begin(), intersections.end()),
        intersections.end());
    return arcs;
}

std::vector<ExactRimArc2> maximal_exact_runs(
    std::vector<ExactRimArc2> arcs,
    const EPoint& center)
{
    const CoordNT exact_x(center.x());
    const CoordNT exact_y(center.y());
    const auto upper = [&exact_x, &exact_y](const GpsPoint& point) {
        const CGAL::Comparison_result y = CGAL::compare(point.y(), exact_y);
        if (y == CGAL::LARGER) {
            return true;
        }
        if (y == CGAL::SMALLER) {
            return false;
        }
        return CGAL::compare(point.x(), exact_x) == CGAL::LARGER;
    };
    std::sort(
        arcs.begin(),
        arcs.end(),
        [&upper](const ExactRimArc2& lhs, const ExactRimArc2& rhs) {
            const bool lhs_upper = upper(lhs.start);
            const bool rhs_upper = upper(rhs.start);
            if (lhs_upper != rhs_upper) {
                return lhs_upper;
            }
            const CGAL::Comparison_result x =
                CGAL::compare(lhs.start.x(), rhs.start.x());
            return lhs_upper ? x == CGAL::LARGER : x == CGAL::SMALLER;
        });
    std::vector<ExactRimArc2> runs;
    for (const ExactRimArc2& arc : arcs) {
        if (!runs.empty() && runs.back().end == arc.start) {
            runs.back().end = arc.end;
            runs.back().reporting_arcs.insert(
                runs.back().reporting_arcs.end(),
                arc.reporting_arcs.begin(),
                arc.reporting_arcs.end());
        } else {
            runs.push_back(arc);
        }
    }
    if (runs.size() > 1 && runs.back().end == runs.front().start) {
        runs.front().start = runs.back().start;
        runs.front().reporting_arcs.insert(
            runs.front().reporting_arcs.begin(),
            runs.back().reporting_arcs.begin(),
            runs.back().reporting_arcs.end());
        runs.pop_back();
    }
    return runs;
}

std::vector<ExactRimArc2> pessimistic_exact_runs(
    const std::vector<ExactRimArc2>& runs,
    const EPoint& center,
    const FT& gap_close_ratio,
    const FT& squared_gap_threshold)
{
    const std::size_t count = runs.size();
    if (count == 0) {
        return {};
    }
    if (count == 1 && runs.front().start == runs.front().end) {
        return runs;
    }
    std::vector<bool> gap_closed(count);
    std::size_t open_count = 0;
    for (std::size_t index = 0; index < count; ++index) {
        gap_closed[index] = !run_exceeds_cap(
            runs[index].end,
            runs[(index + 1) % count].start,
            center,
            gap_close_ratio,
            squared_gap_threshold);
        if (!gap_closed[index]) {
            ++open_count;
        }
    }
    if (open_count == 0) {
        return {{runs.front().start, runs.front().start, {}}};
    }

    std::size_t first = 0;
    while (gap_closed[(first + count - 1) % count]) {
        ++first;
    }
    std::vector<ExactRimArc2> pessimistic;
    GpsPoint start = runs[first].start;
    GpsPoint end = runs[first].end;
    for (std::size_t step = 1; step < count; ++step) {
        const std::size_t index = (first + step) % count;
        if (gap_closed[(index + count - 1) % count]) {
            end = runs[index].end;
        } else {
            pessimistic.push_back({start, end, {}});
            start = runs[index].start;
            end = runs[index].end;
        }
    }
    pessimistic.push_back({start, end, {}});
    return pessimistic;
}

} // namespace

CGAL::Sign audit_sign_mixed_radical_exact(
    const Epeck::FT& a,
    const Epeck::FT& b,
    const Epeck::FT& c,
    const Epeck::FT& d,
    const Epeck::FT& alpha,
    const Epeck::FT& beta)
{
    if (CGAL::sign(alpha) == CGAL::NEGATIVE
        || CGAL::sign(beta) == CGAL::NEGATIVE) {
        throw InvalidMixedRadicalRootError(
            "mixed-radical roots must be nonnegative");
    }
    return sign_mixed_radical_impl(a, b, c, d, alpha, beta);
}

void validate_engagement_input_binary64(
    double cx,
    double cy,
    double tool_radius,
    double cap_chord_ratio,
    double gap_close_ratio)
{
    if (!std::isfinite(cx) || !std::isfinite(cy)
        || !std::isfinite(tool_radius)
        || !std::isfinite(cap_chord_ratio)
        || !std::isfinite(gap_close_ratio)) {
        throw NonFiniteEngagementInputError(
            "engagement station inputs must be finite binary64");
    }
    if (!(tool_radius > 0.0)) {
        throw NonPositiveEngagementToolRadiusError(
            "engagement tool radius must be positive");
    }
    if (!(cap_chord_ratio > 0.0 && cap_chord_ratio <= 4.0)) {
        throw InvalidEngagementCapRatioError(
            "engagement cap chord ratio must lie in (0, 4]");
    }
    if (!(gap_close_ratio >= 0.0 && gap_close_ratio <= 4.0)) {
        throw InvalidEngagementGapRatioError(
            "engagement gap chord ratio must lie in [0, 4]");
    }
}

void validate_segment_certification_input_binary64(
    double x0,
    double y0,
    double x1,
    double y1,
    double tool_radius,
    double cap_radians)
{
    if (!std::isfinite(x0) || !std::isfinite(y0)
        || !std::isfinite(x1) || !std::isfinite(y1)
        || !std::isfinite(tool_radius)
        || !std::isfinite(cap_radians)) {
        throw NonFiniteSegmentCertificationInputError(
            "segment certification inputs must be finite binary64");
    }
    if (!(tool_radius > 0.0)) {
        throw NonPositiveSegmentCertificationToolRadiusError(
            "segment certification tool radius must be positive");
    }
    if (!(cap_radians > 0.0 && cap_radians <= std::numbers::pi)) {
        throw InvalidSegmentCertificationCapError(
            "segment certification cap must lie in (0, pi]");
    }
}

void validate_mixed_radical_input_binary64(
    double a,
    double b,
    double c,
    double d,
    double alpha,
    double beta)
{
    if (!std::isfinite(a) || !std::isfinite(b)
        || !std::isfinite(c) || !std::isfinite(d)
        || !std::isfinite(alpha) || !std::isfinite(beta)) {
        throw NonFiniteMixedRadicalInputError(
            "mixed-radical inputs must be finite binary64");
    }
    if (alpha < 0.0 || beta < 0.0) {
        throw InvalidMixedRadicalRootError(
            "mixed-radical roots must be nonnegative");
    }
}

ExactOneRootPoint2::ExactOneRootPoint2(Epeck::FT root)
    : root_(std::move(root))
{
}

const Epeck::FT& ExactOneRootPoint2::root() const noexcept
{
    return root_;
}

ExactOneRootPoint2 decompose_exact_one_root_point(const GpsPoint& point)
{
    const CoordNT& x = point.x();
    const CoordNT& y = point.y();
    if (x.is_extended() && y.is_extended() && x.root() != y.root()) {
        throw ExactOneRootCoordinateMismatchError(
            "one arrangement point carries mismatched coordinate roots");
    }
    if (x.is_extended()) {
        return ExactOneRootPoint2(x.root());
    }
    if (y.is_extended()) {
        return ExactOneRootPoint2(y.root());
    }
    return ExactOneRootPoint2(Epeck::FT(0));
}

AuditExactStationClassification2::AuditExactStationClassification2(
    AuditExactStationDisposition2 disposition,
    std::vector<GpsPoint> boundary_intersections,
    std::vector<std::vector<std::pair<GpsPoint, GpsPoint>>> true_run_arcs)
    : disposition_(disposition),
      boundary_intersections_(std::move(boundary_intersections)),
      true_run_arcs_(std::move(true_run_arcs))
{
}

const std::vector<std::vector<std::pair<GpsPoint, GpsPoint>>>&
AuditExactStationClassification2::true_run_arcs() const noexcept
{
    return true_run_arcs_;
}

AuditExactStationDisposition2
AuditExactStationClassification2::disposition() const noexcept
{
    return disposition_;
}

const std::vector<GpsPoint>&
AuditExactStationClassification2::boundary_intersections() const noexcept
{
    return boundary_intersections_;
}

AuditExactStationClassification2 classify_audit_unguarded_station_exact(
    const Stock2& stock,
    const EPoint& center,
    const Epeck::FT& tool_radius,
    const Epeck::FT& cap_chord_ratio,
    const Epeck::FT& gap_close_ratio)
{
    if (CGAL::sign(tool_radius) != CGAL::POSITIVE) {
        throw AuditExactStationToolRadiusError(
            "exact station tool radius must be positive");
    }
    if (CGAL::sign(cap_chord_ratio) != CGAL::POSITIVE
        || CGAL::compare(cap_chord_ratio, Epeck::FT(4)) == CGAL::LARGER) {
        throw AuditExactStationCapRatioError(
            "exact station cap chord ratio must lie in (0, 4]");
    }
    if (CGAL::sign(gap_close_ratio) == CGAL::NEGATIVE
        || CGAL::compare(gap_close_ratio, Epeck::FT(4)) == CGAL::LARGER) {
        throw AuditExactStationGapRatioError(
            "exact station gap chord ratio must lie in [0, 4]");
    }
    std::vector<GpsPoint> intersections;
    const std::vector<ExactRimArc2> runs = maximal_exact_runs(
        exact_engaged_arcs(stock, center, tool_radius, intersections),
        center);
    std::vector<std::vector<std::pair<GpsPoint, GpsPoint>>> true_run_arcs;
    true_run_arcs.reserve(runs.size());
    for (const ExactRimArc2& run : runs) {
        true_run_arcs.push_back(run.reporting_arcs);
    }
    const Epeck::FT squared_threshold =
        cap_chord_ratio * tool_radius * tool_radius;
    const Epeck::FT squared_gap_threshold =
        gap_close_ratio * tool_radius * tool_radius;
    for (const ExactRimArc2& run : pessimistic_exact_runs(
             runs, center, gap_close_ratio, squared_gap_threshold)) {
        if (run_exceeds_cap(
                run.start,
                run.end,
                center,
                cap_chord_ratio,
                squared_threshold)) {
            return AuditExactStationClassification2(
                AuditExactStationDisposition2::CAP_EXCEEDED,
                std::move(intersections),
                std::move(true_run_arcs));
        }
    }
    return AuditExactStationClassification2(
        AuditExactStationDisposition2::WITHIN_CAP,
        std::move(intersections),
        std::move(true_run_arcs));
}

AuditExactStationDisposition2 replay_audit_unguarded_station_exact(
    const Stock2& stock,
    const EPoint& center,
    const Epeck::FT& tool_radius,
    const Epeck::FT& cap_chord_ratio)
{
    return classify_audit_unguarded_station_exact(
               stock,
               center,
               tool_radius,
               cap_chord_ratio)
        .disposition();
}
