#include "coverage_certificate_2/target_decision_2.h"

#include "coverage_certificate_2/errors_2.h"

namespace compas_cgal::coverage_certificate {

namespace {
RectSide rect_side(const RectRegion& rect, const LineTable& lines, const ExactPoint& w, const Interval& wx,
                   const Interval& wy, const FilterMode mode, DecisionCounters& counters)
{
    bool boundary = false;
    for (std::size_t s = 0; s < 4; ++s) {
        const CGAL::Sign side = line_side(lines[rect.lines[s]], w, wx, wy, mode, counters);
        if (side == CGAL::ZERO) {
            boundary = true;
        } else if (side != rect.inward[s]) {
            return RectSide::Exterior;
        }
    }
    return boundary ? RectSide::Boundary : RectSide::Interior;
}

/// The VERTEX_DISK side at w: NEGATIVE inside, ZERO on, POSITIVE outside;
/// std::nullopt when w is one-root and the enclosure cannot decide (the
/// cross-field sign is not implemented: VERTEX_DISK pieces are a later phase).
std::optional<CGAL::Sign> vertex_disk_side(const VertexDisk& disk, const ExactPoint& w, const Interval& wx,
                                           const Interval& wy, const MachinableTarget& target)
{
    const Interval dx = wx - disk.x_i, dy = wy - disk.y_i, d2 = dx * dx + dy * dy;
    if (d2.inf() > target.tool_radius_squared_i.sup()) return CGAL::POSITIVE;
    if (d2.sup() < target.tool_radius_squared_i.inf()) return CGAL::NEGATIVE;
    if (w.x.is_extended() || w.y.is_extended()) return std::nullopt;
    return vertex_disk_sign(disk, Point_2(w.x.a0(), w.y.a0()), target.tool_radius_squared);
}
}  // namespace

TargetSide decide_checked_side(const TargetModel& model, const TargetCandidates& candidates,
                               const CheckedSupport& support, const ExactPoint& w, const FilterMode mode,
                               DecisionCounters& counters)
{
    ++counters.target_decisions;
    const Interval wx = interval_of(w.x), wy = interval_of(w.y);
    const LineTable& lines = *model.lines;
    const MachinableTarget& target = *model.target;
    std::vector<std::uint32_t> boundary_rects;
    std::vector<std::uint32_t> boundary_circles;
    // 1. Interiors.
    for (const std::uint32_t r : candidates.rects) {
        const RectSide side = rect_side((*model.rects)[r], lines, w, wx, wy, mode, counters);
        if (side == RectSide::Interior) return TargetSide::Inside;
        if (side == RectSide::Boundary) boundary_rects.push_back(r);
    }
    for (const std::uint32_t c : candidates.circles) {
        const BoundaryCircle& T = model.set->circles[c];
        const auto cmp = compare_to_circle(w, wx, wy, T.circle.center(), T.centre_x, T.centre_y, T.circle.squared_radius(),
                                           T.radius_squared_i, mode, counters);
        if (cmp == CGAL::SMALLER) return TargetSide::Inside;
        if (cmp == CGAL::EQUAL) boundary_circles.push_back(c);
    }
    bool undecided = false;
    for (const VertexDisk& disk : target.vertex_disks) {
        const auto side = vertex_disk_side(disk, w, wx, wy, target);
        if (!side || *side == CGAL::ZERO) {
            undecided = true;
        } else if (*side == CGAL::NEGATIVE) {
            return TargetSide::Inside;
        }
    }
    if (in_closed_domain(target.domain, target.vertex_intervals, w, mode, counters)) return TargetSide::Inside;
    // 2. Boundaries along shared supports (P5/P6).
    for (const std::uint32_t r : boundary_rects) {
        const RectRegion& R = (*model.rects)[r];
        for (std::size_t s = 0; s < 4; ++s) {
            if (line_side(lines[R.lines[s]], w) != CGAL::ZERO) continue;
            const int co = !support.is_line ? 0
                           : R.lines[s] == support.id ? 1
                                                      : line_coincidence(lines[support.id], lines[R.lines[s]]);
            if (co == 0) {
                throw CoverageCertificateNativeError("an exposed-edge witness lies on a target rectangle line the "
                                                     "arrangement did not cut");
            }
            const CGAL::Sign relative = co > 0 ? R.inward[s] : CGAL::opposite(R.inward[s]);
            if (relative == support.check_side) return TargetSide::Inside;
        }
    }
    for (const std::uint32_t c : boundary_circles) {
        const BoundaryCircle& T = model.set->circles[c];
        if (support.is_line) {
            throw CoverageCertificateNativeError("an exposed flank witness lies on a target circle the arrangement did not cut");
        }
        const BoundaryCircle& P = model.set->circles[support.id];
        if (!(P.circle.center() == T.circle.center()) || P.circle.squared_radius() != T.circle.squared_radius()) {
            throw CoverageCertificateNativeError("an exposed-edge witness lies on a target circle the arrangement did not cut");
        }
        // The vertex disk lies inside its circle: it covers the checked side iff
        // the piece checks its inside (an inner annulus circle).
        if (P.check_inside()) return TargetSide::Inside;
    }
    return undecided ? TargetSide::Undecided : TargetSide::Outside;
}

bool target_contains(const TargetModel& model, const TargetCandidates& candidates, const Point_2& u,
                     const FilterMode mode, DecisionCounters& counters)
{
    const ExactPoint w = ExactPoint::rational(u);
    const Interval wx = interval_of(u.x()), wy = interval_of(u.y());
    for (const std::uint32_t r : candidates.rects) {
        if (rect_side((*model.rects)[r], *model.lines, w, wx, wy, mode, counters) != RectSide::Exterior) return true;
    }
    for (const std::uint32_t c : candidates.circles) {
        const BoundaryCircle& T = model.set->circles[c];
        if (compare_to_circle(w, wx, wy, T.circle.center(), T.centre_x, T.centre_y, T.circle.squared_radius(),
                              T.radius_squared_i, mode, counters) != CGAL::LARGER) {
            return true;
        }
    }
    for (const VertexDisk& disk : model.target->vertex_disks) {
        ++counters.exact_target_membership;
        if (vertex_disk_sign(disk, u, model.target->tool_radius_squared) != CGAL::POSITIVE) return true;
    }
    return in_closed_domain(model.target->domain, model.target->vertex_intervals, w, mode, counters);
}

}  // namespace compas_cgal::coverage_certificate
