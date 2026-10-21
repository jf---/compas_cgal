#include "coverage_certificate_2/target_2.h"

#include "coverage_certificate_2/errors_2.h"
#include "coverage_certificate_2/lines_2.h"
#include "coverage_certificate_2/tower_2.h"

#include <CGAL/Polygon_2_algorithms.h>

#include <string>
#include <utility>

namespace compas_cgal::coverage_certificate {

namespace {
using Segment_2 = CK::Segment_2;

std::vector<Point_2> inject_loop(const std::vector<std::array<double, 2>>& loop, const std::string& name)
{
    if (loop.size() < 3) throw InvalidCentreDomainError(name + " needs at least three vertices");
    std::vector<Point_2> points;
    points.reserve(loop.size());
    for (const auto& xy : loop) points.emplace_back(inject_binary64(xy[0]), inject_binary64(xy[1]));
    for (std::size_t i = 0; i < points.size(); ++i) {
        if (points[i] == points[(i + 1) % points.size()]) {
            throw InvalidCentreDomainError(name + " has a zero-length edge at vertex " + std::to_string(i));
        }
    }
    return points;
}

Segment_2 edge(const std::vector<Point_2>& loop, const std::size_t i)
{
    return Segment_2(loop[i], loop[(i + 1) % loop.size()]);
}

bool loops_touch(const std::vector<Point_2>& a, const std::vector<Point_2>& b)
{
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < b.size(); ++j) {
            if (CGAL::do_intersect(edge(a, i), edge(b, j))) return true;
        }
    }
    return false;
}

/// Exact closed membership of q in the polygon with holes, by the crossing
/// number of a ray to +x over every loop edge (holes are inside the outer loop
/// and disjoint, so even-odd parity is membership of the interior). Every
/// comparison and orientation is filtered by certified intervals and decided
/// exactly when the enclosure straddles zero.
bool in_closed_polygon(const std::vector<std::vector<Point_2>>& loops,
                       const std::vector<std::vector<std::array<Interval, 2>>>& shadows, const Point_2& q,
                       const FilterMode mode, DecisionCounters& counters)
{
    const Interval qx = interval_of(q.x()), qy = interval_of(q.y());
    const auto compare_coordinate = [&](const FT& a, const Interval& ai, const FT& b, const Interval& bi) {
        if (mode == FilterMode::Lazy) {
            if (ai.sup() < bi.inf()) return CGAL::SMALLER;
            if (ai.inf() > bi.sup()) return CGAL::LARGER;
        }
        ++counters.exact_target_membership;
        return CGAL::compare(a, b);
    };
    bool inside = false;
    for (std::size_t l = 0; l < loops.size(); ++l) {
        const auto& loop = loops[l];
        const auto& shadow = shadows[l];
        for (std::size_t i = 0; i < loop.size(); ++i) {
            const std::size_t j = (i + 1) % loop.size();
            const Point_2& a = loop[i];
            const Point_2& b = loop[j];
            const Interval &ax = shadow[i][0], &ay = shadow[i][1], &bx = shadow[j][0], &by = shadow[j][1];
            const auto ca = compare_coordinate(a.y(), ay, q.y(), qy);
            const auto cb = compare_coordinate(b.y(), by, q.y(), qy);
            if (ca == cb && ca != CGAL::EQUAL) continue;  // strictly on one side of the ray's line
            CGAL::Orientation o = CGAL::COLLINEAR;
            bool decided = false;
            if (mode == FilterMode::Lazy) {
                const Interval det = (bx - ax) * (qy - ay) - (by - ay) * (qx - ax);
                if (det.inf() > 0) { o = CGAL::LEFT_TURN; decided = true; }
                if (det.sup() < 0) { o = CGAL::RIGHT_TURN; decided = true; }
            }
            if (!decided) {
                ++counters.exact_target_membership;
                o = CGAL::orientation(a, b, q);
            }
            if (o == CGAL::COLLINEAR) {
                // q on the supporting line: on the edge iff within its x-range
                // (a straddling collinear edge always contains q).
                const FT& lo = CGAL::compare(a.x(), b.x()) == CGAL::SMALLER ? a.x() : b.x();
                const FT& hi = CGAL::compare(a.x(), b.x()) == CGAL::SMALLER ? b.x() : a.x();
                if (CGAL::compare(lo, q.x()) != CGAL::LARGER && CGAL::compare(q.x(), hi) != CGAL::LARGER) return true;
                continue;
            }
            const bool a_above = ca == CGAL::LARGER;
            const bool b_above = cb == CGAL::LARGER;
            if (a_above == b_above) continue;
            // The edge meets the line y = q.y right of q iff q is left of the
            // upward-directed edge.
            if ((b_above && o == CGAL::LEFT_TURN) || (a_above && o == CGAL::RIGHT_TURN)) inside = !inside;
        }
    }
    return inside;
}

}  // namespace

/// sign((qx - vx)^2 + (qy - vy)^2 - r^2) for a one-root centre v.
///
/// With vx = x0 + x1*sqrt(d), vy = y0 + y1*sqrt(d) and ex = qx - x0, ey = qy - y0
/// the expression is A + B*sqrt(d), A = ex^2 + ey^2 + (x1^2 + y1^2) d - r^2 and
/// B = -2 (ex x1 + ey y1); all rational. A zero radicand or B is folded away
/// with rational compares before any extension is built (exactness checklist
/// question 17).
CGAL::Sign vertex_disk_sign(const VertexDisk& disk, const Point_2& q, const FT& r_sq)
{
    const FT d = disk.x.is_extended() ? disk.x.root() : (disk.y.is_extended() ? disk.y.root() : FT(0));
    const FT x1 = disk.x.is_extended() ? disk.x.a1() : FT(0);
    const FT y1 = disk.y.is_extended() ? disk.y.a1() : FT(0);
    const FT ex = q.x() - disk.x.a0();
    const FT ey = q.y() - disk.y.a0();
    const FT a = ex * ex + ey * ey + (x1 * x1 + y1 * y1) * d - r_sq;
    const FT b = FT(-2) * (ex * x1 + ey * y1);
    if (CGAL::sign(d) == CGAL::ZERO || CGAL::sign(b) == CGAL::ZERO) return CGAL::sign(a);
    return CGAL::sign(OneRoot(a, b, d));
}

bool in_closed_domain(const CentreDomain& domain, const std::vector<std::vector<std::array<Interval, 2>>>& shadows,
                      const ExactPoint& q, const FilterMode mode, DecisionCounters& counters)
{
    const Interval qx = interval_of(q.x), qy = interval_of(q.y);
    const auto compare_coordinate = [&](const FT& a, const Interval& ai, const OneRoot& b, const Interval& bi) {
        if (mode == FilterMode::Lazy) {
            if (ai.sup() < bi.inf()) return CGAL::SMALLER;
            if (ai.inf() > bi.sup()) return CGAL::LARGER;
        }
        ++counters.exact_target_membership;
        return static_cast<CGAL::Comparison_result>(CGAL::opposite(sign_one_root(b - a)));
    };
    bool inside = false;
    for (std::size_t l = 0; l < domain.loops.size(); ++l) {
        const auto& loop = domain.loops[l];
        const auto& shadow = shadows[l];
        for (std::size_t i = 0; i < loop.size(); ++i) {
            const std::size_t j = (i + 1) % loop.size();
            const Point_2& a = loop[i];
            const Point_2& b = loop[j];
            const Interval &ax = shadow[i][0], &ay = shadow[i][1], &bx = shadow[j][0], &by = shadow[j][1];
            const auto ca = compare_coordinate(a.y(), ay, q.y, qy);
            const auto cb = compare_coordinate(b.y(), by, q.y, qy);
            if (ca == cb && ca != CGAL::EQUAL) continue;
            CGAL::Sign o = CGAL::ZERO;
            bool decided = false;
            if (mode == FilterMode::Lazy) {
                const Interval det = (bx - ax) * (qy - ay) - (by - ay) * (qx - ax);
                if (det.inf() > 0) { o = CGAL::POSITIVE; decided = true; }
                if (det.sup() < 0) { o = CGAL::NEGATIVE; decided = true; }
            }
            if (!decided) {
                ++counters.exact_target_membership;
                o = sign_one_root((b.x() - a.x()) * (q.y - a.y()) - (b.y() - a.y()) * (q.x - a.x()));
            }
            if (o == CGAL::ZERO) {
                const bool a_left = CGAL::compare(a.x(), b.x()) == CGAL::SMALLER;
                const FT& lo = a_left ? a.x() : b.x();
                const FT& hi = a_left ? b.x() : a.x();
                if (sign_one_root(q.x - lo) != CGAL::NEGATIVE && sign_one_root(hi - q.x) != CGAL::NEGATIVE) {
                    // Collinear and within the x-range: on the edge unless the edge
                    // is vertical, where the y-range decides.
                    if (lo != hi) return true;
                    const bool a_low = CGAL::compare(a.y(), b.y()) == CGAL::SMALLER;
                    const FT& ylo = a_low ? a.y() : b.y();
                    const FT& yhi = a_low ? b.y() : a.y();
                    if (sign_one_root(q.y - ylo) != CGAL::NEGATIVE && sign_one_root(yhi - q.y) != CGAL::NEGATIVE) return true;
                }
                continue;
            }
            const bool a_above = ca == CGAL::LARGER;
            const bool b_above = cb == CGAL::LARGER;
            if (a_above == b_above) continue;
            if ((b_above && o == CGAL::POSITIVE) || (a_above && o == CGAL::NEGATIVE)) inside = !inside;
        }
    }
    return inside;
}

CentreDomain CentreDomain::build(const std::vector<std::array<double, 2>>& outer,
                                 const std::vector<std::vector<std::array<double, 2>>>& holes)
{
    CentreDomain domain;
    domain.loops.push_back(inject_loop(outer, "outer loop"));
    for (std::size_t h = 0; h < holes.size(); ++h) domain.loops.push_back(inject_loop(holes[h], "hole " + std::to_string(h)));
    const CK kernel;
    for (std::size_t l = 0; l < domain.loops.size(); ++l) {
        const auto& loop = domain.loops[l];
        const std::string name = l == 0 ? "outer loop" : "hole " + std::to_string(l - 1);
        if (!CGAL::is_simple_2(loop.begin(), loop.end(), kernel)) throw InvalidCentreDomainError(name + " is not simple");
        const auto orientation = CGAL::orientation_2(loop.begin(), loop.end(), kernel);
        if (l == 0 && orientation != CGAL::COUNTERCLOCKWISE) throw InvalidCentreDomainError("outer loop must be CCW");
        if (l > 0 && orientation != CGAL::CLOCKWISE) throw InvalidCentreDomainError(name + " must be CW");
    }
    for (std::size_t h = 1; h < domain.loops.size(); ++h) {
        const std::string name = "hole " + std::to_string(h - 1);
        if (loops_touch(domain.loops[0], domain.loops[h])) throw InvalidCentreDomainError(name + " touches the outer loop");
        if (CGAL::bounded_side_2(domain.loops[0].begin(), domain.loops[0].end(), domain.loops[h][0], kernel) !=
            CGAL::ON_BOUNDED_SIDE) {
            throw InvalidCentreDomainError(name + " is not strictly inside the outer loop");
        }
        for (std::size_t g = 1; g < h; ++g) {
            const auto& a = domain.loops[g];
            const auto& b = domain.loops[h];
            if (loops_touch(a, b) ||
                CGAL::bounded_side_2(a.begin(), a.end(), b[0], kernel) != CGAL::ON_UNBOUNDED_SIDE ||
                CGAL::bounded_side_2(b.begin(), b.end(), a[0], kernel) != CGAL::ON_UNBOUNDED_SIDE) {
                throw InvalidCentreDomainError("holes " + std::to_string(g - 1) + " and " + name + " overlap");
            }
        }
    }
    return domain;
}

VertexDisk VertexDisk::build(const FT& x0, const FT& x1, const FT& y0, const FT& y1, const FT& delta)
{
    if (CGAL::sign(delta) == CGAL::NEGATIVE) throw InvalidCertificateInputError("vertex-disk radicand must be >= 0");
    // Exactness checklist q17: a zero radicand, a zero coefficient or a rational
    // square radicand is folded with rational compares before any extension is
    // built, so an extended coordinate always carries an irrational sqrt(delta).
    const OneRoot x = make_reduced_one_root(x0, x1, delta);
    const OneRoot y = make_reduced_one_root(y0, y1, delta);
    return {x, y, interval_of(x), interval_of(y)};
}

MachinableTarget MachinableTarget::build(CentreDomain domain, const double tool_radius, std::vector<VertexDisk> vertex_disks)
{
    const FT r = inject_binary64(tool_radius);
    if (CGAL::sign(r) != CGAL::POSITIVE) throw InvalidCertificateInputError("tool radius must be positive");
    MachinableTarget target{std::move(domain), r * r, interval_of(r * r), std::move(vertex_disks), {}, {}, r};
    for (const auto& loop : target.domain.loops) {
        auto& shadow = target.vertex_intervals.emplace_back();
        for (const auto& p : loop) shadow.push_back({interval_of(p.x()), interval_of(p.y())});
        for (std::size_t i = 0; i < loop.size(); ++i) {
            const auto& a = loop[i];
            const auto& b = loop[(i + 1) % loop.size()];
            const Interval ax = interval_of(a.x()), ay = interval_of(a.y()), bx = interval_of(b.x()), by = interval_of(b.y());
            target.edge_boxes.push_back({Interval(std::min(ax.inf(), bx.inf())), Interval(std::max(ax.sup(), bx.sup())),
                                         Interval(std::min(ay.inf(), by.inf())), Interval(std::max(ay.sup(), by.sup()))});
        }
    }
    return target;
}

bool MachinableTarget::contains(const Point_2& q, const FilterMode mode, DecisionCounters& counters) const
{
    // POLYGON member.
    if (in_closed_polygon(domain.loops, vertex_intervals, q, mode, counters)) return true;
    // EDGE_CAPSULE members: certified box-distance filter, exact on overlap.
    const Interval qx = interval_of(q.x()), qy = interval_of(q.y());
    std::size_t k = 0;
    for (const auto& loop : domain.loops) {
        for (std::size_t i = 0; i < loop.size(); ++i, ++k) {
            if (mode == FilterMode::Lazy) {
                const auto& box = edge_boxes[k];
                const Interval gx = Interval(std::max({0.0, (box[0] - qx).inf(), (qx - box[1]).inf()}));
                const Interval gy = Interval(std::max({0.0, (box[2] - qy).inf(), (qy - box[3]).inf()}));
                // gx, gy are certified lower bounds of the box gaps (outward-rounded inf).
                if ((gx * gx + gy * gy).inf() > tool_radius_squared_i.sup()) continue;
            }
            ++counters.exact_target_membership;
            if (CGAL::compare(CGAL::squared_distance(q, edge(loop, i)), tool_radius_squared) != CGAL::LARGER) return true;
        }
    }
    // VERTEX_DISK members (one-root centres).
    for (const auto& disk : vertex_disks) {
        if (mode == FilterMode::Lazy) {
            const Interval dx = qx - disk.x_i, dy = qy - disk.y_i;
            const Interval d2 = dx * dx + dy * dy;
            if (d2.inf() > tool_radius_squared_i.sup()) continue;
            if (d2.sup() < tool_radius_squared_i.inf()) return true;
        }
        ++counters.exact_target_membership;
        if (vertex_disk_sign(disk, q, tool_radius_squared) != CGAL::POSITIVE) return true;
    }
    return false;
}

}  // namespace compas_cgal::coverage_certificate
