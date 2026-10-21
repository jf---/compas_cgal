#include "coverage_certificate_2/lines_2.h"

#include "coverage_certificate_2/errors_2.h"

#include <algorithm>
#include <optional>
#include <string>
#include <utility>

namespace compas_cgal::coverage_certificate {


LineSupport LineSupport::build(const FT& nx, const FT& ny, const Point_2& through, const FT& k1, const FT& radicand)
{
    LineSupport line;
    line.nx = nx;
    line.ny = ny;
    line.nn = nx * nx + ny * ny;
    if (CGAL::sign(line.nn) == CGAL::ZERO) throw CoverageCertificateNativeError("line with a zero normal");
    const FT k0 = nx * through.x() + ny * through.y();
    std::optional<FT> folded;
    if (CGAL::sign(k1) == CGAL::ZERO || CGAL::sign(radicand) == CGAL::ZERO) {
        folded = FT(0);
    } else {
        folded = rational_sqrt(radicand);
        if (folded) *folded = k1 * *folded;
    }
    if (folded) {
        // Rational line: n.x = k0 + shift, base point through + (shift / |n|^2) n.
        const FT shift = *folded;
        line.k0 = k0 + shift;
        line.k1 = FT(0);
        line.L = FT(0);
        line.b0x = through.x() + shift / line.nn * nx;
        line.b0y = through.y() + shift / line.nn * ny;
        line.b1x = FT(0);
        line.b1y = FT(0);
    } else {
        line.k0 = k0;
        line.k1 = k1;
        line.L = radicand;
        line.b0x = through.x();
        line.b0y = through.y();
        line.b1x = k1 / line.nn * nx;
        line.b1y = k1 / line.nn * ny;
    }
    line.nx_i = interval_of(line.nx);
    line.ny_i = interval_of(line.ny);
    line.k0_i = interval_of(line.k0);
    line.k1_i = interval_of(line.k1);
    line.sqrt_l_i = CGAL::sqrt(interval_of(line.L));
    line.offset_i = line.k0_i + line.k1_i * line.sqrt_l_i;
    line.b0x_i = interval_of(line.b0x);
    line.b0y_i = interval_of(line.b0y);
    line.b1x_i = interval_of(line.b1x);
    line.b1y_i = interval_of(line.b1y);
    line.nn_i = interval_of(line.nn);
    return line;
}

CGAL::Sign line_side(const LineSupport& line, const ExactPoint& p)
{
    const OneRoot value = line.nx * p.x + line.ny * p.y - line.k0;
    if (line.is_rational()) return sign_one_root(value);
    return sign_tower({value, OneRoot(-line.k1), OneRoot(line.L)});
}

CGAL::Sign line_side(const LineSupport& line, const Point_2& q, const Interval& qx, const Interval& qy,
                     const FilterMode mode, DecisionCounters& counters)
{
    if (mode == FilterMode::Lazy) {
        const Interval v = line.nx_i * qx + line.ny_i * qy - line.offset_i;
        if (v.inf() > 0) return CGAL::POSITIVE;
        if (v.sup() < 0) return CGAL::NEGATIVE;
    }
    ++counters.exact_line_side;
    return line_side(line, ExactPoint::rational(q));
}

CGAL::Comparison_result compare_to_circle(const ExactPoint& p, const Interval& px, const Interval& py, const Point_2& c,
                                          const Interval& cx, const Interval& cy, const FT& r2, const Interval& r2i,
                                          const FilterMode mode, DecisionCounters& counters)
{
    if (mode == FilterMode::Lazy) {
        const Interval dx = px - cx, dy = py - cy, d2 = dx * dx + dy * dy;
        if (d2.sup() < r2i.inf()) return CGAL::SMALLER;
        if (d2.inf() > r2i.sup()) return CGAL::LARGER;
    }
    ++counters.exact_membership;
    const OneRoot dx = p.x - c.x();
    const OneRoot dy = p.y - c.y();
    return static_cast<CGAL::Comparison_result>(sign_one_root(dx * dx + dy * dy - r2));
}

CGAL::Sign line_side(const LineSupport& line, const ExactPoint& p, const Interval& px, const Interval& py,
                     const FilterMode mode, DecisionCounters& counters)
{
    if (mode == FilterMode::Lazy) {
        const Interval v = line.nx_i * px + line.ny_i * py - line.offset_i;
        if (v.inf() > 0) return CGAL::POSITIVE;
        if (v.sup() < 0) return CGAL::NEGATIVE;
    }
    ++counters.exact_line_side;
    return line_side(line, p);
}

ExactPoint line_point(const LineSupport& line, const FT& t)
{
    return {make_one_root(line.b0x + t * line.dx(), line.b1x, line.L),
            make_one_root(line.b0y + t * line.dy(), line.b1y, line.L)};
}

int line_coincidence(const LineSupport& a, const LineSupport& b)
{
    if (CGAL::sign(a.nx * b.ny - a.ny * b.nx) != CGAL::ZERO) return 0;
    // b.n = lambda a.n; same line iff lambda (k0a + k1a sqrt(La)) == k0b + k1b sqrt(Lb).
    const FT lambda = (a.nx * b.nx + a.ny * b.ny) / a.nn;
    const Tower lhs{OneRoot(lambda * a.k0 - b.k0), OneRoot(lambda * a.k1), OneRoot(a.L)};
    if (sign_tower_plus(lhs, OneRoot(-b.k1), OneRoot(b.L)) != CGAL::ZERO) return 0;
    return CGAL::sign(lambda) == CGAL::POSITIVE ? 1 : -1;
}

LineTable::LineTable(const CentreDomain& domain, const FT& tool_radius) : tool_radius_(tool_radius)
{
    for (const auto& loop : domain.loops) {
        for (std::size_t i = 0; i < loop.size(); ++i) {
            starts_.push_back(loop[i]);
            ends_.push_back(loop[(i + 1) % loop.size()]);
        }
    }
}

Point_2 LineTable::edge_point(const std::uint32_t edge, const FT& tau) const
{
    const Point_2& a = starts_[edge];
    const Point_2& b = ends_[edge];
    return Point_2(a.x() + tau * (b.x() - a.x()), a.y() + tau * (b.y() - a.y()));
}

std::uint32_t LineTable::flank(const std::uint32_t edge, const int sigma)
{
    const auto key = std::make_pair(edge, sigma);
    if (const auto it = flanks_.find(key); it != flanks_.end()) return it->second;
    const Point_2& a = starts_[edge];
    const Point_2& b = ends_[edge];
    const FT dx = b.x() - a.x();
    const FT dy = b.y() - a.y();
    // n = (-dy, dx): left normal; |n|^2 = L; n.x = n.a + sigma r sqrt(L).
    const auto id = static_cast<std::uint32_t>(lines_.size());
    lines_.push_back(LineSupport::build(-dy, dx, a, FT(sigma) * tool_radius_, dx * dx + dy * dy));
    flanks_.emplace(key, id);
    return id;
}

std::uint32_t LineTable::end_line(const std::uint32_t edge, const FT& tau)
{
    const auto key = std::make_pair(edge, tau);
    if (const auto it = ends_by_tau_.find(key); it != ends_by_tau_.end()) return it->second;
    const Point_2& a = starts_[edge];
    const Point_2& b = ends_[edge];
    const auto id = static_cast<std::uint32_t>(lines_.size());
    lines_.push_back(LineSupport::build(b.x() - a.x(), b.y() - a.y(), edge_point(edge, tau), FT(0), FT(0)));
    ends_by_tau_.emplace(key, id);
    return id;
}

std::uint32_t LineTable::add_rational_line(const FT& nx, const FT& ny, const Point_2& through)
{
    const auto id = static_cast<std::uint32_t>(lines_.size());
    lines_.push_back(LineSupport::build(nx, ny, through, FT(0), FT(0)));
    return id;
}

RectSide rect_membership(const RectRegion& rect, const LineTable& lines, const ExactPoint& p)
{
    bool boundary = false;
    for (std::size_t i = 0; i < 4; ++i) {
        const CGAL::Sign s = line_side(lines[rect.lines[i]], p);
        if (s == CGAL::ZERO) {
            boundary = true;
        } else if (s != rect.inward[i]) {
            return RectSide::Exterior;
        }
    }
    return boundary ? RectSide::Boundary : RectSide::Interior;
}

RectRegion make_rect(LineTable& lines, const std::uint32_t edge, const FT& tau0, const FT& tau1, const bool target,
                     const std::uint32_t report_index)
{
    if (CGAL::compare(tau0, tau1) != CGAL::SMALLER) throw CoverageCertificateNativeError("rectangle with tau0 >= tau1");
    RectRegion rect;
    rect.edge = edge;
    rect.tau0 = tau0;
    rect.tau1 = tau1;
    rect.target = target;
    rect.report_index = report_index;
    rect.lines = {lines.flank(edge, 1), lines.flank(edge, -1), lines.end_line(edge, tau0), lines.end_line(edge, tau1)};
    // Inside: left flank n.x <= ..., right flank n.x >= ..., end lines d.x >= d.p0 and <= d.p1.
    rect.inward = {CGAL::NEGATIVE, CGAL::POSITIVE, CGAL::POSITIVE, CGAL::NEGATIVE};
    // Certified box: the corners p_i +/- r n / |n|, n = left normal, in intervals.
    const LineSupport& flank = lines[rect.lines[0]];
    const Point_2 p0 = lines.edge_point(edge, tau0);
    const Point_2 p1 = lines.edge_point(edge, tau1);
    // The flank offset vector is (shift / |n|^2) n with shift = r sqrt(L) (folded or not):
    // b1 sqrt(L) for an irrational flank, (b0 - edge start) otherwise; use the line's shadows.
    const Interval ox = flank.b0x_i - interval_of(lines.edge_start(edge).x()) + flank.b1x_i * flank.sqrt_l_i;
    const Interval oy = flank.b0y_i - interval_of(lines.edge_start(edge).y()) + flank.b1y_i * flank.sqrt_l_i;
    double xmin = 0, xmax = 0, ymin = 0, ymax = 0;
    bool first = true;
    for (const Point_2* p : {&p0, &p1}) {
        const Interval px = interval_of(p->x()), py = interval_of(p->y());
        for (const int s : {1, -1}) {
            const Interval cx = px + Interval(s) * ox;
            const Interval cy = py + Interval(s) * oy;
            if (first) {
                xmin = cx.inf(), xmax = cx.sup(), ymin = cy.inf(), ymax = cy.sup();
                first = false;
            } else {
                xmin = std::min(xmin, cx.inf()), xmax = std::max(xmax, cx.sup());
                ymin = std::min(ymin, cy.inf()), ymax = std::max(ymax, cy.sup());
            }
        }
    }
    rect.bbox = CGAL::Bbox_2(xmin, ymin, xmax, ymax);
    return rect;
}

std::uint32_t global_edge(const CentreDomain& domain, const std::uint32_t loop, const std::uint32_t segment)
{
    if (loop >= domain.loops.size()) throw InvalidCertificateInputError("connector names a loop the domain lacks");
    if (segment >= domain.loops[loop].size()) throw InvalidCertificateInputError("connector names a segment the loop lacks");
    std::uint32_t offset = 0;
    for (std::uint32_t l = 0; l < loop; ++l) offset += static_cast<std::uint32_t>(domain.loops[l].size());
    return offset + segment;
}

ConnectorCapsules ConnectorCapsules::build(const CentreDomain& domain, LineTable& lines,
                                           const std::vector<ConnectorSpec>& connectors)
{
    // Straight pieces per global edge: (tau0, tau1, connector index).
    std::vector<std::vector<std::tuple<FT, FT, std::uint32_t>>> per_edge(lines.edge_count());
    ConnectorCapsules out;
    const FT zero(0), one(1);
    for (std::uint32_t c = 0; c < connectors.size(); ++c) {
        const ConnectorSpec& spec = connectors[c];
        if (spec.loop >= domain.loops.size()) throw InvalidCertificateInputError("connector names a loop the domain lacks");
        const auto n = static_cast<std::uint32_t>(domain.loops[spec.loop].size());
        if (spec.segment0 >= n || spec.segment1 >= n) throw InvalidCertificateInputError("connector names a segment the loop lacks");
        for (const FT* tau : {&spec.tau0, &spec.tau1}) {
            if (CGAL::sign(*tau) == CGAL::NEGATIVE || CGAL::compare(*tau, one) != CGAL::SMALLER) {
                throw InvalidCertificateInputError("connector site parameter must lie in [0, 1)");
            }
        }
        const auto add = [&](const std::uint32_t segment, const FT& a, const FT& b) {
            if (CGAL::compare(a, b) == CGAL::SMALLER) {
                per_edge[global_edge(domain, spec.loop, segment)].emplace_back(a, b, c);
                ++out.piece_count;
            }
        };
        if (spec.segment0 == spec.segment1 && CGAL::compare(spec.tau0, spec.tau1) == CGAL::SMALLER) {
            add(spec.segment0, spec.tau0, spec.tau1);
            continue;
        }
        if (spec.segment0 == spec.segment1 && spec.tau0 == spec.tau1) {
            throw InvalidCertificateInputError("connector " + std::to_string(c) + " makes no forward progress");
        }
        // Forward from (segment0, tau0) to the end of that segment, whole segments,
        // then the start of segment1 up to tau1 (wrapping once around the loop).
        add(spec.segment0, spec.tau0, one);
        for (std::uint32_t s = (spec.segment0 + 1) % n; s != spec.segment1; s = (s + 1) % n) add(s, zero, one);
        add(spec.segment1, zero, spec.tau1);
    }
    std::vector<std::pair<Point_2, std::uint32_t>> ends;
    for (std::uint32_t e = 0; e < per_edge.size(); ++e) {
        auto& pieces = per_edge[e];
        if (pieces.empty()) continue;
        std::sort(pieces.begin(), pieces.end(), [](const auto& p, const auto& q) {
            const auto c = CGAL::compare(std::get<0>(p), std::get<0>(q));
            if (c != CGAL::EQUAL) return c == CGAL::SMALLER;
            return std::get<2>(p) < std::get<2>(q);
        });
        FT lo = std::get<0>(pieces[0]);
        FT hi = std::get<1>(pieces[0]);
        std::uint32_t lowest = std::get<2>(pieces[0]);
        std::uint32_t lo_owner = lowest;
        std::uint32_t hi_owner = lowest;
        const auto flush = [&]() {
            out.rects.push_back(make_rect(lines, e, lo, hi, false, lowest));
            ends.emplace_back(lines.edge_point(e, lo), lo_owner);
            ends.emplace_back(lines.edge_point(e, hi), hi_owner);
        };
        for (std::size_t i = 1; i < pieces.size(); ++i) {
            const auto& [a, b, c] = pieces[i];
            if (CGAL::compare(a, hi) != CGAL::LARGER) {  // overlaps or abuts: exact union
                lowest = std::min(lowest, c);
                if (CGAL::compare(b, hi) == CGAL::LARGER) {
                    hi = b;
                    hi_owner = c;
                } else if (b == hi) {
                    hi_owner = std::min(hi_owner, c);
                }
                continue;
            }
            flush();
            lo = a;
            hi = b;
            lowest = c;
            lo_owner = c;
            hi_owner = c;
        }
        flush();
    }
    std::sort(ends.begin(), ends.end(), [](const auto& p, const auto& q) {
        const auto c = CGAL::compare_xy(p.first, q.first);
        if (c != CGAL::EQUAL) return c == CGAL::SMALLER;
        return p.second < q.second;
    });
    for (std::size_t i = 0; i < ends.size(); ++i) {
        if (i > 0 && ends[i].first == ends[i - 1].first) continue;  // lowest index kept (sorted)
        out.disk_centres.push_back(ends[i].first);
        out.disk_report_index.push_back(ends[i].second);
    }
    return out;
}

}  // namespace compas_cgal::coverage_certificate
