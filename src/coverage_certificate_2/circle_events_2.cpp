#include "coverage_certificate_2/circle_events_2.h"

#include "coverage_certificate_2/errors_2.h"
#include "coverage_certificate_2/separation_2.h"

namespace compas_cgal::coverage_certificate {

namespace {
double mid(const Interval& i) { return 0.5 * (i.inf() + i.sup()); }
}  // namespace

void CircleEvents::add_line(const std::uint32_t k, const std::uint32_t line, DecisionCounters& counters)
{
    const BoundaryCircle& circle = set_->circles[k];
    CircleLineRecord record = make_circle_line_record(k, circle, line, (*lines_)[line], mode(), counters);
    if (record.count == 0) return;
    const auto rid = static_cast<std::uint32_t>(records_.size());
    const std::uint8_t multiplicity = record.count == 2 ? 1 : 2;
    for (std::uint8_t b = 0; b < (record.count == 2 ? 2 : 1); ++b) {
        table_->events[k].push_back({rid, b, multiplicity, line, record.x[b], record.y[b], 1});
    }
    records_.push_back(std::move(record));
}

void CircleEvents::add_vertex(const std::uint32_t k, const std::uint32_t vertex, DecisionCounters& counters)
{
    const BoundaryCircle& circle = set_->circles[k];
    CircleVertexRecord record = make_circle_vertex_record(k, circle, vertex, (*vertices_)[vertex], mode(), counters);
    if (record.count == 0) return;
    const auto rid = static_cast<std::uint32_t>(vertex_records_.size());
    const std::uint8_t multiplicity = record.count == 2 ? 1 : 2;
    for (std::uint8_t b = 0; b < (record.count == 2 ? 2 : 1); ++b) {
        table_->events[k].push_back({rid, b, multiplicity, vertex, record.x[b], record.y[b], 2});
    }
    vertex_records_.push_back(std::move(record));
}

const CircleVertexExact& CircleEvents::exact_vertex(const CrossingEvent& event, DecisionCounters& counters)
{
    CircleVertexRecord& record = vertex_records_[event.pair];
    return exact_of(record, set_->circles[record.circle], (*vertices_)[record.vertex], counters);
}

const CircleLineExact& CircleEvents::exact_line(const CrossingEvent& event, DecisionCounters& counters)
{
    CircleLineRecord& record = records_[event.pair];
    return exact_of(record, set_->circles[record.circle], (*lines_)[record.line], counters);
}

Tower CircleEvents::exact_x(const CrossingEvent& event, DecisionCounters& counters)
{
    if (event.source == 0) return table_->exact_x(event, counters);
    if (event.source == 2) return exact_vertex(event, counters).x[event.branch];
    return exact_line(event, counters).x[event.branch];
}

std::pair<Tower, Tower> CircleEvents::exact_xy(const CrossingEvent& event, DecisionCounters& counters)
{
    if (event.source == 0) return {table_->exact_x(event, counters), table_->exact_y(event, counters)};
    if (event.source == 2) {
        const CircleVertexExact& v = exact_vertex(event, counters);
        return {v.x[event.branch], v.y[event.branch]};
    }
    const CircleLineExact& e = exact_line(event, counters);
    return {e.x[event.branch], e.y[event.branch]};
}

int CircleEvents::half(const CrossingEvent& event, const BoundaryCircle& circle, DecisionCounters& counters)
{
    if (event.source == 0) return table_->half(event, circle, counters);
    if (mode() == FilterMode::Lazy) {
        if (event.y.sup() < circle.centre_y.inf()) return 0;
        if (event.y.inf() > circle.centre_y.sup()) return 1;
    }
    ++counters.exact_half;
    const auto [ex, ey] = exact_xy(event, counters);
    const auto& c = circle.circle.center();
    const CGAL::Sign sy = sign_difference(ey, Tower::rational(c.y()));
    if (sy == CGAL::NEGATIVE) return 0;
    if (sy == CGAL::POSITIVE) return 1;
    return sign_difference(ex, Tower::rational(c.x())) == CGAL::NEGATIVE ? 0 : 1;
}

CGAL::Comparison_result CircleEvents::compare_x(const CrossingEvent& a, const CrossingEvent& b,
                                                DecisionCounters& counters)
{
    if (a.source == b.source && a.pair == b.pair && a.branch == b.branch) return CGAL::EQUAL;
    if (a.source == 0 && b.source == 0) return table_->compare_x(a, b, counters);
    if (mode() == FilterMode::Lazy) {
        if (a.x.sup() < b.x.inf()) return CGAL::SMALLER;
        if (a.x.inf() > b.x.sup()) return CGAL::LARGER;
    }
    ++counters.exact_tower_compares;
    const Tower ta = exact_x(a, counters);
    const Tower tb = exact_x(b, counters);
    if (!compatible_bases(ta, tb)) return separate_cross_base(a, b, counters);
    return static_cast<CGAL::Comparison_result>(sign_difference(ta, tb));
}

CGAL::Comparison_result CircleEvents::compare_x(const CrossingEvent& a, const Point_2& q, const Interval& qx,
                                                DecisionCounters& counters)
{
    if (a.source == 0) return table_->compare_x(a, q, qx, counters);
    if (mode() == FilterMode::Lazy) {
        if (a.x.sup() < qx.inf()) return CGAL::SMALLER;
        if (a.x.inf() > qx.sup()) return CGAL::LARGER;
    }
    ++counters.exact_compare_x;
    return static_cast<CGAL::Comparison_result>(sign_difference(exact_x(a, counters), Tower::rational(q.x())));
}

std::pair<double, double> CircleEvents::proposal_xy(const CrossingEvent& event, DecisionCounters& counters)
{
    if (event.source == 0) {
        const auto& p = table_->exact_point(event, counters);
        return {mid(Interval(CGAL::to_interval(p.x()))), mid(Interval(CGAL::to_interval(p.y())))};
    }
    const auto [ex, ey] = exact_xy(event, counters);
    return {mid(interval_of(ex)), mid(interval_of(ey))};
}

bool CircleEvents::same_point_by_provenance(const CrossingEvent& a, const CrossingEvent& b, DecisionCounters& counters)
{
    const auto [xa, ya] = exact_xy(a, counters);
    if (b.source == 2) {
        const CircleVertexRecord& rb = vertex_records_[b.pair];
        const CircleVertexExact& eb = exact_vertex(b, counters);
        return is_circle_vertex_crossing(xa, ya, eb, b.branch, (*vertices_)[rb.vertex]);
    }
    const CircleLineRecord& rb = records_[b.pair];
    const LineSupport& lb = (*lines_)[rb.line];
    // a's point on b's line, and on b's branch.
    const Tower on_line{lb.nx * xa.m + lb.ny * ya.m - lb.k0, lb.nx * xa.n + lb.ny * ya.n, xa.d};
    if (sign_tower_plus(on_line, OneRoot(-lb.k1), OneRoot(lb.L)) != CGAL::ZERO) return false;
    if (rb.count == 1) return true;
    // The chord midpoint of line b projects onto d_b at d_b . c; the + branch
    // lies beyond it.
    const auto& c = set_->circles[rb.circle].circle.center();
    const Tower along{lb.dx() * (xa.m - c.x()) + lb.dy() * (ya.m - c.y()), lb.dx() * xa.n + lb.dy() * ya.n, xa.d};
    return sign_tower(along) == (b.branch == 0 ? CGAL::POSITIVE : CGAL::NEGATIVE);
}

CGAL::Comparison_result CircleEvents::separate_cross_base(const CrossingEvent& a, const CrossingEvent& b,
                                                          DecisionCounters& counters)
{
    ++counters.cross_base_separations;
    if (a.source == 0 || b.source == 0) {
        throw CoverageCertificateNativeError("cross-base separation requested for a circle-circle event");
    }
    if (a.source == 2 && b.source == 2 && a.pair == b.pair) {
        throw CoverageCertificateNativeError("cross-base separation of two crossings of one vertex circle");
    }
    // Equality by provenance: a's point on b's support and on b's branch.
    if (same_point_by_provenance(a, b, counters)) return CGAL::EQUAL;
    const Tower xa = exact_x(a, counters);
    const Tower xb = exact_x(b, counters);
    // Distinct points on one half: separate their x by exact rational bisection.
    const Interval ia = interval_of(xa);
    const Interval ib = interval_of(xb);
    FT lo(std::min(ia.inf(), ib.inf()));
    FT hi(std::max(ia.sup(), ib.sup()));
    // The x's differ (equality was decided above), so the bisection ends once
    // hi - lo < |xa - xb|: within the root-separation bound of xa - xb.
    SeparationBudget budget("cross-base separation");
    const double width = std::max(ia.sup(), ib.sup()) - std::min(ia.inf(), ib.inf());
    for (;;) {
        budget.step([&] {
            return halving_step_bound(width, separation_bits(plus(SeparationMeasure::of(xa), SeparationMeasure::of(xb)),
                                                             tower_radicals({&xa, &xb})));
        });
        const FT m = (lo + hi) / FT(2);
        const CGAL::Sign sa = sign_difference(xa, Tower::rational(m));
        const CGAL::Sign sb = sign_difference(xb, Tower::rational(m));
        if (sa != sb) return static_cast<int>(sa) < static_cast<int>(sb) ? CGAL::SMALLER : CGAL::LARGER;
        if (sa == CGAL::ZERO) throw CoverageCertificateNativeError("distinct cross-base events share an x coordinate");
        if (sa == CGAL::POSITIVE) {
            lo = m;
        } else {
            hi = m;
        }
    }
}

}  // namespace compas_cgal::coverage_certificate
