#include "coverage_certificate_2/line_sweep_2.h"

#include "coverage_certificate_2/errors_2.h"
#include "coverage_certificate_2/separation_2.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>

namespace compas_cgal::coverage_certificate {

namespace {
/// Witness proposals are rounded to dyadics of PROPOSAL_BITS_STEP, 2*STEP, ...
/// up to PROPOSAL_BITS_MAX fractional bits (binary64 has 52). Proposal policy
/// only: shorter dyadics give smaller exact witnesses; acceptance is exact.
constexpr int PROPOSAL_BITS_STEP = 4;
constexpr int PROPOSAL_BITS_MAX = 52;

/// bound < q: decided by the bound's certified enclosure when it excludes q,
/// exactly otherwise.
bool strictly_below(const ChartBound& bound, const FT& q, const Interval& qi)
{
    const Interval b = bound.interval();
    if (b.sup() < qi.inf()) return true;
    if (b.inf() > qi.sup()) return false;
    return bound.minus(q) == CGAL::NEGATIVE;
}

bool strictly_above(const ChartBound& bound, const FT& q, const Interval& qi)
{
    const Interval b = bound.interval();
    if (b.inf() > qi.sup()) return true;
    if (b.sup() < qi.inf()) return false;
    return bound.minus(q) == CGAL::POSITIVE;
}

/// Separation measure and radical count of one chart value.
std::pair<SeparationMeasure, int> measure_of(const ChartBound& bound)
{
    if (bound.vertex != nullptr) return {bound.vertex->measure(), LineVertexExact::RADICALS};
    const Tower t = bound.tower();
    return {SeparationMeasure::of(t), tower_radicals({&t})};
}

double mid(const Interval& i) { return 0.5 * (i.inf() + i.sup()); }

/// Root-separation step bound of a bisection separating two chart values.
long long chart_bisection_bound(const ChartBound& lo, const ChartBound& hi, const double width)
{
    const auto [ml, rl] = measure_of(lo);
    const auto [mh, rh] = measure_of(hi);
    int radicals = rl + rh;
    if (lo.vertex == nullptr && hi.vertex == nullptr) {
        const Tower tl = lo.tower();
        const Tower th = hi.tower();
        radicals = tower_radicals({&tl, &th});
    }
    return halving_step_bound(width, separation_bits(plus(mh, ml), radicals));
}

}  // namespace

Side membership(const CircleRegion& region, const ExactPoint& p)
{
    DecisionCounters scratch;
    const Interval px = interval_of(p.x), py = interval_of(p.y);
    const auto outer = compare_to_circle(p, px, py, region.centre, region.centre_x, region.centre_y, region.outer_squared,
                                         region.outer_squared_i, FilterMode::ExactOnly, scratch);
    if (outer == CGAL::LARGER) return Side::Exterior;
    if (outer == CGAL::EQUAL) return Side::Boundary;
    if (!region.annulus) return Side::Interior;
    const auto inner = compare_to_circle(p, px, py, region.centre, region.centre_x, region.centre_y, region.inner_squared,
                                         region.inner_squared_i, FilterMode::ExactOnly, scratch);
    if (inner == CGAL::SMALLER) return Side::Exterior;
    if (inner == CGAL::EQUAL) return Side::Boundary;
    return Side::Interior;
}

CGAL::Sign ChartBound::minus(const FT& value) const
{
    if (vertex != nullptr) return vertex->sign_minus(branch, value);
    if (rational) return static_cast<CGAL::Sign>(CGAL::compare(q, value));
    return sign_difference(t, Tower::rational(value));
}

const LineVertexExact& LineSweep::exact_vertex(const LineEvent& event, DecisionCounters& counters)
{
    LineVertexRecord& r = vertex_records[event.record];
    return exact_of(r, (*lines)[piece->line], (*vertices)[r.vertex], counters);
}

CGAL::Sign LineSweep::minus_rational(const LineEvent& event, const FT& q, DecisionCounters& counters)
{
    if (event.source == 2) return exact_vertex(event, counters).sign_minus(event.branch, q);
    return sign_difference(exact_t(event, counters), Tower::rational(q));
}

CGAL::Comparison_result LineSweep::compare_events(const LineEvent& a, const LineEvent& b, DecisionCounters& counters)
{
    if (a.source == b.source && a.record == b.record && a.branch == b.branch) return CGAL::EQUAL;
    if (a.source != 2 && b.source != 2) {
        return static_cast<CGAL::Comparison_result>(sign_difference(exact_t(a, counters), exact_t(b, counters)));
    }
    ++counters.exact_target_events;
    if (a.source == 2 && b.source == 2 && a.record == b.record) {
        // Two roots of one vertex circle: the + root (branch 0) is the larger.
        return a.branch == 0 ? CGAL::LARGER : CGAL::SMALLER;
    }
    // A vertex root against a tower event: one exact sign in K = Q(sqrt(L))(sqrt(delta)).
    if (a.source != 2 || b.source != 2) {
        const bool a_vertex = a.source == 2;
        const LineEvent& v = a_vertex ? a : b;
        const LineEvent& o = a_vertex ? b : a;
        const CGAL::Sign s = exact_vertex(v, counters).sign_minus(v.branch, exact_t(o, counters));
        const auto c = static_cast<CGAL::Comparison_result>(a_vertex ? s : CGAL::opposite(s));
        return c;
    }
    // Two different vertex circles (two radicands): separate by exact rational
    // bisection inside the enclosures (a common point of a flank and two vertex
    // circles is a three-support coincidence; the budget turns it into a loud error).
    const auto bound_of = [&](const LineEvent& e) {
        if (e.source == 2) return ChartBound{false, FT(0), Tower::rational(FT(0)), &exact_vertex(e, counters), e.branch, e.t};
        return ChartBound{false, FT(0), exact_t(e, counters)};
    };
    const ChartBound ba = bound_of(a);
    const ChartBound bb = bound_of(b);
    const double lo_d = std::min(a.t.inf(), b.t.inf());
    const double hi_d = std::max(a.t.sup(), b.t.sup());
    FT lo(lo_d);
    FT hi(hi_d);
    SeparationBudget budget("line-event separation");
    for (;;) {
        budget.step([&] { return chart_bisection_bound(ba, bb, hi_d - lo_d); });
        const FT m = (lo + hi) / FT(2);
        const CGAL::Sign sa = ba.minus(m);
        const CGAL::Sign sb = bb.minus(m);
        if (sa != sb) return static_cast<int>(sa) < static_cast<int>(sb) ? CGAL::SMALLER : CGAL::LARGER;
        if (sa == CGAL::ZERO) throw CoverageCertificateNativeError("distinct line events share a chart parameter");
        if (sa == CGAL::POSITIVE) {
            lo = m;
        } else {
            hi = m;
        }
    }
}

Tower LineSweep::exact_t(const LineEvent& event, DecisionCounters& counters)
{
    if (event.source == 2) throw CoverageCertificateNativeError("a vertex-circle root has no tower parameter");
    const LineSupport& support = (*lines)[piece->line];
    if (event.source == 0) {
        CircleLineRecord& r = circle_records[event.record];
        return exact_of(r, set->circles[r.circle], support, counters).t[event.branch];
    }
    LineLineRecord& r = line_records[event.record];
    return exact_t_of(r, support, (*lines)[r.other], counters);
}

std::pair<ChartBound, ChartBound> sub_edge_bounds(LineSweep& sweep, const std::size_t edge, DecisionCounters& counters)
{
    const std::size_t vertices = sweep.vertex_count();
    ChartBound lo{true, sweep.piece->t_lo, Tower::rational(FT(0))};
    ChartBound hi{true, sweep.piece->t_hi, Tower::rational(FT(0))};
    const auto bound_of = [&](const LineEvent& e) {
        if (e.source == 2) {
            return ChartBound{false, FT(0), Tower::rational(FT(0)), &sweep.exact_vertex(e, counters), e.branch, e.t};
        }
        return ChartBound{false, FT(0), sweep.exact_t(e, counters)};
    };
    if (edge > 0) lo = bound_of(sweep.events[sweep.group_starts[edge - 1]]);
    if (edge < vertices) hi = bound_of(sweep.events[sweep.group_starts[edge]]);
    return {lo, hi};
}

FT rational_between(const ChartBound& lo, const ChartBound& hi, DecisionCounters& counters)
{
    const Interval ilo = lo.interval();
    const Interval ihi = hi.interval();
    const double seed = 0.5 * (mid(ilo) + mid(ihi));
    if (std::isfinite(seed)) {
        for (int bits = PROPOSAL_BITS_STEP; bits <= PROPOSAL_BITS_MAX; bits += PROPOSAL_BITS_STEP) {
            const double q = std::ldexp(std::round(std::ldexp(seed, bits)), -bits);
            if (!std::isfinite(q)) continue;
            const FT candidate(q);
            const Interval ci(q);
            if (strictly_below(lo, candidate, ci) && strictly_above(hi, candidate, ci)) return candidate;
        }
    }
    // Exact bisection inside the certified enclosures: a <= lo < hi <= b. Each
    // step halves b - a, and a midpoint lies strictly between lo and hi once
    // b - a < hi - lo, so the loop ends within the root-separation bound of
    // hi - lo (separation_2.h), computed once the lazy threshold is passed.
    FT a(ilo.inf());
    FT b(ihi.sup());
    SeparationBudget budget("P3 chart bisection");
    const double width = ihi.sup() - ilo.inf();
    for (;;) {
        budget.step([&] { return chart_bisection_bound(lo, hi, width); });
        ++counters.witness_bisections;
        const FT m = (a + b) / FT(2);
        const Interval mi = interval_of(m);
        if (!strictly_below(lo, m, mi)) {
            a = m;
        } else if (!strictly_above(hi, m, mi)) {
            b = m;
        } else {
            return m;
        }
    }
}

LineSweep sweep_line(const LinePiece& piece, const CircleRegionSet& set, const LineTable& lines,
                     const std::vector<RectRegion>& rects, const LineCandidates& candidates, const FilterMode mode,
                     DecisionCounters& counters, const std::vector<VertexCircle>* vertex_circles)
{
    ++counters.line_pieces;
    const LineSupport& support = lines[piece.line];
    LineSweep sweep;
    sweep.set = &set;
    sweep.lines = &lines;
    sweep.piece = &piece;
    sweep.vertices = vertex_circles;
    sweep.regions = candidates.regions;

    // Rectangles: coincident lines (P5/P6) first.
    std::vector<std::uint32_t> line_ids;
    for (const std::uint32_t r : candidates.rects) {
        const RectRegion& R = rects[r];
        std::uint8_t coincident_slot = 4;
        bool exclude = false;
        for (std::uint8_t s = 0; s < 4; ++s) {
            const int co = R.lines[s] == piece.line ? 1 : line_coincidence(support, lines[R.lines[s]]);
            if (co == 0) continue;
            const CGAL::Sign relative = co > 0 ? R.inward[s] : CGAL::opposite(R.inward[s]);
            if (relative == piece.check_side) {
                coincident_slot = s;
            } else {
                exclude = true;
            }
        }
        if (coincident_slot != 4 || exclude) ++counters.coincident_line_supports;
        if (exclude) {
            sweep.excluded_rects.push_back(r);
            continue;
        }
        sweep.rects.push_back(r);
        sweep.rect_coincident_slot.push_back(coincident_slot);
        for (std::uint8_t s = 0; s < 4; ++s) {
            if (s != coincident_slot) line_ids.push_back(R.lines[s]);
        }
    }
    for (const std::uint32_t q : candidates.target_lines) {
        if (q != piece.line) line_ids.push_back(q);
    }
    std::sort(line_ids.begin(), line_ids.end());
    line_ids.erase(std::unique(line_ids.begin(), line_ids.end()), line_ids.end());

    // Events from circles and lines, restricted to (t_lo, t_hi).
    std::vector<std::uint32_t> circle_ids;
    for (const std::uint32_t j : sweep.regions) {
        for (std::uint32_t c = set.first_circle[j]; c < set.first_circle[j + 1]; ++c) circle_ids.push_back(c);
    }
    for (const std::uint32_t c : candidates.target_circles) circle_ids.push_back(c);
    const Interval lo_i = interval_of(piece.t_lo), hi_i = interval_of(piece.t_hi);
    std::vector<LineEvent> raw;
    for (const std::uint32_t c : circle_ids) {
        CircleLineRecord record = make_circle_line_record(c, set.circles[c], piece.line, support, mode, counters);
        if (record.count == 0) continue;
        const auto rid = static_cast<std::uint32_t>(sweep.circle_records.size());
        const std::uint8_t multiplicity = record.count == 2 ? 1 : 2;
        for (std::uint8_t b = 0; b < (record.count == 2 ? 2 : 1); ++b) raw.push_back({0, rid, b, multiplicity, c, record.t[b]});
        sweep.circle_records.push_back(std::move(record));
    }
    for (const std::uint32_t q : line_ids) {
        LineLineRecord record = make_line_line_record(support, q, lines[q], false, mode, counters);
        if (record.relation != LineRelation::Crossing) continue;
        const auto rid = static_cast<std::uint32_t>(sweep.line_records.size());
        raw.push_back({1, rid, 0, 1, q, record.t});
        sweep.line_records.push_back(std::move(record));
    }
    for (const std::uint32_t v : candidates.target_vertices) {
        LineVertexRecord record = make_line_vertex_record(support, v, (*vertex_circles)[v], mode, counters);
        if (record.count == 0) continue;
        const auto rid = static_cast<std::uint32_t>(sweep.vertex_records.size());
        const std::uint8_t multiplicity = record.count == 2 ? 1 : 2;
        for (std::uint8_t b = 0; b < (record.count == 2 ? 2 : 1); ++b) raw.push_back({2, rid, b, multiplicity, v, record.t[b]});
        sweep.vertex_records.push_back(std::move(record));
    }
    const auto inside = [&](const LineEvent& e) {
        CGAL::Sign above_lo = CGAL::ZERO, below_hi = CGAL::ZERO;
        bool decided_lo = false, decided_hi = false;
        if (mode == FilterMode::Lazy) {
            if (e.t.sup() < lo_i.inf() || e.t.inf() > hi_i.sup()) return false;
            if (e.t.inf() > lo_i.sup()) above_lo = CGAL::POSITIVE, decided_lo = true;
            if (e.t.sup() < hi_i.inf()) below_hi = CGAL::POSITIVE, decided_hi = true;
        }
        if (!decided_lo || !decided_hi) {
            ++counters.exact_tower_compares;
            if (!decided_lo) above_lo = sweep.minus_rational(e, piece.t_lo, counters);
            if (!decided_hi) below_hi = CGAL::opposite(sweep.minus_rational(e, piece.t_hi, counters));
        }
        return above_lo == CGAL::POSITIVE && below_hi == CGAL::POSITIVE;
    };
    for (const LineEvent& e : raw) {
        if (inside(e)) sweep.events.push_back(e);
    }
    const auto compare = [&](const LineEvent& a, const LineEvent& b) {
        if (a.source == b.source && a.record == b.record && a.branch == b.branch) return CGAL::EQUAL;
        if (mode == FilterMode::Lazy) {
            if (a.t.sup() < b.t.inf()) return CGAL::SMALLER;
            if (a.t.inf() > b.t.sup()) return CGAL::LARGER;
        }
        ++counters.exact_tower_compares;
        return sweep.compare_events(a, b, counters);
    };
    std::sort(sweep.events.begin(), sweep.events.end(),
              [&](const LineEvent& a, const LineEvent& b) { return compare(a, b) == CGAL::SMALLER; });
    for (std::uint32_t i = 0; i < sweep.events.size(); ++i) {
        if (i == 0 || compare(sweep.events[i], sweep.events[i - 1]) != CGAL::EQUAL) sweep.group_starts.push_back(i);
    }
    if (!sweep.events.empty()) sweep.group_starts.push_back(static_cast<std::uint32_t>(sweep.events.size()));

    // Start state at a rational parameter inside the first sub-edge.
    const auto [first_lo, first_hi] = sub_edge_bounds(sweep, 0, counters);
    const FT t_start = rational_between(first_lo, first_hi, counters);
    const ExactPoint w = line_point(support, t_start);
    const Interval wx = interval_of(w.x), wy = interval_of(w.y);
    std::vector<std::uint8_t> in_outer(sweep.regions.size(), 0), in_inner(sweep.regions.size(), 0);
    for (std::size_t t = 0; t < sweep.regions.size(); ++t) {
        const CircleRegion& Y = set.regions[sweep.regions[t]];
        const auto co = compare_to_circle(w, wx, wy, Y.centre, Y.centre_x, Y.centre_y, Y.outer_squared, Y.outer_squared_i,
                                          mode, counters);
        const auto ci = Y.annulus ? compare_to_circle(w, wx, wy, Y.centre, Y.centre_x, Y.centre_y, Y.inner_squared,
                                                      Y.inner_squared_i, mode, counters)
                                  : CGAL::LARGER;
        if (co == CGAL::EQUAL || ci == CGAL::EQUAL) {
            throw CoverageCertificateNativeError("line start witness lies on a circle the crossings did not predict");
        }
        in_outer[t] = co == CGAL::SMALLER;
        in_inner[t] = Y.annulus && ci == CGAL::SMALLER;
    }
    std::vector<CGAL::Sign> line_state(line_ids.size(), CGAL::ZERO);
    for (std::size_t t = 0; t < line_ids.size(); ++t) {
        const LineSupport& Q = lines[line_ids[t]];
        const Interval v = Q.nx_i * wx + Q.ny_i * wy - Q.offset_i;
        if (mode == FilterMode::Lazy && (v.inf() > 0 || v.sup() < 0)) {
            line_state[t] = v.inf() > 0 ? CGAL::POSITIVE : CGAL::NEGATIVE;
        } else {
            ++counters.exact_line_side;
            line_state[t] = line_side(Q, w);
        }
        if (line_state[t] == CGAL::ZERO) {
            // Only a support coincident with the piece passes through w; it has
            // no events and never enters a rectangle's count (handled above).
            if (line_coincidence(support, Q) == 0) {
                throw CoverageCertificateNativeError("line start witness lies on a line the crossings did not predict");
            }
        }
    }
    const auto slot_of = [](const std::vector<std::uint32_t>& list, const std::uint32_t id) -> std::size_t {
        const auto it = std::lower_bound(list.begin(), list.end(), id);
        if (it == list.end() || *it != id) {
            throw CoverageCertificateNativeError("a line event names a support outside the candidate set");
        }
        return static_cast<std::size_t>(it - list.begin());
    };
    std::vector<std::vector<std::pair<std::uint32_t, std::uint8_t>>> line_rects(line_ids.size());
    std::vector<int> rect_good(sweep.rects.size(), 0);
    for (std::uint32_t r = 0; r < sweep.rects.size(); ++r) {
        const RectRegion& R = rects[sweep.rects[r]];
        for (std::uint8_t s = 0; s < 4; ++s) {
            if (s == sweep.rect_coincident_slot[r]) {
                ++rect_good[r];
                continue;
            }
            const std::size_t ls = slot_of(line_ids, R.lines[s]);
            line_rects[ls].emplace_back(r, s);
            rect_good[r] += line_state[ls] == R.inward[s] ? 1 : 0;
        }
    }
    const auto rect_covers = [&](const std::uint32_t r) { return rect_good[r] == 4; };
    const auto region_covers = [&](const std::size_t t) { return in_outer[t] != 0U && in_inner[t] == 0U; };
    int covered = 0;
    for (std::size_t t = 0; t < sweep.regions.size(); ++t) covered += region_covers(t) ? 1 : 0;
    for (std::uint32_t r = 0; r < sweep.rects.size(); ++r) covered += rect_covers(r) ? 1 : 0;

    const std::size_t vertices = sweep.vertex_count();
    sweep.edge_exposed.assign(vertices + 1, 0U);
    sweep.edge_exposed[0] = covered == 0 ? 1U : 0U;
    ++counters.line_sub_edges;
    for (std::size_t v = 0; v < vertices; ++v) {
        for (std::uint32_t i = sweep.group_starts[v]; i < sweep.group_starts[v + 1]; ++i) {
            const LineEvent& e = sweep.events[i];
            if (e.multiplicity != 1 || e.source == 2) continue;
            if (e.source == 1) {
                const std::size_t ls = slot_of(line_ids, e.other);
                for (const auto& [r, s] : line_rects[ls]) {
                    const RectRegion& R = rects[sweep.rects[r]];
                    const int before = rect_covers(r) ? 1 : 0;
                    rect_good[r] += (CGAL::opposite(line_state[ls]) == R.inward[s] ? 1 : 0) -
                                    (line_state[ls] == R.inward[s] ? 1 : 0);
                    covered += (rect_covers(r) ? 1 : 0) - before;
                }
                line_state[ls] = CGAL::opposite(line_state[ls]);
                continue;
            }
            const BoundaryCircle& O = set.circles[e.other];
            if (O.is_target()) continue;
            const std::size_t t = slot_of(sweep.regions, O.owner);
            const int before = region_covers(t) ? 1 : 0;
            if (O.is_inner()) {
                in_inner[t] ^= 1U;
            } else {
                in_outer[t] ^= 1U;
            }
            covered += (region_covers(t) ? 1 : 0) - before;
        }
        sweep.edge_exposed[v + 1] = covered == 0 ? 1U : 0U;
        ++counters.line_sub_edges;
    }
    return sweep;
}

}  // namespace compas_cgal::coverage_certificate
