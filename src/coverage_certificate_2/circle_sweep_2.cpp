#include "coverage_certificate_2/circle_sweep_2.h"

#include "coverage_certificate_2/errors_2.h"

#include <algorithm>
#include <numeric>
#include <functional>
#include <string>

namespace compas_cgal::coverage_certificate {

CGAL::Comparison_result order_event_point(CircleEvents& events, const CrossingEvent& e, const std::uint8_t e_half,
                                          const Point_2& q, const Interval& qx, const std::uint8_t q_half,
                                          DecisionCounters& counters)
{
    if (e_half != q_half) return e_half < q_half ? CGAL::SMALLER : CGAL::LARGER;
    const auto c = events.compare_x(e, q, qx, counters);
    if (c == CGAL::EQUAL) return CGAL::EQUAL;
    // Lower half ascends in x, upper half descends.
    const bool e_first = e_half == 0 ? c == CGAL::SMALLER : c == CGAL::LARGER;
    return e_first ? CGAL::SMALLER : CGAL::LARGER;
}

std::uint32_t start_witness_attempt_bound(const std::size_t counted, const std::size_t target_circles,
                                          const std::size_t lines)
{
    return static_cast<std::uint32_t>(1 + 2 * (2 * counted + target_circles + lines));
}

StartWitness find_start_witness(const BoundaryCircle& C, const std::uint32_t bound,
                                const std::function<bool(const Point_2&)>& on_support, DecisionCounters& counters)
{
    for (std::uint32_t attempt = 0; attempt < bound; ++attempt) {
        StartWitness w{Point_2(), 0, attempt};
        if (attempt == 0) {
            const auto& c = C.circle.center();
            w.point = Point_2(c.x() - C.radius, c.y());
        } else {
            const FT s = FT(static_cast<long long>(2 * attempt - 1)) / FT(2);
            w.point = stereographic_point(C, s);
            w.half = half_of_stereographic(s);
        }
        if (!on_support(w.point)) return w;
        ++counters.start_retries;
    }
    throw StartWitnessSearchError("start-witness search on a boundary circle exceeded its theorem bound " +
                                  std::to_string(bound) +
                                  ": every witness lies on a support, so a support coincides with the circle "
                                  "without being recognised as coincident");
}

CircleSweep sweep_circle(const std::uint32_t k, const CircleRegionSet& set, const LineTable& lines,
                         const std::vector<RectRegion>& rects, const CircleCandidates& candidates, CircleEvents& events,
                         DecisionCounters& counters)
{
    const BoundaryCircle& C = set.circles[k];
    CrossingTable& table = events.table();
    const FilterMode mode = events.mode();
    CircleSweep sweep;
    sweep.rects = candidates.rects;
    sweep.lines = candidates.lines;
    sweep.vertices = candidates.vertices;

    // P5/P6: coincident supports. Target circles are never counted here; their
    // coincidence with a swept circle is the target decision's business.
    for (const std::uint32_t other : table.coincident[k]) {
        const BoundaryCircle& O = set.circles[other];
        if (O.is_target()) continue;
        if (O.region_inside() == C.check_inside()) {
            sweep.covered_by_coincident_support = true;
        } else {
            sweep.excluded.push_back(O.owner);
        }
    }
    std::sort(sweep.excluded.begin(), sweep.excluded.end());
    sweep.excluded.erase(std::unique(sweep.excluded.begin(), sweep.excluded.end()), sweep.excluded.end());
    for (const std::uint32_t j : table.candidates[k]) {
        if (!std::binary_search(sweep.excluded.begin(), sweep.excluded.end(), j)) sweep.counted.push_back(j);
    }
    if (sweep.covered_by_coincident_support) {
        ++counters.covered_by_coincident_support;
        return sweep;
    }

    // CCW order of the events.
    std::vector<CrossingEvent>& raw = table.events[k];
    std::vector<std::uint8_t> raw_half(raw.size());
    for (std::size_t i = 0; i < raw.size(); ++i) raw_half[i] = static_cast<std::uint8_t>(events.half(raw[i], C, counters));
    std::vector<std::uint32_t> order(raw.size());
    std::iota(order.begin(), order.end(), 0U);
    std::sort(order.begin(), order.end(), [&](const std::uint32_t i, const std::uint32_t j) {
        if (raw_half[i] != raw_half[j]) return raw_half[i] < raw_half[j];
        const auto c = events.compare_x(raw[i], raw[j], counters);
        if (c == CGAL::EQUAL) return false;
        return raw_half[i] == 0 ? c == CGAL::SMALLER : c == CGAL::LARGER;
    });
    sweep.events.reserve(raw.size());
    sweep.halves.reserve(raw.size());
    for (const auto i : order) {
        sweep.events.push_back(raw[i]);
        sweep.halves.push_back(raw_half[i]);
    }
    std::vector<CrossingEvent>().swap(raw);
    const auto& evs = sweep.events;
    for (std::uint32_t i = 0; i < evs.size(); ++i) {
        if (i == 0 || sweep.halves[i] != sweep.halves[i - 1] || events.compare_x(evs[i], evs[i - 1], counters) != CGAL::EQUAL) {
            sweep.group_starts.push_back(i);
        }
    }
    if (!evs.empty()) sweep.group_starts.push_back(static_cast<std::uint32_t>(evs.size()));

    // Local slots of every counted candidate region and line.
    const auto slot_of = [](const std::vector<std::uint32_t>& list, const std::uint32_t id) -> std::size_t {
        const auto it = std::lower_bound(list.begin(), list.end(), id);
        if (it == list.end() || *it != id) {
            throw CoverageCertificateNativeError("a crossing event names a support outside the candidate set");
        }
        return static_cast<std::size_t>(it - list.begin());
    };
    std::vector<std::uint8_t> in_outer(sweep.counted.size(), 0);
    std::vector<std::uint8_t> in_inner(sweep.counted.size(), 0);
    std::vector<CGAL::Sign> line_state(sweep.lines.size(), CGAL::ZERO);
    // For each local line, the local rectangles that use it (with the slot).
    std::vector<std::vector<std::pair<std::uint32_t, std::uint8_t>>> line_rects(sweep.lines.size());
    for (std::uint32_t r = 0; r < sweep.rects.size(); ++r) {
        const RectRegion& R = rects[sweep.rects[r]];
        for (std::uint8_t s = 0; s < 4; ++s) line_rects[slot_of(sweep.lines, R.lines[s])].emplace_back(r, s);
    }
    std::vector<int> rect_good(sweep.rects.size(), 0);
    const auto rect_covers = [&](const std::uint32_t r) { return rect_good[r] == 4; };
    const auto region_covers = [&](const std::size_t t) { return in_outer[t] != 0U && in_inner[t] == 0U; };

    // Start state at a rational witness on no support (find_start_witness),
    // with the classification of the accepted witness left in the state arrays.
    const auto& same = table.coincident[k];
    std::size_t cut_circles = 0;
    for (const std::uint32_t tc : table.target_circles[k]) {
        if (std::find(same.begin(), same.end(), tc) == same.end()) ++cut_circles;
    }
    const auto on_support = [&](const Point_2& w) {
        const Interval wx = interval_of(w.x()), wy = interval_of(w.y());
        for (std::size_t t = 0; t < sweep.counted.size(); ++t) {
            const CircleRegion& Y = set.regions[sweep.counted[t]];
            const Interval ddx = wx - Y.centre_x, ddy = wy - Y.centre_y, d2 = ddx * ddx + ddy * ddy;
            const auto cmp = [&](const Interval& r2i, const FT& r2) {
                if (mode == FilterMode::Lazy) {
                    if (d2.sup() < r2i.inf()) return CGAL::SMALLER;
                    if (d2.inf() > r2i.sup()) return CGAL::LARGER;
                }
                ++counters.exact_membership;
                return CGAL::compare(CGAL::squared_distance(w, Y.centre), r2);
            };
            const auto co = cmp(Y.outer_squared_i, Y.outer_squared);
            const auto ci = Y.annulus ? cmp(Y.inner_squared_i, Y.inner_squared) : CGAL::LARGER;
            if (co == CGAL::EQUAL || ci == CGAL::EQUAL) return true;
            in_outer[t] = co == CGAL::SMALLER;
            in_inner[t] = Y.annulus && ci == CGAL::SMALLER;
        }
        for (const std::uint32_t tc : table.target_circles[k]) {
            // A target circle identical to C (a capsule cap at a domain vertex)
            // contains every witness; it has no events on C and is the target
            // decision's coincident support, not a cut.
            if (std::find(same.begin(), same.end(), tc) != same.end()) continue;
            const BoundaryCircle& T = set.circles[tc];
            if (CGAL::compare(CGAL::squared_distance(w, T.circle.center()), T.circle.squared_radius()) == CGAL::EQUAL) {
                return true;
            }
        }
        for (std::size_t t = 0; t < sweep.lines.size(); ++t) {
            line_state[t] = line_side(lines[sweep.lines[t]], w, wx, wy, mode, counters);
            if (line_state[t] == CGAL::ZERO) return true;
        }
        for (const std::uint32_t v : sweep.vertices) {
            const VertexCircle& V = events.vertex(v);
            const Interval dx = wx - V.x_i, dy = wy - V.y_i, d2 = dx * dx + dy * dy;
            if (mode == FilterMode::Lazy && (d2.sup() < V.r_sq_i.inf() || d2.inf() > V.r_sq_i.sup())) continue;
            ++counters.exact_target_events;
            if (vertex_circle_side(Tower::rational(w.x()), Tower::rational(w.y()), V) == CGAL::ZERO) return true;
        }
        return false;
    };
    const StartWitness start = find_start_witness(
        C, start_witness_attempt_bound(sweep.counted.size(), cut_circles + sweep.vertices.size(), sweep.lines.size()),
        on_support, counters);
    if (start.attempt > 0) {
        // Roll back to the wrap edge: undo every transversal flip strictly
        // before w in the CCW order.
        const Point_2& w = start.point;
        const Interval wx = interval_of(w.x());
        for (std::uint32_t i = 0; i < evs.size(); ++i) {
            const auto pos = order_event_point(events, evs[i], sweep.halves[i], w, wx, start.half, counters);
            if (pos == CGAL::EQUAL) throw CoverageCertificateNativeError("start witness coincides with an event");
            if (pos != CGAL::SMALLER || evs[i].multiplicity != 1 || evs[i].source == 2) continue;
            if (evs[i].source == 1) {
                auto& state = line_state[slot_of(sweep.lines, evs[i].other)];
                state = CGAL::opposite(state);
                continue;
            }
            const BoundaryCircle& O = set.circles[evs[i].other];
            if (O.is_target()) continue;
            const std::size_t t = slot_of(sweep.counted, O.owner);
            if (O.is_inner()) {
                in_inner[t] ^= 1U;
            } else {
                in_outer[t] ^= 1U;
            }
        }
    }
    for (std::uint32_t r = 0; r < sweep.rects.size(); ++r) {
        const RectRegion& R = rects[sweep.rects[r]];
        for (std::uint8_t s = 0; s < 4; ++s) rect_good[r] += line_state[slot_of(sweep.lines, R.lines[s])] == R.inward[s] ? 1 : 0;
    }

    int covered = 0;
    for (std::size_t t = 0; t < sweep.counted.size(); ++t) covered += region_covers(t) ? 1 : 0;
    for (std::uint32_t r = 0; r < sweep.rects.size(); ++r) covered += rect_covers(r) ? 1 : 0;
    const std::size_t vertices = sweep.vertex_count();
    if (vertices == 0) {
        ++counters.sub_arcs;
        sweep.edge_exposed.assign(1, covered == 0 ? 1U : 0U);
        return sweep;
    }
    sweep.edge_exposed.assign(vertices, 0U);
    const int covered_at_start = covered;
    for (std::size_t v = 0; v < vertices; ++v) {
        for (std::uint32_t i = sweep.group_starts[v]; i < sweep.group_starts[v + 1]; ++i) {
            if (evs[i].multiplicity != 1 || evs[i].source == 2) continue;
            if (evs[i].source == 1) {
                const std::size_t ls = slot_of(sweep.lines, evs[i].other);
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
            const BoundaryCircle& O = set.circles[evs[i].other];
            if (O.is_target()) continue;
            const std::size_t t = slot_of(sweep.counted, O.owner);
            const int before = region_covers(t) ? 1 : 0;
            if (O.is_inner()) {
                in_inner[t] ^= 1U;
            } else {
                in_outer[t] ^= 1U;
            }
            covered += (region_covers(t) ? 1 : 0) - before;
        }
        sweep.edge_exposed[v] = covered == 0 ? 1U : 0U;
        ++counters.sub_arcs;
    }
    if (covered != covered_at_start) {
        throw CoverageCertificateNativeError("cover counters fail parity around boundary circle " + std::to_string(k));
    }
    return sweep;
}

}  // namespace compas_cgal::coverage_certificate
