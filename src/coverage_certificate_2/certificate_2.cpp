#include "coverage_certificate_2/certificate_2.h"

#include "coverage_certificate_2/circle_events_2.h"
#include "coverage_certificate_2/circle_sweep_2.h"
#include "coverage_certificate_2/driver_support_2.h"
#include "coverage_certificate_2/errors_2.h"
#include "coverage_certificate_2/lazy_crossings_2.h"
#include "coverage_certificate_2/line_sweep_2.h"
#include "coverage_certificate_2/line_witness_2.h"
#include "coverage_certificate_2/target_decision_2.h"
#include "coverage_certificate_2/witness_2.h"

#include <CGAL/box_intersection_d.h>

#include <algorithm>
#include <chrono>
#include <utility>

namespace compas_cgal::coverage_certificate {

namespace {
using Clock = std::chrono::steady_clock;

double seconds(const Clock::time_point a, const Clock::time_point b) { return std::chrono::duration<double>(b - a).count(); }

std::vector<std::uint32_t> vertex_circles(const CircleSweep& sweep, const std::size_t vertex)
{
    std::vector<std::uint32_t> ids;
    for (std::uint32_t i = sweep.group_starts[vertex]; i < sweep.group_starts[vertex + 1]; ++i) {
        ids.push_back(sweep.events[i].other);
    }
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    return ids;
}

}  // namespace

CircleCoreResult run_circle_core(const CircleRegionSet& set, const MachinableTarget& target, const FilterMode mode,
                                 const bool audit)
{
    const auto start = Clock::now();
    CircleCoreResult result;
    result.region_count = set.regions.size();
    result.piece_count = set.circles.size();
    const LineTable lines(target.domain, target.tool_radius);
    const std::vector<RectRegion> rects;
    CrossingTable table = CrossingTable::build(set, mode, result.counters);
    table.verify_enclosures(audit);
    CircleEvents events(set, lines, table);
    const CircleCandidates none;
    CircleWitnessContext context{set, lines, rects, events};
    const auto accept_all = [](const Point_2&) { return true; };
    result.crossings_s = seconds(start, Clock::now());
    for (std::uint32_t k = 0; k < set.circles.size(); ++k) {
        const auto t0 = Clock::now();
        const CircleSweep sweep = sweep_circle(k, set, lines, rects, none, events, result.counters);
        result.sweep_s += seconds(t0, Clock::now());
        if (sweep.covered_by_coincident_support) continue;
        const std::size_t vertices = sweep.vertex_count();
        for (std::size_t edge = 0; edge < sweep.edge_exposed.size(); ++edge) {
            if (sweep.edge_exposed[edge] == 0U) {
                if (audit) audit_covered_edge(k, edge, sweep, context, result.counters);
                continue;
            }
            const auto t1 = Clock::now();
            const auto [s, w] = confirm_exposed_edge(k, edge, sweep, context, result.counters);
            const Point_2 u = uncovered_beside(k, w, sweep, context, accept_all, result.counters);
            const auto t2 = Clock::now();
            const bool in_target = target.contains(u, mode, result.counters);
            result.witness_s += seconds(t1, t2);
            result.target_s += seconds(t2, Clock::now());
            ExposedEdge record{k,
                               set.circles[k].owner,
                               set.circles[k].kind,
                               static_cast<std::uint32_t>(edge),
                               vertices == 0,
                               vertices == 0 ? std::vector<std::uint32_t>{} : vertex_circles(sweep, edge),
                               vertices == 0 ? std::vector<std::uint32_t>{} : vertex_circles(sweep, (edge + 1) % vertices),
                               s,
                               w,
                               u,
                               in_target};
            result.exposed.push_back(std::move(record));
        }
    }
    result.elapsed_s = seconds(start, Clock::now());
    return result;
}

CertificateResult run_certificate(CircleRegionSet set, const MachinableTarget& target,
                                  const std::vector<ConnectorSpec>& connectors, const FilterMode mode, const bool audit)
{
    const auto start = Clock::now();
    CertificateResult result;
    DecisionCounters& counters = result.counters;
    const auto n_full = static_cast<std::uint32_t>(set.regions.size());
    const auto n_conn = static_cast<std::uint32_t>(connectors.size());
    const FT& r = target.tool_radius;
    const FT r_sq = r * r;

    // Regions: cuts, capsule end disks, capsule rectangles; target members.
    LineTable lines(target.domain, r);
    const ConnectorCapsules caps = ConnectorCapsules::build(target.domain, lines, connectors);
    std::vector<std::uint32_t> cap_report;
    cap_report.reserve(caps.disk_report_index.size());
    for (const std::uint32_t c : caps.disk_report_index) cap_report.push_back(n_full + c);
    set.add_capsule_caps(caps.disk_centres, cap_report, r);
    std::vector<RectRegion> target_rects;
    std::vector<Point_2> vertices;
    std::vector<std::uint32_t> members;
    for (std::uint32_t e = 0; e < lines.edge_count(); ++e) {
        target_rects.push_back(make_rect(lines, e, FT(0), FT(1), true, 1 + e));
        vertices.push_back(lines.edge_start(e));
        members.push_back(1 + e);
    }
    set.add_target_circles(vertices, members, r);
    const std::vector<RectRegion>& rects = caps.rects;
    result.region_count = set.regions.size() + rects.size();

    // Flank pieces: two per swept rectangle (checked outward), two per target
    // rectangle (checked toward the edge).
    std::vector<LinePiece> pieces;
    for (std::uint32_t i = 0; i < rects.size(); ++i) {
        const RectRegion& R = rects[i];
        for (const int sigma : {1, -1}) {
            const std::uint32_t id = lines.flank(R.edge, sigma);
            pieces.push_back({id, R.tau0, R.tau1, false, i, R.report_index, sigma > 0 ? CGAL::POSITIVE : CGAL::NEGATIVE,
                              flank_box(lines[id], R.tau0, R.tau1)});
        }
    }
    for (std::uint32_t i = 0; i < target_rects.size(); ++i) {
        const RectRegion& R = target_rects[i];
        for (const int sigma : {1, -1}) {
            const std::uint32_t id = lines.flank(R.edge, sigma);
            pieces.push_back({id, FT(0), FT(1), true, i, R.report_index, sigma > 0 ? CGAL::NEGATIVE : CGAL::POSITIVE,
                              flank_box(lines[id], FT(0), FT(1))});
        }
    }
    const auto n_circles = static_cast<std::uint32_t>(set.circles.size());
    result.piece_count = n_circles + pieces.size();

    // Box index: circle-circle pairs inside the CK table, then the bipartite
    // overlaps of pieces with rectangles and circles.
    CrossingTable table = CrossingTable::build(set, mode, counters);
    table.verify_enclosures(audit);
    std::vector<CGAL::Bbox_2> circle_boxes, rect_boxes, target_boxes, piece_boxes;
    for (const auto& c : set.circles) circle_boxes.push_back(c.bbox);
    for (const auto& R : rects) rect_boxes.push_back(R.bbox);
    for (const auto& R : target_rects) target_boxes.push_back(R.bbox);
    for (const auto& p : pieces) piece_boxes.push_back(p.box);
    const auto circle_rects = overlaps(circle_boxes, rect_boxes);
    const auto circle_targets = overlaps(circle_boxes, target_boxes);
    const auto piece_circles = overlaps(piece_boxes, circle_boxes);
    const auto piece_rects = overlaps(piece_boxes, rect_boxes);
    const auto piece_targets = overlaps(piece_boxes, target_boxes);

    CircleEvents events(set, lines, table);
    std::vector<CircleCandidates> circle_candidates(n_circles);
    std::vector<TargetCandidates> circle_target_candidates(n_circles);
    for (std::uint32_t k = 0; k < n_circles; ++k) {
        CircleCandidates& cc = circle_candidates[k];
        cc.rects = circle_rects[k];
        for (const std::uint32_t i : cc.rects) cc.lines.insert(cc.lines.end(), rects[i].lines.begin(), rects[i].lines.end());
        if (!set.circles[k].is_target()) {
            for (const std::uint32_t i : circle_targets[k]) {
                cc.lines.push_back(target_rects[i].lines[0]);
                cc.lines.push_back(target_rects[i].lines[1]);
            }
            circle_target_candidates[k] = {circle_targets[k], table.target_circles[k]};
        }
        std::sort(cc.lines.begin(), cc.lines.end());
        cc.lines.erase(std::unique(cc.lines.begin(), cc.lines.end()), cc.lines.end());
        for (const std::uint32_t l : cc.lines) events.add_line(k, l, counters);
    }
    std::vector<LineCandidates> line_candidates(pieces.size());
    std::vector<TargetCandidates> line_target_candidates(pieces.size());
    for (std::uint32_t p = 0; p < pieces.size(); ++p) {
        const LinePiece& piece = pieces[p];
        LineCandidates& lc = line_candidates[p];
        TargetCandidates& tc = line_target_candidates[p];
        for (const std::uint32_t c : piece_circles[p]) {
            if (set.circles[c].is_target()) {
                if (!piece.target) tc.circles.push_back(c);
            } else {
                lc.regions.push_back(set.circles[c].owner);
            }
        }
        std::sort(lc.regions.begin(), lc.regions.end());
        lc.regions.erase(std::unique(lc.regions.begin(), lc.regions.end()), lc.regions.end());
        for (const std::uint32_t i : piece_rects[p]) {
            if (!piece.target && i == piece.owner) continue;
            lc.rects.push_back(i);
        }
        if (!piece.target) {
            lc.target_circles = tc.circles;
            tc.rects = piece_targets[p];
            for (const std::uint32_t i : tc.rects) {
                lc.target_lines.push_back(target_rects[i].lines[0]);
                lc.target_lines.push_back(target_rects[i].lines[1]);
            }
            std::sort(lc.target_lines.begin(), lc.target_lines.end());
            lc.target_lines.erase(std::unique(lc.target_lines.begin(), lc.target_lines.end()), lc.target_lines.end());
        }
    }
    result.crossings_s = seconds(start, Clock::now());

    const TargetModel model{&target, &target_rects, &lines, &set};
    CircleWitnessContext circle_context{set, lines, rects, events};
    LineWitnessContext line_context{set, lines, rects, mode};

    // Circle pieces.
    for (std::uint32_t k = 0; k < n_circles; ++k) {
        const BoundaryCircle& C = set.circles[k];
        const auto t0 = Clock::now();
        const CircleSweep sweep = sweep_circle(k, set, lines, rects, circle_candidates[k], events, counters);
        const auto t1 = Clock::now();
        result.circle_sweep_s += seconds(t0, t1);
        if (C.is_target()) ++counters.target_pieces;
        if (sweep.covered_by_coincident_support) {
            ++result.edge_count;
            continue;
        }
        result.edge_count += sweep.edge_exposed.size();
        for (std::size_t edge = 0; edge < sweep.edge_exposed.size(); ++edge) {
            if (sweep.edge_exposed[edge] == 0U) {
                if (audit) audit_covered_edge(k, edge, sweep, circle_context, counters);
                continue;
            }
            // P3, then the target side (cheap rejection of edges outside T),
            // then P4 on the edges that matter.
            const FT s = edge_parameter(k, edge, sweep, circle_context, counters);
            const Point_2 w = stereographic_point(C, s);
            CertificatePieceKind kind = circle_kind(C.kind);
            std::uint32_t region_index = 0;
            std::function<bool(const Point_2&)> accept;
            if (C.is_target()) {
                region_index = n_full + n_conn + C.owner;
                accept = [&](const Point_2& u) {
                    return CGAL::compare(CGAL::squared_distance(u, C.circle.center()), r_sq) != CGAL::LARGER;
                };
            } else {
                const TargetSide side = decide_checked_side(model, circle_target_candidates[k], {false, k, CGAL::ZERO},
                                                            ExactPoint::rational(w), mode, counters);
                if (side == TargetSide::Outside) {
                    ++counters.sweep_exposed_outside_target;
                    continue;
                }
                if (side == TargetSide::Undecided) {
                    ++result.undecided;
                    continue;
                }
                region_index = set.regions[C.owner].report_index;
                accept = [&](const Point_2& u) {
                    return target_contains(model, circle_target_candidates[k], u, mode, counters);
                };
            }
            confirm_circle_witness(k, w, sweep, circle_context, counters);
            const Point_2 u = uncovered_beside(k, w, sweep, circle_context, accept, counters);
            result.exposed.push_back({kind, region_index, k, static_cast<std::uint32_t>(edge), s, u});
        }
        result.witness_s += seconds(t1, Clock::now());
    }

    // Flank pieces.
    for (std::uint32_t p = 0; p < pieces.size(); ++p) {
        const LinePiece& piece = pieces[p];
        const auto t0 = Clock::now();
        LineSweep sweep = sweep_line(piece, set, lines, rects, line_candidates[p], mode, counters);
        const auto t1 = Clock::now();
        result.line_sweep_s += seconds(t0, t1);
        if (piece.target) ++counters.target_pieces;
        result.edge_count += sweep.edge_exposed.size();
        for (std::size_t edge = 0; edge < sweep.edge_exposed.size(); ++edge) {
            if (sweep.edge_exposed[edge] == 0U) {
                if (audit) audit_covered_line_edge(sweep, edge, line_context, counters);
                continue;
            }
            const FT t = line_edge_parameter(sweep, edge, counters);
            CertificatePieceKind kind = CertificatePieceKind::CapsuleFlank;
            std::uint32_t region_index = 0;
            std::function<bool(const Point_2&)> accept;
            if (piece.target) {
                kind = CertificatePieceKind::TargetBoundary;
                region_index = n_full + n_conn + piece.report_index;
                const RectRegion& own = target_rects[piece.owner];
                accept = [&lines, &own](const Point_2& u) {
                    return rect_membership(own, lines, ExactPoint::rational(u)) != RectSide::Exterior;
                };
            } else {
                const ExactPoint w = line_point(lines[piece.line], t);
                const TargetSide side = decide_checked_side(model, line_target_candidates[p],
                                                            {true, piece.line, piece.check_side}, w, mode, counters);
                if (side == TargetSide::Outside) {
                    ++counters.sweep_exposed_outside_target;
                    continue;
                }
                if (side == TargetSide::Undecided) {
                    ++result.undecided;
                    continue;
                }
                region_index = n_full + piece.report_index;
                accept = [&](const Point_2& u) {
                    return target_contains(model, line_target_candidates[p], u, mode, counters);
                };
            }
            confirm_line_witness(sweep, t, line_context, counters);
            const Point_2 u = uncovered_beside_line(sweep, t, line_context, accept, counters);
            result.exposed.push_back({kind, region_index, n_circles + p, static_cast<std::uint32_t>(edge), t, u});
        }
        result.witness_s += seconds(t1, Clock::now());
    }
    result.elapsed_s = seconds(start, Clock::now());
    return result;
}

}  // namespace compas_cgal::coverage_certificate
