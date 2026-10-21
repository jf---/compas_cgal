#include "coverage_certificate_2/exact_certificate_2.h"

#include "coverage_certificate_2/circle_events_2.h"
#include "coverage_certificate_2/circle_sweep_2.h"
#include "coverage_certificate_2/driver_support_2.h"
#include "coverage_certificate_2/errors_2.h"
#include "coverage_certificate_2/lazy_crossings_2.h"
#include "coverage_certificate_2/line_sweep_2.h"
#include "coverage_certificate_2/line_witness_2.h"
#include "coverage_certificate_2/vertex_events_2.h"
#include "coverage_certificate_2/witness_2.h"

#include <algorithm>
#include <chrono>
#include <functional>
#include <string>

namespace compas_cgal::coverage_certificate {

namespace {
using Clock = std::chrono::steady_clock;

double seconds(const Clock::time_point a, const Clock::time_point b) { return std::chrono::duration<double>(b - a).count(); }

void sort_unique(std::vector<std::uint32_t>& ids)
{
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
}

/// The vertex circles with an irrational centre (rational ones are ordinary
/// target circles).
VertexCircle vertex_circle(const VertexDisk& disk, const std::uint32_t index, const FT& r)
{
    VertexCircle V;
    V.x = disk.x;
    V.y = disk.y;
    V.delta = disk.x.is_extended() ? disk.x.root() : disk.y.root();
    V.r_sq = r * r;
    V.x_i = disk.x_i;
    V.y_i = disk.y_i;
    V.r_i = interval_of(r);
    V.r_sq_i = interval_of(V.r_sq);
    V.bbox = CGAL::Bbox_2((V.x_i - V.r_i).inf(), (V.y_i - V.r_i).inf(), (V.x_i + V.r_i).sup(), (V.y_i + V.r_i).sup());
    V.disk = index;
    return V;
}

/// Exact closed membership of the rational point p in the union of the swept
/// regions (pure exact P4, no filter).
bool swept(const Point_2& p, const CircleRegionSet& set, const std::vector<RectRegion>& rects, const LineTable& lines,
           DecisionCounters& counters)
{
    const Interval px = interval_of(p.x()), py = interval_of(p.y());
    for (const auto& region : set.regions) {
        if (membership(region, p, px, py, FilterMode::Lazy, counters) != Side::Exterior) return true;
    }
    for (const auto& rect : rects) {
        if (rect_membership(rect, lines, ExactPoint::rational(p)) != RectSide::Exterior) return true;
    }
    return false;
}
}  // namespace

CertificateResult run_exact_target_certificate(CircleRegionSet set, const CentreDomain& path_domain,
                                               const std::vector<ConnectorSpec>& connectors, const ExactTarget& target,
                                               const FilterMode mode, const bool audit)
{
    const auto start = Clock::now();
    CertificateResult result;
    DecisionCounters& counters = result.counters;
    const auto n_full = static_cast<std::uint32_t>(set.regions.size());
    const auto n_conn = static_cast<std::uint32_t>(connectors.size());
    const FT& r = target.tool_radius;

    // Swept regions: cuts, connector end disks and rectangles (on the path domain).
    LineTable lines(path_domain, r);
    const ConnectorCapsules caps = ConnectorCapsules::build(path_domain, lines, connectors);
    std::vector<std::uint32_t> cap_report;
    cap_report.reserve(caps.disk_report_index.size());
    for (const std::uint32_t c : caps.disk_report_index) cap_report.push_back(n_full + c);
    set.add_capsule_caps(caps.disk_centres, cap_report, r);
    const std::vector<RectRegion>& rects = caps.rects;
    result.region_count = set.regions.size() + rects.size();
    const auto n_swept_circles = static_cast<std::uint32_t>(set.circles.size());

    // Target cut supports: declared lines and circles, vertex circles.
    std::vector<std::uint32_t> segment_lines;
    std::vector<CGAL::Bbox_2> segment_boxes;
    for (const auto& s : target.segments) {
        segment_lines.push_back(lines.add_rational_line(-s.uy, s.ux, s.a));
        segment_boxes.push_back(s.bbox);
    }
    const std::uint32_t member = n_full + n_conn;
    for (const auto& arc : target.arcs) (void)set.add_target_circle(arc.c, arc.rho, member);
    std::vector<VertexCircle> vertices;
    std::vector<CGAL::Bbox_2> vertex_boxes;
    for (std::uint32_t i = 0; i < target.vertex_disks.size(); ++i) {
        const VertexDisk& disk = target.vertex_disks[i];
        if (!disk.x.is_extended() && !disk.y.is_extended()) {
            (void)set.add_target_circle(Point_2(disk.x.a0(), disk.y.a0()), r, member);
            continue;
        }
        vertices.push_back(vertex_circle(disk, i, r));
        vertex_boxes.push_back(vertices.back().bbox);
    }

    // Flank pieces: two per swept rectangle, checked outward.
    std::vector<LinePiece> pieces;
    for (std::uint32_t i = 0; i < rects.size(); ++i) {
        const RectRegion& R = rects[i];
        for (const int sigma : {1, -1}) {
            const std::uint32_t id = lines.flank(R.edge, sigma);
            pieces.push_back({id, R.tau0, R.tau1, false, i, R.report_index, sigma > 0 ? CGAL::POSITIVE : CGAL::NEGATIVE,
                              flank_box(lines[id], R.tau0, R.tau1)});
        }
    }
    result.piece_count = n_swept_circles + pieces.size();

    CrossingTable table = CrossingTable::build(set, mode, counters);
    table.verify_enclosures(audit);
    std::vector<CGAL::Bbox_2> circle_boxes, rect_boxes, piece_boxes;
    for (const auto& c : set.circles) circle_boxes.push_back(c.bbox);
    for (const auto& R : rects) rect_boxes.push_back(R.bbox);
    for (const auto& p : pieces) piece_boxes.push_back(p.box);
    const auto circle_rects = overlaps(circle_boxes, rect_boxes);
    const auto circle_segments = overlaps(circle_boxes, segment_boxes);
    const auto circle_vertices = overlaps(circle_boxes, vertex_boxes);
    const auto piece_circles = overlaps(piece_boxes, circle_boxes);
    const auto piece_rects = overlaps(piece_boxes, rect_boxes);
    const auto piece_segments = overlaps(piece_boxes, segment_boxes);
    const auto piece_vertices = overlaps(piece_boxes, vertex_boxes);

    CircleEvents events(set, lines, table, &vertices);
    std::vector<CircleCandidates> circle_candidates(n_swept_circles);
    for (std::uint32_t k = 0; k < n_swept_circles; ++k) {
        CircleCandidates& cc = circle_candidates[k];
        cc.rects = circle_rects[k];
        for (const std::uint32_t i : cc.rects) cc.lines.insert(cc.lines.end(), rects[i].lines.begin(), rects[i].lines.end());
        for (const std::uint32_t i : circle_segments[k]) cc.lines.push_back(segment_lines[i]);
        sort_unique(cc.lines);
        for (const std::uint32_t l : cc.lines) events.add_line(k, l, counters);
        cc.vertices = circle_vertices[k];
        for (const std::uint32_t v : cc.vertices) events.add_vertex(k, v, counters);
    }
    std::vector<LineCandidates> line_candidates(pieces.size());
    for (std::uint32_t p = 0; p < pieces.size(); ++p) {
        const LinePiece& piece = pieces[p];
        LineCandidates& lc = line_candidates[p];
        for (const std::uint32_t c : piece_circles[p]) {
            if (set.circles[c].is_target()) {
                lc.target_circles.push_back(c);
            } else {
                lc.regions.push_back(set.circles[c].owner);
            }
        }
        sort_unique(lc.regions);
        for (const std::uint32_t i : piece_rects[p]) {
            if (i != piece.owner) lc.rects.push_back(i);
        }
        for (const std::uint32_t i : piece_segments[p]) lc.target_lines.push_back(segment_lines[i]);
        sort_unique(lc.target_lines);
        lc.target_vertices = piece_vertices[p];
        for (const std::uint32_t l : lc.target_lines) {
            if (line_coincidence(lines[piece.line], lines[l]) != 0) {
                throw CoincidentTargetSupportError("a connector flank lies on a declared boundary line of the exact target");
            }
        }
    }
    result.crossings_s = seconds(start, Clock::now());

    CircleWitnessContext circle_context{set, lines, rects, events};
    LineWitnessContext line_context{set, lines, rects, mode};
    const auto in_target = [&](const Point_2& u) { return target.contains(ExactPoint::rational(u), mode, counters); };

    // Swept circles (target circles are cut supports only).
    for (std::uint32_t k = 0; k < n_swept_circles; ++k) {
        const BoundaryCircle& C = set.circles[k];
        for (const std::uint32_t other : table.coincident[k]) {
            if (set.circles[other].is_target()) {
                throw CoincidentTargetSupportError("swept circle " + std::to_string(k) +
                                                   " lies on a declared boundary circle or vertex circle of the exact target");
            }
        }
        const auto t0 = Clock::now();
        const CircleSweep sweep = sweep_circle(k, set, lines, rects, circle_candidates[k], events, counters);
        const auto t1 = Clock::now();
        result.circle_sweep_s += seconds(t0, t1);
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
            const FT s = edge_parameter(k, edge, sweep, circle_context, counters);
            const Point_2 w = stereographic_point(C, s);
            ++counters.target_decisions;
            if (!in_target(w)) {
                ++counters.sweep_exposed_outside_target;
                continue;
            }
            confirm_circle_witness(k, w, sweep, circle_context, counters);
            const Point_2 u = uncovered_beside(k, w, sweep, circle_context, in_target, counters);
            result.exposed.push_back(
                {circle_kind(C.kind), set.regions[C.owner].report_index, k, static_cast<std::uint32_t>(edge), s, u});
        }
        result.witness_s += seconds(t1, Clock::now());
    }

    // Connector flanks.
    for (std::uint32_t p = 0; p < pieces.size(); ++p) {
        const LinePiece& piece = pieces[p];
        const auto t0 = Clock::now();
        LineSweep sweep = sweep_line(piece, set, lines, rects, line_candidates[p], mode, counters, &vertices);
        const auto t1 = Clock::now();
        result.line_sweep_s += seconds(t0, t1);
        result.edge_count += sweep.edge_exposed.size();
        for (std::size_t edge = 0; edge < sweep.edge_exposed.size(); ++edge) {
            if (sweep.edge_exposed[edge] == 0U) {
                if (audit) audit_covered_line_edge(sweep, edge, line_context, counters);
                continue;
            }
            const FT t = line_edge_parameter(sweep, edge, counters);
            ++counters.target_decisions;
            if (!target.contains(line_point(lines[piece.line], t), mode, counters)) {
                ++counters.sweep_exposed_outside_target;
                continue;
            }
            confirm_line_witness(sweep, t, line_context, counters);
            const Point_2 u = uncovered_beside_line(sweep, t, line_context, in_target, counters);
            result.exposed.push_back({CertificatePieceKind::CapsuleFlank, n_full + piece.report_index, n_swept_circles + p,
                                      static_cast<std::uint32_t>(edge), t, u});
        }
        result.witness_s += seconds(t1, Clock::now());
    }

    // The interior anchor: int T is connected, so with no exposed swept edge
    // T is covered iff one interior point is.
    const Point_2 anchor = target.interior_anchor(counters);
    if (!swept(anchor, set, rects, lines, counters)) {
        result.exposed.push_back({CertificatePieceKind::TargetInterior, member, static_cast<std::uint32_t>(result.piece_count), 0,
                                  FT(0), anchor});
    }
    result.elapsed_s = seconds(start, Clock::now());
    return result;
}

}  // namespace compas_cgal::coverage_certificate
