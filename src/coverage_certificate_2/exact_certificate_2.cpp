#include "coverage_certificate_2/exact_certificate_2.h"

#include "coverage_certificate_2/circle_events_2.h"
#include "coverage_certificate_2/circle_sweep_2.h"
#include "coverage_certificate_2/driver_support_2.h"
#include "coverage_certificate_2/errors_2.h"
#include "coverage_certificate_2/exact_geometry_2.h"
#include "coverage_certificate_2/lazy_crossings_2.h"
#include "coverage_certificate_2/line_sweep_2.h"
#include "coverage_certificate_2/line_witness_2.h"
#include "coverage_certificate_2/side_probe_2.h"
#include "coverage_certificate_2/vertex_events_2.h"
#include "coverage_certificate_2/witness_2.h"

#include <algorithm>
#include <chrono>
#include <functional>
#include <map>
#include <string>
#include <utility>

namespace compas_cgal::coverage_certificate {

namespace {
using Clock = std::chrono::steady_clock;

double seconds(const Clock::time_point a, const Clock::time_point b) { return std::chrono::duration<double>(b - a).count(); }

/// piece_connector entry of a polygon capsule flank (it bounds no connector region).
constexpr std::int32_t NO_CONNECTOR_REGION = -1;

void sort_unique(std::vector<std::uint32_t>& ids)
{
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
}

bool boxes_meet(const CGAL::Bbox_2& a, const CGAL::Bbox_2& b)
{
    return !(a.xmax() < b.xmin() || b.xmax() < a.xmin() || a.ymax() < b.ymin() || b.ymax() < a.ymin());
}

/// Certified box of the probe segment [w, w + dir].
CGAL::Bbox_2 segment_box(const Point_2& w, const FT& dx, const FT& dy)
{
    const Interval x0 = interval_of(w.x()), y0 = interval_of(w.y());
    const Interval x1 = x0 + interval_of(dx), y1 = y0 + interval_of(dy);
    return CGAL::Bbox_2(std::min(x0.inf(), x1.inf()), std::min(y0.inf(), y1.inf()), std::max(x0.sup(), x1.sup()),
                        std::max(y0.sup(), y1.sup()));
}

/// The exact target's cut supports, for coincidence tests and side probes.
struct TargetSupports {
    std::vector<std::uint32_t> lines;     ///< declared segment lines
    std::vector<CGAL::Bbox_2> line_boxes; ///< certified boxes of the declared segments
    std::vector<std::uint32_t> circles;   ///< target circle ids in the set (declared arcs, rational vertex disks)
    std::vector<std::uint32_t> vertices;  ///< target vertex circles (irrational centres)
};

/// Exact closed membership of the rational point p in the union of the swept
/// regions (pure exact P4, no filter): cuts, connector capsules, exact connectors.
bool swept(const Point_2& p, const CircleRegionSet& set, const std::vector<RectRegion>& rects, const LineTable& lines,
           const ConnectorSet& connectors, const FT& r, DecisionCounters& counters)
{
    const Interval px = interval_of(p.x()), py = interval_of(p.y());
    for (const auto& region : set.regions) {
        if (!region.covers) continue;
        if (membership(region, p, px, py, FilterMode::Lazy, counters) != Side::Exterior) return true;
    }
    for (const auto& rect : rects) {
        if (rect_membership(rect, lines, ExactPoint::rational(p)) != RectSide::Exterior) return true;
    }
    for (const auto& region : connectors.regions) {
        if (connector_membership(region, ExactPoint::rational(p), r, FilterMode::Lazy, counters) != Side::Exterior) return true;
    }
    return false;
}

/// Whether edge `edge` of a connector support circle is certified (by interval
/// enclosures) to lie outside the closed wedge of every sub-piece it could bound.
/// The support circle is cut at every wedge boundary point (the tangency with the
/// sub-piece's cap circle on the radial line), so an open edge lies wholly inside
/// an open wedge or wholly outside the closed wedge, and one endpoint decides it.
/// A certified decision only: inconclusive enclosures fall through to the exact
/// witness test (and EXACT_ONLY never takes this path).
bool certified_out_of_range(const CircleSweep& sweep, const std::size_t edge, const std::vector<std::uint32_t>& owners,
                            const ConnectorSet& connectors)
{
    const std::size_t vertices = sweep.vertex_count();
    if (vertices == 0) return false;
    const CrossingEvent& a = sweep.events[sweep.group_starts[edge]];
    const CrossingEvent& b = sweep.events[sweep.group_starts[(edge + 1) % vertices]];
    for (const std::uint32_t i : owners) {
        const ConnectorRegion& R = connectors.regions[i];
        if (R.degenerate) continue;
        if (R.full) return false;
        const auto outside = [&](const CrossingEvent& e) {
            const Interval c0 = cross_about_i(R.cx_i, R.cy_i, R.w0.x_i, R.w0.y_i, e.x, e.y);
            const Interval c1 = cross_about_i(R.cx_i, R.cy_i, R.w1.x_i, R.w1.y_i, e.x, e.y);
            // o0 = sense c0, o1 = -sense c1; "certainly negative" per sign.
            const bool o0_neg = R.sense > 0 ? c0.sup() < 0 : c0.inf() > 0;
            const bool o1_neg = R.sense > 0 ? c1.inf() > 0 : c1.sup() < 0;
            switch (R.span) {
                case ArcSpan::Minor: return o0_neg || o1_neg;
                case ArcSpan::Half: return o0_neg;
                case ArcSpan::Major: return o0_neg && o1_neg;
            }
            return false;
        };
        if (!outside(a) && !outside(b)) return false;
    }
    return true;
}

/// Decisions shared by the circle and flank loops.
class ExactDecisions {
public:
    ExactDecisions(const CircleRegionSet& set, const LineTable& lines, const std::vector<VertexCircle>& vertices,
                   const ConnectorSet& connectors, const ExactTarget& target, TargetSupports supports, const FilterMode mode)
        : set_(set), lines_(lines), vertices_(vertices), connectors_(connectors), target_(target),
          supports_(std::move(supports)), mode_(mode)
    {
    }

    /// The supports of connector region `i` (flank lines, support circles, caps).
    [[nodiscard]] ProbeSupports connector_supports(const std::uint32_t i) const
    {
        const ConnectorRegion& R = connectors_.regions[i];
        ProbeSupports s;
        if (!R.degenerate) {
            if (R.is_line) {
                for (const std::uint32_t l : R.flank) s.lines.push_back(&lines_[l]);
            } else {
                const CircleRegion& S = set_.regions[R.support_region];
                s.circles.emplace_back(S.centre, S.outer_squared);
                if (S.annulus) s.circles.emplace_back(S.centre, S.inner_squared);
            }
        }
        for (std::size_t e = 0; e < 2; ++e) {
            if (R.cap_vertex[e] >= 0) s.vertex_circles.push_back(&vertices_[static_cast<std::size_t>(R.cap_vertex[e])]);
            if (R.cap_region[e] >= 0) {
                const CircleRegion& D = set_.regions[static_cast<std::size_t>(R.cap_region[e])];
                s.circles.emplace_back(D.centre, D.outer_squared);
            }
        }
        return s;
    }

    /// Whether a connector region among `candidates` covers the checked side of
    /// an edge at its witness w (dir points to that side): w in its interior, or
    /// w on its boundary (a coincident support) and the probe beside w inside.
    [[nodiscard]] bool connectors_cover(const ExactPoint& w, const FT& dx, const FT& dy,
                                        const std::vector<std::uint32_t>& candidates, DecisionCounters& counters) const
    {
        const Interval wx = interval_of(w.x), wy = interval_of(w.y);
        for (const std::uint32_t i : candidates) {
            if (!may_contain(i, wx, wy)) continue;
            const Side side = connector_membership(connectors_.regions[i], w, target_.tool_radius, mode_, counters);
            if (side == Side::Interior) return true;
            if (side == Side::Exterior) continue;
            const Point_2 q = probe_point(rational(w), dx, dy, connector_supports(i), mode_, counters);
            const Side beside = connector_membership(connectors_.regions[i], ExactPoint::rational(q), target_.tool_radius,
                                                     mode_, counters);
            if (beside == Side::Boundary) {
                throw CoverageCertificateNativeError("a side probe landed on a support of its connector region");
            }
            if (beside == Side::Interior) return true;
        }
        return false;
    }

    /// Whether the point u lies outside every candidate connector region.
    [[nodiscard]] bool outside_connectors(const ExactPoint& u, const std::vector<std::uint32_t>& candidates,
                                          DecisionCounters& counters) const
    {
        const Interval ux = interval_of(u.x), uy = interval_of(u.y);
        for (const std::uint32_t i : candidates) {
            if (!may_contain(i, ux, uy)) continue;
            if (connector_membership(connectors_.regions[i], u, target_.tool_radius, mode_, counters) != Side::Exterior) {
                return false;
            }
        }
        return true;
    }

    /// Whether the checked side of an edge, at its rational witness w (dir to that
    /// side), lies in T. Off every target cut support this is closed membership
    /// at w (w is then not on dT); on one it is membership at the side probe.
    [[nodiscard]] bool checked_side_in_target(const Point_2& w, const FT& dx, const FT& dy, DecisionCounters& counters) const
    {
        const ExactPoint ew = ExactPoint::rational(w);
        const CGAL::Bbox_2 near = segment_box(w, dx, dy);
        const Interval wx = interval_of(w.x()), wy = interval_of(w.y());
        ProbeSupports s;
        bool on_support = false;
        // dT lies on the declared SEGMENTS (and arcs and vertex circles), not on
        // their unrestricted lines: a segment meeting the probe segment has a box
        // meeting the probe segment's box, and off every segment w is not on dT.
        for (std::size_t i = 0; i < supports_.lines.size(); ++i) {
            if (!boxes_meet(supports_.line_boxes[i], near)) continue;
            const LineSupport& line = lines_[supports_.lines[i]];
            if (line_side(line, w, wx, wy, mode_, counters) == CGAL::ZERO) on_support = true;
            s.lines.push_back(&line);
        }
        // A circle meeting the probe segment has a box meeting the segment's box.
        for (const std::uint32_t c : supports_.circles) {
            const BoundaryCircle& C = set_.circles[c];
            if (!boxes_meet(C.bbox, near)) continue;
            if (compare_to_circle(ew, wx, wy, C.circle.center(), C.centre_x, C.centre_y, C.circle.squared_radius(),
                                  C.radius_squared_i, mode_, counters) == CGAL::EQUAL) {
                on_support = true;
            }
            s.circles.emplace_back(C.circle.center(), C.circle.squared_radius());
        }
        for (const std::uint32_t v : supports_.vertices) {
            const VertexCircle& V = vertices_[v];
            if (!boxes_meet(V.bbox, near)) continue;
            const Interval vx = wx - V.x_i, vy = wy - V.y_i, d2 = vx * vx + vy * vy;
            const bool off = mode_ == FilterMode::Lazy && (d2.sup() < V.r_sq_i.inf() || d2.inf() > V.r_sq_i.sup());
            if (!off) {
                ++counters.exact_target_events;
                if (vertex_circle_side(Tower::rational(w.x()), Tower::rational(w.y()), V) == CGAL::ZERO) on_support = true;
            }
            s.vertex_circles.push_back(&V);
        }
        if (!on_support) return target_.contains(ew, mode_, counters);
        ++counters.coincident_target_decisions;
        const Point_2 q = probe_point(w, dx, dy, s, mode_, counters);
        return target_.contains(ExactPoint::rational(q), mode_, counters);
    }

private:
    /// Superset filter: a region lies in its certified box, so a point whose
    /// enclosure misses the box is outside it (and not on its boundary).
    [[nodiscard]] bool may_contain(const std::uint32_t i, const Interval& x, const Interval& y) const
    {
        const CGAL::Bbox_2& b = connectors_.regions[i].bbox;
        return !(x.sup() < b.xmin() || x.inf() > b.xmax() || y.sup() < b.ymin() || y.inf() > b.ymax());
    }

    static Point_2 rational(const ExactPoint& p)
    {
        if (p.x.is_extended() || p.y.is_extended()) {
            throw CoincidentTargetSupportError(
                "a witness on a coincident connector support has irrational coordinates; side probes need a rational witness");
        }
        return Point_2(p.x.a0(), p.y.a0());
    }

    const CircleRegionSet& set_;
    const LineTable& lines_;
    const std::vector<VertexCircle>& vertices_;
    const ConnectorSet& connectors_;
    const ExactTarget& target_;
    TargetSupports supports_;
    FilterMode mode_;
};
}  // namespace

CertificateResult run_exact_target_certificate(CircleRegionSet set, const CentreDomain& path_domain,
                                               const std::vector<ConnectorSpec>& connectors,
                                               const std::vector<ExactConnectorSpec>& exact_connectors,
                                               const ExactTarget& target, const FilterMode mode, const bool audit)
{
    const auto start = Clock::now();
    CertificateResult result;
    DecisionCounters& counters = result.counters;
    const auto n_full = static_cast<std::uint32_t>(set.regions.size());
    const auto n_conn = static_cast<std::uint32_t>(connectors.size());
    const auto n_exact = static_cast<std::uint32_t>(exact_connectors.size());
    const FT& r = target.tool_radius;

    // Swept regions: cuts, connector end disks and rectangles (on the path domain).
    LineTable lines(path_domain, r);
    const ConnectorCapsules caps = ConnectorCapsules::build(path_domain, lines, connectors);
    std::vector<std::uint32_t> cap_report;
    cap_report.reserve(caps.disk_report_index.size());
    for (const std::uint32_t c : caps.disk_report_index) cap_report.push_back(n_full + c);
    set.add_capsule_caps(caps.disk_centres, cap_report, r);
    const std::vector<RectRegion>& rects = caps.rects;

    // Target vertex circles (irrational centres; rational ones are target circles below).
    std::vector<VertexCircle> vertices;
    std::vector<std::uint32_t> rational_disks;
    for (std::uint32_t i = 0; i < target.vertex_disks.size(); ++i) {
        const VertexDisk& disk = target.vertex_disks[i];
        if (!disk.x.is_extended() && !disk.y.is_extended()) {
            rational_disks.push_back(i);
            continue;
        }
        vertices.push_back(make_vertex_circle(disk.x, disk.y, r, i));
    }
    const auto n_target_vertices = static_cast<std::uint32_t>(vertices.size());

    // Exact dOmega* connectors: composite regions; their supports join the set,
    // the line table and the vertex circles (exact_connectors_2.h).
    const ConnectorSet exact = ConnectorSet::build(exact_connectors, target, set, lines, vertices, mode);
    const std::uint32_t exact_base = n_full + n_conn;
    for (std::uint32_t j = 0; j < set.regions.size(); ++j) {
        if (!set.regions[j].covers) set.regions[j].report_index = exact_base + exact.support_reports[j];
    }
    result.region_count = n_full + caps.disk_centres.size() + rects.size() + exact.regions.size();
    const auto n_swept_circles = static_cast<std::uint32_t>(set.circles.size());

    // Target cut supports: declared lines and circles, vertex circles.
    TargetSupports tsupports;
    std::vector<CGAL::Bbox_2> segment_boxes;
    for (const auto& s : target.segments) {
        tsupports.lines.push_back(lines.add_rational_line(-s.uy, s.ux, s.a));
        tsupports.line_boxes.push_back(s.bbox);
        segment_boxes.push_back(s.bbox);
    }
    const std::uint32_t member = exact_base + n_exact;
    for (const auto& arc : target.arcs) tsupports.circles.push_back(set.add_target_circle(arc.c, arc.rho, member));
    for (const std::uint32_t i : rational_disks) {
        const VertexDisk& disk = target.vertex_disks[i];
        tsupports.circles.push_back(set.add_target_circle(Point_2(disk.x.a0(), disk.y.a0()), r, member));
    }
    for (std::uint32_t v = 0; v < n_target_vertices; ++v) tsupports.vertices.push_back(v);
    std::vector<CGAL::Bbox_2> vertex_boxes;
    for (const auto& V : vertices) vertex_boxes.push_back(V.bbox);

    // Flank pieces: two per swept rectangle, and two per line connector sub-piece.
    std::vector<LinePiece> pieces;
    std::vector<std::int32_t> piece_connector;
    for (std::uint32_t i = 0; i < rects.size(); ++i) {
        const RectRegion& R = rects[i];
        for (const int sigma : {1, -1}) {
            const std::uint32_t id = lines.flank(R.edge, sigma);
            pieces.push_back({id, R.tau0, R.tau1, false, i, R.report_index, sigma > 0 ? CGAL::POSITIVE : CGAL::NEGATIVE,
                              flank_box(lines[id], R.tau0, R.tau1)});
            piece_connector.push_back(NO_CONNECTOR_REGION);
        }
    }
    for (std::uint32_t i = 0; i < exact.regions.size(); ++i) {
        const ConnectorRegion& R = exact.regions[i];
        if (!R.is_line || R.degenerate) continue;
        for (std::size_t f = 0; f < 2; ++f) {
            LinePiece piece{R.flank[f], R.t_lo[f], R.t_hi[f], false, i, n_conn + R.connector,
                            f == 0 ? CGAL::POSITIVE : CGAL::NEGATIVE, flank_box(lines[R.flank[f]], R.t_lo[f], R.t_hi[f])};
            piece.owner_is_rect = false;
            pieces.push_back(piece);
            piece_connector.push_back(static_cast<std::int32_t>(i));
        }
    }
    // Reporting only: the lowest connector whose own sub-piece contains a witness.
    std::vector<std::vector<std::uint32_t>> piece_subs(target.piece_bands.size());
    for (std::uint32_t j = 0; j < exact.sub_regions.size(); ++j) piece_subs[exact.sub_regions[j].piece].push_back(j);
    const auto report_connector = [&](const ConnectorRegion& R, const ExactPoint& w) {
        for (const std::uint32_t j : piece_subs[R.piece]) {
            if (connector_in_range(exact.sub_regions[j], w, mode)) return exact.sub_regions[j].connector;
        }
        return R.connector;
    };
    // Reporting only: the line connector regions on each (flank line, side), in
    // index order, so an exposed flank edge is attributed to the lowest one.
    std::map<std::pair<std::uint32_t, std::size_t>, std::vector<std::uint32_t>> flank_regions;
    for (std::uint32_t i = 0; i < exact.regions.size(); ++i) {
        const ConnectorRegion& R = exact.regions[i];
        if (!R.is_line || R.degenerate) continue;
        for (std::size_t f = 0; f < 2; ++f) flank_regions[{R.flank[f], f}].push_back(i);
    }
    // Connector support circles: the pieces each one bounds.
    std::vector<std::vector<std::uint32_t>> circle_connectors(n_swept_circles);
    for (std::uint32_t i = 0; i < exact.regions.size(); ++i) {
        const ConnectorRegion& R = exact.regions[i];
        if (R.is_line || R.degenerate) continue;
        for (std::uint32_t c = set.first_circle[R.support_region]; c < set.first_circle[R.support_region + 1]; ++c) {
            circle_connectors[c].push_back(i);
        }
    }
    std::size_t checked_circles = 0;
    for (std::uint32_t k = 0; k < n_swept_circles; ++k) checked_circles += exact.cap_circle[k] == 0U ? 1 : 0;
    result.piece_count = checked_circles + pieces.size();

    CrossingTable table = CrossingTable::build(set, mode, counters);
    table.verify_enclosures(audit);
    std::vector<CGAL::Bbox_2> circle_boxes, rect_boxes, piece_boxes, connector_boxes;
    for (const auto& c : set.circles) circle_boxes.push_back(c.bbox);
    for (const auto& R : rects) rect_boxes.push_back(R.bbox);
    for (const auto& p : pieces) piece_boxes.push_back(p.box);
    for (const auto& R : exact.regions) connector_boxes.push_back(R.bbox);
    const auto circle_rects = overlaps(circle_boxes, rect_boxes);
    const auto circle_segments = overlaps(circle_boxes, segment_boxes);
    const auto circle_vertices = overlaps(circle_boxes, vertex_boxes);
    const auto circle_exact = overlaps(circle_boxes, connector_boxes);
    const auto piece_circles = overlaps(piece_boxes, circle_boxes);
    const auto piece_rects = overlaps(piece_boxes, rect_boxes);
    const auto piece_segments = overlaps(piece_boxes, segment_boxes);
    const auto piece_vertices = overlaps(piece_boxes, vertex_boxes);
    const auto piece_exact = overlaps(piece_boxes, connector_boxes);

    CircleEvents events(set, lines, table, &vertices);
    std::vector<CircleCandidates> circle_candidates(n_swept_circles);
    for (std::uint32_t k = 0; k < n_swept_circles; ++k) {
        if (exact.cap_circle[k] != 0U) continue;
        CircleCandidates& cc = circle_candidates[k];
        cc.rects = circle_rects[k];
        for (const std::uint32_t i : cc.rects) cc.lines.insert(cc.lines.end(), rects[i].lines.begin(), rects[i].lines.end());
        for (const std::uint32_t i : circle_segments[k]) cc.lines.push_back(tsupports.lines[i]);
        for (const std::uint32_t i : circle_exact[k]) {
            const ConnectorRegion& R = exact.regions[i];
            if (R.is_line && !R.degenerate) cc.lines.insert(cc.lines.end(), R.flank.begin(), R.flank.end());
        }
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
            const BoundaryCircle& C = set.circles[c];
            if (C.is_target() || !set.regions[C.owner].covers) {
                lc.target_circles.push_back(c);
            } else {
                lc.regions.push_back(C.owner);
            }
        }
        sort_unique(lc.regions);
        sort_unique(lc.target_circles);
        for (const std::uint32_t i : piece_rects[p]) {
            if (!piece.owner_is_rect || i != piece.owner) lc.rects.push_back(i);
        }
        for (const std::uint32_t i : piece_segments[p]) lc.target_lines.push_back(tsupports.lines[i]);
        for (const std::uint32_t i : piece_exact[p]) {
            const ConnectorRegion& R = exact.regions[i];
            if (R.is_line && !R.degenerate) lc.target_lines.insert(lc.target_lines.end(), R.flank.begin(), R.flank.end());
        }
        sort_unique(lc.target_lines);
        lc.target_vertices = piece_vertices[p];
    }
    result.crossings_s = seconds(start, Clock::now());

    const ExactDecisions decide(set, lines, vertices, exact, target, std::move(tsupports), mode);
    CircleWitnessContext circle_context{set, lines, rects, events};
    LineWitnessContext line_context{set, lines, rects, mode};

    // Swept circles. Target circles are cut supports only; rational connector cap
    // circles are not pieces (exact_connectors_2.h).
    for (std::uint32_t k = 0; k < n_swept_circles; ++k) {
        if (exact.cap_circle[k] != 0U) continue;
        const BoundaryCircle& C = set.circles[k];
        const auto t0 = Clock::now();
        const CircleSweep sweep = sweep_circle(k, set, lines, rects, circle_candidates[k], events, counters);
        const auto t1 = Clock::now();
        result.circle_sweep_s += seconds(t0, t1);
        const bool connector_circle = !set.regions[C.owner].covers;
        if (connector_circle) ++counters.connector_pieces;
        if (sweep.covered_by_coincident_support) {
            ++result.edge_count;
            continue;
        }
        result.edge_count += sweep.edge_exposed.size();
        const Point_2& c = C.circle.center();
        for (std::size_t edge = 0; edge < sweep.edge_exposed.size(); ++edge) {
            if (sweep.edge_exposed[edge] == 0U) {
                if (audit) audit_covered_edge(k, edge, sweep, circle_context, counters);
                continue;
            }
            if (connector_circle && mode == FilterMode::Lazy && certified_out_of_range(sweep, edge, circle_connectors[k], exact)) {
                ++counters.connector_edges_out_of_range;
                continue;
            }
            const FT s = edge_parameter(k, edge, sweep, circle_context, counters);
            const Point_2 w = stereographic_point(C, s);
            const ExactPoint ew = ExactPoint::rational(w);
            std::uint32_t report = set.regions[C.owner].report_index;
            CertificatePieceKind kind = circle_kind(C.kind);
            if (connector_circle) {
                // A support circle bounds a sub-piece only inside its range.
                const auto& owners = circle_connectors[k];
                const auto in = std::find_if(owners.begin(), owners.end(),
                                             [&](const std::uint32_t i) { return connector_in_range(exact.regions[i], ew, mode); });
                if (in == owners.end()) {
                    ++counters.connector_edges_out_of_range;
                    continue;
                }
                report = exact_base + report_connector(exact.regions[*in], ew);
                kind = CertificatePieceKind::ConnectorArc;
            }
            const FT ox = C.check_inside() ? FT(c.x() - w.x()) : FT(w.x() - c.x());
            const FT oy = C.check_inside() ? FT(c.y() - w.y()) : FT(w.y() - c.y());
            const auto [dx, dy] = probe_direction(ox, oy);
            if (decide.connectors_cover(ew, dx, dy, circle_exact[k], counters)) {
                ++counters.connector_edges_covered;
                continue;
            }
            ++counters.target_decisions;
            if (!decide.checked_side_in_target(w, dx, dy, counters)) {
                ++counters.sweep_exposed_outside_target;
                continue;
            }
            confirm_circle_witness(k, w, sweep, circle_context, counters);
            const auto accept = [&](const Point_2& u) {
                return target.contains(ExactPoint::rational(u), mode, counters) &&
                       decide.outside_connectors(ExactPoint::rational(u), circle_exact[k], counters);
            };
            const Point_2 u = uncovered_beside(k, w, sweep, circle_context, accept, counters);
            result.exposed.push_back({kind, report, k, static_cast<std::uint32_t>(edge), s, u});
        }
        result.witness_s += seconds(t1, Clock::now());
    }

    // Flanks: polygon connector capsule flanks and exact connector flanks.
    for (std::uint32_t p = 0; p < pieces.size(); ++p) {
        const LinePiece& piece = pieces[p];
        const std::int32_t own = piece_connector[p];
        if (own != NO_CONNECTOR_REGION) ++counters.connector_pieces;
        const auto t0 = Clock::now();
        LineSweep sweep = sweep_line(piece, set, lines, rects, line_candidates[p], mode, counters, &vertices);
        const auto t1 = Clock::now();
        result.line_sweep_s += seconds(t0, t1);
        result.edge_count += sweep.edge_exposed.size();
        const LineSupport& support = lines[piece.line];
        for (std::size_t edge = 0; edge < sweep.edge_exposed.size(); ++edge) {
            if (sweep.edge_exposed[edge] == 0U) {
                if (audit) audit_covered_line_edge(sweep, edge, line_context, counters);
                continue;
            }
            const FT t = line_edge_parameter(sweep, edge, counters);
            const ExactPoint ew = line_point(support, t);
            if (own != NO_CONNECTOR_REGION) {
                // A flank bounds its sub-piece only inside the range, and an edge
                // is reported once: by the lowest sub-piece with this flank side.
                const auto i = static_cast<std::uint32_t>(own);
                if (!connector_in_range(exact.regions[i], ew, mode)) {
                    ++counters.connector_edges_out_of_range;
                    continue;
                }
                const std::size_t f = piece.check_side == CGAL::POSITIVE ? 0 : 1;
                bool earlier = false;
                for (const std::uint32_t j : flank_regions.at({piece.line, f})) {
                    if (j >= i || earlier) break;
                    earlier = connector_in_range(exact.regions[j], ew, mode);
                }
                if (earlier) continue;
            }
            if (support.is_rational()) {
                const FT sign = piece.check_side == CGAL::POSITIVE ? FT(1) : FT(-1);
                const auto [dx, dy] = probe_direction(sign * support.nx, sign * support.ny);
                if (decide.connectors_cover(ew, dx, dy, piece_exact[p], counters)) {
                    ++counters.connector_edges_covered;
                    continue;
                }
                ++counters.target_decisions;
                if (!decide.checked_side_in_target(Point_2(ew.x.a0(), ew.y.a0()), dx, dy, counters)) {
                    ++counters.sweep_exposed_outside_target;
                    continue;
                }
            } else {
                // An irrational polygon-connector flank: its witness is a one-root
                // point, which lies on no rational target or connector support
                // unless it shares the support exactly (then it is not modelled).
                bool covered = false;
                for (const std::uint32_t i : piece_exact[p]) {
                    const Side side = connector_membership(exact.regions[i], ew, r, mode, counters);
                    if (side == Side::Boundary) {
                        throw CoincidentTargetSupportError("an irrational connector flank lies on an exact connector support");
                    }
                    covered = covered || side == Side::Interior;
                }
                if (covered) {
                    ++counters.connector_edges_covered;
                    continue;
                }
                ++counters.target_decisions;
                if (!target.contains(ew, mode, counters)) {
                    ++counters.sweep_exposed_outside_target;
                    continue;
                }
            }
            confirm_line_witness(sweep, t, line_context, counters);
            const auto accept = [&](const Point_2& u) {
                return target.contains(ExactPoint::rational(u), mode, counters) &&
                       decide.outside_connectors(ExactPoint::rational(u), piece_exact[p], counters);
            };
            const Point_2 u = uncovered_beside_line(sweep, t, line_context, accept, counters);
            const std::uint32_t flank_report =
                own != NO_CONNECTOR_REGION ? exact_base + report_connector(exact.regions[static_cast<std::uint32_t>(own)], ew)
                                           : n_full + piece.report_index;
            result.exposed.push_back({CertificatePieceKind::CapsuleFlank, flank_report, n_swept_circles + p,
                                      static_cast<std::uint32_t>(edge), t, u});
        }
        result.witness_s += seconds(t1, Clock::now());
    }

    // The interior anchor: int T is connected, so with no exposed swept edge
    // T is covered iff one interior point is.
    const Point_2 anchor = target.interior_anchor(counters);
    if (!swept(anchor, set, rects, lines, exact, r, counters)) {
        result.exposed.push_back({CertificatePieceKind::TargetInterior, member, static_cast<std::uint32_t>(result.piece_count), 0,
                                  FT(0), anchor});
    }
    result.elapsed_s = seconds(start, Clock::now());
    return result;
}

}  // namespace compas_cgal::coverage_certificate
