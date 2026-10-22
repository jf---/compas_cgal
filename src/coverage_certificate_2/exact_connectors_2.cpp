#include "coverage_certificate_2/exact_connectors_2.h"

#include "coverage_certificate_2/errors_2.h"
#include "coverage_certificate_2/exact_geometry_2.h"

#include <CGAL/box_intersection_d.h>

#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <utility>

namespace compas_cgal::coverage_certificate {

namespace {

/// VertexCircle::disk of a cap that is not a target vertex disk.
constexpr std::uint32_t NO_TARGET_DISK = std::numeric_limits<std::uint32_t>::max();

/// Identical representation (same parts, same root): equal values, decided
/// without arithmetic (the converse need not hold; the exact test follows).
bool identical(const OneRoot& a, const OneRoot& b)
{
    if (a.is_extended() != b.is_extended() || a.a0() != b.a0()) return false;
    return !a.is_extended() || (a.a1() == b.a1() && a.root() == b.root());
}

bool identical(const OneRootPoint& a, const OneRootPoint& b) { return identical(a.x, b.x) && identical(a.y, b.y); }

bool boxes_overlap(const Interval& ax, const Interval& ay, const Interval& bx, const Interval& by)
{
    return !(ax.sup() < bx.inf() || bx.sup() < ax.inf() || ay.sup() < by.inf() || by.sup() < ay.inf());
}

/// Exact point equality: identical representations, else (LAZY) disjoint
/// certified boxes, else the exact one-root signs.
bool same_point(const OneRootPoint& a, const OneRootPoint& b, const FilterMode mode)
{
    if (identical(a, b)) return true;
    if (mode == FilterMode::Lazy && !boxes_overlap(a.x_i, a.y_i, b.x_i, b.y_i)) return false;
    return sign_difference_one_roots(a.x, b.x) == CGAL::ZERO && sign_difference_one_roots(a.y, b.y) == CGAL::ZERO;
}

ExactPoint exact_of(const OneRootPoint& p) { return {p.x, p.y}; }

/// Exact point equality: identical representations, else (LAZY) disjoint
/// certified boxes (equal points have overlapping enclosures), else exact.
bool same_point_filtered(const OneRoot& ax, const OneRoot& ay, const Interval& ax_i, const Interval& ay_i, const OneRoot& bx,
                         const OneRoot& by, const Interval& bx_i, const Interval& by_i, const FilterMode mode)
{
    if (identical(ax, bx) && identical(ay, by)) return true;
    if (mode == FilterMode::Lazy && !boxes_overlap(ax_i, ay_i, bx_i, by_i)) return false;
    return sign_difference_one_roots(ax, bx) == CGAL::ZERO && sign_difference_one_roots(ay, by) == CGAL::ZERO;
}

/// The travel tangent of dOmega* piece `piece` at its end point p: u for an
/// offset line, sense perp(p - c) for a circle piece (perp(x, y) = (-y, x)).
/// Its coordinates lie in p's field.
std::array<OneRoot, 2> travel_tangent(const ExactTarget& target, const std::uint32_t piece, const OneRootPoint& p)
{
    const auto [is_line, band_index] = target.piece_bands[piece];
    if (is_line) {
        const LineBand& band = target.line_bands[band_index];
        return {OneRoot(band.ux), OneRoot(band.uy)};
    }
    const CircleBand& band = target.circle_bands[band_index];
    const OneRoot ex = p.x - band.c.x(), ey = p.y - band.c.y();
    if (band.sense > 0) return {-ey, ex};
    return {ey, -ex};
}

/// The precondition of the cap-circle theorem that the declared pieces carry
/// (exact_connectors_2.h): they form closed chains -- every piece ends exactly
/// where exactly one piece starts, and every piece start is reached exactly once
/// -- and no junction is reflex. Omega* lies on the left of travel (an offset
/// line's interior normal n is u turned left), so a junction from tangent a to
/// tangent b is convex when a x b > 0, G1 when a x b == 0 and a . b > 0, and
/// reflex (or a reversal) otherwise. Every sign is exact (cross_product towers
/// over the junction point's fields), behind certified interval filters in LAZY.
///
/// Raises:
///     BoundaryChainError: a chain is open or branches, or a junction is reflex.
void verify_boundary_chains(const ExactTarget& target, const FilterMode mode)
{
    constexpr std::uint32_t NONE = std::numeric_limits<std::uint32_t>::max();
    const std::size_t n = target.piece_bands.size();
    std::vector<std::uint32_t> next(n, NONE);
    std::vector<std::uint32_t> reached(n, 0);
    // Candidate (end i, start j) pairs: in LAZY those whose certified boxes
    // overlap (CGAL::box_intersection_d, closed boxes: a superset of the equal
    // pairs); in EXACT_ONLY every pair. Visited in (i, j) order either way.
    std::vector<std::vector<std::uint32_t>> candidates(n);
    if (mode == FilterMode::Lazy) {
        using PointBox = CGAL::Box_intersection_d::Box_with_info_d<double, 2, std::uint32_t, CGAL::Box_intersection_d::ID_EXPLICIT>;
        const auto box_of = [](const OneRootPoint& p, const std::uint32_t index) {
            double lo[2] = {p.x_i.inf(), p.y_i.inf()};
            double hi[2] = {p.x_i.sup(), p.y_i.sup()};
            return PointBox(lo, hi, index);
        };
        std::vector<PointBox> ends, starts;
        for (std::uint32_t i = 0; i < n; ++i) {
            ends.push_back(box_of(target.piece_targets[i], i));
            starts.push_back(box_of(target.piece_sources[i], i));
        }
        CGAL::box_intersection_d(ends.begin(), ends.end(), starts.begin(), starts.end(),
                                 [&](const PointBox& end, const PointBox& start) { candidates[end.info()].push_back(start.info()); });
        for (std::vector<std::uint32_t>& list : candidates) std::sort(list.begin(), list.end());
    } else {
        for (std::vector<std::uint32_t>& list : candidates) {
            list.resize(n);
            for (std::uint32_t j = 0; j < n; ++j) list[j] = j;
        }
    }
    for (std::uint32_t i = 0; i < n; ++i) {
        const OneRootPoint& end = target.piece_targets[i];
        for (const std::uint32_t j : candidates[i]) {
            const OneRootPoint& start = target.piece_sources[j];
            if (!same_point_filtered(end.x, end.y, end.x_i, end.y_i, start.x, start.y, start.x_i, start.y_i, mode)) continue;
            if (next[i] != NONE) {
                throw BoundaryChainError("dOmega* piece " + std::to_string(i) + " ends where two pieces start");
            }
            next[i] = j;
            ++reached[j];
        }
        if (next[i] == NONE) throw BoundaryChainError("dOmega* piece " + std::to_string(i) + " ends where no piece starts");
    }
    for (std::uint32_t j = 0; j < n; ++j) {
        if (reached[j] != 1) {
            throw BoundaryChainError("dOmega* piece " + std::to_string(j) + " starts where " + std::to_string(reached[j]) +
                                     " pieces end (closed chains need exactly one)");
        }
    }
    const auto decide = [&](const Interval& enclosure, auto&& exact) {
        if (mode == FilterMode::Lazy) {
            if (enclosure.inf() > 0) return CGAL::POSITIVE;
            if (enclosure.sup() < 0) return CGAL::NEGATIVE;
        }
        return static_cast<CGAL::Sign>(exact());
    };
    for (std::uint32_t i = 0; i < n; ++i) {
        const std::array<OneRoot, 2> a = travel_tangent(target, i, target.piece_targets[i]);
        const std::array<OneRoot, 2> b = travel_tangent(target, next[i], target.piece_sources[next[i]]);
        const Interval ax = interval_of(a[0]), ay = interval_of(a[1]), bx = interval_of(b[0]), by = interval_of(b[1]);
        const CGAL::Sign turn = decide(ax * by - ay * bx, [&] {
            return sign_tower(tower_sub(cross_product(a[0], b[1]), cross_product(a[1], b[0])));
        });
        if (turn == CGAL::POSITIVE) continue;
        const CGAL::Sign along = decide(ax * bx + ay * by, [&] {
            return sign_tower(tower_add(cross_product(a[0], b[0]), cross_product(a[1], b[1])));
        });
        if (turn == CGAL::ZERO && along == CGAL::POSITIVE) continue;
        throw BoundaryChainError("the junction from dOmega* piece " + std::to_string(i) + " to piece " +
                                 std::to_string(next[i]) + (turn == CGAL::NEGATIVE ? " is reflex" : " reverses the travel"));
    }
}

CGAL::Bbox_2 disk_box(const Interval& cx, const Interval& cy, const Interval& radius)
{
    return CGAL::Bbox_2((cx - radius).inf(), (cy - radius).inf(), (cx + radius).sup(), (cy + radius).sup());
}

/// Certified outward box of the annular sector {inner <= |x - c| <= outer, x in
/// the closed wedge of the arc about c from w0 to w1 in `sense` with `span`}
/// (R: the arc's own radius, |w - c| = R). The sector is bounded by two radial
/// segments, whose ends are the four corners c + (w - c) rho / R, and two arcs;
/// an arc's coordinate extremes are at its ends (corners) or at an axis point
/// c +/- rho e that lies in the wedge, and the inner arc's axis points are
/// dominated by the outer arc's. So the box of the corners and of every outer
/// axis point not certified outside the closed wedge contains the sector.
CGAL::Bbox_2 sector_box(const Interval& cx, const Interval& cy, const Interval& big_r, const Interval& outer,
                        const Interval& inner, const OneRootPoint& w0, const OneRootPoint& w1, const int sense,
                        const ArcSpan span)
{
    double x0 = std::numeric_limits<double>::infinity(), y0 = x0, x1 = -x0, y1 = -x0;
    const auto add = [&](const Interval& x, const Interval& y) {
        x0 = std::min(x0, x.inf());
        y0 = std::min(y0, y.inf());
        x1 = std::max(x1, x.sup());
        y1 = std::max(y1, y.sup());
    };
    for (const OneRootPoint* w : {&w0, &w1}) {
        for (const Interval& rho : {outer, inner}) {
            const Interval scale = rho / big_r;
            add(cx + (w->x_i - cx) * scale, cy + (w->y_i - cy) * scale);
        }
    }
    const std::array<std::array<int, 2>, 4> axes{{{1, 0}, {-1, 0}, {0, 1}, {0, -1}}};
    for (const auto& e : axes) {
        const Interval ex = cx + Interval(e[0]), ey = cy + Interval(e[1]);
        // o0 = sense (w0 - c) x (e - c), o1 = -sense (w1 - c) x (e - c); "certainly negative".
        const Interval c0 = cross_about_i(cx, cy, w0.x_i, w0.y_i, ex, ey);
        const Interval c1 = cross_about_i(cx, cy, w1.x_i, w1.y_i, ex, ey);
        const bool o0_neg = sense > 0 ? c0.sup() < 0 : c0.inf() > 0;
        const bool o1_neg = sense > 0 ? c1.inf() > 0 : c1.sup() < 0;
        bool outside = false;
        switch (span) {
            case ArcSpan::Minor: outside = o0_neg || o1_neg; break;
            case ArcSpan::Half: outside = o0_neg; break;
            case ArcSpan::Major: outside = o0_neg && o1_neg; break;
        }
        if (!outside) add(cx + outer * Interval(e[0]), cy + outer * Interval(e[1]));
    }
    return CGAL::Bbox_2(x0, y0, x1, y1);
}

/// span_of behind a certified filter on its orientation sign (LAZY); the exact
/// rule decides a straddling enclosure (and the half-turn case).
ArcSpan filtered_span(const Point_2& c, const OneRootPoint& a, const OneRootPoint& b, const int sense, const FilterMode mode)
{
    if (mode == FilterMode::Lazy) {
        const Interval cross = cross_about_i(interval_of(c.x()), interval_of(c.y()), a.x_i, a.y_i, b.x_i, b.y_i);
        if (cross.inf() > 0) return sense > 0 ? ArcSpan::Minor : ArcSpan::Major;
        if (cross.sup() < 0) return sense > 0 ? ArcSpan::Major : ArcSpan::Minor;
    }
    return span_of(c, a.x, a.y, b.x, b.y, sense);
}

/// Position along an offset-line piece: u.p.
OneRoot along(const LineBand& band, const OneRootPoint& p) { return band.ux * p.x + band.uy * p.y; }

/// Strict order of two points on one circle piece, travelling from `source`
/// about c in `sense` (half-turn class, then orientation; exact over mixed
/// radicands, exact_geometry_2.h). Every sign is filtered (LAZY): a certified
/// interval enclosure that excludes 0 decides it, otherwise the exact tower rule
/// does; EXACT_ONLY always runs the exact rule.
class ArcOrder {
public:
    ArcOrder(const Point_2& c, const int sense, const OneRootPoint& source, const FilterMode mode)
        : c_(c), sense_(sense), source_(source), mode_(mode), cx_(interval_of(c.x())), cy_(interval_of(c.y()))
    {
    }

    /// 0 at the source, 1 the first open half turn, 2 the antipode, 3 the second half turn.
    [[nodiscard]] int half(const OneRootPoint& p) const
    {
        if (identical(p, source_)) return 0;
        const CGAL::Sign s = times(sense_, decide(cross_about_i(cx_, cy_, source_.x_i, source_.y_i, p.x_i, p.y_i), [&] {
                                       return sign_tower(cross_about(c_, source_.x, source_.y, exact_of(p)));
                                   }));
        if (s == CGAL::POSITIVE) return 1;
        if (s == CGAL::NEGATIVE) return 3;
        const Interval dot_i = (source_.x_i - cx_) * (p.x_i - cx_) + (source_.y_i - cy_) * (p.y_i - cy_);
        const CGAL::Sign dot = decide(dot_i, [&] {
            return sign_tower(
                tower_add(cross_product(source_.x - c_.x(), p.x - c_.x()), cross_product(source_.y - c_.y(), p.y - c_.y())));
        });
        return dot == CGAL::POSITIVE ? 0 : 2;
    }

    [[nodiscard]] bool less(const OneRootPoint& p, const OneRootPoint& q) const
    {
        if (identical(p, q)) return false;
        const int hp = half(p);
        const int hq = half(q);
        if (hp != hq) return hp < hq;
        if (hp == 0 || hp == 2) return false;
        return times(sense_, decide(cross_about_i(cx_, cy_, p.x_i, p.y_i, q.x_i, q.y_i),
                                    [&] { return sign_tower(cross_about(c_, p.x, p.y, exact_of(q))); })) == CGAL::POSITIVE;
    }

private:
    template <class Exact>
    [[nodiscard]] CGAL::Sign decide(const Interval& enclosure, const Exact& exact) const
    {
        if (mode_ == FilterMode::Lazy) {
            if (enclosure.inf() > 0) return CGAL::POSITIVE;
            if (enclosure.sup() < 0) return CGAL::NEGATIVE;
        }
        return exact();
    }

    Point_2 c_;
    int sense_;
    OneRootPoint source_;
    FilterMode mode_;
    Interval cx_, cy_;
};

/// One covered stretch of a target piece, in the piece's own order. `to_end`
/// marks a stretch that reaches the piece's end (for a closed circle piece, the
/// full turn back to its source).
struct Stretch {
    OneRootPoint lo, hi;
    bool to_end;
};

/// Verify that the stretches cover the piece from `source` to its end.
template <class Less>
bool covers_piece(std::vector<Stretch> stretches, const OneRootPoint& source, const OneRootPoint& target, const bool closed,
                  const Less& less)
{
    std::sort(stretches.begin(), stretches.end(), [&](const Stretch& a, const Stretch& b) { return less(a.lo, b.lo); });
    OneRootPoint reached = source;
    bool at_end = false;
    for (const Stretch& s : stretches) {
        if (at_end) break;
        if (less(reached, s.lo)) return false;  // a gap before s
        if (s.to_end) {
            at_end = true;
        } else if (less(reached, s.hi)) {
            reached = s.hi;
        }
    }
    if (at_end) return true;
    if (closed) return false;
    return !less(reached, target);
}

}  // namespace

VertexCircle make_vertex_circle(const OneRoot& x, const OneRoot& y, const FT& r, const std::uint32_t disk)
{
    VertexCircle V;
    V.x = x;
    V.y = y;
    V.delta = x.is_extended() ? x.root() : y.root();
    V.r_sq = r * r;
    V.x_i = interval_of(x);
    V.y_i = interval_of(y);
    V.r_i = interval_of(r);
    V.r_sq_i = interval_of(V.r_sq);
    V.bbox = disk_box(V.x_i, V.y_i, V.r_i);
    V.disk = disk;
    return V;
}

ConnectorSet ConnectorSet::build(const std::vector<ExactConnectorSpec>& connectors, const ExactTarget& target,
                                 CircleRegionSet& set, LineTable& lines, std::vector<VertexCircle>& vertices,
                                 const FilterMode mode)
{
    ConnectorSet out;
    const FT& r = target.tool_radius;
    const Interval r_i = interval_of(r);
    std::map<std::tuple<FT, FT, FT>, std::uint32_t> flank_ids;      // (nx, ny, level') -> line
    std::map<std::tuple<FT, FT, FT>, std::uint32_t> support_ids;    // (cx, cy, R) -> region
    std::map<std::pair<FT, FT>, std::uint32_t> rational_caps;       // centre -> region
    const auto report_of = [&](const std::uint32_t region, const std::uint32_t connector) {
        if (out.support_reports.size() <= region) out.support_reports.resize(region + 1, NO_TARGET_DISK);
        out.support_reports[region] = std::min(out.support_reports[region], connector);
    };

    const auto flank_line = [&](const FT& nx, const FT& ny, const FT& level) {
        const auto key = std::make_tuple(nx, ny, level);
        const auto found = flank_ids.find(key);
        if (found != flank_ids.end()) return found->second;
        const std::uint32_t id = lines.add_rational_line(nx, ny, Point_2(level * nx, level * ny));
        flank_ids.emplace(key, id);
        out.flank_lines.push_back(id);
        return id;
    };
    // LAZY cap lookup: vertex circles by the lower end of their certified x
    // enclosure, kept in step with `vertices` (which only grows).
    constexpr std::uint32_t NO_VERTEX = std::numeric_limits<std::uint32_t>::max();
    std::multimap<double, std::uint32_t> by_x_lower;
    double widest_x = 0.0;
    std::size_t indexed = 0;
    const auto index_vertices = [&] {
        for (; indexed < vertices.size(); ++indexed) {
            const Interval& x = vertices[indexed].x_i;
            by_x_lower.emplace(x.inf(), static_cast<std::uint32_t>(indexed));
            // Rounded up: the search window only grows.
            widest_x = std::max(widest_x, (Interval(x.sup()) - Interval(x.inf())).sup());
        }
    };
    const auto cap = [&](const OneRootPoint& w, const std::uint32_t connector, std::int32_t& vertex, std::int32_t& region) {
        vertex = -1;
        region = -1;
        if (w.is_rational()) {
            const auto key = std::make_pair(w.x.a0(), w.y.a0());
            auto found = rational_caps.find(key);
            if (found == rational_caps.end()) {
                found = rational_caps.emplace(key, set.add_support_region(Point_2(key.first, key.second), r, FT(0), connector)).first;
            }
            region = static_cast<std::int32_t>(found->second);
            report_of(found->second, connector);
            return;
        }
        // Shared by VALUE with an identical vertex circle already there (a target
        // vertex disk or an earlier cap): exact equality, the lowest index first.
        // LAZY visits only the circles whose certified x enclosure can meet w's
        // (an index ordered by lower end, searched back by the widest enclosure);
        // EXACT_ONLY visits every circle.
        std::uint32_t found = NO_VERTEX;
        const auto visit = [&](const std::uint32_t v) {
            const VertexCircle& V = vertices[v];
            if (v < found && same_point_filtered(w.x, w.y, w.x_i, w.y_i, V.x, V.y, V.x_i, V.y_i, mode)) found = v;
        };
        if (mode == FilterMode::Lazy) {
            index_vertices();
            const auto last = by_x_lower.upper_bound(w.x_i.sup());
            // Rounded down: a circle whose enclosure meets w's has lower end >= w.inf - widest.
            const double first = (Interval(w.x_i.inf()) - Interval(widest_x)).inf();
            for (auto it = by_x_lower.lower_bound(first); it != last; ++it) visit(it->second);
        } else {
            for (std::uint32_t v = 0; v < vertices.size(); ++v) visit(v);
        }
        if (found != NO_VERTEX) {
            vertex = static_cast<std::int32_t>(found);
            return;
        }
        vertices.push_back(make_vertex_circle(w.x, w.y, r, NO_TARGET_DISK));
        vertex = static_cast<std::int32_t>(vertices.size() - 1);
    };

    // Covered stretches per target piece, for the dOmega* cover check.
    std::vector<std::vector<Stretch>> stretches(target.piece_bands.size());
    // `sub`: a connector's own sub-piece (validated, recorded for the cover check
    // and for reporting); otherwise a whole dOmega* piece, the region actually
    // swept, whose supports and end caps are registered.
    const auto make = [&](const ExactConnectorPieceSpec& spec, const std::uint32_t k, const bool sub) {
            if (spec.piece >= target.piece_bands.size()) {
                throw InvalidCertificateInputError("an exact connector names dOmega* piece " + std::to_string(spec.piece) +
                                                   " of " + std::to_string(target.piece_bands.size()));
            }
            const auto [is_line, band_index] = target.piece_bands[spec.piece];
            ConnectorRegion R;
            R.connector = k;
            R.piece = spec.piece;
            R.is_line = is_line;
            R.w0 = spec.start;
            R.w1 = spec.end;
            R.full = spec.full;
            R.degenerate = !spec.full && same_point(spec.start, spec.end, mode);
            R.sense = 1;
            R.span = ArcSpan::Major;
            R.support_region = 0;
            R.flank = {0, 0};
            if (is_line) {
                const LineBand& band = target.line_bands[band_index];
                for (const OneRootPoint* w : {&spec.start, &spec.end}) {
                    if (sign_one_root(band.nx * w->x + band.ny * w->y - band.level) != CGAL::ZERO) {
                        throw InvalidCertificateInputError("an exact connector end is not on its offset line");
                    }
                }
                if (spec.full) throw InvalidCertificateInputError("a line sub-piece cannot be a full turn");
                R.nx = band.nx;
                R.ny = band.ny;
                R.level = band.level;
                R.ux = band.ux;
                R.uy = band.uy;
                R.s0 = along(band, spec.start);
                R.s1 = along(band, spec.end);
                R.ux_i = interval_of(R.ux);
                R.uy_i = interval_of(R.uy);
                R.nx_i = interval_of(R.nx);
                R.ny_i = interval_of(R.ny);
                R.level_i = interval_of(R.level);
                R.s0_i = interval_of(R.s0);
                R.s1_i = interval_of(R.s1);
                if (sign_difference_one_roots(R.s1, R.s0) == CGAL::NEGATIVE ||
                    sign_difference_one_roots(R.s0, band.s0) == CGAL::NEGATIVE ||
                    sign_difference_one_roots(band.s1, R.s1) == CGAL::NEGATIVE) {
                    throw InvalidCertificateInputError("an exact line sub-piece runs backwards or leaves its piece");
                }
                if (!sub) R.flank = {flank_line(band.nx, band.ny, band.level + r), flank_line(band.nx, band.ny, band.level - r)};
                // On a flank n.x = level' (|n| = 1) the chart is x(t) = level' n + t u, so t = u.x.
                for (std::size_t f = 0; f < 2; ++f) {
                    R.t_lo[f] = FT(interval_of(R.s0).inf());
                    R.t_hi[f] = FT(interval_of(R.s1).sup());
                }
                const Interval nx = interval_of(band.nx), ny = interval_of(band.ny), ux = interval_of(band.ux),
                               uy = interval_of(band.uy);
                const Interval l_hi = interval_of(band.level) + r_i, l_lo = interval_of(band.level) - r_i;
                double x0 = std::numeric_limits<double>::infinity(), y0 = x0, x1 = -x0, y1 = -x0;
                for (const Interval& s : {interval_of(R.s0), interval_of(R.s1)}) {
                    for (const Interval& l : {l_lo, l_hi}) {
                        const Interval x = s * ux + l * nx, y = s * uy + l * ny;
                        x0 = std::min(x0, x.inf());
                        y0 = std::min(y0, y.inf());
                        x1 = std::max(x1, x.sup());
                        y1 = std::max(y1, y.sup());
                    }
                }
                R.bbox = CGAL::Bbox_2(x0, y0, x1, y1);
                if (sub) stretches[spec.piece].push_back({spec.start, spec.end, false});
            } else {
                const CircleBand& band = target.circle_bands[band_index];
                const FT radius_sq = band.radius * band.radius;
                for (const OneRootPoint* w : {&spec.start, &spec.end}) {
                    const OneRoot ex = w->x - band.c.x(), ey = w->y - band.c.y();
                    if (sign_one_root(ex * ex + ey * ey - radius_sq) != CGAL::ZERO) {
                        throw InvalidCertificateInputError("an exact connector end is not on its offset circle");
                    }
                }
                R.c = band.c;
                R.radius = band.radius;
                R.sense = band.sense;
                const bool closed_piece = same_point(band.w0, band.w1, mode);
                const bool degenerate = same_point(spec.start, spec.end, mode);
                if (spec.full && !(closed_piece && degenerate)) {
                    throw InvalidCertificateInputError("a full-turn sub-piece must start and end at its closed piece's source");
                }
                R.span = (spec.full || degenerate) ? ArcSpan::Major
                                                   : filtered_span(band.c, spec.start, spec.end, band.sense, mode);
                const ArcOrder order(band.c, band.sense, band.w0, mode);
                // On an open piece both ends lie in [w0, w1] (every point of a
                // closed piece does): checked like a line sub-piece's range.
                if (!closed_piece) {
                    for (const OneRootPoint* w : {&spec.start, &spec.end}) {
                        if (order.less(band.w1, *w)) throw InvalidCertificateInputError("an exact arc sub-piece leaves its piece");
                    }
                }
                if (!spec.full && !degenerate && !closed_piece && order.less(spec.end, spec.start)) {
                    throw InvalidCertificateInputError("an exact arc sub-piece runs backwards");
                }
                const FT outer = band.radius + r;
                const FT inner = CGAL::compare(band.radius, r) == CGAL::LARGER ? FT(band.radius - r) : FT(0);
                R.annulus = CGAL::compare(band.radius, r) == CGAL::LARGER;
                R.outer_sq = outer * outer;
                R.inner_sq = inner * inner;
                R.cx_i = interval_of(band.c.x());
                R.cy_i = interval_of(band.c.y());
                R.outer_sq_i = interval_of(R.outer_sq);
                R.inner_sq_i = interval_of(R.inner_sq);
                if (!sub) {
                    const auto key = std::make_tuple(band.c.x(), band.c.y(), band.radius);
                    auto found = support_ids.find(key);
                    if (found == support_ids.end()) found = support_ids.emplace(key, set.add_support_region(band.c, outer, inner, k)).first;
                    R.support_region = found->second;
                    report_of(found->second, k);
                }
                R.bbox = (spec.full || degenerate)
                             ? disk_box(interval_of(band.c.x()), interval_of(band.c.y()), interval_of(outer))
                             : sector_box(R.cx_i, R.cy_i, interval_of(band.radius), interval_of(outer), interval_of(inner),
                                          spec.start, spec.end, band.sense, R.span);
                if (!sub) {
                    // A whole piece records no stretch.
                } else if (spec.full) {
                    stretches[spec.piece].push_back({band.w0, band.w0, true});
                } else if (degenerate) {
                    stretches[spec.piece].push_back({spec.start, spec.end, false});
                } else if (closed_piece && order.less(spec.end, spec.start)) {
                    stretches[spec.piece].push_back({spec.start, spec.start, true});
                    stretches[spec.piece].push_back({band.w0, spec.end, false});
                } else {
                    const bool to_end = closed_piece && same_point(spec.end, band.w0, mode);
                    stretches[spec.piece].push_back({spec.start, spec.end, to_end});
                }
            }
            R.cap_vertex = {-1, -1};
            R.cap_region = {-1, -1};
            // A full turn's end disk lies inside its own annulus: no cap.
            if (!sub && !spec.full) {
                cap(spec.start, k, R.cap_vertex[0], R.cap_region[0]);
                cap(spec.end, k, R.cap_vertex[1], R.cap_region[1]);
            }
            for (const OneRootPoint* w : {&spec.start, &spec.end}) R.bbox += disk_box(w->x_i, w->y_i, r_i);
            return R;
    };
    for (std::uint32_t k = 0; k < connectors.size(); ++k) {
        for (const ExactConnectorPieceSpec& spec : connectors[k].pieces) out.sub_regions.push_back(make(spec, k, true));
    }

    // The cap-circle theorem needs K = dOmega*: the declared pieces form closed,
    // non-reflex chains, and the sub-pieces cover every one of them.
    if (!out.sub_regions.empty()) verify_boundary_chains(target, mode);
    if (!out.sub_regions.empty()) {
        for (std::uint32_t piece = 0; piece < target.piece_bands.size(); ++piece) {
            const auto [is_line, band_index] = target.piece_bands[piece];
            bool covered = false;
            if (is_line) {
                const LineBand& band = target.line_bands[band_index];
                // Ordered on the line's own coordinate u.x, from band.s0 to band.s1.
                const Interval ux = interval_of(band.ux), uy = interval_of(band.uy);
                const auto less = [&](const OneRootPoint& a, const OneRootPoint& b) {
                    if (identical(a, b)) return false;
                    if (mode == FilterMode::Lazy) {
                        const Interval gap = (ux * a.x_i + uy * a.y_i) - (ux * b.x_i + uy * b.y_i);
                        if (gap.sup() < 0) return true;
                        if (gap.inf() > 0) return false;
                    }
                    return sign_difference_one_roots(along(band, a), along(band, b)) == CGAL::NEGATIVE;
                };
                std::vector<Stretch>& list = stretches[piece];
                std::sort(list.begin(), list.end(), [&](const Stretch& a, const Stretch& b) { return less(a.lo, b.lo); });
                OneRoot reached = band.s0;
                for (const Stretch& s : list) {
                    if (sign_difference_one_roots(along(band, s.lo), reached) == CGAL::POSITIVE) break;
                    const OneRoot hi = along(band, s.hi);
                    if (sign_difference_one_roots(hi, reached) == CGAL::POSITIVE) reached = hi;
                }
                covered = sign_difference_one_roots(reached, band.s1) != CGAL::NEGATIVE;
            } else {
                const CircleBand& band = target.circle_bands[band_index];
                const bool closed_piece = same_point(band.w0, band.w1, mode);
                const ArcOrder order(band.c, band.sense, band.w0, mode);
                covered = covers_piece(stretches[piece], band.w0, band.w1, closed_piece,
                                       [&](const OneRootPoint& a, const OneRootPoint& b) { return order.less(a, b); });
            }
            if (!covered) {
                throw ConnectorBoundaryCoverageError(
                    "the exact connectors leave part of dOmega* piece " + std::to_string(piece) +
                    " unswept; the one-root cap circles would have to be checked as pieces, which is not modelled");
            }
        }
    }

    // K = dOmega*, so the swept set of the connectors is exactly dOmega* (+) D_r:
    // the union of piece (+) D_r over the WHOLE pieces, whatever the sub-piece
    // split (each sub-piece (+) D_r lies in its piece's, and the pieces' are
    // covered by the sub-pieces'). One region per piece: its end caps are the
    // dOmega* vertices only, never the sites inside a piece.
    if (!out.sub_regions.empty()) {
        std::vector<std::uint32_t> lowest(target.piece_bands.size(), NO_TARGET_DISK);
        for (const ConnectorRegion& R : out.sub_regions) lowest[R.piece] = std::min(lowest[R.piece], R.connector);
        for (std::uint32_t piece = 0; piece < target.piece_bands.size(); ++piece) {
            const auto [is_line, band_index] = target.piece_bands[piece];
            const OneRootPoint& source = is_line ? target.piece_sources[piece] : target.circle_bands[band_index].w0;
            const OneRootPoint& end = is_line ? target.piece_targets[piece] : target.circle_bands[band_index].w1;
            const bool full = !is_line && same_point(source, end, mode);
            out.regions.push_back(make({piece, full, source, end}, lowest[piece], false));
        }
    }
    std::sort(out.flank_lines.begin(), out.flank_lines.end());

    out.cap_circle.assign(set.circles.size(), 0U);
    for (const auto& [centre, region] : rational_caps) {
        (void)centre;
        for (std::uint32_t c = set.first_circle[region]; c < set.first_circle[region + 1]; ++c) out.cap_circle[c] = 1U;
    }
    return out;
}

Side connector_membership(const ConnectorRegion& region, const ExactPoint& p, const FT& r, const FilterMode mode,
                          DecisionCounters& counters)
{
    const Interval px = interval_of(p.x), py = interval_of(p.y);
    // The region lies in its certified box: a point whose enclosure misses it is outside.
    if (mode == FilterMode::Lazy && (px.sup() < region.bbox.xmin() || px.inf() > region.bbox.xmax() ||
                                     py.sup() < region.bbox.ymin() || py.inf() > region.bbox.ymax())) {
        return Side::Exterior;
    }
    // Interval filter, then the exact sign (EXACT_ONLY: always exact).
    const auto decide = [&](const Interval& enclosure, auto&& exact) {
        if (mode == FilterMode::Lazy) {
            if (enclosure.inf() > 0) return CGAL::POSITIVE;
            if (enclosure.sup() < 0) return CGAL::NEGATIVE;
        }
        ++counters.exact_connector_membership;
        return static_cast<CGAL::Sign>(exact());
    };
    const FT r2 = r * r;
    const Interval r_i = interval_of(r), r2_i = interval_of(r2);
    bool boundary = false;
    for (const OneRootPoint* w : {&region.w0, &region.w1}) {
        const Interval dx = px - w->x_i, dy = py - w->y_i;
        const CGAL::Sign s = decide(dx * dx + dy * dy - r2_i, [&] { return sign_tower(squared_distance_minus(p, w->x, w->y, r2)); });
        if (s == CGAL::NEGATIVE) return Side::Interior;
        if (s == CGAL::ZERO) boundary = true;
    }
    const Side outside = boundary ? Side::Boundary : Side::Exterior;
    if (region.degenerate) return outside;
    if (region.is_line) {
        const Interval u_i = region.ux_i * px + region.uy_i * py;
        const auto u = [&] { return OneRoot(region.ux * p.x + region.uy * p.y); };
        const CGAL::Sign a0 = decide(u_i - region.s0_i, [&] { return sign_difference_one_roots(u(), region.s0); });
        if (a0 == CGAL::NEGATIVE) return outside;
        const CGAL::Sign a1 = decide(u_i - region.s1_i, [&] { return sign_difference_one_roots(u(), region.s1); });
        if (a1 == CGAL::POSITIVE) return outside;
        const Interval level_i = region.nx_i * px + region.ny_i * py - region.level_i;
        const auto level = [&] { return OneRoot(region.nx * p.x + region.ny * p.y - region.level); };
        const CGAL::Sign hi = decide(level_i - r_i, [&] { return sign_one_root(level() - r); });
        if (hi == CGAL::POSITIVE) return outside;
        const CGAL::Sign lo = decide(level_i + r_i, [&] { return sign_one_root(level() + r); });
        if (lo == CGAL::NEGATIVE) return outside;
        if (hi == CGAL::ZERO || lo == CGAL::ZERO || a0 == CGAL::ZERO || a1 == CGAL::ZERO) return Side::Boundary;
        return Side::Interior;
    }
    const Interval& cx = region.cx_i;
    const Interval& cy = region.cy_i;
    const Interval ex_i = px - cx, ey_i = py - cy, d2_i = ex_i * ex_i + ey_i * ey_i;
    const auto d2 = [&] {
        const OneRoot ex = p.x - region.c.x(), ey = p.y - region.c.y();
        return OneRoot(ex * ex + ey * ey);
    };
    const CGAL::Sign so = decide(d2_i - region.outer_sq_i, [&] { return sign_one_root(d2() - region.outer_sq); });
    if (so == CGAL::POSITIVE) return outside;
    CGAL::Sign si = CGAL::POSITIVE;
    if (region.annulus) {
        si = decide(d2_i - region.inner_sq_i, [&] { return sign_one_root(d2() - region.inner_sq); });
        if (si == CGAL::NEGATIVE) return outside;
    }
    if (region.full) return (so == CGAL::ZERO || si == CGAL::ZERO) ? Side::Boundary : Side::Interior;
    const auto side_of = [&](const OneRootPoint& w) {
        return decide(cross_about_i(cx, cy, w.x_i, w.y_i, px, py), [&] { return sign_tower(cross_about(region.c, w.x, w.y, p)); });
    };
    const CGAL::Sign o0 = times(region.sense, side_of(region.w0));
    const CGAL::Sign o1 = CGAL::opposite(times(region.sense, side_of(region.w1)));
    if (!in_wedge(region.span, o0, o1)) return outside;
    if (so == CGAL::ZERO || si == CGAL::ZERO || !in_open_wedge(region.span, o0, o1)) return Side::Boundary;
    return Side::Interior;
}

bool connector_in_range(const ConnectorRegion& region, const ExactPoint& p, const FilterMode mode)
{
    if (region.degenerate) return false;
    const Interval px = interval_of(p.x), py = interval_of(p.y);
    const auto decide = [&](const Interval& enclosure, auto&& exact) {
        if (mode == FilterMode::Lazy) {
            if (enclosure.inf() > 0) return CGAL::POSITIVE;
            if (enclosure.sup() < 0) return CGAL::NEGATIVE;
        }
        return static_cast<CGAL::Sign>(exact());
    };
    if (region.is_line) {
        const Interval u_i = region.ux_i * px + region.uy_i * py;
        const auto u = [&] { return OneRoot(region.ux * p.x + region.uy * p.y); };
        return decide(u_i - region.s0_i, [&] { return sign_difference_one_roots(u(), region.s0); }) == CGAL::POSITIVE &&
               decide(u_i - region.s1_i, [&] { return sign_difference_one_roots(u(), region.s1); }) == CGAL::NEGATIVE;
    }
    if (region.full) return true;
    const auto side_of = [&](const OneRootPoint& w) {
        return decide(cross_about_i(region.cx_i, region.cy_i, w.x_i, w.y_i, px, py),
                      [&] { return sign_tower(cross_about(region.c, w.x, w.y, p)); });
    };
    const CGAL::Sign o0 = times(region.sense, side_of(region.w0));
    const CGAL::Sign o1 = CGAL::opposite(times(region.sense, side_of(region.w1)));
    return in_open_wedge(region.span, o0, o1);
}

}  // namespace compas_cgal::coverage_certificate
