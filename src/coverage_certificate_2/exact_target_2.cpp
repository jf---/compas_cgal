#include "coverage_certificate_2/exact_target_2.h"

#include "coverage_certificate_2/errors_2.h"
#include "coverage_certificate_2/line_sweep_2.h"

#include <algorithm>
#include <string>
#include <utility>

namespace compas_cgal::coverage_certificate {

namespace {

/// Interval filter, then the exact sign (counted as a target-membership
/// fallback).
template <class ExactFn>
CGAL::Sign filtered(const Interval& enclosure, const FilterMode mode, DecisionCounters& counters, ExactFn&& exact)
{
    if (mode == FilterMode::Lazy) {
        if (enclosure.inf() > 0) return CGAL::POSITIVE;
        if (enclosure.sup() < 0) return CGAL::NEGATIVE;
    }
    ++counters.exact_target_membership;
    return exact();
}

/// sign(a - b) for one-roots over possibly different radicands.
CGAL::Sign sign_difference_one_roots(const OneRoot& a, const OneRoot& b)
{
    const FT ba = base_of(a);
    const FT bb = base_of(b);
    if (CGAL::sign(ba) == CGAL::ZERO || CGAL::sign(bb) == CGAL::ZERO || ba == bb) return sign_one_root(a - b);
    return sign_rational_two(a.a0() - b.a0(), a.a1(), ba, -b.a1(), bb);
}

/// (w - c) x (p - c) for w over one radicand and p over another, as a tower over
/// p's field (tower_2.h, cross_product).
Tower cross_about(const Point_2& c, const OneRoot& wx, const OneRoot& wy, const ExactPoint& p)
{
    return tower_sub(cross_product(wx - c.x(), p.y - c.y()), cross_product(wy - c.y(), p.x - c.x()));
}

Interval cross_about_i(const Interval& cx, const Interval& cy, const Interval& wx, const Interval& wy, const Interval& px,
                       const Interval& py)
{
    return (wx - cx) * (py - cy) - (wy - cy) * (px - cx);
}

/// |p - w|^2 - r2 for a one-root centre w, as a tower over p's field.
Tower squared_distance_minus(const ExactPoint& p, const OneRoot& wx, const OneRoot& wy, const FT& r2)
{
    const FT delta = wx.is_extended() ? wx.root() : (wy.is_extended() ? wy.root() : FT(0));
    const FT x1 = wx.is_extended() ? wx.a1() : FT(0);
    const FT y1 = wy.is_extended() ? wy.a1() : FT(0);
    const OneRoot ex = p.x - wx.a0();
    const OneRoot ey = p.y - wy.a0();
    // (ex - x1 s)^2 + (ey - y1 s)^2 - r2 with s = sqrt(delta).
    const OneRoot m = ex * ex + ey * ey + OneRoot((x1 * x1 + y1 * y1) * delta - r2);
    const OneRoot n = FT(-2) * (ex * x1 + ey * y1);
    return {m, n, OneRoot(delta)};
}

/// Wedge membership of p for an arc about c from w0 to w1 in `sense` with span
/// `span` (closed; the apex c is inside). o0, o1 are the signs of
/// sense (w0 - c) x (p - c) and sense (p - c) x (w1 - c).
bool in_wedge(const ArcSpan span, const CGAL::Sign o0, const CGAL::Sign o1)
{
    switch (span) {
        case ArcSpan::Minor: return o0 != CGAL::NEGATIVE && o1 != CGAL::NEGATIVE;
        case ArcSpan::Half: return o0 != CGAL::NEGATIVE;
        case ArcSpan::Major: return o0 != CGAL::NEGATIVE || o1 != CGAL::NEGATIVE;
    }
    throw CoverageCertificateNativeError("unknown arc span");
}

CGAL::Sign times(const int sense, const CGAL::Sign s) { return sense > 0 ? s : CGAL::opposite(s); }

ArcSpan span_of(const Point_2& c, const OneRoot& ax, const OneRoot& ay, const OneRoot& bx, const OneRoot& by,
                const int sense)
{
    const ExactPoint b{bx, by};
    const CGAL::Sign s = times(sense, sign_tower(cross_about(c, ax, ay, b)));
    if (s == CGAL::POSITIVE) return ArcSpan::Minor;
    if (s == CGAL::NEGATIVE) return ArcSpan::Major;
    if (sign_difference_one_roots(ax, bx) == CGAL::ZERO && sign_difference_one_roots(ay, by) == CGAL::ZERO) {
        throw InvalidCertificateInputError("an arc of the exact target starts and ends at one point");
    }
    return ArcSpan::Half;
}

CGAL::Bbox_2 box_of(const Interval& x0, const Interval& x1, const Interval& y0, const Interval& y1)
{
    return CGAL::Bbox_2(std::min(x0.inf(), x1.inf()), std::min(y0.inf(), y1.inf()), std::max(x0.sup(), x1.sup()),
                        std::max(y0.sup(), y1.sup()));
}

CGAL::Bbox_2 disk_box(const Interval& cx, const Interval& cy, const Interval& radius)
{
    return box_of(cx - radius, cx + radius, cy - radius, cy + radius);
}

bool overlaps(const CGAL::Bbox_2& box, const Interval& x, const Interval& y)
{
    return !(x.sup() < box.xmin() || x.inf() > box.xmax() || y.sup() < box.ymin() || y.inf() > box.ymax());
}

CGAL::Sign orientation_sign(const Point_2& c, const Point_2& a, const Point_2& x)
{
    const CGAL::Orientation o = CGAL::orientation(c, a, x);
    return o == CGAL::LEFT_TURN ? CGAL::POSITIVE : (o == CGAL::RIGHT_TURN ? CGAL::NEGATIVE : CGAL::ZERO);
}

/// Whether p comes strictly before q travelling from a about c in `sense`
/// (exact; the arc-prefix idiom: half-turn class, then orientation).
bool before_on_arc(const Point_2& c, const int sense, const Point_2& a, const Point_2& p, const Point_2& q)
{
    const auto half = [&](const Point_2& x) {
        const CGAL::Sign s = times(sense, orientation_sign(c, a, x));
        if (s == CGAL::POSITIVE) return 1;
        if (s == CGAL::NEGATIVE) return 3;
        // On the line through a: at a (0) or at the antipode (2).
        return CGAL::sign((x.x() - c.x()) * (a.x() - c.x()) + (x.y() - c.y()) * (a.y() - c.y())) == CGAL::POSITIVE ? 0 : 2;
    };
    const int hp = half(p);
    const int hq = half(q);
    if (hp != hq) return hp < hq;
    return times(sense, orientation_sign(c, p, q)) == CGAL::POSITIVE;
}

DeclaredCircleArc make_arc(const DeclaredPrimitiveSpec& spec)
{
    DeclaredCircleArc arc;
    arc.c = spec.centre;
    arc.a = spec.start;
    arc.b = spec.end;
    arc.rho = spec.radius;
    arc.sense = spec.sense;
    if (CGAL::sign(arc.rho) != CGAL::POSITIVE || (arc.sense != 1 && arc.sense != -1)) {
        throw InvalidCertificateInputError("declared arc needs a positive radius and sense +/-1");
    }
    if (CGAL::squared_distance(arc.a, arc.c) != arc.rho * arc.rho || CGAL::squared_distance(arc.b, arc.c) != arc.rho * arc.rho) {
        throw InvalidCertificateInputError("declared arc endpoints are not on its circle");
    }
    arc.span = span_of(arc.c, OneRoot(arc.a.x()), OneRoot(arc.a.y()), OneRoot(arc.b.x()), OneRoot(arc.b.y()), arc.sense);
    arc.cx_i = interval_of(arc.c.x());
    arc.cy_i = interval_of(arc.c.y());
    arc.rho_i = interval_of(arc.rho);
    arc.bbox = disk_box(arc.cx_i, arc.cy_i, arc.rho_i);
    // y-monotone split at the top and bottom points inside the open arc.
    const Point_2 top(arc.c.x(), arc.c.y() + arc.rho);
    const Point_2 bottom(arc.c.x(), arc.c.y() - arc.rho);
    std::vector<Point_2> cuts;
    for (const Point_2& e : {top, bottom}) {
        if (e == arc.a || e == arc.b) continue;
        if (before_on_arc(arc.c, arc.sense, arc.a, e, arc.b)) cuts.push_back(e);
    }
    if (cuts.size() == 2 && !before_on_arc(arc.c, arc.sense, arc.a, cuts[0], cuts[1])) std::swap(cuts[0], cuts[1]);
    std::vector<Point_2> chain{arc.a};
    chain.insert(chain.end(), cuts.begin(), cuts.end());
    chain.push_back(arc.b);
    for (std::size_t i = 0; i + 1 < chain.size(); ++i) {
        const Point_2& p = chain[i];
        const Point_2& q = chain[i + 1];
        int side = 0;
        for (const Point_2* e : {&p, &q}) {
            const auto cx = CGAL::compare(e->x(), arc.c.x());
            if (cx != CGAL::EQUAL) side = cx == CGAL::LARGER ? 1 : -1;
        }
        if (side == 0) {
            // A half circle from bottom to top (or back): CCW from the bottom runs
            // through the right half.
            const bool upward = CGAL::compare(q.y(), p.y()) == CGAL::LARGER;
            side = (upward == (arc.sense > 0)) ? 1 : -1;
        }
        arc.monotone.push_back({p, q, side});
    }
    return arc;
}

}  // namespace

OneRootPoint OneRootPoint::build(const FT& x0, const FT& x1, const FT& y0, const FT& y1, const FT& delta)
{
    if (CGAL::sign(delta) == CGAL::NEGATIVE) throw InvalidCertificateInputError("one-root point radicand must be >= 0");
    OneRootPoint p{make_reduced_one_root(x0, x1, delta), make_reduced_one_root(y0, y1, delta), FT(0), {}, {}};
    p.delta = p.x.is_extended() ? p.x.root() : (p.y.is_extended() ? p.y.root() : FT(0));
    p.x_i = interval_of(p.x);
    p.y_i = interval_of(p.y);
    return p;
}

ExactTarget ExactTarget::build(const double tool_radius, const std::vector<DeclaredPrimitiveSpec>& primitives,
                               const std::vector<ExactPieceSpec>& pieces, std::vector<VertexDisk> vertex_disks)
{
    ExactTarget target;
    target.tool_radius = inject_binary64(tool_radius);
    if (CGAL::sign(target.tool_radius) != CGAL::POSITIVE) throw InvalidCertificateInputError("tool radius must be positive");
    target.tool_radius_squared = target.tool_radius * target.tool_radius;
    target.tool_radius_squared_i = interval_of(target.tool_radius_squared);
    if (primitives.empty()) throw InvalidCertificateInputError("the exact target needs a declared boundary");
    for (const auto& spec : primitives) {
        if (!spec.is_line) {
            target.arcs.push_back(make_arc(spec));
            continue;
        }
        if (CGAL::sign(spec.length) != CGAL::POSITIVE || spec.ux * spec.ux + spec.uy * spec.uy != FT(1)) {
            throw InvalidCertificateInputError("declared line needs a positive length and a unit direction");
        }
        if (Point_2(spec.start.x() + spec.length * spec.ux, spec.start.y() + spec.length * spec.uy) != spec.end) {
            throw InvalidCertificateInputError("declared line end is not start + length u");
        }
        DeclaredSegment s{spec.start, spec.end, spec.ux, spec.uy, spec.length, interval_of(spec.start.x()),
                          interval_of(spec.start.y()), interval_of(spec.end.x()), interval_of(spec.end.y()), {}};
        s.bbox = box_of(s.ax_i, s.bx_i, s.ay_i, s.by_i);
        target.segments.push_back(std::move(s));
    }
    const Interval r_i = interval_of(target.tool_radius);
    for (const auto& piece : pieces) {
        if (piece.is_line) {
            LineBand band{piece.nx, piece.ny, piece.level, piece.ux, piece.uy, {}, {}, {}};
            if (piece.nx * piece.nx + piece.ny * piece.ny != FT(1) || piece.ux != piece.ny || piece.uy != -piece.nx) {
                throw InvalidCertificateInputError("offset-line piece needs a unit normal and u = (n_y, -n_x)");
            }
            band.s0 = piece.ux * piece.source.x + piece.uy * piece.source.y;
            band.s1 = piece.ux * piece.target.x + piece.uy * piece.target.y;
            if (sign_difference_one_roots(band.s1, band.s0) != CGAL::POSITIVE) {
                throw InvalidCertificateInputError("offset-line piece does not travel forward");
            }
            const Interval s0 = interval_of(band.s0), s1 = interval_of(band.s1);
            const Interval l1 = interval_of(band.level), l0 = l1 - r_i;
            const Interval ux = interval_of(band.ux), uy = interval_of(band.uy), nx = interval_of(band.nx),
                           ny = interval_of(band.ny);
            Interval xs[4], ys[4];
            int k = 0;
            for (const Interval& s : {s0, s1}) {
                for (const Interval& l : {l0, l1}) {
                    xs[k] = s * ux + l * nx;
                    ys[k] = s * uy + l * ny;
                    ++k;
                }
            }
            band.bbox = box_of(xs[0], xs[1], ys[0], ys[1]) + box_of(xs[2], xs[3], ys[2], ys[3]);
            target.line_bands.push_back(std::move(band));
            continue;
        }
        if (piece.sense != 1 && piece.sense != -1) throw InvalidCertificateInputError("circle piece sense must be +/-1");
        const FT& big_r = piece.radius;
        const FT& rho = piece.boundary_radius;
        if (CGAL::sign(big_r) != CGAL::POSITIVE || CGAL::sign(rho) == CGAL::NEGATIVE) {
            throw InvalidCertificateInputError("circle piece needs R > 0 and rho >= 0");
        }
        CircleBand band;
        band.c = piece.centre;
        const bool r_outer = CGAL::compare(big_r, rho) == CGAL::LARGER;
        band.inner_sq = r_outer ? rho * rho : big_r * big_r;
        band.outer_sq = r_outer ? big_r * big_r : rho * rho;
        band.sense = piece.sense;
        band.w0 = piece.source;
        band.w1 = piece.target;
        band.span = (sign_difference_one_roots(band.w0.x, band.w1.x) == CGAL::ZERO &&
                     sign_difference_one_roots(band.w0.y, band.w1.y) == CGAL::ZERO)
                        ? ArcSpan::Major
                        : span_of(band.c, band.w0.x, band.w0.y, band.w1.x, band.w1.y, band.sense);
        band.cx_i = interval_of(band.c.x());
        band.cy_i = interval_of(band.c.y());
        band.inner_sq_i = interval_of(band.inner_sq);
        band.outer_sq_i = interval_of(band.outer_sq);
        band.bbox = disk_box(band.cx_i, band.cy_i, CGAL::sqrt(band.outer_sq_i));
        target.circle_bands.push_back(std::move(band));
    }
    target.vertex_disks = std::move(vertex_disks);
    for (const auto& disk : target.vertex_disks) target.vertex_boxes.push_back(disk_box(disk.x_i, disk.y_i, r_i));
    return target;
}

bool ExactTarget::in_centre_domain(const ExactPoint& p, const FilterMode mode, DecisionCounters& counters) const
{
    const Interval px = interval_of(p.x), py = interval_of(p.y);
    const Interval r2 = tool_radius_squared_i;
    // dist(p, dP) >= r against every primitive (box filter: a primitive whose
    // r-widened box misses p is farther than r).
    const Interval r_i = CGAL::sqrt(r2);
    const auto squared_distance_sign = [&](const Point_2& q, const Interval& qx, const Interval& qy) {
        // sign(|p - q|^2 - r^2)
        const Interval dx = px - qx, dy = py - qy;
        return filtered(dx * dx + dy * dy - r2, mode, counters, [&] {
            const OneRoot ex = p.x - q.x(), ey = p.y - q.y();
            return sign_one_root(ex * ex + ey * ey - tool_radius_squared);
        });
    };
    for (const auto& s : segments) {
        const CGAL::Bbox_2 wide(s.bbox.xmin() - r_i.sup(), s.bbox.ymin() - r_i.sup(), s.bbox.xmax() + r_i.sup(),
                                s.bbox.ymax() + r_i.sup());
        if (!overlaps(wide, px, py)) continue;
        const Interval ux = interval_of(s.ux), uy = interval_of(s.uy);
        const Interval along = ux * (px - s.ax_i) + uy * (py - s.ay_i);
        const CGAL::Sign before = filtered(along, mode, counters, [&] {
            return sign_one_root(s.ux * (p.x - s.a.x()) + s.uy * (p.y - s.a.y()));
        });
        if (before == CGAL::NEGATIVE) {
            if (squared_distance_sign(s.a, s.ax_i, s.ay_i) == CGAL::NEGATIVE) return false;
            continue;
        }
        const CGAL::Sign after = filtered(along - interval_of(s.length), mode, counters, [&] {
            return sign_one_root(s.ux * (p.x - s.a.x()) + s.uy * (p.y - s.a.y()) - s.length);
        });
        if (after == CGAL::POSITIVE) {
            if (squared_distance_sign(s.b, s.bx_i, s.by_i) == CGAL::NEGATIVE) return false;
            continue;
        }
        // Perpendicular distance: (n.(p - a))^2 < r^2 with n = (-u_y, u_x).
        const Interval normal = uy * (s.ax_i - px) + ux * (py - s.ay_i);
        const CGAL::Sign near = filtered(normal * normal - r2, mode, counters, [&] {
            const OneRoot v = s.ux * (p.y - s.a.y()) - s.uy * (p.x - s.a.x());
            return sign_one_root(v * v - tool_radius_squared);
        });
        if (near == CGAL::NEGATIVE) return false;
    }
    for (const auto& arc : arcs) {
        const CGAL::Bbox_2 wide(arc.bbox.xmin() - r_i.sup(), arc.bbox.ymin() - r_i.sup(), arc.bbox.xmax() + r_i.sup(),
                                arc.bbox.ymax() + r_i.sup());
        if (!overlaps(wide, px, py)) continue;
        const Interval ax = interval_of(arc.a.x()), ay = interval_of(arc.a.y()), bx = interval_of(arc.b.x()),
                       by = interval_of(arc.b.y());
        const CGAL::Sign o0 = times(arc.sense, filtered(cross_about_i(arc.cx_i, arc.cy_i, ax, ay, px, py), mode, counters, [&] {
                                        return sign_tower(cross_about(arc.c, OneRoot(arc.a.x()), OneRoot(arc.a.y()), p));
                                    }));
        const CGAL::Sign o1 = CGAL::opposite(times(
            arc.sense, filtered(cross_about_i(arc.cx_i, arc.cy_i, bx, by, px, py), mode, counters, [&] {
                return sign_tower(cross_about(arc.c, OneRoot(arc.b.x()), OneRoot(arc.b.y()), p));
            })));
        if (in_wedge(arc.span, o0, o1)) {
            // (|p - c| - rho)^2 < r^2  <=>  (rho - r)^2 < |p - c|^2 < (rho + r)^2 (lower bound only if rho > r).
            const Interval dx = px - arc.cx_i, dy = py - arc.cy_i, d2 = dx * dx + dy * dy;
            const FT hi = (arc.rho + tool_radius) * (arc.rho + tool_radius);
            const CGAL::Sign below_hi = filtered(d2 - interval_of(hi), mode, counters, [&] {
                const OneRoot ex = p.x - arc.c.x(), ey = p.y - arc.c.y();
                return sign_one_root(ex * ex + ey * ey - hi);
            });
            if (below_hi != CGAL::NEGATIVE) continue;
            if (CGAL::compare(arc.rho, tool_radius) != CGAL::LARGER) return false;
            const FT lo = (arc.rho - tool_radius) * (arc.rho - tool_radius);
            const CGAL::Sign above_lo = filtered(d2 - interval_of(lo), mode, counters, [&] {
                const OneRoot ex = p.x - arc.c.x(), ey = p.y - arc.c.y();
                return sign_one_root(ex * ex + ey * ey - lo);
            });
            if (above_lo == CGAL::POSITIVE) return false;
            continue;
        }
        if (squared_distance_sign(arc.a, ax, ay) == CGAL::NEGATIVE) return false;
        if (squared_distance_sign(arc.b, bx, by) == CGAL::NEGATIVE) return false;
    }
    // p in P: crossing number of the ray to +x over the y-monotone pieces of dP.
    const auto y_compare = [&](const Point_2& q) {
        return filtered(interval_of(q.y()) - py, mode, counters, [&] { return sign_one_root(OneRoot(q.y()) - p.y); });
    };
    const auto is_point = [&](const Point_2& q) {
        return sign_one_root(p.x - q.x()) == CGAL::ZERO && sign_one_root(p.y - q.y()) == CGAL::ZERO;
    };
    bool inside = false;
    for (const auto& s : segments) {
        if (py.sup() < s.bbox.ymin() || py.inf() > s.bbox.ymax() || px.inf() > s.bbox.xmax()) {
            // Entirely above, below, or left of p: no crossing to the right, p not on it.
            continue;
        }
        const CGAL::Sign ca = y_compare(s.a);  // sign(a.y - p.y)
        const CGAL::Sign cb = y_compare(s.b);
        if (ca == cb && ca != CGAL::ZERO) continue;
        const Interval det = (s.bx_i - s.ax_i) * (py - s.ay_i) - (s.by_i - s.ay_i) * (px - s.ax_i);
        const CGAL::Sign o = filtered(det, mode, counters, [&] {
            return sign_one_root((s.b.x() - s.a.x()) * (p.y - s.a.y()) - (s.b.y() - s.a.y()) * (p.x - s.a.x()));
        });
        if (o == CGAL::ZERO) {
            // Collinear: on the segment iff between its ends along u.
            const OneRoot along = s.ux * (p.x - s.a.x()) + s.uy * (p.y - s.a.y());
            if (sign_one_root(along) != CGAL::NEGATIVE && sign_one_root(along - s.length) != CGAL::POSITIVE) return true;
            continue;
        }
        const bool a_above = ca == CGAL::POSITIVE;
        const bool b_above = cb == CGAL::POSITIVE;
        if (a_above == b_above) continue;
        if ((b_above && o == CGAL::POSITIVE) || (a_above && o == CGAL::NEGATIVE)) inside = !inside;
    }
    for (const auto& arc : arcs) {
        if (py.sup() < arc.bbox.ymin() || py.inf() > arc.bbox.ymax() || px.inf() > arc.bbox.xmax()) continue;
        for (const auto& piece : arc.monotone) {
            if (is_point(piece.p) || is_point(piece.q)) return true;
            const CGAL::Sign cp = y_compare(piece.p);
            const CGAL::Sign cq = y_compare(piece.q);
            // Half-open in y: the lower end counts, the upper end does not.
            const bool p_low = CGAL::compare(piece.p.y(), piece.q.y()) == CGAL::SMALLER;
            const CGAL::Sign c_lo = p_low ? cp : cq;
            const CGAL::Sign c_hi = p_low ? cq : cp;
            if (c_lo == CGAL::POSITIVE || c_hi == CGAL::NEGATIVE) continue;  // p.y below lo or above hi
            // The piece's point at height p.y: x = c_x + side sqrt(rho^2 - (p.y - c_y)^2).
            const Interval ey = py - arc.cy_i;
            const Interval disc = arc.rho_i * arc.rho_i - ey * ey;
            const Interval xc = arc.cx_i + Interval(piece.side) * CGAL::sqrt(Interval(std::max(0.0, disc.inf()), std::max(0.0, disc.sup())));
            const CGAL::Sign right = filtered(xc - px, mode, counters, [&] {
                const OneRoot e = p.y - arc.c.y();
                return sign_tower({OneRoot(arc.c.x()) - p.x, OneRoot(FT(piece.side)), OneRoot(arc.rho * arc.rho) - e * e});
            });
            if (right == CGAL::ZERO) return true;
            if (c_hi == CGAL::ZERO) continue;  // the upper end does not count
            if (right == CGAL::POSITIVE) inside = !inside;
        }
    }
    return inside;
}

bool ExactTarget::contains(const ExactPoint& p, const FilterMode mode, DecisionCounters& counters) const
{
    const Interval px = interval_of(p.x), py = interval_of(p.y);
    for (std::size_t i = 0; i < vertex_disks.size(); ++i) {
        if (!overlaps(vertex_boxes[i], px, py)) continue;
        const VertexDisk& disk = vertex_disks[i];
        const Interval dx = px - disk.x_i, dy = py - disk.y_i;
        const CGAL::Sign s = filtered(dx * dx + dy * dy - tool_radius_squared_i, mode, counters, [&] {
            return sign_tower(squared_distance_minus(p, disk.x, disk.y, tool_radius_squared));
        });
        if (s != CGAL::POSITIVE) return true;
    }
    for (const auto& band : line_bands) {
        if (!overlaps(band.bbox, px, py)) continue;
        const Interval nx = interval_of(band.nx), ny = interval_of(band.ny), ux = interval_of(band.ux),
                       uy = interval_of(band.uy);
        const Interval lv = nx * px + ny * py - interval_of(band.level);
        const auto level = [&](const FT& shift) {
            return filtered(lv + Interval(CGAL::to_interval(shift)), mode, counters, [&] {
                return sign_one_root(band.nx * p.x + band.ny * p.y - band.level + shift);
            });
        };
        if (level(FT(0)) == CGAL::POSITIVE) continue;
        if (level(tool_radius) == CGAL::NEGATIVE) continue;
        const Interval along = ux * px + uy * py;
        const auto along_minus = [&](const OneRoot& bound) {
            return sign_difference_one_roots(band.ux * p.x + band.uy * p.y, bound);
        };
        if (filtered(along - interval_of(band.s0), mode, counters, [&] { return along_minus(band.s0); }) == CGAL::NEGATIVE) {
            continue;
        }
        if (filtered(along - interval_of(band.s1), mode, counters, [&] { return along_minus(band.s1); }) == CGAL::POSITIVE) {
            continue;
        }
        return true;
    }
    for (const auto& band : circle_bands) {
        if (!overlaps(band.bbox, px, py)) continue;
        const Interval dx = px - band.cx_i, dy = py - band.cy_i, d2 = dx * dx + dy * dy;
        const auto d2_minus = [&](const FT& r2) {
            const OneRoot ex = p.x - band.c.x(), ey = p.y - band.c.y();
            return sign_one_root(ex * ex + ey * ey - r2);
        };
        if (filtered(d2 - band.outer_sq_i, mode, counters, [&] { return d2_minus(band.outer_sq); }) == CGAL::POSITIVE) continue;
        if (filtered(d2 - band.inner_sq_i, mode, counters, [&] { return d2_minus(band.inner_sq); }) == CGAL::NEGATIVE) continue;
        const CGAL::Sign o0 = times(band.sense, filtered(cross_about_i(band.cx_i, band.cy_i, band.w0.x_i, band.w0.y_i, px, py),
                                                         mode, counters, [&] {
                                                             return sign_tower(cross_about(band.c, band.w0.x, band.w0.y, p));
                                                         }));
        const CGAL::Sign o1 = CGAL::opposite(
            times(band.sense, filtered(cross_about_i(band.cx_i, band.cy_i, band.w1.x_i, band.w1.y_i, px, py), mode, counters,
                                       [&] { return sign_tower(cross_about(band.c, band.w1.x, band.w1.y, p)); })));
        if (in_wedge(band.span, o0, o1)) return true;
    }
    return in_centre_domain(p, mode, counters);
}

Point_2 ExactTarget::interior_anchor(DecisionCounters& counters) const
{
    if (line_bands.empty()) {
        throw InvalidCertificateInputError("the exact target has no offset-line piece to anchor an interior point on");
    }
    const LineBand& band = line_bands.front();
    const auto as_tower = [](const OneRoot& v) {
        if (!v.is_extended()) return Tower::rational(v.a0());
        return Tower{OneRoot(v.a0()), OneRoot(v.a1()), OneRoot(v.root())};
    };
    const FT tau = rational_between({false, FT(0), as_tower(band.s0)}, {false, FT(0), as_tower(band.s1)}, counters);
    const FT level = band.level - tool_radius / FT(2);
    return Point_2(tau * band.ux + level * band.nx, tau * band.uy + level * band.ny);
}

}  // namespace compas_cgal::coverage_certificate
