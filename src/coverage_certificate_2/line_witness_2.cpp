#include "coverage_certificate_2/line_witness_2.h"

#include "coverage_certificate_2/errors_2.h"

#include <string>

namespace compas_cgal::coverage_certificate {

namespace {
/// Contract bound on P7 halvings (see witness_2.cpp): termination is a theorem,
/// the bound turns a logic error into a loud failure.
constexpr int UNCOVERED_HALVING_CONTRACT_BOUND = 4096;

/// Contract bound on widening the initial rational bracket of the flank offset
/// (each step widens by 1 mm; the certified enclosure is already a bracket, so a
/// correct run never widens at all).
constexpr int BRACKET_WIDENING_CONTRACT_BOUND = 64;

/// Whether a counted rectangle covers the witness w on the piece: every line not
/// coincident with the piece strictly on its inward side (the coincident one is
/// satisfied by construction). A non-coincident line through w is a breach.
bool rect_covers_at(const LineSweep& sweep, const std::size_t r, const RectRegion& rect, const LineTable& lines,
                    const ExactPoint& w)
{
    bool covers = true;
    for (std::uint8_t s = 0; s < 4; ++s) {
        if (s == sweep.rect_coincident_slot[r]) continue;
        const CGAL::Sign side = line_side(lines[rect.lines[s]], w);
        if (side == CGAL::ZERO) {
            throw CoverageCertificateNativeError("a flank witness lies on a rectangle line the crossings did not predict");
        }
        covers = covers && side == rect.inward[s];
    }
    return covers;
}

bool outside_rect(const RectRegion& rect, const LineTable& lines, const Point_2& q, const Interval& qx, const Interval& qy,
                  const FilterMode mode, DecisionCounters& counters)
{
    for (std::size_t s = 0; s < 4; ++s) {
        const CGAL::Sign side = line_side(lines[rect.lines[s]], q, qx, qy, mode, counters);
        if (side != CGAL::ZERO && side != rect.inward[s]) return true;
    }
    return false;
}

/// sign(q |n|^2 - k1 sqrt(L)) = sign(q - alpha) for the flank offset alpha.
CGAL::Sign offset_sign(const LineSupport& line, const FT& q)
{
    return sign_one_root(make_one_root(q * line.nn, -line.k1, line.L));
}
}  // namespace

FT line_edge_parameter(LineSweep& sweep, const std::size_t edge, DecisionCounters& counters)
{
    const auto [lo, hi] = sub_edge_bounds(sweep, edge, counters);
    return rational_between(lo, hi, counters);
}

void confirm_line_witness(LineSweep& sweep, const FT& t, LineWitnessContext& context, DecisionCounters& counters)
{
    const ExactPoint w = line_point(context.lines[sweep.piece->line], t);
    const Interval wx = interval_of(w.x), wy = interval_of(w.y);
    for (const std::uint32_t j : sweep.regions) {
        const CircleRegion& Y = context.set.regions[j];
        const auto outer = compare_to_circle(w, wx, wy, Y.centre, Y.centre_x, Y.centre_y, Y.outer_squared, Y.outer_squared_i,
                                             context.mode, counters);
        const bool exterior = outer == CGAL::LARGER ||
                              (Y.annulus && compare_to_circle(w, wx, wy, Y.centre, Y.centre_x, Y.centre_y, Y.inner_squared,
                                                              Y.inner_squared_i, context.mode, counters) == CGAL::SMALLER);
        if (!exterior) {
            throw CoverageCertificateNativeError("exposed flank witness lies in or on circle region " + std::to_string(j) +
                                                 ": the cover counters disagree with P4");
        }
    }
    for (std::size_t r = 0; r < sweep.rects.size(); ++r) {
        if (rect_covers_at(sweep, r, context.rects[sweep.rects[r]], context.lines, w)) {
            throw CoverageCertificateNativeError("exposed flank witness is covered by rectangle " +
                                                 std::to_string(sweep.rects[r]) + ": the cover counters disagree with P4");
        }
    }
}

void audit_covered_line_edge(LineSweep& sweep, const std::size_t edge, LineWitnessContext& context,
                             DecisionCounters& counters)
{
    const auto [lo, hi] = sub_edge_bounds(sweep, edge, counters);
    const FT t = rational_between(lo, hi, counters);
    const ExactPoint w = line_point(context.lines[sweep.piece->line], t);
    bool covered = false;
    for (const std::uint32_t j : sweep.regions) {
        const Side side = membership(context.set.regions[j], w);
        if (side == Side::Boundary) {
            throw CoverageCertificateNativeError("audit flank witness lies on circle region " + std::to_string(j) +
                                                 ", which the crossings did not predict");
        }
        covered = covered || side == Side::Interior;
    }
    for (std::size_t r = 0; r < sweep.rects.size(); ++r) {
        covered = rect_covers_at(sweep, r, context.rects[sweep.rects[r]], context.lines, w) || covered;
    }
    if (!covered) {
        throw CoverageCertificateNativeError("audit: flank sub-edge " + std::to_string(edge) +
                                             " is counted covered but its exact witness lies in no candidate");
    }
}

Point_2 uncovered_beside_line(LineSweep& sweep, const FT& t, LineWitnessContext& context,
                              const std::function<bool(const Point_2&)>& accept, DecisionCounters& counters)
{
    const LinePiece& piece = *sweep.piece;
    const LineSupport& line = context.lines[piece.line];
    const bool positive = piece.check_side == CGAL::POSITIVE;
    const FT px = line.b0x + t * line.dx();
    const FT py = line.b0y + t * line.dy();
    // Exact rational bracket lo < alpha < hi of the normal offset alpha.
    FT lo(0), hi(0);
    if (!line.is_rational()) {
        const Interval alpha = line.k1_i * line.sqrt_l_i / line.nn_i;
        lo = FT(alpha.inf());
        hi = FT(alpha.sup());
        for (int i = 0; offset_sign(line, lo) != CGAL::NEGATIVE; ++i) {
            if (i == BRACKET_WIDENING_CONTRACT_BOUND) throw CoverageCertificateNativeError("flank offset bracket failed");
            lo -= FT(1);
        }
        for (int i = 0; offset_sign(line, hi) != CGAL::POSITIVE; ++i) {
            if (i == BRACKET_WIDENING_CONTRACT_BOUND) throw CoverageCertificateNativeError("flank offset bracket failed");
            hi += FT(1);
        }
    }
    const FT xmin(piece.box.xmin()), xmax(piece.box.xmax()), ymin(piece.box.ymin()), ymax(piece.box.ymax());
    FT h(1);
    for (int i = 0; i < UNCOVERED_HALVING_CONTRACT_BOUND; ++i) {
        ++counters.uncovered_halvings;
        h /= FT(2);
        FT lambda;
        if (line.is_rational()) {
            lambda = positive ? h : FT(-h);
        } else {
            while (CGAL::compare(hi - lo, h) != CGAL::SMALLER) {
                const FT m = (lo + hi) / FT(2);
                if (offset_sign(line, m) == CGAL::NEGATIVE) {
                    lo = m;
                } else {
                    hi = m;
                }
            }
            lambda = positive ? hi : lo;
        }
        const Point_2 u(px + lambda * line.nx, py + lambda * line.ny);
        if (CGAL::compare(u.x(), xmin) == CGAL::SMALLER || CGAL::compare(u.x(), xmax) == CGAL::LARGER ||
            CGAL::compare(u.y(), ymin) == CGAL::SMALLER || CGAL::compare(u.y(), ymax) == CGAL::LARGER) {
            continue;
        }
        const Interval ux = interval_of(u.x()), uy = interval_of(u.y());
        bool clear = true;
        for (std::size_t s = 0; clear && s < sweep.regions.size(); ++s) {
            clear = membership(context.set.regions[sweep.regions[s]], u, ux, uy, context.mode, counters) == Side::Exterior;
        }
        for (std::size_t s = 0; clear && s < sweep.rects.size(); ++s) {
            clear = outside_rect(context.rects[sweep.rects[s]], context.lines, u, ux, uy, context.mode, counters);
        }
        for (std::size_t s = 0; clear && s < sweep.excluded_rects.size(); ++s) {
            clear = outside_rect(context.rects[sweep.excluded_rects[s]], context.lines, u, ux, uy, context.mode, counters);
        }
        if (clear && !piece.target && piece.owner_is_rect) {
            clear = outside_rect(context.rects[piece.owner], context.lines, u, ux, uy, context.mode, counters);
        }
        if (clear && accept(u)) return u;
    }
    throw CoverageCertificateNativeError("P7 uncovered-point search exceeded its contract bound on a flank");
}

}  // namespace compas_cgal::coverage_certificate
