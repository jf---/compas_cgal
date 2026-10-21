#include "coverage_certificate_2/line_crossings_2.h"

#include "coverage_certificate_2/errors_2.h"

#include <algorithm>
#include <cmath>
#include <utility>

namespace compas_cgal::coverage_certificate {

namespace {
bool overlaps(const Interval& a, const Interval& b) { return !(a.sup() < b.inf() || b.sup() < a.inf()); }

/// Filter policy (never a decision): a lazy root enclosure wider than this,
/// relative to 1 + |t|, is replaced by the exact towers' enclosure at once
/// (same constant and reasoning as vertex_events_2.cpp).
constexpr double ROOT_RELATIVE_WIDTH = 0x1p-40;

bool too_wide(const Interval& t)
{
    const double scale = 1.0 + std::max(std::abs(t.inf()), std::abs(t.sup()));
    return t.sup() - t.inf() > ROOT_RELATIVE_WIDTH * scale;
}

/// The crossing parameter on P of a non-parallel line Q (den = n_Q . d_P != 0).
Tower line_crossing_tower(const LineSupport& piece, const LineSupport& other, const FT& den)
{
    return {make_one_root((other.k0 - other.nx * piece.b0x - other.ny * piece.b0y) / den,
                          -(other.nx * piece.b1x + other.ny * piece.b1y) / den, piece.L),
            OneRoot(CGAL::sign(other.L) == CGAL::ZERO ? FT(0) : other.k1 / den), OneRoot(other.L)};
}
}  // namespace

CircleLineExact exact_circle_line(const BoundaryCircle& circle, const LineSupport& line)
{
    const Point_2& c = circle.circle.center();
    const FT& rho_sq = circle.circle.squared_radius();
    const FT& ld = line.nn;
    const FT ex = line.b0x - c.x();
    const FT ey = line.b0y - c.y();
    const FT beta = line.dx() * ex + line.dy() * ey;
    const FT g0 = ex * ex + ey * ey + (line.b1x * line.b1x + line.b1y * line.b1y) * line.L - rho_sq;
    const FT g1 = FT(2) * (ex * line.b1x + ey * line.b1y);
    const OneRoot disc = make_one_root(beta * beta - ld * g0, -ld * g1, line.L);
    CircleLineExact out;
    const CGAL::Sign s = sign_one_root(disc);
    if (s == CGAL::NEGATIVE) return out;
    const FT t_mid = -beta / ld;
    const OneRoot mx = make_one_root(line.b0x + t_mid * line.dx(), line.b1x, line.L);
    const OneRoot my = make_one_root(line.b0y + t_mid * line.dy(), line.b1y, line.L);
    if (s == CGAL::ZERO) {
        out.count = 1;
        for (std::size_t b = 0; b < 2; ++b) {
            out.t[b] = Tower::rational(t_mid);
            out.x[b] = Tower::of(mx);
            out.y[b] = Tower::of(my);
        }
        return out;
    }
    out.count = 2;
    for (std::size_t b = 0; b < 2; ++b) {
        const FT unit = (b == 0 ? FT(1) : FT(-1)) / ld;
        out.t[b] = {OneRoot(t_mid), OneRoot(unit), disc};
        out.x[b] = {mx, OneRoot(unit * line.dx()), disc};
        out.y[b] = {my, OneRoot(unit * line.dy()), disc};
    }
    return out;
}

CircleLineFilter filter_circle_line(const BoundaryCircle& circle, const LineSupport& line)
{
    const Interval ex = line.b0x_i - circle.centre_x;
    const Interval ey = line.b0y_i - circle.centre_y;
    const Interval dx = line.ny_i;
    const Interval dy = -line.nx_i;
    const Interval& ld = line.nn_i;
    const Interval beta = dx * ex + dy * ey;
    const Interval l = line.sqrt_l_i * line.sqrt_l_i;
    const Interval g0 = ex * ex + ey * ey + (line.b1x_i * line.b1x_i + line.b1y_i * line.b1y_i) * l - circle.radius_squared_i;
    const Interval g1 = Interval(2) * (ex * line.b1x_i + ey * line.b1y_i);
    const Interval disc = beta * beta - ld * g0 - ld * g1 * line.sqrt_l_i;
    CircleLineFilter out;
    if (disc.sup() < 0) {
        out.count = 0;
        return out;
    }
    if (!(disc.inf() > 0)) return out;
    out.count = 2;
    const Interval root = CGAL::sqrt(disc);
    for (std::size_t b = 0; b < 2; ++b) {
        const Interval t = (b == 0 ? (-beta + root) : (-beta - root)) / ld;
        out.t[b] = t;
        out.x[b] = line.b0x_i + line.b1x_i * line.sqrt_l_i + t * dx;
        out.y[b] = line.b0y_i + line.b1y_i * line.sqrt_l_i + t * dy;
    }
    return out;
}

CircleLineRecord make_circle_line_record(const std::uint32_t circle_id, const BoundaryCircle& circle,
                                         const std::uint32_t line_id, const LineSupport& line, const FilterMode mode,
                                         DecisionCounters& counters)
{
    ++counters.circle_line_pairs;
    CircleLineRecord record{circle_id, line_id, 0, {}, {}, {}, nullptr};
    std::optional<CircleLineFilter> lazy;
    if (mode == FilterMode::Lazy) {
        const CircleLineFilter f = filter_circle_line(circle, line);
        if (f.count) {
            record.count = *f.count;
            record.t = f.t;
            record.x = f.x;
            record.y = f.y;
            if (record.count == 0 || (!too_wide(f.t[0]) && !too_wide(f.t[1]))) return record;
            // Near tangency the double discriminant has lost its digits: build the
            // exact towers now, whose enclosures are tight (the count stays the
            // filter's certified one).
            lazy = f;
        }
    }
    ++counters.circle_line_exact;
    CircleLineExact exact = exact_circle_line(circle, line);
    if (lazy && exact.count != *lazy->count) {
        throw CoverageCertificateNativeError("a filter-certified circle-line count disagrees with the exact count");
    }
    record.count = exact.count;
    if (exact.count == 1) ++counters.line_tangencies;
    for (std::size_t b = 0; b < exact.count && b < 2; ++b) {
        record.t[b] = lazy ? tight_interval(exact.t[b]) : interval_of(exact.t[b]);
        record.x[b] = lazy ? tight_interval(exact.x[b]) : interval_of(exact.x[b]);
        record.y[b] = lazy ? tight_interval(exact.y[b]) : interval_of(exact.y[b]);
    }
    if (exact.count == 1) {
        record.t[1] = record.t[0];
        record.x[1] = record.x[0];
        record.y[1] = record.y[0];
    }
    if (lazy) {
        // Both are certified enclosures: keep the intersection.
        for (std::size_t b = 0; b < 2; ++b) {
            for (auto [mine, theirs] : {std::pair{&record.t[b], &lazy->t[b]}, std::pair{&record.x[b], &lazy->x[b]},
                                        std::pair{&record.y[b], &lazy->y[b]}}) {
                if (!overlaps(*mine, *theirs)) {
                    throw CoverageCertificateNativeError("an exact circle-line crossing lies outside its certified enclosure");
                }
                *mine = Interval(std::max(mine->inf(), theirs->inf()), std::min(mine->sup(), theirs->sup()));
            }
        }
    }
    record.exact = std::make_unique<CircleLineExact>(std::move(exact));
    return record;
}

const CircleLineExact& exact_of(CircleLineRecord& record, const BoundaryCircle& circle, const LineSupport& line,
                                DecisionCounters& counters)
{
    if (record.exact) return *record.exact;
    ++counters.materialised_line_events;
    CircleLineExact exact = exact_circle_line(circle, line);
    if (exact.count != record.count) {
        throw CoverageCertificateNativeError("a filter-certified circle-line count disagrees with the exact count");
    }
    for (std::size_t b = 0; b < exact.count; ++b) {
        if (!overlaps(interval_of(exact.t[b]), record.t[b]) || !overlaps(interval_of(exact.x[b]), record.x[b]) ||
            !overlaps(interval_of(exact.y[b]), record.y[b])) {
            throw CoverageCertificateNativeError("an exact circle-line crossing lies outside its certified enclosure");
        }
    }
    record.exact = std::make_unique<CircleLineExact>(std::move(exact));
    return *record.exact;
}

LineLineRecord make_line_line_record(const LineSupport& piece, const std::uint32_t other_id, const LineSupport& other,
                                     const bool same_support, const FilterMode mode, DecisionCounters& counters)
{
    LineLineRecord record{other_id, LineRelation::Crossing, 0, CGAL::ZERO, Interval(0), nullptr};
    if (same_support) {
        record.relation = LineRelation::Coincident;
        record.coincidence = 1;
        return record;
    }
    const Interval den_i = other.nx_i * piece.ny_i + other.ny_i * (-piece.nx_i);
    const Interval proj_i = other.nx_i * (piece.b0x_i + piece.b1x_i * piece.sqrt_l_i) +
                            other.ny_i * (piece.b0y_i + piece.b1y_i * piece.sqrt_l_i);
    if (mode == FilterMode::Lazy && (den_i.inf() > 0 || den_i.sup() < 0)) {
        record.t = (other.offset_i - proj_i) / den_i;
        if (!too_wide(record.t)) return record;
        // Nearly parallel: build the exact tower now for a tight enclosure.
        const FT den = other.nx * piece.dx() + other.ny * piece.dy();
        record.exact_t = std::make_unique<Tower>(line_crossing_tower(piece, other, den));
        const Interval exact_t = tight_interval(*record.exact_t);
        if (!overlaps(exact_t, record.t)) {
            throw CoverageCertificateNativeError("an exact line crossing lies outside its certified enclosure");
        }
        record.t = Interval(std::max(exact_t.inf(), record.t.inf()), std::min(exact_t.sup(), record.t.sup()));
        return record;
    }
    ++counters.line_line_exact;
    const FT den = other.nx * piece.dx() + other.ny * piece.dy();
    if (CGAL::sign(den) != CGAL::ZERO) {
        record.exact_t = std::make_unique<Tower>(line_crossing_tower(piece, other, den));
        record.t = interval_of(*record.exact_t);
        return record;
    }
    // Parallel: P's side relative to Q is constant; zero means one line.
    const OneRoot proj = make_one_root(other.nx * piece.b0x + other.ny * piece.b0y - other.k0,
                                       other.nx * piece.b1x + other.ny * piece.b1y, piece.L);
    record.parallel_side = sign_tower({proj, OneRoot(-other.k1), OneRoot(other.L)});
    if (record.parallel_side == CGAL::ZERO) {
        record.relation = LineRelation::Coincident;
        record.coincidence = CGAL::sign(other.nx * piece.nx + other.ny * piece.ny) == CGAL::POSITIVE ? 1 : -1;
    } else {
        record.relation = LineRelation::Parallel;
    }
    return record;
}

const Tower& exact_t_of(LineLineRecord& record, const LineSupport& piece, const LineSupport& other,
                        DecisionCounters& counters)
{
    if (record.relation != LineRelation::Crossing) throw CoverageCertificateNativeError("exact crossing of non-crossing lines");
    if (record.exact_t) return *record.exact_t;
    ++counters.materialised_line_events;
    const FT den = other.nx * piece.dx() + other.ny * piece.dy();
    if (CGAL::sign(den) == CGAL::ZERO) throw CoverageCertificateNativeError("a filter-certified crossing pair is parallel");
    const Tower t = line_crossing_tower(piece, other, den);
    if (!overlaps(interval_of(t), record.t)) {
        throw CoverageCertificateNativeError("an exact line crossing lies outside its certified enclosure");
    }
    record.exact_t = std::make_unique<Tower>(t);
    return *record.exact_t;
}

}  // namespace compas_cgal::coverage_certificate
