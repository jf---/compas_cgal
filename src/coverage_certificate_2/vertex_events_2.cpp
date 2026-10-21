#include "coverage_certificate_2/vertex_events_2.h"

#include "coverage_certificate_2/errors_2.h"

#include <algorithm>
#include <cmath>

namespace compas_cgal::coverage_certificate {

namespace {
bool overlaps(const Interval& a, const Interval& b) { return !(a.sup() < b.inf() || b.sup() < a.inf()); }

Interval sqrt_nonnegative(const Interval& v)
{
    return CGAL::sqrt(Interval(v.inf() < 0 ? 0.0 : v.inf(), v.sup() < 0 ? 0.0 : v.sup()));
}

FT part0(const OneRoot& v) { return v.a0(); }
FT part1(const OneRoot& v) { return v.is_extended() ? v.a1() : FT(0); }

/// a / b for one-roots over one radicand: a conj(b) / N(b), N(b) = b0^2 - b1^2 delta.
OneRoot divide_one_roots(const OneRoot& a, const OneRoot& b, const FT& delta)
{
    const FT b0 = part0(b), b1 = part1(b);
    const FT norm = b0 * b0 - b1 * b1 * delta;
    if (CGAL::sign(norm) == CGAL::ZERO) throw CoverageCertificateNativeError("one-root division by a zero-norm value");
    const OneRoot conj = make_one_root(b0, -b1, delta);
    return divide(a * conj, norm);
}

void require_same_delta(const FT& a, const FT& b)
{
    if (CGAL::sign(a) != CGAL::ZERO && CGAL::sign(b) != CGAL::ZERO && a != b) {
        throw CoverageCertificateNativeError("tower-plus-root arithmetic across two radicands");
    }
}
}  // namespace

TowerPlusRoot TowerPlusRoot::of(const Tower& t, const FT& delta) { return {t, Tower::rational(FT(0)), delta}; }

TowerPlusRoot TowerPlusRoot::of(const OneRoot& v, const FT& delta)
{
    if (v.is_extended() && CGAL::sign(v.root()) != CGAL::ZERO && v.root() != delta) {
        throw CoverageCertificateNativeError("one-root over a radicand other than the vertex's");
    }
    return {Tower::rational(part0(v)), Tower::rational(part1(v)), delta};
}

TowerPlusRoot operator+(const TowerPlusRoot& a, const TowerPlusRoot& b)
{
    require_same_delta(a.delta, b.delta);
    return {tower_add(a.P, b.P), tower_add(a.Q, b.Q), a.delta};
}

TowerPlusRoot operator-(const TowerPlusRoot& a, const TowerPlusRoot& b)
{
    require_same_delta(a.delta, b.delta);
    return {tower_sub(a.P, b.P), tower_sub(a.Q, b.Q), a.delta};
}

TowerPlusRoot operator*(const TowerPlusRoot& a, const TowerPlusRoot& b)
{
    require_same_delta(a.delta, b.delta);
    const Tower qq = tower_scale(tower_mul(a.Q, b.Q), OneRoot(a.delta));
    return {tower_add(tower_mul(a.P, b.P), qq), tower_add(tower_mul(a.P, b.Q), tower_mul(a.Q, b.P)), a.delta};
}

CGAL::Sign vertex_circle_side(const Tower& x, const Tower& y, const VertexCircle& V)
{
    const TowerPlusRoot dx = TowerPlusRoot::of(x, V.delta) - TowerPlusRoot::of(V.x, V.delta);
    const TowerPlusRoot dy = TowerPlusRoot::of(y, V.delta) - TowerPlusRoot::of(V.y, V.delta);
    return (dx * dx + dy * dy - TowerPlusRoot::of(OneRoot(V.r_sq), V.delta)).sign();
}

// ---------------------------------------------------------------------------
// On a rational circle
// ---------------------------------------------------------------------------

CircleVertexExact exact_circle_vertex(const BoundaryCircle& circle, const VertexCircle& V)
{
    const Point_2& c = circle.circle.center();
    const FT& rho_sq = circle.circle.squared_radius();
    CircleVertexExact out;
    out.ex = V.x - c.x();
    out.ey = V.y - c.y();
    const OneRoot D = out.ex * out.ex + out.ey * out.ey;
    if (sign_one_root(D) == CGAL::ZERO) throw CoverageCertificateNativeError("a vertex circle is concentric with a rational circle");
    const OneRoot q = D + OneRoot(rho_sq - V.r_sq);
    const OneRoot disc = FT(4) * rho_sq * D - q * q;
    const CGAL::Sign s = sign_one_root(disc);
    if (s == CGAL::NEGATIVE) return out;
    const OneRoot two_d = FT(2) * D;
    const OneRoot l = divide_one_roots(q, two_d, V.delta);
    out.mx = l * out.ex + c.x();
    out.my = l * out.ey + c.y();
    if (s == CGAL::ZERO) {
        out.count = 1;
        for (std::size_t b = 0; b < 2; ++b) {
            out.x[b] = Tower::of(out.mx);
            out.y[b] = Tower::of(out.my);
        }
        return out;
    }
    out.count = 2;
    const OneRoot nx = divide_one_roots(-out.ey, two_d, V.delta);
    const OneRoot ny = divide_one_roots(out.ex, two_d, V.delta);
    for (std::size_t b = 0; b < 2; ++b) {
        const FT sigma = b == 0 ? FT(1) : FT(-1);
        out.x[b] = {out.mx, sigma * nx, disc};
        out.y[b] = {out.my, sigma * ny, disc};
    }
    return out;
}

CircleVertexRecord make_circle_vertex_record(const std::uint32_t circle_id, const BoundaryCircle& circle,
                                             const std::uint32_t vertex_id, const VertexCircle& V, const FilterMode mode,
                                             DecisionCounters& counters)
{
    ++counters.vertex_circle_pairs;
    CircleVertexRecord record{circle_id, vertex_id, 0, {}, {}, nullptr};
    if (mode == FilterMode::Lazy) {
        const Interval ex = V.x_i - circle.centre_x, ey = V.y_i - circle.centre_y;
        const Interval D = ex * ex + ey * ey;
        const Interval q = D + circle.radius_squared_i - V.r_sq_i;
        const Interval disc = Interval(4) * D * circle.radius_squared_i - q * q;
        if (disc.sup() < 0) return record;
        if (disc.inf() > 0) {
            const Interval l = q / (Interval(2) * D);
            const Interval s = CGAL::sqrt(disc) / (Interval(2) * D);
            const Interval mx = circle.centre_x + l * ex, my = circle.centre_y + l * ey;
            record.count = 2;
            record.x = {mx - s * ey, mx + s * ey};
            record.y = {my + s * ex, my - s * ex};
            return record;
        }
    }
    ++counters.vertex_circle_exact;
    CircleVertexExact exact = exact_circle_vertex(circle, V);
    record.count = exact.count;
    if (exact.count == 1) ++counters.vertex_circle_tangencies;
    for (std::size_t b = 0; b < 2 && exact.count > 0; ++b) {
        record.x[b] = interval_of(exact.x[b]);
        record.y[b] = interval_of(exact.y[b]);
    }
    record.exact = std::make_unique<CircleVertexExact>(std::move(exact));
    return record;
}

const CircleVertexExact& exact_of(CircleVertexRecord& record, const BoundaryCircle& circle, const VertexCircle& V,
                                  DecisionCounters& counters)
{
    if (record.exact) return *record.exact;
    ++counters.materialised_vertex_events;
    CircleVertexExact exact = exact_circle_vertex(circle, V);
    if (exact.count != record.count) {
        throw CoverageCertificateNativeError("a filter-certified circle-vertex count disagrees with the exact count");
    }
    for (std::size_t b = 0; b < exact.count; ++b) {
        if (!overlaps(interval_of(exact.x[b]), record.x[b]) || !overlaps(interval_of(exact.y[b]), record.y[b])) {
            throw CoverageCertificateNativeError("an exact circle-vertex crossing lies outside its certified enclosure");
        }
    }
    record.exact = std::make_unique<CircleVertexExact>(std::move(exact));
    return *record.exact;
}

bool is_circle_vertex_crossing(const Tower& x, const Tower& y, const CircleVertexExact& exact, const std::uint8_t branch,
                               const VertexCircle& V)
{
    if (exact.count == 0 || vertex_circle_side(x, y, V) != CGAL::ZERO) return false;
    if (exact.count == 1) return true;
    // (X - m) . perp(e) is +s |e|^2 on the + branch, -s |e|^2 on the other.
    const TowerPlusRoot side = (TowerPlusRoot::of(x, V.delta) - TowerPlusRoot::of(exact.mx, V.delta)) *
                                   TowerPlusRoot::of(-exact.ey, V.delta) +
                               (TowerPlusRoot::of(y, V.delta) - TowerPlusRoot::of(exact.my, V.delta)) *
                                   TowerPlusRoot::of(exact.ex, V.delta);
    return side.sign() == (branch == 0 ? CGAL::POSITIVE : CGAL::NEGATIVE);
}

// ---------------------------------------------------------------------------
// On a line support
// ---------------------------------------------------------------------------

CGAL::Sign LineVertexExact::sign_minus(const std::uint8_t branch, const FT& q) const
{
    if (count == 0) throw CoverageCertificateNativeError("sign of a vertex root on a line it does not cross");
    const FT ma = part0(mu) - q * ld;
    const FT mb = part1(mu);
    const OneRoot M = make_one_root(ma, mb, delta);
    const CGAL::Sign sm = sign_one_root(M);
    const Interval disc_i = interval_of(disc_a) + interval_of(disc_b) * CGAL::sqrt(interval_of(delta));
    if (count == 1) return require_sign_agrees(sm, interval_of(M), "line-vertex tangent root");
    const CGAL::Sign sigma = branch == 0 ? CGAL::POSITIVE : CGAL::NEGATIVE;
    CGAL::Sign exact = sigma;
    if (sm != CGAL::ZERO && sm != sigma) {
        // compare M^2 with disc: (ma^2 + mb^2 delta - disc_a) + (2 ma mb - disc_b) sqrt(delta).
        const Tower diff{OneRoot(ma * ma + mb * mb * delta) - disc_a, OneRoot(FT(2) * ma * mb) - disc_b, OneRoot(delta)};
        const CGAL::Sign c = sign_tower(diff);
        exact = c == CGAL::POSITIVE ? sm : (c == CGAL::NEGATIVE ? sigma : CGAL::ZERO);
    }
    const Interval root = sqrt_nonnegative(disc_i);
    return require_sign_agrees(exact, interval_of(M) + (branch == 0 ? root : -root), "line-vertex root");
}

Interval KNumber::interval() const { return interval_of(a) + interval_of(b) * CGAL::sqrt(interval_of(delta)); }

namespace {
void require_same_k(const KNumber& x, const KNumber& y)
{
    if (x.delta != y.delta) throw CoverageCertificateNativeError("K arithmetic across two vertex radicands");
}
}  // namespace

KNumber operator+(const KNumber& x, const KNumber& y)
{
    require_same_k(x, y);
    require_same_base(x.a, y.a);
    return {x.a + y.a, x.b + y.b, x.delta};
}

KNumber operator-(const KNumber& x, const KNumber& y)
{
    require_same_k(x, y);
    require_same_base(x.a, y.a);
    return {x.a - y.a, x.b - y.b, x.delta};
}

KNumber operator*(const KNumber& x, const KNumber& y)
{
    require_same_k(x, y);
    require_same_base(x.a, y.a);
    return {x.a * y.a + x.b * y.b * x.delta, x.a * y.b + x.b * y.a, x.delta};
}

CGAL::Sign sign_k_root(const KNumber& A, const KNumber& B, const KNumber& D)
{
    const CGAL::Sign sd = D.sign();
    if (sd == CGAL::NEGATIVE) throw CoverageCertificateNativeError("K-root radicand is negative");
    const CGAL::Sign sa = A.sign();
    if (sd == CGAL::ZERO) return sa;
    const CGAL::Sign sb = B.sign();
    CGAL::Sign exact = CGAL::ZERO;
    if (sb == CGAL::ZERO) {
        exact = sa;
    } else if (sa == CGAL::ZERO || sa == sb) {
        exact = sb;
    } else {
        const CGAL::Sign c = (A * A - B * B * D).sign();
        exact = c == CGAL::POSITIVE ? sa : (c == CGAL::NEGATIVE ? sb : CGAL::ZERO);
    }
    return require_sign_agrees(exact, A.interval() + B.interval() * sqrt_nonnegative(D.interval()), "K root");
}

CGAL::Sign LineVertexExact::sign_minus(const std::uint8_t branch, const Tower& t) const
{
    if (count == 0) throw CoverageCertificateNativeError("sign of a vertex root on a line it does not cross");
    // ld (t_V - t) = (mu - ld m) + sigma sqrt(disc) - ld n sqrt(d), ld > 0.
    const KNumber A{OneRoot(part0(mu)) - ld * t.m, OneRoot(part1(mu)), delta};
    const KNumber C{-ld * t.n, OneRoot(FT(0)), delta};
    const KNumber d{t.d, OneRoot(FT(0)), delta};
    const KNumber disc{disc_a, disc_b, delta};
    const KNumber sigma{OneRoot(FT(branch == 0 ? 1 : -1)), OneRoot(FT(0)), delta};
    if (count == 1) return sign_k_root(A, C, d);
    const CGAL::Sign s1 = sign_k_root(A, sigma, disc);
    const CGAL::Sign sc = C.sign();
    if (sc == CGAL::ZERO || d.sign() == CGAL::ZERO) return s1;
    if (s1 == CGAL::ZERO) return sc;
    if (s1 == sc) return s1;
    // compare (A + sigma sqrt(disc))^2 with C^2 d.
    const CGAL::Sign c = sign_k_root(A * A + disc - C * C * d, sigma * A + sigma * A, disc);
    return c == CGAL::POSITIVE ? s1 : (c == CGAL::NEGATIVE ? sc : CGAL::ZERO);
}

SeparationMeasure LineVertexExact::measure() const
{
    const SeparationMeasure disc = plus(SeparationMeasure::of(disc_a), times(SeparationMeasure::of(disc_b),
                                                                             square_root(SeparationMeasure::of(delta))));
    return over(plus(SeparationMeasure::of(mu), square_root(disc)), SeparationMeasure::of(ld));
}

LineVertexExact exact_line_vertex(const LineSupport& line, const VertexCircle& V)
{
    LineVertexExact out;
    out.delta = V.delta;
    out.ld = line.nn;
    const FT dx = line.dx(), dy = line.dy();
    const FT bwx = line.b0x - part0(V.x), bwy = line.b0y - part0(V.y);
    const FT w1x = part1(V.x), w1y = part1(V.y);
    const FT beta_a = dx * bwx + dy * bwy;
    const FT beta_b = -(dx * w1x + dy * w1y);
    out.mu = make_one_root(-beta_a, -beta_b, V.delta);
    const FT b1b1 = line.b1x * line.b1x + line.b1y * line.b1y;
    const OneRoot gamma_a = make_one_root(bwx * bwx + bwy * bwy + b1b1 * line.L + (w1x * w1x + w1y * w1y) * V.delta - V.r_sq,
                                          FT(2) * (bwx * line.b1x + bwy * line.b1y), line.L);
    const OneRoot gamma_b = make_one_root(FT(-2) * (bwx * w1x + bwy * w1y), FT(-2) * (line.b1x * w1x + line.b1y * w1y), line.L);
    out.disc_a = OneRoot(beta_a * beta_a + beta_b * beta_b * V.delta) - out.ld * gamma_a;
    out.disc_b = OneRoot(FT(2) * beta_a * beta_b) - out.ld * gamma_b;
    const CGAL::Sign s = sign_tower({out.disc_a, out.disc_b, OneRoot(V.delta)});
    out.count = s == CGAL::NEGATIVE ? 0 : (s == CGAL::ZERO ? 1 : 2);
    return out;
}

namespace {
/// Filter policy (never a decision): a lazy root enclosure wider than this,
/// relative to 1 + |t|, is replaced by the exact-parts enclosure at once. 2^-48
/// is above the ~2^-52 width of a well-conditioned root and far below any
/// gap the sweep needs to order distinct events cheaply.
constexpr double VERTEX_ROOT_RELATIVE_WIDTH = 0x1p-40;

bool too_wide(const Interval& t)
{
    const double scale = 1.0 + std::max(std::abs(t.inf()), std::abs(t.sup()));
    return t.sup() - t.inf() > VERTEX_ROOT_RELATIVE_WIDTH * scale;
}

/// The same for a + b sqrt(delta) in K, a, b over sqrt(L).
Interval tight_interval_k(const OneRoot& a, const OneRoot& b, const FT& delta)
{
    const Interval ia = tight_interval(a);
    const Interval ib = tight_interval(b);
    const Interval root = CGAL::sqrt(interval_of(delta));
    const bool cancel = (ia.inf() > 0 && ib.sup() < 0) || (ia.sup() < 0 && ib.inf() > 0);
    if (!cancel) return ia + ib * root;
    const OneRoot norm = a * a - b * b * delta;
    return tight_interval(norm) / (ia - ib * root);
}

std::array<Interval, 2> root_intervals(const LineVertexExact& e)
{
    const Interval disc = tight_interval_k(e.disc_a, e.disc_b, e.delta);
    const Interval root = sqrt_nonnegative(disc);
    const Interval mu = tight_interval(e.mu), ld = interval_of(e.ld);
    if (e.count == 1) return {mu / ld, mu / ld};
    return {(mu + root) / ld, (mu - root) / ld};
}
}  // namespace

LineVertexRecord make_line_vertex_record(const LineSupport& line, const std::uint32_t vertex_id, const VertexCircle& V,
                                         const FilterMode mode, DecisionCounters& counters)
{
    ++counters.vertex_line_pairs;
    LineVertexRecord record{vertex_id, 0, {}, nullptr};
    if (mode == FilterMode::Lazy) {
        const Interval x0 = line.b0x_i + line.b1x_i * line.sqrt_l_i - V.x_i;
        const Interval y0 = line.b0y_i + line.b1y_i * line.sqrt_l_i - V.y_i;
        const Interval dx = line.ny_i, dy = -line.nx_i;
        const Interval beta = dx * x0 + dy * y0;
        const Interval gamma = x0 * x0 + y0 * y0 - V.r_sq_i;
        const Interval disc = beta * beta - line.nn_i * gamma;
        if (disc.sup() < 0) return record;
        if (disc.inf() > 0) {
            const Interval root = CGAL::sqrt(disc);
            record.count = 2;
            record.t = {(-beta + root) / line.nn_i, (-beta - root) / line.nn_i};
            if (!too_wide(record.t[0]) && !too_wide(record.t[1])) return record;
            // Near tangency the double discriminant has lost its digits; the
            // exact parts cancel exactly, so materialise now for a tight enclosure.
        }
    }
    ++counters.vertex_line_exact;
    LineVertexExact exact = exact_line_vertex(line, V);
    if (record.count == 2 && exact.count != 2) {
        throw CoverageCertificateNativeError("a filter-certified line-vertex count disagrees with the exact count");
    }
    const auto lazy_t = record.t;
    const bool had_lazy = record.count == 2;
    record.count = exact.count;
    if (exact.count == 1) ++counters.vertex_circle_tangencies;
    if (exact.count > 0) record.t = root_intervals(exact);
    for (std::size_t b = 0; had_lazy && b < 2; ++b) {
        if (!overlaps(lazy_t[b], record.t[b])) {
            throw CoverageCertificateNativeError("an exact line-vertex root lies outside its certified enclosure");
        }
        record.t[b] = Interval(std::max(lazy_t[b].inf(), record.t[b].inf()), std::min(lazy_t[b].sup(), record.t[b].sup()));
    }
    record.exact = std::make_unique<LineVertexExact>(std::move(exact));
    return record;
}

const LineVertexExact& exact_of(LineVertexRecord& record, const LineSupport& line, const VertexCircle& V,
                                DecisionCounters& counters)
{
    if (record.exact) return *record.exact;
    ++counters.materialised_vertex_events;
    LineVertexExact exact = exact_line_vertex(line, V);
    if (exact.count != record.count) {
        throw CoverageCertificateNativeError("a filter-certified line-vertex count disagrees with the exact count");
    }
    const auto t = root_intervals(exact);
    for (std::size_t b = 0; b < exact.count; ++b) {
        if (!overlaps(t[b], record.t[b])) {
            throw CoverageCertificateNativeError("an exact line-vertex root lies outside its certified enclosure");
        }
        // Both are certified enclosures of the same root: keep their intersection
        // (the exact parts have the cancellation of the discriminant done exactly).
        record.t[b] = Interval(std::max(t[b].inf(), record.t[b].inf()), std::min(t[b].sup(), record.t[b].sup()));
    }
    record.exact = std::make_unique<LineVertexExact>(std::move(exact));
    return *record.exact;
}

bool is_line_vertex_root(const Tower& t, const LineSupport& line, const LineVertexExact& exact, const std::uint8_t branch,
                         const VertexCircle& V)
{
    if (exact.count == 0) return false;
    const OneRoot bx = make_one_root(line.b0x, line.b1x, line.L);
    const OneRoot by = make_one_root(line.b0y, line.b1y, line.L);
    const Tower x{bx + t.m * line.dx(), t.n * line.dx(), t.d};
    const Tower y{by + t.m * line.dy(), t.n * line.dy(), t.d};
    if (vertex_circle_side(x, y, V) != CGAL::ZERO) return false;
    if (exact.count == 1) return true;
    // The + root lies beyond the chord midpoint mu / Ld.
    const Tower shifted{t.m - part0(exact.mu) / exact.ld, t.n, t.d};
    const CGAL::Sign side = sign_tower_plus(shifted, OneRoot(-part1(exact.mu) / exact.ld), OneRoot(exact.delta));
    return side == (branch == 0 ? CGAL::POSITIVE : CGAL::NEGATIVE);
}

}  // namespace compas_cgal::coverage_certificate
