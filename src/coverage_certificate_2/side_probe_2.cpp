#include "coverage_certificate_2/side_probe_2.h"

#include "coverage_certificate_2/errors_2.h"
#include "coverage_certificate_2/separation_2.h"
#include "coverage_certificate_2/tower_2.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace compas_cgal::coverage_certificate {

namespace {

/// The two roots of a tau^2 + b tau + c0 = 0 (a > 0 rational, b and c0 one-roots
/// over one radicand) as towers (-b / 2a) +/- (1 / 2a) sqrt(b^2 - 4 a c0), or none
/// when the discriminant is negative.
void quadratic_roots(const FT& a, const OneRoot& b, const OneRoot& c0, std::vector<Tower>& roots)
{
    const OneRoot disc = b * b - FT(4) * a * c0;
    if (sign_one_root(disc) == CGAL::NEGATIVE) return;
    const FT two_a = FT(2) * a;
    const OneRoot m = divide(-b, two_a);
    const FT half = FT(1) / two_a;
    for (const FT& n : {half, FT(-half)}) roots.push_back({m, OneRoot(n), disc});
}

}  // namespace

std::pair<FT, FT> probe_direction(const FT& dx, const FT& dy)
{
    if (CGAL::sign(dx) == CGAL::ZERO && CGAL::sign(dy) == CGAL::ZERO) {
        throw CoverageCertificateNativeError("side probe with a zero direction");
    }
    const double big = std::max(std::fabs(CGAL::to_double(dx)), std::fabs(CGAL::to_double(dy)));
    int exponent = 0;
    (void)std::frexp(big, &exponent);  // big in [2^(e-1), 2^e)
    const int shift = PROBE_DIRECTION_LOG2 - exponent;
    // 2^|shift| must be a finite binary64: fail loud rather than clamp (only a
    // direction beyond ~2^1000 mm or below ~2^-1000 mm reaches this).
    const double power = std::ldexp(1.0, shift >= 0 ? shift : -shift);
    if (!std::isfinite(power)) throw CoverageCertificateNativeError("side probe direction scale 2^" + std::to_string(shift) + " is out of binary64 range");
    const FT scale = shift >= 0 ? FT(power) : FT(1) / FT(power);
    return {dx * scale, dy * scale};
}

/// LAZY pre-pass of one support: its crossing parameters with the ray as
/// certified intervals. True when every one is decided by its enclosure (no
/// real root: discriminant enclosure below 0; a root below 0: imposes nothing; a
/// root above 0: `lower` takes its lower end); false when an enclosure meets 0,
/// and then the support takes the exact path.
bool interval_quadratic(const Interval& a, const Interval& b, const Interval& c0, double& lower)
{
    const Interval disc = b * b - Interval(4) * a * c0;
    if (disc.sup() < 0) return true;
    if (disc.inf() < 0) return false;
    const Interval root = CGAL::sqrt(disc);
    for (const Interval& tau : {(-b + root) / (Interval(2) * a), (-b - root) / (Interval(2) * a)}) {
        if (tau.sup() < 0) continue;
        if (tau.inf() > 0) {
            lower = std::min(lower, tau.inf());
            continue;
        }
        return false;
    }
    return true;
}

FT probe_step(const Point_2& w, const FT& dx, const FT& dy, const ProbeSupports& supports, const FilterMode mode,
              DecisionCounters& counters)
{
    ++counters.side_probes;
    const bool lazy = mode == FilterMode::Lazy;
    // LAZY: a support whose crossing parameters all have enclosures away from 0
    // is decided without exact work: a positive one is cleared by any h below
    // its enclosure's lower end (certified: tau >= lower > h). Only the other
    // supports form exact roots. EXACT_ONLY: every support takes the exact path.
    double lower = std::numeric_limits<double>::infinity();
    const Interval wx = interval_of(w.x()), wy = interval_of(w.y()), dx_i = interval_of(dx), dy_i = interval_of(dy);
    const Interval a_i = dx_i * dx_i + dy_i * dy_i;
    std::vector<Tower> roots;
    const FT a = dx * dx + dy * dy;
    for (const LineSupport* line : supports.lines) {
        // n.(w + tau d) = k0 + k1 sqrt(L)
        const FT nd = line->nx * dx + line->ny * dy;
        if (CGAL::sign(nd) == CGAL::ZERO) continue;  // parallel: meets the ray nowhere or along w's own line
        if (lazy) {
            const Interval nd_i = interval_of(nd);
            if (nd_i.inf() > 0 || nd_i.sup() < 0) {
                const Interval tau = (line->k0_i + line->k1_i * line->sqrt_l_i - (line->nx_i * wx + line->ny_i * wy)) / nd_i;
                if (tau.sup() < 0) continue;
                if (tau.inf() > 0) {
                    lower = std::min(lower, tau.inf());
                    continue;
                }
            }
        }
        ++counters.exact_probe_roots;
        const FT nw = line->nx * w.x() + line->ny * w.y();
        roots.push_back(Tower::of(make_one_root((line->k0 - nw) / nd, line->k1 / nd, line->L)));
    }
    for (const auto& [c, r2] : supports.circles) {
        const FT ex = w.x() - c.x(), ey = w.y() - c.y();
        if (lazy) {
            const Interval ex_i = wx - interval_of(c.x()), ey_i = wy - interval_of(c.y());
            const Interval b_i = Interval(2) * (dx_i * ex_i + dy_i * ey_i);
            if (interval_quadratic(a_i, b_i, ex_i * ex_i + ey_i * ey_i - interval_of(r2), lower)) continue;
        }
        ++counters.exact_probe_roots;
        const FT b = FT(2) * (dx * ex + dy * ey), c0 = ex * ex + ey * ey - r2;
        if (CGAL::sign(c0) == CGAL::ZERO) {
            // w on the circle: the roots are 0 (imposes nothing) and -b / a, rational.
            roots.push_back(Tower::rational(-b / a));
            continue;
        }
        quadratic_roots(a, OneRoot(b), OneRoot(c0), roots);
    }
    for (const VertexCircle* V : supports.vertex_circles) {
        if (lazy) {
            const Interval ex_i = wx - V->x_i, ey_i = wy - V->y_i;
            const Interval b_i = Interval(2) * (dx_i * ex_i + dy_i * ey_i);
            if (interval_quadratic(a_i, b_i, ex_i * ex_i + ey_i * ey_i - V->r_sq_i, lower)) continue;
        }
        ++counters.exact_probe_roots;
        const OneRoot ex = OneRoot(w.x()) - V->x, ey = OneRoot(w.y()) - V->y;
        const OneRoot b = FT(2) * (dx * ex + dy * ey), c0 = ex * ex + ey * ey - V->r_sq;
        if (sign_one_root(c0) == CGAL::ZERO) {
            // w on the vertex circle: the roots are 0 and -b / a, a one-root.
            roots.push_back(Tower::of(divide(-b, a)));
            continue;
        }
        quadratic_roots(a, b, c0, roots);
    }
    // The largest power of two (at most 1) strictly below every certified lower
    // end: a lower end of exactly 1 (a support crossed at tau = 1) must halve too.
    FT h(1);
    if (lower <= 1.0) {
        int exponent = 0;
        (void)std::frexp(lower, &exponent);  // lower in [2^(e-1), 2^e)
        const double below = std::ldexp(1.0, exponent - 2);
        if (!(below > 0.0)) {
            throw CoverageCertificateNativeError("side probe: a certified crossing parameter is below the binary64 range");
        }
        h = FT(below);
    }
    std::vector<const Tower*> positive;
    for (const Tower& t : roots) {
        if (sign_tower(t) == CGAL::POSITIVE) positive.push_back(&t);
    }
    SeparationBudget budget("side probe");
    for (;;) {
        const bool clear = std::all_of(positive.begin(), positive.end(), [&](const Tower* t) {
            return sign_tower({t->m - h, t->n, t->d}) == CGAL::POSITIVE;
        });
        if (clear) return h;
        budget.step([&] {
            double bits = 0.0;
            for (const Tower* t : positive) bits = std::max(bits, separation_bits(SeparationMeasure::of(*t), tower_radicals({t})));
            return halving_step_bound(1.0, bits);
        });
        h /= FT(2);
    }
}

Point_2 probe_point(const Point_2& w, const FT& dx, const FT& dy, const ProbeSupports& supports, const FilterMode mode,
                    DecisionCounters& counters)
{
    const FT h = probe_step(w, dx, dy, supports, mode, counters);
    return Point_2(w.x() + h * dx, w.y() + h * dy);
}

}  // namespace compas_cgal::coverage_certificate
