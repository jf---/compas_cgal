#include "coverage_certificate_2/witness_2.h"

#include "coverage_certificate_2/errors_2.h"
#include "coverage_certificate_2/separation_2.h"

#include <cmath>
#include <numbers>
#include <optional>
#include <string>

namespace compas_cgal::coverage_certificate {

namespace {
/// Contract bound on P7 halvings. Termination is a theorem (the uncovered face
/// is open and touches the edge at w), so this never limits a correct run; it
/// turns a logic error into a loud failure instead of a hang. 4096 halvings is
/// h = 2^-4096, far below any separation a binary64-built arrangement can have.
constexpr int UNCOVERED_HALVING_CONTRACT_BOUND = 4096;

/// Witness proposals are rounded to dyadics of PROPOSAL_BITS_STEP, 2*STEP, ...
/// up to PROPOSAL_BITS_MAX fractional bits (binary64 has 52). Proposal policy
/// only: shorter dyadics give smaller exact witnesses; acceptance is exact.
constexpr int PROPOSAL_BITS_STEP = 4;
constexpr int PROPOSAL_BITS_MAX = 52;

bool is_extreme(const FT& s)
{
    return CGAL::sign(s) == CGAL::ZERO || s == FT(1) || s == FT(-1);
}

struct Endpoint {
    const CrossingEvent* event;
    std::uint8_t half;
    double angle;  ///< proposal only, in [-pi, pi]
};

/// A double angle for proposing witnesses; never used to decide anything.
double proposal_angle(const std::pair<double, double>& p, const std::uint8_t half, const BoundaryCircle& C)
{
    const double cx = CGAL::to_double(C.circle.center().x());
    const double cy = CGAL::to_double(C.circle.center().y());
    double angle = std::atan2(p.second - cy, p.first - cx);
    if (half == 0 && angle > 0.5 * std::numbers::pi) angle -= 2 * std::numbers::pi;
    if (half == 1 && angle < -0.5 * std::numbers::pi) angle += 2 * std::numbers::pi;
    return angle;
}

class EdgeSearch {
public:
    EdgeSearch(const BoundaryCircle& circle, CircleEvents& events, DecisionCounters& counters, const Endpoint& a,
               const Endpoint& b, const bool wrap)
        : circle_(circle), events_(events), counters_(counters), a_(a), b_(b), wrap_(wrap)
    {
    }

    /// -1: Q(s) at or before a; +1: at or after b; 0: strictly inside. For a
    /// wrap edge (a = last vertex, b = first vertex, through the leftmost
    /// point) only 0 and -1 are reported.
    [[nodiscard]] int classify(const FT& s) const
    {
        const Point_2 q = stereographic_point(circle_, s);
        const Interval qx = interval_of(q.x());
        const std::uint8_t qh = half_of_stereographic(s);
        const auto pa = order_event_point(events_, *a_.event, a_.half, q, qx, qh, counters_);
        if (wrap_) {
            if (pa == CGAL::SMALLER) return 0;
            const auto pb = order_event_point(events_, *b_.event, b_.half, q, qx, qh, counters_);
            return pb == CGAL::LARGER ? 0 : -1;
        }
        if (pa != CGAL::SMALLER) return -1;
        const auto pb = order_event_point(events_, *b_.event, b_.half, q, qx, qh, counters_);
        return pb == CGAL::LARGER ? 0 : 1;
    }

    [[nodiscard]] FT find() const
    {
        const double seed = std::tan(0.5 * proposal_theta());
        std::optional<FT> found;
        if (std::isfinite(seed)) {
            // Short dyadics first: they keep the witness rationals small. Each
            // proposal is accepted only by the exact classification.
            for (int bits = PROPOSAL_BITS_STEP; bits <= PROPOSAL_BITS_MAX && !found; bits += PROPOSAL_BITS_STEP) {
                const double q = std::ldexp(std::round(std::ldexp(seed, bits)), -bits);
                if (!std::isfinite(q)) continue;
                const FT s(q);
                if (!is_extreme(s) && classify(s) == 0) found = s;
            }
        }
        if (!found) found = wrap_ ? search_wrap() : search_between(std::isfinite(seed) ? FT(seed) : FT(0));
        if (is_extreme(*found)) found = away_from_extreme(*found);
        return *found;
    }

private:
    /// Proposal angle (doubles, never decisive): the midpoint of the longest
    /// piece of (theta_a, theta_b) between consecutive extreme angles k*pi/2,
    /// i.e. the point of the edge farthest from both its ends and the extreme
    /// points, where P7's outward offset stays inside the circle's box longest.
    [[nodiscard]] double proposal_theta() const
    {
        const double lo = a_.angle;
        const double hi = b_.angle + (wrap_ ? 2 * std::numbers::pi : 0.0);
        const double quarter = 0.5 * std::numbers::pi;
        double best = 0.5 * (lo + hi);
        double best_half_width = -1.0;
        double left = lo;
        for (double k = std::floor(lo / quarter) + 1; left < hi; k += 1) {
            const double right = std::min(hi, k * quarter);
            if (right > left && 0.5 * (right - left) > best_half_width) {
                best_half_width = 0.5 * (right - left);
                best = 0.5 * (left + right);
            }
            left = std::max(left, right);
        }
        return best;
    }

    /// One exact search step, within the root-separation bound of the edge.
    void step(SeparationBudget& budget, const FT* inner = nullptr) const
    {
        ++counters_.witness_bisections;
        budget.step([&] { return search_bound(inner); });
    }

    /// The stereographic parameter s = Y / (X + rho) of an endpoint, (X, Y) = P - c,
    /// as a separation measure; nullopt at the leftmost point (s infinite).
    struct EndpointS {
        SeparationMeasure s;
        Tower x, y;
        double log2_abs_bound;  ///< log2 |s| <= this
    };
    [[nodiscard]] std::optional<EndpointS> endpoint_s(const Endpoint& e) const
    {
        const auto [x, y] = events_.exact_xy(*e.event, counters_);
        const auto& c = circle_.circle.center();
        const Tower shifted{x.m - c.x() + circle_.radius, x.n, x.d};
        if (sign_tower(shifted) == CGAL::ZERO) return std::nullopt;
        const Tower dy{y.m - c.y(), y.n, y.d};
        const SeparationMeasure den = SeparationMeasure::of(shifted);
        // |Y| <= rho and X + rho >= its separation, so |s| <= rho 2^bits.
        const double log2_abs =
            SeparationMeasure::of(circle_.radius).log2_u + separation_bits(den, tower_radicals({&x}));
        return EndpointS{over(SeparationMeasure::of(dy), den), x, y, log2_abs};
    }

    /// Root-separation bound of the whole search (walks, bisection, extreme
    /// avoidance): every phase halves or doubles toward s-values whose magnitude
    /// is at most 2^log2S, and stops within the separation of the edge's
    /// s-length (both ends finite) or of the distance to the finite end.
    /// With `inner` (a parameter strictly inside the edge), also within the
    /// separation of inner from each finite end (the extreme-avoidance phase).
    [[nodiscard]] long long search_bound(const FT* inner) const
    {
        const auto ea = endpoint_s(a_);
        const auto eb = endpoint_s(b_);
        double log2_s = 1.0;
        for (const auto* e : {&ea, &eb}) {
            if (*e) log2_s = std::max(log2_s, (*e)->log2_abs_bound);
        }
        double gap_bits = log2_s;
        if (ea && eb) {
            gap_bits = separation_bits(plus(ea->s, eb->s), tower_radicals({&ea->x, &ea->y, &eb->x, &eb->y}));
        }
        if (inner != nullptr) {
            for (const auto* e : {&ea, &eb}) {
                if (!*e) continue;
                const SeparationMeasure d = plus((*e)->s, SeparationMeasure::of(*inner));
                gap_bits = std::max(gap_bits, separation_bits(d, tower_radicals({&(*e)->x, &(*e)->y})));
            }
        }
        // Two doubling walks and two halving phases, each within log2S + gap bits.
        return halving_step_bound(std::ldexp(1.0, 0), 4.0 * (log2_s + 2.0) + 2.0 * gap_bits);
    }

    [[nodiscard]] FT search_wrap() const
    {
        // Q(s) tends to the leftmost point from the upper half as s grows, so it
        // passes the last vertex: the doubling terminates.
        SeparationBudget budget("P3 circle witness search");
        for (FT s(2);; s *= FT(2)) {
            step(budget);
            if (classify(s) == 0) return s;
        }
    }

    [[nodiscard]] FT search_between(const FT& start) const
    {
        SeparationBudget budget("P3 circle witness search");
        const int c = classify(start);
        if (c == 0) return start;
        FT lo = start;
        FT hi = start;
        FT width(1);
        if (c < 0) {
            // Walk up until past a (b is a later vertex, finite in s).
            for (;; width *= FT(2)) {
                step(budget);
                hi = lo + width;
                const int ch = classify(hi);
                if (ch == 0) return hi;
                if (ch > 0) break;
                lo = hi;
            }
        } else {
            // Walk down until before b (when a is the leftmost point every
            // finite s is after it, and the walk ends inside the edge).
            for (;; width *= FT(2)) {
                step(budget);
                lo = hi - width;
                const int cl = classify(lo);
                if (cl == 0) return lo;
                if (cl < 0) break;
                hi = lo;
            }
        }
        // classify(lo) == -1 and classify(hi) == +1: bisect. Once hi - lo is
        // below twice the edge's length in s, the midpoint is strictly inside.
        for (;;) {
            step(budget);
            const FT mid = (lo + hi) / FT(2);
            const int m = classify(mid);
            if (m == 0) return mid;
            if (m < 0) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
    }

    [[nodiscard]] FT away_from_extreme(const FT& s) const
    {
        SeparationBudget budget("P3 circle witness search");
        FT h(1);
        for (;;) {
            step(budget, &s);
            h /= FT(2);
            for (const FT candidate : {FT(s + h), FT(s - h)}) {
                if (!is_extreme(candidate) && classify(candidate) == 0) return candidate;
            }
        }
    }

    const BoundaryCircle& circle_;
    CircleEvents& events_;
    DecisionCounters& counters_;
    Endpoint a_;
    Endpoint b_;
    bool wrap_;
};

/// Exact side of the rational point q relative to a rectangle (no filter).
RectSide rect_side_exact(const RectRegion& rect, const LineTable& lines, const Point_2& q)
{
    return rect_membership(rect, lines, ExactPoint::rational(q));
}

/// Filtered: true when q lies outside the closed rectangle.
bool outside_rect(const RectRegion& rect, const LineTable& lines, const Point_2& q, const Interval& qx, const Interval& qy,
                  const FilterMode mode, DecisionCounters& counters)
{
    for (std::size_t s = 0; s < 4; ++s) {
        const CGAL::Sign side = line_side(lines[rect.lines[s]], q, qx, qy, mode, counters);
        if (side != CGAL::ZERO && side != rect.inward[s]) return true;
    }
    return false;
}
}  // namespace

FT edge_parameter(const std::uint32_t k, const std::size_t edge, const CircleSweep& sweep, CircleWitnessContext& context,
                  DecisionCounters& counters)
{
    const BoundaryCircle& C = context.set.circles[k];
    const std::size_t vertices = sweep.vertex_count();
    if (vertices == 0) return FT(1, 2);
    const std::size_t ia = sweep.group_starts[edge];
    const std::size_t ib = sweep.group_starts[(edge + 1) % vertices];
    const Endpoint a{&sweep.events[ia], sweep.halves[ia],
                     proposal_angle(context.events.proposal_xy(sweep.events[ia], counters), sweep.halves[ia], C)};
    const Endpoint b{&sweep.events[ib], sweep.halves[ib],
                     proposal_angle(context.events.proposal_xy(sweep.events[ib], counters), sweep.halves[ib], C)};
    return EdgeSearch(C, context.events, counters, a, b, edge + 1 == vertices).find();
}

void audit_covered_edge(const std::uint32_t k, const std::size_t edge, const CircleSweep& sweep,
                        CircleWitnessContext& context, DecisionCounters& counters)
{
    const FT s = edge_parameter(k, edge, sweep, context, counters);
    const Point_2 w = stereographic_point(context.set.circles[k], s);
    const Interval wx = interval_of(w.x()), wy = interval_of(w.y());
    bool covered = false;
    for (const std::uint32_t j : sweep.counted) {
        const Side side = membership(context.set.regions[j], w, wx, wy, FilterMode::ExactOnly, counters);
        if (side == Side::Boundary) {
            throw CoverageCertificateNativeError("audit witness on circle " + std::to_string(k) + " lies on region " +
                                                 std::to_string(j) + ", which the crossings did not predict");
        }
        covered = covered || side == Side::Interior;
    }
    for (const std::uint32_t r : sweep.rects) {
        const RectSide side = rect_side_exact(context.rects[r], context.lines, w);
        if (side == RectSide::Boundary) {
            throw CoverageCertificateNativeError("audit witness on circle " + std::to_string(k) + " lies on rectangle " +
                                                 std::to_string(r) + ", which the crossings did not predict");
        }
        covered = covered || side == RectSide::Interior;
    }
    if (!covered) {
        throw CoverageCertificateNativeError("audit: edge " + std::to_string(edge) + " of circle " + std::to_string(k) +
                                             " is counted covered but its exact witness lies in no candidate region");
    }
}

std::pair<FT, Point_2> confirm_exposed_edge(const std::uint32_t k, const std::size_t edge, const CircleSweep& sweep,
                                            CircleWitnessContext& context, DecisionCounters& counters)
{
    const FT s = edge_parameter(k, edge, sweep, context, counters);
    const Point_2 w = stereographic_point(context.set.circles[k], s);
    confirm_circle_witness(k, w, sweep, context, counters);
    return {s, w};
}

void confirm_circle_witness(const std::uint32_t k, const Point_2& w, const CircleSweep& sweep,
                            CircleWitnessContext& context, DecisionCounters& counters)
{
    const FilterMode mode = context.events.mode();
    // P4: the decision.
    const Interval wx = interval_of(w.x()), wy = interval_of(w.y());
    for (const std::uint32_t j : sweep.counted) {
        const Side side = membership(context.set.regions[j], w, wx, wy, mode, counters);
        if (side != Side::Exterior) {
            throw CoverageCertificateNativeError("exposed-edge witness on circle " + std::to_string(k) +
                                                 (side == Side::Interior ? " lies inside" : " lies on") + " region " +
                                                 std::to_string(j) + ": the cover counters disagree with P4");
        }
    }
    for (const std::uint32_t r : sweep.rects) {
        if (!outside_rect(context.rects[r], context.lines, w, wx, wy, mode, counters)) {
            throw CoverageCertificateNativeError("exposed-edge witness on circle " + std::to_string(k) +
                                                 " lies in rectangle " + std::to_string(r) +
                                                 ": the cover counters disagree with P4");
        }
    }
}

Point_2 uncovered_beside(const std::uint32_t k, const Point_2& w, const CircleSweep& sweep, CircleWitnessContext& context,
                         const std::function<bool(const Point_2&)>& accept, DecisionCounters& counters)
{
    const BoundaryCircle& C = context.set.circles[k];
    const FilterMode mode = context.events.mode();
    const auto& c = C.circle.center();
    const FT one(1);
    FT h(1);
    for (int i = 0; i < UNCOVERED_HALVING_CONTRACT_BOUND; ++i) {
        ++counters.uncovered_halvings;
        h /= FT(2);
        const FT f = C.check_inside() ? FT(one - h) : FT(one + h);
        const Point_2 u(c.x() + f * (w.x() - c.x()), c.y() + f * (w.y() - c.y()));
        if (!C.check_inside()) {
            const FT ax = CGAL::abs(u.x() - c.x());
            const FT ay = CGAL::abs(u.y() - c.y());
            if (CGAL::compare(ax, C.radius) == CGAL::LARGER || CGAL::compare(ay, C.radius) == CGAL::LARGER) continue;
        }
        const Interval ux = interval_of(u.x()), uy = interval_of(u.y());
        bool clear = C.is_target() || membership(context.set.regions[C.owner], u, ux, uy, mode, counters) == Side::Exterior;
        for (std::size_t t = 0; clear && t < sweep.counted.size(); ++t) {
            clear = membership(context.set.regions[sweep.counted[t]], u, ux, uy, mode, counters) == Side::Exterior;
        }
        for (std::size_t t = 0; clear && t < sweep.excluded.size(); ++t) {
            clear = membership(context.set.regions[sweep.excluded[t]], u, ux, uy, mode, counters) == Side::Exterior;
        }
        for (std::size_t t = 0; clear && t < sweep.rects.size(); ++t) {
            clear = outside_rect(context.rects[sweep.rects[t]], context.lines, u, ux, uy, mode, counters);
        }
        if (clear && accept(u)) return u;
    }
    throw CoverageCertificateNativeError("P7 uncovered-point search exceeded its contract bound on circle " +
                                         std::to_string(k));
}

}  // namespace compas_cgal::coverage_certificate
