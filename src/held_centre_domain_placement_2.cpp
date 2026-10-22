#include "held_centre_domain_placement_2.h"

#include "held_centre_domain_2.h"

#include <CGAL/Fraction_traits.h>
#include <CGAL/Uncertain.h>

#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <tuple>
#include <utility>

namespace compas_cgal::held_centre_domain {

namespace {

using Root = PlacementRoot;
using Box = PlacementInterval;

constexpr double INFINITY_BINARY64 = std::numeric_limits<double>::infinity();

Root root_of(const ExactRational& a0, const ExactRational& a1, const ExactRational& root)
{
    if (CGAL::is_zero(a1) || CGAL::is_zero(root)) return Root(a0);
    return Root(a0, a1, root);
}

Root root_of(const OneRootT<ExactRational>& value) { return root_of(value.a, value.b, value.r); }

/// The value with a tight certified enclosure (held_centre_domain_query_2.h
/// `enclosure`: the root bracketed by doubles whose squares are compared
/// exactly, so cancellation between a0 and a1 sqrt(r) does not widen it).
EnclosedRoot enclosed(const Root& value)
{
    if (!value.is_extended()) return {value, Box(CGAL::to_interval(value.a0()))};
    return {value, Box(enclosure(OneRoot{FT(value.a0()), FT(value.a1()), FT(value.root())}))};
}

Box box_of(const ExactRational& value) { return Box(CGAL::to_interval(value)); }

ExactRational injected(double value)
{
    if (!std::isfinite(value)) throw InvalidCentreDomainInputError("placement inputs must be finite binary64");
    return ExactRational(value);  // exact: every finite binary64 is a dyadic rational
}

bool even(double value) { return (std::bit_cast<std::uint64_t>(value) & 1U) == 0U; }

double positive_zero(double value) { return value == 0.0 ? 0.0 : value; }

using Integer = CGAL::Fraction_traits<ExactRational>::Numerator_type;

/// m * 2^e: a binary64, or the midpoint of two neighbouring ones, as an exact
/// integer pair (no rational normalisation).
struct Dyadic {
    Integer m;
    int e;
};

/// x = mantissa * 2^exponent with an integer mantissa of at most 53 bits (exact).
Dyadic dyadic_of(double x)
{
    if (x == 0.0) return {Integer(0), 0};
    int exponent = 0;
    const double fraction = std::frexp(x, &exponent);  // x = fraction * 2^exponent, |fraction| in [0.5, 1)
    constexpr int MANTISSA_BITS = 53;
    return {Integer(static_cast<std::int64_t>(std::ldexp(fraction, MANTISSA_BITS))), exponent - MANTISSA_BITS};
}

/// (a + b) / 2 of two binary64 values, exactly.
Dyadic midpoint_of(double a, double b)
{
    const Dyadic x = dyadic_of(a), y = dyadic_of(b);
    const int low = std::min(x.e, y.e);
    return {Integer((x.m << (x.e - low)) + (y.m << (y.e - low))), low - 1};
}

/// n / d with d > 0, NOT reduced: the per-placement arithmetic of the circle
/// sites multiplies a few of these and compares by cross-multiplication, so no
/// gcd is ever taken (the normalisation of the exact rational type dominates
/// otherwise). Built from exact rationals (CGAL's Fraction_traits) and dyadics.
struct IntegerRatio {
    Integer n;
    Integer d;
};

IntegerRatio ratio_of(const ExactRational& value)
{
    IntegerRatio out;
    CGAL::Fraction_traits<ExactRational>::Decompose()(value, out.n, out.d);
    return out;
}

/// a - q for a rational a and a dyadic q.
IntegerRatio minus(const IntegerRatio& a, const Dyadic& q)
{
    if (q.e >= 0) return {Integer(a.n - ((q.m * a.d) << q.e)), a.d};
    return {Integer((a.n << -q.e) - q.m * a.d), Integer(a.d << -q.e)};
}

IntegerRatio negated(const IntegerRatio& a) { return {Integer(-a.n), a.d}; }

IntegerRatio times(const IntegerRatio& a, const IntegerRatio& b) { return {Integer(a.n * b.n), Integer(a.d * b.d)}; }

IntegerRatio plus(const IntegerRatio& a, const IntegerRatio& b) { return {Integer(a.n * b.d + b.n * a.d), Integer(a.d * b.d)}; }

/// a <=> b, cross-multiplied (both denominators positive).
CGAL::Comparison_result compare_ratios(const IntegerRatio& a, const IntegerRatio& b)
{
    return CGAL::compare(Integer(a.n * b.d), Integer(b.n * a.d));
}

int sign_of_ratio(const IntegerRatio& a) { return static_cast<int>(CGAL::sign(a.n)); }

std::pair<Integer, Integer> parts_of(const ExactRational& value)
{
    const IntegerRatio ratio = ratio_of(value);
    return {ratio.n, ratio.d};
}

/// Exact comparison of a fixed value a + b sqrt(r) (b = 0 for a rational)
/// with dyadic rationals m: sign((a - m) + b sqrt(r)) is sign(a - m) when b = 0,
/// sign(b) when a - m is zero or of b's sign, else sign(a - m) times the sign of
/// (a - m)^2 - b^2 r -- the decomposition Sqrt_extension::compare makes. The
/// search compares one value with several midpoints, so a = an / ad and
/// b^2 r = bn / bd are kept as integer pairs and every comparison is
/// cross-multiplied in integers: no gcd, no division.
class ExactVersus {
public:
    ExactVersus(const ExactRational& a, int b_sign, const ExactRational& b_squared_r, double estimate)
        : b_sign_(b_sign), estimate_(estimate)
    {
        std::tie(an_, ad_) = parts_of(a);
        std::tie(bn_, bd_) = parts_of(b_squared_r);
    }

    ExactVersus(const IntegerRatio& a, int b_sign, const IntegerRatio& b_squared_r, double estimate)
        : an_(a.n), ad_(a.d), bn_(b_squared_r.n), bd_(b_squared_r.d), b_sign_(b_sign), estimate_(estimate)
    {
    }

    [[nodiscard]] CGAL::Comparison_result versus(const Dyadic& m) const
    {
        // a - m = U / den with den > 0.
        Integer U, den;
        if (m.e >= 0) {
            U = an_ - ((m.m * ad_) << m.e);
            den = ad_;
        } else {
            U = (an_ << -m.e) - m.m * ad_;
            den = ad_ << -m.e;
        }
        const int u_sign = static_cast<int>(CGAL::sign(U));
        if (b_sign_ == 0) return static_cast<CGAL::Comparison_result>(u_sign);
        if (u_sign == 0 || u_sign == b_sign_) return static_cast<CGAL::Comparison_result>(b_sign_);
        const int squares = static_cast<int>(CGAL::compare(Integer(U * U * bd_), Integer(bn_ * den * den)));
        return static_cast<CGAL::Comparison_result>(u_sign * squares);
    }

    [[nodiscard]] double estimate() const { return estimate_; }

private:
    Integer an_, ad_, bn_, bd_;
    int b_sign_;
    double estimate_;
};

ExactVersus versus_of(const Root& value, std::optional<double> estimate)
{
    const double seed = estimate ? *estimate : [&] {
        const std::pair<double, double> box = value.to_interval();
        return box.first + (box.second - box.first) / 2;
    }();
    if (!value.is_extended()) return ExactVersus(value.a0(), 0, ExactRational(0), seed);
    return ExactVersus(value.a0(), static_cast<int>(CGAL::sign(value.a1())), ExactRational(value.a1() * value.a1() * value.root()), seed);
}

/// The binary64 nearest an exact value, ties to even. Starting from the
/// estimate, a candidate d is accepted when the value lies between the
/// midpoints of d with its neighbours; otherwise the neighbour on the value's
/// side is tried. The neighbours only enumerate candidates; every acceptance
/// is an exact comparison.
double nearest_by(const ExactVersus& value)
{
    double candidate = value.estimate();
    for (;;) {
        const double up = std::nextafter(candidate, INFINITY_BINARY64);
        const double down = std::nextafter(candidate, -INFINITY_BINARY64);
        if (!std::isfinite(up) || !std::isfinite(down)) {
            throw InvalidCentreDomainInputError("an exact placement value is outside the binary64 range");
        }
        const CGAL::Comparison_result above = value.versus(midpoint_of(candidate, up));
        if (above == CGAL::LARGER) {
            candidate = up;
            continue;
        }
        if (above == CGAL::EQUAL) return positive_zero(even(candidate) ? candidate : up);
        const CGAL::Comparison_result below = value.versus(midpoint_of(down, candidate));
        if (below == CGAL::SMALLER) {
            candidate = down;
            continue;
        }
        if (below == CGAL::EQUAL) return positive_zero(even(candidate) ? candidate : down);
        return positive_zero(candidate);
    }
}

double round_down_by(const ExactVersus& value)
{
    const double nearest = nearest_by(value);
    return value.versus(dyadic_of(nearest)) == CGAL::SMALLER ? std::nextafter(nearest, -INFINITY_BINARY64) : nearest;
}

double nearest_of(const Root& value, std::optional<double> estimate = std::nullopt) { return nearest_by(versus_of(value, estimate)); }

double round_down_of(const Root& value, std::optional<double> estimate = std::nullopt) { return round_down_by(versus_of(value, estimate)); }

int sign_of(const Root& value) { return static_cast<int>(CGAL::sign(value)); }

/// Filtered exact sign: the enclosure's sign when certain, else the exact one.
template <typename Exact>
int filtered_sign(const Box& box, const Exact& exact)
{
    const CGAL::Uncertain<CGAL::Sign> certain = CGAL::sign(box);
    if (CGAL::is_certain(certain)) return static_cast<int>(CGAL::get_certain(certain));
    return exact();
}

/// Slope key y / x of a normal direction; `nullopt` for x = 0.
std::optional<ExactRational> slope_key(const ExactRational& x, const ExactRational& y)
{
    if (CGAL::is_zero(x)) return std::nullopt;
    return ExactRational(y / x);
}

}  // namespace

double nearest_binary64(const ExactRational& value)
{
    return nearest_by(ExactVersus(value, 0, ExactRational(0), CGAL::to_double(value)));
}

double round_down_binary64(const ExactRational& value)
{
    return round_down_by(ExactVersus(value, 0, ExactRational(0), CGAL::to_double(value)));
}

std::vector<PlacementFrame> placement_frames(const std::vector<QueryPiece>& pieces,
                                             const std::vector<Kernel::Vector_2>& declared_normals)
{
    std::map<ExactRational, std::pair<ExactRational, ExactRational>> by_slope;
    std::optional<std::pair<ExactRational, ExactRational>> vertical;
    for (const Kernel::Vector_2& normal : declared_normals) {
        const ExactRational nx(CGAL::exact(normal.x()));
        const ExactRational ny(CGAL::exact(normal.y()));
        const std::optional<ExactRational> key = slope_key(nx, ny);
        if (key) {
            by_slope.emplace(*key, std::make_pair(nx, ny));
        } else {
            vertical = std::make_pair(nx, ny);
        }
    }
    std::vector<PlacementFrame> frames(pieces.size());
    for (std::size_t index = 0; index < pieces.size(); ++index) {
        const QueryPiece& piece = pieces[index];
        const PieceNumbers<ExactRational>& n = piece.rational;
        PlacementFrame& frame = frames[index];
        frame.circular = piece.circular;
        frame.vertex_x = root_of(n.sx);
        frame.vertex_y = root_of(n.sy);
        if (piece.circular) {
            frame.source_x = enclosed(frame.vertex_x - n.cx);
            frame.source_y = enclosed(frame.vertex_y - n.cy);
            frame.target_x = enclosed(root_of(n.tx) - n.cx);
            frame.target_y = enclosed(root_of(n.ty) - n.cy);
            frame.centre_x_box = box_of(n.cx);
            frame.centre_y_box = box_of(n.cy);
            frame.centre_x_view = CGAL::to_double(n.cx);
            frame.centre_y_view = CGAL::to_double(n.cy);
            frame.radius_view = CGAL::to_double(n.radius);
            std::tie(frame.centre_x_n, frame.centre_x_d) = parts_of(n.cx);
            std::tie(frame.centre_y_n, frame.centre_y_d) = parts_of(n.cy);
            std::tie(frame.radius_n, frame.radius_d) = parts_of(n.radius);
            continue;
        }
        // (a, b) is a normal of the piece's line.
        const std::optional<ExactRational> key = slope_key(n.a, n.b);
        const std::pair<ExactRational, ExactRational>* found = nullptr;
        if (key) {
            const auto entry = by_slope.find(*key);
            if (entry != by_slope.end()) found = &entry->second;
        } else if (vertical) {
            found = &*vertical;
        }
        if (found == nullptr) continue;
        // Inward: to the left of the travel direction (dx, dy), cross(d, n) > 0.
        const bool left = CGAL::sign(ExactRational(n.dx * found->second - n.dy * found->first)) == CGAL::POSITIVE;
        frame.has_normal = true;
        frame.nx = left ? found->first : ExactRational(-found->first);
        frame.ny = left ? found->second : ExactRational(-found->second);
        // a x + b y + c = 0 with (a, b) = lambda n: n.x = -c / lambda on the support.
        frame.level = ExactRational(-n.c / ExactRational(n.a * frame.nx + n.b * frame.ny));
        // Travel unit u = (n_y, -n_x): n is its left normal.
        frame.ux_box = box_of(frame.ny);
        frame.uy_box = box_of(ExactRational(-frame.nx));
        frame.source_along = enclosed(frame.vertex_x * frame.ny - frame.vertex_y * frame.nx);
        frame.target_along = enclosed(root_of(n.tx) * frame.ny - root_of(n.ty) * frame.nx);
    }
    return frames;
}

PlacedCircle place_on_normal(const QueryPiece& piece, const PlacementFrame& frame, double qx, double qy, bool vertex,
                             double guide_radius)
{
    if (!std::isfinite(guide_radius) || !(guide_radius > 0.0)) {
        throw InvalidCentreDomainInputError("a guide radius must be a positive finite binary64");
    }
    const PieceNumbers<ExactRational>& n = piece.rational;
    const ExactRational px = injected(qx);
    const ExactRational py = injected(qy);
    PlacedCircle out;
    double distance = guide_radius;
    if (piece.circular && piece.sense == 1 && CGAL::compare(n.radius, ExactRational(guide_radius)) == CGAL::SMALLER) {
        // Inside a convex arc of radius rho' the tangent circles q + t*n stay
        // inside its circle only for t <= rho' (t = rho' is the circle itself).
        distance = round_down_binary64(n.radius);
        out.limited = true;
    }
    const ExactRational step(distance);

    if (!piece.circular) {
        if (!frame.has_normal) throw NoRationalUnitNormalError("a line piece of the centre domain has no rational unit normal");
        out.normal_x = nearest_binary64(frame.nx);
        out.normal_y = nearest_binary64(frame.ny);
        if (!vertex) {
            // Position u.q along the travel direction against u.s and u.t.
            const Box along_box = frame.ux_box * Box(qx) + frame.uy_box * Box(qy);
            std::optional<ExactRational> along;
            auto exact_along = [&]() -> const ExactRational& {
                if (!along) along = ExactRational(frame.ny * px - frame.nx * py);
                return *along;
            };
            const int versus_source = filtered_sign(along_box - frame.source_along.box,
                                                    [&] { return -static_cast<int>(frame.source_along.value.compare(exact_along())); });
            if (versus_source != 0) {
                const int versus_target = filtered_sign(along_box - frame.target_along.box,
                                                        [&] { return -static_cast<int>(frame.target_along.value.compare(exact_along())); });
                if (versus_source < 0 || versus_target >= 0) {
                    out.status = PlacementStatus::OutsidePiece;
                    return out;
                }
                // q + t n = p + (t - h) n with h = n.p - level.
                const ExactRational offset(step - (frame.nx * px + frame.ny * py - frame.level));
                out.placement = SitePlacement::Interior;
                out.cx = nearest_binary64(ExactRational(px + offset * frame.nx));
                out.cy = nearest_binary64(ExactRational(py + offset * frame.ny));
            } else {
                vertex = true;
            }
        }
        if (vertex) {
            out.placement = SitePlacement::Source;
            out.cx = nearest_of(frame.vertex_x + ExactRational(step * frame.nx));
            out.cy = nearest_of(frame.vertex_y + ExactRational(step * frame.ny));
        }
        const ExactRational clearance(frame.nx * ExactRational(out.cx) + frame.ny * ExactRational(out.cy) - frame.level);
        out.radius = CGAL::compare(step, clearance) != CGAL::LARGER ? distance : round_down_binary64(clearance);
        return out;
    }

    const int sense = piece.sense;
    const ExactRational exact_sense(sense);
    const IntegerRatio centre_x{frame.centre_x_n, frame.centre_x_d};
    const IntegerRatio centre_y{frame.centre_y_n, frame.centre_y_d};
    const IntegerRatio radius{frame.radius_n, frame.radius_d};
    if (!vertex) {
        // w = q - C as unreduced integer ratios; the exact rationals only on the
        // exact fallbacks of the filtered signs.
        const IntegerRatio wx_ratio = negated(minus(centre_x, dyadic_of(qx)));
        const IntegerRatio wy_ratio = negated(minus(centre_y, dyadic_of(qy)));
        if (sign_of_ratio(wx_ratio) == 0 && sign_of_ratio(wy_ratio) == 0) {
            out.status = PlacementStatus::ProbeAtCircleCentre;
            return out;
        }
        const Box wx_box = Box(qx) - frame.centre_x_box;
        const Box wy_box = Box(qy) - frame.centre_y_box;
        // Half turn of w from the source vector s - C: by the sign of the cross
        // product, and at a collinear w by the sign of the dot product.
        const EnclosedRoot& ux = frame.source_x;
        const EnclosedRoot& uy = frame.source_y;
        const ExactRational wx(px - n.cx);
        const ExactRational wy(py - n.cy);
        const int cross = filtered_sign(ux.box * wy_box - uy.box * wx_box, [&] { return sign_of(ux.value * wy - uy.value * wx); });
        const int side = sense * cross;
        int half = side > 0 ? 0 : 1;
        if (side == 0) {
            half = filtered_sign(ux.box * wx_box + uy.box * wy_box, [&] { return sign_of(ux.value * wx + uy.value * wy); }) > 0 ? 0 : 1;
        }
        if (half == 0 && cross == 0) {
            vertex = true;
        } else {
            bool before = piece.full || half < piece.target_half;
            if (!piece.full && half == piece.target_half) {
                const EnclosedRoot& vx = frame.target_x;
                const EnclosedRoot& vy = frame.target_y;
                before = sense * -filtered_sign(vx.box * wy_box - vy.box * wx_box, [&] { return sign_of(vx.value * wy - vy.value * wx); }) > 0;
            }
            if (!before) {
                out.status = PlacementStatus::OutsidePiece;
                return out;
            }
            // C + P w / sqrt(D), P = rho' - sigma t, D = |w|^2: per coordinate
            // a = C_x, b = P w_x / D, so b^2 D = P^2 w_x^2 / D (one root, sqrt(D)).
            const IntegerRatio radial = plus(times(wx_ratio, wx_ratio), times(wy_ratio, wy_ratio));
            out.placement = SitePlacement::Interior;
            out.half = half;
            // Binary64 views: the normal view as the reference reports it, and
            // the centre estimate that seeds the exact nearest search.
            const IntegerRatio zero{Integer(0), Integer(1)};
            const double wx_view = nearest_by(ExactVersus(wx_ratio, 0, zero, qx - frame.centre_x_view));
            const double wy_view = nearest_by(ExactVersus(wy_ratio, 0, zero, qy - frame.centre_y_view));
            const double length = std::sqrt(nearest_by(ExactVersus(radial, 0, zero, wx_view * wx_view + wy_view * wy_view)));
            out.normal_x = (-sense * wx_view) / length;
            out.normal_y = (-sense * wy_view) / length;
            const double reach = frame.radius_view - sense * distance;
            const Dyadic sigma_t = dyadic_of(sense * distance);
            const IntegerRatio along = minus(radius, sigma_t);
            const IntegerRatio along_squared_over_d{Integer(along.n * along.n * radial.d), Integer(along.d * along.d * radial.n)};
            const int along_sign = sign_of_ratio(along);
            out.cx = nearest_by(ExactVersus(centre_x, along_sign * sign_of_ratio(wx_ratio), times(along_squared_over_d, times(wx_ratio, wx_ratio)),
                                            frame.centre_x_view + reach * wx_view / length));
            out.cy = nearest_by(ExactVersus(centre_y, along_sign * sign_of_ratio(wy_ratio), times(along_squared_over_d, times(wy_ratio, wy_ratio)),
                                            frame.centre_y_view + reach * wy_view / length));
        }
    }
    if (vertex) {
        // q + t * (-sigma / rho') (q - C) = (1 - k) q + k C with k = sigma t / rho'.
        const ExactRational k(exact_sense * step / n.radius);
        const ExactRational keep(1 - k);
        out.placement = SitePlacement::Source;
        out.half = 0;
        out.cx = nearest_of(frame.vertex_x * keep + ExactRational(k * n.cx));
        out.cy = nearest_of(frame.vertex_y * keep + ExactRational(k * n.cy));
        const ExactRational inward(ExactRational(-sense) / n.radius);
        out.normal_x = nearest_of(frame.source_x.value * inward);
        out.normal_y = nearest_of(frame.source_y.value * inward);
    }
    // Clearance of the rounded centre from the support circle: sqrt(D) - rho'
    // (Omega* outside, sense -1) or rho' - sqrt(D) (inside, sense +1).
    const IntegerRatio dx = negated(minus(centre_x, dyadic_of(out.cx)));
    const IntegerRatio dy = negated(minus(centre_y, dyadic_of(out.cy)));
    const IntegerRatio squared = plus(times(dx, dx), times(dy, dy));
    // t <= sqrt(D) - rho'  <=>  t + rho' <= 0 or (t + rho')^2 <= D;
    // t <= rho' - sqrt(D)  <=>  rho' - t >= 0 and D <= (rho' - t)^2.
    const IntegerRatio reach = sense < 0 ? negated(minus(negated(radius), dyadic_of(distance))) : minus(radius, dyadic_of(distance));
    const bool fits = sense < 0 ? (sign_of_ratio(reach) <= 0 || compare_ratios(times(reach, reach), squared) != CGAL::LARGER)
                                : (sign_of_ratio(reach) >= 0 && compare_ratios(squared, times(reach, reach)) != CGAL::LARGER);
    if (fits) {
        out.radius = distance;
        return out;
    }
    // Estimate without cancellation: sqrt(D) - rho' = (D - rho'^2) / (sqrt(D) + rho'),
    // the numerator rounded from its exact value (the estimate only seeds the
    // search; a few ulps off costs a few exact comparisons, a cancelled one
    // costs one per ulp of rho').
    const double root_view = std::hypot(out.cx - frame.centre_x_view, out.cy - frame.centre_y_view);
    const ExactRational gap(ExactRational(ExactRational(squared.n) / ExactRational(squared.d)) - n.radius * n.radius);
    const double bound_view = CGAL::to_double(gap) / (root_view + frame.radius_view);
    const IntegerRatio offset = sense < 0 ? negated(radius) : radius;
    out.radius = round_down_by(ExactVersus(offset, -sense, squared, sense < 0 ? bound_view : -bound_view));
    return out;
}

}  // namespace compas_cgal::held_centre_domain
