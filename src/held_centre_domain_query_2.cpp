#include "held_centre_domain_query_2.h"

#include <CGAL/FPU.h>
#include <CGAL/Uncertain.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

namespace compas_cgal::held_centre_domain {

namespace {

// --- number-type adapters: an interval sign that cannot be decided throws
// --- CGAL::Uncertain_conversion_exception (the static-filter idiom).

int sign_of(const FT& value) { return static_cast<int>(CGAL::sign(value)); }
int sign_of(const ExactRational& value) { return static_cast<int>(CGAL::sign(value)); }
int sign_of(const Interval& value)
{
    const CGAL::Sign certain = CGAL::sign(value);
    return static_cast<int>(certain);
}
bool zero(const FT& value) { return CGAL::is_zero(value); }
bool zero(const ExactRational& value) { return CGAL::is_zero(value); }
bool zero(const Interval& value)
{
    const bool certain = CGAL::is_zero(value);
    return certain;
}

template <typename NT>
bool is_rational(const OneRootT<NT>& value)
{
    return zero(value.b) || zero(value.r);
}

template <typename NT>
const NT& root_of(const OneRootT<NT>& x, const OneRootT<NT>& y)
{
    return is_rational(x) ? y.r : x.r;
}

template <typename NT>
int sign_one_root(const OneRootT<NT>& value)
{
    const int a = sign_of(value.a);
    const int b = zero(value.r) ? 0 : sign_of(value.b);
    if (b == 0) return a;
    if (a == 0 || a == b) return b;
    return a * sign_of(value.a * value.a - value.b * value.b * value.r);
}

/// Exact sign of a + b*sqrt(alpha) + c*sqrt(beta): sign(u) and sign(c) with
/// u = a + b*sqrt(alpha), and when they disagree, the sign of u^2 - c^2*beta,
/// which lies in Q(sqrt(alpha)). Valid for any alpha, beta >= 0.
template <typename NT>
int sign_two_roots(const std::type_identity_t<NT>& a, const std::type_identity_t<NT>& b, const NT& alpha,
                   const std::type_identity_t<NT>& c, const NT& beta)
{
    const int first = sign_one_root(OneRootT<NT>{a, b, alpha});
    const int second = zero(beta) ? 0 : sign_of(c);
    if (second == 0) return first;
    if (first == 0 || first == second) return second;
    const int squared = sign_one_root(OneRootT<NT>{a * a + b * b * alpha - c * c * beta, 2 * a * b, alpha});
    return first > 0 ? squared : -squared;
}

template <typename NT>
int compare_values(const OneRootT<NT>& first, const OneRootT<NT>& second)
{
    const NT difference = first.a - second.a;
    if (is_rational(second)) return is_rational(first) ? sign_of(difference) : sign_one_root(OneRootT<NT>{difference, first.b, first.r});
    if (is_rational(first)) return sign_one_root(OneRootT<NT>{difference, -second.b, second.r});
    return sign_two_roots(difference, first.b, first.r, -second.b, second.r);
}

/// k0 + kx*x + ky*y of a one-root point (the coordinates share one root).
template <typename NT>
OneRootT<NT> linear(const std::type_identity_t<NT>& k0, const std::type_identity_t<NT>& kx,
                    const std::type_identity_t<NT>& ky, const OneRootT<NT>& x, const OneRootT<NT>& y)
{
    return {k0 + kx * x.a + ky * y.a, kx * x.b + ky * y.b, root_of(x, y)};
}

/// |v - q|^2 for a one-root point v and a rational query q.
template <typename NT>
OneRootT<NT> squared_distance(const OneRootT<NT>& x, const OneRootT<NT>& y, const NT& qx, const NT& qy)
{
    const NT dx = x.a - qx;
    const NT dy = y.a - qy;
    const NT& root = root_of(x, y);
    return {dx * dx + x.b * x.b * root + dy * dy + y.b * y.b * root, 2 * (dx * x.b + dy * y.b), root};
}

/// Half-turn class of direction q - c from s - c in the travel sense: 0 for a
/// swept angle in (0, pi), 1 for [pi, 2pi); `at_start` for angle 0.
struct HalfTurn {
    int half;
    bool at_start;
};

template <typename NT>
HalfTurn half_turn(const PieceNumbers<NT>& n, int sense, const OneRootT<NT>& px, const OneRootT<NT>& py,
                   const NT& qx, const NT& qy)
{
    // cross(s - c, q - c) and dot(s - c, q - c) for a RATIONAL q, linear in the point s = (px, py).
    const NT ux = qx - n.cx;
    const NT uy = qy - n.cy;
    const int side = sense * sign_one_root(linear(-n.cx * uy + n.cy * ux, uy, -ux, px, py));
    if (side > 0) return {0, false};
    if (side < 0) return {1, false};
    return sign_one_root(linear(-n.cx * ux - n.cy * uy, ux, uy, px, py)) > 0 ? HalfTurn{0, true} : HalfTurn{1, false};
}

/// Sign of s1*t1 + combine * s2*t2 for s in Q(sqrt alpha), t in Q(sqrt beta).
template <typename NT>
int sign_of_products(const OneRootT<NT>& s1, const OneRootT<NT>& t1, const OneRootT<NT>& s2, const OneRootT<NT>& t2,
                     int combine, const NT& alpha, const NT& beta)
{
    // A + B sqrt(alpha) + C sqrt(beta) + D sqrt(alpha beta), grouped as u + sqrt(beta) w.
    const NT A = s1.a * t1.a + combine * s2.a * t2.a;
    const NT B = s1.b * t1.a + combine * s2.b * t2.a;
    const NT C = s1.a * t1.b + combine * s2.a * t2.b;
    const NT D = s1.b * t1.b + combine * s2.b * t2.b;
    const int u = sign_one_root(OneRootT<NT>{A, B, alpha});
    const int w = zero(beta) ? 0 : sign_one_root(OneRootT<NT>{C, D, alpha});
    if (w == 0) return u;
    if (u == 0 || u == w) return w;
    const int squared =
        sign_one_root(OneRootT<NT>{A * A + B * B * alpha - beta * (C * C + D * D * alpha), 2 * A * B - 2 * beta * C * D, alpha});
    return u > 0 ? squared : -squared;
}

template <typename NT>
const PieceNumbers<NT>& numbers(const QueryPiece& piece);
template <>
const PieceNumbers<FT>& numbers<FT>(const QueryPiece& piece) { return piece.exact; }
template <>
const PieceNumbers<Interval>& numbers<Interval>(const QueryPiece& piece) { return piece.approximate; }
template <>
const PieceNumbers<ExactRational>& numbers<ExactRational>(const QueryPiece& piece) { return piece.rational; }

/// Rigorous lower bound of the squared distance from a query enclosed by qx x qy
/// to a certified outward piece box (interval arithmetic: call under rounding
/// protection).
double box_lower_bound(const QueryPiece& piece, const Interval& qx, const Interval& qy)
{
    const Interval x = Interval(piece.xmin, piece.xmax) - qx;
    const Interval y = Interval(piece.ymin, piece.ymax) - qy;
    const Interval gx = x.inf() > 0 ? x : x.sup() < 0 ? -x : Interval(0);
    const Interval gy = y.inf() > 0 ? y : y.sup() < 0 ? -y : Interval(0);
    return (CGAL::square(gx) + CGAL::square(gy)).inf();
}

/// Whether the query projects strictly inside the piece (not onto an end).
template <typename NT>
bool strictly_inside(const QueryPiece& piece, const NT& qx, const NT& qy)
{
    const PieceNumbers<NT>& n = numbers<NT>(piece);
    if (!piece.circular) {
        const NT along = n.dx * qx + n.dy * qy;
        return sign_one_root(linear(along, -n.dx, -n.dy, n.sx, n.sy)) > 0 &&
               sign_one_root(linear(-along, n.dx, n.dy, n.tx, n.ty)) > 0;
    }
    const HalfTurn at_query = half_turn(n, piece.sense, n.sx, n.sy, qx, qy);
    if (at_query.at_start) return false;  // the radial projection is the source vertex
    if (piece.full) return true;
    if (at_query.half != piece.target_half) return at_query.half < piece.target_half;
    // sense * cross(q - c, t - c) > 0: q precedes t within the same half turn.
    const NT ux = qx - n.cx;
    const NT uy = qy - n.cy;
    return piece.sense * sign_one_root(linear(-ux * n.cy + uy * n.cx, ux, -uy, n.ty, n.tx)) > 0;
}

template <typename NT>
OneRootT<NT> interior_squared_distance(const QueryPiece& piece, const NT& qx, const NT& qy)
{
    const PieceNumbers<NT>& n = numbers<NT>(piece);
    if (!piece.circular) {
        const NT height = n.a * qx + n.b * qy + n.c;
        return {height * height / n.norm_squared, NT(0), NT(0)};
    }
    const NT ux = qx - n.cx;
    const NT uy = qy - n.cy;
    const NT radial = ux * ux + uy * uy;
    return {radial + n.radius * n.radius, -2 * n.radius, radial};
}

}  // namespace

int sign(const OneRoot& value) { return sign_one_root(value); }

int compare(const OneRoot& first, const OneRoot& second) { return compare_values(first, second); }


namespace {

/// Tight certified enclosure of an exact rational (from its exact value, so a
/// lazy construction history cannot widen it).
Interval tight(const FT& value) { return Interval(CGAL::to_interval(CGAL::exact(value))); }

/// Tight certified enclosure of a + b*sqrt(r): sqrt(r) is bracketed by two
/// doubles whose squares are compared with r exactly, and a + b*s is evaluated
/// exactly at both ends, so cancellation between a and b*sqrt(r) costs nothing.
Interval tight(const OneRoot& value)
{
    if (CGAL::is_zero(value.b) || CGAL::is_zero(value.r)) return tight(value.a);
    double low = std::sqrt(CGAL::to_double(value.r));
    double high = low;
    while (CGAL::compare(FT(low) * FT(low), value.r) == CGAL::LARGER) low = std::nextafter(low, 0.0);
    while (CGAL::compare(FT(high) * FT(high), value.r) == CGAL::SMALLER) high = std::nextafter(high, std::numeric_limits<double>::infinity());
    const Interval first = tight(FT(value.a + value.b * FT(low)));
    const Interval second = tight(FT(value.a + value.b * FT(high)));
    return Interval(std::min(first.inf(), second.inf()), std::max(first.sup(), second.sup()));
}

/// A one-root value as an enclosure (b = r = 0 in the interval stage): sound,
/// because the interval stage only ever needs certified enclosures.
OneRootT<Interval> enclosed(const OneRoot& value) { return {tight(value), Interval(0), Interval(0)}; }

}  // namespace

std::pair<double, double> enclosure(const OneRoot& value)
{
    const Interval result = tight(value);
    return {result.inf(), result.sup()};
}

void finish(QueryPiece& piece)
{
    const PieceNumbers<FT>& e = piece.exact;
    auto exact = [](const FT& value) { return ExactRational(CGAL::exact(value)); };
    auto exact_root = [&](const OneRoot& value) {
        return OneRootT<ExactRational>{exact(value.a), exact(value.b), exact(value.r)};
    };
    piece.rational = {exact(e.a),  exact(e.b),  exact(e.c),  exact(e.norm_squared), exact(e.dx),
                      exact(e.dy), exact(e.cx), exact(e.cy), exact(e.radius),       exact_root(e.sx),
                      exact_root(e.sy), exact_root(e.tx), exact_root(e.ty)};
    piece.approximate = {tight(e.a),  tight(e.b),  tight(e.c),      tight(e.norm_squared), tight(e.dx),
                         tight(e.dy), tight(e.cx), tight(e.cy),     tight(e.radius),       enclosed(e.sx),
                         enclosed(e.sy), enclosed(e.tx), enclosed(e.ty)};
    if (!piece.circular || piece.full) return;
    // cross(s - c, t - c) mixes the roots of s and t: decided over the rational parts.
    const OneRoot sx{e.sx.a - e.cx, e.sx.b, e.sx.r};
    const OneRoot sy{e.sy.a - e.cy, e.sy.b, e.sy.r};
    const OneRoot tx{e.tx.a - e.cx, e.tx.b, e.tx.r};
    const OneRoot ty{e.ty.a - e.cy, e.ty.b, e.ty.r};
    const FT& alpha = root_of(sx, sy);
    const FT& beta = root_of(tx, ty);
    const int side = piece.sense * sign_of_products(sx, ty, sy, tx, -1, alpha, beta);
    if (side != 0) {
        piece.target_half = side > 0 ? 0 : 1;
        return;
    }
    // Collinear with the start: the target is the antipode (a distinct point of the circle).
    piece.target_half = 1;
}

namespace {

/// One candidate nearest point with a certified enclosure of its squared distance.
struct Candidate {
    NearestKind kind;
    std::size_t piece;
    bool ambiguous;  // every point of the piece is nearest (the query is a circle's centre)
    Interval squared_distance;
};

/// Round-to-nearest for the exact stage, inside the interval stage's upward rounding.
using NearestRounding = CGAL::Protect_FPU_rounding<true>;

/// The query's exact coordinates, converted once per query on first need.
struct ExactQuery {
    const RationalPoint& point;
    bool ready = false;
    ExactRational x, y;
    const ExactRational& qx() { prepare(); return x; }
    const ExactRational& qy() { prepare(); return y; }
    void prepare()
    {
        if (ready) return;
        x = CGAL::exact(point.x());
        y = CGAL::exact(point.y());
        ready = true;
    }
};

OneRootT<ExactRational> exact_squared_distance(const std::vector<QueryPiece>& pieces, const Candidate& candidate,
                                               ExactQuery& query)
{
    const NearestRounding nearest(CGAL_FE_TONEAREST);
    const PieceNumbers<ExactRational>& n = pieces[candidate.piece].rational;
    if (candidate.kind == NearestKind::Vertex) return squared_distance(n.sx, n.sy, query.qx(), query.qy());
    if (candidate.ambiguous) return {n.radius * n.radius, ExactRational(0), ExactRational(0)};
    return interior_squared_distance(pieces[candidate.piece], query.qx(), query.qy());
}

/// Filtered three-way comparison of two candidates' squared distances.
int compare_candidates(const std::vector<QueryPiece>& pieces, const Candidate& first, const Candidate& second,
                       ExactQuery& query, QueryWorkspace& workspace)
{
    if (first.squared_distance.sup() < second.squared_distance.inf()) return -1;
    if (first.squared_distance.inf() > second.squared_distance.sup()) return 1;
    ++workspace.exact_fallbacks;
    const OneRootT<ExactRational> left = exact_squared_distance(pieces, first, query);
    const OneRootT<ExactRational> right = exact_squared_distance(pieces, second, query);
    const NearestRounding nearest(CGAL_FE_TONEAREST);
    return compare_values(left, right);
}

/// Filtered: does the query project strictly inside the piece?
bool filtered_strictly_inside(const QueryPiece& piece, const Interval& qx, const Interval& qy, ExactQuery& query,
                              QueryWorkspace& workspace)
{
    try {
        return strictly_inside(piece, qx, qy);
    } catch (const CGAL::Uncertain_conversion_exception&) {
        ++workspace.exact_fallbacks;
        const NearestRounding nearest(CGAL_FE_TONEAREST);
        return strictly_inside(piece, query.qx(), query.qy());
    }
}

/// Filtered: is the query exactly the centre of the piece's circle?
bool filtered_at_centre(const QueryPiece& piece, const Interval& qx, const Interval& qy, ExactQuery& query)
{
    const Interval radial = CGAL::square(qx - piece.approximate.cx) + CGAL::square(qy - piece.approximate.cy);
    if (radial.inf() > 0) return false;
    const NearestRounding nearest(CGAL_FE_TONEAREST);
    return query.qx() == piece.rational.cx && query.qy() == piece.rational.cy;
}

Interval vertex_interval(const QueryPiece& piece, const Interval& qx, const Interval& qy)
{
    return CGAL::square(piece.approximate.sx.a - qx) + CGAL::square(piece.approximate.sy.a - qy);
}

/// Enclosure of the interior distance without cancellation: (sqrt(D) - R)^2 for a circle.
Interval interior_interval(const QueryPiece& piece, const Interval& qx, const Interval& qy)
{
    const PieceNumbers<Interval>& n = piece.approximate;
    if (!piece.circular) return CGAL::square(n.a * qx + n.b * qy + n.c) / n.norm_squared;
    return CGAL::square(CGAL::sqrt(CGAL::square(qx - n.cx) + CGAL::square(qy - n.cy)) - n.radius);
}

}  // namespace

NearestDecision decide_nearest(const std::vector<QueryPiece>& pieces, const RationalPoint& query, const FT* gap_squared,
                               QueryWorkspace& workspace)
{
    std::vector<double>& lower = workspace.lower;
    std::vector<std::size_t>& seen = workspace.seen_vertices;
    lower.resize(pieces.size());
    const CGAL::Protect_FPU_rounding<true> protect;
    const Interval qx(CGAL::to_interval(query.x()));
    const Interval qy(CGAL::to_interval(query.y()));
    std::size_t closest = 0;
    for (std::size_t i = 0; i < pieces.size(); ++i) {
        lower[i] = box_lower_bound(pieces[i], qx, qy);
        if (lower[i] < lower[closest]) closest = i;
    }
    ExactQuery exact_query{query};
    Candidate best{NearestKind::Vertex, closest, false, vertex_interval(pieces[closest], qx, qy)};
    bool ambiguous = false;
    seen.assign(1, closest);
    auto offer = [&](const Candidate& candidate) {
        const int order = compare_candidates(pieces, candidate, best, exact_query, workspace);
        if (order < 0) {
            best = candidate;
            ambiguous = candidate.ambiguous;
        } else if (order == 0) {
            ambiguous = true;  // distinct candidates: an interior point is never a vertex or another piece's point
        }
    };
    auto offer_vertex = [&](std::size_t vertex) {
        if (std::find(seen.begin(), seen.end(), vertex) != seen.end()) return;
        seen.push_back(vertex);
        offer({NearestKind::Vertex, vertex, false, vertex_interval(pieces[vertex], qx, qy)});
    };
    for (std::size_t i = 0; i < pieces.size(); ++i) {
        if (lower[i] > best.squared_distance.sup()) continue;
        const QueryPiece& piece = pieces[i];
        offer_vertex(i);
        offer_vertex(piece.target_vertex);
        if (piece.circular && filtered_at_centre(piece, qx, qy, exact_query)) {
            offer({NearestKind::Interior, i, true, CGAL::square(piece.approximate.radius)});
            continue;
        }
        if (filtered_strictly_inside(piece, qx, qy, exact_query, workspace)) {
            offer({NearestKind::Interior, i, false, interior_interval(piece, qx, qy)});
        }
    }
    bool beyond = false;
    if (gap_squared != nullptr) {
        const Interval gap(CGAL::to_interval(*gap_squared));
        if (best.squared_distance.inf() > gap.sup()) {
            beyond = true;
        } else if (best.squared_distance.sup() > gap.inf()) {
            ++workspace.exact_fallbacks;
            const OneRootT<ExactRational> distance = exact_squared_distance(pieces, best, exact_query);
            const NearestRounding nearest(CGAL_FE_TONEAREST);
            beyond = compare_values(distance, OneRootT<ExactRational>{CGAL::exact(*gap_squared), ExactRational(0), ExactRational(0)}) > 0;
        }
    }
    return {ambiguous, best.kind, best.piece, beyond};
}

NearestValues nearest_values(const std::vector<QueryPiece>& pieces, const RationalPoint& query, const NearestDecision& decision)
{
    const QueryPiece& piece = pieces[decision.piece];
    const PieceNumbers<FT>& n = piece.exact;
    const FT& qx = query.x();
    const FT& qy = query.y();
    if (decision.kind == NearestKind::Vertex) return {squared_distance(n.sx, n.sy, qx, qy), n.sx, n.sy};
    if (!piece.circular) {
        const FT height = n.a * qx + n.b * qy + n.c;
        const FT scale = height / n.norm_squared;
        return {{height * scale, FT(0), FT(0)}, {qx - scale * n.a, FT(0), FT(0)}, {qy - scale * n.b, FT(0), FT(0)}};
    }
    const FT ux = qx - n.cx;
    const FT uy = qy - n.cy;
    const FT radial = ux * ux + uy * uy;
    if (CGAL::is_zero(radial)) return {{n.radius * n.radius, FT(0), FT(0)}, n.sx, n.sy};
    const FT scale = n.radius / radial;
    return {{radial + n.radius * n.radius, -2 * n.radius, radial}, {n.cx, scale * ux, radial}, {n.cy, scale * uy, radial}};
}

}  // namespace compas_cgal::held_centre_domain
