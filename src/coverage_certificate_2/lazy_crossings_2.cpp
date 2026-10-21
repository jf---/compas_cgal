#include "coverage_certificate_2/lazy_crossings_2.h"

#include "coverage_certificate_2/errors_2.h"

#include <CGAL/box_intersection_d.h>

#include <algorithm>
#include <iterator>
#include <utility>
#include <variant>

namespace compas_cgal::coverage_certificate {

namespace {
using Box = CGAL::Box_intersection_d::Box_with_info_d<double, 2, std::uint32_t>;
using IntersectionResult = CGAL::CK2_Intersection_traits<CK, Circle_2, Circle_2>::type;

bool overlaps(const Interval& a, const Interval& b) { return !(a.sup() < b.inf() || b.sup() < a.inf()); }

/// Order the two exact crossing points of circles A, B as (+ branch, - branch).
/// See the header: p+ - p- = 2 s (-d_y, d_x) with s > 0, so rational signs of
/// d_x, d_y plus one exact coordinate comparison decide the labels.
std::array<CircularPoint_2, 2> label_branches(const Point_2& ca, const Point_2& cb, const CircularPoint_2& p,
                                              const CircularPoint_2& q)
{
    const CK kernel;
    const FT dx = cb.x() - ca.x();
    const FT dy = cb.y() - ca.y();
    bool p_is_plus = false;
    if (CGAL::sign(dy) != CGAL::ZERO) {
        const auto c = kernel.compare_x_2_object()(p, q);
        if (c == CGAL::EQUAL) throw CoverageCertificateNativeError("two crossings with d_y != 0 share an x coordinate");
        p_is_plus = (CGAL::sign(dy) == CGAL::POSITIVE) == (c == CGAL::SMALLER);
    } else {
        const auto c = kernel.compare_y_2_object()(p, q);
        if (c == CGAL::EQUAL) throw CoverageCertificateNativeError("two crossings with d_y == 0 share a y coordinate");
        p_is_plus = (CGAL::sign(dx) == CGAL::POSITIVE) == (c == CGAL::LARGER);
    }
    return p_is_plus ? std::array<CircularPoint_2, 2>{p, q} : std::array<CircularPoint_2, 2>{q, p};
}

struct ExactCrossings {
    std::array<CircularPoint_2, 2> points;
    std::uint8_t count = 0;
    std::uint8_t multiplicity = 1;
};

ExactCrossings exact_crossings(const BoundaryCircle& a, const BoundaryCircle& b)
{
    std::vector<IntersectionResult> result;
    CK().intersect_2_object()(a.circle, b.circle, std::back_inserter(result));
    ExactCrossings out;
    std::array<CircularPoint_2, 2> raw;
    for (const auto& item : result) {
        const auto* point = std::get_if<std::pair<CircularPoint_2, unsigned>>(&item);
        if (point == nullptr) throw CoverageCertificateNativeError("non-concentric circles reported as coincident");
        if (out.count == 2) throw CoverageCertificateNativeError("two circles reported more than two crossings");
        raw[out.count++] = point->first;
        out.multiplicity = static_cast<std::uint8_t>(point->second);
    }
    if (out.count == 2) {
        if (out.multiplicity != 1) throw CoverageCertificateNativeError("two crossings reported with multiplicity 2");
        out.points = label_branches(a.circle.center(), b.circle.center(), raw[0], raw[1]);
    } else if (out.count == 1) {
        if (out.multiplicity != 2) throw CoverageCertificateNativeError("single crossing reported as transversal");
        out.points = {raw[0], raw[0]};
    }
    return out;
}
}  // namespace

CrossingTable CrossingTable::build(const CircleRegionSet& set, const FilterMode mode, DecisionCounters& counters)
{
    CrossingTable table;
    table.set_ = &set;
    table.mode_ = mode;
    const auto& circles = set.circles;
    const std::size_t n = circles.size();
    counters.boundary_circles = n;
    table.events.resize(n);
    table.candidates.resize(n);
    table.coincident.resize(n);
    table.target_circles.resize(n);

    std::vector<Box> boxes;
    boxes.reserve(n);
    for (std::uint32_t k = 0; k < n; ++k) boxes.emplace_back(circles[k].bbox, k);
    std::vector<std::pair<std::uint32_t, std::uint32_t>> candidate_pairs;
    CGAL::box_self_intersection_d(boxes.begin(), boxes.end(), [&](const Box& a, const Box& b) {
        const BoundaryCircle& A = circles[a.info()];
        const BoundaryCircle& B = circles[b.info()];
        if (A.is_target() && B.is_target()) return;  // the target-edge rule needs sweeps only
        if (!A.is_target() && !B.is_target() && A.owner == B.owner) return;
        candidate_pairs.emplace_back(std::min(a.info(), b.info()), std::max(a.info(), b.info()));
    });
    // box_self_intersection_d reports pairs in an implementation-defined order;
    // sort so that pair indices, and with them every downstream tie, are a
    // function of the input alone.
    std::sort(candidate_pairs.begin(), candidate_pairs.end());
    counters.candidate_pairs = candidate_pairs.size();
    for (const auto& [a, b] : candidate_pairs) {
        if (circles[b].is_target()) {
            table.target_circles[a].push_back(b);
        } else {
            table.candidates[a].push_back(circles[b].owner);
        }
        if (circles[a].is_target()) {
            table.target_circles[b].push_back(a);
        } else {
            table.candidates[b].push_back(circles[a].owner);
        }
    }
    for (auto* lists : {&table.candidates, &table.target_circles}) {
        for (auto& list : *lists) {
            std::sort(list.begin(), list.end());
            list.erase(std::unique(list.begin(), list.end()), list.end());
        }
    }

    const auto push_events = [&](const std::uint32_t pid, const std::uint32_t a, const std::uint32_t b,
                                 const std::uint8_t count, const std::uint8_t multiplicity) {
        const auto& pair = table.pairs_[pid];
        for (std::uint8_t branch = 0; branch < count; ++branch) {
            table.events[a].push_back({pid, branch, multiplicity, b, pair.box_x[branch], pair.box_y[branch]});
            table.events[b].push_back({pid, branch, multiplicity, a, pair.box_x[branch], pair.box_y[branch]});
        }
    };

    for (const auto& [a, b] : candidate_pairs) {
        const auto& A = circles[a];
        const auto& B = circles[b];
        const Interval dx = B.centre_x - A.centre_x;
        const Interval dy = B.centre_y - A.centre_y;
        const bool maybe_concentric = !(dx.inf() > 0 || dx.sup() < 0 || dy.inf() > 0 || dy.sup() < 0);
        if (mode == FilterMode::ExactOnly || maybe_concentric) {
            if (A.circle.center() == B.circle.center()) {
                if (A.circle.squared_radius() == B.circle.squared_radius()) {
                    ++counters.coincident_pairs;
                    table.coincident[a].push_back(b);
                    table.coincident[b].push_back(a);
                } else {
                    ++counters.pairs_exact;
                }
                continue;
            }
        }
        if (mode == FilterMode::Lazy) {
            const Interval d2 = dx * dx + dy * dy;
            const Interval q = d2 + A.radius_squared_i - B.radius_squared_i;
            const Interval disc = 4 * d2 * A.radius_squared_i - q * q;
            if (disc.sup() < 0) {
                ++counters.pairs_filtered_disjoint;
                continue;
            }
            if (disc.inf() > 0) {
                ++counters.pairs_filtered_transversal;
                const Interval l = q / (2 * d2);
                const Interval s = CGAL::sqrt(disc) / (2 * d2);
                const Interval mx = A.centre_x + l * dx;
                const Interval my = A.centre_y + l * dy;
                const auto pid = static_cast<std::uint32_t>(table.pairs_.size());
                table.pairs_.push_back({a, b, true, {mx - s * dy, mx + s * dy}, {my + s * dx, my - s * dx}, std::nullopt});
                push_events(pid, a, b, 2, 1);
                continue;
            }
        }
        ++counters.pairs_exact;
        const ExactCrossings crossings = exact_crossings(A, B);
        if (crossings.count == 0) continue;
        if (crossings.multiplicity == 2) ++counters.tangencies;
        Pair pair{a, b, false, {}, {}, crossings.points};
        for (std::uint8_t branch = 0; branch < 2; ++branch) {
            pair.box_x[branch] = Interval(CGAL::to_interval(crossings.points[branch].x()));
            pair.box_y[branch] = Interval(CGAL::to_interval(crossings.points[branch].y()));
        }
        const auto pid = static_cast<std::uint32_t>(table.pairs_.size());
        table.pairs_.push_back(std::move(pair));
        push_events(pid, a, b, crossings.count, crossings.multiplicity);
    }
    return table;
}

void CrossingTable::materialise(Pair& pair, DecisionCounters& counters)
{
    if (pair.exact) return;
    ++counters.materialised_pairs;
    const auto& circles = set_->circles;
    const ExactCrossings crossings = exact_crossings(circles[pair.a], circles[pair.b]);
    if (crossings.count != 2) {
        throw CoverageCertificateNativeError("a certified-transversal pair did not give two exact crossings");
    }
    for (std::uint8_t branch = 0; branch < 2; ++branch) {
        const Interval ex(CGAL::to_interval(crossings.points[branch].x()));
        const Interval ey(CGAL::to_interval(crossings.points[branch].y()));
        if (!overlaps(ex, pair.box_x[branch]) || !overlaps(ey, pair.box_y[branch])) {
            throw CoverageCertificateNativeError("an exact crossing lies outside its certified branch enclosure");
        }
    }
    if (overlaps(pair.box_x[0], pair.box_x[1]) && overlaps(pair.box_y[0], pair.box_y[1])) {
        ++counters.exact_branch_matches;
    }
    pair.exact = crossings.points;
}

const CircularPoint_2& CrossingTable::exact_point(const CrossingEvent& event, DecisionCounters& counters)
{
    Pair& pair = pairs_[event.pair];
    materialise(pair, counters);
    return (*pair.exact)[event.branch];
}

const CrossingTable::PairTowers& CrossingTable::towers_of(const std::uint32_t pid, DecisionCounters& counters)
{
    if (towers_.size() != pairs_.size()) towers_.resize(pairs_.size());
    if (towers_[pid]) return *towers_[pid];
    ++counters.pair_towers;
    const Pair& pair = pairs_[pid];
    if (pair.exact) {
        // The exact CK points already carry their coordinates as one-root
        // numbers a0 + a1 sqrt(root): read them, build nothing.
        auto read = std::make_unique<PairTowers>();
        for (std::size_t b = 0; b < 2; ++b) {
            read->x[b] = tower_of_root((*pair.exact)[b].x());
            read->y[b] = tower_of_root((*pair.exact)[b].y());
        }
        towers_[pid] = std::move(read);
        return *towers_[pid];
    }
    const auto& circles = set_->circles;
    const auto& ca = circles[pair.a].circle.center();
    const auto& cb = circles[pair.b].circle.center();
    const FT dx = cb.x() - ca.x();
    const FT dy = cb.y() - ca.y();
    const FT d2 = dx * dx + dy * dy;
    const FT ra2 = circles[pair.a].circle.squared_radius();
    const FT q = d2 + ra2 - circles[pair.b].circle.squared_radius();
    const FT disc = FT(4) * d2 * ra2 - q * q;
    if (CGAL::sign(disc) == CGAL::NEGATIVE) throw CoverageCertificateNativeError("a crossing pair has a negative discriminant");
    const FT two_d2 = FT(2) * d2;
    const FT mx = ca.x() + q * dx / two_d2;
    const FT my = ca.y() + q * dy / two_d2;
    auto towers = std::make_unique<PairTowers>();
    for (std::size_t b = 0; b < 2; ++b) {
        const FT sigma = b == 0 ? FT(1) : FT(-1);
        towers->x[b] = {OneRoot(mx), OneRoot(sigma * (-dy) / two_d2), OneRoot(disc)};
        towers->y[b] = {OneRoot(my), OneRoot(sigma * dx / two_d2), OneRoot(disc)};
        // Contract on the filter: the certified enclosure holds the exact value.
        if (verify_enclosures_ &&
            (!overlaps(interval_of(towers->x[b]), pair.box_x[b]) || !overlaps(interval_of(towers->y[b]), pair.box_y[b]))) {
            throw CoverageCertificateNativeError("an exact crossing lies outside its certified branch enclosure");
        }
    }
    towers_[pid] = std::move(towers);
    return *towers_[pid];
}

const Tower& CrossingTable::exact_x(const CrossingEvent& event, DecisionCounters& counters)
{
    return towers_of(event.pair, counters).x[event.branch];
}

const Tower& CrossingTable::exact_y(const CrossingEvent& event, DecisionCounters& counters)
{
    return towers_of(event.pair, counters).y[event.branch];
}

int CrossingTable::half(const CrossingEvent& event, const BoundaryCircle& circle, DecisionCounters& counters)
{
    if (mode_ == FilterMode::Lazy) {
        if (event.y.sup() < circle.centre_y.inf()) return 0;
        if (event.y.inf() > circle.centre_y.sup()) return 1;
    }
    ++counters.exact_half;
    const auto& c = circle.circle.center();
    const CGAL::Sign sy = sign_difference(exact_y(event, counters), Tower::rational(c.y()));
    if (sy == CGAL::NEGATIVE) return 0;
    if (sy == CGAL::POSITIVE) return 1;
    return sign_difference(exact_x(event, counters), Tower::rational(c.x())) == CGAL::NEGATIVE ? 0 : 1;
}

CGAL::Comparison_result CrossingTable::compare_x(const CrossingEvent& a, const CrossingEvent& b,
                                                 DecisionCounters& counters)
{
    if (a.pair == b.pair && a.branch == b.branch) return CGAL::EQUAL;
    if (mode_ == FilterMode::Lazy) {
        if (a.x.sup() < b.x.inf()) return CGAL::SMALLER;
        if (a.x.inf() > b.x.sup()) return CGAL::LARGER;
    }
    ++counters.exact_compare_x;
    return static_cast<CGAL::Comparison_result>(sign_difference(exact_x(a, counters), exact_x(b, counters)));
}

CGAL::Comparison_result CrossingTable::compare_x(const CrossingEvent& a, const Point_2& q, const Interval& qx,
                                                 DecisionCounters& counters)
{
    if (mode_ == FilterMode::Lazy) {
        if (a.x.sup() < qx.inf()) return CGAL::SMALLER;
        if (a.x.inf() > qx.sup()) return CGAL::LARGER;
    }
    ++counters.exact_compare_x;
    return static_cast<CGAL::Comparison_result>(sign_difference(exact_x(a, counters), Tower::rational(q.x())));
}

}  // namespace compas_cgal::coverage_certificate
