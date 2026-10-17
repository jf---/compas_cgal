#include "held_boundary_clearance_2.h"

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/Interval_nt.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <optional>
#include <utility>

namespace {
using namespace held_boundary_clearance;
using Interval = CGAL::Interval_nt<>;
using ApproximateKernel = CGAL::Simple_cartesian<CGAL::Interval_nt_advanced>;
// Epick's own exact ring type for the end of its filter chain. Both predicates
// below are division-free, so a ring type decides them exactly.
using ExactKernel = CGAL::Simple_cartesian<CGAL::internal::Exact_ring_selector<double>::Type>;
using ToApproximate = CGAL::Cartesian_converter<Kernel, ApproximateKernel>;
using ToExact = CGAL::Cartesian_converter<Kernel, ExactKernel>;

// Exact squared distance from a station to one site as numerator/denominator,
// ring operations only. A vertex site (interior == false) is |p - a|² / 1; an
// interior foot of segment ab is det(b - a, p - a)² / |b - a|².
template <class K>
std::pair<typename K::FT, typename K::FT> site_squared_distance(
    const typename K::Point_2& station, const typename K::Point_2& source, const typename K::Point_2& target, bool interior)
{
    using FT = typename K::FT;
    const auto offset = station - source;
    if (!interior) return {offset.squared_length(), FT(1)};
    const auto edge = target - source;
    return {CGAL::square(CGAL::determinant(edge, offset)), edge.squared_length()};
}

// compare(d(p, first site), d(p, second site)), cross-multiplied by the
// positive denominators.
template <class K>
struct Compare_site_distance_2 {
    using result_type = typename K::Comparison_result;
    result_type operator()(const typename K::Point_2& station,
                           const typename K::Point_2& first_source, const typename K::Point_2& first_target, bool first_interior,
                           const typename K::Point_2& second_source, const typename K::Point_2& second_target, bool second_interior) const
    {
        const auto [first_numerator, first_denominator] = site_squared_distance<K>(station, first_source, first_target, first_interior);
        const auto [second_numerator, second_denominator] = site_squared_distance<K>(station, second_source, second_target, second_interior);
        return CGAL::compare(first_numerator * second_denominator, second_numerator * first_denominator);
    }
};

// sqrt(c) - sqrt(m) <= gap for c = d(p, candidate), m = d(p, minimum), gap >= 0,
// without roots: with r = c - m - gap², the claim is c <= m or r <= 0 or
// r² <= 4·m·gap² — the Fraction path's _within_distance_gap as a disjunction.
// With c = Nc/Dc, m = Nm/Dm and R = r·Dc·Dm (Dc, Dm > 0):
//   c <= m       iff Nc·Dm <= Nm·Dc
//   r <= 0       iff R <= 0
//   r² <= 4mg²   iff R² <= 4·Nm·gap²·Dc²·Dm
// The disjunction is evaluated without branching so the interval stage never
// converts an uncertain comparison (no exception on filter failure).
template <class K>
struct Within_distance_gap_2 {
    using result_type = typename K::Boolean;
    result_type operator()(const typename K::Point_2& station,
                           const typename K::Point_2& candidate_source, const typename K::Point_2& candidate_target, bool candidate_interior,
                           const typename K::Point_2& minimum_source, const typename K::Point_2& minimum_target, bool minimum_interior,
                           const typename K::FT& gap) const
    {
        using FT = typename K::FT;
        const auto [nc, dc] = site_squared_distance<K>(station, candidate_source, candidate_target, candidate_interior);
        const auto [nm, dm] = site_squared_distance<K>(station, minimum_source, minimum_target, minimum_interior);
        const FT candidate_scaled = nc * dm;
        const FT minimum_scaled = nm * dc;
        const FT gap_squared = CGAL::square(gap);
        const FT remainder = candidate_scaled - minimum_scaled - gap_squared * dc * dm;
        const FT bound = FT(4) * nm * gap_squared * CGAL::square(dc) * dm;
        return (CGAL::compare(candidate_scaled, minimum_scaled) != CGAL::LARGER)
            | (CGAL::sign(remainder) != CGAL::POSITIVE)
            | (CGAL::compare(CGAL::square(remainder), bound) != CGAL::LARGER);
    }
};

using CompareSiteDistance = CGAL::Filtered_predicate<Compare_site_distance_2<ExactKernel>,
    Compare_site_distance_2<ApproximateKernel>, ToExact, ToApproximate>;
using WithinDistanceGap = CGAL::Filtered_predicate<Within_distance_gap_2<ExactKernel>,
    Within_distance_gap_2<ApproximateKernel>, ToExact, ToApproximate>;

enum class Foot { Source, Interior, Target };

// Clamp classification of the projection parameter t = (p-a)·(b-a)/|b-a|²,
// decided by the exact Epick angle predicate: t <= 0 iff angle(p, a, b) is not
// acute; t >= 1 iff angle(p, b, a) is not acute.
Foot classify_foot(const Point& station, const Segment& segment)
{
    if (CGAL::angle(station, segment.source(), segment.target()) != CGAL::ACUTE) return Foot::Source;
    if (CGAL::angle(station, segment.target(), segment.source()) != CGAL::ACUTE) return Foot::Target;
    return Foot::Interior;
}

// Certified upper bound of the exact squared station-segment distance.
double squared_distance_upper(const Point& station, const Segment& segment, Foot foot)
{
    const Interval px(station.x()), py(station.y());
    const Interval ax(segment.source().x()), ay(segment.source().y());
    const Interval bx(segment.target().x()), by(segment.target().y());
    if (foot == Foot::Source) return (CGAL::square(px - ax) + CGAL::square(py - ay)).sup();
    if (foot == Foot::Target) return (CGAL::square(px - bx) + CGAL::square(py - by)).sup();
    const Interval ex = bx - ax, ey = by - ay;
    const Interval cross = ex * (py - ay) - ey * (px - ax);
    return (CGAL::square(cross) / (CGAL::square(ex) + CGAL::square(ey))).sup();
}

// Exact injection: value = mantissa · 2^exponent with an integer mantissa.
struct Dyadic {
    std::int64_t mantissa;
    int exponent;
};

Dyadic dyadic(double value)
{
    if (value == 0.0) return {0, 0};
    int exponent = 0;
    const double fraction = std::frexp(value, &exponent);
    constexpr int mantissa_bits = 53;  // IEEE-754 binary64 significand, exact in int64
    return {static_cast<std::int64_t>(std::ldexp(fraction, mantissa_bits)), exponent - mantissa_bits};
}

Integer scaled(const Dyadic& value, int exponent)
{
    // Shift the magnitude; |mantissa| < 2^53, so negation cannot overflow.
    Integer result(value.mantissa < 0 ? -value.mantissa : value.mantissa);
    result <<= static_cast<unsigned>(value.exponent - exponent);
    return value.mantissa < 0 ? Integer(-result) : result;
}

// Exact closed form of one site; the vertex case is parameter zero at `source`.
AdmissibleSite exact_site(std::size_t segment_id, const Point& station, const Segment& segment, bool vertex)
{
    const std::array<Dyadic, 6> inputs{dyadic(station.x()), dyadic(station.y()),
        dyadic(segment.source().x()), dyadic(segment.source().y()),
        dyadic(segment.target().x()), dyadic(segment.target().y())};
    std::optional<int> common;
    for (const auto& input : inputs) {
        if (input.mantissa != 0) common = common ? std::min(*common, input.exponent) : input.exponent;
    }
    const int exponent = common.value_or(0);
    const Integer px = scaled(inputs[0], exponent), py = scaled(inputs[1], exponent);
    const Integer ax = scaled(inputs[2], exponent), ay = scaled(inputs[3], exponent);
    const Integer wx = px - ax, wy = py - ay;
    if (vertex) {
        return {segment_id, Integer(0), Integer(1), ax, ay, wx * wx + wy * wy, exponent};
    }
    const Integer ex = scaled(inputs[4], exponent) - ax, ey = scaled(inputs[5], exponent) - ay;
    const Integer along = wx * ex + wy * ey;
    const Integer length_squared = ex * ex + ey * ey;
    return {segment_id, along, length_squared, ax * length_squared + along * ex, ay * length_squared + along * ey,
            (wx * wx + wy * wy) * length_squared - along * along, exponent};
}

void require_finite(const HeldBoundaryClearance2::XY& point, const char* message)
{
    if (!std::isfinite(point[0]) || !std::isfinite(point[1])) {
        throw InvalidHeldBoundaryClearanceInputError(message);
    }
}

// One candidate site: a vertex (represented by the segment it starts) or the
// interior foot of one segment. Vertices reached from both incident segments
// collapse to one identity before any distance is compared.
struct Candidate {
    std::size_t segment;
    bool vertex;
};

} // namespace

HeldBoundaryClearance2::HeldBoundaryClearance2(const std::vector<XY>& ring)
{
    if (ring.size() < 2) {
        throw InvalidHeldBoundaryClearanceInputError("Boundary clearance requires a closed ring of at least two vertices.");
    }
    for (const auto& point : ring) require_finite(point, "Boundary clearance ring coordinates must be finite.");
    segments_.reserve(ring.size());
    hint_segments_.reserve(ring.size());
    for (std::size_t index = 0; index < ring.size(); ++index) {
        const XY& start = ring[index];
        const XY& end = ring[(index + 1) % ring.size()];
        const Point source(start[0], start[1]);
        const Point target(end[0], end[1]);
        if (source == target) {
            throw DegenerateHeldBoundarySegmentError("Boundary clearance ring contains a zero-length segment.");
        }
        segments_.emplace_back(source, target);
        hint_segments_.emplace_back(HintKernel::Point_2(start[0], start[1]), HintKernel::Point_2(end[0], end[1]));
    }
    tree_.insert(segments_.cbegin(), segments_.cend());
    tree_.build();
    hint_tree_.insert(hint_segments_.cbegin(), hint_segments_.cend());
    hint_tree_.build();
    hint_tree_.accelerate_distance_queries();
}

std::vector<AdmissibleSite> HeldBoundaryClearance2::station_sites(const XY& station, double gap) const
{
    const Point point(station[0], station[1]);
    const std::size_t count = segments_.size();
    // Pruning only: any segment's certified distance bound bounds the minimum.
    const auto hint = hint_tree_.closest_point_and_primitive(HintKernel::Point_2(station[0], station[1]));
    const auto hint_index = static_cast<std::size_t>(std::distance(hint_segments_.cbegin(), hint.second));
    const Segment& hint_segment = segments_[hint_index];
    const double upper = squared_distance_upper(point, hint_segment, classify_foot(point, hint_segment));
    const Interval reach = CGAL::sqrt(Interval(upper)) + Interval(gap);
    const Interval x(station[0]), y(station[1]);
    const std::array<double, 4> corners{(x - reach).inf(), (y - reach).inf(), (x + reach).sup(), (y + reach).sup()};

    std::vector<FilterTree::Primitive_id> hits;
    if (std::all_of(corners.begin(), corners.end(), [](double value) { return std::isfinite(value); })) {
        tree_.all_intersected_primitives(Kernel::Iso_rectangle_2(corners[0], corners[1], corners[2], corners[3]),
                                         std::back_inserter(hits));
    } else {
        // An overflowing bound (extreme coordinates) prunes nothing: every
        // segment stays a candidate and is decided exactly below.
        for (auto segment = segments_.cbegin(); segment != segments_.cend(); ++segment) hits.push_back(segment);
    }
    if (hits.empty()) {
        throw BrokenHeldBoundaryClearancePruningError("Boundary clearance pruning box lost the hint segment.");
    }

    std::vector<Candidate> candidates;
    candidates.reserve(hits.size());
    for (const auto& hit : hits) {
        const auto index = static_cast<std::size_t>(std::distance(segments_.cbegin(), hit));
        switch (classify_foot(point, segments_[index])) {
        case Foot::Source: candidates.push_back({index, true}); break;
        case Foot::Target: candidates.push_back({(index + 1) % count, true}); break;
        case Foot::Interior: candidates.push_back({index, false}); break;
        }
    }
    // Canonical order (segment, parameter): a vertex site has parameter zero
    // and precedes the interior site of the segment it starts.
    std::sort(candidates.begin(), candidates.end(), [](const Candidate& first, const Candidate& second) {
        if (first.segment != second.segment) return first.segment < second.segment;
        return first.vertex && !second.vertex;
    });
    candidates.erase(std::unique(candidates.begin(), candidates.end(), [](const Candidate& first, const Candidate& second) {
        return first.segment == second.segment && first.vertex == second.vertex;
    }), candidates.end());

    const CompareSiteDistance compare_distance;
    const WithinDistanceGap within_gap;
    const auto source = [&](const Candidate& candidate) { return segments_[candidate.segment].source(); };
    const auto target = [&](const Candidate& candidate) { return segments_[candidate.segment].target(); };
    const Candidate* minimum = &candidates.front();
    for (const auto& candidate : candidates) {
        if (compare_distance(point, source(candidate), target(candidate), !candidate.vertex,
                             source(*minimum), target(*minimum), !minimum->vertex) == CGAL::SMALLER) {
            minimum = &candidate;
        }
    }
    std::vector<AdmissibleSite> sites;
    for (const auto& candidate : candidates) {
        if (&candidate != minimum
            && !within_gap(point, source(candidate), target(candidate), !candidate.vertex,
                           source(*minimum), target(*minimum), !minimum->vertex, gap)) {
            continue;
        }
        sites.push_back(exact_site(candidate.segment, point, segments_[candidate.segment], candidate.vertex));
    }
    return sites;
}

std::vector<std::vector<AdmissibleSite>> HeldBoundaryClearance2::distance_admissible_sites(
    const std::vector<XY>& stations, double gap) const
{
    if (!std::isfinite(gap) || gap < 0) {
        throw InvalidHeldBoundaryClearanceInputError("Boundary clearance distance gap must be finite and non-negative.");
    }
    for (const auto& station : stations) require_finite(station, "Boundary clearance station coordinates must be finite.");
    std::vector<std::vector<AdmissibleSite>> result;
    result.reserve(stations.size());
    for (const auto& station : stations) result.push_back(station_sites(station, gap));
    return result;
}
