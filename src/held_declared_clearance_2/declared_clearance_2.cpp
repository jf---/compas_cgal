#include "held_declared_clearance_2/declared_clearance_2.h"

#include "held_declared_clearance_2/boxes_2.h"
#include "held_declared_clearance_2/candidates_2.h"
#include "held_declared_clearance_2/crossings_2.h"
#include "held_declared_clearance_2/errors_2.h"
#include "held_declared_clearance_2/owner_lemma_2.h"
#include "held_declared_clearance_2/reporting_2.h"

#include <algorithm>
#include <utility>

namespace compas_cgal::held_declared_clearance {

namespace {

/// One decided (non-CLEAR) primitive: exact sign, reported distance.
struct Decided {
    int sign;
    double distance;
    std::size_t index;
};

/// The Python oracle's selection: worst sign, then (reported distance, index).
ClearanceVerdict select(const std::vector<Decided>& decided, double tool_radius)
{
    if (decided.empty()) return {Clearance::Clear, -1, 0.0};
    int worst = 1;
    for (const Decided& item : decided) worst = std::min(worst, item.sign);
    const Decided* best = nullptr;
    for (const Decided& item : decided) {
        if (item.sign != worst) continue;
        if (best == nullptr || item.distance < best->distance ||
            (item.distance == best->distance && item.index < best->index)) {
            best = &item;
        }
    }
    const double depth = worst < 0 ? tool_radius - best->distance : 0.0;
    return {static_cast<Clearance>(worst), static_cast<std::ptrdiff_t>(best->index), depth};
}

/// Python's min(candidate.approx() for candidate in candidates).
double reported_minimum(const std::vector<Distance>& candidates)
{
    double result = candidates.front().approx();
    for (std::size_t i = 1; i < candidates.size(); ++i) {
        const double value = candidates[i].approx();
        if (value < result) result = value;
    }
    return result;
}

/// Python's max(a, b, c): the first maximal argument.
double reported_maximum(double first, double second, double third)
{
    double result = first;
    if (second > result) result = second;
    if (third > result) result = third;
    return result;
}

int minimum_sign(const std::vector<Distance>& candidates, const FT& threshold)
{
    int result = 1;
    for (const Distance& candidate : candidates) result = std::min(result, candidate.sign_minus(threshold));
    return result;
}

bool listed(const std::vector<std::size_t>& owners, std::size_t index)
{
    return std::find(owners.begin(), owners.end(), index) != owners.end();
}

}  // namespace

std::unique_ptr<DeclaredClearance2> DeclaredClearance2::build(std::vector<Primitive> primitives, const FT& tool_radius)
{
    if (!CGAL::is_positive(tool_radius)) throw NonPositiveClearanceRadiusError("the tool radius is not positive");
    const double reach = binary64(tool_radius);
    if (FT(reach) != tool_radius) throw NonPositiveClearanceRadiusError("the tool radius is not a binary64");
    if (primitives.empty()) throw InvalidDeclaredPrimitiveError("a declared boundary needs primitives");
    for (std::size_t index = 0; index < primitives.size(); ++index) {
        if (primitives[index].index != index) throw InvalidDeclaredPrimitiveError("primitive indices must be 0..n-1 in order");
    }
    return std::unique_ptr<DeclaredClearance2>(new DeclaredClearance2(std::move(primitives), tool_radius, reach));
}

namespace {

std::vector<Box> boxes_of(const std::vector<Primitive>& primitives)
{
    std::vector<Box> boxes;
    boxes.reserve(primitives.size());
    for (const Primitive& primitive : primitives) boxes.push_back(primitive_box(primitive));
    return boxes;
}

}  // namespace

DeclaredClearance2::DeclaredClearance2(std::vector<Primitive> primitives, const FT& tool_radius, double reach)
    : primitives_(std::move(primitives)), tool_radius_(tool_radius), reach_(reach), index_(boxes_of(primitives_))
{
}

ClearanceVerdict DeclaredClearance2::circle(const CircleElement& element, ClearanceWorkspace& workspace) const
{
    if (!CGAL::is_positive(element.radius)) throw NonPositiveClearanceRadiusError("a circle radius is not positive");
    ++workspace.counters.elements;
    const FT outer = element.radius + tool_radius_;
    const FT inner = element.radius - tool_radius_;
    const double radius = binary64(element.radius);
    std::vector<Decided> decided;
    index_.near(disk_box(element.centre, element.radius), reach_, workspace.near);
    for (const std::size_t index : workspace.near) {
        ++workspace.counters.pairs;
        const Primitive& primitive = primitives_[index];
        const Distance low_form = nearest(primitive, element.centre);
        const int low_sign = low_form.sign_minus(outer);
        if (low_sign > 0) continue;
        const double low = low_form.approx();
        if (low_sign == 0) {
            // lo == R + r forces hi >= lo > R - r: a touch from the outside.
            decided.push_back({0, std::max(low - radius, 0.0), index});
            continue;
        }
        // lo < R + r: clear only when the whole primitive is inside the hole (hi < R - r).
        const std::vector<Distance> highs = farthest(primitive, element.centre);
        int high_sign = -1;
        for (const Distance& form : highs) high_sign = std::max(high_sign, form.sign_minus(inner));
        double high = highs.front().approx();
        for (std::size_t i = 1; i < highs.size(); ++i) {
            const double value = highs[i].approx();
            if (value > high) high = value;
        }
        const int sign = std::max({low_sign, -high_sign, -1});
        if (sign > 0) continue;
        decided.push_back({sign, reported_maximum(low - radius, radius - high, 0.0), index});
    }
    return select(decided, reach_);
}

ClearanceVerdict DeclaredClearance2::segment(const SegmentElement& element, const std::vector<std::size_t>& owners,
                                             ClearanceWorkspace& workspace) const
{
    ++workspace.counters.elements;
    std::vector<Decided> decided;
    index_.near(segment_box(element), reach_, workspace.near);
    for (const std::size_t index : workspace.near) {
        ++workspace.counters.pairs;
        const Primitive& primitive = primitives_[index];
        if (listed(owners, index) && line_owner_touches(element, primitive, tool_radius_)) {
            ++workspace.counters.owner_touches;
            decided.push_back({0, reported_minimum(segment_candidates(element, primitive)), index});
            continue;
        }
        if (segment_crosses(element, primitive)) {
            ++workspace.counters.crossings;
            decided.push_back({-1, 0.0, index});
            continue;
        }
        const std::vector<Distance> candidates = segment_candidates(element, primitive);
        const int sign = minimum_sign(candidates, tool_radius_);
        if (sign > 0) continue;
        decided.push_back({sign, reported_minimum(candidates), index});
    }
    return select(decided, reach_);
}

ClearanceVerdict DeclaredClearance2::arc(const ArcElement& element, const std::vector<std::size_t>& owners,
                                         ClearanceWorkspace& workspace) const
{
    ++workspace.counters.elements;
    std::vector<Decided> decided;
    index_.near(arc_box(element), reach_, workspace.near);
    for (const std::size_t index : workspace.near) {
        ++workspace.counters.pairs;
        const Primitive& primitive = primitives_[index];
        if (listed(owners, index) && arc_owner_touches(element, primitive, tool_radius_)) {
            ++workspace.counters.owner_touches;
            decided.push_back({0, reported_minimum(arc_candidates(element, primitive)), index});
            continue;
        }
        if (arc_crosses(element, primitive)) {
            ++workspace.counters.crossings;
            decided.push_back({-1, 0.0, index});
            continue;
        }
        const std::vector<Distance> candidates = arc_candidates(element, primitive);
        const int sign = minimum_sign(candidates, tool_radius_);
        if (sign > 0) continue;
        decided.push_back({sign, reported_minimum(candidates), index});
    }
    return select(decided, reach_);
}

}  // namespace compas_cgal::held_declared_clearance
