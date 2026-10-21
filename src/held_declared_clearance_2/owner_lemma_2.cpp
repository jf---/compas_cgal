#include "held_declared_clearance_2/owner_lemma_2.h"

namespace compas_cgal::held_declared_clearance {

namespace {

/// n.(p - S) == r and u.(p - S) in [0, Lambda], exactly.
bool on_owner_offset(const OneRootPoint& point, const Primitive& line, const FT& tool_radius)
{
    const Root dx = shifted(point.x(), -line.start.x());
    const Root dy = shifted(point.y(), -line.start.y());
    if (sign(shifted(plus(scaled(line.normal.x(), dx), scaled(line.normal.y(), dy)), -tool_radius)) != 0) return false;
    const Root along = plus(scaled(line.unit.x(), dx), scaled(line.unit.y(), dy));
    return sign(along) >= 0 && sign(shifted(along, -line.length)) <= 0;
}

}  // namespace

bool line_owner_touches(const SegmentElement& element, const Primitive& owner, const FT& tool_radius)
{
    if (owner.is_arc()) return false;
    return on_owner_offset(element.start, owner, tool_radius) && on_owner_offset(element.end, owner, tool_radius);
}

bool arc_owner_touches(const ArcElement& element, const Primitive& owner, const FT& tool_radius)
{
    if (!owner.is_arc() || element.full || element.centre != owner.centre) return false;
    if (CGAL::compare(CGAL::square(element.radius - owner.radius), CGAL::square(tool_radius)) != CGAL::EQUAL) return false;
    // The element's sweep, read in the owner's sense, runs from `first` to `last`.
    const bool same_sense = element.sense == owner.sense;
    const OneRootPoint& first = same_sense ? element.start : element.end;
    const OneRootPoint& last = same_sense ? element.end : element.start;
    const OneRootPoint start = embed(owner.start);
    const OneRootPoint end = embed(owner.end);
    return compare_on_arc(owner.centre, owner.sense, start, first, end) <= 0 &&
           compare_on_arc(owner.centre, owner.sense, start, last, end) <= 0 &&
           compare_on_arc(owner.centre, owner.sense, start, first, last) <= 0;
}

}  // namespace compas_cgal::held_declared_clearance
