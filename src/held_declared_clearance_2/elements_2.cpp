#include "held_declared_clearance_2/elements_2.h"

#include "held_declared_clearance_2/errors_2.h"
#include "held_declared_clearance_2/primitives_2.h"

namespace compas_cgal::held_declared_clearance {

SegmentElement SegmentElement::build(const OneRootPoint& start, const OneRootPoint& end, const Vector& direction,
                                     const Vector& normal, const FT& level)
{
    if (direction == CGAL::NULL_VECTOR) throw DegenerateClearanceElementError("a segment element needs a direction");
    if (normal == CGAL::NULL_VECTOR) throw DegenerateClearanceElementError("a segment element needs a normal");
    if (!CGAL::is_zero(normal * direction)) {
        throw DegenerateClearanceElementError("a segment element's normal is not orthogonal to its direction");
    }
    return {start, end, direction, normal, level};
}

SegmentElement SegmentElement::between(const RationalPoint& start, const RationalPoint& end)
{
    if (start == end) throw DegenerateClearanceElementError("a segment element needs two distinct points");
    const Vector direction = end - start;
    const Vector normal = direction.perpendicular(CGAL::COUNTERCLOCKWISE);
    return {embed(start), embed(end), direction, normal, normal * (start - CGAL::ORIGIN)};
}

Root SegmentElement::position(const OneRootPoint& point) const
{
    return plus(scaled(direction.x(), point.x()), scaled(direction.y(), point.y()));
}

bool SegmentElement::contains_parameter_of(const RationalPoint& point) const
{
    const Root here = rational(direction * (point - CGAL::ORIGIN));
    return compare(position(start), here) <= 0 && compare(here, position(end)) <= 0;
}

Distance SegmentElement::distance_to_point(const RationalPoint& point) const
{
    if (contains_parameter_of(point)) {
        const FT offset = normal * (point - CGAL::ORIGIN) - level;
        return Distance::root(rational(offset * offset / normal.squared_length()));
    }
    return Distance::root(smaller(squared_distance(start, point), squared_distance(end, point)));
}

ArcElement ArcElement::build(const RationalPoint& centre, const FT& radius, int sense, const OneRootPoint& start,
                             const OneRootPoint& end)
{
    if (!CGAL::is_positive(radius)) throw DegenerateClearanceElementError("an arc element needs a positive radius");
    if (sense != 1 && sense != -1) throw DegenerateClearanceElementError("an arc element's sense is not +1 or -1");
    return {centre, radius, sense, start, end, compare_xy(start, end) == 0};
}

bool ArcElement::spans(const OneRootPoint& probe) const
{
    return full || compare_on_arc(centre, sense, start, probe, end) <= 0;
}

Distance ArcElement::distance_to_point(const RationalPoint& point) const
{
    const FT square = CGAL::squared_distance(point, centre);
    if (CGAL::is_zero(square)) return Distance::absolute(rational(radius));
    if (spans(embed(point))) return Distance::shifted_root(rational(square), -radius);
    return Distance::root(smaller(squared_distance(start, point), squared_distance(end, point)));
}

}  // namespace compas_cgal::held_declared_clearance
