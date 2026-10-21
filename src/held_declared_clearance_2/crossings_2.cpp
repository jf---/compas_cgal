#include "held_declared_clearance_2/crossings_2.h"

#include <utility>
#include <vector>

namespace compas_cgal::held_declared_clearance {

namespace {

struct LinePoint {
    Root along;  // the point is anchor + along * direction
    OneRootPoint point;
};

/// Points anchor + lambda * direction on the circle (centre, radius), with their lambda.
std::vector<LinePoint> line_circle_points(const RationalPoint& anchor, const Vector& direction,
                                          const RationalPoint& centre, const FT& radius)
{
    const Vector offset = anchor - centre;
    const FT a = direction.squared_length();
    const FT b = 2 * (direction * offset);
    const FT c = offset.squared_length() - radius * radius;
    const FT discriminant = b * b - 4 * a * c;
    std::vector<LinePoint> found;
    if (CGAL::is_negative(discriminant)) return found;
    const FT twice = 2 * a;
    for (const int side : {1, -1}) {
        if (side == -1 && CGAL::is_zero(discriminant)) break;
        Root along = make_deciding_root(-b / twice, FT(side) / twice, discriminant);
        OneRootPoint point = make_point(shifted(scaled(direction.x(), along), anchor.x()),
                                        shifted(scaled(direction.y(), along), anchor.y()));
        found.push_back({std::move(along), std::move(point)});
    }
    return found;
}

/// Intersection points of two non-concentric circles.
std::vector<OneRootPoint> circle_circle_points(const RationalPoint& first, const FT& first_radius,
                                               const RationalPoint& second, const FT& second_radius)
{
    std::vector<OneRootPoint> found;
    const Vector between = second - first;
    const FT distance = between.squared_length();
    if (CGAL::is_zero(distance)) return found;
    const FT first_squared = first_radius * first_radius;
    const FT along = (first_squared - second_radius * second_radius + distance) / (2 * distance);
    const FT height = first_squared / distance - along * along;
    if (CGAL::is_negative(height)) return found;
    const RationalPoint base = first + along * between;
    for (const int side : {1, -1}) {
        if (side == -1 && CGAL::is_zero(height)) break;
        found.push_back(make_point(make_deciding_root(base.x(), -side * between.y(), height),
                                   make_deciding_root(base.y(), side * between.x(), height)));
    }
    return found;
}

/// low <= value <= high, exactly across extensions.
bool within(const Root& low, const Root& value, const Root& high)
{
    return compare(low, value) <= 0 && compare(value, high) <= 0;
}

}  // namespace

bool segment_crosses(const SegmentElement& element, const Primitive& primitive)
{
    const Root low = element.position(element.start);
    const Root high = element.position(element.end);
    if (!primitive.is_arc()) {
        const FT start_side = element.normal * (primitive.start - CGAL::ORIGIN) - element.level;
        const FT end_side = element.normal * (primitive.end - CGAL::ORIGIN) - element.level;
        if (sign(start_side) * sign(end_side) > 0) return false;
        if (CGAL::is_zero(start_side) && CGAL::is_zero(end_side)) {
            FT first = element.direction * (primitive.start - CGAL::ORIGIN);
            FT second = element.direction * (primitive.end - CGAL::ORIGIN);
            if (second < first) std::swap(first, second);
            return !(compare(low, rational(second)) > 0 || compare(high, rational(first)) < 0);
        }
        const FT fraction = start_side / (start_side - end_side);
        const RationalPoint crossing = primitive.start + fraction * (primitive.end - primitive.start);
        return within(low, rational(element.direction * (crossing - CGAL::ORIGIN)), high);
    }
    const Vector anchor = (element.level / element.normal.squared_length()) * element.normal;
    const FT base = element.direction * anchor;
    const FT scale = element.direction.squared_length();
    for (const LinePoint& found :
         line_circle_points(CGAL::ORIGIN + anchor, element.direction, primitive.centre, primitive.radius)) {
        if (within(low, shifted(scaled(scale, found.along), base), high) && spans(primitive, found.point)) return true;
    }
    return false;
}

bool arc_crosses(const ArcElement& element, const Primitive& primitive)
{
    if (!primitive.is_arc()) {
        for (const LinePoint& found : line_circle_points(primitive.start, primitive.unit, element.centre, element.radius)) {
            if (sign(found.along) >= 0 && compare(found.along, rational(primitive.length)) <= 0 &&
                element.spans(found.point)) {
                return true;
            }
        }
        return false;
    }
    for (const OneRootPoint& point : circle_circle_points(element.centre, element.radius, primitive.centre, primitive.radius)) {
        if (element.spans(point) && spans(primitive, point)) return true;
    }
    return false;
}

}  // namespace compas_cgal::held_declared_clearance
