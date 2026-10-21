#include "held_declared_clearance_2/primitives_2.h"

#include "held_declared_clearance_2/errors_2.h"

#include <string>

namespace compas_cgal::held_declared_clearance {

namespace {

/// Exact u.(p - start) of a line, in p's extension.
Root along(const Primitive& line, const OneRootPoint& point)
{
    return plus(scaled(line.unit.x(), shifted(point.x(), -line.start.x())),
                scaled(line.unit.y(), shifted(point.y(), -line.start.y())));
}

/// Exact n.(p - start) of a line, in p's extension.
Root height(const Primitive& line, const OneRootPoint& point)
{
    return plus(scaled(line.normal.x(), shifted(point.x(), -line.start.x())),
                scaled(line.normal.y(), shifted(point.y(), -line.start.y())));
}

std::string where(std::size_t index) { return "declared primitive " + std::to_string(index) + ": "; }

}  // namespace

Root squared_distance(const OneRootPoint& point, const RationalPoint& origin)
{
    const Root dx = shifted(point.x(), -origin.x());
    const Root dy = shifted(point.y(), -origin.y());
    return plus(times(dx, dx), times(dy, dy));
}

const Root& smaller(const Root& first, const Root& second) { return compare(second, first) < 0 ? second : first; }

Primitive Primitive::line(std::size_t index, const RationalPoint& start, const RationalPoint& end, const Vector& unit,
                          const FT& length)
{
    if (unit.squared_length() != FT(1)) throw InvalidDeclaredPrimitiveError(where(index) + "direction is not exactly unit");
    if (!CGAL::is_positive(length)) throw InvalidDeclaredPrimitiveError(where(index) + "length is not positive");
    if (start + length * unit != end) throw InvalidDeclaredPrimitiveError(where(index) + "end is not start + length * u");
    Primitive primitive;
    primitive.index = index;
    primitive.kind = PrimitiveKind::Line;
    primitive.start = start;
    primitive.end = end;
    primitive.unit = unit;
    primitive.normal = unit.perpendicular(CGAL::COUNTERCLOCKWISE);
    primitive.length = length;
    return primitive;
}

Primitive Primitive::arc(std::size_t index, const RationalPoint& start, const RationalPoint& end,
                         const RationalPoint& centre, const FT& radius, int sense)
{
    if (!CGAL::is_positive(radius)) throw InvalidDeclaredPrimitiveError(where(index) + "radius is not positive");
    if (sense != 1 && sense != -1) throw InvalidDeclaredPrimitiveError(where(index) + "sense is not +1 or -1");
    const FT squared = radius * radius;
    if (CGAL::compare_squared_distance(start, centre, squared) != CGAL::EQUAL ||
        CGAL::compare_squared_distance(end, centre, squared) != CGAL::EQUAL) {
        throw InvalidDeclaredPrimitiveError(where(index) + "an end is not exactly on the circle");
    }
    Primitive primitive;
    primitive.index = index;
    primitive.kind = PrimitiveKind::Arc;
    primitive.start = start;
    primitive.end = end;
    primitive.centre = centre;
    primitive.radius = radius;
    primitive.sense = sense;
    return primitive;
}

bool spans_direction(const Primitive& arc, const Vector& direction)
{
    return compare_on_arc(arc.sense, arc.start - arc.centre, direction, arc.end - arc.centre) <= 0;
}

bool spans(const Primitive& arc, const OneRootPoint& probe)
{
    if (is_rational(probe)) return spans_direction(arc, rational_point(probe) - arc.centre);
    return compare_on_arc(arc.centre, arc.sense, embed(arc.start), probe, embed(arc.end)) <= 0;
}

Distance distance_to(const Primitive& primitive, const OneRootPoint& point)
{
    if (!primitive.is_arc()) {
        const Root position = along(primitive, point);
        if (sign(position) >= 0 && sign(shifted(position, -primitive.length)) <= 0) {
            return Distance::absolute(height(primitive, point));
        }
        return Distance::root(smaller(squared_distance(point, primitive.start), squared_distance(point, primitive.end)));
    }
    Root square = squared_distance(point, primitive.centre);
    if (sign(square) == 0) return Distance::absolute(rational(primitive.radius));
    if (spans(primitive, point)) return Distance::shifted_root(std::move(square), -primitive.radius);
    return Distance::root(smaller(squared_distance(point, primitive.start), squared_distance(point, primitive.end)));
}

Distance nearest(const Primitive& primitive, const RationalPoint& centre)
{
    if (!primitive.is_arc()) {
        const Vector offset = centre - primitive.start;
        const FT position = primitive.unit * offset;
        if (!CGAL::is_negative(position) && position <= primitive.length) {
            return Distance::absolute(rational(primitive.normal * offset));
        }
        const RationalPoint& end = CGAL::is_negative(position) ? primitive.start : primitive.end;
        return Distance::root(rational(CGAL::squared_distance(centre, end)));
    }
    const Vector offset = centre - primitive.centre;
    const FT square = offset.squared_length();
    if (CGAL::is_zero(square)) return Distance::absolute(rational(primitive.radius));
    if (spans_direction(primitive, offset)) return Distance::shifted_root(rational(square), -primitive.radius);
    const FT to_start = CGAL::squared_distance(centre, primitive.start);
    const FT to_end = CGAL::squared_distance(centre, primitive.end);
    return Distance::root(rational(to_end < to_start ? to_end : to_start));
}

std::vector<Distance> farthest(const Primitive& primitive, const RationalPoint& centre)
{
    std::vector<Distance> found;
    if (primitive.is_arc() && centre == primitive.centre) {
        found.push_back(Distance::absolute(rational(primitive.radius)));
        return found;
    }
    found.push_back(Distance::root(rational(CGAL::squared_distance(centre, primitive.start))));
    found.push_back(Distance::root(rational(CGAL::squared_distance(centre, primitive.end))));
    if (!primitive.is_arc()) return found;
    const Vector offset = centre - primitive.centre;
    if (spans_direction(primitive, -offset)) {
        found.push_back(Distance::shifted_root(rational(offset.squared_length()), primitive.radius));
    }
    return found;
}

}  // namespace compas_cgal::held_declared_clearance
