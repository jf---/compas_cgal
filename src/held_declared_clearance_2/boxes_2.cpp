#include "held_declared_clearance_2/boxes_2.h"

#include <algorithm>
#include <utility>

namespace compas_cgal::held_declared_clearance {

namespace {

Box hull(const std::pair<double, double>& x0, const std::pair<double, double>& y0, const std::pair<double, double>& x1,
         const std::pair<double, double>& y1)
{
    return {std::min(x0.first, x1.first), std::min(y0.first, y1.first), std::max(x0.second, x1.second),
            std::max(y0.second, y1.second)};
}

void include(Box& box, const std::pair<double, double>& x, const std::pair<double, double>& y)
{
    box.xmin = std::min(box.xmin, x.first);
    box.ymin = std::min(box.ymin, y.first);
    box.xmax = std::max(box.xmax, x.second);
    box.ymax = std::max(box.ymax, y.second);
}

/// The ends' box grown by every axis extreme centre +- radius e_i whose
/// direction the sweep spans (`spans(direction)`).
template <typename Spans>
Box swept_box(Box box, const RationalPoint& centre, const FT& radius, Spans spans)
{
    const FT zero(0);
    if (spans(Vector(radius, zero))) include(box, CGAL::to_interval(centre.x() + radius), CGAL::to_interval(centre.y()));
    if (spans(Vector(-radius, zero))) include(box, CGAL::to_interval(centre.x() - radius), CGAL::to_interval(centre.y()));
    if (spans(Vector(zero, radius))) include(box, CGAL::to_interval(centre.x()), CGAL::to_interval(centre.y() + radius));
    if (spans(Vector(zero, -radius))) include(box, CGAL::to_interval(centre.x()), CGAL::to_interval(centre.y() - radius));
    return box;
}

}  // namespace

Box disk_box(const RationalPoint& centre, const FT& radius)
{
    return hull(CGAL::to_interval(centre.x() - radius), CGAL::to_interval(centre.y() - radius),
                CGAL::to_interval(centre.x() + radius), CGAL::to_interval(centre.y() + radius));
}

Box primitive_box(const Primitive& primitive)
{
    const Box ends = hull(CGAL::to_interval(primitive.start.x()), CGAL::to_interval(primitive.start.y()),
                          CGAL::to_interval(primitive.end.x()), CGAL::to_interval(primitive.end.y()));
    if (!primitive.is_arc()) return ends;
    return swept_box(ends, primitive.centre, primitive.radius,
                     [&](const Vector& direction) { return spans_direction(primitive, direction); });
}

Box arc_box(const ArcElement& element)
{
    if (element.full) return disk_box(element.centre, element.radius);
    const Box ends = hull(CGAL::to_interval(element.start.x()), CGAL::to_interval(element.start.y()),
                          CGAL::to_interval(element.end.x()), CGAL::to_interval(element.end.y()));
    return swept_box(ends, element.centre, element.radius,
                     [&](const Vector& direction) { return element.spans(embed(element.centre + direction)); });
}

Box segment_box(const SegmentElement& element)
{
    return hull(CGAL::to_interval(element.start.x()), CGAL::to_interval(element.start.y()),
                CGAL::to_interval(element.end.x()), CGAL::to_interval(element.end.y()));
}

}  // namespace compas_cgal::held_declared_clearance
