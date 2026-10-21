#include "held_declared_clearance_2/candidates_2.h"

namespace compas_cgal::held_declared_clearance {

std::vector<Distance> segment_candidates(const SegmentElement& element, const Primitive& primitive)
{
    std::vector<Distance> candidates;
    candidates.reserve(6);
    candidates.push_back(distance_to(primitive, element.start));
    candidates.push_back(distance_to(primitive, element.end));
    candidates.push_back(element.distance_to_point(primitive.start));
    candidates.push_back(element.distance_to_point(primitive.end));
    if (!primitive.is_arc() || !element.contains_parameter_of(primitive.centre)) return candidates;
    // Arc points on the line's normal through the arc centre: C -+ rho * n_hat.
    const FT norm_squared = element.normal.squared_length();
    const FT height = element.normal * (primitive.centre - CGAL::ORIGIN) - element.level;
    for (const int side : {1, -1}) {
        const RationalPoint probe = primitive.centre - FT(side) * element.normal;
        if (!spans(primitive, embed(probe))) continue;
        // Signed distance of C - side * rho * n_hat to the line: h / |n| - side * rho.
        candidates.push_back(Distance::absolute(make_root(-side * primitive.radius, height / norm_squared, norm_squared)));
    }
    return candidates;
}

std::vector<Distance> arc_candidates(const ArcElement& element, const Primitive& primitive)
{
    std::vector<Distance> candidates;
    candidates.reserve(8);
    candidates.push_back(element.distance_to_point(primitive.start));
    candidates.push_back(element.distance_to_point(primitive.end));
    if (!element.full) {
        candidates.push_back(distance_to(primitive, element.start));
        candidates.push_back(distance_to(primitive, element.end));
    }
    if (!primitive.is_arc()) {
        const Vector offset = element.centre - primitive.start;
        const FT position = primitive.unit * offset;
        if (CGAL::is_negative(position) || position > primitive.length) return candidates;
        const FT height = primitive.normal * offset;
        for (const int side : {1, -1}) {
            const RationalPoint probe = element.centre - FT(side) * primitive.normal;
            if (element.spans(embed(probe))) candidates.push_back(Distance::absolute(rational(height - side * element.radius)));
        }
        return candidates;
    }
    const Vector between = primitive.centre - element.centre;
    const FT square = between.squared_length();
    if (CGAL::is_zero(square)) return candidates;
    for (const int own : {1, -1}) {
        if (!element.spans(embed(element.centre + FT(own) * between))) continue;
        for (const int other : {1, -1}) {
            // x_E = C' + own rho' v_hat, x_X = C + other rho v_hat (v = C - C'):
            // |x_X - x_E| = |sqrt(V) + other rho - own rho'|.
            if (spans(primitive, embed(primitive.centre + FT(other) * between))) {
                candidates.push_back(
                    Distance::shifted_root(rational(square), other * primitive.radius - own * element.radius));
            }
        }
    }
    return candidates;
}

}  // namespace compas_cgal::held_declared_clearance
