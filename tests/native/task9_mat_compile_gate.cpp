#include "segment_site_mat.h"

#include <cstdlib>
#include <optional>

namespace {

bool parameter_domains_are_exact()
{
    using Point = MatTraits::Point_2;
    const MatTraits::Line_2 line(Point(0, 0), Point(1, 0));
    const MatTraits::Ray_2 ray(Point(0, 0), Point(1, 0));
    const MatTraits::Segment_2 segment(
        Point(0, 0),
        Point(1, 0));
    const SegmentSiteParabola2 parabola(
        Point(0, 1),
        MatTraits::Line_2(0, 1, 1),
        Point(-2, 1),
        Point(2, 1));

    const MatParameterDomain2 line_domain =
        exact_parameter_domain(line);
    const MatParameterDomain2 ray_domain =
        exact_parameter_domain(ray);
    const MatParameterDomain2 segment_domain =
        exact_parameter_domain(segment);
    const MatParameterDomain2 parabola_domain =
        exact_parameter_domain(parabola);
    return !line_domain.lower.has_value()
        && !line_domain.upper.has_value()
        && ray_domain.lower == MatTraits::FT(0)
        && !ray_domain.upper.has_value()
        && segment_domain.lower == MatTraits::FT(0)
        && segment_domain.upper == MatTraits::FT(1)
        && parabola_domain.lower.has_value()
        && parabola_domain.upper.has_value()
        && CGAL::compare(
               *parabola_domain.lower,
               *parabola_domain.upper)
            != CGAL::LARGER;
}

bool clearance_boundaries_are_exact()
{
    const RationalPrimitiveParameterization2 line{
        {CORE::BigRat(0), CORE::BigRat(1)},
        {CORE::BigRat(0)},
        std::nullopt,
        std::nullopt,
    };
    const RationalPrimitiveParameterization2 ray{
        line.x_coefficients,
        line.y_coefficients,
        CORE::BigRat(0),
        std::nullopt,
    };
    const RationalPrimitiveParameterization2 segment{
        {CORE::BigRat(0), CORE::BigRat(2)},
        {CORE::BigRat(0)},
        CORE::BigRat(0),
        CORE::BigRat(1),
    };
    const RationalPrimitiveParameterization2 parabola{
        {CORE::BigRat(0), CORE::BigRat(1)},
        {CORE::BigRat(0), CORE::BigRat(0), CORE::BigRat(1)},
        CORE::BigRat(-1),
        CORE::BigRat(1),
    };

    const ClearanceRootBoundary2 line_boundary =
        point_clearance_boundary(line, 0, 0, 1);
    const ClearanceRootBoundary2 ray_boundary =
        point_clearance_boundary(ray, 0, 0, 1);
    const ClearanceRootBoundary2 segment_boundary =
        point_clearance_boundary(segment, 0, 0, 1);
    const ClearanceRootBoundary2 parabola_boundary =
        point_clearance_boundary(parabola, 0, 0, 1);
    const ClearanceRootBoundary2 zero =
        point_clearance_boundary(
            {{CORE::BigRat(0)}, {CORE::BigRat(0)}, 0, 1},
            0,
            0,
            0);
    const ClearanceRootBoundary2 positive =
        point_clearance_boundary(
            {{CORE::BigRat(1)}, {CORE::BigRat(0)}, 0, 1},
            0,
            0,
            0);
    const ClearanceRootBoundary2 negative =
        point_clearance_boundary(
            {{CORE::BigRat(0)}, {CORE::BigRat(0)}, 0, 1},
            0,
            0,
            1);

    return line_boundary.roots.size() == 2
        && ray_boundary.roots.size() == 1
        && segment_boundary.roots.size() == 1
        && parabola_boundary.roots.size() == 2
        && zero.constant_sign == CGAL::ZERO
        && positive.constant_sign == CGAL::POSITIVE
        && negative.constant_sign == CGAL::NEGATIVE;
}

} // namespace

int main()
{
    const SegmentSiteMatCompileEvidence2 evidence =
        segment_site_mat_compile_spike();
    return evidence.delaunay_valid
            && evidence.assigned_dual_primitives > 0
            && evidence.exact_clearance_roots == 2
            && evidence.matched_generator_sites
                == evidence.delaunay_vertices
            && parameter_domains_are_exact()
            && clearance_boundaries_are_exact()
        ? EXIT_SUCCESS
        : EXIT_FAILURE;
}
