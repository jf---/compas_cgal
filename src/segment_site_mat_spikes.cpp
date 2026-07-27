#include "segment_site_mat.h"

#include <iterator>
#include <vector>

#include <CGAL/Object.h>
#include <CGAL/Polynomial_traits_d.h>

namespace {

std::size_t exact_clearance_root_count()
{
    using Polynomial =
        ExactAlgebraicKernel1::Polynomial_1;
    const std::vector<ExactAlgebraicInteger1>
        coefficients{-1, 0, 1};
    const Polynomial clearance =
        typename CGAL::Polynomial_traits_d<
            Polynomial>::Construct_polynomial()(
            coefficients.begin(),
            coefficients.end());
    std::vector<
        ExactAlgebraicKernel1::Algebraic_real_1>
        roots;
    ExactAlgebraicKernel1 kernel;
    kernel.solve_1_object()(
        clearance,
        true,
        std::back_inserter(roots));
    return roots.size();
}

std::size_t assign_dual_primitive(
    const CGAL::Object& dual)
{
    MatTraits::Line_2 line;
    MatTraits::Ray_2 ray;
    MatTraits::Segment_2 segment;
    SegmentSiteParabola2 parabola;
    if (CGAL::assign(line, dual)
        || CGAL::assign(ray, dual)
        || CGAL::assign(segment, dual)) {
        return 1;
    }
    if (!CGAL::assign(parabola, dual)) {
        return 0;
    }

    // `t` and `f` are the exact algebraic parameterization API. Drawing
    // helpers (`compute_k`, `generate_points`, streaming) are forbidden.
    const MatTraits::FT first_parameter =
        parabola.t(parabola.p1);
    const MatTraits::Point_2 reconstructed =
        parabola.f(first_parameter);
    return reconstructed == parabola.p1 ? 1 : 0;
}

} // namespace

SegmentSiteMatCompileEvidence2
segment_site_mat_compile_spike()
{
    using Point = MatTraits::Point_2;
    using Site = MatTraits::Site_2;
    const Point lower_left(0, 0);
    const Point lower_right(4, 0);
    const Point upper_left(0, 3);
    const Point upper_right(4, 3);
    const Point isolated_point(6, 1);
    const std::vector<GeneratorSite2> generators{
        {"lower-left", Site::construct_site_2(lower_left)},
        {"lower-right", Site::construct_site_2(lower_right)},
        {"upper-left", Site::construct_site_2(upper_left)},
        {"upper-right", Site::construct_site_2(upper_right)},
        {"isolated", Site::construct_site_2(isolated_point)},
        {"lower-segment",
         Site::construct_site_2(lower_left, lower_right)},
        {"upper-segment",
         Site::construct_site_2(upper_left, upper_right)},
    };

    SegmentSiteDelaunay2 delaunay;
    delaunay.insert(lower_left, lower_right);
    delaunay.insert(upper_left, upper_right);
    delaunay.insert(isolated_point);

    SegmentSiteVoronoi2 voronoi(delaunay);
    static_cast<void>(voronoi);

    std::size_t assigned = 0;
    for (auto edge = delaunay.finite_edges_begin();
         edge != delaunay.finite_edges_end();
         ++edge) {
        assigned += assign_dual_primitive(
            delaunay.primal(edge));
    }
    const std::size_t delaunay_vertices = static_cast<std::size_t>(
        std::distance(
            delaunay.finite_vertices_begin(),
            delaunay.finite_vertices_end()));
    return {
        delaunay.is_valid(),
        assigned,
        exact_clearance_root_count(),
        matched_generator_site_count(delaunay, generators),
        delaunay_vertices,
    };
}
