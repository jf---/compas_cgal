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
    const RationalPrimitiveParameterization2 clipped_line{
        line.x_coefficients,
        line.y_coefficients,
        CORE::BigRat(-2),
        CORE::BigRat(2),
    };
    const ClearanceRootBoundary2 clipped_line_boundary =
        point_clearance_boundary(clipped_line, 0, 0, 1);

    const std::vector<MatAdmissibleComponent2> line_components =
        maximal_clearance_components(
            "line-dual",
            line,
            line_boundary);
    const std::vector<MatAdmissibleComponent2> ray_components =
        maximal_clearance_components(
            "ray-dual",
            ray,
            ray_boundary);
    const std::vector<MatAdmissibleComponent2> segment_components =
        maximal_clearance_components(
            "segment-dual",
            segment,
            segment_boundary);
    const std::vector<MatAdmissibleComponent2> parabola_components =
        maximal_clearance_components(
            "parabola-dual",
            parabola,
            parabola_boundary);
    const std::vector<MatAdmissibleComponent2>
        clipped_line_components =
            maximal_clearance_components(
                "clipped-line-dual",
                clipped_line,
                clipped_line_boundary);
    const std::vector<MatAdmissibleComponent2> zero_components =
        maximal_clearance_components(
            "zero-dual",
            {{CORE::BigRat(0)}, {CORE::BigRat(0)}, 0, 1},
            zero);
    const std::vector<MatAdmissibleComponent2> positive_components =
        maximal_clearance_components(
            "positive-dual",
            {{CORE::BigRat(1)}, {CORE::BigRat(0)}, 0, 1},
            positive);
    const std::vector<MatAdmissibleComponent2> negative_components =
        maximal_clearance_components(
            "negative-dual",
            {{CORE::BigRat(0)}, {CORE::BigRat(0)}, 0, 1},
            negative);

    MatDomainPolygon2 outer;
    outer.push_back({-3, -2});
    outer.push_back({3, -2});
    outer.push_back({3, 2});
    outer.push_back({-3, 2});
    MatDomainPolygon2 hole;
    hole.push_back({-1, -1});
    hole.push_back({-1, 1});
    hole.push_back({1, 1});
    hole.push_back({1, -1});
    const std::vector<MatDomainPolygon2> holes{hole};
    const MatDomainPolygonWithHoles2 domain(
        outer,
        holes.begin(),
        holes.end());

    const std::vector<MatAdmissibleComponent2>
        clipped_domain_line =
            clip_linear_clearance_components(
                "D-line-dual",
                line,
                line_boundary,
                domain);
    const std::vector<MatAdmissibleComponent2>
        clipped_domain_ray =
            clip_linear_clearance_components(
                "D-ray-dual",
                ray,
                ray_boundary,
                domain);
    const RationalPrimitiveParameterization2
        crossing_segment{
            {CORE::BigRat(-4), CORE::BigRat(8)},
            {CORE::BigRat(0)},
            CORE::BigRat(0),
            CORE::BigRat(1),
        };
    const ClearanceRootBoundary2
        crossing_segment_boundary =
            point_clearance_boundary(
                crossing_segment,
                0,
                10,
                0);
    const std::vector<MatAdmissibleComponent2>
        clipped_domain_segment =
            clip_linear_clearance_components(
                "D-segment-dual",
                crossing_segment,
                crossing_segment_boundary,
                domain);
    const RationalPrimitiveParameterization2
        crossing_parabola{
            {CORE::BigRat(0), CORE::BigRat(1)},
            {
                CORE::BigRat(0),
                CORE::BigRat(0),
                CORE::BigRat(1),
            },
            CORE::BigRat(-2),
            CORE::BigRat(2),
        };
    const ClearanceRootBoundary2
        crossing_parabola_boundary =
            point_clearance_boundary(
                crossing_parabola,
                0,
                0,
                0);
    const std::vector<MatAdmissibleComponent2>
        clipped_domain_parabola =
            clip_parabola_clearance_components(
                "D-parabola-dual",
                crossing_parabola,
                crossing_parabola_boundary,
                domain);
    MatDomainPolygon2 source_outer;
    source_outer.push_back({-4, -1});
    source_outer.push_back({4, -1});
    source_outer.push_back({4, 2});
    source_outer.push_back({-4, 2});
    const MatDomainPolygonWithHoles2 source_domain(
        source_outer,
        holes.begin(),
        holes.end());
    const std::vector<MatAdmissibleComponent2>
        source_parabola_components =
            clip_source_parabola_clearance_components(
                "source-parabola-dual",
                {
                    "irrational-focus",
                    {0, 1},
                    {1, 0},
                    2,
                },
                {
                    "open-directrix",
                    0,
                    1,
                    1,
                },
                CORE::BigRat(-4),
                CORE::BigRat(4),
                {
                    CGAL::POSITIVE,
                    {},
                    {},
                },
                source_domain);

    return line_boundary.roots.size() == 2
        && ray_boundary.roots.size() == 1
        && segment_boundary.roots.size() == 1
        && parabola_boundary.roots.size() == 2
        && zero.constant_sign == CGAL::ZERO
        && positive.constant_sign == CGAL::POSITIVE
        && negative.constant_sign == CGAL::NEGATIVE
        && line_components.size() == 2
        && !line_components.front().lower.parameter.has_value()
        && !line_components.back().upper.parameter.has_value()
        && ray_components.size() == 1
        && segment_components.size() == 1
        && parabola_components.size() == 2
        && clipped_line_components.size() == 2
        && clipped_line_components.front().lower.parameter.has_value()
        && clipped_line_components.back().upper.parameter.has_value()
        && clipped_line_components.front().component_id
            == "clipped-line-dual/component-0"
        && clipped_line_components.front()
            .upper.provenance_ids.front()
            == "clipped-line-dual/clearance-root-0"
        && clipped_line_components.back()
            .lower.provenance_ids.front()
            == "clipped-line-dual/clearance-root-1"
        && zero_components.size() == 1
        && positive_components.size() == 1
        && negative_components.empty()
        && clipped_domain_line.size() == 2
        && clipped_domain_line.front().lower.parameter.has_value()
        && clipped_domain_line.front().upper.parameter.has_value()
        && clipped_domain_line.back().lower.parameter.has_value()
        && clipped_domain_line.back().upper.parameter.has_value()
        && clipped_domain_ray.size() == 1
        && clipped_domain_segment.size() == 2
        && clipped_domain_segment.front()
            .lower.provenance_ids.front()
            == "D-segment-dual/D-outer/edge-3"
        && clipped_domain_segment.back()
            .upper.provenance_ids.front()
            == "D-segment-dual/D-outer/edge-1"
        && clipped_domain_parabola.size() == 2
        && clipped_domain_parabola.front()
            .lower.provenance_ids.front()
            == "D-parabola-dual/D-outer/edge-2"
               "/algebraic-solution-0"
        && clipped_domain_parabola.back()
            .upper.provenance_ids.front()
            == "D-parabola-dual/D-outer/edge-2"
               "/algebraic-solution-1"
        && source_parabola_components.size() == 2
        && source_parabola_components.front()
            .lower.provenance_ids.front()
            .find("source-algebraic-solution")
            != std::string::npos
        && source_parabola_components.back()
            .upper.provenance_ids.front()
            .find("source-algebraic-solution")
            != std::string::npos;
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
