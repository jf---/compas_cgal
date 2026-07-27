#include "segment_site_mat.h"

#include <algorithm>
#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <CGAL/Object.h>

namespace {

bool same_unoriented_segment(
    const MatTraits::Segment_2& lhs,
    const MatTraits::Segment_2& rhs)
{
    return (lhs.source() == rhs.source()
            && lhs.target() == rhs.target())
        || (lhs.source() == rhs.target()
            && lhs.target() == rhs.source());
}

bool nonparallel_segment_segment_producer_contract()
{
    using Point = MatTraits::Point_2;
    using Site = MatTraits::Site_2;
    const Point lower_left(-20, 0);
    const Point lower_right(8, 0);
    const Point diagonal_lower(5, -4);
    const Point diagonal_upper(20, 11);
    const std::vector<GeneratorSite2> generators{
        {"lower-left", Site::construct_site_2(lower_left)},
        {"lower-right", Site::construct_site_2(lower_right)},
        {
            "diagonal-lower",
            Site::construct_site_2(diagonal_lower),
        },
        {
            "diagonal-upper",
            Site::construct_site_2(diagonal_upper),
        },
        {
            "lower-segment",
            Site::construct_site_2(lower_left, lower_right),
        },
        {
            "diagonal-segment",
            Site::construct_site_2(
                diagonal_lower,
                diagonal_upper),
        },
    };
    SegmentSiteDelaunay2 delaunay;
    delaunay.insert(lower_left, lower_right);
    delaunay.insert(diagonal_lower, diagonal_upper);
    if (require_generator_site_bijection(
            delaunay,
            generators)
        != 6) {
        return false;
    }

    const CORE::Expr sqrt_two =
        CORE::sqrt(CORE::Expr(2));
    const Point upper_transition(
        8,
        CORE::Expr(1) + sqrt_two);
    const Point lower_transition(
        8,
        CORE::Expr(1) - sqrt_two);
    const Point upper_feature(
        CORE::Expr(9)
            - CORE::Expr(11) * sqrt_two,
        CORE::Expr(22)
            + CORE::Expr(11) * sqrt_two);
    const Point lower_feature(
        CORE::Expr(9)
            - CORE::Expr(4) * sqrt_two,
        CORE::Expr(-8)
            + CORE::Expr(4) * sqrt_two);
    const MatTraits::Segment_2 expected_upper(
        upper_feature,
        upper_transition);
    const MatTraits::Segment_2 expected_lower(
        lower_transition,
        lower_feature);
    const MatTraits::Segment_2 transition_span(
        lower_transition,
        upper_transition);

    SegmentSiteVoronoi2 voronoi(delaunay);
    std::size_t raw_segment_segment_edges = 0;
    std::size_t rejected_transition_edges = 0;
    std::size_t rejected_transition_segments = 0;
    std::size_t rejected_transition_parabolas = 0;
    bool distinct_parent_transition_is_parabola = false;
    bool self_transition_is_segment = false;
    for (auto edge = voronoi.dual().finite_edges_begin();
         edge != voronoi.dual().finite_edges_end();
         ++edge) {
        const auto raw = *edge;
        const auto face = raw.first;
        const int index = raw.second;
        const std::string first_id =
            stable_generator_site_id(
                face->vertex(
                        voronoi.dual().ccw(index))
                    ->site(),
                generators);
        const std::string second_id =
            stable_generator_site_id(
                face->vertex(
                        voronoi.dual().cw(index))
                    ->site(),
                generators);
        const CGAL::Object primal =
            voronoi.dual().primal(raw);
        const MatTraits::Segment_2* segment =
            CGAL::object_cast<MatTraits::Segment_2>(
                &primal);
        const SegmentSiteParabola2* parabola =
            CGAL::object_cast<SegmentSiteParabola2>(
                &primal);
        const std::vector<std::string> pair_ids =
            ordered_generator_site_ids(
                first_id,
                second_id);
        const bool rejected =
            voronoi.edge_rejector()(
                voronoi.dual(),
                raw);
        if (pair_ids
            == std::vector<std::string>{
                "diagonal-segment",
                "lower-segment",
            }) {
            if (rejected || segment == nullptr
                || (!same_unoriented_segment(
                         *segment,
                         expected_upper)
                    && !same_unoriented_segment(
                        *segment,
                        expected_lower))) {
                return false;
            }
            ++raw_segment_segment_edges;
            continue;
        }
        if (!rejected) {
            continue;
        }
        const bool expected_transition_pair =
            pair_ids
                == std::vector<std::string>{
                    "diagonal-segment",
                    "lower-right",
                }
            || pair_ids
                == std::vector<std::string>{
                    "lower-right",
                    "lower-segment",
                };
        const bool segment_matches =
            segment != nullptr
            && same_unoriented_segment(
                *segment,
                transition_span);
        const bool parabola_matches =
            parabola != nullptr
            && ((parabola->p1 == lower_transition
                 && parabola->p2 == upper_transition)
                || (parabola->p1 == upper_transition
                    && parabola->p2
                        == lower_transition));
        if (!expected_transition_pair
            || (!segment_matches
                && !parabola_matches)) {
            return false;
        }
        if (pair_ids
            == std::vector<std::string>{
                "diagonal-segment",
                "lower-right",
            }) {
            if (!parabola_matches || segment_matches) {
                return false;
            }
            distinct_parent_transition_is_parabola = true;
        } else {
            if (!segment_matches || parabola_matches) {
                return false;
            }
            self_transition_is_segment = true;
        }
        ++rejected_transition_edges;
        rejected_transition_segments += segment_matches;
        rejected_transition_parabolas += parabola_matches;
    }

    std::size_t normalized_segment_segment_edges = 0;
    std::size_t matching_normalized_primitives = 0;
    std::size_t composite_normalized_primitives = 0;
    for (auto halfedge = voronoi.halfedges_begin();
         halfedge != voronoi.halfedges_end();
         ++halfedge) {
        const std::string up_id =
            stable_generator_site_id(
                halfedge->up()->site(),
                generators);
        const std::string down_id =
            stable_generator_site_id(
                halfedge->down()->site(),
                generators);
        if (ordered_generator_site_ids(up_id, down_id)
            != std::vector<std::string>{
                "diagonal-segment",
                "lower-segment",
            }
            || up_id != "diagonal-segment") {
            continue;
        }
        if (!halfedge->has_source()
            || !halfedge->has_target()) {
            return false;
        }
        MatTraits::Segment_2 primal_segment;
        if (!CGAL::assign(
                primal_segment,
                voronoi.dual().primal(
                    halfedge->dual()))) {
            return false;
        }
        const MatTraits::Segment_2 adaptor_segment(
            halfedge->source()->point(),
            halfedge->target()->point());
        const std::string left_id =
            stable_generator_site_id(
                halfedge->left()->site(),
                generators);
        const std::string right_id =
            stable_generator_site_id(
                halfedge->right()->site(),
                generators);
        if (same_unoriented_segment(
                primal_segment,
                expected_upper)) {
            if (!same_unoriented_segment(
                    adaptor_segment,
                    expected_upper)
                || left_id != "diagonal-upper"
                || right_id != "lower-right") {
                return false;
            }
            ++matching_normalized_primitives;
        } else if (same_unoriented_segment(
                       primal_segment,
                       expected_lower)) {
            if (!same_unoriented_segment(
                    adaptor_segment,
                    MatTraits::Segment_2(
                        lower_feature,
                        upper_transition))
                || left_id != "lower-right"
                || right_id != "diagonal-lower") {
                return false;
            }
            ++composite_normalized_primitives;
        } else {
            return false;
        }
        ++normalized_segment_segment_edges;
    }
    return raw_segment_segment_edges == 2
        && rejected_transition_edges == 2
        && rejected_transition_segments == 1
        && rejected_transition_parabolas == 1
        && distinct_parent_transition_is_parabola
        && self_transition_is_segment
        && normalized_segment_segment_edges == 2
        && matching_normalized_primitives == 1
        && composite_normalized_primitives == 1;
}

bool nonparallel_charts_equal(
    const NonparallelSegmentBisectorParameterization2& lhs,
    const NonparallelSegmentBisectorParameterization2& rhs)
{
    return lhs.first_segment_id == rhs.first_segment_id
        && lhs.second_segment_id == rhs.second_segment_id
        && lhs.branch_sign == rhs.branch_sign
        && lhs.x_rational == rhs.x_rational
        && lhs.x_radical == rhs.x_radical
        && lhs.y_rational == rhs.y_rational
        && lhs.y_radical == rhs.y_radical
        && lhs.radicand == rhs.radicand;
}

bool nonparallel_segment_charts_are_exact()
{
    const MatExactOpenSegmentSource2 diagonal =
        canonical_open_segment_source(
            "diagonal-segment",
            5,
            -4,
            20,
            11);
    const MatExactOpenSegmentSource2 lower =
        canonical_open_segment_source(
            "lower-segment",
            -20,
            0,
            8,
            0);
    const CORE::Expr sqrt_two =
        CORE::sqrt(CORE::Expr(2));
    const MatTraits::Point_2 upper_feature(
        CORE::Expr(9)
            - CORE::Expr(11) * sqrt_two,
        CORE::Expr(22)
            + CORE::Expr(11) * sqrt_two);
    const MatTraits::Point_2 upper_transition(
        8,
        CORE::Expr(1) + sqrt_two);
    const MatTraits::Point_2 lower_transition(
        8,
        CORE::Expr(1) - sqrt_two);
    const MatTraits::Point_2 lower_feature(
        CORE::Expr(9)
            - CORE::Expr(4) * sqrt_two,
        CORE::Expr(-8)
            + CORE::Expr(4) * sqrt_two);
    const NonparallelSegmentBisectorParameterization2
        upper =
            nonparallel_segment_bisector_parameterization(
                lower,
                diagonal,
                {
                    upper_feature,
                    upper_transition,
                });
    const NonparallelSegmentBisectorParameterization2
        upper_reversed =
            nonparallel_segment_bisector_parameterization(
                canonical_open_segment_source(
                    "diagonal-segment",
                    20,
                    11,
                    5,
                    -4),
                canonical_open_segment_source(
                    "lower-segment",
                    8,
                    0,
                    -20,
                    0),
                {
                    upper_transition,
                    upper_feature,
                });
    const NonparallelSegmentBisectorParameterization2
        lower_branch =
            nonparallel_segment_bisector_parameterization(
                diagonal,
                lower,
                {
                    lower_transition,
                    lower_feature,
                });
    if (upper.first_segment_id != "diagonal-segment"
        || upper.second_segment_id != "lower-segment"
        || upper.branch_sign != 1
        || upper.radicand != 2
        || upper.x_rational
            != std::vector<CORE::BigRat>{9, -1}
        || upper.x_radical
            != std::vector<CORE::BigRat>{0, 1}
        || upper.y_rational
            != std::vector<CORE::BigRat>{0, -1}
        || upper.y_radical
            != std::vector<CORE::BigRat>{0, 0}
        || !nonparallel_charts_equal(
            upper,
            upper_reversed)
        || lower_branch.branch_sign != -1
        || lower_branch.radicand != 2
        || lower_branch.x_rational
            != std::vector<CORE::BigRat>{9, 1}
        || lower_branch.x_radical
            != std::vector<CORE::BigRat>{0, 1}
        || lower_branch.y_rational
            != std::vector<CORE::BigRat>{0, 1}
        || lower_branch.y_radical
            != std::vector<CORE::BigRat>{0, 0}) {
        return false;
    }

    const MatExactOpenSegmentSource2 horizontal =
        canonical_open_segment_source(
            "horizontal",
            -4,
            0,
            4,
            0);
    const MatExactOpenSegmentSource2 vertical =
        canonical_open_segment_source(
            "vertical",
            5,
            -4,
            5,
            4);
    const NonparallelSegmentBisectorParameterization2
        rational =
            nonparallel_segment_bisector_parameterization(
                horizontal,
                vertical,
                {
                    {1, 4},
                    {4, 1},
                });
    return rational.branch_sign == 1
        && rational.radicand == 1
        && rational.x_rational
            == std::vector<CORE::BigRat>{5, 1}
        && rational.x_radical
            == std::vector<CORE::BigRat>{0, 0}
        && rational.y_rational
            == std::vector<CORE::BigRat>{0, -1}
        && rational.y_radical
            == std::vector<CORE::BigRat>{0, 0};
}

bool nonparallel_feature_parameters_are_exact()
{
    const MatExactOpenSegmentSource2 diagonal =
        canonical_open_segment_source(
            "diagonal-segment",
            5,
            -4,
            20,
            11);
    const MatExactOpenSegmentSource2 lower =
        canonical_open_segment_source(
            "lower-segment",
            -20,
            0,
            8,
            0);
    const CORE::Expr sqrt_two =
        CORE::sqrt(CORE::Expr(2));
    const MatTraits::Point_2 upper_feature(
        CORE::Expr(9)
            - CORE::Expr(11) * sqrt_two,
        CORE::Expr(22)
            + CORE::Expr(11) * sqrt_two);
    const MatTraits::Point_2 upper_transition(
        8,
        CORE::Expr(1) + sqrt_two);
    const MatTraits::Point_2 lower_transition(
        8,
        CORE::Expr(1) - sqrt_two);
    const MatTraits::Point_2 lower_feature(
        CORE::Expr(9)
            - CORE::Expr(4) * sqrt_two,
        CORE::Expr(-8)
            + CORE::Expr(4) * sqrt_two);
    const NonparallelSegmentBisectorParameterization2 upper =
        nonparallel_segment_bisector_parameterization(
            diagonal,
            lower,
            {
                upper_feature,
                upper_transition,
            });
    const NonparallelSegmentBisectorParameterization2 lower_branch =
        nonparallel_segment_bisector_parameterization(
            diagonal,
            lower,
            {
                lower_transition,
                lower_feature,
            });
    const auto equals =
        [](const MatQuadraticFieldValue2& value,
           const CORE::BigRat& rational,
           const CORE::BigRat& radical) {
            return value.rational == rational
                && value.radical == radical;
        };
    if (!equals(
            nonparallel_segment_tangent_parameter(
                upper,
                diagonal,
                20,
                11),
            -22,
            -11)
        || !equals(
            nonparallel_segment_tangent_parameter(
                upper,
                diagonal,
                5,
                -4),
            8,
            4)
        || !equals(
            nonparallel_segment_tangent_parameter(
                upper,
                lower,
                -20,
                0),
            -29,
            -29)
        || !equals(
            nonparallel_segment_tangent_parameter(
                upper,
                lower,
                8,
                0),
            -1,
            -1)
        || !equals(
            nonparallel_segment_tangent_parameter(
                lower_branch,
                diagonal,
                5,
                -4),
            -8,
            4)
        || !equals(
            nonparallel_segment_tangent_parameter(
                lower_branch,
                diagonal,
                20,
                11),
            22,
            -11)
        || !equals(
            nonparallel_segment_tangent_parameter(
                lower_branch,
                lower,
                -20,
                0),
            29,
            -29)
        || !equals(
            nonparallel_segment_tangent_parameter(
                lower_branch,
                lower,
                8,
                0),
            1,
            -1)) {
        return false;
    }

    const MatExactOpenSegmentSource2 horizontal =
        canonical_open_segment_source(
            "horizontal",
            -4,
            0,
            4,
            0);
    const MatExactOpenSegmentSource2 vertical =
        canonical_open_segment_source(
            "vertical",
            5,
            -4,
            5,
            4);
    const NonparallelSegmentBisectorParameterization2 rational =
        nonparallel_segment_bisector_parameterization(
            horizontal,
            vertical,
            {
                {1, 4},
                {4, 1},
            });
    return equals(
               nonparallel_segment_tangent_parameter(
                   rational,
                   horizontal,
                   -4,
                   0),
               -9,
               0)
        && equals(
               nonparallel_segment_tangent_parameter(
                   rational,
                   horizontal,
                   4,
                   0),
               -1,
               0)
        && equals(
               nonparallel_segment_tangent_parameter(
                   rational,
                   vertical,
                   5,
                   -4),
               4,
               0)
        && equals(
               nonparallel_segment_tangent_parameter(
                   rational,
                   vertical,
                   5,
                   4),
               -4,
               0);
}

bool nonparallel_feature_domains_are_exact()
{
    const NonparallelSegmentFeatureDomain2 upper =
        intersect_nonparallel_segment_feature_domains(
            {
                {8, 4},
                {"diagonal-lower"},
            },
            {
                {-22, -11},
                {"diagonal-upper"},
            },
            {
                {-29, -29},
                {"lower-left"},
            },
            {
                {-1, -1},
                {"lower-right"},
            },
            2);
    const NonparallelSegmentFeatureDomain2 upper_reversed =
        intersect_nonparallel_segment_feature_domains(
            {
                {-22, -11},
                {"diagonal-upper"},
            },
            {
                {8, 4},
                {"diagonal-lower"},
            },
            {
                {-1, -1},
                {"lower-right"},
            },
            {
                {-29, -29},
                {"lower-left"},
            },
            2);
    const NonparallelSegmentFeatureDomain2 lower =
        intersect_nonparallel_segment_feature_domains(
            {
                {-8, 4},
                {"diagonal-lower"},
            },
            {
                {22, -11},
                {"diagonal-upper"},
            },
            {
                {29, -29},
                {"lower-left"},
            },
            {
                {1, -1},
                {"lower-right"},
            },
            2);
    const auto same_boundary =
        [](const MatQuadraticFieldDomainBoundary2& lhs,
           const MatQuadraticFieldDomainBoundary2& rhs) {
            return lhs.parameter.rational
                    == rhs.parameter.rational
                && lhs.parameter.radical
                    == rhs.parameter.radical
                && lhs.provenance_ids
                    == rhs.provenance_ids;
        };
    if (upper.radicand != 2
        || upper.lower.parameter.rational != -22
        || upper.lower.parameter.radical != -11
        || upper.lower.provenance_ids
            != std::vector<std::string>{
                "diagonal-upper",
            }
        || upper.upper.parameter.rational != -1
        || upper.upper.parameter.radical != -1
        || upper.upper.provenance_ids
            != std::vector<std::string>{
                "lower-right",
            }
        || upper_reversed.radicand != upper.radicand
        || !same_boundary(
            upper_reversed.lower,
            upper.lower)
        || !same_boundary(
            upper_reversed.upper,
            upper.upper)
        || lower.lower.parameter.rational != -8
        || lower.lower.parameter.radical != 4
        || lower.lower.provenance_ids
            != std::vector<std::string>{
                "diagonal-lower",
            }
        || lower.upper.parameter.rational != 1
        || lower.upper.parameter.radical != -1
        || lower.upper.provenance_ids
            != std::vector<std::string>{
                "lower-right",
            }) {
        return false;
    }

    const NonparallelSegmentFeatureDomain2 tied =
        intersect_nonparallel_segment_feature_domains(
            {
                {0, 0},
                {"first-lower"},
            },
            {
                {2, 0},
                {"first-upper"},
            },
            {
                {0, 0},
                {"second-lower"},
            },
            {
                {3, 0},
                {"second-upper"},
            },
            2);
    return tied.lower.parameter.rational == 0
        && tied.lower.parameter.radical == 0
        && tied.lower.provenance_ids
            == std::vector<std::string>{
                "first-lower",
                "second-lower",
            }
        && tied.upper.parameter.rational == 2
        && tied.upper.parameter.radical == 0
        && tied.upper.provenance_ids
            == std::vector<std::string>{
                "first-upper",
            };
}

bool quadratic_feature_parameters_order_and_embed()
{
    const CORE::BigRat radicand = 2;
    const MatQuadraticFieldValue2 upper_lower{
        -22,
        -11,
    };
    const MatQuadraticFieldValue2 upper_upper{
        -1,
        -1,
    };
    const MatQuadraticFieldValue2 lower_lower{
        -8,
        4,
    };
    const MatQuadraticFieldValue2 lower_upper{
        1,
        -1,
    };
    if (quadratic_field_compare(
            upper_lower,
            upper_upper,
            radicand)
            != CGAL::SMALLER
        || quadratic_field_compare(
               lower_lower,
               lower_upper,
               radicand)
            != CGAL::SMALLER
        || quadratic_field_compare(
               upper_lower,
               upper_lower,
               radicand)
            != CGAL::EQUAL) {
        return false;
    }

    ExactAlgebraicKernel1 kernel;
    const auto compare = kernel.compare_1_object();
    const auto construct =
        kernel.construct_algebraic_real_1_object();
    const auto upper_lower_algebraic =
        quadratic_field_algebraic_real(
            upper_lower,
            radicand);
    const auto upper_upper_algebraic =
        quadratic_field_algebraic_real(
            upper_upper,
            radicand);
    const auto lower_lower_algebraic =
        quadratic_field_algebraic_real(
            lower_lower,
            radicand);
    const auto lower_upper_algebraic =
        quadratic_field_algebraic_real(
            lower_upper,
            radicand);
    return compare(
               construct(-38),
               upper_lower_algebraic)
            == CGAL::SMALLER
        && compare(
               upper_lower_algebraic,
               construct(-37))
            == CGAL::SMALLER
        && compare(
               construct(-3),
               upper_upper_algebraic)
            == CGAL::SMALLER
        && compare(
               upper_upper_algebraic,
               construct(-2))
            == CGAL::SMALLER
        && compare(
               upper_lower_algebraic,
               upper_upper_algebraic)
            == CGAL::SMALLER
        && compare(
               lower_lower_algebraic,
               lower_upper_algebraic)
            == CGAL::SMALLER
        && algebraic_root_identity_v1(
               upper_lower_algebraic)
            == algebraic_root_id_v1(
                {242, 44, 1},
                0)
        && algebraic_root_identity_v1(
               upper_upper_algebraic)
            == algebraic_root_id_v1(
                {-1, 2, 1},
                0)
        && algebraic_root_identity_v1(
               lower_lower_algebraic)
            == algebraic_root_id_v1(
                {32, 16, 1},
                1)
        && algebraic_root_identity_v1(
               lower_upper_algebraic)
            == algebraic_root_id_v1(
                {-1, -2, 1},
                0)
        && compare(
               quadratic_field_algebraic_real(
                   {-9, 0},
                   1),
               construct(-9))
            == CGAL::EQUAL;
}

bool nonparallel_clearance_is_rational_quadratic()
{
    const MatExactOpenSegmentSource2 diagonal =
        canonical_open_segment_source(
            "diagonal-segment",
            5,
            -4,
            20,
            11);
    const MatExactOpenSegmentSource2 lower =
        canonical_open_segment_source(
            "lower-segment",
            -20,
            0,
            8,
            0);
    const CORE::Expr sqrt_two =
        CORE::sqrt(CORE::Expr(2));
    const MatTraits::Point_2 upper_feature(
        CORE::Expr(9)
            - CORE::Expr(11) * sqrt_two,
        CORE::Expr(22)
            + CORE::Expr(11) * sqrt_two);
    const MatTraits::Point_2 upper_transition(
        8,
        CORE::Expr(1) + sqrt_two);
    const MatTraits::Point_2 lower_transition(
        8,
        CORE::Expr(1) - sqrt_two);
    const MatTraits::Point_2 lower_feature(
        CORE::Expr(9)
            - CORE::Expr(4) * sqrt_two,
        CORE::Expr(-8)
            + CORE::Expr(4) * sqrt_two);
    const NonparallelSegmentBisectorParameterization2 upper =
        nonparallel_segment_bisector_parameterization(
            diagonal,
            lower,
            {
                upper_feature,
                upper_transition,
            });
    const NonparallelSegmentBisectorParameterization2 lower_branch =
        nonparallel_segment_bisector_parameterization(
            diagonal,
            lower,
            {
                lower_transition,
                lower_feature,
            });
    const ClearanceRootBoundary2 upper_boundary =
        nonparallel_segment_clearance_boundary(
            upper,
            diagonal,
            lower,
            4);
    const ClearanceRootBoundary2 lower_boundary =
        nonparallel_segment_clearance_boundary(
            lower_branch,
            lower,
            diagonal,
            4);
    ExactAlgebraicKernel1 kernel;
    const auto compare = kernel.compare_1_object();
    const auto construct =
        kernel.construct_algebraic_real_1_object();
    if (upper_boundary.constant_sign.has_value()
        || lower_boundary.constant_sign.has_value()
        || upper_boundary.primitive_coefficients
            != std::vector<ExactAlgebraicInteger1>{
                -4,
                0,
                1,
            }
        || lower_boundary.primitive_coefficients
            != upper_boundary.primitive_coefficients
        || upper_boundary.roots.size() != 2
        || lower_boundary.roots.size() != 2) {
        return false;
    }
    for (std::size_t index = 0; index < 2; ++index) {
        if (compare(
                upper_boundary.roots[index].parameter,
                lower_boundary.roots[index].parameter)
                != CGAL::EQUAL
            || upper_boundary.roots[index].event_id
                != algebraic_root_id_v1(
                    {-4, 0, 1},
                    index)
            || lower_boundary.roots[index].event_id
                != upper_boundary.roots[index].event_id) {
            return false;
        }
    }
    if (compare(
            upper_boundary.roots[0].parameter,
            construct(-2))
            != CGAL::EQUAL
        || compare(
               upper_boundary.roots[1].parameter,
               construct(2))
            != CGAL::EQUAL) {
        return false;
    }

    const MatExactOpenSegmentSource2 steep =
        canonical_open_segment_source(
            "first",
            -1,
            2,
            1,
            -2);
    const MatExactOpenSegmentSource2 crossing =
        canonical_open_segment_source(
            "second",
            -2,
            -1,
            2,
            1);
    const NonparallelSegmentBisectorParameterization2 scaled =
        nonparallel_segment_bisector_parameterization(
            steep,
            crossing,
            {
                {0, 0},
                {1, 3},
            });
    const ClearanceRootBoundary2 scaled_boundary =
        nonparallel_segment_clearance_boundary(
            scaled,
            steep,
            crossing,
            5);
    return !scaled_boundary.constant_sign.has_value()
        && scaled_boundary.primitive_coefficients
            == std::vector<ExactAlgebraicInteger1>{
                -1,
                0,
                25,
            }
        && scaled_boundary.roots.size() == 2
        && compare(
               scaled_boundary.roots[0].parameter,
               construct(CORE::BigRat(-1, 5)))
            == CGAL::EQUAL
        && compare(
               scaled_boundary.roots[1].parameter,
               construct(CORE::BigRat(1, 5)))
            == CGAL::EQUAL
        && scaled_boundary.roots[0].event_id
            == algebraic_root_id_v1(
                {-1, 0, 25},
                0)
        && scaled_boundary.roots[1].event_id
            == algebraic_root_id_v1(
                {-1, 0, 25},
                1);
}

bool nonparallel_feature_clearance_components_are_exact()
{
    const MatExactOpenSegmentSource2 diagonal =
        canonical_open_segment_source(
            "diagonal-segment",
            5,
            -4,
            20,
            11);
    const MatExactOpenSegmentSource2 lower =
        canonical_open_segment_source(
            "lower-segment",
            -20,
            0,
            8,
            0);
    const CORE::Expr sqrt_two =
        CORE::sqrt(CORE::Expr(2));
    const MatTraits::Point_2 upper_feature(
        CORE::Expr(9)
            - CORE::Expr(11) * sqrt_two,
        CORE::Expr(22)
            + CORE::Expr(11) * sqrt_two);
    const MatTraits::Point_2 upper_transition(
        8,
        CORE::Expr(1) + sqrt_two);
    const MatTraits::Point_2 lower_transition(
        8,
        CORE::Expr(1) - sqrt_two);
    const MatTraits::Point_2 lower_feature(
        CORE::Expr(9)
            - CORE::Expr(4) * sqrt_two,
        CORE::Expr(-8)
            + CORE::Expr(4) * sqrt_two);
    const NonparallelSegmentBisectorParameterization2 upper =
        nonparallel_segment_bisector_parameterization(
            diagonal,
            lower,
            {
                upper_feature,
                upper_transition,
            });
    const NonparallelSegmentBisectorParameterization2 lower_branch =
        nonparallel_segment_bisector_parameterization(
            diagonal,
            lower,
            {
                lower_transition,
                lower_feature,
            });
    const NonparallelSegmentFeatureDomain2 upper_feature_domain =
        intersect_nonparallel_segment_feature_domains(
            {
                {8, 4},
                {"diagonal-lower"},
            },
            {
                {-22, -11},
                {"diagonal-upper"},
            },
            {
                {-29, -29},
                {"lower-left"},
            },
            {
                {-1, -1},
                {"lower-right"},
            },
            2);
    const NonparallelSegmentFeatureDomain2 lower_feature_domain =
        intersect_nonparallel_segment_feature_domains(
            {
                {-8, 4},
                {"diagonal-lower"},
            },
            {
                {22, -11},
                {"diagonal-upper"},
            },
            {
                {29, -29},
                {"lower-left"},
            },
            {
                {1, -1},
                {"lower-right"},
            },
            2);
    const std::vector<MatAdmissibleComponent2>
        upper_radius_three =
            maximal_nonparallel_segment_clearance_components(
                "nonparallel-upper",
                upper_feature_domain,
                nonparallel_segment_clearance_boundary(
                    upper,
                    diagonal,
                    lower,
                    9));
    const std::vector<MatAdmissibleComponent2>
        upper_radius_two =
            maximal_nonparallel_segment_clearance_components(
                "nonparallel-upper",
                upper_feature_domain,
                nonparallel_segment_clearance_boundary(
                    upper,
                    diagonal,
                    lower,
                    4));
    const std::vector<MatAdmissibleComponent2>
        upper_radius_forty =
            maximal_nonparallel_segment_clearance_components(
                "nonparallel-upper",
                upper_feature_domain,
                nonparallel_segment_clearance_boundary(
                    upper,
                    diagonal,
                    lower,
                    1600));
    const std::vector<MatAdmissibleComponent2>
        lower_radius_two =
            maximal_nonparallel_segment_clearance_components(
                "nonparallel-lower",
                lower_feature_domain,
                nonparallel_segment_clearance_boundary(
                    lower_branch,
                    diagonal,
                    lower,
                    4));
    ExactAlgebraicKernel1 kernel;
    const auto compare = kernel.compare_1_object();
    const auto construct =
        kernel.construct_algebraic_real_1_object();
    const auto endpoint_is =
        [&compare](
            const MatParameterEndpoint2& endpoint,
            const ExactAlgebraicKernel1::Algebraic_real_1&
                parameter,
            const std::vector<std::string>& provenance) {
            return endpoint.parameter.has_value()
                && compare(
                       *endpoint.parameter,
                       parameter)
                    == CGAL::EQUAL
                && endpoint.provenance_ids == provenance;
        };
    if (upper_radius_three.size() != 1
        || upper_radius_three[0].component_id
            != "nonparallel-upper/component-0"
        || !endpoint_is(
            upper_radius_three[0].lower,
            quadratic_field_algebraic_real(
                {-22, -11},
                2),
            {"diagonal-upper"})
        || !endpoint_is(
            upper_radius_three[0].upper,
            construct(-3),
            {
                algebraic_root_id_v1(
                    {-9, 0, 1},
                    0),
            })
        || upper_radius_two.size() != 1
        || !endpoint_is(
            upper_radius_two[0].lower,
            quadratic_field_algebraic_real(
                {-22, -11},
                2),
            {"diagonal-upper"})
        || !endpoint_is(
            upper_radius_two[0].upper,
            quadratic_field_algebraic_real(
                {-1, -1},
                2),
            {"lower-right"})
        || !upper_radius_forty.empty()) {
        return false;
    }
    return lower_radius_two.size() == 1
        && lower_radius_two[0].component_id
            == "nonparallel-lower/component-0"
        && endpoint_is(
            lower_radius_two[0].lower,
            quadratic_field_algebraic_real(
                {-8, 4},
                2),
            {"diagonal-lower"})
        && endpoint_is(
            lower_radius_two[0].upper,
            construct(-2),
            {
                algebraic_root_id_v1(
                    {-4, 0, 1},
                    0),
            });
}

bool nonparallel_domain_clipping_is_exact()
{
    const MatExactOpenSegmentSource2 diagonal =
        canonical_open_segment_source(
            "diagonal-segment",
            5,
            -4,
            20,
            11);
    const MatExactOpenSegmentSource2 lower =
        canonical_open_segment_source(
            "lower-segment",
            -20,
            0,
            8,
            0);
    const CORE::Expr sqrt_two =
        CORE::sqrt(CORE::Expr(2));
    const MatTraits::Point_2 upper_feature(
        CORE::Expr(9)
            - CORE::Expr(11) * sqrt_two,
        CORE::Expr(22)
            + CORE::Expr(11) * sqrt_two);
    const MatTraits::Point_2 upper_transition(
        8,
        CORE::Expr(1) + sqrt_two);
    const NonparallelSegmentBisectorParameterization2 upper =
        nonparallel_segment_bisector_parameterization(
            diagonal,
            lower,
            {
                upper_feature,
                upper_transition,
            });
    const NonparallelSegmentFeatureDomain2 feature =
        intersect_nonparallel_segment_feature_domains(
            {
                {8, 4},
                {"diagonal-lower"},
            },
            {
                {-22, -11},
                {"diagonal-upper"},
            },
            {
                {-29, -29},
                {"lower-left"},
            },
            {
                {-1, -1},
                {"lower-right"},
            },
            2);
    MatDomainPolygon2 outer;
    outer.push_back({0, 4});
    outer.push_back({10, 4});
    outer.push_back({10, 10});
    outer.push_back({0, 10});
    MatDomainPolygon2 hole;
    hole.push_back({5, 6});
    hole.push_back({5, 8});
    hole.push_back({8, 8});
    hole.push_back({8, 6});
    const std::vector<MatDomainPolygon2> holes{hole};
    const MatDomainPolygonWithHoles2 domain(
        outer,
        holes.begin(),
        holes.end());
    const std::vector<MatAdmissibleComponent2> components =
        clip_nonparallel_segment_clearance_components(
            "nonparallel-D",
            upper,
            feature,
            nonparallel_segment_clearance_boundary(
                upper,
                diagonal,
                lower,
                9),
            domain);
    if (components.size() != 2) {
        return false;
    }
    ExactAlgebraicKernel1 kernel;
    const auto compare = kernel.compare_1_object();
    const auto construct =
        kernel.construct_algebraic_real_1_object();
    const auto endpoint_has =
        [&compare](
            const MatParameterEndpoint2& endpoint,
            const CORE::BigRat& parameter,
            const std::string& provenance) {
            return endpoint.parameter.has_value()
                && compare(
                       *endpoint.parameter,
                       parameter)
                    == CGAL::EQUAL
                && std::find(
                       endpoint.provenance_ids.begin(),
                       endpoint.provenance_ids.end(),
                       provenance)
                    != endpoint.provenance_ids.end()
                && std::find(
                       endpoint.provenance_ids.begin(),
                       endpoint.provenance_ids.end(),
                       algebraic_root_identity_v1(
                           ExactAlgebraicKernel1()
                               .construct_algebraic_real_1_object()(
                                   parameter)))
                    != endpoint.provenance_ids.end();
        };
    return components[0].component_id
            == "nonparallel-D/component-0"
        && components[1].component_id
            == "nonparallel-D/component-1"
        && endpoint_has(
            components[0].lower,
            -10,
            "nonparallel-D/D-outer/edge-2")
        && endpoint_has(
            components[0].upper,
            -8,
            "nonparallel-D/D-hole-0/edge-1")
        && endpoint_has(
            components[1].lower,
            -6,
            "nonparallel-D/D-hole-0/edge-3")
        && endpoint_has(
            components[1].upper,
            -4,
            "nonparallel-D/D-outer/edge-0")
        && compare(
               *components[0].lower.parameter,
               construct(-10))
            == CGAL::EQUAL;
}

bool collinear_domain_intersections_are_exact()
{
    const MatExactOpenSegmentSource2 horizontal =
        canonical_open_segment_source(
            "horizontal",
            -4,
            0,
            4,
            0);
    const MatExactOpenSegmentSource2 vertical =
        canonical_open_segment_source(
            "vertical",
            5,
            -4,
            5,
            4);
    const NonparallelSegmentBisectorParameterization2 primitive =
        nonparallel_segment_bisector_parameterization(
            horizontal,
            vertical,
            {
                {1, 4},
                {4, 1},
            });
    const NonparallelSegmentFeatureDomain2 feature =
        intersect_nonparallel_segment_feature_domains(
            {
                {-9, 0},
                {"horizontal-left"},
            },
            {
                {-1, 0},
                {"horizontal-right"},
            },
            {
                {4, 0},
                {"vertical-lower"},
            },
            {
                {-4, 0},
                {"vertical-upper"},
            },
            1);
    const ClearanceRootBoundary2 clearance =
        nonparallel_segment_clearance_boundary(
            primitive,
            horizontal,
            vertical,
            0);

    MatDomainPolygon2 disjoint;
    disjoint.push_back({6, -1});
    disjoint.push_back({7, -2});
    disjoint.push_back({8, 0});
    if (!clip_nonparallel_segment_clearance_components(
             "collinear-disjoint",
             primitive,
             feature,
             clearance,
             MatDomainPolygonWithHoles2(disjoint))
             .empty()) {
        return false;
    }

    MatDomainPolygon2 touching;
    touching.push_back({0, 3});
    touching.push_back({1, 4});
    touching.push_back({0, 5});
    const std::vector<MatAdmissibleComponent2> components =
        clip_nonparallel_segment_clearance_components(
            "collinear-touch",
            primitive,
            feature,
            clearance,
            MatDomainPolygonWithHoles2(touching));
    if (components.size() != 1
        || !components[0].lower.parameter.has_value()
        || !components[0].upper.parameter.has_value()) {
        return false;
    }
    ExactAlgebraicKernel1 kernel;
    const auto compare = kernel.compare_1_object();
    const auto construct =
        kernel.construct_algebraic_real_1_object();
    const auto has =
        [](const MatParameterEndpoint2& endpoint,
           const std::string& provenance) {
            return std::find(
                       endpoint.provenance_ids.begin(),
                       endpoint.provenance_ids.end(),
                       provenance)
                != endpoint.provenance_ids.end();
        };
    return compare(
               *components[0].lower.parameter,
               construct(-4))
            == CGAL::EQUAL
        && compare(
               *components[0].upper.parameter,
               construct(-4))
            == CGAL::EQUAL
        && has(
            components[0].lower,
            "vertical-upper")
        && has(
            components[0].lower,
            "collinear-touch/D-outer/edge-1")
        && has(
            components[0].lower,
            algebraic_root_identity_v1(
                construct(-4)));
}

bool unsupported_nonparallel_charts_fail_loudly()
{
    const MatExactOpenSegmentSource2 horizontal =
        canonical_open_segment_source(
            "horizontal",
            -4,
            0,
            4,
            0);
    const MatExactOpenSegmentSource2 parallel =
        canonical_open_segment_source(
            "parallel",
            -4,
            3,
            4,
            3);
    const MatExactOpenSegmentSource2 diagonal =
        canonical_open_segment_source(
            "diagonal",
            5,
            -4,
            20,
            11);
    const CORE::Expr sqrt_two =
        CORE::sqrt(CORE::Expr(2));
    const NonparallelSegmentBisectorParameterization2 valid =
        nonparallel_segment_bisector_parameterization(
            horizontal,
            diagonal,
            {
                {9, 0},
                {
                    8,
                    CORE::Expr(1) + sqrt_two,
                },
            });
    bool parallel_rejected = false;
    bool degenerate_primitive_rejected = false;
    bool unbound_branch_rejected = false;
    bool mismatched_source_rejected = false;
    bool off_support_endpoint_rejected = false;
    bool invalid_radicand_rejected = false;
    bool mismatched_clearance_sources_rejected = false;
    bool nonrational_clearance_rejected = false;
    bool negative_clearance_radius_rejected = false;
    bool empty_feature_domain_rejected = false;
    bool missing_feature_provenance_rejected = false;
    bool mismatched_feature_domain_rejected = false;
    bool overlapping_domain_boundary_rejected = false;
    try {
        static_cast<void>(
            nonparallel_segment_bisector_parameterization(
                horizontal,
                parallel,
                {
                    {0, CORE::BigRat(3, 2)},
                    {1, CORE::BigRat(3, 2)},
                }));
    } catch (const ParallelSegmentSupportsError&) {
        parallel_rejected = true;
    }
    try {
        static_cast<void>(
            nonparallel_segment_bisector_parameterization(
                horizontal,
                diagonal,
                {
                    {9, 0},
                    {9, 0},
                }));
    } catch (const DegenerateLiveSegmentPrimitiveError&) {
        degenerate_primitive_rejected = true;
    }
    try {
        static_cast<void>(
            nonparallel_segment_bisector_parameterization(
                horizontal,
                diagonal,
                {
                    {0, 0},
                    {1, 0},
                }));
    } catch (
        const UnboundNonparallelSegmentBranchError&) {
        unbound_branch_rejected = true;
    }
    try {
        static_cast<void>(
            nonparallel_segment_tangent_parameter(
                valid,
                parallel,
                0,
                3));
    } catch (
        const MismatchedNonparallelSegmentSourceError&) {
        mismatched_source_rejected = true;
    }
    try {
        static_cast<void>(
            nonparallel_segment_tangent_parameter(
                valid,
                horizontal,
                0,
                1));
    } catch (const OffSupportSegmentEndpointError&) {
        off_support_endpoint_rejected = true;
    }
    try {
        static_cast<void>(
            quadratic_field_compare(
                {0, 0},
                {1, 0},
                0));
    } catch (
        const InvalidQuadraticFieldRadicandError&) {
        invalid_radicand_rejected = true;
    }
    try {
        static_cast<void>(
            nonparallel_segment_clearance_boundary(
                valid,
                horizontal,
                parallel,
                1));
    } catch (
        const MismatchedNonparallelSegmentClearanceError&) {
        mismatched_clearance_sources_rejected = true;
    }
    try {
        static_cast<void>(
            nonparallel_segment_clearance_boundary(
                valid,
                horizontal,
                canonical_open_segment_source(
                    "diagonal",
                    5,
                    -4,
                    5,
                    4),
                1));
    } catch (
        const NonrationalNonparallelSegmentClearanceError&) {
        nonrational_clearance_rejected = true;
    }
    try {
        static_cast<void>(
            nonparallel_segment_clearance_boundary(
                valid,
                horizontal,
                diagonal,
                -1));
    } catch (
        const NegativeClearanceRadiusSquaredError&) {
        negative_clearance_radius_rejected = true;
    }
    try {
        static_cast<void>(
            intersect_nonparallel_segment_feature_domains(
                {
                    {0, 0},
                    {"first-lower"},
                },
                {
                    {1, 0},
                    {"first-upper"},
                },
                {
                    {1, 0},
                    {"second-lower"},
                },
                {
                    {2, 0},
                    {"second-upper"},
                },
                2));
    } catch (
        const EmptyNonparallelSegmentFeatureDomainError&) {
        empty_feature_domain_rejected = true;
    }
    try {
        static_cast<void>(
            intersect_nonparallel_segment_feature_domains(
                {
                    {0, 0},
                    {},
                },
                {
                    {2, 0},
                    {"first-upper"},
                },
                {
                    {0, 0},
                    {"second-lower"},
                },
                {
                    {3, 0},
                    {"second-upper"},
                },
                2));
    } catch (
        const MissingNonparallelSegmentFeatureProvenanceError&) {
        missing_feature_provenance_rejected = true;
    }
    try {
        const NonparallelSegmentFeatureDomain2 wrong_field =
            intersect_nonparallel_segment_feature_domains(
                {
                    {0, 0},
                    {"first-lower"},
                },
                {
                    {2, 0},
                    {"first-upper"},
                },
                {
                    {0, 0},
                    {"second-lower"},
                },
                {
                    {3, 0},
                    {"second-upper"},
                },
                1);
        MatDomainPolygon2 outer;
        outer.push_back({-10, -10});
        outer.push_back({10, -10});
        outer.push_back({10, 10});
        outer.push_back({-10, 10});
        static_cast<void>(
            clip_nonparallel_segment_clearance_components(
                "wrong-field",
                valid,
                wrong_field,
                nonparallel_segment_clearance_boundary(
                    valid,
                    horizontal,
                    diagonal,
                    1),
                MatDomainPolygonWithHoles2(outer)));
    } catch (
        const MismatchedNonparallelSegmentFeatureDomainError&) {
        mismatched_feature_domain_rejected = true;
    }
    try {
        const MatExactOpenSegmentSource2 vertical =
            canonical_open_segment_source(
                "vertical",
                5,
                -4,
                5,
                4);
        const NonparallelSegmentBisectorParameterization2 rational =
            nonparallel_segment_bisector_parameterization(
                horizontal,
                vertical,
                {
                    {1, 4},
                    {4, 1},
                });
        const NonparallelSegmentFeatureDomain2 rational_feature =
            intersect_nonparallel_segment_feature_domains(
                {
                    {-9, 0},
                    {"horizontal-left"},
                },
                {
                    {-1, 0},
                    {"horizontal-right"},
                },
                {
                    {4, 0},
                    {"vertical-lower"},
                },
                {
                    {-4, 0},
                    {"vertical-upper"},
                },
                1);
        MatDomainPolygon2 overlapping;
        overlapping.push_back({1, 4});
        overlapping.push_back({4, 1});
        overlapping.push_back({6, 6});
        static_cast<void>(
            clip_nonparallel_segment_clearance_components(
                "overlapping-D",
                rational,
                rational_feature,
                nonparallel_segment_clearance_boundary(
                    rational,
                    horizontal,
                    vertical,
                    0),
                MatDomainPolygonWithHoles2(
                    overlapping)));
    } catch (const OverlappingDomainBoundaryError&) {
        overlapping_domain_boundary_rejected = true;
    }
    return parallel_rejected
        && degenerate_primitive_rejected
        && unbound_branch_rejected
        && mismatched_source_rejected
        && off_support_endpoint_rejected
        && invalid_radicand_rejected
        && mismatched_clearance_sources_rejected
        && nonrational_clearance_rejected
        && negative_clearance_radius_rejected
        && empty_feature_domain_rejected
        && missing_feature_provenance_rejected
        && mismatched_feature_domain_rejected
        && overlapping_domain_boundary_rejected;
}

bool has_provenance(
    const MatParameterEndpoint2& endpoint,
    const std::string& expected)
{
    return std::find(
               endpoint.provenance_ids.begin(),
               endpoint.provenance_ids.end(),
               expected)
        != endpoint.provenance_ids.end();
}

bool endpoints_equal(
    const MatParameterEndpoint2& lhs,
    const MatParameterEndpoint2& rhs)
{
    if (lhs.parameter.has_value()
        != rhs.parameter.has_value()
        || lhs.provenance_ids != rhs.provenance_ids) {
        return false;
    }
    if (!lhs.parameter.has_value()) {
        return true;
    }
    ExactAlgebraicKernel1 kernel;
    return kernel.compare_1_object()(
               *lhs.parameter,
               *rhs.parameter)
        == CGAL::EQUAL;
}

bool graphs_equal(
    const MatExactGraph2& lhs,
    const MatExactGraph2& rhs)
{
    if (lhs.nodes.size() != rhs.nodes.size()
        || lhs.edges.size() != rhs.edges.size()
        || lhs.rejected_incident_transitions
            != rhs.rejected_incident_transitions
        || lhs.matched_generator_sites
            != rhs.matched_generator_sites) {
        return false;
    }
    for (std::size_t index = 0;
         index < lhs.nodes.size();
         ++index) {
        if (lhs.nodes[index].node_id
                != rhs.nodes[index].node_id
            || lhs.nodes[index].provenance_ids
                != rhs.nodes[index].provenance_ids
            || lhs.nodes[index].generator_site_ids
                != rhs.nodes[index].generator_site_ids
            || lhs.nodes[index].parent_site_ids
                != rhs.nodes[index].parent_site_ids) {
            return false;
        }
    }
    for (std::size_t index = 0;
         index < lhs.edges.size();
         ++index) {
        if (lhs.edges[index].edge_id
                != rhs.edges[index].edge_id
            || lhs.edges[index].primitive_kind
                != rhs.edges[index].primitive_kind
            || lhs.edges[index].source_node_id
                != rhs.edges[index].source_node_id
            || lhs.edges[index].target_node_id
                != rhs.edges[index].target_node_id
            || lhs.edges[index].generator_site_ids
                != rhs.edges[index].generator_site_ids
            || lhs.edges[index].parent_site_ids
                != rhs.edges[index].parent_site_ids
            || lhs.edges[index].original_dual_id
                != rhs.edges[index].original_dual_id
            || !endpoints_equal(
                lhs.edges[index].source_endpoint,
                rhs.edges[index].source_endpoint)
            || !endpoints_equal(
                lhs.edges[index].target_endpoint,
                rhs.edges[index].target_endpoint)) {
            return false;
        }
    }
    return true;
}

bool parallel_chart_and_clearance_are_exact()
{
    const MatExactOpenSegmentSource2 lower =
        canonical_open_segment_source(
            "lower-segment",
            0,
            0,
            4,
            0);
    const MatExactOpenSegmentSource2 upper =
        canonical_open_segment_source(
            "upper-segment",
            4,
            3,
            0,
            3);
    const RationalPrimitiveParameterization2 primitive =
        parallel_segment_bisector_parameterization(
            upper,
            lower);
    if (primitive.x_coefficients
            != std::vector<CORE::BigRat>{0, 1}
        || primitive.y_coefficients
            != std::vector<CORE::BigRat>{
                CORE::BigRat(3, 2),
                0,
            }
        || primitive.domain_lower.has_value()
        || primitive.domain_upper.has_value()
        || parallel_segment_tangent_parameter(
               primitive,
               4,
               3)
            != 4) {
        return false;
    }
    const ClearanceRootBoundary2 positive =
        parallel_segment_clearance_boundary(
            primitive,
            lower,
            upper,
            2);
    const ClearanceRootBoundary2 plateau =
        parallel_segment_clearance_boundary(
            primitive,
            lower,
            upper,
            CORE::BigRat(9, 4));
    const ClearanceRootBoundary2 negative =
        parallel_segment_clearance_boundary(
            primitive,
            lower,
            upper,
            CORE::BigRat(5, 2));
    const MatExactOpenSegmentSource2 vertical_left =
        canonical_open_segment_source(
            "vertical-left",
            0,
            3,
            0,
            -2);
    const MatExactOpenSegmentSource2 vertical_right =
        canonical_open_segment_source(
            "vertical-right",
            4,
            -1,
            4,
            5);
    const RationalPrimitiveParameterization2 vertical =
        parallel_segment_bisector_parameterization(
            vertical_right,
            vertical_left);
    const MatExactOpenSegmentSource2 diagonal_lower =
        canonical_open_segment_source(
            "diagonal-lower",
            4,
            4,
            0,
            0);
    const MatExactOpenSegmentSource2 diagonal_upper =
        canonical_open_segment_source(
            "diagonal-upper",
            0,
            2,
            4,
            6);
    const RationalPrimitiveParameterization2 diagonal =
        parallel_segment_bisector_parameterization(
            diagonal_lower,
            diagonal_upper);
    return positive.constant_sign == CGAL::POSITIVE
        && plateau.constant_sign == CGAL::ZERO
        && negative.constant_sign == CGAL::NEGATIVE
        && positive.roots.empty()
        && plateau.roots.empty()
        && negative.roots.empty()
        && vertical.x_coefficients
            == std::vector<CORE::BigRat>{2, 0}
        && vertical.y_coefficients
            == std::vector<CORE::BigRat>{0, 1}
        && diagonal.x_coefficients
            == std::vector<CORE::BigRat>{
                CORE::BigRat(-1, 2),
                1,
            }
        && diagonal.y_coefficients
            == std::vector<CORE::BigRat>{
                CORE::BigRat(1, 2),
                1,
            };
}

bool unsupported_parallel_inputs_fail_loudly()
{
    const MatExactOpenSegmentSource2 horizontal =
        canonical_open_segment_source(
            "horizontal",
            0,
            0,
            4,
            0);
    const MatExactOpenSegmentSource2 vertical =
        canonical_open_segment_source(
            "vertical",
            0,
            0,
            0,
            4);
    const MatExactOpenSegmentSource2 coincident =
        canonical_open_segment_source(
            "coincident",
            -2,
            0,
            8,
            0);
    bool nonparallel_rejected = false;
    bool coincident_rejected = false;
    bool duplicate_identity_rejected = false;
    bool mismatched_clearance_rejected = false;
    bool rescaled_chart_rejected = false;
    bool negative_radius_rejected = false;
    bool empty_feature_rejected = false;
    bool external_limiter_rejected = false;
    bool mismatched_owned_domain_rejected = false;
    bool ownerless_domain_rejected = false;
    bool bounded_outside_feature_rejected = false;
    bool bounded_ownerless_domain_rejected = false;
    try {
        static_cast<void>(
            parallel_segment_bisector_parameterization(
                horizontal,
                vertical));
    } catch (const NonparallelSegmentSupportsError&) {
        nonparallel_rejected = true;
    }
    try {
        static_cast<void>(
            parallel_segment_bisector_parameterization(
                horizontal,
                coincident));
    } catch (const CoincidentSegmentSupportsError&) {
        coincident_rejected = true;
    }
    try {
        static_cast<void>(
            parallel_segment_bisector_parameterization(
                horizontal,
                canonical_open_segment_source(
                    "horizontal",
                    0,
                    3,
                    4,
                    3)));
    } catch (
        const DuplicateOpenSegmentSourceIdentityError&) {
        duplicate_identity_rejected = true;
    }
    const MatExactOpenSegmentSource2 upper =
        canonical_open_segment_source(
            "upper",
            0,
            3,
            4,
            3);
    try {
        static_cast<void>(
            parallel_segment_clearance_boundary(
                {
                    {0, 1},
                    {1, 0},
                    std::nullopt,
                    std::nullopt,
                },
                horizontal,
                upper,
                0));
    } catch (
        const MismatchedParallelSegmentClearanceError&) {
        mismatched_clearance_rejected = true;
    }
    try {
        static_cast<void>(
            parallel_segment_clearance_boundary(
                {
                    {0, 2},
                    {
                        CORE::BigRat(3, 2),
                        0,
                    },
                    std::nullopt,
                    std::nullopt,
                },
                horizontal,
                upper,
                0));
    } catch (
        const MismatchedParallelSegmentClearanceError&) {
        rescaled_chart_rejected = true;
    }
    try {
        const RationalPrimitiveParameterization2 primitive =
            parallel_segment_bisector_parameterization(
                horizontal,
                upper);
        static_cast<void>(
            parallel_segment_clearance_boundary(
                primitive,
                horizontal,
                upper,
                -1));
    } catch (const NegativeClearanceRadiusSquaredError&) {
        negative_radius_rejected = true;
    }
    try {
        static_cast<void>(
            segment_site_disjoint_parallel_segment_graph_spike());
    } catch (
        const EmptyParallelSegmentFeatureDomainError&) {
        empty_feature_rejected = true;
    }
    try {
        static_cast<void>(
            segment_site_external_limited_parallel_segment_graph_spike());
    } catch (
        const UnsupportedSegmentSegmentLimiterError&) {
        external_limiter_rejected = true;
    }
    try {
        RationalPrimitiveParameterization2 primitive =
            parallel_segment_bisector_parameterization(
                horizontal,
                upper);
        primitive.domain_lower = 0;
        primitive.domain_upper = 4;
        ExactAlgebraicKernel1 kernel;
        MatDomainPolygon2 outer;
        outer.push_back({-1, -1});
        outer.push_back({5, -1});
        outer.push_back({5, 4});
        outer.push_back({-1, 4});
        static_cast<void>(
            clip_owned_linear_clearance_components(
                "mismatched-owned-domain",
                primitive,
                {
                    kernel.construct_algebraic_real_1_object()(1),
                    {"wrong-lower"},
                },
                {
                    kernel.construct_algebraic_real_1_object()(4),
                    {"upper"},
                },
                parallel_segment_clearance_boundary(
                    primitive,
                    horizontal,
                    upper,
                    2),
                MatDomainPolygonWithHoles2(outer)));
    } catch (const MismatchedLinearCellDomainError&) {
        mismatched_owned_domain_rejected = true;
    }
    try {
        RationalPrimitiveParameterization2 primitive =
            parallel_segment_bisector_parameterization(
                horizontal,
                upper);
        primitive.domain_lower = 0;
        primitive.domain_upper = 4;
        ExactAlgebraicKernel1 kernel;
        MatDomainPolygon2 outer;
        outer.push_back({-1, -1});
        outer.push_back({5, -1});
        outer.push_back({5, 4});
        outer.push_back({-1, 4});
        static_cast<void>(
            clip_owned_linear_clearance_components(
                "ownerless-domain",
                primitive,
                {
                    kernel.construct_algebraic_real_1_object()(0),
                    {},
                },
                {
                    kernel.construct_algebraic_real_1_object()(4),
                    {"upper"},
                },
                parallel_segment_clearance_boundary(
                    primitive,
                    horizontal,
                    upper,
                    2),
                MatDomainPolygonWithHoles2(outer)));
    } catch (
        const MissingOwnedLinearCellProvenanceError&) {
        ownerless_domain_rejected = true;
    }
    try {
        RationalPrimitiveParameterization2 primitive =
            parallel_segment_bisector_parameterization(
                horizontal,
                upper);
        primitive.domain_lower = 0;
        primitive.domain_upper = 4;
        ExactAlgebraicKernel1 kernel;
        MatDomainPolygon2 outer;
        outer.push_back({-2, -1});
        outer.push_back({5, -1});
        outer.push_back({5, 4});
        outer.push_back({-2, 4});
        static_cast<void>(
            clip_bounded_linear_clearance_components(
                "outside-bounded-domain",
                primitive,
                {
                    kernel.construct_algebraic_real_1_object()(-1),
                    {"external-lower"},
                },
                {
                    kernel.construct_algebraic_real_1_object()(3),
                    {"external-upper"},
                },
                parallel_segment_clearance_boundary(
                    primitive,
                    horizontal,
                    upper,
                    2),
                MatDomainPolygonWithHoles2(outer)));
    } catch (const MismatchedLinearCellDomainError&) {
        bounded_outside_feature_rejected = true;
    }
    try {
        RationalPrimitiveParameterization2 primitive =
            parallel_segment_bisector_parameterization(
                horizontal,
                upper);
        primitive.domain_lower = 0;
        primitive.domain_upper = 4;
        ExactAlgebraicKernel1 kernel;
        MatDomainPolygon2 outer;
        outer.push_back({-1, -1});
        outer.push_back({5, -1});
        outer.push_back({5, 4});
        outer.push_back({-1, 4});
        static_cast<void>(
            clip_bounded_linear_clearance_components(
                "ownerless-bounded-domain",
                primitive,
                {
                    kernel.construct_algebraic_real_1_object()(0),
                    {},
                },
                {
                    kernel.construct_algebraic_real_1_object()(3),
                    {"external-upper"},
                },
                parallel_segment_clearance_boundary(
                    primitive,
                    horizontal,
                    upper,
                    2),
                MatDomainPolygonWithHoles2(outer)));
    } catch (
        const MissingOwnedLinearCellProvenanceError&) {
        bounded_ownerless_domain_rejected = true;
    }
    return nonparallel_rejected
        && coincident_rejected
        && duplicate_identity_rejected
        && mismatched_clearance_rejected
        && rescaled_chart_rejected
        && negative_radius_rejected
        && empty_feature_rejected
        && external_limiter_rejected
        && mismatched_owned_domain_rejected
        && ownerless_domain_rejected
        && bounded_outside_feature_rejected
        && bounded_ownerless_domain_rejected;
}

bool segment_segment_production_graph_gate()
{
    const MatExactGraph2 graph =
        segment_site_segment_segment_graph_spike();
    const MatExactGraph2 reversed =
        segment_site_reversed_segment_segment_graph_spike();
    const MatExactGraph2 positive =
        segment_site_segment_segment_graph_spike(2);
    const MatExactGraph2 negative =
        segment_site_segment_segment_graph_spike(
            CORE::BigRat(5, 2));
    const MatExactGraph2 domain_coincident =
        segment_site_domain_coincident_parallel_segment_graph_spike();
    if (graph.edges.size() != 1
        || graph.nodes.size() != 2
        || graph.rejected_incident_transitions != 0
        || graph.matched_generator_sites != 7
        || !graphs_equal(graph, reversed)
        || !graphs_equal(graph, positive)
        || !negative.edges.empty()
        || !negative.nodes.empty()
        || negative.matched_generator_sites != 7
        || domain_coincident.edges.size() != 1
        || domain_coincident.nodes.size() != 2) {
        return false;
    }
    const MatExactGraphEdge2& edge = graph.edges.front();
    if (edge.primitive_kind != "LINE"
        || edge.original_dual_id.empty()
        || edge.generator_site_ids
            != std::vector<std::string>{
                "lower-segment",
                "upper-segment",
            }
        || edge.parent_site_ids
            != edge.generator_site_ids
        || !edge.source_endpoint.parameter.has_value()
        || !edge.target_endpoint.parameter.has_value()) {
        return false;
    }

    ExactAlgebraicKernel1 kernel;
    const auto zero =
        kernel.construct_algebraic_real_1_object()(0);
    const auto four =
        kernel.construct_algebraic_real_1_object()(4);
    const std::string dual_id =
        stable_dual_identity_v1(
            "segment-segment",
            {
                "lower-segment",
                "upper-segment",
            });
    if (edge.original_dual_id != dual_id) {
        return false;
    }
    std::vector<std::string> expected_source_provenance{
        "upper-left",
        algebraic_root_identity_v1(zero),
        dual_id + "/D-outer/edge-3",
    };
    std::vector<std::string> expected_target_provenance{
        "upper-right",
        algebraic_root_identity_v1(four),
        dual_id + "/D-outer/edge-1",
    };
    union_stable_ids(expected_source_provenance, {});
    union_stable_ids(expected_target_provenance, {});
    const MatExactGraphEdge2& coincident_edge =
        domain_coincident.edges.front();
    return kernel.compare_1_object()(
               *edge.source_endpoint.parameter,
               zero)
            == CGAL::EQUAL
        && kernel.compare_1_object()(
               *edge.target_endpoint.parameter,
               four)
            == CGAL::EQUAL
        && has_provenance(
            edge.source_endpoint,
            "upper-left")
        && has_provenance(
            edge.target_endpoint,
            "upper-right")
        && has_provenance(
            edge.source_endpoint,
            algebraic_root_identity_v1(zero))
        && has_provenance(
            edge.target_endpoint,
            algebraic_root_identity_v1(four))
        && coincident_edge.source_endpoint.provenance_ids
            == expected_source_provenance
        && coincident_edge.target_endpoint.provenance_ids
            == expected_target_provenance;
}

bool point_limited_parallel_segment_segment_production_graph_gate()
{
    const MatExactGraph2 graph =
        segment_site_point_limited_parallel_segment_graph_spike();
    const MatExactGraph2 repeated =
        segment_site_point_limited_parallel_segment_graph_spike();
    const MatExactGraph2 reversed =
        segment_site_reversed_point_limited_parallel_segment_graph_spike();
    const MatExactGraph2 positive =
        segment_site_point_limited_parallel_segment_graph_spike(2);
    const MatExactGraph2 negative =
        segment_site_point_limited_parallel_segment_graph_spike(
            CORE::BigRat(5, 2));
    if (graph.edges.size() != 1
        || graph.nodes.size() != 2
        || graph.rejected_incident_transitions != 0
        || graph.matched_generator_sites != 7
        || !graphs_equal(graph, repeated)
        || !graphs_equal(graph, reversed)
        || !graphs_equal(graph, positive)
        || !negative.edges.empty()
        || !negative.nodes.empty()
        || negative.matched_generator_sites != 7) {
        return false;
    }

    const std::vector<std::string> parent_sites{
        "lower-segment",
        "upper-segment",
    };
    const std::string dual_id =
        stable_dual_identity_v1(
            "segment-segment",
            parent_sites);
    const std::string limiter_root_id =
        algebraic_root_id_v1(
            {-7, 2},
            0);
    const MatExactGraphEdge2& edge =
        graph.edges.front();
    if (edge.primitive_kind != "LINE"
        || edge.original_dual_id != dual_id
        || edge.generator_site_ids != parent_sites
        || edge.parent_site_ids != parent_sites
        || !edge.source_endpoint.parameter.has_value()
        || !edge.target_endpoint.parameter.has_value()
        || !has_provenance(
            edge.source_endpoint,
            "upper-left")
        || !has_provenance(
            edge.target_endpoint,
            "external-limiter")
        || !has_provenance(
            edge.target_endpoint,
            limiter_root_id)
        || !has_provenance(
            edge.target_endpoint,
            "parallel-point-limiter/"
            "equation-factor-multiplicity/1")) {
        return false;
    }
    ExactAlgebraicKernel1 kernel;
    return kernel.compare_1_object()(
               *edge.source_endpoint.parameter,
               kernel.construct_algebraic_real_1_object()(0))
            == CGAL::EQUAL
        && kernel.compare_1_object()(
               *edge.target_endpoint.parameter,
               kernel.construct_algebraic_real_1_object()(
                   CORE::BigRat(7, 2)))
            == CGAL::EQUAL;
}

bool segment_limited_parallel_segment_segment_production_graph_gate()
{
    bool negative_radius_rejected = false;
    try {
        static_cast<void>(
            segment_site_segment_limited_parallel_segment_graph_spike(-1));
    } catch (
        const NegativeClearanceRadiusSquaredError&) {
        negative_radius_rejected = true;
    }
    const MatExactGraph2 graph =
        segment_site_segment_limited_parallel_segment_graph_spike();
    const MatExactGraph2 repeated =
        segment_site_segment_limited_parallel_segment_graph_spike();
    const MatExactGraph2 reversed =
        segment_site_reversed_segment_limited_parallel_segment_graph_spike();
    const MatExactGraph2 positive =
        segment_site_segment_limited_parallel_segment_graph_spike(2);
    const MatExactGraph2 negative =
        segment_site_segment_limited_parallel_segment_graph_spike(
            CORE::BigRat(5, 2));
    if (graph.edges.size() != 1
        || graph.nodes.size() != 2
        || graph.rejected_incident_transitions != 0
        || graph.matched_generator_sites != 9
        || !negative_radius_rejected
        || !graphs_equal(graph, repeated)
        || !graphs_equal(graph, reversed)
        || !graphs_equal(graph, positive)
        || !negative.edges.empty()
        || !negative.nodes.empty()
        || negative.matched_generator_sites != 9) {
        return false;
    }

    const std::vector<std::string> parent_sites{
        "lower-segment",
        "upper-segment",
    };
    const std::string dual_id =
        stable_dual_identity_v1(
            "segment-segment",
            parent_sites);
    const std::string limiter_root_id =
        algebraic_root_id_v1(
            {-7, 2},
            0);
    const MatExactGraphEdge2& edge =
        graph.edges.front();
    if (edge.primitive_kind != "LINE"
        || edge.original_dual_id != dual_id
        || edge.generator_site_ids != parent_sites
        || edge.parent_site_ids != parent_sites
        || !edge.source_endpoint.parameter.has_value()
        || !edge.target_endpoint.parameter.has_value()
        || !has_provenance(
            edge.source_endpoint,
            "upper-left")
        || !has_provenance(
            edge.target_endpoint,
            "external-segment-limiter")
        || !has_provenance(
            edge.target_endpoint,
            limiter_root_id)
        || !has_provenance(
            edge.target_endpoint,
            "parallel-segment-limiter/"
            "equation-factor-multiplicity/1")) {
        return false;
    }
    ExactAlgebraicKernel1 kernel;
    return kernel.compare_1_object()(
               *edge.source_endpoint.parameter,
               kernel.construct_algebraic_real_1_object()(0))
            == CGAL::EQUAL
        && kernel.compare_1_object()(
               *edge.target_endpoint.parameter,
               kernel.construct_algebraic_real_1_object()(
                   CORE::BigRat(7, 2)))
            == CGAL::EQUAL;
}

bool rectangle_central_parallel_production_graph_gate()
{
    bool negative_radius_rejected = false;
    try {
        static_cast<void>(
            segment_site_rectangle_central_parallel_graph_spike(-1));
    } catch (
        const NegativeClearanceRadiusSquaredError&) {
        negative_radius_rejected = true;
    }
    const MatExactGraph2 graph =
        segment_site_rectangle_central_parallel_graph_spike();
    const MatExactGraph2 repeated =
        segment_site_rectangle_central_parallel_graph_spike();
    const MatExactGraph2 reversed =
        segment_site_reversed_rectangle_central_parallel_graph_spike();
    const MatExactGraph2 positive =
        segment_site_rectangle_central_parallel_graph_spike(3);
    const MatExactGraph2 negative =
        segment_site_rectangle_central_parallel_graph_spike(5);
    if (graph.edges.size() != 1
        || graph.nodes.size() != 2
        || graph.rejected_incident_transitions != 0
        || graph.matched_generator_sites != 8
        || !negative_radius_rejected
        || !graphs_equal(graph, repeated)
        || !graphs_equal(graph, reversed)
        || !graphs_equal(graph, positive)
        || !negative.edges.empty()
        || !negative.nodes.empty()
        || negative.matched_generator_sites != 8) {
        return false;
    }

    const std::vector<std::string> parent_sites{
        "bottom-segment",
        "top-segment",
    };
    const MatExactGraphEdge2& edge =
        graph.edges.front();
    if (edge.primitive_kind != "LINE"
        || edge.original_dual_id
            != stable_dual_identity_v1(
                "segment-segment",
                parent_sites)
        || edge.generator_site_ids != parent_sites
        || edge.parent_site_ids != parent_sites
        || !edge.source_endpoint.parameter.has_value()
        || !edge.target_endpoint.parameter.has_value()
        || !has_provenance(
            edge.source_endpoint,
            "left-segment")
        || !has_provenance(
            edge.source_endpoint,
            algebraic_root_id_v1(
                {2, 1},
                0))
        || !has_provenance(
            edge.target_endpoint,
            "right-segment")
        || !has_provenance(
            edge.target_endpoint,
            algebraic_root_id_v1(
                {-2, 1},
                0))
        || !has_provenance(
            edge.source_endpoint,
            "parallel-segment-limiter/"
            "equation-factor-multiplicity/1")
        || !has_provenance(
            edge.target_endpoint,
            "parallel-segment-limiter/"
            "equation-factor-multiplicity/1")) {
        return false;
    }
    ExactAlgebraicKernel1 kernel;
    return kernel.compare_1_object()(
               *edge.source_endpoint.parameter,
               kernel.construct_algebraic_real_1_object()(-2))
            == CGAL::EQUAL
        && kernel.compare_1_object()(
               *edge.target_endpoint.parameter,
               kernel.construct_algebraic_real_1_object()(2))
            == CGAL::EQUAL;
}

bool nonparallel_segment_segment_production_graph_gate()
{
    bool negative_radius_rejected = false;
    try {
        static_cast<void>(
            segment_site_nonparallel_segment_segment_graph_spike(-1));
    } catch (const NegativeClearanceRadiusSquaredError&) {
        negative_radius_rejected = true;
    }
    const MatExactGraph2 graph =
        segment_site_nonparallel_segment_segment_graph_spike();
    const MatExactGraph2 repeated =
        segment_site_nonparallel_segment_segment_graph_spike();
    const MatExactGraph2 reversed =
        segment_site_reversed_nonparallel_segment_segment_graph_spike();
    const MatExactGraph2 radius_one =
        segment_site_nonparallel_segment_segment_graph_spike(1);
    const std::vector<std::string> parent_sites{
        "diagonal-segment",
        "lower-segment",
    };
    const std::vector<std::string> segment_features{
        "diagonal-segment",
        "lower-segment",
    };
    const std::vector<std::string> parabola_features{
        "diagonal-segment",
        "lower-right",
    };
    const std::vector<std::string> transition_features{
        "diagonal-segment",
        "lower-right",
        "lower-segment",
    };
    if (graph.edges.size() != 3
        || graph.nodes.size() != 4
        || graph.rejected_incident_transitions != 1
        || graph.matched_generator_sites != 6
        || !negative_radius_rejected
        || !graphs_equal(graph, repeated)
        || !graphs_equal(graph, reversed)) {
        return false;
    }

    std::map<std::string, std::vector<std::string>>
        primitive_kinds_by_dual;
    std::map<std::string, std::size_t> node_degree;
    std::map<std::string, std::size_t> edge_id_counts;
    std::size_t parabola_count = 0;
    std::size_t line_count = 0;
    for (const MatExactGraphEdge2& edge : graph.edges) {
        if (edge.parent_site_ids != parent_sites
            || edge.original_dual_id.empty()
            || ++edge_id_counts[edge.edge_id] != 1) {
            return false;
        }
        if (edge.primitive_kind == "LINE") {
            if (edge.generator_site_ids
                != segment_features) {
                return false;
            }
            ++line_count;
        } else if (edge.primitive_kind == "PARABOLA") {
            if (edge.generator_site_ids
                != parabola_features) {
                return false;
            }
            ++parabola_count;
        } else {
            return false;
        }
        primitive_kinds_by_dual[edge.original_dual_id]
            .push_back(edge.primitive_kind);
        ++node_degree[edge.source_node_id];
        ++node_degree[edge.target_node_id];
    }
    if (line_count != 2 || parabola_count != 1
        || primitive_kinds_by_dual.size() != 2) {
        return false;
    }

    bool single_line_dual = false;
    bool composite_line_parabola_dual = false;
    for (auto& [dual_id, primitive_kinds] :
         primitive_kinds_by_dual) {
        if (dual_id.empty()) {
            return false;
        }
        std::sort(
            primitive_kinds.begin(),
            primitive_kinds.end());
        if (primitive_kinds
            == std::vector<std::string>{"LINE"}) {
            single_line_dual = true;
        } else if (primitive_kinds
                   == std::vector<std::string>{
                       "LINE",
                       "PARABOLA",
                   }) {
            composite_line_parabola_dual = true;
        } else {
            return false;
        }
    }
    if (!single_line_dual
        || !composite_line_parabola_dual) {
        return false;
    }

    std::size_t internal_nodes = 0;
    std::size_t terminal_nodes = 0;
    for (const MatExactGraphNode2& node : graph.nodes) {
        const auto degree = node_degree.find(node.node_id);
        if (degree == node_degree.end()
            || node.parent_site_ids != parent_sites) {
            return false;
        }
        if (degree->second == 2) {
            if (node.generator_site_ids
                != transition_features) {
                return false;
            }
            ++internal_nodes;
        } else if (degree->second == 1) {
            if (node.generator_site_ids
                != segment_features) {
                return false;
            }
            ++terminal_nodes;
        } else {
            return false;
        }
    }
    if (internal_nodes != 2
        || terminal_nodes != 2
        || radius_one.edges.size() != 3
        || radius_one.nodes.size() != 5
        || radius_one.rejected_incident_transitions != 1
        || radius_one.matched_generator_sites != 6) {
        return false;
    }
    std::map<std::string, std::size_t> clipped_degree;
    for (const MatExactGraphEdge2& edge :
         radius_one.edges) {
        if (edge.parent_site_ids != parent_sites) {
            return false;
        }
        ++clipped_degree[edge.source_node_id];
        ++clipped_degree[edge.target_node_id];
    }
    std::size_t clipped_internal = 0;
    std::size_t clipped_terminal = 0;
    for (const auto& [node_id, degree] :
         clipped_degree) {
        if (node_id.empty()) {
            return false;
        }
        clipped_internal += degree == 2;
        clipped_terminal += degree == 1;
    }
    return clipped_internal == 1
        && clipped_terminal == 4;
}

bool point_limited_parallel_failures_are_named()
{
    using Point = MatTraits::Point_2;
    using Site = MatTraits::Site_2;
    const Point lower_left(-2, 0);
    const Point lower_right(6, 0);
    const Point upper_left(0, 3);
    const Point upper_right(4, 3);
    const Point limiter_point(
        5,
        CORE::BigRat(3, 2));
    const std::vector<GeneratorSite2> generators{
        {"lower-left", Site::construct_site_2(lower_left)},
        {"lower-right", Site::construct_site_2(lower_right)},
        {"upper-left", Site::construct_site_2(upper_left)},
        {"upper-right", Site::construct_site_2(upper_right)},
        {
            "external-limiter",
            Site::construct_site_2(limiter_point),
        },
        {
            "lower-segment",
            Site::construct_site_2(
                lower_left,
                lower_right),
        },
        {
            "upper-segment",
            Site::construct_site_2(
                upper_left,
                upper_right),
        },
    };
    SegmentSiteDelaunay2 delaunay;
    delaunay.insert(lower_left, lower_right);
    delaunay.insert(upper_left, upper_right);
    delaunay.insert(limiter_point);
    require_generator_site_bijection(
        delaunay,
        generators);
    SegmentSiteVoronoi2 voronoi(delaunay);
    const std::vector<std::string> generator_ids{
        "lower-segment",
        "upper-segment",
    };
    std::vector<SegmentSiteVoronoi2::Halfedge_handle>
        matching;
    for (auto halfedge = voronoi.halfedges_begin();
         halfedge != voronoi.halfedges_end();
         ++halfedge) {
        const MatTraits::Site_2 up =
            halfedge->up()->site();
        const MatTraits::Site_2 down =
            halfedge->down()->site();
        if (up.is_segment()
            && down.is_segment()
            && ordered_generator_site_ids(
                   stable_generator_site_id(
                       up,
                       generators),
                   stable_generator_site_id(
                       down,
                       generators))
                == generator_ids
            && stable_generator_site_id(
                   up,
                   generators)
                == generator_ids.front()) {
            matching.push_back(halfedge);
        }
    }
    if (matching.size() != 1) {
        return false;
    }

    const MatExactOpenSegmentSource2 lower_segment =
        canonical_open_segment_source(
            "lower-segment",
            -2,
            0,
            6,
            0);
    const MatExactOpenSegmentSource2 upper_segment =
        canonical_open_segment_source(
            "upper-segment",
            0,
            3,
            4,
            3);
    RationalPrimitiveParameterization2 primitive =
        parallel_segment_bisector_parameterization(
            lower_segment,
            upper_segment);
    primitive.domain_lower = 0;
    primitive.domain_upper = 4;
    const RationalDomainRoot2 lower{
        0,
        {"upper-left"},
    };
    const RationalDomainRoot2 upper{
        4,
        {"upper-right"},
    };
    const MatExactPointSiteSource2 limiter{
        "external-limiter",
        {5, 0},
        {
            CORE::BigRat(3, 2),
            0,
        },
        1,
    };

    bool missing_rejected = false;
    bool duplicate_rejected = false;
    bool nonrational_rejected = false;
    bool mismatch_rejected = false;
    bool chart_rejected = false;
    try {
        static_cast<void>(
            bind_parallel_segment_segment_cell_endpoints(
                primitive,
                lower,
                upper,
                lower_segment,
                upper_segment,
                {},
                generator_ids,
                generators,
                voronoi,
                matching.front()));
    } catch (
        const UnsupportedSegmentSegmentLimiterError&) {
        missing_rejected = true;
    }
    try {
        static_cast<void>(
            bind_parallel_segment_segment_cell_endpoints(
                primitive,
                lower,
                upper,
                lower_segment,
                upper_segment,
                {limiter, limiter},
                generator_ids,
                generators,
                voronoi,
                matching.front()));
    } catch (
        const AmbiguousParallelSegmentPointLimiterError&) {
        duplicate_rejected = true;
    }
    try {
        static_cast<void>(
            bind_parallel_segment_segment_cell_endpoints(
                primitive,
                lower,
                upper,
                lower_segment,
                upper_segment,
                {
                    {
                        "external-limiter",
                        {5, 1},
                        {
                            CORE::BigRat(3, 2),
                            0,
                        },
                        2,
                    },
                },
                generator_ids,
                generators,
                voronoi,
                matching.front()));
    } catch (
        const NonRationalParallelSegmentPointLimiterError&) {
        nonrational_rejected = true;
    }
    try {
        static_cast<void>(
            bind_parallel_segment_segment_cell_endpoints(
                primitive,
                lower,
                upper,
                lower_segment,
                upper_segment,
                {
                    {
                        "external-limiter",
                        {4, 0},
                        {
                            CORE::BigRat(3, 2),
                            0,
                        },
                        1,
                    },
                },
                generator_ids,
                generators,
                voronoi,
                matching.front()));
    } catch (
        const MismatchedLiveSegmentSegmentBridgeError&) {
        mismatch_rejected = true;
    }
    try {
        RationalPrimitiveParameterization2 mutated =
            primitive;
        mutated.x_coefficients[1] *= 2;
        static_cast<void>(
            bind_parallel_segment_segment_cell_endpoints(
                mutated,
                lower,
                upper,
                lower_segment,
                upper_segment,
                {limiter},
                generator_ids,
                generators,
                voronoi,
                matching.front()));
    } catch (
        const MismatchedLiveSegmentSegmentBridgeError&) {
        chart_rejected = true;
    }
    return missing_rejected
        && duplicate_rejected
        && nonrational_rejected
        && mismatch_rejected
        && chart_rejected;
}

bool segment_limited_parallel_failures_are_named()
{
    using Point = MatTraits::Point_2;
    using Site = MatTraits::Site_2;
    const Point lower_left(-2, 0);
    const Point lower_right(6, 0);
    const Point upper_left(0, 3);
    const Point upper_right(4, 3);
    const Point limiter_lower(5, 1);
    const Point limiter_upper(5, 2);
    const std::vector<GeneratorSite2> generators{
        {"lower-left", Site::construct_site_2(lower_left)},
        {"lower-right", Site::construct_site_2(lower_right)},
        {"upper-left", Site::construct_site_2(upper_left)},
        {"upper-right", Site::construct_site_2(upper_right)},
        {"external-limiter-lower", Site::construct_site_2(limiter_lower)},
        {"external-limiter-upper", Site::construct_site_2(limiter_upper)},
        {
            "lower-segment",
            Site::construct_site_2(lower_left, lower_right),
        },
        {
            "upper-segment",
            Site::construct_site_2(upper_left, upper_right),
        },
        {
            "external-segment-limiter",
            Site::construct_site_2(limiter_lower, limiter_upper),
        },
    };
    SegmentSiteDelaunay2 delaunay;
    delaunay.insert(lower_left, lower_right);
    delaunay.insert(upper_left, upper_right);
    delaunay.insert(limiter_lower, limiter_upper);
    require_generator_site_bijection(
        delaunay,
        generators);
    SegmentSiteVoronoi2 voronoi(delaunay);
    const std::vector<std::string> generator_ids{
        "lower-segment",
        "upper-segment",
    };
    std::vector<SegmentSiteVoronoi2::Halfedge_handle>
        matching;
    for (auto halfedge = voronoi.halfedges_begin();
         halfedge != voronoi.halfedges_end();
         ++halfedge) {
        const MatTraits::Site_2 up =
            halfedge->up()->site();
        const MatTraits::Site_2 down =
            halfedge->down()->site();
        if (up.is_segment()
            && down.is_segment()
            && ordered_generator_site_ids(
                   stable_generator_site_id(
                       up,
                       generators),
                   stable_generator_site_id(
                       down,
                       generators))
                == generator_ids
            && stable_generator_site_id(
                   up,
                   generators)
                == generator_ids.front()) {
            matching.push_back(halfedge);
        }
    }
    if (matching.size() != 1) {
        return false;
    }

    const MatExactOpenSegmentSource2 lower_segment =
        canonical_open_segment_source(
            "lower-segment",
            -2,
            0,
            6,
            0);
    const MatExactOpenSegmentSource2 upper_segment =
        canonical_open_segment_source(
            "upper-segment",
            0,
            3,
            4,
            3);
    const MatExactOpenSegmentSource2 limiter =
        canonical_open_segment_source(
            "external-segment-limiter",
            5,
            1,
            5,
            2);
    RationalPrimitiveParameterization2 primitive =
        parallel_segment_bisector_parameterization(
            lower_segment,
            upper_segment);
    primitive.domain_lower = 0;
    primitive.domain_upper = 4;
    const RationalDomainRoot2 lower{
        0,
        {"upper-left"},
    };
    const RationalDomainRoot2 upper{
        4,
        {"upper-right"},
    };

    bool missing_rejected = false;
    bool duplicate_rejected = false;
    bool mismatch_rejected = false;
    try {
        static_cast<void>(
            bind_parallel_segment_segment_cell_endpoints(
                primitive,
                lower,
                upper,
                lower_segment,
                upper_segment,
                {},
                {},
                generator_ids,
                generators,
                voronoi,
                matching.front()));
    } catch (
        const UnsupportedSegmentSegmentLimiterError&) {
        missing_rejected = true;
    }
    try {
        static_cast<void>(
            bind_parallel_segment_segment_cell_endpoints(
                primitive,
                lower,
                upper,
                lower_segment,
                upper_segment,
                {},
                {limiter, limiter},
                generator_ids,
                generators,
                voronoi,
                matching.front()));
    } catch (
        const AmbiguousParallelSegmentOpenLimiterError&) {
        duplicate_rejected = true;
    }
    try {
        static_cast<void>(
            bind_parallel_segment_segment_cell_endpoints(
                primitive,
                lower,
                upper,
                lower_segment,
                upper_segment,
                {},
                {
                    canonical_open_segment_source(
                        "external-segment-limiter",
                        4,
                        1,
                        4,
                        2),
                },
                generator_ids,
                generators,
                voronoi,
                matching.front()));
    } catch (
        const MismatchedLiveSegmentSegmentBridgeError&) {
        mismatch_rejected = true;
    }
    return missing_rejected
        && duplicate_rejected
        && mismatch_rejected;
}

bool rectangle_adaptor_characterization_gate()
{
    using Point = MatTraits::Point_2;
    using Site = MatTraits::Site_2;
    const Point lower_left(-4, -2);
    const Point lower_right(4, -2);
    const Point upper_right(4, 2);
    const Point upper_left(-4, 2);
    const std::vector<GeneratorSite2> generators{
        {"lower-left", Site::construct_site_2(lower_left)},
        {"lower-right", Site::construct_site_2(lower_right)},
        {"upper-right", Site::construct_site_2(upper_right)},
        {"upper-left", Site::construct_site_2(upper_left)},
        {
            "bottom-segment",
            Site::construct_site_2(
                lower_left,
                lower_right),
        },
        {
            "right-segment",
            Site::construct_site_2(
                lower_right,
                upper_right),
        },
        {
            "top-segment",
            Site::construct_site_2(
                upper_right,
                upper_left),
        },
        {
            "left-segment",
            Site::construct_site_2(
                upper_left,
                lower_left),
        },
    };
    SegmentSiteDelaunay2 delaunay;
    delaunay.insert(lower_left, lower_right);
    delaunay.insert(lower_right, upper_right);
    delaunay.insert(upper_right, upper_left);
    delaunay.insert(upper_left, lower_left);
    require_generator_site_bijection(
        delaunay,
        generators);
    SegmentSiteVoronoi2 voronoi(delaunay);
    const std::set<std::vector<std::string>>
        expected_incident_point_segment_pairs{
            {"bottom-segment", "lower-left"},
            {"bottom-segment", "lower-right"},
            {"left-segment", "lower-left"},
            {"left-segment", "upper-left"},
            {"lower-right", "right-segment"},
            {"right-segment", "upper-right"},
            {"top-segment", "upper-left"},
            {"top-segment", "upper-right"},
        };
    const std::map<
        std::vector<std::string>,
        std::vector<std::string>>
        expected_segment_segment_owners{
            {
                {"bottom-segment", "left-segment"},
                {"top-segment", "lower-left"},
            },
            {
                {"bottom-segment", "right-segment"},
                {"lower-right", "top-segment"},
            },
            {
                {"bottom-segment", "top-segment"},
                {"right-segment", "left-segment"},
            },
            {
                {"left-segment", "top-segment"},
                {"bottom-segment", "upper-left"},
            },
            {
                {"right-segment", "top-segment"},
                {"upper-right", "bottom-segment"},
            },
        };
    std::set<std::vector<std::string>>
        incident_point_segment_pairs;
    std::map<
        std::vector<std::string>,
        std::vector<std::string>>
        segment_segment_owners;
    for (auto halfedge = voronoi.halfedges_begin();
         halfedge != voronoi.halfedges_end();
         ++halfedge) {
        const MatTraits::Site_2 up =
            halfedge->up()->site();
        const MatTraits::Site_2 down =
            halfedge->down()->site();
        const std::string up_id =
            stable_generator_site_id(
                up,
                generators);
        const std::string down_id =
            stable_generator_site_id(
                down,
                generators);
        const std::vector<std::string> pair =
            ordered_generator_site_ids(
                up_id,
                down_id);
        if (up_id != pair.front()) {
            continue;
        }
        const CGAL::Object primal =
            voronoi.dual().primal(
                halfedge->dual());
        MatTraits::Ray_2 ray;
        MatTraits::Segment_2 segment;
        if (up.is_segment() && down.is_segment()) {
            if (!halfedge->has_source()
                || !halfedge->has_target()
                || !CGAL::assign(segment, primal)) {
                return false;
            }
            segment_segment_owners.emplace(
                pair,
                std::vector<std::string>{
                    stable_generator_site_id(
                        halfedge->left()->site(),
                        generators),
                    stable_generator_site_id(
                        halfedge->right()->site(),
                        generators),
                });
            continue;
        }
        if (up.is_point() == down.is_point()
            || halfedge->has_source()
                == halfedge->has_target()
            || !CGAL::assign(ray, primal)) {
            return false;
        }
        incident_point_segment_pairs.insert(pair);
    }
    return incident_point_segment_pairs
            == expected_incident_point_segment_pairs
        && segment_segment_owners
            == expected_segment_segment_owners;
}

} // namespace

bool segment_segment_producer_gate()
{
    using Point = MatTraits::Point_2;
    using Site = MatTraits::Site_2;
    const Point lower_left(0, 0);
    const Point lower_right(4, 0);
    const Point upper_left(0, 3);
    const Point upper_right(4, 3);
    const Point limiter(6, 1);
    const std::vector<GeneratorSite2> generators{
        {"lower-left", Site::construct_site_2(lower_left)},
        {"lower-right", Site::construct_site_2(lower_right)},
        {"upper-left", Site::construct_site_2(upper_left)},
        {"upper-right", Site::construct_site_2(upper_right)},
        {"limiter", Site::construct_site_2(limiter)},
        {
            "lower-segment",
            Site::construct_site_2(lower_left, lower_right),
        },
        {
            "upper-segment",
            Site::construct_site_2(upper_left, upper_right),
        },
    };

    SegmentSiteDelaunay2 delaunay;
    delaunay.insert(lower_left, lower_right);
    delaunay.insert(upper_left, upper_right);
    delaunay.insert(limiter);
    require_generator_site_bijection(delaunay, generators);
    SegmentSiteVoronoi2 voronoi(delaunay);

    const std::vector<std::string> expected_generators{
        "lower-segment",
        "upper-segment",
    };
    std::size_t matched = 0;
    for (auto halfedge = voronoi.halfedges_begin();
         halfedge != voronoi.halfedges_end();
         ++halfedge) {
        const MatTraits::Site_2 up = halfedge->up()->site();
        const MatTraits::Site_2 down = halfedge->down()->site();
        if (!up.is_segment() || !down.is_segment()) {
            continue;
        }
        const std::string up_id =
            stable_generator_site_id(up, generators);
        const std::string down_id =
            stable_generator_site_id(down, generators);
        if (ordered_generator_site_ids(up_id, down_id)
                != expected_generators
            || up_id != expected_generators.front()) {
            continue;
        }
        if (!halfedge->has_source()
            || !halfedge->has_target()) {
            return false;
        }

        MatTraits::Segment_2 primal;
        if (!CGAL::assign(
                primal,
                voronoi.dual().primal(
                    halfedge->dual()))) {
            return false;
        }
        const MatTraits::Segment_2 adaptor_segment(
            halfedge->source()->point(),
            halfedge->target()->point());
        if (!same_unoriented_segment(
                primal,
                adaptor_segment)) {
            return false;
        }

        const std::string left_id =
            stable_generator_site_id(
                halfedge->left()->site(),
                generators);
        const std::string right_id =
            stable_generator_site_id(
                halfedge->right()->site(),
                generators);
        if (left_id != "upper-right"
            || right_id != "upper-left"
            || halfedge->source()->point()
                != MatTraits::Point_2(
                    4,
                    CORE::BigRat(3, 2))
            || halfedge->target()->point()
                != MatTraits::Point_2(
                    0,
                    CORE::BigRat(3, 2))) {
            return false;
        }

        RationalPrimitiveParameterization2 primitive =
            parallel_segment_bisector_parameterization(
                canonical_open_segment_source(
                    "lower-segment",
                    0,
                    0,
                    4,
                    0),
                canonical_open_segment_source(
                    "upper-segment",
                    0,
                    3,
                    4,
                    3));
        primitive.domain_lower = 0;
        primitive.domain_upper = 4;
        bool owner_mutation_rejected = false;
        try {
            static_cast<void>(
                bind_parallel_segment_segment_cell_endpoints(
                    primitive,
                    {
                        0,
                        {
                            "lower-left",
                            "upper-left",
                        },
                    },
                    {
                        4,
                        {
                            "lower-right",
                        },
                    },
                    expected_generators,
                    generators,
                    voronoi,
                    halfedge));
        } catch (
            const UnsupportedSegmentSegmentLimiterError&) {
            owner_mutation_rejected = true;
        }
        if (!owner_mutation_rejected) {
            return false;
        }
        ++matched;
    }
    if (matched != 1) {
        return false;
    }
    if (!parallel_chart_and_clearance_are_exact()) {
        return false;
    }
    if (!unsupported_parallel_inputs_fail_loudly()) {
        return false;
    }
    return segment_segment_production_graph_gate()
        && rectangle_adaptor_characterization_gate()
        && point_limited_parallel_segment_segment_production_graph_gate()
        && segment_limited_parallel_segment_segment_production_graph_gate()
        && rectangle_central_parallel_production_graph_gate()
        && point_limited_parallel_failures_are_named()
        && segment_limited_parallel_failures_are_named()
        && nonparallel_segment_segment_production_graph_gate()
        && nonparallel_segment_segment_producer_contract()
        && nonparallel_segment_charts_are_exact()
        && nonparallel_feature_parameters_are_exact()
        && nonparallel_feature_domains_are_exact()
        && quadratic_feature_parameters_order_and_embed()
        && nonparallel_clearance_is_rational_quadratic()
        && nonparallel_feature_clearance_components_are_exact()
        && nonparallel_domain_clipping_is_exact()
        && collinear_domain_intersections_are_exact()
        && unsupported_nonparallel_charts_fail_loudly();
}
