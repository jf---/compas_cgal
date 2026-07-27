#include "segment_site_mat.h"

#include <algorithm>
#include <cstddef>
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
    bool parallel_rejected = false;
    bool degenerate_primitive_rejected = false;
    bool unbound_branch_rejected = false;
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
    return parallel_rejected
        && degenerate_primitive_rejected
        && unbound_branch_rejected;
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
                != rhs.nodes[index].generator_site_ids) {
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
    return nonparallel_rejected
        && coincident_rejected
        && duplicate_identity_rejected
        && mismatched_clearance_rejected
        && rescaled_chart_rejected
        && negative_radius_rejected
        && empty_feature_rejected
        && external_limiter_rejected
        && mismatched_owned_domain_rejected
        && ownerless_domain_rejected;
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
        || edge.generator_site_ids
            != std::vector<std::string>{
                "lower-segment",
                "upper-segment",
            }
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
        && nonparallel_segment_segment_producer_contract()
        && nonparallel_segment_charts_are_exact()
        && unsupported_nonparallel_charts_fail_loudly();
}
