#include "segment_site_clipping.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

#include <CGAL/Polynomial_traits_d.h>

namespace {

using AlgebraicParameter =
    ExactAlgebraicKernel1::Algebraic_real_1;

CGAL::Sign clearance_sign_at(
    const std::vector<ExactAlgebraicInteger1>& coefficients,
    const AlgebraicParameter& parameter,
    const ExactAlgebraicKernel1& kernel)
{
    using Polynomial = ExactAlgebraicKernel1::Polynomial_1;
    const Polynomial polynomial =
        typename CGAL::Polynomial_traits_d<
            Polynomial>::Construct_polynomial()(
            coefficients.begin(),
            coefficients.end());
    return kernel.sign_at_1_object()(polynomial, parameter);
}

CGAL::Sign clearance_sign_on_open_cell(
    const std::vector<ExactAlgebraicInteger1>& coefficients,
    const MatParameterEndpoint2& lower,
    const MatParameterEndpoint2& upper,
    const ExactAlgebraicKernel1& kernel)
{
    CORE::BigRat witness;
    if (!lower.parameter.has_value()
        && !upper.parameter.has_value()) {
        witness = 0;
    } else if (!lower.parameter.has_value()) {
        witness = upper.parameter->low() - 1;
    } else if (!upper.parameter.has_value()) {
        witness = lower.parameter->high() + 1;
    } else {
        witness = kernel.bound_between_1_object()(
            *lower.parameter,
            *upper.parameter);
    }

    CORE::BigRat value = 0;
    for (auto coefficient = coefficients.rbegin();
         coefficient != coefficients.rend();
         ++coefficient) {
        value *= witness;
        value += *coefficient;
    }
    return CGAL::sign(value);
}

CORE::BigRat open_cell_witness(
    const MatParameterEndpoint2& lower,
    const MatParameterEndpoint2& upper,
    const ExactAlgebraicKernel1& kernel)
{
    if (!lower.parameter.has_value()
        && !upper.parameter.has_value()) {
        return 0;
    }
    if (!lower.parameter.has_value()) {
        return upper.parameter->low() - 1;
    }
    if (!upper.parameter.has_value()) {
        return lower.parameter->high() + 1;
    }
    return kernel.bound_between_1_object()(
        *lower.parameter,
        *upper.parameter);
}

void append_polygon_intersections(
    const MatDomainPolygon2& polygon,
    const std::string& ring_id,
    const RationalPrimitiveParameterization2& primitive,
    std::vector<RationalDomainRoot2>& roots)
{
    const CORE::BigRat x0 = primitive.x_coefficients.front();
    const CORE::BigRat y0 = primitive.y_coefficients.front();
    const CORE::BigRat vx =
        primitive.x_coefficients.size() == 1
        ? CORE::BigRat(0)
        : primitive.x_coefficients[1];
    const CORE::BigRat vy =
        primitive.y_coefficients.size() == 1
        ? CORE::BigRat(0)
        : primitive.y_coefficients[1];

    std::size_t edge_index = 0;
    for (auto edge = polygon.edges_begin();
         edge != polygon.edges_end();
         ++edge, ++edge_index) {
        const CORE::BigRat ax = edge->source().x();
        const CORE::BigRat ay = edge->source().y();
        const CORE::BigRat ex =
            edge->target().x() - ax;
        const CORE::BigRat ey =
            edge->target().y() - ay;
        const CORE::BigRat wx = ax - x0;
        const CORE::BigRat wy = ay - y0;
        const CORE::BigRat denominator =
            vx * ey - vy * ex;
        const CORE::BigRat collinearity =
            wx * vy - wy * vx;
        if (denominator == 0) {
            if (collinearity == 0) {
                throw OverlappingDomainBoundaryError(
                    "linear primitive overlaps domain boundary");
            }
            continue;
        }

        const CORE::BigRat parameter =
            (wx * ey - wy * ex) / denominator;
        const CORE::BigRat edge_parameter =
            collinearity / denominator;
        if (edge_parameter < 0 || edge_parameter > 1) {
            continue;
        }
        if (primitive.domain_lower.has_value()
            && parameter < *primitive.domain_lower) {
            continue;
        }
        if (primitive.domain_upper.has_value()
            && parameter > *primitive.domain_upper) {
            continue;
        }
        roots.push_back(
            {
                parameter,
                {
                    ring_id + "/edge-"
                    + std::to_string(edge_index),
                },
            });
    }
}

std::vector<RationalDomainRoot2> linear_domain_roots(
    const std::string& dual_id,
    const RationalPrimitiveParameterization2& primitive,
    const MatDomainPolygonWithHoles2& domain)
{
    if (primitive.x_coefficients.size() > 2
        || primitive.y_coefficients.size() > 2) {
        throw InvalidRationalPrimitiveError(
            "D clipping requires a linear primitive");
    }
    const CORE::BigRat vx =
        primitive.x_coefficients.size() == 1
        ? CORE::BigRat(0)
        : primitive.x_coefficients[1];
    const CORE::BigRat vy =
        primitive.y_coefficients.size() == 1
        ? CORE::BigRat(0)
        : primitive.y_coefficients[1];
    if (vx == 0 && vy == 0) {
        throw InvalidRationalPrimitiveError(
            "linear primitive direction is zero");
    }
    if (domain.is_unbounded()) {
        throw InvalidRationalPrimitiveError(
            "D clipping requires a bounded outer polygon");
    }

    std::vector<RationalDomainRoot2> roots;
    append_polygon_intersections(
        domain.outer_boundary(),
        dual_id + "/D-outer",
        primitive,
        roots);
    std::size_t hole_index = 0;
    for (auto hole = domain.holes_begin();
         hole != domain.holes_end();
         ++hole, ++hole_index) {
        append_polygon_intersections(
            *hole,
            dual_id + "/D-hole-"
                + std::to_string(hole_index),
            primitive,
            roots);
    }
    std::sort(
        roots.begin(),
        roots.end(),
        [](const RationalDomainRoot2& lhs,
           const RationalDomainRoot2& rhs) {
            return lhs.parameter < rhs.parameter;
        });
    std::vector<RationalDomainRoot2> unique;
    for (RationalDomainRoot2& root : roots) {
        if (!unique.empty()
            && unique.back().parameter == root.parameter) {
            unique.back().provenance_ids.insert(
                unique.back().provenance_ids.end(),
                root.provenance_ids.begin(),
                root.provenance_ids.end());
        } else {
            unique.push_back(std::move(root));
        }
    }
    return unique;
}

std::vector<AlgebraicDomainRoot2>
algebraic_linear_domain_roots(
    const std::string& original_dual_id,
    const RationalPrimitiveParameterization2& primitive,
    const MatDomainPolygonWithHoles2& domain)
{
    ExactAlgebraicKernel1 kernel;
    std::vector<AlgebraicDomainRoot2> roots;
    for (const RationalDomainRoot2& root :
         linear_domain_roots(
             original_dual_id,
             primitive,
             domain)) {
        const auto parameter =
            kernel.construct_algebraic_real_1_object()(
                root.parameter);
        std::vector<std::string> provenance =
            root.provenance_ids;
        provenance.push_back(
            algebraic_root_identity_v1(parameter));
        roots.push_back(
            {
                parameter,
                std::move(provenance),
            });
    }
    return roots;
}

void append_source_parabola_intersections(
    const MatDomainPolygon2& polygon,
    const std::string& ring_id,
    const SourceParabolaParameterization2& primitive,
    const MatParameterEndpoint2& domain_lower,
    const MatParameterEndpoint2& domain_upper,
    std::vector<AlgebraicDomainRoot2>& roots)
{
    ExactAlgebraicKernel1 kernel1;
    ExactAlgebraicKernel2 kernel2;
    const auto compare = kernel1.compare_1_object();
    std::size_t edge_index = 0;
    for (auto edge = polygon.edges_begin();
         edge != polygon.edges_end();
         ++edge, ++edge_index) {
        const CORE::BigRat ax = edge->source().x();
        const CORE::BigRat ay = edge->source().y();
        const CORE::BigRat ex = edge->target().x() - ax;
        const CORE::BigRat ey = edge->target().y() - ay;
        const RadicalEquation2 x_equation =
            radical_equation(
                primitive.x_rational,
                primitive.x_radical,
                ax,
                ex);
        const RadicalEquation2 y_equation =
            radical_equation(
                primitive.y_rational,
                primitive.y_radical,
                ay,
                ey);
        const auto x_norm = radical_norm(
            x_equation,
            primitive.radicand);
        const auto y_norm = radical_norm(
            y_equation,
            primitive.radicand);
        if (CGAL::is_zero(x_norm)
            || CGAL::is_zero(y_norm)
            || !kernel2.is_coprime_2_object()(
                x_norm,
                y_norm)) {
            throw OverlappingDomainBoundaryError(
                "source parabola boundary intersection is not isolated");
        }

        std::vector<
            std::pair<
                ExactAlgebraicKernel2::Algebraic_real_2,
                ExactAlgebraicKernel2::Multiplicity_type>>
            solutions;
        kernel2.solve_2_object()(
            x_norm,
            y_norm,
            std::back_inserter(solutions));
        for (const auto& [solution, multiplicity] : solutions) {
            static_cast<void>(multiplicity);
            const auto parameter =
                kernel2.compute_x_2_object()(solution);
            const auto edge_parameter =
                kernel2.compute_y_2_object()(solution);
            if (!radical_equation_holds(
                    x_equation,
                    solution)
                || !radical_equation_holds(
                    y_equation,
                    solution)
                || compare(edge_parameter, CORE::BigRat(0))
                    == CGAL::SMALLER
                || compare(edge_parameter, CORE::BigRat(1))
                    == CGAL::LARGER
                || (domain_lower.parameter.has_value()
                    && compare(
                           parameter,
                           *domain_lower.parameter)
                        == CGAL::SMALLER)
                || (domain_upper.parameter.has_value()
                    && compare(
                           parameter,
                           *domain_upper.parameter)
                        == CGAL::LARGER)) {
                continue;
            }
            roots.push_back(
                {
                    parameter,
                    {
                        ring_id + "/edge-"
                        + std::to_string(edge_index),
                        algebraic_root_identity_v1(parameter),
                    },
                });
        }
    }
}

std::vector<AlgebraicDomainRoot2>
source_parabola_domain_roots(
    const std::string& dual_id,
    const SourceParabolaParameterization2& primitive,
    const MatParameterEndpoint2& domain_lower,
    const MatParameterEndpoint2& domain_upper,
    const MatDomainPolygonWithHoles2& domain)
{
    std::vector<AlgebraicDomainRoot2> roots;
    append_source_parabola_intersections(
        domain.outer_boundary(),
        dual_id + "/D-outer",
        primitive,
        domain_lower,
        domain_upper,
        roots);
    std::size_t hole_index = 0;
    for (auto hole = domain.holes_begin();
         hole != domain.holes_end();
         ++hole, ++hole_index) {
        append_source_parabola_intersections(
            *hole,
            dual_id + "/D-hole-"
                + std::to_string(hole_index),
            primitive,
            domain_lower,
            domain_upper,
            roots);
    }
    ExactAlgebraicKernel1 kernel;
    const auto compare = kernel.compare_1_object();
    std::sort(
        roots.begin(),
        roots.end(),
        [&compare](const AlgebraicDomainRoot2& lhs,
                   const AlgebraicDomainRoot2& rhs) {
            return compare(lhs.parameter, rhs.parameter)
                == CGAL::SMALLER;
        });
    std::vector<AlgebraicDomainRoot2> unique;
    for (AlgebraicDomainRoot2& root : roots) {
        if (!unique.empty()
            && compare(
                   unique.back().parameter,
                   root.parameter)
                == CGAL::EQUAL) {
            unique.back().provenance_ids.insert(
                unique.back().provenance_ids.end(),
                root.provenance_ids.begin(),
                root.provenance_ids.end());
        } else {
            unique.push_back(std::move(root));
        }
    }
    return unique;
}

std::vector<ExactAlgebraicInteger1>
canonical_polynomial_coefficients(
    const ExactAlgebraicKernel1::Polynomial_1& polynomial)
{
    const int degree = CGAL::degree(polynomial);
    std::vector<ExactAlgebraicInteger1> coefficients;
    coefficients.reserve(
        static_cast<std::size_t>(degree + 1));
    for (int index = 0; index <= degree; ++index) {
        coefficients.push_back(polynomial[index]);
    }
    return coefficients;
}

ClearanceRootBoundary2 clearance_boundary_from_polynomial(
    RationalPolynomial clearance,
    const RationalPrimitiveParameterization2* rational_domain)
{
    trim(clearance);
    if (clearance.size() == 1) {
        return {
            CGAL::sign(clearance.front()),
            {},
            {},
        };
    }

    const std::vector<ExactAlgebraicInteger1> coefficients =
        primitive_integer_coefficients(clearance);
    using Polynomial = ExactAlgebraicKernel1::Polynomial_1;
    const Polynomial polynomial =
        typename CGAL::Polynomial_traits_d<
            Polynomial>::Construct_polynomial()(
            coefficients.begin(),
            coefficients.end());
    ExactAlgebraicKernel1 kernel;
    const Polynomial square_free =
        typename CGAL::Polynomial_traits_d<
            Polynomial>::Canonicalize()(
            kernel.make_square_free_1_object()(
                polynomial));
    const std::vector<ExactAlgebraicInteger1>
        square_free_coefficients =
            canonical_polynomial_coefficients(
                square_free);
    std::vector<ExactAlgebraicKernel1::Algebraic_real_1>
        isolated_roots;
    kernel.solve_1_object()(
        square_free,
        true,
        std::back_inserter(isolated_roots));

    std::vector<ClearanceRootEvent2> roots;
    for (std::size_t ordinal = 0;
         ordinal < isolated_roots.size();
         ++ordinal) {
        if (rational_domain != nullptr
            && !root_is_in_domain(
                isolated_roots[ordinal],
                *rational_domain,
                kernel)) {
            continue;
        }
        roots.push_back(
            {
                isolated_roots[ordinal],
                algebraic_root_id_v1(
                    square_free_coefficients,
                    ordinal),
            });
    }
    return {
        std::nullopt,
        coefficients,
        roots,
    };
}

} // namespace

ClearanceRootBoundary2 point_clearance_boundary(
    const RationalPrimitiveParameterization2& primitive,
    const CORE::BigRat& site_x,
    const CORE::BigRat& site_y,
    const CORE::BigRat& radius_squared)
{
    if (primitive.x_coefficients.empty()
        || primitive.y_coefficients.empty()) {
        throw InvalidRationalPrimitiveError(
            "primitive coordinate polynomial is empty");
    }
    if (primitive.domain_lower.has_value()
        && primitive.domain_upper.has_value()
        && *primitive.domain_lower > *primitive.domain_upper) {
        throw InvalidRationalPrimitiveError(
            "primitive parameter domain is reversed");
    }

    RationalPolynomial dx = primitive.x_coefficients;
    RationalPolynomial dy = primitive.y_coefficients;
    dx.front() -= site_x;
    dy.front() -= site_y;
    trim(dx);
    trim(dy);

    RationalPolynomial clearance = multiply(dx, dx);
    add_in_place(clearance, multiply(dy, dy));
    clearance.front() -= radius_squared;
    return clearance_boundary_from_polynomial(
        std::move(clearance),
        &primitive);
}

ClearanceRootBoundary2
source_parabola_clearance_boundary(
    const MatExactPointSiteSource2& point_site,
    const MatExactOpenSegmentSource2& segment_site,
    const CORE::BigRat& radius_squared)
{
    if (radius_squared < 0) {
        throw NegativeClearanceRadiusSquaredError(
            "source parabola squared clearance radius is negative");
    }
    const SourceParabolaParameterization2 source =
        source_parameterization(
            point_site,
            segment_site);
    const auto has_radical_coefficient =
        [](const std::vector<CORE::BigRat>& coefficients) {
            return std::any_of(
                coefficients.begin(),
                coefficients.end(),
                [](const CORE::BigRat& coefficient) {
                    return coefficient != 0;
                });
        };
    if (point_site.x.radical != 0
        || point_site.y.radical != 0
        || has_radical_coefficient(source.x_radical)
        || has_radical_coefficient(source.y_radical)) {
        throw NonRationalParabolaClearanceError(
            "source parabola clearance requires a rational point site");
    }
    return point_clearance_boundary(
        {
            source.x_rational,
            source.y_rational,
            std::nullopt,
            std::nullopt,
        },
        point_site.x.rational,
        point_site.y.rational,
        radius_squared);
}

ClearanceRootBoundary2
parallel_segment_clearance_boundary(
    const RationalPrimitiveParameterization2& primitive,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const CORE::BigRat& radius_squared)
{
    if (radius_squared < 0) {
        throw NegativeClearanceRadiusSquaredError(
            "parallel S-S squared clearance radius is negative");
    }
    if (primitive.x_coefficients.size() != 2
        || primitive.y_coefficients.size() != 2) {
        throw InvalidRationalPrimitiveError(
            "parallel S-S clearance requires an affine chart");
    }
    const RationalPrimitiveParameterization2 canonical =
        parallel_segment_bisector_parameterization(
            first_segment,
            second_segment);
    if (primitive.x_coefficients
            != canonical.x_coefficients
        || primitive.y_coefficients
            != canonical.y_coefficients) {
        throw MismatchedParallelSegmentClearanceError(
            "parallel S-S clearance chart is not canonical");
    }
    const auto squared_distance_to_support =
        [&primitive](
            const MatExactOpenSegmentSource2& segment)
            -> CORE::BigRat {
            const CORE::BigRat line_constant =
                segment.line_a
                    * primitive.x_coefficients[0]
                + segment.line_b
                    * primitive.y_coefficients[0]
                + segment.line_c;
            const CORE::BigRat line_direction =
                segment.line_a
                    * primitive.x_coefficients[1]
                + segment.line_b
                    * primitive.y_coefficients[1];
            if (line_direction != 0) {
                throw NonconstantParallelSegmentClearanceError(
                    "parallel S-S chart is not parallel to a source support");
            }
            const CORE::BigRat line_norm =
                segment.line_a * segment.line_a
                + segment.line_b * segment.line_b;
            if (line_norm == 0) {
                throw InvalidRationalPrimitiveError(
                    "parallel S-S source support has zero normal");
            }
            return line_constant * line_constant / line_norm;
        };
    const CORE::BigRat first_distance =
        squared_distance_to_support(first_segment);
    const CORE::BigRat second_distance =
        squared_distance_to_support(second_segment);
    if (first_distance != second_distance) {
        throw MismatchedParallelSegmentClearanceError(
            "parallel S-S chart is not equidistant from both sources");
    }
    const CORE::BigRat clearance =
        first_distance - radius_squared;
    return {
        CGAL::sign(clearance),
        {},
        {},
    };
}

ClearanceRootBoundary2
nonparallel_segment_clearance_boundary(
    const NonparallelSegmentBisectorParameterization2& primitive,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const CORE::BigRat& radius_squared)
{
    if (radius_squared < 0) {
        throw NegativeClearanceRadiusSquaredError(
            "nonparallel S-S squared clearance radius is negative");
    }
    if (primitive.x_rational.size() != 2
        || primitive.x_radical.size() != 2
        || primitive.y_rational.size() != 2
        || primitive.y_radical.size() != 2) {
        throw InvalidRationalPrimitiveError(
            "nonparallel S-S clearance requires an affine field chart");
    }
    if (primitive.radicand <= 0) {
        throw InvalidQuadraticFieldRadicandError(
            "nonparallel S-S clearance radicand is not positive");
    }

    const MatExactOpenSegmentSource2* ordered_first =
        &first_segment;
    const MatExactOpenSegmentSource2* ordered_second =
        &second_segment;
    if (ordered_second->stable_site_id
        < ordered_first->stable_site_id) {
        std::swap(ordered_first, ordered_second);
    }
    if (ordered_first->stable_site_id
            != primitive.first_segment_id
        || ordered_second->stable_site_id
            != primitive.second_segment_id) {
        throw MismatchedNonparallelSegmentClearanceError(
            "nonparallel S-S clearance sources do not match the chart");
    }

    const auto squared_distance_polynomials =
        [&primitive](
            const MatExactOpenSegmentSource2& segment) {
            const CORE::BigRat normal_squared =
                segment.line_a * segment.line_a
                + segment.line_b * segment.line_b;
            if (normal_squared == 0) {
                throw DegenerateNonparallelSegmentClearanceError(
                    "nonparallel S-S support has zero normal");
            }
            const CORE::BigRat constant_rational =
                segment.line_a * primitive.x_rational[0]
                + segment.line_b * primitive.y_rational[0]
                + segment.line_c;
            const CORE::BigRat constant_radical =
                segment.line_a * primitive.x_radical[0]
                + segment.line_b * primitive.y_radical[0];
            const CORE::BigRat direction_rational =
                segment.line_a * primitive.x_rational[1]
                + segment.line_b * primitive.y_rational[1];
            const CORE::BigRat direction_radical =
                segment.line_a * primitive.x_radical[1]
                + segment.line_b * primitive.y_radical[1];
            RationalPolynomial rational{
                (
                    constant_rational * constant_rational
                    + primitive.radicand
                        * constant_radical
                        * constant_radical)
                    / normal_squared,
                CORE::BigRat(2)
                    * (
                        constant_rational
                            * direction_rational
                        + primitive.radicand
                            * constant_radical
                            * direction_radical)
                    / normal_squared,
                (
                    direction_rational * direction_rational
                    + primitive.radicand
                        * direction_radical
                        * direction_radical)
                    / normal_squared,
            };
            RationalPolynomial radical{
                CORE::BigRat(2)
                    * constant_rational
                    * constant_radical
                    / normal_squared,
                CORE::BigRat(2)
                    * (
                        constant_rational
                            * direction_radical
                        + constant_radical
                            * direction_rational)
                    / normal_squared,
                CORE::BigRat(2)
                    * direction_rational
                    * direction_radical
                    / normal_squared,
            };
            trim(rational);
            trim(radical);
            return std::pair<
                RationalPolynomial,
                RationalPolynomial>{
                std::move(rational),
                std::move(radical),
            };
        };
    auto [first_rational, first_radical] =
        squared_distance_polynomials(*ordered_first);
    auto [second_rational, second_radical] =
        squared_distance_polynomials(*ordered_second);
    const auto has_nonzero =
        [](const RationalPolynomial& polynomial) {
            return std::any_of(
                polynomial.begin(),
                polynomial.end(),
                [](const CORE::BigRat& coefficient) {
                    return coefficient != 0;
                });
        };
    if (has_nonzero(first_radical)
        || has_nonzero(second_radical)) {
        throw NonrationalNonparallelSegmentClearanceError(
            "nonparallel S-S squared clearance retained a radical");
    }
    if (first_rational != second_rational
        || first_rational.size() != 3
        || first_rational[0] != 0
        || first_rational[1] != 0) {
        throw MismatchedNonparallelSegmentClearanceError(
            "nonparallel S-S chart is not equidistant from both sources");
    }
    if (first_rational[2] <= 0) {
        throw DegenerateNonparallelSegmentClearanceError(
            "nonparallel S-S squared-clearance coefficient is not positive");
    }

    const CORE::BigRat determinant =
        ordered_first->line_a * ordered_second->line_b
        - ordered_second->line_a * ordered_first->line_b;
    const CORE::BigRat second_norm =
        ordered_second->line_a * ordered_second->line_a
        + ordered_second->line_b * ordered_second->line_b;
    const CORE::BigRat expected_quadratic =
        second_norm * determinant * determinant;
    if (first_rational[2] != expected_quadratic) {
        throw MismatchedNonparallelSegmentClearanceError(
            "nonparallel S-S clearance scale is not canonical");
    }

    first_rational[0] -= radius_squared;
    return clearance_boundary_from_polynomial(
        std::move(first_rational),
        nullptr);
}

std::vector<MatAdmissibleComponent2>
maximal_clearance_components(
    const std::string& original_dual_id,
    const RationalPrimitiveParameterization2& primitive,
    const ClearanceRootBoundary2& boundary)
{
    if (original_dual_id.empty()) {
        throw InvalidRationalPrimitiveError(
            "original dual identity is empty");
    }

    ExactAlgebraicKernel1 kernel;
    TaggedEndpoint2 lower{
        domain_endpoint(
            original_dual_id,
            "lower",
            primitive.domain_lower,
            kernel),
        false,
    };
    TaggedEndpoint2 upper{
        domain_endpoint(
            original_dual_id,
            "upper",
            primitive.domain_upper,
            kernel),
        false,
    };

    if (boundary.constant_sign.has_value()) {
        if (*boundary.constant_sign == CGAL::NEGATIVE) {
            return {};
        }
        return {
            {
                original_dual_id + "/component-0",
                lower.endpoint,
                upper.endpoint,
            },
        };
    }
    if (boundary.primitive_coefficients.empty()) {
        throw InvalidRationalPrimitiveError(
            "nonconstant clearance boundary has no polynomial");
    }

    std::vector<TaggedEndpoint2> endpoints{lower};
    for (std::size_t index = 0;
         index < boundary.roots.size();
         ++index) {
        append_root_endpoint(
            endpoints,
            upper,
            boundary.roots[index].parameter,
            boundary.roots[index].event_id,
            kernel);
    }
    if (endpoints.front().endpoint.parameter.has_value()
        && upper.endpoint.parameter.has_value()
        && kernel.compare_1_object()(
               *endpoints.front().endpoint.parameter,
               *upper.endpoint.parameter)
            == CGAL::EQUAL) {
        endpoints.front().endpoint.provenance_ids.insert(
            endpoints.front().endpoint.provenance_ids.end(),
            upper.endpoint.provenance_ids.begin(),
            upper.endpoint.provenance_ids.end());
    } else {
        endpoints.push_back(upper);
    }

    if (endpoints.size() == 1) {
        if (clearance_sign_at(
                boundary.primitive_coefficients,
                *endpoints.front().endpoint.parameter,
                kernel)
            == CGAL::NEGATIVE) {
            return {};
        }
        return {
            {
                original_dual_id + "/component-0",
                endpoints.front().endpoint,
                endpoints.front().endpoint,
            },
        };
    }

    std::vector<bool> retained_cells;
    retained_cells.reserve(endpoints.size() - 1);
    for (std::size_t index = 0;
         index + 1 < endpoints.size();
         ++index) {
        retained_cells.push_back(
            clearance_sign_on_open_cell(
                boundary.primitive_coefficients,
                endpoints[index].endpoint,
                endpoints[index + 1].endpoint,
                kernel)
            == CGAL::POSITIVE);
    }

    struct OrderedComponent2 {
        std::size_t order;
        MatParameterEndpoint2 lower;
        MatParameterEndpoint2 upper;
    };
    std::vector<OrderedComponent2> ordered;
    for (std::size_t index = 0;
         index < retained_cells.size();) {
        if (!retained_cells[index]) {
            ++index;
            continue;
        }
        const std::size_t first = index;
        while (index + 1 < retained_cells.size()
               && retained_cells[index + 1]) {
            ++index;
        }
        ordered.push_back(
            {
                2 * first,
                endpoints[first].endpoint,
                endpoints[index + 1].endpoint,
            });
        ++index;
    }
    for (std::size_t index = 0;
         index < endpoints.size();
         ++index) {
        if (!endpoints[index].clearance_root) {
            continue;
        }
        const bool retained_left =
            index > 0 && retained_cells[index - 1];
        const bool retained_right =
            index < retained_cells.size()
            && retained_cells[index];
        if (!retained_left && !retained_right) {
            ordered.push_back(
                {
                    2 * index + 1,
                    endpoints[index].endpoint,
                    endpoints[index].endpoint,
                });
        }
    }
    std::sort(
        ordered.begin(),
        ordered.end(),
        [](const OrderedComponent2& lhs,
           const OrderedComponent2& rhs) {
            return lhs.order < rhs.order;
        });

    std::vector<MatAdmissibleComponent2> components;
    components.reserve(ordered.size());
    for (std::size_t index = 0;
         index < ordered.size();
         ++index) {
        components.push_back(
            {
                original_dual_id + "/component-"
                    + std::to_string(index),
                ordered[index].lower,
                ordered[index].upper,
            });
    }
    return components;
}

std::vector<MatAdmissibleComponent2>
clip_clearance_components_with_domain_roots(
    const std::string& original_dual_id,
    const MatParameterEndpoint2& domain_lower,
    const MatParameterEndpoint2& domain_upper,
    const ClearanceRootBoundary2& boundary,
    const std::function<bool(const CORE::BigRat&)>&
        domain_cell_contains,
    const std::vector<AlgebraicDomainRoot2>& domain_roots)
{
    if (original_dual_id.empty()) {
        throw InvalidRationalPrimitiveError(
            "original dual identity is empty");
    }
    if (!boundary.constant_sign.has_value()
        && boundary.primitive_coefficients.empty()) {
        throw InvalidRationalPrimitiveError(
            "nonconstant clearance boundary has no polynomial");
    }

    struct CombinedEndpoint2 {
        MatParameterEndpoint2 endpoint;
        bool clearance_root;
        bool domain_boundary;
    };
    ExactAlgebraicKernel1 kernel;
    CombinedEndpoint2 lower{
        domain_lower,
        false,
        false,
    };
    CombinedEndpoint2 upper{
        domain_upper,
        false,
        false,
    };
    const auto compare = kernel.compare_1_object();
    const auto is_in_domain =
        [&lower, &upper, &compare](
            const ExactAlgebraicKernel1::Algebraic_real_1&
                parameter) {
            return
                (!lower.endpoint.parameter.has_value()
                 || compare(
                        parameter,
                        *lower.endpoint.parameter)
                     != CGAL::SMALLER)
                && (!upper.endpoint.parameter.has_value()
                    || compare(
                           parameter,
                           *upper.endpoint.parameter)
                        != CGAL::LARGER);
        };
    std::vector<CombinedEndpoint2> candidates;
    for (std::size_t index = 0;
         index < boundary.roots.size();
         ++index) {
        if (!is_in_domain(
                boundary.roots[index].parameter)) {
            continue;
        }
        candidates.push_back(
            {
                {
                    boundary.roots[index].parameter,
                    {
                        boundary.roots[index].event_id,
                    },
                },
                true,
                false,
            });
    }
    for (const AlgebraicDomainRoot2& root : domain_roots) {
        if (!is_in_domain(root.parameter)) {
            continue;
        }
        candidates.push_back(
            {
                {
                    root.parameter,
                    root.provenance_ids,
                },
                false,
                true,
            });
    }
    std::sort(
        candidates.begin(),
        candidates.end(),
        [&compare](const CombinedEndpoint2& lhs,
                   const CombinedEndpoint2& rhs) {
            return compare(
                       *lhs.endpoint.parameter,
                       *rhs.endpoint.parameter)
                == CGAL::SMALLER;
        });

    std::vector<CombinedEndpoint2> endpoints{lower};
    for (CombinedEndpoint2& candidate : candidates) {
        CombinedEndpoint2* equal_endpoint = nullptr;
        if (endpoints.back().endpoint.parameter.has_value()
            && compare(
                   *endpoints.back().endpoint.parameter,
                   *candidate.endpoint.parameter)
                == CGAL::EQUAL) {
            equal_endpoint = &endpoints.back();
        } else if (upper.endpoint.parameter.has_value()
                   && compare(
                          *upper.endpoint.parameter,
                          *candidate.endpoint.parameter)
                       == CGAL::EQUAL) {
            equal_endpoint = &upper;
        }
        if (equal_endpoint == nullptr) {
            endpoints.push_back(std::move(candidate));
            continue;
        }
        equal_endpoint->clearance_root =
            equal_endpoint->clearance_root
            || candidate.clearance_root;
        equal_endpoint->domain_boundary =
            equal_endpoint->domain_boundary
            || candidate.domain_boundary;
        equal_endpoint->endpoint.provenance_ids.insert(
            equal_endpoint->endpoint.provenance_ids.end(),
            candidate.endpoint.provenance_ids.begin(),
            candidate.endpoint.provenance_ids.end());
    }
    if (endpoints.front().endpoint.parameter.has_value()
        && upper.endpoint.parameter.has_value()
        && compare(
               *endpoints.front().endpoint.parameter,
               *upper.endpoint.parameter)
            == CGAL::EQUAL) {
        endpoints.front().endpoint.provenance_ids.insert(
            endpoints.front().endpoint.provenance_ids.end(),
            upper.endpoint.provenance_ids.begin(),
            upper.endpoint.provenance_ids.end());
    } else {
        endpoints.push_back(upper);
    }

    std::vector<bool> inside_cells;
    std::vector<bool> retained_cells;
    inside_cells.reserve(endpoints.size() - 1);
    retained_cells.reserve(endpoints.size() - 1);
    for (std::size_t index = 0;
         index + 1 < endpoints.size();
         ++index) {
        const CORE::BigRat witness = open_cell_witness(
            endpoints[index].endpoint,
            endpoints[index + 1].endpoint,
            kernel);
        const bool inside = domain_cell_contains(witness);
        const bool clearance_admissible =
            boundary.constant_sign.has_value()
            ? *boundary.constant_sign != CGAL::NEGATIVE
            : clearance_sign_on_open_cell(
                  boundary.primitive_coefficients,
                  endpoints[index].endpoint,
                  endpoints[index + 1].endpoint,
                  kernel)
                == CGAL::POSITIVE;
        inside_cells.push_back(inside);
        retained_cells.push_back(
            inside && clearance_admissible);
    }

    struct OrderedComponent2 {
        std::size_t order;
        MatParameterEndpoint2 lower;
        MatParameterEndpoint2 upper;
    };
    std::vector<OrderedComponent2> ordered;
    for (std::size_t index = 0;
         index < retained_cells.size();) {
        if (!retained_cells[index]) {
            ++index;
            continue;
        }
        const std::size_t first = index;
        while (index + 1 < retained_cells.size()
               && retained_cells[index + 1]) {
            ++index;
        }
        ordered.push_back(
            {
                2 * first,
                endpoints[first].endpoint,
                endpoints[index + 1].endpoint,
            });
        ++index;
    }
    for (std::size_t index = 0;
         index < endpoints.size();
         ++index) {
        const bool retained_left =
            index > 0 && retained_cells[index - 1];
        const bool retained_right =
            index < retained_cells.size()
            && retained_cells[index];
        if (retained_left || retained_right) {
            continue;
        }
        const bool inside =
            endpoints[index].domain_boundary
            || (index > 0 && inside_cells[index - 1])
            || (index < inside_cells.size()
                && inside_cells[index]);
        if (!inside) {
            continue;
        }
        const bool clearance_admissible =
            boundary.constant_sign.has_value()
            ? *boundary.constant_sign != CGAL::NEGATIVE
            : clearance_sign_at(
                  boundary.primitive_coefficients,
                  *endpoints[index].endpoint.parameter,
                  kernel)
                != CGAL::NEGATIVE;
        if (inside && clearance_admissible) {
            ordered.push_back(
                {
                    2 * index + 1,
                    endpoints[index].endpoint,
                    endpoints[index].endpoint,
                });
        }
    }
    std::sort(
        ordered.begin(),
        ordered.end(),
        [](const OrderedComponent2& lhs,
           const OrderedComponent2& rhs) {
            return lhs.order < rhs.order;
        });

    std::vector<MatAdmissibleComponent2> components;
    components.reserve(ordered.size());
    for (std::size_t index = 0;
         index < ordered.size();
         ++index) {
        components.push_back(
            {
                original_dual_id + "/component-"
                    + std::to_string(index),
                ordered[index].lower,
                ordered[index].upper,
            });
    }
    return components;
}

std::vector<MatAdmissibleComponent2>
maximal_nonparallel_segment_clearance_components(
    const std::string& original_dual_id,
    const NonparallelSegmentFeatureDomain2& feature_domain,
    const ClearanceRootBoundary2& boundary)
{
    const MatParameterEndpoint2 lower{
        quadratic_field_algebraic_real(
            feature_domain.lower.parameter,
            feature_domain.radicand),
        feature_domain.lower.provenance_ids,
    };
    const MatParameterEndpoint2 upper{
        quadratic_field_algebraic_real(
            feature_domain.upper.parameter,
            feature_domain.radicand),
        feature_domain.upper.provenance_ids,
    };
    ExactAlgebraicKernel1 kernel;
    if (kernel.compare_1_object()(
            *lower.parameter,
            *upper.parameter)
        != CGAL::SMALLER) {
        throw EmptyNonparallelSegmentFeatureDomainError(
            "nonparallel S-S algebraic feature domain is not increasing");
    }
    return clip_clearance_components_with_domain_roots(
        original_dual_id,
        lower,
        upper,
        boundary,
        [](const CORE::BigRat&) {
            return true;
        },
        {});
}

std::vector<MatAdmissibleComponent2>
clip_nonparallel_segment_clearance_components(
    const std::string& original_dual_id,
    const NonparallelSegmentBisectorParameterization2& primitive,
    const NonparallelSegmentFeatureDomain2& feature_domain,
    const ClearanceRootBoundary2& boundary,
    const MatDomainPolygonWithHoles2& domain)
{
    if (primitive.radicand != feature_domain.radicand) {
        throw MismatchedNonparallelSegmentFeatureDomainError(
            "nonparallel S-S chart and feature domain use different fields");
    }
    const MatParameterEndpoint2 lower{
        quadratic_field_algebraic_real(
            feature_domain.lower.parameter,
            feature_domain.radicand),
        feature_domain.lower.provenance_ids,
    };
    const MatParameterEndpoint2 upper{
        quadratic_field_algebraic_real(
            feature_domain.upper.parameter,
            feature_domain.radicand),
        feature_domain.upper.provenance_ids,
    };
    return clip_clearance_components_with_domain_roots(
        original_dual_id,
        lower,
        upper,
        boundary,
        [&domain, &primitive](const CORE::BigRat& parameter) {
            return nonparallel_segment_domain_contains(
                domain,
                primitive,
                parameter);
        },
        nonparallel_segment_domain_roots(
            original_dual_id,
            primitive,
            feature_domain,
            domain));
}

std::vector<MatAdmissibleComponent2>
clip_linear_clearance_components(
    const std::string& original_dual_id,
    const RationalPrimitiveParameterization2& primitive,
    const ClearanceRootBoundary2& boundary,
    const MatDomainPolygonWithHoles2& domain)
{
    ExactAlgebraicKernel1 kernel;
    return clip_clearance_components_with_domain_roots(
        original_dual_id,
        domain_endpoint(
            original_dual_id,
            "lower",
            primitive.domain_lower,
            kernel),
        domain_endpoint(
            original_dual_id,
            "upper",
            primitive.domain_upper,
            kernel),
        boundary,
        [&domain, &primitive](const CORE::BigRat& parameter) {
            return domain_contains(
                domain,
                evaluate_rational_polynomial(
                    primitive.x_coefficients,
                    parameter),
                evaluate_rational_polynomial(
                    primitive.y_coefficients,
                    parameter));
        },
        algebraic_linear_domain_roots(
            original_dual_id,
            primitive,
            domain));
}

std::vector<MatAdmissibleComponent2>
clip_owned_linear_clearance_components(
    const std::string& original_dual_id,
    const RationalPrimitiveParameterization2& primitive,
    const MatParameterEndpoint2& domain_lower,
    const MatParameterEndpoint2& domain_upper,
    const ClearanceRootBoundary2& boundary,
    const MatDomainPolygonWithHoles2& domain)
{
    if (!primitive.domain_lower.has_value()
        || !primitive.domain_upper.has_value()
        || !domain_lower.parameter.has_value()
        || !domain_upper.parameter.has_value()) {
        throw MismatchedLinearCellDomainError(
            "owned linear cell requires bounded matching endpoints");
    }
    if (domain_lower.provenance_ids.empty()
        || domain_upper.provenance_ids.empty()) {
        throw MissingOwnedLinearCellProvenanceError(
            "owned linear cell endpoint provenance is empty");
    }
    ExactAlgebraicKernel1 kernel;
    const auto construct =
        kernel.construct_algebraic_real_1_object();
    const auto compare = kernel.compare_1_object();
    if (compare(
            *domain_lower.parameter,
            construct(*primitive.domain_lower))
            != CGAL::EQUAL
        || compare(
               *domain_upper.parameter,
               construct(*primitive.domain_upper))
            != CGAL::EQUAL
        || compare(
               *domain_lower.parameter,
               *domain_upper.parameter)
            != CGAL::SMALLER) {
        throw MismatchedLinearCellDomainError(
            "owned linear endpoints differ from primitive domain");
    }
    return clip_clearance_components_with_domain_roots(
        original_dual_id,
        domain_lower,
        domain_upper,
        boundary,
        [&domain, &primitive](const CORE::BigRat& parameter) {
            return domain_contains(
                domain,
                evaluate_rational_polynomial(
                    primitive.x_coefficients,
                    parameter),
                evaluate_rational_polynomial(
                    primitive.y_coefficients,
                    parameter));
        },
        algebraic_linear_domain_roots(
            original_dual_id,
            primitive,
            domain));
}

std::vector<MatAdmissibleComponent2>
clip_parabola_clearance_components(
    const std::string& original_dual_id,
    const RationalPrimitiveParameterization2& primitive,
    const ClearanceRootBoundary2& boundary,
    const MatDomainPolygonWithHoles2& domain)
{
    ExactAlgebraicKernel1 kernel;
    return clip_clearance_components_with_domain_roots(
        original_dual_id,
        domain_endpoint(
            original_dual_id,
            "lower",
            primitive.domain_lower,
            kernel),
        domain_endpoint(
            original_dual_id,
            "upper",
            primitive.domain_upper,
            kernel),
        boundary,
        [&domain, &primitive](const CORE::BigRat& parameter) {
            return domain_contains(
                domain,
                evaluate_rational_polynomial(
                    primitive.x_coefficients,
                    parameter),
                evaluate_rational_polynomial(
                    primitive.y_coefficients,
                    parameter));
        },
        parabola_domain_roots(
            original_dual_id,
            primitive,
            domain));
}

std::vector<MatAdmissibleComponent2>
clip_source_parabola_clearance_components(
    const std::string& original_dual_id,
    const MatExactPointSiteSource2& point_site,
    const MatExactOpenSegmentSource2& segment_site,
    const MatParameterEndpoint2& domain_lower,
    const MatParameterEndpoint2& domain_upper,
    const ClearanceRootBoundary2& boundary,
    const MatDomainPolygonWithHoles2& domain)
{
    if (!domain_lower.parameter.has_value()
        || !domain_upper.parameter.has_value()) {
        throw InvalidParabolaCellDomainError(
            "source parabola cell requires two bounded endpoints");
    }
    ExactAlgebraicKernel1 kernel;
    if (kernel.compare_1_object()(
            *domain_lower.parameter,
            *domain_upper.parameter)
        != CGAL::SMALLER) {
        throw InvalidParabolaCellDomainError(
            "source parabola cell endpoints are not strictly ordered");
    }
    if (domain_lower.provenance_ids.empty()
        || domain_upper.provenance_ids.empty()) {
        throw InvalidParabolaCellDomainError(
            "source parabola cell endpoint ownership is empty");
    }
    const SourceParabolaParameterization2 source =
        source_parameterization(
            point_site,
            segment_site);
    return clip_clearance_components_with_domain_roots(
        original_dual_id,
        domain_lower,
        domain_upper,
        boundary,
        [&domain, &source](const CORE::BigRat& parameter) {
            return source_domain_contains(
                domain,
                source,
                parameter);
        },
        source_parabola_domain_roots(
            original_dual_id,
            source,
            domain_lower,
            domain_upper,
            domain));
}
