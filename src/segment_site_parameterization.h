#pragma once

#include "exact_algebraic_1.h"

#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Cartesian.h>
#include <CGAL/CORE/BigRat.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Parabola_segment_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Segment_Delaunay_graph_traits_2.h>

using MatKernel =
    CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;
using MatTraits =
    CGAL::Segment_Delaunay_graph_traits_without_intersections_2<
        MatKernel,
        CGAL::Field_with_sqrt_tag>;
using SegmentSiteParabola2 =
    CGAL::Parabola_segment_2<MatTraits>;
using MatDomainKernel2 = CGAL::Cartesian<CORE::BigRat>;
using MatDomainPolygon2 = CGAL::Polygon_2<MatDomainKernel2>;
using MatDomainPolygonWithHoles2 =
    CGAL::Polygon_with_holes_2<MatDomainKernel2>;

struct MatParameterDomain2 {
    std::optional<MatTraits::FT> lower;
    std::optional<MatTraits::FT> upper;
};

MatParameterDomain2
exact_parameter_domain(const MatTraits::Line_2& line);

MatParameterDomain2
exact_parameter_domain(const MatTraits::Ray_2& ray);

MatParameterDomain2
exact_parameter_domain(const MatTraits::Segment_2& segment);

MatParameterDomain2
exact_parameter_domain(const SegmentSiteParabola2& parabola);

struct RationalPrimitiveParameterization2 {
    std::vector<CORE::BigRat> x_coefficients;
    std::vector<CORE::BigRat> y_coefficients;
    std::optional<CORE::BigRat> domain_lower;
    std::optional<CORE::BigRat> domain_upper;
};

struct MatQuadraticFieldValue2 {
    CORE::BigRat rational;
    CORE::BigRat radical;
};

class InvalidQuadraticFieldRadicandError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class InvalidQuadraticFieldEmbeddingError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

CGAL::Sign quadratic_field_sign(
    const MatQuadraticFieldValue2& value,
    const CORE::BigRat& radicand);

CGAL::Comparison_result quadratic_field_compare(
    const MatQuadraticFieldValue2& lhs,
    const MatQuadraticFieldValue2& rhs,
    const CORE::BigRat& radicand);

ExactAlgebraicKernel1::Algebraic_real_1
quadratic_field_algebraic_real(
    const MatQuadraticFieldValue2& value,
    const CORE::BigRat& radicand);

struct MatExactPointSiteSource2 {
    std::string stable_site_id;
    MatQuadraticFieldValue2 x;
    MatQuadraticFieldValue2 y;
    CORE::BigRat radicand;
};

class MatExactOpenSegmentSource2 {
public:
    const std::string stable_site_id;
    const CORE::BigRat line_a;
    const CORE::BigRat line_b;
    const CORE::BigRat line_c;

private:
    MatExactOpenSegmentSource2(
        std::string stable_site_id,
        CORE::BigRat line_a,
        CORE::BigRat line_b,
        CORE::BigRat line_c)
        : stable_site_id(std::move(stable_site_id)),
          line_a(std::move(line_a)),
          line_b(std::move(line_b)),
          line_c(std::move(line_c))
    {
    }

    friend MatExactOpenSegmentSource2
    canonical_open_segment_source(
        std::string stable_site_id,
        const CORE::BigRat& source_x,
        const CORE::BigRat& source_y,
        const CORE::BigRat& target_x,
        const CORE::BigRat& target_y);
};

class EmptyOpenSegmentSourceIdentityError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class DegenerateOpenSegmentSourceError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class InvalidRationalPrimitiveError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class OverlappingDomainBoundaryError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class DuplicateOpenSegmentSourceIdentityError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class NonparallelSegmentSupportsError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class CoincidentSegmentSupportsError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class EmptyParallelSegmentFeatureDomainError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class ParallelSegmentSupportsError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class DegenerateLiveSegmentPrimitiveError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class UnboundNonparallelSegmentBranchError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class MismatchedNonparallelSegmentSourceError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class OffSupportSegmentEndpointError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class DegenerateNonparallelFeatureProjectionError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class MissingNonparallelSegmentFeatureProvenanceError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class EmptyNonparallelSegmentFeatureDomainError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

MatExactOpenSegmentSource2 canonical_open_segment_source(
    std::string stable_site_id,
    const CORE::BigRat& source_x,
    const CORE::BigRat& source_y,
    const CORE::BigRat& target_x,
    const CORE::BigRat& target_y);

RationalPrimitiveParameterization2
parallel_segment_bisector_parameterization(
    const MatExactOpenSegmentSource2& first,
    const MatExactOpenSegmentSource2& second);

CORE::BigRat parallel_segment_tangent_parameter(
    const RationalPrimitiveParameterization2& primitive,
    const CORE::BigRat& x,
    const CORE::BigRat& y);

class NonparallelSegmentBisectorParameterization2 {
public:
    const std::string first_segment_id;
    const std::string second_segment_id;
    const int branch_sign;
    const std::vector<CORE::BigRat> x_rational;
    const std::vector<CORE::BigRat> x_radical;
    const std::vector<CORE::BigRat> y_rational;
    const std::vector<CORE::BigRat> y_radical;
    const CORE::BigRat radicand;

private:
    NonparallelSegmentBisectorParameterization2(
        std::string first_segment_id,
        std::string second_segment_id,
        int branch_sign,
        std::vector<CORE::BigRat> x_rational,
        std::vector<CORE::BigRat> x_radical,
        std::vector<CORE::BigRat> y_rational,
        std::vector<CORE::BigRat> y_radical,
        CORE::BigRat radicand)
        : first_segment_id(std::move(first_segment_id)),
          second_segment_id(std::move(second_segment_id)),
          branch_sign(branch_sign),
          x_rational(std::move(x_rational)),
          x_radical(std::move(x_radical)),
          y_rational(std::move(y_rational)),
          y_radical(std::move(y_radical)),
          radicand(std::move(radicand))
    {
    }

    friend NonparallelSegmentBisectorParameterization2
    nonparallel_segment_bisector_parameterization(
        const MatExactOpenSegmentSource2& first,
        const MatExactOpenSegmentSource2& second,
        const MatTraits::Segment_2& live_primitive);
};

NonparallelSegmentBisectorParameterization2
nonparallel_segment_bisector_parameterization(
    const MatExactOpenSegmentSource2& first,
    const MatExactOpenSegmentSource2& second,
    const MatTraits::Segment_2& live_primitive);

MatQuadraticFieldValue2
nonparallel_segment_tangent_parameter(
    const NonparallelSegmentBisectorParameterization2& primitive,
    const MatExactOpenSegmentSource2& segment,
    const CORE::BigRat& endpoint_x,
    const CORE::BigRat& endpoint_y);

struct MatQuadraticFieldDomainBoundary2 {
    MatQuadraticFieldValue2 parameter;
    std::vector<std::string> provenance_ids;
};

class NonparallelSegmentFeatureDomain2 {
public:
    const MatQuadraticFieldDomainBoundary2 lower;
    const MatQuadraticFieldDomainBoundary2 upper;
    const CORE::BigRat radicand;

private:
    NonparallelSegmentFeatureDomain2(
        MatQuadraticFieldDomainBoundary2 lower,
        MatQuadraticFieldDomainBoundary2 upper,
        CORE::BigRat radicand)
        : lower(std::move(lower)),
          upper(std::move(upper)),
          radicand(std::move(radicand))
    {
    }

    friend NonparallelSegmentFeatureDomain2
    intersect_nonparallel_segment_feature_domains(
        MatQuadraticFieldDomainBoundary2 first_source,
        MatQuadraticFieldDomainBoundary2 first_target,
        MatQuadraticFieldDomainBoundary2 second_source,
        MatQuadraticFieldDomainBoundary2 second_target,
        const CORE::BigRat& radicand);
};

NonparallelSegmentFeatureDomain2
intersect_nonparallel_segment_feature_domains(
    MatQuadraticFieldDomainBoundary2 first_source,
    MatQuadraticFieldDomainBoundary2 first_target,
    MatQuadraticFieldDomainBoundary2 second_source,
    MatQuadraticFieldDomainBoundary2 second_target,
    const CORE::BigRat& radicand);

using RationalPolynomial = std::vector<CORE::BigRat>;

void trim(RationalPolynomial& polynomial);
RationalPolynomial multiply(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs);
void add_in_place(
    RationalPolynomial& target,
    const RationalPolynomial& source);
std::vector<ExactAlgebraicInteger1>
primitive_integer_coefficients(
    const RationalPolynomial& polynomial);
bool root_is_in_domain(
    const ExactAlgebraicKernel1::Algebraic_real_1& root,
    const RationalPrimitiveParameterization2& primitive,
    const ExactAlgebraicKernel1& kernel);
CORE::BigRat evaluate_rational_polynomial(
    const std::vector<CORE::BigRat>& coefficients,
    const CORE::BigRat& parameter);

struct SourceParabolaParameterization2 {
    std::vector<CORE::BigRat> x_rational;
    std::vector<CORE::BigRat> x_radical;
    std::vector<CORE::BigRat> y_rational;
    std::vector<CORE::BigRat> y_radical;
    CORE::BigRat radicand;
};

SourceParabolaParameterization2 source_parameterization(
    const MatExactPointSiteSource2& point,
    const MatExactOpenSegmentSource2& segment);

struct RadicalEquation2 {
    ExactAlgebraicKernel2::Polynomial_2 rational;
    ExactAlgebraicKernel2::Polynomial_2 radical;
};

RadicalEquation2 radical_equation(
    const std::vector<CORE::BigRat>& rational,
    const std::vector<CORE::BigRat>& radical,
    const CORE::BigRat& edge_origin,
    const CORE::BigRat& edge_direction);
ExactAlgebraicKernel2::Polynomial_2 radical_norm(
    const RadicalEquation2& equation,
    const CORE::BigRat& radicand);
bool radical_equation_holds(
    const RadicalEquation2& equation,
    const ExactAlgebraicKernel2::Algebraic_real_2& point,
    const ExactAlgebraicKernel2& kernel);
bool source_domain_contains(
    const MatDomainPolygonWithHoles2& domain,
    const SourceParabolaParameterization2& primitive,
    const CORE::BigRat& parameter);
bool domain_contains(
    const MatDomainPolygonWithHoles2& domain,
    const CORE::BigRat& x,
    const CORE::BigRat& y);
