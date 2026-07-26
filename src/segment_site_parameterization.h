#pragma once

#include "exact_algebraic_1.h"

#include <optional>
#include <stdexcept>
#include <string>
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

struct MatExactPointSiteSource2 {
    std::string stable_site_id;
    MatQuadraticFieldValue2 x;
    MatQuadraticFieldValue2 y;
    CORE::BigRat radicand;
};

struct MatExactOpenSegmentSource2 {
    std::string stable_site_id;
    CORE::BigRat line_a;
    CORE::BigRat line_b;
    CORE::BigRat line_c;
};

class InvalidRationalPrimitiveError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class OverlappingDomainBoundaryError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

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
