#pragma once

#include "segment_site_parameterization.h"
#include "segment_site_provenance.h"

#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

struct ClearanceRootEvent2 {
    ExactAlgebraicKernel1::Algebraic_real_1 parameter;
    std::string event_id;
};

struct ClearanceRootBoundary2 {
    std::optional<CGAL::Sign> constant_sign;
    std::vector<ExactAlgebraicInteger1> primitive_coefficients;
    std::vector<ClearanceRootEvent2> roots;
};

class InvalidParabolaCellDomainError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class NegativeClearanceRadiusSquaredError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class NonRationalParabolaClearanceError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class NonconstantParallelSegmentClearanceError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class MismatchedParallelSegmentClearanceError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class MismatchedNonparallelSegmentClearanceError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class NonrationalNonparallelSegmentClearanceError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class DegenerateNonparallelSegmentClearanceError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class MismatchedLinearCellDomainError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class MissingOwnedLinearCellProvenanceError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

ClearanceRootBoundary2 point_clearance_boundary(
    const RationalPrimitiveParameterization2& primitive,
    const CORE::BigRat& site_x,
    const CORE::BigRat& site_y,
    const CORE::BigRat& radius_squared);

ClearanceRootBoundary2
source_parabola_clearance_boundary(
    const MatExactPointSiteSource2& point_site,
    const MatExactOpenSegmentSource2& segment_site,
    const CORE::BigRat& radius_squared);

ClearanceRootBoundary2
parallel_segment_clearance_boundary(
    const RationalPrimitiveParameterization2& primitive,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const CORE::BigRat& radius_squared);

ClearanceRootBoundary2
nonparallel_segment_clearance_boundary(
    const NonparallelSegmentBisectorParameterization2& primitive,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const CORE::BigRat& radius_squared);

std::vector<MatAdmissibleComponent2>
maximal_clearance_components(
    const std::string& original_dual_id,
    const RationalPrimitiveParameterization2& primitive,
    const ClearanceRootBoundary2& boundary);

std::vector<MatAdmissibleComponent2>
maximal_nonparallel_segment_clearance_components(
    const std::string& original_dual_id,
    const NonparallelSegmentFeatureDomain2& feature_domain,
    const ClearanceRootBoundary2& boundary);

std::vector<MatAdmissibleComponent2>
clip_nonparallel_segment_clearance_components(
    const std::string& original_dual_id,
    const NonparallelSegmentBisectorParameterization2& primitive,
    const NonparallelSegmentFeatureDomain2& feature_domain,
    const ClearanceRootBoundary2& boundary,
    const MatDomainPolygonWithHoles2& domain);

std::vector<MatAdmissibleComponent2>
clip_linear_clearance_components(
    const std::string& original_dual_id,
    const RationalPrimitiveParameterization2& primitive,
    const ClearanceRootBoundary2& boundary,
    const MatDomainPolygonWithHoles2& domain);

std::vector<MatAdmissibleComponent2>
clip_owned_linear_clearance_components(
    const std::string& original_dual_id,
    const RationalPrimitiveParameterization2& primitive,
    const MatParameterEndpoint2& domain_lower,
    const MatParameterEndpoint2& domain_upper,
    const ClearanceRootBoundary2& boundary,
    const MatDomainPolygonWithHoles2& domain);

std::vector<MatAdmissibleComponent2>
clip_parabola_clearance_components(
    const std::string& original_dual_id,
    const RationalPrimitiveParameterization2& primitive,
    const ClearanceRootBoundary2& boundary,
    const MatDomainPolygonWithHoles2& domain);

std::vector<MatAdmissibleComponent2>
clip_source_parabola_clearance_components(
    const std::string& original_dual_id,
    const MatExactPointSiteSource2& point_site,
    const MatExactOpenSegmentSource2& segment_site,
    const MatParameterEndpoint2& domain_lower,
    const MatParameterEndpoint2& domain_upper,
    const ClearanceRootBoundary2& boundary,
    const MatDomainPolygonWithHoles2& domain);
