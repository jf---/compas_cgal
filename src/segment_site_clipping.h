#pragma once

#include "segment_site_parameterization.h"
#include "segment_site_provenance.h"

#include <optional>
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

ClearanceRootBoundary2 point_clearance_boundary(
    const RationalPrimitiveParameterization2& primitive,
    const CORE::BigRat& site_x,
    const CORE::BigRat& site_y,
    const CORE::BigRat& radius_squared);

std::vector<MatAdmissibleComponent2>
maximal_clearance_components(
    const std::string& original_dual_id,
    const RationalPrimitiveParameterization2& primitive,
    const ClearanceRootBoundary2& boundary);

std::vector<MatAdmissibleComponent2>
clip_linear_clearance_components(
    const std::string& original_dual_id,
    const RationalPrimitiveParameterization2& primitive,
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
    const std::optional<CORE::BigRat>& domain_lower,
    const std::optional<CORE::BigRat>& domain_upper,
    const ClearanceRootBoundary2& boundary,
    const MatDomainPolygonWithHoles2& domain);
