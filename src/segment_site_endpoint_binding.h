#pragma once

#include "segment_site_provenance.h"

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

class UnboundLiveParabolaEndpointError : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class MismatchedLiveParabolaBridgeError : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class AmbiguousLiveParabolaEndpointError : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class IdenticallyZeroSegmentLimiterEquationError : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class NonCanonicalQuadraticFieldError : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class CoincidentSegmentEndpointsError : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class ZeroSegmentLineNormalError : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class MismatchedLiveSegmentSegmentBridgeError
    : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class UnboundLiveSegmentSegmentEndpointError
    : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class UnsupportedSegmentSegmentLimiterError
    : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

bool exact_open_segment_feature_contains(
    const SourceParabolaParameterization2& parabola,
    const MatExactPointSiteSource2& segment_source,
    const MatExactPointSiteSource2& segment_target,
    const ExactAlgebraicKernel1::Algebraic_real_1& parameter);

void require_distinct_live_parabola_endpoints(
    const SegmentSiteParabola2& parabola);

MatParameterEndpoint2 bind_point_limiter_parabola_endpoint(
    const MatExactPointSiteSource2& focus,
    const MatExactOpenSegmentSource2& segment,
    const MatExactPointSiteSource2& segment_source,
    const MatExactPointSiteSource2& segment_target,
    const MatExactPointSiteSource2& limiter,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge);

MatParameterEndpoint2 bind_segment_limiter_parabola_endpoint(
    const MatExactPointSiteSource2& focus,
    const MatExactOpenSegmentSource2& segment,
    const MatExactPointSiteSource2& segment_source,
    const MatExactPointSiteSource2& segment_target,
    const MatExactOpenSegmentSource2& limiter,
    const MatExactPointSiteSource2& limiter_source,
    const MatExactPointSiteSource2& limiter_target,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge);

MatParameterEndpoint2 bind_segment_limiter_parabola_endpoint(
    const MatExactPointSiteSource2& focus,
    const MatExactOpenSegmentSource2& segment,
    const MatExactPointSiteSource2& segment_source,
    const MatExactPointSiteSource2& segment_target,
    const MatExactOpenSegmentSource2& limiter,
    const MatExactPointSiteSource2& limiter_source,
    const MatExactPointSiteSource2& limiter_target,
    const SegmentSiteParabola2& live_parabola,
    const MatTraits::Point_2& live_point);

std::pair<MatParameterEndpoint2, MatParameterEndpoint2>
bind_parallel_segment_segment_cell_endpoints(
    const RationalPrimitiveParameterization2& primitive,
    const RationalDomainRoot2& lower,
    const RationalDomainRoot2& upper,
    const std::vector<std::string>& generator_ids,
    const std::vector<GeneratorSite2>& generators,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge);
