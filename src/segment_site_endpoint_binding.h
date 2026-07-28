#pragma once

#include "segment_site_provenance.h"

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

class CanonicalMatSiteGeometryIndex2;
class CanonicalMatRationalSources2;

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

class AmbiguousParallelSegmentPointLimiterError
    : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class NonRationalParallelSegmentPointLimiterError
    : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class AmbiguousParallelSegmentOpenLimiterError
    : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class MismatchedLiveNonparallelSegmentBridgeError
    : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class UnboundLiveNonparallelSegmentEndpointError
    : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class UnsupportedNonparallelSegmentLimiterError
    : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class AmbiguousNonparallelSegmentOpenLimiterError
    : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class IdenticallyZeroNonparallelSegmentLimiterEquationError
    : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class UnknownCanonicalMatParabolaSourceError
    : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class IncidentCanonicalMatParabolaSourceError
    : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class InvalidCanonicalMatParabolaFeatureDomainError
    : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class UnknownCanonicalMatPointRaySourceError
    : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class MismatchedLivePointPointRayError
    : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class UnboundLivePointPointRayEndpointError
    : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class NonRationalLivePointPointRayEndpointError
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
bind_point_segment_cell_endpoints(
    const CanonicalMatRationalSources2& sources,
    const std::string& focus_id,
    const std::string& segment_id,
    const CanonicalMatSiteGeometryIndex2& site_index,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge);

RationalPrimitiveParameterization2
bind_point_point_ray_parameterization(
    const CanonicalMatRationalSources2& sources,
    const std::vector<std::string>& generator_ids,
    const CanonicalMatSiteGeometryIndex2& site_index,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge);

std::pair<MatParameterEndpoint2, MatParameterEndpoint2>
bind_parallel_segment_segment_cell_endpoints(
    const RationalPrimitiveParameterization2& primitive,
    const RationalDomainRoot2& lower,
    const RationalDomainRoot2& upper,
    const std::vector<std::string>& generator_ids,
    const std::vector<GeneratorSite2>& generators,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge);

std::pair<MatParameterEndpoint2, MatParameterEndpoint2>
bind_nonparallel_segment_segment_cell_endpoints(
    const NonparallelSegmentBisectorParameterization2& primitive,
    const NonparallelSegmentFeatureDomain2& feature_domain,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const std::vector<MatExactOpenSegmentSource2>& segment_limiters,
    const std::vector<std::string>& generator_ids,
    const std::vector<GeneratorSite2>& generators,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge);

std::pair<MatParameterEndpoint2, MatParameterEndpoint2>
bind_nonparallel_segment_segment_cell_endpoints(
    const NonparallelSegmentBisectorParameterization2& primitive,
    const NonparallelSegmentFeatureDomain2& feature_domain,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const std::vector<MatExactOpenSegmentSource2>& segment_limiters,
    const std::vector<std::string>& generator_ids,
    const CanonicalMatSiteGeometryIndex2& site_index,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge);

std::pair<MatParameterEndpoint2, MatParameterEndpoint2>
bind_parallel_segment_segment_cell_endpoints(
    const RationalPrimitiveParameterization2& primitive,
    const RationalDomainRoot2& lower,
    const RationalDomainRoot2& upper,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const std::vector<MatExactPointSiteSource2>& point_limiters,
    const std::vector<MatExactOpenSegmentSource2>& segment_limiters,
    const std::vector<std::string>& generator_ids,
    const std::vector<GeneratorSite2>& generators,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge);

std::pair<MatParameterEndpoint2, MatParameterEndpoint2>
bind_parallel_segment_segment_cell_endpoints(
    const RationalPrimitiveParameterization2& primitive,
    const RationalDomainRoot2& lower,
    const RationalDomainRoot2& upper,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const std::vector<MatExactPointSiteSource2>& point_limiters,
    const std::vector<MatExactOpenSegmentSource2>& segment_limiters,
    const std::vector<std::string>& generator_ids,
    const CanonicalMatSiteGeometryIndex2& site_index,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge);

std::pair<MatParameterEndpoint2, MatParameterEndpoint2>
bind_parallel_segment_segment_cell_endpoints(
    const RationalPrimitiveParameterization2& primitive,
    const RationalDomainRoot2& lower,
    const RationalDomainRoot2& upper,
    const MatExactOpenSegmentSource2& first_segment,
    const MatExactOpenSegmentSource2& second_segment,
    const std::vector<MatExactPointSiteSource2>& point_limiters,
    const std::vector<std::string>& generator_ids,
    const std::vector<GeneratorSite2>& generators,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge);
