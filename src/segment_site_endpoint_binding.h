#pragma once

#include "segment_site_provenance.h"

#include <stdexcept>

struct MatLiveParabolaEndpoint2 {
  MatParameterEndpoint2 endpoint;
  bool matches_live_p1_or_p2;
};

class UnboundLiveParabolaEndpointError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

MatLiveParabolaEndpoint2 bind_point_limiter_parabola_endpoint(
    const MatExactPointSiteSource2 &focus,
    const MatExactOpenSegmentSource2 &segment,
    const MatExactPointSiteSource2 &segment_source,
    const MatExactPointSiteSource2 &segment_target,
    const MatExactPointSiteSource2 &limiter,
    const MatTraits::Point_2 &live_point,
    const SegmentSiteParabola2 &live_parabola);
