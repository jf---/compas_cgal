#pragma once

// The stations of one machining circle, shared by the arrangement engine
// (engagement_circle_prefixes, engagement_circle_exceeded) and the rim engine
// (AnnulusView2), so both decide the same exact centres and prefix extents.

#include "audit_certification_2.h"
#include "engagement_2.h"
#include "exact_motion_2.h"
#include "swept_prefix_engagement_2.h"

#include <cstddef>
#include <vector>

// One station of a swept-prefix motion: its exact centre and prefix extent,
// placed BEFORE any classification so that a decision-only scan raises on
// exactly the placement errors the full prefix call raises on.
struct ArcPrefixStation2 {
    EPoint centre;
    ArcPrefixExtent2 extent;
};

struct CirclePrefixStations2 {
    EPoint guide;
    EPoint start;
    bool clockwise;
    std::vector<ArcPrefixStation2> stations;
};

// Validates and places every station exactly as engagement_circle_prefixes does:
// parameter 0 is the authored entry, parameters > 0 lie on the authored radius.
CirclePrefixStations2 circle_prefix_stations(
    double guide_cx, double guide_cy, double guide_radius,
    double entry_x, double entry_y, bool clockwise,
    const std::vector<double>& parameters, double tool_radius,
    double cap_chord_ratio, double gap_close_ratio);

// The probing order of the decision-only queries (engagement_2.cpp).
std::vector<std::size_t> station_probe_order(std::size_t count);

// The reporting view of an exact classification at an exact station centre.
EngagementSample engagement_sample_at_station(
    const AuditExactStationClassification2& classification,
    const EPoint& centre);

// An exact rational unit vector steered by the double direction (x, y) != 0. The
// double is a construction parameter, never a decision.
EVector exact_unit_along(double x, double y);

// That unit's construction: ((1 - t^2), 2 t) / (1 + t^2), negated when reflected.
// The double t IS the exact rational.
struct RationalUnitSteering2 {
    double t;
    bool reflected;
};
RationalUnitSteering2 rational_unit_steering(double x, double y);

// The exact rotation a circle's station at 0 < parameter < 1 applies to the motion's
// phase: the unit vector that phase (1, 0) rotates to.
EVector circle_station_rotation(double parameter, bool clockwise);
