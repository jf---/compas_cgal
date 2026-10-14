#pragma once

// In-process engagement: a station INSIDE a motion is classified against the
// stock minus that motion's own swept prefix, not against the static stock the
// motion starts from. The static engaged rim arcs are harvested once
// (classify_audit_unguarded_station_exact) and clipped exactly to the part the
// prefix has not removed:
//
//   line  -- the prefix is the capsule from the motion start to the station; on
//            the station's rim it covers exactly the closed rear half-disk
//            (offset . direction <= 0).
//   arc   -- the prefix is the start disk union the annular sector swept from
//            the start ray to the station ray. On the station's rim a point
//            outside the sector is covered iff it lies in the start disk (the
//            squared distance to a guide point is unimodal in the guide angle,
//            so its minimum over the prefix is at an end of the prefix), and a
//            point inside the sector is always covered.
//
// Every covering decision is an exact predicate on one-root rim points; cells
// between exact event points are decided at a rational witness strictly inside
// the cell (constructed by exact-comparison interval refinement, never by a
// search over engagement). The query never mutates the stock: depletion is the
// replay's job.

#include "audit_certification_2.h"

#include <stdexcept>

// Raised when a swept-prefix rim cell cannot be decided: distinct cell
// endpoints sharing an x coordinate on an x-monotone rim piece, a degenerate
// radial direction, or a witness that fails its own exact containment check.
// Each is an invariant breach of the clipping construction, never an input
// condition.
class SweptPrefixRimCellError : public AuditCertificationError {
public:
    using AuditCertificationError::AuditCertificationError;
};

// Raised when an arc station lies exactly on the start ray strictly inside the
// motion, so the exact orientation cannot tell a zero prefix from a full turn
// (the binary64 angle seam rounded the station onto the seam).
class ArcPrefixExtentAmbiguousError : public AuditCertificationError {
public:
    using AuditCertificationError::AuditCertificationError;
};

// How far an arc's prefix has turned at a station, decided exactly from points:
// ZERO and FULL_TURN are structural (parameter 0; parameter 1 of a closed
// circle), the two interior extents come from classify_arc_prefix_extent_exact.
// At exactly a half turn the intersection and union sector tests coincide, so
// AT_MOST_HALF_TURN includes it without a convention.
enum class ArcPrefixExtent2 {
    ZERO,
    AT_MOST_HALF_TURN,
    BEYOND_HALF_TURN,
    FULL_TURN,
};

// Exact extent of the prefix from motion_start to station around guide_center
// in the travel direction: CGAL::orientation of (guide, start, station), negated
// for clockwise travel; on the start line, CGAL::angle decides the exact half
// turn and the start ray itself raises ArcPrefixExtentAmbiguousError.
ArcPrefixExtent2 classify_arc_prefix_extent_exact(
    const EPoint& guide_center,
    const EPoint& motion_start,
    const EPoint& station_center,
    bool clockwise);

AuditExactStationClassification2 classify_line_swept_prefix_station_exact(
    const Stock2& stock,
    const EPoint& motion_start,
    const EPoint& motion_end,
    const Epeck::FT& parameter,
    const Epeck::FT& tool_radius,
    const Epeck::FT& cap_chord_ratio,
    const Epeck::FT& gap_close_ratio);

AuditExactStationClassification2 classify_arc_swept_prefix_station_exact(
    const Stock2& stock,
    const EPoint& guide_center,
    const EPoint& motion_start,
    const EPoint& station_center,
    bool clockwise,
    ArcPrefixExtent2 extent,
    const Epeck::FT& tool_radius,
    const Epeck::FT& cap_chord_ratio,
    const Epeck::FT& gap_close_ratio);
