#pragma once

#include "exact_region_2.h"

#include <vector>

ReachPolygon reach_disk_polygon(
    const ReachKernelPoint& center,
    const ReachFT& radius);
std::vector<ReachPolygon> reach_capsule_parts(
    const ReachKernelPoint& start,
    const ReachKernelPoint& end,
    const ReachFT& radius);
std::vector<ReachPolygon> reach_arc_sweep_parts(
    const ReachXCurve& guide_arc,
    const ReachFT& tool_radius);
// The `_into` forms build the sweep in a set the caller already owns, and the
// by-value forms delegate to them. A ReachSet has no move constructor, so a set
// returned by value can only reach long-lived storage through CGAL's Gps copy
// constructor -- which leaves the copy's arrangement reading the traits of the
// returned object, and that object then dies. A caller that must keep the sweep
// alive (ExactRegion2 storage, or a sweep a region's arrangement may end up
// borrowing traits from) builds it with the `_into` form instead.
void reach_join_parts_into(
    ReachSet& target,
    const std::vector<ReachPolygon>& polygons,
    const std::vector<ReachPolygonWithHoles>& polygons_with_holes);
ReachSet reach_join_parts(
    const std::vector<ReachPolygon>& polygons,
    const std::vector<ReachPolygonWithHoles>& polygons_with_holes);
void reach_full_circle_sweep_into(
    ReachSet& target,
    const ReachKernelPoint& center,
    const ReachKernelVector& phase_vector,
    const ReachFT& tool_radius);
ReachSet reach_full_circle_sweep(
    const ReachKernelPoint& center,
    const ReachKernelVector& phase_vector,
    const ReachFT& tool_radius);
