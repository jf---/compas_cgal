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
ReachSet reach_join_parts(
    const std::vector<ReachPolygon>& polygons,
    const std::vector<ReachPolygonWithHoles>& polygons_with_holes);
ReachSet reach_full_circle_sweep(
    const ReachKernelPoint& center,
    const ReachKernelVector& phase_vector,
    const ReachFT& tool_radius);
