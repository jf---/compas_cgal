#pragma once

#include "exact_region_2.h"
#include "reachable_input_2.h"

// World-XY millimetres: circles=(cx, cy, guide radius), segments=(x0,y0,x1,y1),
// disks=(cx,cy). All motions use the separately injected cutter radius in mm.
// Final coverage only; no chronological engagement or continuity assertion.
ExactRegion2 remaining_material(
    const ExactRegion2& target,
    Eigen::Ref<const compas::RowMatrixXd> circles,
    Eigen::Ref<const compas::RowMatrixXd> segments,
    Eigen::Ref<const compas::RowMatrixXd> disks,
    double tool_radius);
