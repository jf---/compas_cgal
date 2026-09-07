#pragma once

#include "reachable_arrangement_2.h"

class BoundaryContactConstructionError : public std::runtime_error {
public: using std::runtime_error::runtime_error;
};

ReachKernelPoint reachable_kernel_point(const ReachPoint& point);
ReachPoint sample_reachable_boundary(const ReachableBoundaryCurve2& primitive, double parameter);
