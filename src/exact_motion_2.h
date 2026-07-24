#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

using Epeck = CGAL::Exact_predicates_exact_constructions_kernel;
using EPoint = Epeck::Point_2;
using EVector = Epeck::Vector_2;
using ECircle = Epeck::Circle_2;
using ESegment = Epeck::Segment_2;

struct ExactSegmentMotion2 {
    EPoint start;
    EPoint end;
};

struct ExactCircleMotion2 {
    EPoint center;
    EVector phase_vector;
    bool clockwise;
};
