#pragma once

#include "exact_algebraic_1.h"

#include <cstddef>

#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Parabola_segment_2.h>
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_policies_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_traits_2.h>
#include <CGAL/Voronoi_diagram_2.h>

using MatKernel =
    CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;
using MatTraits =
    CGAL::Segment_Delaunay_graph_traits_without_intersections_2<
        MatKernel,
        CGAL::Field_with_sqrt_tag>;
using SegmentSiteDelaunay2 =
    CGAL::Segment_Delaunay_graph_2<MatTraits>;
using SegmentSiteAdaptationTraits2 =
    CGAL::Segment_Delaunay_graph_adaptation_traits_2<
        SegmentSiteDelaunay2>;
using SegmentSiteAdaptationPolicy2 =
    CGAL::Segment_Delaunay_graph_degeneracy_removal_policy_2<
        SegmentSiteDelaunay2>;
using SegmentSiteVoronoi2 = CGAL::Voronoi_diagram_2<
    SegmentSiteDelaunay2,
    SegmentSiteAdaptationTraits2,
    SegmentSiteAdaptationPolicy2>;
using SegmentSiteParabola2 =
    CGAL::Parabola_segment_2<MatTraits>;

struct SegmentSiteMatCompileEvidence2 {
    bool delaunay_valid;
    std::size_t assigned_dual_primitives;
    std::size_t exact_clearance_roots;
};

SegmentSiteMatCompileEvidence2
segment_site_mat_compile_spike();
