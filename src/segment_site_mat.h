#pragma once

#include "segment_site_clipping.h"

#include <cstddef>
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_policies_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_traits_2.h>
#include <CGAL/Voronoi_diagram_2.h>

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

struct SegmentSiteMatCompileEvidence2 {
    bool delaunay_valid;
    std::size_t assigned_dual_primitives;
    std::size_t exact_clearance_roots;
    std::size_t matched_generator_sites;
    std::size_t delaunay_vertices;
};

SegmentSiteMatCompileEvidence2
segment_site_mat_compile_spike();
