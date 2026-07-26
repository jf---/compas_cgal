#pragma once

#include "segment_site_clipping.h"
#include "segment_site_endpoint_binding.h"

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

struct MatExactGraphNode2 {
    std::string node_id;
    MatParameterEndpoint2 endpoint;
    std::vector<std::string> generator_site_ids;
};

struct MatExactGraphEdge2 {
    std::string edge_id;
    std::string primitive_kind;
    std::string source_node_id;
    std::string target_node_id;
    std::vector<std::string> generator_site_ids;
};

struct MatExactGraph2 {
    std::vector<MatExactGraphNode2> nodes;
    std::vector<MatExactGraphEdge2> edges;
    std::size_t rejected_incident_transitions;
};

SegmentSiteMatCompileEvidence2
segment_site_mat_compile_spike();

MatExactGraph2 segment_site_live_graph_spike();

MatExactGraph2 segment_site_generic_graph_spike();
