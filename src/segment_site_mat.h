#pragma once

#include "segment_site_clipping.h"
#include "segment_site_endpoint_binding.h"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

struct NormalizedPointSource2 {
    std::string stable_site_id;
    CORE::BigRat x;
    CORE::BigRat y;
};

class InsufficientPointSitesError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class EmptyPointSiteIdentityError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class DuplicatePointSiteIdentityError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class CoincidentPointSitesError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class NegativeSquaredRadiusError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class DegeneratePointSiteTopologyError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class AmbiguousCompositeSiteOwnerError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class IncompleteCompositeSegmentChainError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class UnsupportedCompositeSegmentPrimitiveError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct SegmentSiteMatCompileEvidence2 {
    bool delaunay_valid;
    std::size_t assigned_dual_primitives;
    std::size_t exact_clearance_roots;
    std::size_t matched_generator_sites;
    std::size_t delaunay_vertices;
};

struct MatExactGraphNode2 {
    std::string node_id;
    std::vector<std::string> provenance_ids;
    std::vector<std::string> generator_site_ids;
    std::vector<std::string> parent_site_ids;
};

struct MatExactGraphEdge2 {
    std::string edge_id;
    std::string primitive_kind;
    std::string original_dual_id;
    std::string source_node_id;
    std::string target_node_id;
    MatParameterEndpoint2 source_endpoint;
    MatParameterEndpoint2 target_endpoint;
    std::vector<std::string> generator_site_ids;
    std::vector<std::string> parent_site_ids;
};

struct MatExactGraph2 {
    std::vector<MatExactGraphNode2> nodes;
    std::vector<MatExactGraphEdge2> edges;
    std::size_t rejected_incident_transitions;
    std::size_t matched_generator_sites;
};

SegmentSiteMatCompileEvidence2
segment_site_mat_compile_spike();

MatExactGraph2 exact_point_site_graph(
    const std::vector<NormalizedPointSource2>& points,
    const MatDomainPolygonWithHoles2& domain,
    const CORE::BigRat& radius_squared);

MatExactGraph2 segment_site_live_graph_spike();

MatExactGraph2 segment_site_generic_graph_spike();

MatExactGraph2
segment_site_reversed_source_graph_spike();

MatExactGraph2
segment_site_true_radius_graph_spike();

MatExactGraph2
segment_site_segment_segment_graph_spike();

MatExactGraph2
segment_site_segment_segment_graph_spike(
    const CORE::BigRat& radius_squared);

MatExactGraph2
segment_site_reversed_segment_segment_graph_spike();

MatExactGraph2
segment_site_nonparallel_segment_segment_graph_spike();

MatExactGraph2
segment_site_nonparallel_segment_segment_graph_spike(
    const CORE::BigRat& radius_squared);

MatExactGraph2
segment_site_reversed_nonparallel_segment_segment_graph_spike();

MatExactGraph2
segment_site_disjoint_parallel_segment_graph_spike();

MatExactGraph2
segment_site_external_limited_parallel_segment_graph_spike();

MatExactGraph2
segment_site_domain_coincident_parallel_segment_graph_spike();
