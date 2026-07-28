#pragma once

#include "segment_site_mat.h"

#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

class CanonicalMatSiteCatalog2;

class NonCanonicalMatGraphNodeOrderError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class NonCanonicalMatNodeSiteOrderError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class MatNodeSiteCsrOverflowError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class NonCanonicalMatGraphEdgeOrderError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class UnknownMatGraphNodeIdentityError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class InvalidMatGraphEdgeRecordError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class InvalidMatEndpointEvidenceError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class UnsupportedMatGraphCurveKindError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class MatNumericGraphTableOverflowError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

inline constexpr std::int64_t
    MAT_ENDPOINT_ORIGINAL_VORONOI_VERTEX = 1;
inline constexpr std::int64_t
    MAT_ENDPOINT_DOMAIN_BOUNDARY = 2;
inline constexpr std::int64_t
    MAT_ENDPOINT_CLEARANCE_ROOT = 4;

inline constexpr std::int64_t MAT_DUAL_POINT_POINT = 0;
inline constexpr std::int64_t MAT_DUAL_POINT_SEGMENT = 1;
inline constexpr std::int64_t MAT_DUAL_SEGMENT_SEGMENT = 2;

struct MatNodeSiteCsr2 {
    std::vector<std::int64_t> node_site_offsets;
    std::vector<std::string> node_site_ids;
};

struct MatNumericNodeSiteCsr2 {
    std::vector<std::int64_t> node_site_offsets;
    std::vector<std::int64_t> node_site_ids;
};

struct MatNumericGraphTable2 {
    std::vector<std::string> node_ids;
    std::vector<std::string> original_dual_ids;
    std::vector<std::array<std::int64_t, 8>> edges;
    std::vector<std::int64_t> node_site_offsets;
    std::vector<std::int64_t> node_site_ids;
    std::vector<std::array<std::int64_t, 3>>
        site_provenance;
    std::vector<std::array<std::int64_t, 2>>
        edge_endpoint_provenance_flags;
    std::vector<std::int64_t>
        endpoint_feature_offsets;
    std::vector<std::array<std::int64_t, 5>>
        endpoint_features;
    std::vector<std::array<std::int64_t, 3>>
        edge_exact_flags;

    bool operator==(
        const MatNumericGraphTable2&) const = default;
};

MatNodeSiteCsr2
node_site_csr(const MatExactGraph2& graph);

MatNumericNodeSiteCsr2
numeric_node_site_csr(
    const MatExactGraph2& graph,
    const CanonicalMatSiteCatalog2& catalog);

MatNumericGraphTable2 numeric_graph_table(
    const MatExactGraph2& graph,
    const CanonicalMatSiteCatalog2& catalog);
