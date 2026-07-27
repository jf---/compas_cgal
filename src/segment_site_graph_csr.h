#pragma once

#include "segment_site_mat.h"

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

struct MatNodeSiteCsr2 {
    std::vector<std::int64_t> node_site_offsets;
    std::vector<std::string> node_site_ids;
};

struct MatNumericNodeSiteCsr2 {
    std::vector<std::int64_t> node_site_offsets;
    std::vector<std::int64_t> node_site_ids;
};

MatNodeSiteCsr2
node_site_csr(const MatExactGraph2& graph);

MatNumericNodeSiteCsr2
numeric_node_site_csr(
    const MatExactGraph2& graph,
    const CanonicalMatSiteCatalog2& catalog);
