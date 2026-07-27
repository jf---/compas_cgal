#pragma once

#include "segment_site_mat.h"

#include <cstddef>
#include <map>
#include <string>
#include <vector>

void append_exact_graph_components(
    const std::string& dual_id,
    const std::string& primitive_kind,
    const std::vector<std::string>& generator_ids,
    const std::vector<MatAdmissibleComponent2>& components,
    MatExactGraph2& graph,
    std::map<std::string, std::size_t>& node_indices);

void append_exact_graph_components(
    const std::string& piece_dual_id,
    const std::string& original_dual_id,
    const std::string& primitive_kind,
    const std::vector<std::string>& generator_ids,
    const std::vector<std::string>& parent_site_ids,
    const std::vector<MatAdmissibleComponent2>& components,
    MatExactGraph2& graph,
    std::map<std::string, std::size_t>& node_indices);
