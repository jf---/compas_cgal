#pragma once

#include "segment_site_mat.h"

#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

class InvalidMatNeckCutGraphError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class UnknownMatNeckCutTargetError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class InvalidMatNeckPlateauError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class NonSeparatingMatNeckCandidateError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class MatNeckSeparatingCut2 {
public:
    const std::vector<std::vector<std::string>>&
    edge_partitions() const noexcept;

    bool operator==(
        const MatNeckSeparatingCut2&) const = default;

private:
    explicit MatNeckSeparatingCut2(
        std::vector<std::vector<std::string>>
            edge_partitions);

    std::vector<std::vector<std::string>>
        edge_partitions_;

    friend MatNeckSeparatingCut2
    strict_edge_neck_separating_cut(
        const MatExactGraph2& graph,
        const std::string& edge_id);
    friend MatNeckSeparatingCut2
    vertex_neck_separating_cut(
        const MatExactGraph2& graph,
        const std::string& node_id);
    friend MatNeckSeparatingCut2
    plateau_neck_separating_cut(
        const MatExactGraph2& graph,
        const std::vector<std::string>&
            plateau_node_ids,
        const std::vector<std::string>&
            plateau_edge_ids);
    friend class MatNeckCutIndex2;
};

class MatNeckCutIndex2 {
public:
    static MatNeckCutIndex2 build(
        const MatExactGraph2& graph);
    static MatNeckCutIndex2 build(
        MatExactGraph2&& graph) = delete;

    MatNeckSeparatingCut2 strict_edge_cut(
        const std::string& edge_id) const;
    MatNeckSeparatingCut2 vertex_cut(
        const std::string& node_id) const;
    MatNeckSeparatingCut2 plateau_cut(
        const std::vector<std::string>&
            plateau_node_ids,
        const std::vector<std::string>&
            plateau_edge_ids) const;

private:
    MatNeckCutIndex2(
        const MatExactGraph2& graph,
        std::vector<
            std::array<std::size_t, 2>>
            edge_nodes,
        std::vector<std::vector<std::size_t>>
            incident_edges);

    std::size_t edge_index(
        const std::string& edge_id) const;
    std::size_t node_index(
        const std::string& node_id) const;
    std::size_t other_node(
        std::size_t edge,
        std::size_t node) const;
    std::vector<bool> connected_component(
        std::size_t start) const;
    MatNeckSeparatingCut2 separating_cut(
        const std::vector<bool>& active_nodes,
        const std::vector<bool>& removed_nodes,
        const std::vector<bool>& removed_edges) const;

    const MatExactGraph2* graph_;
    std::vector<std::array<std::size_t, 2>>
        edge_nodes_;
    std::vector<std::vector<std::size_t>>
        incident_edges_;
};

MatNeckSeparatingCut2
strict_edge_neck_separating_cut(
    const MatExactGraph2& graph,
    const std::string& edge_id);

MatNeckSeparatingCut2
vertex_neck_separating_cut(
    const MatExactGraph2& graph,
    const std::string& node_id);

MatNeckSeparatingCut2
plateau_neck_separating_cut(
    const MatExactGraph2& graph,
    const std::vector<std::string>&
        plateau_node_ids,
    const std::vector<std::string>&
        plateau_edge_ids);
