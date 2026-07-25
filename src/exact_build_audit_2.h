#pragma once

#include <cstddef>

struct ReachableDomainBuildAudit2 {
    std::size_t geometry_passes = 0;
    std::size_t provenance_arrangements = 0;
    std::size_t center_extractions = 0;
    std::size_t material_batch_unions = 0;
    std::size_t subset_decisions = 0;
    std::size_t residual_differences = 0;
    std::size_t source_geometric_rematches = 0;
    std::size_t input_vertex_count = 0;
    std::size_t ring_rotation_comparisons = 0;
    std::size_t cycle_element_count = 0;
    std::size_t cycle_rotation_comparisons = 0;
    std::size_t selected_face_count = 0;
    std::size_t selected_adjacency_count = 0;
};

struct CoverageTransitionAudit2 {
    std::size_t sweep_constructions = 0;
    std::size_t accumulated_unions = 0;
    std::size_t residual_differences = 0;
    std::size_t replay_constructions = 0;
    std::size_t immutable_region_deep_copies = 0;
};
