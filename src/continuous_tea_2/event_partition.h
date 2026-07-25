#pragma once

#include "boundary_events.h"
#include "partition_certificate.h"

#include <vector>

EventPartitionCertificate2 partition_projections(
    const std::vector<ProjectionInput2>& projections);

EventPartitionCertificate2 partition_pullback_overlap(
    const ProjectionRecord2& projection,
    const std::vector<PartitionEvent2>& events);

std::vector<TrimmedLineBranch2> solve_trimmed_line_branches(
    const std::vector<std::string>& line_support,
    const std::vector<std::string>& trim_start,
    const std::vector<std::string>& trim_end,
    const std::vector<std::string>& segment_motion,
    const std::string& cutter_radius,
    const std::string& rim_chart);

ProjectedRegularizationVertex2 project_regularization_vertex(
    const Stock2& stock,
    std::size_t first_index,
    std::size_t second_index,
    const std::string& vertex_id);
