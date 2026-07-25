#pragma once

#include "partition_certificate.h"

#include <string>
#include <vector>

std::vector<ParameterChart2> parameter_charts();

VerifiedEventPartition2 verify_chart_coverage(
    const std::vector<ParameterChart2>& charts);

ProjectionRecord2 construct_pullback(
    const std::string& motion_kind,
    const std::vector<std::string>& motion_data,
    const std::string& support_kind,
    const std::vector<std::string>& support_data,
    const std::string& cutter_radius,
    const std::string& center_chart,
    const std::string& rim_chart);

ProjectionRecord2 projection_from_grid(
    const std::string& projection_id,
    const std::vector<std::vector<std::string>>& coefficient_rows,
    const std::string& degree_bound_id);
