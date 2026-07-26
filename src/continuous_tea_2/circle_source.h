#pragma once

#include "../stock_2.h"
#include "boundary_events.h"
#include "partition_certificate.h"

#include <optional>
#include <string>
#include <vector>

std::string full_circle_stock_identity(
    const std::vector<BoundaryFeatureRecord2>& records);

std::optional<std::string>
full_circle_uniform_disposition(
    const Stock2& stock,
    const std::vector<BoundaryFeatureRecord2>& records,
    double center_x,
    double center_y,
    double phase_dx,
    double phase_dy,
    double tool_radius);

bool has_material_rational_chart_witness(
    const Stock2& stock,
    const Epeck::FT& center_x,
    const Epeck::FT& center_y,
    const Epeck::FT& phase_x,
    const Epeck::FT& phase_y,
    const Epeck::FT& tool_radius);

EventPartitionCertificate2
construct_full_circle_uniform_partition(
    const std::string& stock_identity,
    const std::string& disposition);

EventPartitionCertificate2
construct_full_circle_boundary_pullback_partition(
    const std::string& stock_identity,
    const std::vector<std::string>& motion_data,
    const std::string& cutter_radius,
    const std::string& cap_chord_ratio,
    const std::vector<std::string>& line_sources,
    const std::vector<std::string>& circle_sources);
