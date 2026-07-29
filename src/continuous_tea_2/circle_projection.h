#pragma once

#include "../stock_2.h"
#include "partition_certificate.h"

#include <optional>
#include <string>
#include <vector>

std::string exact_rational_text(
    const Epeck::FT& value);

std::optional<std::string> rational_coordinate_text(
    const GpsPoint::CoordNT& coordinate);

EventPartitionCertificate2
partition_full_circle_boundary_geometry(
    const std::vector<std::string>& motion_data,
    const std::string& cutter_radius,
    const std::string& cap_chord_ratio,
    const std::vector<std::string>& line_sources,
    const std::vector<std::string>& circle_sources,
    const std::vector<std::string>& pair_requests);
