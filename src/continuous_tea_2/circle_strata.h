#pragma once

#include "../exact_motion_2.h"
#include "partition_certificate.h"

#include <string>
#include <optional>
#include <vector>

class Stock2;

struct FullCircleAuthorityParameter2 {
    int chart;
    Epeck::FT parameter;
};

struct FullCircleCellAuthority2 {
    ContinuousTeaVerdict verdict;
    std::string whole_rim_disposition;
    std::string canonical_bytes;
    std::string canonical_digest;
    std::optional<FullCircleAuthorityParameter2> violating_parameter;
};

std::vector<std::string>
derive_full_circle_pair_requests(
    const Stock2& stock,
    const EventPartitionCertificate2& topology_partition,
    const Epeck::FT& center_x,
    const Epeck::FT& center_y,
    const Epeck::FT& phase_dx,
    const Epeck::FT& phase_dy,
    const Epeck::FT& tool_radius,
    const Epeck::FT& cap_chord_ratio);

EventPartitionCertificate2
construct_full_circle_pair_closed_partition(
    const Stock2& stock,
    const std::string& stock_identity,
    const std::vector<std::string>& motion_data,
    const std::string& cutter_radius,
    const std::string& cap_chord_ratio,
    const std::vector<std::string>& line_sources,
    const std::vector<std::string>& circle_sources,
    const Epeck::FT& center_x,
    const Epeck::FT& center_y,
    const Epeck::FT& phase_dx,
    const Epeck::FT& phase_dy,
    const Epeck::FT& tool_radius,
    const Epeck::FT& cap_chord_ratio_value);

FullCircleCellAuthority2
construct_full_circle_cell_authority(
    const Stock2& stock,
    const VerifiedEventPartition2& verified_partition,
    const Epeck::FT& center_x,
    const Epeck::FT& center_y,
    const Epeck::FT& phase_dx,
    const Epeck::FT& phase_dy,
    const Epeck::FT& tool_radius,
    const Epeck::FT& cap_chord_ratio);

bool verify_full_circle_cell_authority(
    const Stock2& stock,
    const VerifiedEventPartition2& verified_partition,
    const Epeck::FT& center_x,
    const Epeck::FT& center_y,
    const Epeck::FT& phase_dx,
    const Epeck::FT& phase_dy,
    const Epeck::FT& tool_radius,
    const Epeck::FT& cap_chord_ratio,
    const FullCircleCellAuthority2& candidate);

class IncompleteFullCircleCellAuthorityError
    : public EventSubstrateError {
public:
    using EventSubstrateError::EventSubstrateError;
};
