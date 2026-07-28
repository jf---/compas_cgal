#pragma once

#include "partition_certificate.h"

#include <string>

class Stock2;

struct FullCircleCellAuthority2 {
    ContinuousTeaVerdict verdict;
    std::string whole_rim_disposition;
    std::string canonical_bytes;
    std::string canonical_digest;
};

FullCircleCellAuthority2
construct_full_circle_cell_authority(
    const Stock2& stock,
    const VerifiedEventPartition2& verified_partition,
    double center_x,
    double center_y,
    double phase_dx,
    double phase_dy,
    double tool_radius,
    double cap_chord_ratio);

bool verify_full_circle_cell_authority(
    const Stock2& stock,
    const VerifiedEventPartition2& verified_partition,
    double center_x,
    double center_y,
    double phase_dx,
    double phase_dy,
    double tool_radius,
    double cap_chord_ratio,
    const FullCircleCellAuthority2& candidate);

class IncompleteFullCircleCellAuthorityError
    : public EventSubstrateError {
public:
    using EventSubstrateError::EventSubstrateError;
};
