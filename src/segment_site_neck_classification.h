#pragma once

#include "segment_site_neck_evidence_bytes.h"

#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

class InvalidMatNeckWidthBoundariesError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

std::pair<std::vector<std::int64_t>, std::vector<std::string>>
validate_and_classify_necks_v1(
    const MatClearanceProfileGraph2 &bundle,
    const std::vector<std::string> &neck_evidence,
    const std::vector<std::string> &squared_width_boundaries);
