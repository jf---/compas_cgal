#pragma once

#include "partition_certificate.h"

#include <string>
#include <vector>

struct CcwOrientation2 {
    std::string disposition;
    std::string determinant_sign;
    std::string dot_sign;
};

CcwOrientation2 classify_ccw_orientation(
    const std::vector<std::string>& first_numerator,
    const std::vector<std::string>& first_denominator,
    const std::vector<std::string>& second_numerator,
    const std::vector<std::string>& second_denominator,
    const std::string& motion_parameter);

EventPartitionCertificate2 partition_cap_crossings(
    const std::vector<std::string>& first_numerator,
    const std::vector<std::string>& first_denominator,
    const std::vector<std::string>& second_numerator,
    const std::vector<std::string>& second_denominator,
    const std::string& cap_numerator,
    const std::string& cap_denominator,
    const PartitionEvent2& event);
