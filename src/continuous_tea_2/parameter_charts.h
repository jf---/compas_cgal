#pragma once

#include "partition_certificate.h"

#include <string>
#include <vector>

std::vector<ParameterChart2> parameter_charts();

VerifiedEventPartition2 verify_chart_coverage(
    const std::vector<ParameterChart2>& charts);
