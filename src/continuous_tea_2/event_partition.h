#pragma once

#include "partition_certificate.h"

#include <vector>

EventPartitionCertificate2 partition_projections(
    const std::vector<ProjectionInput2>& projections);
