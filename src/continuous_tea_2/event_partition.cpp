#include "event_partition.h"

#include "../exact_algebraic_1.h"
#include "parameter_charts.h"

EventPartitionCertificate2 partition_projections(
    const std::vector<ProjectionInput2>& projections)
{
    EventPartitionCertificate2 certificate =
        partition_integer_projections(projections);
    certificate.charts = parameter_charts();
    return certificate;
}
