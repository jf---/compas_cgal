#pragma once

#include "continuous_tea_2/partition_certificate.h"

#include <cstddef>
#include <string>
#include <vector>

#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Algebraic_kernel_d_2.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>
#include <CGAL/CORE/BigInt.h>

using ExactAlgebraicInteger1 = CORE::BigInt;
using ExactAlgebraicKernel1 =
    CGAL::Algebraic_kernel_d_1<ExactAlgebraicInteger1>;
using ExactAlgebraicKernel2 =
    CGAL::Algebraic_kernel_d_2<ExactAlgebraicInteger1>;
using ExactAlgebraicTraits2 =
    CGAL::Arr_algebraic_segment_traits_2<ExactAlgebraicInteger1>;

AlgebraicBackendEvidence2 exact_algebraic_backend_evidence();

std::string algebraic_root_id_v1(
    const std::vector<ExactAlgebraicInteger1>& primitive_factor,
    std::size_t root_ordinal);

EventPartitionCertificate2 partition_integer_projections(
    const std::vector<ProjectionInput2>& projections);
