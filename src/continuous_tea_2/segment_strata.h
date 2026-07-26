#pragma once

#include "boundary_events.h"
#include "../exact_algebraic_1.h"
#include "segment_partition.h"

#include <vector>

struct SegmentBranchState2 {
    SegmentBoundaryBranch2 branch;
    ExactAlgebraicKernel1::Algebraic_real_1 rim_parameter;
    ExactAlgebraicKernel1::Polynomial_1 rim_factor;
    std::size_t rim_order;
};

std::vector<SegmentBranchState2> segment_branches_at(
    const std::vector<BoundaryFeatureRecord2>& records,
    const SegmentEventSource2& source,
    const std::string& witness_numerator,
    const std::string& witness_denominator);

std::vector<BranchPairDisposition2>
segment_branch_pair_dispositions(
    const std::vector<SegmentBranchState2>& branches,
    const SegmentEventSource2& source);

SegmentCellStratum2 construct_segment_cell_stratum(
    const Stock2& stock,
    const SegmentEventSource2& source,
    const std::string& witness_numerator,
    const std::string& witness_denominator);

std::vector<std::string>
segment_pair_literal_signs();

std::vector<std::string>
segment_rational_square_root_cases();
