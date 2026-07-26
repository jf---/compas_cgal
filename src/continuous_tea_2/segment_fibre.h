#pragma once

#include "boundary_events.h"
#include "partition_certificate.h"
#include "segment_source.h"
#include "segment_strata.h"

#include <vector>

struct SegmentFibreEvaluation2 {
    std::vector<SegmentBranchState2> branches;
    std::vector<std::string> active_branch_ids;
    std::vector<PartitionEvent2> events;
    std::vector<BranchPairDisposition2> pair_dispositions;
    bool algebraic_root_evaluated;
    bool original_equations_rechecked;
    bool orientation_rechecked;
    bool trim_predicates_rechecked;
};

SegmentFibreEvaluation2 evaluate_segment_fibre(
    const std::vector<BoundaryFeatureRecord2>& records,
    const SegmentEventSource2& source,
    const std::vector<ProjectionRecord2>& pullbacks,
    const std::vector<ProjectionInput2>& event_projections,
    const AlgebraicRootRecord2& root,
    const std::vector<PartitionEvent2>& events,
    const std::vector<SegmentBranchState2>& left_states,
    const std::vector<SegmentBranchState2>& right_states);
