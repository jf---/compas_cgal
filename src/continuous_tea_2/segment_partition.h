#pragma once

#include "boundary_events.h"
#include "partition_certificate.h"
#include "segment_source.h"

#include <string>
#include <vector>

struct SegmentBoundaryBranch2 {
    std::string branch_id;
    std::string feature_id;
    std::string support_id;
    std::string support_kind;
    std::string trim_id;
    std::string vertex_id;
    std::string material_side;
    std::string trim_disposition;
    std::string rim_chart_id;
    std::size_t rim_sheet_ordinal;
    std::string rim_root_id;
    std::vector<std::string> rim_factor_coefficients;
    std::size_t rim_root_ordinal;
};

struct BranchPairDisposition2 {
    std::string pair_sheet_id;
    std::string first_branch_id;
    std::string second_branch_id;
    std::string orientation_disposition;
    std::string cap_disposition;
};

struct SegmentEventStratum2 {
    std::string kind;
    std::string root_id;
    std::string local_root_id;
    std::string global_fibre_id;
    std::string chart_id;
    std::string witness_numerator;
    std::string witness_denominator;
    std::vector<std::string> root_factor_coefficients;
    std::size_t root_ordinal;
    std::vector<std::string> active_branch_ids;
    std::vector<PartitionEvent2> events;
    std::vector<BranchPairDisposition2>
        left_pair_dispositions;
    std::vector<BranchPairDisposition2> pair_dispositions;
    std::vector<BranchPairDisposition2>
        right_pair_dispositions;
    bool algebraic_root_evaluated;
    bool original_equations_rechecked;
    bool orientation_rechecked;
    bool trim_predicates_rechecked;
};

struct SegmentCellStratum2 {
    std::vector<SegmentBoundaryBranch2> branches;
    SegmentEventStratum2 stratum;
};

struct SegmentEventPartition2 {
    SegmentEventSource2 source;
    std::vector<std::string> boundary_feature_ids;
    std::vector<ProjectionRecord2> projections;
    std::vector<SegmentBoundaryBranch2> branches;
    std::vector<SegmentEventStratum2> strata;
    EventPartitionCertificate2 certificate;
    std::string canonical_bytes;
    std::string canonical_digest;
};

struct VerifiedSegmentEventPartition2 {
    ContinuousTeaVerdict verdict;
    SegmentEventPartition2 partition;
};

SegmentEventPartition2 construct_segment_event_partition(
    const Stock2& stock,
    const SegmentEventSource2& source);

SegmentEventPartition2 construct_segment_event_partition(
    const Stock2& stock,
    double x0,
    double y0,
    double x1,
    double y1,
    double tool_radius,
    double cap_chord_ratio);

VerifiedSegmentEventPartition2 verify_segment_event_partition(
    const Stock2& stock,
    const SegmentEventSource2& source,
    const SegmentEventPartition2& candidate);

std::string canonical_segment_cell_stratum(
    const SegmentCellStratum2& cell);

SegmentEventPartition2 mutate_segment_event_partition(
    const SegmentEventPartition2& partition,
    const std::string& mutation);

class IncompleteSegmentPartitionError
    : public EventSubstrateError {
public:
    using EventSubstrateError::EventSubstrateError;
};
