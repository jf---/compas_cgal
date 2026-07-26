#pragma once

#include "result.h"

#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct AlgebraicBackendEvidence2 {
    std::string cgal_version;
    std::string integer_backend;
    std::string algebraic_kernel_1;
    std::string algebraic_kernel_2;
    std::string arrangement_traits;
    std::vector<std::string> compile_definitions;
};

struct ParameterChart2 {
    std::string chart_id;
    std::string family;
    std::string domain_low;
    std::string domain_high;
    std::string orientation;
    std::string start_seam_id;
    std::string end_seam_id;
    bool owns_start_seam;
    bool owns_end_seam;
};

struct PartitionEvent2 {
    std::string kind;
    std::string feature_id;
    std::string support_id;
    std::string trim_id;
    std::string vertex_id;
    std::string branch_id;
    std::string disposition;
    std::size_t left_active_count = 0;
    std::size_t right_active_count = 0;
    bool incidence_permutation_rechecked = false;
    bool original_equations_rechecked = false;
    bool orientation_rechecked = false;
    std::string trim_disposition = "unchecked";
    std::string pair_sheet_id;
    std::string first_feature_id;
    std::string second_feature_id;
    std::string first_chart_id;
    std::string second_chart_id;
    std::string first_branch_id;
    std::string second_branch_id;

    PartitionEvent2() = default;
    PartitionEvent2(
        std::string kind_,
        std::string feature_id_,
        std::string support_id_,
        std::string trim_id_,
        std::string vertex_id_,
        std::string branch_id_,
        std::string disposition_)
        : kind(std::move(kind_)),
          feature_id(std::move(feature_id_)),
          support_id(std::move(support_id_)),
          trim_id(std::move(trim_id_)),
          vertex_id(std::move(vertex_id_)),
          branch_id(std::move(branch_id_)),
          disposition(std::move(disposition_))
    {
    }
};

struct ProjectionInput2 {
    std::string projection_id;
    std::vector<std::string> coefficients;
    std::vector<PartitionEvent2> events;

    ProjectionInput2() = default;
    ProjectionInput2(
        std::string projection_id_,
        std::vector<std::string> coefficients_,
        std::vector<PartitionEvent2> events_)
        : projection_id(std::move(projection_id_)),
          coefficients(std::move(coefficients_)),
          events(std::move(events_))
    {
    }
};

struct AlgebraicRootRecord2 {
    std::string root_id;
    std::vector<std::string> factor_coefficients;
    std::size_t root_ordinal;
    unsigned int multiplicity;
    std::string interval_low;
    std::string interval_high;
};

struct ParameterCell2 {
    std::string lower_root_id;
    std::string upper_root_id;
    std::string witness_numerator;
    std::string witness_denominator;
    std::string disposition;
};

struct ActiveBoundaryBranch2 {
    std::string branch_id;
    std::string support_id;
};

struct EventFibre2 {
    std::string root_id;
    std::vector<PartitionEvent2> events;
    std::vector<ActiveBoundaryBranch2>
        left_active_branches;
    std::vector<ActiveBoundaryBranch2>
        right_active_branches;
    std::string ccw_direction;
    std::string cw_direction;
};

struct RationalTrimInterval2 {
    std::string rim_parameter;
    std::string motion_domain_low;
    std::string motion_domain_high;
};

struct TrimmedLineBranch2 {
    std::string rim_parameter;
    std::string motion_domain_low;
    std::string motion_domain_high;
    std::optional<RationalTrimInterval2>
        rational_convenience;
    AlgebraicRootRecord2 rim_root;
    std::vector<std::vector<std::string>>
        lower_trim_predicate_rows;
    std::vector<std::vector<std::string>>
        upper_trim_predicate_rows;
    bool rational_convenience_available;
    bool domain_nonempty_rechecked;
    bool complementarity_rechecked;
    std::string trim_disposition;
    bool rejected_outside_closed_domain;
    std::string feature_id;
    std::string local_support_id;
    std::string local_trimmed_feature_id;
    std::string trim_id;
    std::string branch_id;
};

struct ProjectedRegularizationVertex2 {
    AlgebraicRootRecord2 root;
    std::string vertex_id;
    std::string first_trim_disposition;
    std::string second_trim_disposition;
    std::string conjugate_disposition;
};

struct ProjectionRecord2 {
    std::string projection_id;
    std::vector<std::vector<std::string>> coefficient_rows;
    std::vector<std::vector<std::string>> factor_coefficients;
    int actual_motion_degree = -1;
    int actual_rim_degree = -1;
    int bound_motion_degree = -1;
    int bound_rim_degree = -1;
    std::string degree_bound_id;
    std::string normalized_coefficient_bytes;
};

struct OverlapInterval2 {
    std::string kind;
    std::string domain_low;
    std::string domain_high;
    std::string witness_numerator;
    std::string witness_denominator;
    std::string orientation_disposition;
    std::string feature_id;
    std::string support_id;
    std::string trim_id;
    std::string branch_id;
};

struct ChartSeam2 {
    std::string seam_id;
    std::string owner_chart_id;
};

struct EventPartitionCertificate2 {
    AlgebraicBackendEvidence2 build_evidence;
    std::vector<ParameterChart2> charts;
    std::vector<ProjectionRecord2> projections;
    std::vector<AlgebraicRootRecord2> roots;
    std::vector<ParameterCell2> cells;
    std::vector<EventFibre2> fibres;
    std::vector<OverlapInterval2> overlaps;
    std::vector<ChartSeam2> seams;
    std::string source_kind;
    std::string source_payload;
    std::string canonical_bytes;
    std::string canonical_digest;
};

struct VerifiedEventPartition2 {
    ContinuousTeaVerdict verdict;
    EventPartitionCertificate2 partition;
};

class EventSubstrateError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class AlgebraicBackendError : public EventSubstrateError {
public:
    using EventSubstrateError::EventSubstrateError;
};

class InvalidAlgebraicPolynomialError : public EventSubstrateError {
public:
    using EventSubstrateError::EventSubstrateError;
};

class UnsupportedAlgebraicDegreeError : public EventSubstrateError {
public:
    using EventSubstrateError::EventSubstrateError;
};

class AlgebraicRootIsolationError : public EventSubstrateError {
public:
    using EventSubstrateError::EventSubstrateError;
};

class ChartCoverageError : public EventSubstrateError {
public:
    using EventSubstrateError::EventSubstrateError;
};

class ProjectionDegreeBoundError : public EventSubstrateError {
public:
    using EventSubstrateError::EventSubstrateError;
};

class TrimFilterError : public EventSubstrateError {
public:
    using EventSubstrateError::EventSubstrateError;
};

class TrimEndpointOffSupportError : public TrimFilterError {
public:
    using TrimFilterError::TrimFilterError;
};

class DegenerateTrimError : public TrimFilterError {
public:
    using TrimFilterError::TrimFilterError;
};

class UnsupportedLineMotionError : public TrimFilterError {
public:
    using TrimFilterError::TrimFilterError;
};

class EventPartitionVerificationError : public EventSubstrateError {
public:
    using EventSubstrateError::EventSubstrateError;
};
