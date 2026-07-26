#include "compas.h"
#include "continuous_tea_2/boundary_events.h"
#include "continuous_tea_2/cap_partition.h"
#include "continuous_tea_2/event_certificate.h"
#include "continuous_tea_2/event_partition.h"
#include "continuous_tea_2/parameter_charts.h"
#include "continuous_tea_2/partition_certificate.h"
#include "continuous_tea_2/result.h"
#include "continuous_tea_2/segment_partition.h"
#include "continuous_tea_2/segment_oracle.h"
#include "continuous_tea_2/segment_source.h"
#include "continuous_tea_2/segment_strata.h"
#include "continuous_tea_2/sha256.h"
#include "exact_algebraic_1.h"

#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <new>
#include <utility>

namespace {

nb::bytes as_bytes(const std::string& value)
{
    return nb::bytes(value.data(), value.size());
}

std::string from_bytes(const nb::bytes& value)
{
    char* data = nullptr;
    Py_ssize_t size = 0;
    if (PyBytes_AsStringAndSize(
            value.ptr(),
            &data,
            &size)
        != 0) {
        throw nb::python_error();
    }
    return std::string(
        data,
        static_cast<std::size_t>(size));
}

nb::list coefficient_bytes(
    const std::vector<std::string>& values)
{
    nb::list result;
    for (const std::string& value : values) {
        result.append(as_bytes(value));
    }
    return result;
}

nb::object bytes_tuple(
    const std::vector<std::string>& values)
{
    return nb::module_::import_("builtins")
        .attr("tuple")(coefficient_bytes(values));
}

nb::object string_tuple(
    const std::vector<std::string>& values)
{
    nb::list result;
    for (const std::string& value : values) {
        result.append(value);
    }
    return nb::module_::import_("builtins")
        .attr("tuple")(result);
}

nb::object string_matrix_tuple(
    const std::vector<std::vector<std::string>>& rows)
{
    nb::list result;
    for (const auto& row : rows) {
        result.append(string_tuple(row));
    }
    return nb::module_::import_("builtins")
        .attr("tuple")(result);
}

} // namespace

NB_MODULE(_continuous_tea_2, m)
{
    nb::exception<BoundaryExtractionError> extraction_error(
        m,
        "BoundaryExtractionError");
    nb::exception<DegenerateBoundarySupportError>(
        m,
        "DegenerateBoundarySupportError",
        extraction_error.ptr());
    nb::exception<MissingBoundaryEndpointError>(
        m,
        "MissingBoundaryEndpointError",
        extraction_error.ptr());
    nb::exception<MissingBoundaryIntersectionError>(
        m,
        "MissingBoundaryIntersectionError",
        extraction_error.ptr());
    nb::exception<BoundaryFeatureIndexError>(
        m,
        "BoundaryFeatureIndexError",
        extraction_error.ptr());

    nb::exception<EventSubstrateError> substrate_error(
        m,
        "EventSubstrateError");
    nb::exception<AlgebraicBackendError>(
        m,
        "AlgebraicBackendError",
        substrate_error.ptr());
    nb::exception<InvalidAlgebraicPolynomialError>(
        m,
        "InvalidAlgebraicPolynomialError",
        substrate_error.ptr());
    nb::exception<UnsupportedAlgebraicDegreeError>(
        m,
        "UnsupportedAlgebraicDegreeError",
        substrate_error.ptr());
    nb::exception<AlgebraicRootIsolationError>(
        m,
        "AlgebraicRootIsolationError",
        substrate_error.ptr());
    nb::exception<ChartCoverageError>(
        m,
        "ChartCoverageError",
        substrate_error.ptr());
    nb::exception<ProjectionDegreeBoundError>(
        m,
        "ProjectionDegreeBoundError",
        substrate_error.ptr());
    nb::exception<TrimFilterError> trim_filter_error(
        m,
        "TrimFilterError",
        substrate_error.ptr());
    nb::exception<TrimEndpointOffSupportError>(
        m,
        "TrimEndpointOffSupportError",
        trim_filter_error.ptr());
    nb::exception<DegenerateTrimError>(
        m,
        "DegenerateTrimError",
        trim_filter_error.ptr());
    nb::exception<UnsupportedLineMotionError>(
        m,
        "UnsupportedLineMotionError",
        trim_filter_error.ptr());
    nb::exception<EventPartitionVerificationError>(
        m,
        "EventPartitionVerificationError",
        substrate_error.ptr());
    nb::exception<NonFiniteSegmentInputError>(
        m,
        "NonFiniteSegmentInputError",
        substrate_error.ptr());
    nb::exception<ZeroLengthSegmentMotionError>(
        m,
        "ZeroLengthSegmentMotionError",
        substrate_error.ptr());
    nb::exception<NonPositiveToolRadiusError>(
        m,
        "NonPositiveToolRadiusError",
        substrate_error.ptr());
    nb::exception<InvalidCapChordRatioError>(
        m,
        "InvalidCapChordRatioError",
        substrate_error.ptr());
    nb::exception<UnsupportedAlgebraicVertexProjectionError>(
        m,
        "UnsupportedAlgebraicVertexProjectionError",
        substrate_error.ptr());
    nb::exception<IncompleteSegmentPartitionError>(
        m,
        "IncompleteSegmentPartitionError",
        substrate_error.ptr());

    nb::enum_<ContinuousTeaVerdict>(
        m,
        "ContinuousTeaVerdict")
        .value(
            "CERTIFIED",
            ContinuousTeaVerdict::CERTIFIED)
        .value(
            "CAP_EXCEEDED",
            ContinuousTeaVerdict::CAP_EXCEEDED)
        .value(
            "UNRESOLVED_DEGENERACY",
            ContinuousTeaVerdict::UNRESOLVED_DEGENERACY);

    nb::class_<AlgebraicBackendEvidence2>(
        m,
        "AlgebraicBackendEvidence2")
        .def_ro(
            "cgal_version",
            &AlgebraicBackendEvidence2::cgal_version)
        .def_ro(
            "integer_backend",
            &AlgebraicBackendEvidence2::integer_backend)
        .def_ro(
            "algebraic_kernel_1",
            &AlgebraicBackendEvidence2::algebraic_kernel_1)
        .def_ro(
            "algebraic_kernel_2",
            &AlgebraicBackendEvidence2::algebraic_kernel_2)
        .def_ro(
            "arrangement_traits",
            &AlgebraicBackendEvidence2::arrangement_traits)
        .def_prop_ro(
            "compile_definitions",
            [](const AlgebraicBackendEvidence2& evidence) {
                return string_tuple(
                    evidence.compile_definitions);
            });

    nb::class_<ExactBinary64Rational2>(
        m,
        "ExactBinary64Rational2")
        .def_prop_ro(
            "numerator",
            &ExactBinary64Rational2::numerator)
        .def_prop_ro(
            "denominator",
            &ExactBinary64Rational2::denominator)
        .def_prop_ro(
            "text",
            &ExactBinary64Rational2::text)
        .def_prop_ro(
            "canonical_bytes",
            [](const ExactBinary64Rational2& value) {
                return as_bytes(
                    value.canonical_bytes());
            });

    nb::class_<SegmentEventSource2>(
        m,
        "SegmentEventSource2")
        .def_static(
            "from_binary64",
            &SegmentEventSource2::from_binary64,
            "x0"_a,
            "y0"_a,
            "x1"_a,
            "y1"_a,
            "tool_radius"_a,
            "cap_chord_ratio"_a)
        .def_prop_ro(
            "x0",
            &SegmentEventSource2::x0,
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "y0",
            &SegmentEventSource2::y0,
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "x1",
            &SegmentEventSource2::x1,
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "y1",
            &SegmentEventSource2::y1,
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "tool_radius",
            &SegmentEventSource2::tool_radius,
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "cap_chord_ratio",
            &SegmentEventSource2::cap_chord_ratio,
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "motion_data",
            [](const SegmentEventSource2& source) {
                return string_tuple(
                    source.motion_data());
            })
        .def_prop_ro(
            "canonical_bytes",
            [](const SegmentEventSource2& source) {
                return as_bytes(
                    source.canonical_bytes());
            })
        .def_prop_ro(
            "canonical_digest",
            [](const SegmentEventSource2& source) {
                return as_bytes(
                    source.canonical_digest());
            });

    nb::class_<SegmentBoundaryBranch2>(
        m,
        "SegmentBoundaryBranch2")
        .def_prop_ro(
            "branch_id",
            [](const SegmentBoundaryBranch2& branch) {
                return as_bytes(branch.branch_id);
            })
        .def_prop_ro(
            "feature_id",
            [](const SegmentBoundaryBranch2& branch) {
                return as_bytes(branch.feature_id);
            })
        .def_prop_ro(
            "support_id",
            [](const SegmentBoundaryBranch2& branch) {
                return as_bytes(branch.support_id);
            })
        .def_ro(
            "support_kind",
            &SegmentBoundaryBranch2::support_kind)
        .def_prop_ro(
            "trim_id",
            [](const SegmentBoundaryBranch2& branch) {
                return as_bytes(branch.trim_id);
            })
        .def_prop_ro(
            "vertex_id",
            [](const SegmentBoundaryBranch2& branch) {
                return as_bytes(branch.vertex_id);
            })
        .def_ro(
            "material_side",
            &SegmentBoundaryBranch2::material_side)
        .def_ro(
            "trim_disposition",
            &SegmentBoundaryBranch2::trim_disposition)
        .def_ro(
            "rim_chart_id",
            &SegmentBoundaryBranch2::rim_chart_id)
        .def_ro(
            "rim_sheet_ordinal",
            &SegmentBoundaryBranch2::
                rim_sheet_ordinal)
        .def_prop_ro(
            "rim_root_id",
            [](const SegmentBoundaryBranch2& branch) {
                return as_bytes(branch.rim_root_id);
            })
        .def_prop_ro(
            "rim_factor_coefficients",
            [](const SegmentBoundaryBranch2& branch) {
                return string_tuple(
                    branch.rim_factor_coefficients);
            })
        .def_ro(
            "rim_root_ordinal",
            &SegmentBoundaryBranch2::
                rim_root_ordinal);

    nb::class_<BranchPairDisposition2>(
        m,
        "BranchPairDisposition2")
        .def_prop_ro(
            "pair_sheet_id",
            [](const BranchPairDisposition2& pair) {
                return as_bytes(
                    pair.pair_sheet_id);
            })
        .def_prop_ro(
            "first_branch_id",
            [](const BranchPairDisposition2& pair) {
                return as_bytes(
                    pair.first_branch_id);
            })
        .def_prop_ro(
            "second_branch_id",
            [](const BranchPairDisposition2& pair) {
                return as_bytes(
                    pair.second_branch_id);
            })
        .def_ro(
            "orientation_disposition",
            &BranchPairDisposition2::
                orientation_disposition)
        .def_ro(
            "cap_disposition",
            &BranchPairDisposition2::cap_disposition);

    nb::class_<SegmentEventStratum2>(
        m,
        "SegmentEventStratum2")
        .def_ro("kind", &SegmentEventStratum2::kind)
        .def_prop_ro(
            "root_id",
            [](const SegmentEventStratum2& stratum) {
                return as_bytes(stratum.root_id);
            })
        .def_prop_ro(
            "local_root_id",
            [](const SegmentEventStratum2& stratum) {
                return as_bytes(
                    stratum.local_root_id);
            })
        .def_prop_ro(
            "global_fibre_id",
            [](const SegmentEventStratum2& stratum) {
                return as_bytes(
                    stratum.global_fibre_id);
            })
        .def_ro(
            "chart_id",
            &SegmentEventStratum2::chart_id)
        .def_ro(
            "witness_numerator",
            &SegmentEventStratum2::
                witness_numerator)
        .def_ro(
            "witness_denominator",
            &SegmentEventStratum2::
                witness_denominator)
        .def_prop_ro(
            "root_factor_coefficients",
            [](const SegmentEventStratum2& stratum) {
                return string_tuple(
                    stratum.root_factor_coefficients);
            })
        .def_ro(
            "root_ordinal",
            &SegmentEventStratum2::root_ordinal)
        .def_prop_ro(
            "active_branch_ids",
            [](const SegmentEventStratum2& stratum) {
                return bytes_tuple(
                    stratum.active_branch_ids);
            })
        .def_ro(
            "events",
            &SegmentEventStratum2::events)
        .def_ro(
            "left_pair_dispositions",
            &SegmentEventStratum2::
                left_pair_dispositions)
        .def_ro(
            "pair_dispositions",
            &SegmentEventStratum2::
                pair_dispositions)
        .def_ro(
            "right_pair_dispositions",
            &SegmentEventStratum2::
                right_pair_dispositions)
        .def_ro(
            "algebraic_root_evaluated",
            &SegmentEventStratum2::
                algebraic_root_evaluated)
        .def_ro(
            "original_equations_rechecked",
            &SegmentEventStratum2::
                original_equations_rechecked)
        .def_ro(
            "orientation_rechecked",
            &SegmentEventStratum2::
                orientation_rechecked)
        .def_ro(
            "trim_predicates_rechecked",
            &SegmentEventStratum2::
                trim_predicates_rechecked);

    nb::class_<SegmentCellStratum2>(
        m,
        "SegmentCellStratum2")
        .def_ro(
            "branches",
            &SegmentCellStratum2::branches)
        .def_ro(
            "stratum",
            &SegmentCellStratum2::stratum);

    nb::class_<SegmentEventPartition2>(
        m,
        "SegmentEventPartition2")
        .def_prop_ro(
            "source",
            [](const SegmentEventPartition2& partition)
                -> const SegmentEventSource2& {
                return partition.source;
            },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "boundary_feature_ids",
            [](const SegmentEventPartition2& partition) {
                return bytes_tuple(
                    partition.boundary_feature_ids);
            })
        .def_ro(
            "projections",
            &SegmentEventPartition2::projections)
        .def_ro(
            "branches",
            &SegmentEventPartition2::branches)
        .def_ro(
            "strata",
            &SegmentEventPartition2::strata)
        .def_ro(
            "certificate",
            &SegmentEventPartition2::certificate)
        .def_prop_ro(
            "canonical_bytes",
            [](const SegmentEventPartition2& partition) {
                return as_bytes(
                    partition.canonical_bytes);
            })
        .def_prop_ro(
            "canonical_digest",
            [](const SegmentEventPartition2& partition) {
                return as_bytes(
                    partition.canonical_digest);
            });

    nb::class_<VerifiedSegmentEventPartition2>(
        m,
        "VerifiedSegmentEventPartition2")
        .def_ro(
            "verdict",
            &VerifiedSegmentEventPartition2::verdict)
        .def_ro(
            "partition",
            &VerifiedSegmentEventPartition2::partition);

    nb::class_<ParameterChart2>(m, "ParameterChart2")
        .def_ro("chart_id", &ParameterChart2::chart_id)
        .def_ro("family", &ParameterChart2::family)
        .def_ro("domain_low", &ParameterChart2::domain_low)
        .def_ro(
            "domain_high",
            &ParameterChart2::domain_high)
        .def_ro(
            "orientation",
            &ParameterChart2::orientation)
        .def_prop_ro(
            "start_seam_id",
            [](const ParameterChart2& chart) {
                return as_bytes(chart.start_seam_id);
            })
        .def_prop_ro(
            "end_seam_id",
            [](const ParameterChart2& chart) {
                return as_bytes(chart.end_seam_id);
            })
        .def_ro(
            "owns_start_seam",
            &ParameterChart2::owns_start_seam)
        .def_ro(
            "owns_end_seam",
            &ParameterChart2::owns_end_seam);

    nb::class_<PartitionEvent2>(m, "PartitionEvent2")
        .def(
            "__init__",
            [](PartitionEvent2* event,
               std::string kind,
               nb::bytes feature_id,
               nb::bytes support_id,
               nb::bytes trim_id,
               nb::bytes vertex_id,
               nb::bytes branch_id,
               std::string disposition) {
                new (event) PartitionEvent2(
                    std::move(kind),
                    from_bytes(feature_id),
                    from_bytes(support_id),
                    from_bytes(trim_id),
                    from_bytes(vertex_id),
                    from_bytes(branch_id),
                    std::move(disposition));
            },
            "kind"_a,
            "feature_id"_a,
            "support_id"_a,
            "trim_id"_a,
            "vertex_id"_a,
            "branch_id"_a,
            "disposition"_a)
        .def_ro("kind", &PartitionEvent2::kind)
        .def_prop_ro(
            "feature_id",
            [](const PartitionEvent2& event) {
                return as_bytes(event.feature_id);
            })
        .def_prop_ro(
            "support_id",
            [](const PartitionEvent2& event) {
                return as_bytes(event.support_id);
            })
        .def_prop_ro(
            "trim_id",
            [](const PartitionEvent2& event) {
                return as_bytes(event.trim_id);
            })
        .def_prop_ro(
            "vertex_id",
            [](const PartitionEvent2& event) {
                return as_bytes(event.vertex_id);
            })
        .def_prop_ro(
            "branch_id",
            [](const PartitionEvent2& event) {
                return as_bytes(event.branch_id);
            })
        .def_ro(
            "disposition",
            &PartitionEvent2::disposition)
        .def_ro(
            "left_active_count",
            &PartitionEvent2::left_active_count)
        .def_ro(
            "right_active_count",
            &PartitionEvent2::right_active_count)
        .def_ro(
            "incidence_permutation_rechecked",
            &PartitionEvent2::
                incidence_permutation_rechecked)
        .def_ro(
            "original_equations_rechecked",
            &PartitionEvent2::original_equations_rechecked)
        .def_ro(
            "orientation_rechecked",
            &PartitionEvent2::orientation_rechecked)
        .def_ro(
            "trim_disposition",
            &PartitionEvent2::trim_disposition)
        .def_prop_ro(
            "pair_sheet_id",
            [](const PartitionEvent2& event) {
                return as_bytes(event.pair_sheet_id);
            })
        .def_prop_ro(
            "first_feature_id",
            [](const PartitionEvent2& event) {
                return as_bytes(event.first_feature_id);
            })
        .def_prop_ro(
            "second_feature_id",
            [](const PartitionEvent2& event) {
                return as_bytes(event.second_feature_id);
            })
        .def_ro(
            "first_chart_id",
            &PartitionEvent2::first_chart_id)
        .def_ro(
            "second_chart_id",
            &PartitionEvent2::second_chart_id)
        .def_prop_ro(
            "first_branch_id",
            [](const PartitionEvent2& event) {
                return as_bytes(event.first_branch_id);
            })
        .def_prop_ro(
            "second_branch_id",
            [](const PartitionEvent2& event) {
                return as_bytes(event.second_branch_id);
            });

    nb::class_<ProjectionInput2>(m, "ProjectionInput2")
        .def(
            nb::init<
                std::string,
                std::vector<std::string>,
                std::vector<PartitionEvent2>>(),
            "projection_id"_a,
            "coefficients"_a,
            "events"_a);

    nb::class_<AlgebraicRootRecord2>(
        m,
        "AlgebraicRootRecord2")
        .def(
            "__init__",
            [](AlgebraicRootRecord2* self,
               const nb::bytes& root_id,
               std::vector<std::string> factor_coefficients,
               std::size_t root_ordinal,
               unsigned int multiplicity,
               std::string interval_low,
               std::string interval_high) {
                new (self) AlgebraicRootRecord2{
                    from_bytes(root_id),
                    std::move(factor_coefficients),
                    root_ordinal,
                    multiplicity,
                    std::move(interval_low),
                    std::move(interval_high),
                };
            },
            "root_id"_a,
            "factor_coefficients"_a,
            "root_ordinal"_a,
            "multiplicity"_a,
            "interval_low"_a,
            "interval_high"_a)
        .def_prop_ro(
            "root_id",
            [](const AlgebraicRootRecord2& root) {
                return as_bytes(root.root_id);
            })
        .def_prop_ro(
            "factor_coefficients",
            [](const AlgebraicRootRecord2& root) {
                return string_tuple(
                    root.factor_coefficients);
            })
        .def_ro(
            "root_ordinal",
            &AlgebraicRootRecord2::root_ordinal)
        .def_ro(
            "multiplicity",
            &AlgebraicRootRecord2::multiplicity)
        .def_ro(
            "interval_low",
            &AlgebraicRootRecord2::interval_low)
        .def_ro(
            "interval_high",
            &AlgebraicRootRecord2::interval_high);

    nb::class_<ParameterCell2>(m, "ParameterCell2")
        .def_prop_ro(
            "lower_root_id",
            [](const ParameterCell2& cell) {
                return as_bytes(cell.lower_root_id);
            })
        .def_prop_ro(
            "upper_root_id",
            [](const ParameterCell2& cell) {
                return as_bytes(cell.upper_root_id);
            })
        .def_ro(
            "witness_numerator",
            &ParameterCell2::witness_numerator)
        .def_ro(
            "witness_denominator",
            &ParameterCell2::witness_denominator)
        .def_ro(
            "disposition",
            &ParameterCell2::disposition);

    nb::class_<EventFibre2>(m, "EventFibre2")
        .def_prop_ro(
            "root_id",
            [](const EventFibre2& fibre) {
                return as_bytes(fibre.root_id);
            })
        .def_ro("events", &EventFibre2::events);

    nb::class_<RationalTrimInterval2>(
        m,
        "RationalTrimInterval2")
        .def_ro(
            "rim_parameter",
            &RationalTrimInterval2::rim_parameter)
        .def_ro(
            "motion_domain_low",
            &RationalTrimInterval2::motion_domain_low)
        .def_ro(
            "motion_domain_high",
            &RationalTrimInterval2::motion_domain_high);

    nb::class_<TrimmedLineBranch2>(m, "TrimmedLineBranch2")
        .def_ro(
            "rim_parameter",
            &TrimmedLineBranch2::rim_parameter)
        .def_ro(
            "motion_domain_low",
            &TrimmedLineBranch2::motion_domain_low)
        .def_ro(
            "motion_domain_high",
            &TrimmedLineBranch2::motion_domain_high)
        .def_ro(
            "rational_convenience",
            &TrimmedLineBranch2::rational_convenience)
        .def_ro(
            "rim_root",
            &TrimmedLineBranch2::rim_root)
        .def_prop_ro(
            "lower_trim_predicate_rows",
            [](const TrimmedLineBranch2& branch) {
                return string_matrix_tuple(
                    branch.lower_trim_predicate_rows);
            })
        .def_prop_ro(
            "upper_trim_predicate_rows",
            [](const TrimmedLineBranch2& branch) {
                return string_matrix_tuple(
                    branch.upper_trim_predicate_rows);
            })
        .def_ro(
            "rational_convenience_available",
            &TrimmedLineBranch2::
                rational_convenience_available)
        .def_ro(
            "domain_nonempty_rechecked",
            &TrimmedLineBranch2::
                domain_nonempty_rechecked)
        .def_ro(
            "complementarity_rechecked",
            &TrimmedLineBranch2::
                complementarity_rechecked)
        .def_ro(
            "trim_disposition",
            &TrimmedLineBranch2::trim_disposition)
        .def_ro(
            "rejected_outside_closed_domain",
            &TrimmedLineBranch2::
                rejected_outside_closed_domain)
        .def_prop_ro(
            "feature_id",
            [](const TrimmedLineBranch2& branch) {
                return as_bytes(branch.feature_id);
            })
        .def_prop_ro(
            "local_support_id",
            [](const TrimmedLineBranch2& branch) {
                return as_bytes(
                    branch.local_support_id);
            })
        .def_prop_ro(
            "local_trimmed_feature_id",
            [](const TrimmedLineBranch2& branch) {
                return as_bytes(
                    branch.local_trimmed_feature_id);
            })
        .def_prop_ro(
            "trim_id",
            [](const TrimmedLineBranch2& branch) {
                return as_bytes(branch.trim_id);
            })
        .def_prop_ro(
            "branch_id",
            [](const TrimmedLineBranch2& branch) {
                return as_bytes(branch.branch_id);
            });

    nb::class_<ProjectedRegularizationVertex2>(
        m,
        "ProjectedRegularizationVertex2")
        .def_ro(
            "root",
            &ProjectedRegularizationVertex2::root)
        .def_prop_ro(
            "vertex_id",
            [](const ProjectedRegularizationVertex2& projected) {
                return as_bytes(projected.vertex_id);
            })
        .def_ro(
            "first_trim_disposition",
            &ProjectedRegularizationVertex2::
                first_trim_disposition)
        .def_ro(
            "second_trim_disposition",
            &ProjectedRegularizationVertex2::
                second_trim_disposition)
        .def_ro(
            "conjugate_disposition",
            &ProjectedRegularizationVertex2::
                conjugate_disposition);

    nb::class_<CcwOrientation2>(m, "CcwOrientation2")
        .def_ro(
            "disposition",
            &CcwOrientation2::disposition)
        .def_ro(
            "determinant_sign",
            &CcwOrientation2::determinant_sign)
        .def_ro(
            "dot_sign",
            &CcwOrientation2::dot_sign);

    nb::class_<ProjectionRecord2>(m, "ProjectionRecord2")
        .def_ro(
            "projection_id",
            &ProjectionRecord2::projection_id)
        .def_prop_ro(
            "coefficient_rows",
            [](const ProjectionRecord2& projection) {
                return string_matrix_tuple(
                    projection.coefficient_rows);
            })
        .def_prop_ro(
            "factor_coefficients",
            [](const ProjectionRecord2& projection) {
                return string_matrix_tuple(
                    projection.factor_coefficients);
            })
        .def_prop_ro(
            "actual_degree",
            [](const ProjectionRecord2& projection) {
                return nb::make_tuple(
                    projection.actual_motion_degree,
                    projection.actual_rim_degree);
            })
        .def_prop_ro(
            "degree_bound",
            [](const ProjectionRecord2& projection) {
                return nb::make_tuple(
                    projection.bound_motion_degree,
                    projection.bound_rim_degree);
            })
        .def_ro(
            "degree_bound_id",
            &ProjectionRecord2::degree_bound_id)
        .def_prop_ro(
            "normalized_coefficient_bytes",
            [](const ProjectionRecord2& projection) {
                return as_bytes(
                    projection.normalized_coefficient_bytes);
            });

    nb::class_<OverlapInterval2>(m, "OverlapInterval2")
        .def_ro("kind", &OverlapInterval2::kind)
        .def_ro(
            "domain_low",
            &OverlapInterval2::domain_low)
        .def_ro(
            "domain_high",
            &OverlapInterval2::domain_high)
        .def_ro(
            "witness_numerator",
            &OverlapInterval2::witness_numerator)
        .def_ro(
            "witness_denominator",
            &OverlapInterval2::witness_denominator)
        .def_ro(
            "orientation_disposition",
            &OverlapInterval2::orientation_disposition);

    nb::class_<ChartSeam2>(m, "ChartSeam2")
        .def_prop_ro(
            "seam_id",
            [](const ChartSeam2& seam) {
                return as_bytes(seam.seam_id);
            })
        .def_ro(
            "owner_chart_id",
            &ChartSeam2::owner_chart_id);

    nb::class_<EventPartitionCertificate2>(
        m,
        "EventPartitionCertificate2")
        .def(
            "__init__",
            [](EventPartitionCertificate2* self,
               AlgebraicBackendEvidence2 build_evidence,
               std::vector<ParameterChart2> charts,
               std::vector<ProjectionRecord2> projections,
               std::vector<AlgebraicRootRecord2> roots,
               std::vector<ParameterCell2> cells,
               std::vector<EventFibre2> fibres,
               std::vector<OverlapInterval2> overlaps,
               std::vector<ChartSeam2> seams,
               std::string source_kind,
               const nb::bytes& source_payload) {
                EventPartitionCertificate2 certificate{
                    std::move(build_evidence),
                    std::move(charts),
                    std::move(projections),
                    std::move(roots),
                    std::move(cells),
                    std::move(fibres),
                    std::move(overlaps),
                    std::move(seams),
                    std::move(source_kind),
                    from_bytes(source_payload),
                    {},
                    {},
                };
                finalize_event_partition(certificate);
                new (self) EventPartitionCertificate2(
                    std::move(certificate));
            },
            "build_evidence"_a,
            "charts"_a,
            "projections"_a,
            "roots"_a,
            "cells"_a,
            "fibres"_a,
            "overlaps"_a,
            "seams"_a,
            "source_kind"_a,
            "source_payload"_a)
        .def_ro(
            "build_evidence",
            &EventPartitionCertificate2::build_evidence)
        .def_ro(
            "charts",
            &EventPartitionCertificate2::charts)
        .def_ro(
            "projections",
            &EventPartitionCertificate2::projections)
        .def_ro(
            "roots",
            &EventPartitionCertificate2::roots)
        .def_ro(
            "cells",
            &EventPartitionCertificate2::cells)
        .def_ro(
            "fibres",
            &EventPartitionCertificate2::fibres)
        .def_ro(
            "overlaps",
            &EventPartitionCertificate2::overlaps)
        .def_ro(
            "seams",
            &EventPartitionCertificate2::seams)
        .def_ro(
            "source_kind",
            &EventPartitionCertificate2::source_kind)
        .def_prop_ro(
            "source_payload",
            [](const EventPartitionCertificate2& certificate) {
                return as_bytes(
                    certificate.source_payload);
            })
        .def_prop_ro(
            "canonical_bytes",
            [](const EventPartitionCertificate2& certificate) {
                return as_bytes(
                    certificate.canonical_bytes);
            })
        .def_prop_ro(
            "canonical_digest",
            [](const EventPartitionCertificate2& certificate) {
                return as_bytes(
                    certificate.canonical_digest);
            });

    nb::class_<VerifiedEventPartition2>(
        m,
        "VerifiedEventPartition2")
        .def_ro(
            "verdict",
            &VerifiedEventPartition2::verdict)
        .def_ro(
            "partition",
            &VerifiedEventPartition2::partition);

    m.attr("BOUNDARY_EVENT_KINDS") = nb::make_tuple(
        "transverse",
        "tangent",
        "vertex",
        "overlap",
        "seam");

    nb::class_<BoundaryFeatureRecord2>(
        m,
        "BoundaryFeatureRecord2")
        .def_prop_ro(
            "support_kind",
            [](const BoundaryFeatureRecord2& record) {
                return record.support_kind;
            })
        .def_prop_ro(
            "support_coefficients",
            [](const BoundaryFeatureRecord2& record) {
                return coefficient_bytes(
                    record.support_coefficients);
            })
        .def_ro(
            "primitive_coefficients",
            &BoundaryFeatureRecord2::primitive_coefficients)
        .def_prop_ro(
            "support_id",
            [](const BoundaryFeatureRecord2& record) {
                return as_bytes(record.support_id);
            })
        .def_prop_ro(
            "source_exact",
            [](const BoundaryFeatureRecord2& record) {
                return as_bytes(record.source_exact);
            })
        .def_prop_ro(
            "target_exact",
            [](const BoundaryFeatureRecord2& record) {
                return as_bytes(record.target_exact);
            })
        .def_prop_ro(
            "source_vertex_id",
            [](const BoundaryFeatureRecord2& record) {
                return as_bytes(record.source_vertex_id);
            })
        .def_prop_ro(
            "target_vertex_id",
            [](const BoundaryFeatureRecord2& record) {
                return as_bytes(record.target_vertex_id);
            })
        .def_prop_ro(
            "source_incidence",
            [](const BoundaryFeatureRecord2& record) {
                return as_bytes(record.source_incidence);
            })
        .def_prop_ro(
            "target_incidence",
            [](const BoundaryFeatureRecord2& record) {
                return as_bytes(record.target_incidence);
            })
        .def_prop_ro(
            "material_side",
            [](const BoundaryFeatureRecord2& record) {
                return record.material_side;
            })
        .def_prop_ro(
            "trim_predicate",
            [](const BoundaryFeatureRecord2& record) {
                return as_bytes(record.trim_predicate);
            })
        .def_prop_ro(
            "feature_id",
            [](const BoundaryFeatureRecord2& record) {
                return as_bytes(record.feature_id);
            })
        .def_ro(
            "overlap_multiplicity",
            &BoundaryFeatureRecord2::overlap_multiplicity);

    nb::class_<BoundaryEvent2>(m, "BoundaryEvent2")
        .def_ro("kind", &BoundaryEvent2::kind)
        .def_prop_ro(
            "first_feature_id",
            [](const BoundaryEvent2& event) {
                return as_bytes(event.first_feature_id);
            })
        .def_prop_ro(
            "second_feature_id",
            [](const BoundaryEvent2& event) {
                return as_bytes(event.second_feature_id);
            })
        .def_prop_ro(
            "vertex_id",
            [](const BoundaryEvent2& event) {
                return as_bytes(event.vertex_id);
            })
        .def_prop_ro(
            "overlap_id",
            [](const BoundaryEvent2& event) {
                return as_bytes(event.overlap_id);
            })
        .def_prop_ro(
            "exact_overlap_record",
            [](const BoundaryEvent2& event) {
                return as_bytes(
                    event.exact_overlap_record);
            })
        .def_ro(
            "multiplicity",
            &BoundaryEvent2::multiplicity);

    m.def(
        "extract_boundary_records",
        &extract_boundary_records,
        "stock"_a);
    m.def(
        "classify_boundary_pair",
        &classify_boundary_pair,
        "stock"_a,
        "first_index"_a,
        "second_index"_a);
    m.def(
        "extract_boundary_events",
        &extract_boundary_events,
        "stock"_a);
    m.def(
        "exact_algebraic_backend_evidence",
        &exact_algebraic_backend_evidence);
    m.def(
        "parameter_charts",
        &parameter_charts);
    m.def(
        "verify_chart_coverage",
        &verify_chart_coverage,
        "charts"_a);
    m.def(
        "construct_pullback",
        &construct_pullback,
        "motion_kind"_a,
        "motion_data"_a,
        "support_kind"_a,
        "support_data"_a,
        "cutter_radius"_a,
        "center_chart"_a,
        "rim_chart"_a);
    m.def(
        "projection_from_grid",
        &projection_from_grid,
        "projection_id"_a,
        "coefficient_rows"_a,
        "degree_bound_id"_a);
    m.def(
        "partition_pullback_overlap",
        &partition_pullback_overlap,
        "projection"_a,
        "events"_a);
    m.def(
        "solve_trimmed_line_branches",
        &solve_trimmed_line_branches,
        "line_support"_a,
        "trim_start"_a,
        "trim_end"_a,
        "segment_motion"_a,
        "cutter_radius"_a,
        "rim_chart"_a);
    m.def(
        "project_regularization_vertex",
        [](const Stock2& stock,
           std::size_t first_index,
           std::size_t second_index,
           const nb::bytes& vertex_id) {
            return project_regularization_vertex(
                stock,
                first_index,
                second_index,
                from_bytes(vertex_id));
        },
        "stock"_a,
        "first_index"_a,
        "second_index"_a,
        "vertex_id"_a);
    m.def(
        "classify_ccw_orientation",
        &classify_ccw_orientation,
        "first_numerator"_a,
        "first_denominator"_a,
        "second_numerator"_a,
        "second_denominator"_a,
        "motion_parameter"_a);
    m.def(
        "partition_cap_crossings",
        &partition_cap_crossings,
        "first_numerator"_a,
        "first_denominator"_a,
        "second_numerator"_a,
        "second_denominator"_a,
        "cap_numerator"_a,
        "cap_denominator"_a,
        "event"_a);
    m.def(
        "verify_event_partition",
        static_cast<VerifiedEventPartition2 (*)(
            const EventPartitionCertificate2&)>(
            &verify_event_partition),
        "certificate"_a);
    m.def(
        "mutate_certificate_record",
        &mutate_certificate_record,
        "certificate"_a,
        "mutation"_a);
    m.def(
        "sha256_bytes",
        [](const nb::bytes& input) {
            return as_bytes(
                sha256_bytes(from_bytes(input)));
        },
        "input"_a);
    m.def(
        "partition_projections",
        &partition_projections,
        "projections"_a);
    m.def(
        "segment_pair_literal_signs",
        []() {
            return string_tuple(
                segment_pair_literal_signs());
        });
    m.def(
        "segment_rational_square_root_cases",
        []() {
            return string_tuple(
                segment_rational_square_root_cases());
        });
    m.def(
        "construct_segment_cell_stratum",
        &construct_segment_cell_stratum,
        "stock"_a,
        "source"_a,
        "witness_numerator"_a,
        "witness_denominator"_a);
    m.def(
        "construct_segment_event_partition",
        static_cast<SegmentEventPartition2 (*)(
            const Stock2&,
            double,
            double,
            double,
            double,
            double,
            double)>(
            &construct_segment_event_partition),
        "stock"_a,
        "x0"_a,
        "y0"_a,
        "x1"_a,
        "y1"_a,
        "tool_radius"_a,
        "cap_chord_ratio"_a);
    m.def(
        "verify_segment_event_partition",
        &verify_segment_event_partition,
        "stock"_a,
        "source"_a,
        "candidate"_a);
    m.def(
        "mutate_segment_event_partition",
        &mutate_segment_event_partition,
        "partition"_a,
        "mutation"_a);
}
