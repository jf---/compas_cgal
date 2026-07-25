#include "compas.h"
#include "continuous_tea_2/boundary_events.h"
#include "continuous_tea_2/cap_partition.h"
#include "continuous_tea_2/circle_oracle.h"
#include "continuous_tea_2/event_certificate.h"
#include "continuous_tea_2/event_partition.h"
#include "continuous_tea_2/event_trace.h"
#include "continuous_tea_2/parameter_charts.h"
#include "continuous_tea_2/partition_certificate.h"
#include "continuous_tea_2/result.h"
#include "continuous_tea_2/sha256.h"
#include "exact_algebraic_1.h"

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

std::string from_bytes_handle(nb::handle value)
{
    if (!PyBytes_Check(value.ptr())) {
        throw nb::type_error(
            "canonical identity sequence entries must be bytes");
    }
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

std::vector<std::string> from_byte_sequence(
    const nb::sequence& values)
{
    std::vector<std::string> result;
    for (nb::handle value : values) {
        result.push_back(
            from_bytes_handle(value));
    }
    return result;
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

nb::object bytes_tuple(
    const std::vector<std::string>& values)
{
    nb::list result;
    for (const std::string& value : values) {
        result.append(as_bytes(value));
    }
    return nb::module_::import_("builtins")
        .attr("tuple")(result);
}

std::string continuous_tea_verdict_text(
    ContinuousTeaVerdict verdict)
{
    switch (verdict) {
    case ContinuousTeaVerdict::CERTIFIED:
        return "certified";
    case ContinuousTeaVerdict::CAP_EXCEEDED:
        return "cap_exceeded";
    case ContinuousTeaVerdict::UNRESOLVED_DEGENERACY:
        return "unresolved";
    }
    throw EventTraceVerificationError(
        "unknown continuous TEA verdict");
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
    nb::exception<TrimFilterError>(
        m,
        "TrimFilterError",
        substrate_error.ptr());
    nb::exception<EventPartitionVerificationError>(
        m,
        "EventPartitionVerificationError",
        substrate_error.ptr());
    nb::exception<EventTraceVerificationError>(
        m,
        "EventTraceVerificationError",
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
            "original_equations_rechecked",
            &PartitionEvent2::original_equations_rechecked)
        .def_ro(
            "orientation_rechecked",
            &PartitionEvent2::orientation_rechecked)
        .def_ro(
            "trim_disposition",
            &PartitionEvent2::trim_disposition);

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

    nb::class_<EventTraceEvent2>(
        m,
        "EventTraceEvent2")
        .def(
            "__init__",
            [](EventTraceEvent2* self,
               const nb::bytes& root_id,
               const nb::bytes& global_fibre_id,
               std::string kind,
               const nb::sequence& feature_ids,
               const nb::sequence& branch_ids,
               unsigned int multiplicity,
               std::string disposition,
               std::size_t motion_order) {
                new (self) EventTraceEvent2(
                    make_event_trace_event(
                        from_bytes(root_id),
                        from_bytes(global_fibre_id),
                        kind,
                        from_byte_sequence(feature_ids),
                        from_byte_sequence(branch_ids),
                        multiplicity,
                        disposition,
                        motion_order));
            },
            "root_id"_a,
            "global_fibre_id"_a,
            "kind"_a,
            "feature_ids"_a,
            "branch_ids"_a,
            "multiplicity"_a,
            "disposition"_a,
            "motion_order"_a)
        .def_prop_ro(
            "root_id",
            [](const EventTraceEvent2& event) {
                return as_bytes(event.root_id);
            })
        .def_prop_ro(
            "global_fibre_id",
            [](const EventTraceEvent2& event) {
                return as_bytes(event.global_fibre_id);
            })
        .def_ro("kind", &EventTraceEvent2::kind)
        .def_prop_ro(
            "feature_ids",
            [](const EventTraceEvent2& event) {
                return bytes_tuple(event.feature_ids);
            })
        .def_prop_ro(
            "branch_ids",
            [](const EventTraceEvent2& event) {
                return bytes_tuple(event.branch_ids);
            })
        .def_ro(
            "multiplicity",
            &EventTraceEvent2::multiplicity)
        .def_ro(
            "disposition",
            &EventTraceEvent2::disposition)
        .def_ro(
            "motion_order",
            &EventTraceEvent2::motion_order)
        .def_prop_ro(
            "canonical_bytes",
            [](const EventTraceEvent2& event) {
                return as_bytes(event.canonical_bytes);
            })
        .def_prop_ro(
            "canonical_id",
            [](const EventTraceEvent2& event) {
                return as_bytes(event.canonical_id);
            });

    nb::class_<EventTrace2>(m, "EventTrace2")
        .def_prop_ro(
            "exact_verdict",
            [](const EventTrace2& trace) {
                return continuous_tea_verdict_text(
                    trace.verdict);
            })
        .def_ro(
            "partition",
            &EventTrace2::partition)
        .def_ro("events", &EventTrace2::events)
        .def_ro(
            "motion_chart_id",
            &EventTrace2::motion_chart_id)
        .def_prop_ro(
            "motion_identity",
            [](const EventTrace2& trace) {
                return as_bytes(trace.motion_identity);
            })
        .def_prop_ro(
            "effective_cap_bytes",
            [](const EventTrace2& trace) {
                return as_bytes(
                    trace.effective_cap_bytes);
            })
        .def_ro(
            "whole_rim_disposition",
            &EventTrace2::whole_rim_disposition)
        .def_ro(
            "oracle_strategy_version",
            &EventTrace2::oracle_strategy_version)
        .def_ro(
            "event_cell_count",
            &EventTrace2::event_cell_count)
        .def_prop_ro(
            "canonical_bytes",
            [](const EventTrace2& trace) {
                return as_bytes(trace.canonical_bytes);
            })
        .def_prop_ro(
            "canonical_digest",
            [](const EventTrace2& trace) {
                return as_bytes(trace.canonical_digest);
            });

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
        "order_full_circle_events",
        &order_full_circle_events,
        "verified_partition"_a,
        "clockwise"_a,
        "events"_a);
    m.def(
        "build_event_trace",
        [](const EventPartitionCertificate2& partition,
           const std::string& motion_chart_id,
           const nb::bytes& motion_identity,
           const nb::bytes& effective_cap_bytes,
           ContinuousTeaVerdict verdict,
           const std::string& whole_rim_disposition,
           const std::string& oracle_strategy_version,
           std::vector<EventTraceEvent2> events) {
            return build_event_trace(
                partition,
                motion_chart_id,
                from_bytes(motion_identity),
                from_bytes(effective_cap_bytes),
                verdict,
                whole_rim_disposition,
                oracle_strategy_version,
                std::move(events));
        },
        "partition"_a,
        "motion_chart_id"_a,
        "motion_identity"_a,
        "effective_cap_bytes"_a,
        "verdict"_a,
        "whole_rim_disposition"_a,
        "oracle_strategy_version"_a,
        "events"_a);
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
}
