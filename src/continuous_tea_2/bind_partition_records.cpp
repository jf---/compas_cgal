#include "bind_partition_records.h"

#include "partition_certificate.h"

#include <string>
#include <vector>

#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;

namespace {

nb::bytes as_bytes(const std::string& value)
{
    return nb::bytes(value.data(), value.size());
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

void bind_partition_records(nb::module_& module)
{
    nb::class_<ActiveBoundaryBranch2>(
        module,
        "ActiveBoundaryBranch2")
        .def_prop_ro(
            "branch_id",
            [](const ActiveBoundaryBranch2& branch) {
                return as_bytes(branch.branch_id);
            })
        .def_prop_ro(
            "feature_id",
            [](const ActiveBoundaryBranch2& branch) {
                return as_bytes(branch.feature_id);
            })
        .def_prop_ro(
            "support_id",
            [](const ActiveBoundaryBranch2& branch) {
                return as_bytes(branch.support_id);
            })
        .def_prop_ro(
            "trim_id",
            [](const ActiveBoundaryBranch2& branch) {
                return as_bytes(branch.trim_id);
            })
        .def_ro(
            "chart_id",
            &ActiveBoundaryBranch2::chart_id)
        .def_ro(
            "sheet_ordinal",
            &ActiveBoundaryBranch2::sheet_ordinal)
        .def_prop_ro(
            "root_id",
            [](const ActiveBoundaryBranch2& branch) {
                return as_bytes(branch.root_id);
            });

    nb::class_<EventFibre2>(module, "EventFibre2")
        .def_prop_ro(
            "root_id",
            [](const EventFibre2& fibre) {
                return as_bytes(fibre.root_id);
            })
        .def_ro("events", &EventFibre2::events)
        .def_prop_ro(
            "local_root_ids",
            [](const EventFibre2& fibre) {
                return bytes_tuple(
                    fibre.local_root_ids);
            })
        .def_prop_ro(
            "seam_id",
            [](const EventFibre2& fibre) {
                return as_bytes(fibre.seam_id);
            })
        .def_ro(
            "left_active_branches",
            &EventFibre2::left_active_branches)
        .def_ro(
            "right_active_branches",
            &EventFibre2::right_active_branches)
        .def_ro(
            "ccw_direction",
            &EventFibre2::ccw_direction)
        .def_ro(
            "cw_direction",
            &EventFibre2::cw_direction);

    nb::class_<ProjectionRecord2>(
        module,
        "ProjectionRecord2")
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
            })
        .def_prop_ro(
            "signed_predicate_coefficients",
            [](const ProjectionRecord2& projection) {
                return string_tuple(
                    projection
                        .signed_predicate_coefficients);
            });

    nb::class_<OverlapInterval2>(
        module,
        "OverlapInterval2")
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
            &OverlapInterval2::orientation_disposition)
        .def_prop_ro(
            "feature_id",
            [](const OverlapInterval2& overlap) {
                return as_bytes(overlap.feature_id);
            })
        .def_prop_ro(
            "support_id",
            [](const OverlapInterval2& overlap) {
                return as_bytes(overlap.support_id);
            })
        .def_prop_ro(
            "trim_id",
            [](const OverlapInterval2& overlap) {
                return as_bytes(overlap.trim_id);
            })
        .def_prop_ro(
            "branch_id",
            [](const OverlapInterval2& overlap) {
                return as_bytes(overlap.branch_id);
            });
}
