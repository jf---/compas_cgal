#include "reachable_domain_2.h"
#include "reachable_errors_2.h"

#include <string>
#include <vector>

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace {

nb::list bytes_sequence(const std::vector<std::string>& records)
{
    nb::list result;
    for (const std::string& record : records) {
        result.append(nb::bytes(record.data(), record.size()));
    }
    return result;
}

nb::bytes bytes_value(const std::string& value)
{
    return nb::bytes(value.data(), value.size());
}

} // namespace

NB_MODULE(_coverage_2, m)
{
    nb::exception<ReachableDomainConstructionError> reachable_error(
        m,
        "ReachableDomainConstructionError");
    nb::exception<InvalidReachableDomainInputError>(
        m,
        "InvalidReachableDomainInputError",
        reachable_error.ptr());
    nb::exception<ReachableArrangementTopologyError>(
        m,
        "ReachableArrangementTopologyError",
        reachable_error.ptr());
    nb::exception<PocketNotMachinableError>(
        m,
        "PocketNotMachinableError",
        reachable_error.ptr());
    nb::exception<ReachableMaterialContainmentError>(
        m,
        "ReachableMaterialContainmentError",
        reachable_error.ptr());
    nb::class_<ExactRegion2>(m, "ExactRegion2")
        .def("clone", &ExactRegion2::clone)
        .def("contains", &ExactRegion2::contains, "x"_a, "y"_a)
        .def("is_empty", &ExactRegion2::is_empty)
        .def("component_count", &ExactRegion2::component_count)
        .def("is_subset_of", &ExactRegion2::is_subset_of, "other"_a)
        .def("exactly_equals", &ExactRegion2::exactly_equals, "other"_a);

    nb::class_<ReachableDomainCertificate2>(
        m,
        "ReachableDomainCertificate2")
        .def_prop_ro(
            "strategy_version",
            [](const ReachableDomainCertificate2& certificate) {
                return bytes_value(certificate.strategy_version);
            })
        .def_prop_ro(
            "source_curve_records",
            [](const ReachableDomainCertificate2& certificate) {
                return bytes_sequence(
                    certificate.source_curve_records);
            })
        .def_prop_ro(
            "arrangement_vertex_records",
            [](const ReachableDomainCertificate2& certificate) {
                return bytes_sequence(
                    certificate.arrangement_vertex_records);
            })
        .def_prop_ro(
            "selected_cell_records",
            [](const ReachableDomainCertificate2& certificate) {
                return bytes_sequence(
                    certificate.selected_cell_records);
            })
        .def_prop_ro(
            "component_records",
            [](const ReachableDomainCertificate2& certificate) {
                return bytes_sequence(
                    certificate.component_records);
            })
        .def_ro(
            "exact_cell_selection",
            &ReachableDomainCertificate2::exact_cell_selection)
        .def_ro(
            "complete_source_provenance",
            &ReachableDomainCertificate2::complete_source_provenance)
        .def_ro(
            "reachable_subset_of_design",
            &ReachableDomainCertificate2::reachable_subset_of_design)
        .def(
            "matches_exact_inputs",
            &ReachableDomainCertificate2::matches_exact_inputs,
            "boundary"_a,
            "holes"_a,
            "tool_radius"_a);

    nb::class_<ReachableDomain2>(m, "ReachableDomain2")
        .def(
            nb::init<
                Eigen::Ref<const compas::RowMatrixXd>,
                const std::vector<compas::RowMatrixXd>&,
                double>(),
            "design_boundary"_a,
            "holes"_a,
            "tool_radius"_a)
        .def("design_region", &ReachableDomain2::design_region)
        .def("center_domain", &ReachableDomain2::center_domain)
        .def(
            "reachable_material",
            &ReachableDomain2::reachable_material)
        .def(
            "unreachable_residual",
            &ReachableDomain2::unreachable_residual)
        .def("certificate", &ReachableDomain2::certificate);

}
