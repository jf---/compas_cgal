#include "coverage_2.h"
#include "cutter_centre_domain_2.h"
#include "reachable_arrangement_2.h"
#include "reachable_domain_2.h"
#include "reachable_errors_2.h"
#include "reachable_input_2.h"
#include "reachable_material_predicate_2.h"

#include <algorithm>
#include <string>
#include <tuple>
#include <vector>

#include <CGAL/number_utils.h>

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
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

std::tuple<double, double> reporting_point(const ReachPoint& point)
{
    return {
        CGAL::to_double(point.x()),
        CGAL::to_double(point.y()),
    };
}

const ReachableBoundaryCycle2& outer_center_boundary(
    const ReachableArrangement2& arrangement)
{
    const auto found = std::find_if(
        arrangement.center_boundary_cycles.begin(),
        arrangement.center_boundary_cycles.end(),
        [](const ReachableBoundaryCycle2& cycle) {
            return cycle.orientation == CGAL::COUNTERCLOCKWISE;
        });
    if (found == arrangement.center_boundary_cycles.end()) {
        throw ReachableArrangementTopologyError(
            "reachable arrangement exposes no counterclockwise center boundary");
    }
    return *found;
}

std::tuple<double, double> arc_center_reporting_point(
    const ReachableBoundaryCurve2& primitive)
{
    if (!primitive.curve.is_circular()) {
        throw ReachableArrangementTopologyError(
            "linear center-boundary primitive has no arc center");
    }
    const ReachKernelPoint center =
        primitive.curve.supporting_circle().center();
    return {
        CGAL::to_double(center.x()),
        CGAL::to_double(center.y()),
    };
}

double arc_reporting_radius(
    const ReachableBoundaryCurve2& primitive)
{
    if (!primitive.curve.is_circular()) {
        throw ReachableArrangementTopologyError(
            "linear center-boundary primitive has no arc radius");
    }
    return CGAL::to_double(
        CGAL::sqrt(
            primitive.curve.supporting_circle().squared_radius()));
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
    nb::exception<ReachableMaterialPredicateGeometryError>(
        m,
        "ReachableMaterialPredicateGeometryError",
        reachable_error.ptr());

    nb::class_<ReachPoint>(
        m,
        "WorldXYBoundaryPointMm",
        "Native exact boundary contact in world XY, measured in millimetres.")
        .def_prop_ro(
            "reporting_xy_mm",
            &reporting_point,
            "Approximate coordinates for reporting only; never geometric input.")
        .def(
            "__eq__",
            [](const ReachPoint& left, const ReachPoint& right) {
                return ReachTraits().equal_2_object()(left, right);
            },
            nb::is_operator())
        .def(
            "__ne__",
            [](const ReachPoint& left, const ReachPoint& right) {
                return !ReachTraits().equal_2_object()(left, right);
            },
            nb::is_operator());

    nb::class_<ReachableBoundaryCurve2>(
        m,
        "ReachableBoundaryPrimitive2")
        .def_prop_ro(
            "kind",
            [](const ReachableBoundaryCurve2& primitive) {
                return primitive.curve.is_linear()
                    ? std::string("line")
                    : std::string("arc");
            })
        .def_prop_ro(
            "start",
            [](const ReachableBoundaryCurve2& primitive) -> ReachPoint {
                return primitive.curve.source();
            })
        .def_prop_ro(
            "end",
            [](const ReachableBoundaryCurve2& primitive) -> ReachPoint {
                return primitive.curve.target();
            })
        .def_prop_ro(
            "start_mm",
            [](const ReachableBoundaryCurve2& primitive) {
                return reporting_point(primitive.curve.source());
            })
        .def_prop_ro(
            "end_mm",
            [](const ReachableBoundaryCurve2& primitive) {
                return reporting_point(primitive.curve.target());
            })
        .def_prop_ro(
            "source_piece_records",
            [](const ReachableBoundaryCurve2& primitive) {
                return bytes_sequence(primitive.source_piece_ids);
            })
        .def_prop_ro(
            "arc_center_mm",
            &arc_center_reporting_point)
        .def_prop_ro(
            "arc_radius_mm",
            &arc_reporting_radius)
        .def_prop_ro(
            "arc_counterclockwise",
            [](const ReachableBoundaryCurve2& primitive) {
                if (!primitive.curve.is_circular()) {
                    throw ReachableArrangementTopologyError(
                        "linear center-boundary primitive has no arc orientation");
                }
                return primitive.curve.orientation()
                    == CGAL::COUNTERCLOCKWISE;
            });

    nb::class_<ReachableBoundaryCycle2>(
        m,
        "ReachableBoundaryCycle2")
        .def_prop_ro(
            "counterclockwise",
            [](const ReachableBoundaryCycle2& cycle) {
                return cycle.orientation == CGAL::COUNTERCLOCKWISE;
            })
        .def_ro("primitives", &ReachableBoundaryCycle2::curves)
        .def(
            "ccw_transition",
            [](const ReachableBoundaryCycle2& cycle,
               const ReachPoint& start,
               const ReachPoint& end) {
                return reachable_ccw_transition(cycle, start, end).curves;
            },
            "start"_a,
            "end"_a,
            "Traverse positive CCW boundary progress between exact native contacts.");

    m.def(
        "build_center_boundary_cycle",
        [](Eigen::Ref<const compas::RowMatrixXd> design_boundary,
           const std::vector<compas::RowMatrixXd>& holes,
           double tool_radius) {
            ReachableArrangement2 arrangement =
                build_reachable_arrangement(
                    canonical_reach_input(
                        design_boundary,
                        holes,
                        tool_radius));
            return outer_center_boundary(arrangement);
        },
        "design_boundary"_a,
        "holes"_a,
        "tool_radius"_a);

    nb::exception<InvalidCoverageGeometryError>(
        m,
        "InvalidCoverageGeometryError",
        PyExc_RuntimeError);
    nb::exception<CoverageTransitionError>(
        m,
        "CoverageTransitionError",
        PyExc_RuntimeError);
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

    nb::class_<CutterCentreDomain2>(m, "CutterCentreDomain2")
        .def_static(
            "build",
            &CutterCentreDomain2::build,
            "design_boundary"_a,
            "holes"_a,
            "tool_radius"_a)
        .def("contains", &CutterCentreDomain2::contains, "x"_a, "y"_a);

    nb::class_<ReachableMaterialPredicate2>(
        m,
        "ReachableMaterialPredicate2")
        .def_static(
            "build",
            &ReachableMaterialPredicate2::build,
            "design_boundary"_a,
            "holes"_a,
            "tool_radius"_a)
        .def("contains", &ReachableMaterialPredicate2::contains, "x"_a, "y"_a);

    nb::class_<CoverageSweepRecord2>(m, "CoverageSweepRecord2")
        .def_prop_ro(
            "strategy_version",
            [](const CoverageSweepRecord2& record) {
                return bytes_value(record.strategy_version);
            })
        .def_prop_ro(
            "structural_record",
            [](const CoverageSweepRecord2& record) {
                return bytes_value(record.structural_record);
            })
        .def_ro("segment", &CoverageSweepRecord2::segment)
        .def_ro("center_x", &CoverageSweepRecord2::center_x)
        .def_ro("center_y", &CoverageSweepRecord2::center_y)
        .def_ro("first_x", &CoverageSweepRecord2::first_x)
        .def_ro("first_y", &CoverageSweepRecord2::first_y)
        .def_ro("tool_radius", &CoverageSweepRecord2::tool_radius)
        .def(
            "matches_exact_segment",
            &CoverageSweepRecord2::matches_exact_segment,
            "x0"_a,
            "y0"_a,
            "x1"_a,
            "y1"_a,
            "tool_radius"_a)
        .def(
            "matches_exact_full_circle",
            &CoverageSweepRecord2::matches_exact_full_circle,
            "center_x"_a,
            "center_y"_a,
            "phase_x"_a,
            "phase_y"_a,
            "tool_radius"_a);

    nb::class_<Coverage2>(m, "Coverage2")
        .def(
            nb::init<
                const ExactRegion2&,
                double,
                double,
                double>(),
            "reachable_material"_a,
            "precleared_x"_a,
            "precleared_y"_a,
            "precleared_radius"_a)
        .def("clone", &Coverage2::clone)
        .def(
            "add_segment_sweep",
            &Coverage2::add_segment_sweep,
            "x0"_a,
            "y0"_a,
            "x1"_a,
            "y1"_a,
            "tool_radius"_a)
        .def(
            "add_full_circle_sweep",
            &Coverage2::add_full_circle_sweep,
            "center_x"_a,
            "center_y"_a,
            "phase_x"_a,
            "phase_y"_a,
            "tool_radius"_a)
        .def("residual_is_empty", &Coverage2::residual_is_empty)
        .def(
            "residual_component_count",
            &Coverage2::residual_component_count)
        .def("residual", &Coverage2::residual)
        .def(
            "accumulated_sweeps",
            &Coverage2::accumulated_sweeps)
        .def_prop_ro(
            "residual_component_records",
            [](const Coverage2& coverage) {
                return bytes_sequence(
                    coverage.residual_component_records());
            })
        .def_prop_ro(
            "sweep_records",
            [](const Coverage2& coverage) {
                return bytes_sequence(coverage.sweep_records());
            })
        .def(
            "exact_residual_relation",
            &Coverage2::exact_residual_relation)
        .def_prop_ro(
            "strategy_version",
            [](const Coverage2& coverage) {
                return bytes_value(coverage.strategy_version());
            });
}
