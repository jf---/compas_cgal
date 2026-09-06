#include "audit_classification_2.h"
#include "audit_identity_bindings_2.h"
#include "audit_replay_bindings_2.h"
#include "engagement_2_bindings.h"
#include "stock_2.h"
#include "stock_local_2.h"

#include <cmath>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace nb::literals;

void register_held_disk_contour_2(nb::module_& module);

namespace {

std::vector<ExactCenterParameter2> to_exact_center_parameters(
    const std::vector<std::tuple<int, std::size_t, std::size_t>>& parameters)
{
    std::vector<ExactCenterParameter2> result;
    result.reserve(parameters.size());
    for (const auto& [chart, numerator, denominator] : parameters) {
        result.push_back({chart, numerator, denominator});
    }
    return result;
}

} // namespace

NB_MODULE(_stock_2, m)
{
    register_held_disk_contour_2(m);
    register_audit_classification_2(m);
    register_audit_identity_2(m);
    register_audit_replay_2(m);
    nb::exception<ExactDepletionConstructionError> construction_error(
        m,
        "ExactDepletionConstructionError");
    nb::exception<ExactDepletionCenterLimitError>(
        m,
        "ExactDepletionCenterLimitError",
        construction_error.ptr());
    nb::exception<ExactArcDepletionPolicyError>(
        m,
        "ExactArcDepletionPolicyError",
        construction_error.ptr());
    nb::exception<ExactArcForgedTraceError>(
        m,
        "ExactArcForgedTraceError",
        construction_error.ptr());
    nb::exception<NonFiniteExactArcDepletionInputError>(
        m,
        "NonFiniteExactArcDepletionInputError",
        PyExc_ValueError);

    // Argument faults at the double boundary: ValueError-derived so they read the
    // same way as every other malformed-input rejection in this module.
    nb::exception<InvalidAnnulusRadiiError>(m, "InvalidAnnulusRadiiError", PyExc_ValueError);
    nb::exception<NonFiniteAnnulusInputError>(m, "NonFiniteAnnulusInputError", PyExc_ValueError);
    nb::exception<NonFiniteCapsuleInputError>(m, "NonFiniteCapsuleInputError", PyExc_ValueError);

    // Not an argument fault: the quad capsule's exact half-width certificate did
    // not hold. Registered explicitly under RuntimeError -- nanobind's default
    // base is Exception, which would leave "broken invariant" and "malformed
    // input" with no common ancestor a caller could tell apart.
    nb::exception<CapsuleQuadCertificateError>(m, "CapsuleQuadCertificateError", PyExc_RuntimeError);

    // Not an argument fault: a broken internal invariant of the local depletion
    // path, surfaced as a RuntimeError so it can never be confused with -- or
    // caught alongside -- a malformed-input rejection.
    nb::exception<LocalDepletionEscapedError>(m, "LocalDepletionEscapedError");

    nb::class_<DepletionTrace>(m, "DepletionTrace")
        .def_prop_ro("center_count", [](const DepletionTrace& trace) {
            return trace.center_count;
        })
        .def_prop_ro("center_parameters", [](const DepletionTrace& trace) {
            std::vector<std::tuple<int, std::size_t, std::size_t>> result;
            result.reserve(trace.center_parameters.size());
            for (const ExactCenterParameter2& parameter : trace.center_parameters) {
                result.emplace_back(
                    parameter.chart,
                    parameter.numerator,
                    parameter.denominator);
            }
            return result;
        })
        .def(
            "matches_exact_inputs",
            [](const DepletionTrace& trace,
               double expected_tool_radius,
               double expected_max_chord,
               std::size_t expected_center_count_limit) {
                return trace.matches_exact_inputs(
                    Epeck::FT(expected_tool_radius),
                    Epeck::FT(expected_max_chord),
                    expected_center_count_limit);
            },
            "expected_tool_radius"_a,
            "expected_max_chord"_a,
            "expected_center_count_limit"_a)
        .def_prop_ro("strategy_version", [](const DepletionTrace& trace) {
            return nb::bytes(
                trace.strategy_version.data(),
                trace.strategy_version.size());
        })
        .def_ro("cyclic", &DepletionTrace::cyclic)
        .def_ro("exact_incidence", &DepletionTrace::exact_incidence)
        .def_ro("exact_parameters_in_range", &DepletionTrace::exact_parameters_in_range)
        .def_ro("exact_anchors_present", &DepletionTrace::exact_anchors_present)
        .def_ro("exact_removal_radius_valid", &DepletionTrace::exact_removal_radius_valid)
        .def_ro("exact_chord_bound_holds", &DepletionTrace::exact_chord_bound_holds)
        .def_ro("exact_seam_chord_bound_holds", &DepletionTrace::exact_seam_chord_bound_holds);

    nb::class_<ExactArcDepletionTrace2>(m, "ExactArcDepletionTrace2")
        .def_prop_ro("center_count", [](const ExactArcDepletionTrace2& trace) {
            return trace.parameters().size();
        })
        .def(
            "matches_exact_inputs",
            [](const ExactArcDepletionTrace2& trace,
               double tool_radius,
               double max_chord,
               std::int64_t center_count_limit) {
                if (!std::isfinite(tool_radius)
                    || !std::isfinite(max_chord)) {
                    throw NonFiniteExactArcDepletionInputError(
                        "exact arc matcher inputs must be finite");
                }
                if (center_count_limit <= 0) {
                    throw ExactDepletionCenterLimitError(
                        "exact arc center-count limit must be positive");
                }
                return trace.matches_exact_inputs(
                    Epeck::FT(tool_radius),
                    Epeck::FT(max_chord),
                    static_cast<std::size_t>(center_count_limit));
            },
            "tool_radius"_a,
            "max_chord"_a,
            "center_count_limit"_a)
        .def(
            "matches_motion",
            &ExactArcDepletionTrace2::matches_motion,
            "motion"_a)
        .def_prop_ro("canonical_bytes", [](const ExactArcDepletionTrace2& trace) {
            return nb::bytes(
                trace.canonical_bytes().data(),
                trace.canonical_bytes().size());
        })
        .def_prop_ro("digest", [](const ExactArcDepletionTrace2& trace) {
            return nb::bytes(
                trace.digest().bytes().data(),
                trace.digest().bytes().size());
        })
        .def_prop_ro("strategy_version", [](const ExactArcDepletionTrace2& trace) {
            return nb::bytes(
                trace.strategy_version().data(),
                trace.strategy_version().size());
        })
        .def_prop_ro("cyclic", &ExactArcDepletionTrace2::cyclic);

    nb::class_<Stock2>(m, "Stock2")
        .def(nb::init<Eigen::Ref<const compas::RowMatrixXd>,
                      const std::vector<compas::RowMatrixXd>&>(),
             "boundary"_a, "holes"_a)
        .def("contains", &Stock2::contains, "x"_a, "y"_a)
        .def("is_empty", &Stock2::is_empty)
        .def("clone", &Stock2::clone)
        .def("is_subset_of", &Stock2::is_subset_of, "other"_a)
        .def("exactly_equals", &Stock2::exactly_equals, "other"_a)
        .def("subtract_capsule", &Stock2::subtract_capsule,
             "x0"_a, "y0"_a, "x1"_a, "y1"_a, "radius"_a)
        .def("subtract_capsule_quad", &Stock2::subtract_capsule_quad,
             "x0"_a, "y0"_a, "x1"_a, "y1"_a, "radius"_a)
        .def("subtract_arc_sweep", &Stock2::subtract_arc_sweep,
             "cx"_a, "cy"_a, "sx"_a, "sy"_a, "ex"_a, "ey"_a, "cw"_a, "tool_radius"_a)
        .def(
            "subtract_exact_segment",
            [](Stock2& stock,
               double x0,
               double y0,
               double x1,
               double y1,
               double tool_radius,
               double max_chord,
               std::size_t center_count_limit) {
                return stock.subtract_exact_segment(
                    {EPoint(x0, y0), EPoint(x1, y1)},
                    Epeck::FT(tool_radius),
                    Epeck::FT(max_chord),
                    center_count_limit);
            },
            "x0"_a,
            "y0"_a,
            "x1"_a,
            "y1"_a,
            "tool_radius"_a,
            "max_chord"_a,
            "center_count_limit"_a)
        .def(
            "subtract_exact_full_circle",
            [](Stock2& stock,
               double cx,
               double cy,
               double phase_x,
               double phase_y,
               bool clockwise,
               double tool_radius,
               double max_chord,
               std::size_t center_count_limit) {
                return stock.subtract_exact_full_circle(
                    {EPoint(cx, cy), EVector(phase_x, phase_y), clockwise},
                    Epeck::FT(tool_radius),
                    Epeck::FT(max_chord),
                    center_count_limit);
            },
            "cx"_a,
            "cy"_a,
            "phase_x"_a,
            "phase_y"_a,
            "clockwise"_a,
            "tool_radius"_a,
            "max_chord"_a,
            "center_count_limit"_a)
        .def(
            "subtract_exact_arc",
            [](Stock2& stock,
               const AuditArcMotion2& motion,
               double tool_radius,
               double max_chord,
               std::int64_t center_count_limit) {
                if (!std::isfinite(tool_radius)
                    || !std::isfinite(max_chord)) {
                    throw NonFiniteExactArcDepletionInputError(
                        "exact arc tool radius and chord bound must be finite");
                }
                if (center_count_limit <= 0) {
                    throw ExactDepletionCenterLimitError(
                        "exact arc center-count limit must be positive");
                }
                return stock.subtract_exact_arc(
                    motion,
                    Epeck::FT(tool_radius),
                    Epeck::FT(max_chord),
                    static_cast<std::size_t>(center_count_limit));
            },
            "motion"_a,
            "tool_radius"_a,
            "max_chord"_a,
            "center_count_limit"_a)
        .def("subtract_disk", &Stock2::subtract_disk, "cx"_a, "cy"_a, "radius"_a)
        .def("subtract_annulus", &Stock2::subtract_annulus,
             "cx"_a, "cy"_a, "inner_radius"_a, "outer_radius"_a)
        .def("subtract_circle_sweep", &Stock2::subtract_circle_sweep,
             "cx"_a, "cy"_a, "guide_radius"_a, "tool_radius"_a)
        .def("subtract_disk_local", &Stock2::subtract_disk_local,
             "cx"_a, "cy"_a, "radius"_a)
        .def("subtract_annulus_local", &Stock2::subtract_annulus_local,
             "cx"_a, "cy"_a, "inner_radius"_a, "outer_radius"_a)
        .def("subtract_arc_sweep_local", &Stock2::subtract_arc_sweep_local,
             "cx"_a, "cy"_a, "sx"_a, "sy"_a, "ex"_a, "ey"_a, "cw"_a, "tool_radius"_a)
        .def("representation_is_valid", &Stock2::representation_is_valid)
        .def(
            "arrangement_stats",
            [](const Stock2& stock) {
                const Stock2::ArrangementStats stats = stock.arrangement_stats();
                return std::make_tuple(stats.vertices, stats.halfedges, stats.faces);
            })
        .def(
            "coordinate_digits",
            [](const Stock2& stock) {
                const Stock2::CoordinateDigits digits = stock.coordinate_digits();
                return std::make_tuple(digits.max_digits, digits.mean_digits, digits.sampled);
            });

    m.def(
        "exact_segment_point_is_incident",
        [](double x0, double y0, double x1, double y1, double px, double py) {
            return exact_segment_point_is_incident(
                {EPoint(x0, y0), EPoint(x1, y1)},
                EPoint(px, py));
        },
        "x0"_a,
        "y0"_a,
        "x1"_a,
        "y1"_a,
        "px"_a,
        "py"_a);
    m.def(
        "exact_circle_point_is_incident",
        [](double cx,
           double cy,
           double phase_x,
           double phase_y,
           double px,
           double py) {
            return exact_circle_point_is_incident(
                {EPoint(cx, cy), EVector(phase_x, phase_y), false},
                EPoint(px, py));
        },
        "cx"_a,
        "cy"_a,
        "phase_x"_a,
        "phase_y"_a,
        "px"_a,
        "py"_a);
    m.def(
        "exact_segment_structural_density_holds",
        [](double x0,
           double y0,
           double x1,
           double y1,
           double max_chord,
           const std::vector<std::tuple<int, std::size_t, std::size_t>>& parameters) {
            return exact_segment_structural_density_holds(
                {EPoint(x0, y0), EPoint(x1, y1)},
                Epeck::FT(max_chord),
                to_exact_center_parameters(parameters));
        },
        "x0"_a,
        "y0"_a,
        "x1"_a,
        "y1"_a,
        "max_chord"_a,
        "parameters"_a);
    m.def(
        "exact_full_circle_structural_density_holds",
        [](double cx,
           double cy,
           double phase_x,
           double phase_y,
           bool clockwise,
           double max_chord,
           const std::vector<std::tuple<int, std::size_t, std::size_t>>& parameters) {
            return exact_full_circle_structural_density_holds(
                {EPoint(cx, cy), EVector(phase_x, phase_y), clockwise},
                Epeck::FT(max_chord),
                to_exact_center_parameters(parameters));
        },
        "cx"_a,
        "cy"_a,
        "phase_x"_a,
        "phase_y"_a,
        "clockwise"_a,
        "max_chord"_a,
        "parameters"_a);
    m.def(
        "exact_segment_undercover_holds",
        [](double x0,
           double y0,
           double x1,
           double y1,
           double exact_length,
           double tool_radius,
           double max_chord,
           std::size_t center_count_limit) {
            return exact_segment_undercover_holds(
                {EPoint(x0, y0), EPoint(x1, y1)},
                Epeck::FT(exact_length),
                Epeck::FT(tool_radius),
                Epeck::FT(max_chord),
                center_count_limit);
        },
        "x0"_a,
        "y0"_a,
        "x1"_a,
        "y1"_a,
        "exact_length"_a,
        "tool_radius"_a,
        "max_chord"_a,
        "center_count_limit"_a);
    m.def(
        "exact_full_circle_undercover_holds",
        [](double cx,
           double cy,
           double phase_x,
           double phase_y,
           double guide_radius,
           double tool_radius,
           double max_chord,
           std::size_t center_count_limit) {
            return exact_full_circle_undercover_holds(
                {EPoint(cx, cy), EVector(phase_x, phase_y), false},
                Epeck::FT(guide_radius),
                Epeck::FT(tool_radius),
                Epeck::FT(max_chord),
                center_count_limit);
        },
        "cx"_a,
        "cy"_a,
        "phase_x"_a,
        "phase_y"_a,
        "guide_radius"_a,
        "tool_radius"_a,
        "max_chord"_a,
        "center_count_limit"_a);
    m.def(
        "exact_segment_induction_holds",
        [](const Stock2& initial,
           double x0,
           double y0,
           double x1,
           double y1,
           double exact_length,
           double tool_radius,
           double max_chord,
           std::size_t center_count_limit) {
            return exact_segment_induction_holds(
                initial,
                {EPoint(x0, y0), EPoint(x1, y1)},
                Epeck::FT(exact_length),
                Epeck::FT(tool_radius),
                Epeck::FT(max_chord),
                center_count_limit);
        },
        "initial"_a,
        "x0"_a,
        "y0"_a,
        "x1"_a,
        "y1"_a,
        "exact_length"_a,
        "tool_radius"_a,
        "max_chord"_a,
        "center_count_limit"_a);
    m.def(
        "exact_full_circle_induction_holds",
        [](const Stock2& initial,
           double cx,
           double cy,
           double phase_x,
           double phase_y,
           double guide_radius,
           double tool_radius,
           double max_chord,
           std::size_t center_count_limit) {
            return exact_full_circle_induction_holds(
                initial,
                {EPoint(cx, cy), EVector(phase_x, phase_y), false},
                Epeck::FT(guide_radius),
                Epeck::FT(tool_radius),
                Epeck::FT(max_chord),
                center_count_limit);
        },
        "initial"_a,
        "cx"_a,
        "cy"_a,
        "phase_x"_a,
        "phase_y"_a,
        "guide_radius"_a,
        "tool_radius"_a,
        "max_chord"_a,
        "center_count_limit"_a);
    m.def("exact_depletion_strategy_version", []() {
        const std::string& version = exact_depletion_strategy_version();
        return nb::bytes(version.data(), version.size());
    });
    m.def("exact_entry_depletion_strategy_version", []() {
        static const std::string version = "exact-precleared-entry-disk-v1";
        return nb::bytes(version.data(), version.size());
    });

    register_engagement(m);
}
