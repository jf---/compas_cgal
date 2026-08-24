#include "audit_classification_2.h"

#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>

namespace nb = nanobind;
using namespace nb::literals;

void register_audit_classification_2(nb::module_& module)
{
    nb::exception<AuditNonFiniteInputError>(
        module, "AuditNonFiniteInputError", PyExc_ValueError);
    nb::exception<AuditInvalidPlaneError>(
        module, "AuditInvalidPlaneError", PyExc_ValueError);
    nb::exception<AuditUnsupportedGeometryError> unsupported(
        module, "AuditUnsupportedGeometryError", PyExc_ValueError);
    nb::exception<AuditOffPlaneError>(
        module, "AuditOffPlaneError", unsupported.ptr());
    nb::exception<AuditContradictoryRoleError>(
        module, "AuditContradictoryRoleError", PyExc_ValueError);
    nb::exception<AuditContradictoryOrientationError>(
        module, "AuditContradictoryOrientationError", PyExc_ValueError);

    nb::class_<AuditSegmentMotion2>(
        module, "AuditSegmentMotion2", nb::is_final())
        .def_prop_ro("digest", [](const AuditSegmentMotion2& motion) {
            return nb::bytes(
                motion.digest().bytes().data(),
                motion.digest().bytes().size());
        });
    nb::class_<AuditCircleMotion2>(
        module, "AuditCircleMotion2", nb::is_final())
        .def_prop_ro("digest", [](const AuditCircleMotion2& motion) {
            return nb::bytes(
                motion.digest().bytes().data(),
                motion.digest().bytes().size());
        });
    nb::class_<AuditArcMotion2>(
        module, "AuditArcMotion2", nb::is_final())
        .def_prop_ro("digest", [](const AuditArcMotion2& motion) {
            return nb::bytes(
                motion.digest().bytes().data(),
                motion.digest().bytes().size());
        });
    nb::class_<AuditVerticalPlunge2>(
        module, "AuditVerticalPlunge2", nb::is_final())
        .def_prop_ro("digest", [](const AuditVerticalPlunge2& motion) {
            return nb::bytes(
                motion.digest().bytes().data(),
                motion.digest().bytes().size());
        });
    nb::class_<AuditVerticalRetract2>(
        module, "AuditVerticalRetract2", nb::is_final())
        .def_prop_ro("digest", [](const AuditVerticalRetract2& motion) {
            return nb::bytes(
                motion.digest().bytes().data(),
                motion.digest().bytes().size());
        });
    nb::class_<AuditClearanceTransport2>(
        module, "AuditClearanceTransport2", nb::is_final())
        .def_prop_ro("digest", [](const AuditClearanceTransport2& motion) {
            return nb::bytes(
                motion.digest().bytes().data(),
                motion.digest().bytes().size());
        });

    module.def(
        "audit_arc_phase_strategy_version",
        []() {
            return nb::bytes(
                audit_arc_motion_strategy_version().data(),
                audit_arc_motion_strategy_version().size());
        });

    module.def(
        "classify_audit_line", &classify_audit_line,
        "start"_a, "end"_a, "cut_z"_a, "clearance_z"_a,
        "operation_role"_a);
    module.def(
        "classify_audit_circle", &classify_audit_circle,
        "center"_a, "xaxis"_a, "yaxis"_a, "radius"_a,
        "clockwise"_a, "cut_z"_a, "clearance_z"_a,
        "operation_role"_a);
    module.def(
        "classify_audit_arc", &classify_audit_arc,
        "center"_a, "xaxis"_a, "yaxis"_a, "radius"_a,
        "start_angle"_a, "end_angle"_a, "clockwise"_a,
        "cut_z"_a, "clearance_z"_a, "operation_role"_a);
}
