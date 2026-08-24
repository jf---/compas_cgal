#include "audit_identity_bindings_2.h"

#include "audit_classification_2.h"
#include "audit_policy_2.h"
#include "audit_request_identity_2.h"
#include "audit_stock_identity_2.h"
#include "audit_strategy_versions_2.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace {

nb::bytes digest_bytes(const std::string& bytes)
{
    return nb::bytes(bytes.data(), bytes.size());
}

std::size_t positive_center_limit(nb::handle value)
{
    if (!PyLong_CheckExact(value.ptr())) {
        throw AuditPolicyCenterCountLimitError(
            "audit policy center-count limit must be an exact positive integer");
    }
    const long long parsed = PyLong_AsLongLong(value.ptr());
    if (PyErr_Occurred()) {
        PyErr_Clear();
        throw AuditPolicyCenterCountLimitError(
            "audit policy center-count limit exceeds the native integer domain");
    }
    if (parsed <= 0
        || static_cast<unsigned long long>(parsed)
            > std::numeric_limits<std::size_t>::max()) {
        throw AuditPolicyCenterCountLimitError(
            "audit policy center-count limit must be positive and representable");
    }
    return static_cast<std::size_t>(parsed);
}

std::vector<NativeMotionDigest2> closed_motion_digests(const nb::tuple& motions)
{
    if (motions.size() == 0) {
        throw AuditNativeRequestMotionError(
            "native audit request requires at least one opaque motion");
    }
    std::vector<NativeMotionDigest2> result;
    result.reserve(motions.size());
    for (nb::handle motion : motions) {
        if (nb::isinstance<AuditSegmentMotion2>(motion)) {
            result.push_back(nb::cast<const AuditSegmentMotion2&>(motion).digest());
        } else if (nb::isinstance<AuditCircleMotion2>(motion)) {
            result.push_back(nb::cast<const AuditCircleMotion2&>(motion).digest());
        } else if (nb::isinstance<AuditArcMotion2>(motion)) {
            result.push_back(nb::cast<const AuditArcMotion2&>(motion).digest());
        } else if (nb::isinstance<AuditVerticalPlunge2>(motion)) {
            result.push_back(nb::cast<const AuditVerticalPlunge2&>(motion).digest());
        } else if (nb::isinstance<AuditVerticalRetract2>(motion)) {
            result.push_back(nb::cast<const AuditVerticalRetract2&>(motion).digest());
        } else if (nb::isinstance<AuditClearanceTransport2>(motion)) {
            result.push_back(nb::cast<const AuditClearanceTransport2&>(motion).digest());
        } else {
            throw AuditNativeRequestMotionError(
                "native audit request motion is outside the closed opaque domain");
        }
    }
    return result;
}

} // namespace

void register_audit_identity_2(nb::module_& module)
{
    nb::exception<AuditPolicyNonFiniteInputError>(
        module, "AuditPolicyNonFiniteInputError", PyExc_ValueError);
    nb::exception<AuditPolicyEngagementCapRangeError>(
        module, "AuditPolicyEngagementCapRangeError", PyExc_ValueError);
    nb::exception<AuditPolicyCapSurrogateMismatchError>(
        module, "AuditPolicyCapSurrogateMismatchError", PyExc_ValueError);
    nb::exception<AuditPolicyToolRadiusError>(
        module, "AuditPolicyToolRadiusError", PyExc_ValueError);
    nb::exception<AuditPolicyDepletionChordBoundError>(
        module, "AuditPolicyDepletionChordBoundError", PyExc_ValueError);
    nb::exception<AuditPolicyCenterCountLimitError>(
        module, "AuditPolicyCenterCountLimitError", PyExc_ValueError);
    nb::exception<AuditNativeStockIdentityError> stock_error(
        module, "AuditNativeStockIdentityError", PyExc_ValueError);
    nb::exception<AuditNativeStockShapeError>(
        module, "AuditNativeStockShapeError", stock_error.ptr());
    nb::exception<AuditNativeStockNonFiniteInputError>(
        module, "AuditNativeStockNonFiniteInputError", stock_error.ptr());
    nb::exception<AuditNativeStockRingError>(
        module, "AuditNativeStockRingError", stock_error.ptr());
    nb::exception<AuditNativeStockDuplicateHoleError>(
        module, "AuditNativeStockDuplicateHoleError", stock_error.ptr());
    nb::exception<AuditNativeRequestMotionError>(
        module, "AuditNativeRequestMotionError", PyExc_ValueError);

    nb::class_<AuditPolicy2>(module, "AuditPolicy2", nb::is_final())
        .def_prop_ro("digest", [](const AuditPolicy2& policy) {
            return digest_bytes(policy.digest().bytes());
        });
    nb::class_<AuditNativeStockIdentity2>(
        module, "AuditNativeStockIdentity2", nb::is_final())
        .def_prop_ro("canonical_bytes", [](const AuditNativeStockIdentity2& identity) {
            return digest_bytes(identity.canonical_bytes());
        })
        .def_prop_ro("digest", [](const AuditNativeStockIdentity2& identity) {
            return digest_bytes(identity.digest().bytes());
        });
    nb::class_<AuditNativeRequestIdentity2>(
        module, "AuditNativeRequestIdentity2", nb::is_final())
        .def_prop_ro("canonical_bytes", [](const AuditNativeRequestIdentity2& identity) {
            return digest_bytes(identity.canonical_bytes());
        })
        .def_prop_ro("digest", [](const AuditNativeRequestIdentity2& identity) {
            return digest_bytes(identity.digest().bytes());
        });

    module.def(
        "build_audit_policy",
        [](double tool_radius,
           double authored_radians,
           double supplied_chord_ratio,
           double depletion_chord_bound,
           nb::handle center_count_limit) {
            if (!std::isfinite(tool_radius)
                || !std::isfinite(authored_radians)
                || !std::isfinite(supplied_chord_ratio)
                || !std::isfinite(depletion_chord_bound)) {
                throw AuditPolicyNonFiniteInputError(
                    "audit policy scalar inputs must be finite binary64 values");
            }
            const AuditCapObservation2 cap = AuditCapObservation2::build(
                authored_radians, supplied_chord_ratio);
            return AuditPolicy2::build(
                cap,
                Epeck::FT(tool_radius),
                Epeck::FT(depletion_chord_bound),
                positive_center_limit(center_count_limit));
        },
        "tool_radius_mm"_a,
        "engagement_cap_radians"_a,
        "engagement_cap_chord_ratio"_a,
        "depletion_chord_bound_mm"_a,
        "center_count_limit"_a);
    module.def(
        "build_audit_native_stock_identity",
        &AuditNativeStockIdentity2::build,
        "boundary"_a,
        "holes"_a);
    module.def(
        "build_audit_native_request_identity",
        [](Eigen::Ref<const compas::RowMatrixXd> boundary,
           const std::vector<compas::RowMatrixXd>& holes,
           const AuditPolicy2& policy,
           const nb::tuple& motions) {
            return AuditNativeRequestIdentity2::build(
                AuditNativeStockIdentity2::build(boundary, holes),
                policy,
                closed_motion_digests(motions));
        },
        "boundary"_a,
        "holes"_a,
        "policy"_a,
        "motions"_a);
    module.def("audit_native_decision_contract_version", []() {
        return digest_bytes(audit_native_decision_contract_version());
    });
    module.def("audit_native_depletion_contract_version", []() {
        return digest_bytes(audit_native_depletion_contract_version());
    });
}
