#include "audit_replay_bindings_2.h"

#include "audit_depletion_witness_2.h"
#include "audit_replay_2.h"

#include <string>
#include <vector>

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace {

nb::bytes as_bytes(const std::string& value)
{
    return nb::bytes(value.data(), value.size());
}

std::string from_bytes(nb::handle value, const char* field)
{
    if (!PyBytes_Check(value.ptr())) {
        const std::string message = std::string(field) + " must be bytes";
        throw nb::type_error(message.c_str());
    }
    char* data = nullptr;
    Py_ssize_t size = 0;
    if (PyBytes_AsStringAndSize(value.ptr(), &data, &size) != 0) {
        throw nb::python_error();
    }
    return std::string(data, static_cast<std::size_t>(size));
}

std::vector<std::string> operation_digest_bytes(const nb::tuple& values)
{
    std::vector<std::string> result;
    result.reserve(values.size());
    for (nb::handle value : values) {
        result.push_back(from_bytes(value, "authenticated operation digest"));
    }
    return result;
}

AuthenticatedOperationDigest2 authenticated_operation_digest(
    const nb::bytes& value)
{
    return AuthenticatedOperationDigestAuthority2::from_external_bytes(
        from_bytes(value, "authenticated operation digest"));
}

const char* verdict_text(AuditTeaVerdict2 verdict)
{
    switch (verdict) {
    case AuditTeaVerdict2::CERTIFIED:
        return "certified";
    case AuditTeaVerdict2::CAP_EXCEEDED:
        return "cap_exceeded";
    case AuditTeaVerdict2::UNRESOLVED:
        return "unresolved";
    }
    throw AuditResultEvidenceError("unknown audit result verdict");
}

const char* reason_text(AuditNonEngagingReason2 reason)
{
    switch (reason) {
    case AuditNonEngagingReason2::VERTICAL_RETRACT:
        return "vertical_retract";
    case AuditNonEngagingReason2::CLEARANCE_TRANSPORT:
        return "clearance_transport";
    }
    throw AuditResultEvidenceError("unknown non-engaging result reason");
}

} // namespace

void register_audit_replay_2(nb::module_& module)
{
    nb::exception<AuditDigestSizeError>(
        module, "AuditDigestSizeError", PyExc_ValueError);
    nb::exception<AuditReplayError> replay_error(
        module, "AuditReplayError", PyExc_RuntimeError);
    nb::exception<AuditReplayRequestIdentityError>(
        module, "AuditReplayRequestIdentityError", replay_error.ptr());
    nb::exception<AuditReplayCardinalityError>(
        module, "AuditReplayCardinalityError", replay_error.ptr());
    nb::exception<AuditReplayOperationIdentityError>(
        module, "AuditReplayOperationIdentityError", replay_error.ptr());
    nb::exception<AuditReplayMotionIdentityError>(
        module, "AuditReplayMotionIdentityError", replay_error.ptr());
    nb::exception<AuditReplayOperationExhaustedError>(
        module, "AuditReplayOperationExhaustedError", replay_error.ptr());
    nb::exception<AuditReplayIncompleteError>(
        module, "AuditReplayIncompleteError", replay_error.ptr());
    nb::exception<AuditReplayFinalizedError>(
        module, "AuditReplayFinalizedError", replay_error.ptr());
    nb::exception<AuditReplayDecisionEvidenceError>(
        module, "AuditReplayDecisionEvidenceError", replay_error.ptr());
    nb::exception<AuditReplayDepletionEvidenceError>(
        module, "AuditReplayDepletionEvidenceError", replay_error.ptr());
    nb::exception<AuditResultEvidenceError>(
        module, "AuditResultEvidenceError", PyExc_RuntimeError);
    nb::exception<AuditDepletionWitnessError>(
        module, "AuditDepletionWitnessError", PyExc_RuntimeError);

    nb::class_<AuditReplay2>(module, "AuditReplay2", nb::is_final());
    nb::class_<AuditLateralResult2>(
        module, "AuditLateralResult2", nb::is_final())
        .def_prop_ro("verdict", [](const AuditLateralResult2& result) {
            return verdict_text(result.verdict());
        })
        .def_prop_ro("evidence_count", &AuditLateralResult2::evidence_count)
        .def_prop_ro(
            "authenticated_operation_digest",
            [](const AuditLateralResult2& result) {
                return as_bytes(
                    result.authenticated_operation_digest().bytes());
            })
        .def_prop_ro("decision_digest", [](const AuditLateralResult2& result) {
            return as_bytes(result.decision_digest().bytes());
        })
        .def_prop_ro(
            "depletion_witness_digest",
            [](const AuditLateralResult2& result) {
                return as_bytes(result.depletion_witness_digest().bytes());
            })
        .def_prop_ro("pre_lineage", [](const AuditLateralResult2& result) {
            return as_bytes(result.pre_lineage().bytes());
        })
        .def_prop_ro("post_lineage", [](const AuditLateralResult2& result) {
            return as_bytes(result.post_lineage().bytes());
        })
        .def_prop_ro("digest", [](const AuditLateralResult2& result) {
            return as_bytes(result.digest().bytes());
        });
    nb::class_<AuditPlungeResult2>(
        module, "AuditPlungeResult2", nb::is_final())
        .def_prop_ro(
            "authenticated_operation_digest",
            [](const AuditPlungeResult2& result) {
                return as_bytes(
                    result.authenticated_operation_digest().bytes());
            })
        .def_prop_ro(
            "depletion_witness_digest",
            [](const AuditPlungeResult2& result) {
                return as_bytes(result.depletion_witness_digest().bytes());
            })
        .def_prop_ro("pre_lineage", [](const AuditPlungeResult2& result) {
            return as_bytes(result.pre_lineage().bytes());
        })
        .def_prop_ro("post_lineage", [](const AuditPlungeResult2& result) {
            return as_bytes(result.post_lineage().bytes());
        })
        .def_prop_ro("digest", [](const AuditPlungeResult2& result) {
            return as_bytes(result.digest().bytes());
        });
    nb::class_<AuditNonEngagingResult2>(
        module, "AuditNonEngagingResult2", nb::is_final())
        .def_prop_ro("reason", [](const AuditNonEngagingResult2& result) {
            return reason_text(result.reason());
        })
        .def_prop_ro(
            "authenticated_operation_digest",
            [](const AuditNonEngagingResult2& result) {
                return as_bytes(
                    result.authenticated_operation_digest().bytes());
            })
        .def_prop_ro("pre_lineage", [](const AuditNonEngagingResult2& result) {
            return as_bytes(result.pre_lineage().bytes());
        })
        .def_prop_ro("post_lineage", [](const AuditNonEngagingResult2& result) {
            return as_bytes(result.post_lineage().bytes());
        })
        .def_prop_ro("digest", [](const AuditNonEngagingResult2& result) {
            return as_bytes(result.digest().bytes());
        });
    nb::class_<AuditReplayCompletion2>(
        module, "AuditReplayCompletion2", nb::is_final())
        .def_prop_ro("request_digest", [](const AuditReplayCompletion2& result) {
            return as_bytes(result.request_digest().bytes());
        })
        .def_prop_ro("operation_count", &AuditReplayCompletion2::operation_count)
        .def_prop_ro(
            "terminal_lineage",
            [](const AuditReplayCompletion2& result) {
                return as_bytes(result.terminal_lineage().bytes());
            })
        .def_prop_ro("digest", [](const AuditReplayCompletion2& result) {
            return as_bytes(result.digest().bytes());
        });

    module.def(
        "begin_audit_replay",
        [](Eigen::Ref<const compas::RowMatrixXd> boundary,
           const std::vector<compas::RowMatrixXd>& holes,
           nb::bytes input_digest,
           const AuditNativeRequestIdentity2& native_request,
           const AuditPolicy2& policy,
           const AuditDecisionLimits2& decision_limits,
           const nb::tuple& operation_digests) {
            return begin_audit_replay(
                boundary,
                holes,
                from_bytes(input_digest, "audit input digest"),
                native_request,
                policy,
                decision_limits,
                operation_digest_bytes(operation_digests));
        },
        "boundary"_a,
        "holes"_a,
        "input_digest"_a,
        "native_request"_a,
        "policy"_a,
        "decision_limits"_a,
        "authenticated_operation_digests"_a);
    module.def(
        "audit_deplete_segment",
        [](AuditReplay2& replay,
           const AuditSegmentMotion2& motion,
           const nb::bytes& operation_digest) {
            return audit_deplete_segment(
                replay, motion, authenticated_operation_digest(operation_digest));
        },
        "replay"_a,
        "motion"_a,
        "authenticated_operation_digest"_a);
    module.def(
        "audit_deplete_circle",
        [](AuditReplay2& replay,
           const AuditCircleMotion2& motion,
           const nb::bytes& operation_digest) {
            return audit_deplete_circle(
                replay, motion, authenticated_operation_digest(operation_digest));
        },
        "replay"_a,
        "motion"_a,
        "authenticated_operation_digest"_a);
    module.def(
        "audit_deplete_arc",
        [](AuditReplay2& replay,
           const AuditArcMotion2& motion,
           const nb::bytes& operation_digest) {
            return audit_deplete_arc(
                replay, motion, authenticated_operation_digest(operation_digest));
        },
        "replay"_a,
        "motion"_a,
        "authenticated_operation_digest"_a);
    module.def(
        "deplete_audit_plunge",
        [](AuditReplay2& replay,
           const AuditVerticalPlunge2& motion,
           const nb::bytes& operation_digest) {
            return deplete_audit_plunge(
                replay, motion, authenticated_operation_digest(operation_digest));
        },
        "replay"_a,
        "motion"_a,
        "authenticated_operation_digest"_a);
    module.def(
        "record_audit_retract",
        [](AuditReplay2& replay,
           const AuditVerticalRetract2& motion,
           const nb::bytes& operation_digest) {
            return record_audit_retract(
                replay, motion, authenticated_operation_digest(operation_digest));
        },
        "replay"_a,
        "motion"_a,
        "authenticated_operation_digest"_a);
    module.def(
        "record_audit_clearance",
        [](AuditReplay2& replay,
           const AuditClearanceTransport2& motion,
           const nb::bytes& operation_digest) {
            return record_audit_clearance(
                replay, motion, authenticated_operation_digest(operation_digest));
        },
        "replay"_a,
        "motion"_a,
        "authenticated_operation_digest"_a);
    module.def("finish_audit_replay", &finish_audit_replay, "replay"_a);
}
