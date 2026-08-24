#pragma once

#include "audit_classification_core_2.h"
#include "audit_motion_result_2.h"
#include "compas_matrix.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

struct ReplayProgress2;
struct AuditReplayRequestState2;
struct AuditReplayStorage2;
enum class AuditReplayFailurePoint2;
struct AuditReplayInspection2;
struct AuditReplayInstrumentation2;

class AuditReplayError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class AuditReplayRequestIdentityError : public AuditReplayError {
public:
    using AuditReplayError::AuditReplayError;
};

class AuditReplayCardinalityError : public AuditReplayError {
public:
    using AuditReplayError::AuditReplayError;
};

class AuditReplayOperationIdentityError : public AuditReplayError {
public:
    using AuditReplayError::AuditReplayError;
};

class AuditReplayMotionIdentityError : public AuditReplayError {
public:
    using AuditReplayError::AuditReplayError;
};

class AuditReplayOperationExhaustedError : public AuditReplayError {
public:
    using AuditReplayError::AuditReplayError;
};

class AuditReplayIncompleteError : public AuditReplayError {
public:
    using AuditReplayError::AuditReplayError;
};

class AuditReplayFinalizedError : public AuditReplayError {
public:
    using AuditReplayError::AuditReplayError;
};

class AuditReplayInjectedFailure : public AuditReplayError {
public:
    using AuditReplayError::AuditReplayError;
};

class AuditReplayDecisionEvidenceError : public AuditReplayError {
public:
    using AuditReplayError::AuditReplayError;
};

class AuditReplayDepletionEvidenceError : public AuditReplayError {
public:
    using AuditReplayError::AuditReplayError;
};

class AuditReplay2 {
public:
    static AuditReplay2 build(
        Eigen::Ref<const compas::RowMatrixXd> boundary,
        const std::vector<compas::RowMatrixXd>& holes,
        const AuditInputDigest2& input_digest,
        const AuditNativeRequestIdentity2& native_request,
        const AuditPolicy2& policy,
        const AuditDecisionLimits2& decision_limits,
        std::vector<AuthenticatedOperationDigest2>
            authenticated_operation_digests);

    AuditReplay2(AuditReplay2&&) noexcept;
    AuditReplay2& operator=(AuditReplay2&&) noexcept;
    AuditReplay2(const AuditReplay2&) = delete;
    AuditReplay2& operator=(const AuditReplay2&) = delete;
    ~AuditReplay2();

private:
    AuditReplay2(
        AuditReplayRequestState2 request,
        Stock2 stock,
        ReplayProgress2 progress);
    AuditLateralResult2 deplete_segment(
        const AuditSegmentMotion2& motion,
        const AuthenticatedOperationDigest2& operation_digest);
    AuditLateralResult2 deplete_circle(
        const AuditCircleMotion2& motion,
        const AuthenticatedOperationDigest2& operation_digest);
    AuditLateralResult2 deplete_arc(
        const AuditArcMotion2& motion,
        const AuthenticatedOperationDigest2& operation_digest);
    AuditPlungeResult2 deplete_plunge(
        const AuditVerticalPlunge2& motion,
        const AuthenticatedOperationDigest2& operation_digest);
    AuditNonEngagingResult2 record_retract(
        const AuditVerticalRetract2& motion,
        const AuthenticatedOperationDigest2& operation_digest);
    AuditNonEngagingResult2 record_clearance(
        const AuditClearanceTransport2& motion,
        const AuthenticatedOperationDigest2& operation_digest);
    AuditReplayCompletion2 finish();
    template <class Motion, class Apply, class Validate>
    AuditLateralResult2 deplete_lateral(
        const Motion& motion,
        const AuthenticatedOperationDigest2& operation_digest,
        Apply&& apply,
        Validate&& validate);
    template <class Motion>
    AuditNonEngagingResult2 record_non_engaging(
        const Motion& motion,
        const AuthenticatedOperationDigest2& operation_digest,
        AuditNonEngagingReason2 reason);

    std::unique_ptr<AuditReplayStorage2> storage_;

    friend AuditLateralResult2 audit_deplete_segment(
        AuditReplay2&,
        const AuditSegmentMotion2&,
        const AuthenticatedOperationDigest2&);
    friend AuditLateralResult2 audit_deplete_circle(
        AuditReplay2&,
        const AuditCircleMotion2&,
        const AuthenticatedOperationDigest2&);
    friend AuditLateralResult2 audit_deplete_arc(
        AuditReplay2&,
        const AuditArcMotion2&,
        const AuthenticatedOperationDigest2&);
    friend AuditPlungeResult2 deplete_audit_plunge(
        AuditReplay2&,
        const AuditVerticalPlunge2&,
        const AuthenticatedOperationDigest2&);
    friend AuditNonEngagingResult2 record_audit_retract(
        AuditReplay2&,
        const AuditVerticalRetract2&,
        const AuthenticatedOperationDigest2&);
    friend AuditNonEngagingResult2 record_audit_clearance(
        AuditReplay2&,
        const AuditClearanceTransport2&,
        const AuthenticatedOperationDigest2&);
    friend AuditReplayCompletion2 finish_audit_replay(AuditReplay2&);
    friend void inject_audit_replay_failure_for_test(
        AuditReplay2&,
        AuditReplayFailurePoint2) noexcept;
    friend void clear_audit_replay_failure_for_test(AuditReplay2&) noexcept;
    friend AuditReplayInspection2 inspect_audit_replay_for_test(
        const AuditReplay2&);
    friend AuditReplayInstrumentation2 audit_replay_instrumentation_for_test(
        const AuditReplay2&) noexcept;
};

AuditLateralResult2 audit_deplete_segment(
    AuditReplay2& replay,
    const AuditSegmentMotion2& motion,
    const AuthenticatedOperationDigest2& authenticated_operation_digest);
AuditLateralResult2 audit_deplete_circle(
    AuditReplay2& replay,
    const AuditCircleMotion2& motion,
    const AuthenticatedOperationDigest2& authenticated_operation_digest);
AuditLateralResult2 audit_deplete_arc(
    AuditReplay2& replay,
    const AuditArcMotion2& motion,
    const AuthenticatedOperationDigest2& authenticated_operation_digest);
AuditPlungeResult2 deplete_audit_plunge(
    AuditReplay2& replay,
    const AuditVerticalPlunge2& motion,
    const AuthenticatedOperationDigest2& authenticated_operation_digest);
AuditNonEngagingResult2 record_audit_retract(
    AuditReplay2& replay,
    const AuditVerticalRetract2& motion,
    const AuthenticatedOperationDigest2& authenticated_operation_digest);
AuditNonEngagingResult2 record_audit_clearance(
    AuditReplay2& replay,
    const AuditClearanceTransport2& motion,
    const AuthenticatedOperationDigest2& authenticated_operation_digest);
AuditReplayCompletion2 finish_audit_replay(AuditReplay2& replay);

AuditReplay2 begin_audit_replay(
    Eigen::Ref<const compas::RowMatrixXd> boundary,
    const std::vector<compas::RowMatrixXd>& holes,
    std::string input_digest_bytes,
    const AuditNativeRequestIdentity2& native_request,
    const AuditPolicy2& policy,
    const AuditDecisionLimits2& decision_limits,
    std::vector<std::string> authenticated_operation_digest_bytes);
