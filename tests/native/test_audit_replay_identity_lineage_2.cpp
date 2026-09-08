#include "audit_replay_fixtures_2.h"

#include "audit_depletion_witness_2.h"
#include "audit_lineage_2.h"
#include "audit_motion_result_2.h"
#include "canonical_encoding.h"
#include "stock_exact_depletion_2.h"

#include <set>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

using namespace audit_replay_fixtures;

namespace {

template <class Type>
concept HasPublicBuild = requires { &Type::build; };

template <class Type>
concept HasPublicCreate = requires { &Type::create; };

template <class Type>
concept HasPublicDerive = requires { &Type::derive; };

template <class Type>
concept HasPublicMint = requires { &Type::mint; };

template <class Authority>
concept HasPublicCanonicalHash = requires(std::string_view canonical) {
    Authority::hash_canonical(canonical);
};

static_assert(!HasPublicBuild<AuditMotionResult2>);
static_assert(!HasPublicCreate<AuditMotionResult2>);
static_assert(!HasPublicMint<AuditMotionResult2>);
static_assert(!HasPublicBuild<AuditLineage2>);
static_assert(!HasPublicCreate<AuditLineage2>);
static_assert(!HasPublicDerive<AuditLineage2>);
static_assert(!HasPublicMint<AuditLineage2>);
static_assert(!HasPublicBuild<AuditLateralResult2>);
static_assert(!HasPublicBuild<AuditPlungeResult2>);
static_assert(!HasPublicBuild<AuditNonEngagingResult2>);
static_assert(!HasPublicBuild<AuditReplayCompletion2>);
static_assert(!HasPublicCanonicalHash<ExactDepletionTraceDigestAuthority2>);
static_assert(!HasPublicCanonicalHash<DepletionWitnessDigestAuthority2>);
static_assert(!HasPublicCanonicalHash<StockLineageDigestAuthority2>);
static_assert(!HasPublicCanonicalHash<AuditResultDigestAuthority2>);
static_assert(!std::is_aggregate_v<AuditDepletionWitness2>);
static_assert(!std::is_aggregate_v<AuditLineage2>);
static_assert(!std::is_aggregate_v<AuditLateralResult2>);
static_assert(!std::is_aggregate_v<AuditPlungeResult2>);
static_assert(!std::is_aggregate_v<AuditNonEngagingResult2>);
static_assert(!std::is_aggregate_v<AuditReplayCompletion2>);
static_assert(!std::is_default_constructible_v<AuditReplay2>);
static_assert(!std::is_copy_constructible_v<AuditReplay2>);
static_assert(!std::is_default_constructible_v<AuditDepletionWitness2>);
static_assert(!std::is_default_constructible_v<AuditLineage2>);
static_assert(!std::is_default_constructible_v<AuditLateralResult2>);
static_assert(!std::is_default_constructible_v<AuditPlungeResult2>);
static_assert(!std::is_default_constructible_v<AuditNonEngagingResult2>);
static_assert(!std::is_default_constructible_v<AuditReplayCompletion2>);
static_assert(!std::is_constructible_v<ExactDepletionTraceDigest2, std::string>);
static_assert(!std::is_constructible_v<DepletionWitnessDigest2, std::string>);
static_assert(!std::is_constructible_v<StockLineageDigest2, std::string>);
static_assert(!std::is_constructible_v<AuditResultDigest2, std::string>);
static_assert(!std::is_constructible_v<
              AuditLateralResult2,
              AuditTeaVerdict2,
              std::size_t,
              NativeDecisionDigest2,
              DepletionWitnessDigest2,
              StockLineageDigest2,
              StockLineageDigest2>);
static_assert(!std::is_same_v<
              ExactDepletionTraceDigest2,
              DepletionWitnessDigest2>);
static_assert(std::is_nothrow_move_constructible_v<AuditLateralResult2>);
static_assert(std::is_nothrow_move_constructible_v<AuditPlungeResult2>);
static_assert(std::is_nothrow_move_constructible_v<AuditNonEngagingResult2>);
static_assert(std::is_nothrow_move_constructible_v<AuditReplayCompletion2>);
static_assert(noexcept(std::declval<Stock2&>().swap(std::declval<Stock2&>())));
static_assert(std::is_nothrow_swappable_v<ReplayProgress2>);

template <class Error, class Function>
void preflight_rejects(Function&& function, const char* message)
{
    bool threw = false;
    try {
        std::forward<Function>(function)();
    } catch (const Error&) {
        threw = true;
    }
    require(threw, message);
    require(
        audit_replay_construction_instrumentation_for_test().stock_clone_count
            == 0,
        "replay preflight rejection reached stock cloning");
}

void limits_are_request_and_decision_identity_gate()
{
    // Kills hidden replay defaults and decisions made under foreign limits.
    const compas::RowMatrixXd boundary = rectangle(-5.0, -5.0, 5.0, 5.0);
    const AuditPolicy2 audit_policy = policy();
    const AuditSegmentMotion2 motion = segment();
    const AuthenticatedOperationDigest2 operation =
        operation_digest("limits");
    const AuditDecisionLimits2 baseline = decision_limits();
    const std::vector<AuditDecisionLimits2> changed{
        decision_limits(Epeck::FT(1) / Epeck::FT(8192), 8, 256),
        decision_limits(Epeck::FT(1) / Epeck::FT(4096), 7, 256),
        decision_limits(Epeck::FT(1) / Epeck::FT(4096), 8, 255),
    };
    const AuditNativeRequestIdentity2 request = one_motion_request(
        boundary, audit_policy, baseline, motion);
    for (const AuditDecisionLimits2& foreign : changed) {
        require(
            request.digest().bytes()
                != one_motion_request(
                       boundary, audit_policy, foreign, motion)
                       .digest()
                       .bytes(),
            "native request identity ignores one exact decision limit");
        reset_audit_replay_construction_instrumentation_for_test();
        preflight_rejects<AuditReplayRequestIdentityError>(
            [&] {
                static_cast<void>(AuditReplay2::build(
                    boundary,
                    {},
                    input_digest(),
                    request,
                    audit_policy,
                    foreign,
                    {operation}));
            },
            "replay accepted limits foreign to native request");
    }

    Stock2 direct_stock(boundary, {});
    const AuditDecisionWitness2 direct = certify_audit_tea_exact(
        direct_stock, motion, audit_policy, baseline);
    AuditReplay2 replay = one_motion_replay(
        boundary, audit_policy, baseline, motion, operation);
    const AuditLateralResult2 result =
        audit_deplete_segment(replay, motion, operation);
    require(
        result.decision_digest().bytes() == direct.digest().bytes(),
        "replay did not retain the direct decision under bound limits");
    require(
        result.evidence_count() == exact_evidence_count(direct.counters())
            && result.evidence_count() > 0,
        "evidence count is not exact station plus coverage replay count");
}

void request_preflight_gate()
{
    // Kills malformed ingress, foreign policy/rings, and cardinality after clone.
    const compas::RowMatrixXd boundary = rectangle(-5.0, -5.0, 5.0, 5.0);
    const compas::RowMatrixXd hole = square_hole();
    const AuditPolicy2 audit_policy = policy();
    const AuditDecisionLimits2 limits = decision_limits();
    const AuditSegmentMotion2 motion = segment();
    const AuthenticatedOperationDigest2 operation =
        operation_digest("preflight");
    const AuditNativeRequestIdentity2 request =
        AuditNativeRequestIdentity2::build(
            AuditNativeStockIdentity2::build(boundary, {hole}),
            audit_policy,
            limits,
            {motion.digest()});

    reset_audit_replay_construction_instrumentation_for_test();
    preflight_rejects<AuditDigestSizeError>(
        [&] {
            static_cast<void>(begin_audit_replay(
                boundary,
                {hole},
                "short",
                request,
                audit_policy,
                limits,
                {operation.bytes()}));
        },
        "short external input digest was accepted");
    preflight_rejects<AuditDigestSizeError>(
        [&] {
            static_cast<void>(begin_audit_replay(
                boundary,
                {hole},
                input_digest().bytes(),
                request,
                audit_policy,
                limits,
                {"short"}));
        },
        "short authenticated-operation digest was accepted");
    preflight_rejects<AuditReplayCardinalityError>(
        [&] {
            static_cast<void>(AuditReplay2::build(
                boundary,
                {hole},
                input_digest(),
                request,
                audit_policy,
                limits,
                {}));
        },
        "missing authenticated-operation digest was accepted");
    preflight_rejects<AuditReplayCardinalityError>(
        [&] {
            static_cast<void>(AuditReplay2::build(
                boundary,
                {hole},
                input_digest(),
                request,
                audit_policy,
                limits,
                {operation, operation_digest("extra")}));
        },
        "surplus authenticated-operation digest was accepted");
    preflight_rejects<AuditReplayRequestIdentityError>(
        [&] {
            static_cast<void>(AuditReplay2::build(
                boundary,
                {hole},
                input_digest(),
                request,
                policy(2048),
                limits,
                {operation}));
        },
        "foreign audit policy was accepted");
    preflight_rejects<AuditReplayRequestIdentityError>(
        [&] {
            static_cast<void>(AuditReplay2::build(
                rectangle(-6.0, -5.0, 5.0, 5.0),
                {hole},
                input_digest(),
                request,
                audit_policy,
                limits,
                {operation}));
        },
        "changed stock boundary was accepted");
    preflight_rejects<AuditReplayRequestIdentityError>(
        [&] {
            static_cast<void>(AuditReplay2::build(
                boundary,
                {square_hole(-3.0, -3.0, -1.5, -2.0)},
                input_digest(),
                request,
                audit_policy,
                limits,
                {operation}));
        },
        "changed stock hole was accepted");

    AuditReplay2 equivalent = AuditReplay2::build(
        reversed_rotated_ring(boundary),
        {reversed_rotated_ring(hole)},
        input_digest(),
        request,
        audit_policy,
        limits,
        {operation});
    require(
        inspect_audit_replay_for_test(equivalent).cursor == 0,
        "canonical ring rotation/orientation did not preserve request identity");
    require(
        audit_replay_construction_instrumentation_for_test().stock_clone_count
            == 0,
        "canonical request construction cloned stock");
}

void lineage_component_mutation_gate()
{
    // Kills retagging or omission of any seed/transition component.
    const compas::RowMatrixXd boundary = rectangle(-5.0, -5.0, 5.0, 5.0);
    const AuditPolicy2 audit_policy = policy();
    const AuditDecisionLimits2 limits = decision_limits();
    const AuditSegmentMotion2 motion = segment();
    const AuditSegmentMotion2 changed_motion = segment(1.0);
    const AuditNativeRequestIdentity2 request = one_motion_request(
        boundary, audit_policy, limits, motion);
    const AuditNativeRequestIdentity2 changed_request = one_motion_request(
        boundary, audit_policy, limits, changed_motion);
    const AuditInputDigest2 input = input_digest();
    const AuditLineage2 seed = AuditReplayTestAuthority2::seed_lineage(
        input, request);
    const AuditLineage2 changed_input_seed =
        AuditReplayTestAuthority2::seed_lineage(
            input_digest("changed"), request);
    const AuditLineage2 changed_request_seed =
        AuditReplayTestAuthority2::seed_lineage(input, changed_request);
    require(
        seed.digest().bytes() != input.bytes()
            && seed.digest().bytes() != request.digest().bytes()
            && seed.digest().bytes() != changed_input_seed.digest().bytes()
            && seed.digest().bytes() != changed_request_seed.digest().bytes(),
        "lineage seed retags or omits one typed root");

    Stock2 authority(boundary, {});
    const AuditDecisionWitness2 decision = certify_audit_tea_exact(
        authority, motion, audit_policy, limits);
    Stock2 trial = authority.clone();
    const AuditDepletionWitness2 depletion =
        apply_audit_segment_depletion_to_trial(
            authority, trial, motion, audit_policy);
    Stock2 changed_authority(boundary, {});
    const AuditDecisionWitness2 changed_decision = certify_audit_tea_exact(
        changed_authority, changed_motion, audit_policy, limits);
    Stock2 changed_trial = changed_authority.clone();
    const AuditDepletionWitness2 changed_depletion =
        apply_audit_segment_depletion_to_trial(
            changed_authority,
            changed_trial,
            changed_motion,
            audit_policy);
    const AuthenticatedOperationDigest2 operation =
        operation_digest("lineage");
    const AuditLineage2 baseline =
        AuditReplayTestAuthority2::transition_lineage(
            request, 0, seed, operation, motion.digest(), decision, depletion);
    const std::vector<AuditLineage2> mutations{
        AuditReplayTestAuthority2::transition_lineage(
            changed_request, 0, seed, operation, motion.digest(), decision, depletion),
        AuditReplayTestAuthority2::transition_lineage(
            request, 1, seed, operation, motion.digest(), decision, depletion),
        AuditReplayTestAuthority2::transition_lineage(
            request, 0, changed_input_seed, operation, motion.digest(), decision, depletion),
        AuditReplayTestAuthority2::transition_lineage(
            request, 0, seed, operation_digest("changed"), motion.digest(), decision, depletion),
        AuditReplayTestAuthority2::transition_lineage(
            request, 0, seed, operation, changed_motion.digest(), decision, depletion),
        AuditReplayTestAuthority2::transition_lineage(
            request, 0, seed, operation, motion.digest(), changed_decision, depletion),
        AuditReplayTestAuthority2::transition_lineage(
            request, 0, seed, operation, motion.digest(), decision, changed_depletion),
    };
    std::set<std::string> digests{baseline.digest().bytes()};
    for (const AuditLineage2& mutation : mutations) {
        digests.insert(mutation.digest().bytes());
    }
    require(
        digests.size() == mutations.size() + 1,
        "lineage transition aliases a typed transition component");

    const AuditVerticalPlunge2 plunge_motion = plunge();
    const AuditVerticalPlunge2 changed_plunge = std::get<AuditVerticalPlunge2>(
        classify_audit_line(
            {2.0, 3.0, 5.0},
            {2.0, 3.0, 0.0},
            0.0,
            5.0,
            "plunge"));
    const AuditNativeRequestIdentity2 plunge_request = one_motion_request(
        boundary, audit_policy, limits, plunge_motion);
    Stock2 plunge_authority(boundary, {});
    Stock2 plunge_trial = plunge_authority.clone();
    const AuditDepletionWitness2 plunge_depletion =
        apply_audit_plunge_depletion_to_trial(
            plunge_authority,
            plunge_trial,
            plunge_motion,
            audit_policy);
    Stock2 changed_plunge_authority(boundary, {});
    Stock2 changed_plunge_trial = changed_plunge_authority.clone();
    const AuditDepletionWitness2 changed_plunge_depletion =
        apply_audit_plunge_depletion_to_trial(
            changed_plunge_authority,
            changed_plunge_trial,
            changed_plunge,
            audit_policy);
    const AuditLineage2 plunge_baseline =
        AuditReplayTestAuthority2::transition_plunge_lineage(
            plunge_request,
            0,
            seed,
            operation,
            plunge_motion.digest(),
            plunge_depletion);
    const std::vector<AuditLineage2> plunge_mutations{
        AuditReplayTestAuthority2::transition_plunge_lineage(
            changed_request,
            0,
            seed,
            operation,
            plunge_motion.digest(),
            plunge_depletion),
        AuditReplayTestAuthority2::transition_plunge_lineage(
            plunge_request,
            1,
            seed,
            operation,
            plunge_motion.digest(),
            plunge_depletion),
        AuditReplayTestAuthority2::transition_plunge_lineage(
            plunge_request,
            0,
            changed_input_seed,
            operation,
            plunge_motion.digest(),
            plunge_depletion),
        AuditReplayTestAuthority2::transition_plunge_lineage(
            plunge_request,
            0,
            seed,
            operation_digest("changed-plunge"),
            plunge_motion.digest(),
            plunge_depletion),
        AuditReplayTestAuthority2::transition_plunge_lineage(
            plunge_request,
            0,
            seed,
            operation,
            changed_plunge.digest(),
            plunge_depletion),
        AuditReplayTestAuthority2::transition_plunge_lineage(
            plunge_request,
            0,
            seed,
            operation,
            plunge_motion.digest(),
            changed_plunge_depletion),
    };
    std::set<std::string> plunge_digests{plunge_baseline.digest().bytes()};
    for (const AuditLineage2& mutation : plunge_mutations) {
        plunge_digests.insert(mutation.digest().bytes());
    }
    require(
        plunge_digests.size() == plunge_mutations.size() + 1,
        "plunge lineage aliases a typed transition component");
}

void no_material_and_nonengaging_lineage_gate()
{
    // Kills stock-change-derived chronology and non-engaging lineage drift.
    const AuditPolicy2 audit_policy = policy();
    const AuditDecisionLimits2 limits = decision_limits();
    const AuditSegmentMotion2 motion = segment();
    const AuthenticatedOperationDigest2 operation =
        operation_digest("no-material");
    const compas::RowMatrixXd disjoint = rectangle(10.0, 10.0, 12.0, 12.0);
    AuditReplay2 no_material = one_motion_replay(
        disjoint, audit_policy, limits, motion, operation);
    const AuditReplayInspection2 before =
        inspect_audit_replay_for_test(no_material);
    const AuditLateralResult2 no_material_result =
        audit_deplete_segment(no_material, motion, operation);
    const AuditReplayInspection2 after =
        inspect_audit_replay_for_test(no_material);
    require(
        before.stock_digest == after.stock_digest
            && no_material_result.pre_lineage().bytes()
                != no_material_result.post_lineage().bytes(),
        "no-material mutating operation did not advance lineage");

    const compas::RowMatrixXd boundary = rectangle(-5.0, -5.0, 5.0, 5.0);
    const AuditVerticalRetract2 retract_motion = retract();
    const AuditClearanceTransport2 clearance_motion = clearance();
    const AuditSegmentMotion2 second_motion = segment(1.0);
    const std::vector<AuthenticatedOperationDigest2> operations{
        operation_digest("first"),
        operation_digest("retract"),
        operation_digest("clearance"),
        operation_digest("second"),
    };
    const AuditNativeRequestIdentity2 request =
        AuditNativeRequestIdentity2::build(
            AuditNativeStockIdentity2::build(boundary, {}),
            audit_policy,
            limits,
            {
                motion.digest(),
                retract_motion.digest(),
                clearance_motion.digest(),
                second_motion.digest(),
            });
    AuditReplay2 replay = AuditReplay2::build(
        boundary,
        {},
        input_digest(),
        request,
        audit_policy,
        limits,
        operations);
    const AuditLateralResult2 first =
        audit_deplete_segment(replay, motion, operations[0]);
    const AuditNonEngagingResult2 retract_result =
        record_audit_retract(replay, retract_motion, operations[1]);
    const AuditNonEngagingResult2 clearance_result =
        record_audit_clearance(replay, clearance_motion, operations[2]);
    const AuditLateralResult2 second =
        audit_deplete_segment(replay, second_motion, operations[3]);
    require(
        first.post_lineage().bytes() == retract_result.pre_lineage().bytes()
            && retract_result.pre_lineage().bytes()
                == retract_result.post_lineage().bytes()
            && retract_result.post_lineage().bytes()
                == clearance_result.pre_lineage().bytes()
            && clearance_result.pre_lineage().bytes()
                == clearance_result.post_lineage().bytes()
            && clearance_result.post_lineage().bytes()
                == second.pre_lineage().bytes(),
        "non-engaging chain changed or disconnected lineage");
    require(
        retract_result.authenticated_operation_digest().bytes()
                    == operations[1].bytes()
            && clearance_result.authenticated_operation_digest().bytes()
                == operations[2].bytes()
            && first.authenticated_operation_digest().bytes()
                == operations[0].bytes()
            && second.authenticated_operation_digest().bytes()
                == operations[3].bytes()
            &&
        retract_result.reason()
                == AuditNonEngagingReason2::VERTICAL_RETRACT
            && clearance_result.reason()
                == AuditNonEngagingReason2::CLEARANCE_TRANSPORT
            && retract_result.digest().bytes()
                != clearance_result.digest().bytes(),
        "retract and clearance reasons do not produce distinct identities");
}

void result_and_completion_component_gate()
{
    // Kills result/completion hashes that omit a typed constituent.
    const compas::RowMatrixXd boundary = rectangle(-5.0, -5.0, 5.0, 5.0);
    const AuditPolicy2 audit_policy = policy();
    const AuditDecisionLimits2 limits = decision_limits();
    const AuditSegmentMotion2 motion = segment();
    const AuditSegmentMotion2 changed_motion = segment(1.0);
    const AuthenticatedOperationDigest2 operation =
        operation_digest("result");
    const AuditNativeRequestIdentity2 request = one_motion_request(
        boundary, audit_policy, limits, motion);
    const AuditNativeRequestIdentity2 changed_request = one_motion_request(
        boundary, audit_policy, limits, changed_motion);

    Stock2 authority(boundary, {});
    const AuditDecisionWitness2 decision = certify_audit_tea_exact(
        authority, motion, audit_policy, limits);
    Stock2 trial = authority.clone();
    const AuditDepletionWitness2 depletion =
        apply_audit_segment_depletion_to_trial(
            authority, trial, motion, audit_policy);
    Stock2 changed_authority(boundary, {});
    const AuditDecisionWitness2 changed_decision = certify_audit_tea_exact(
        changed_authority, changed_motion, audit_policy, limits);
    Stock2 changed_trial = changed_authority.clone();
    const AuditDepletionWitness2 changed_depletion =
        apply_audit_segment_depletion_to_trial(
            changed_authority,
            changed_trial,
            changed_motion,
            audit_policy);

    AuditReplay2 replay = AuditReplay2::build(
        boundary,
        {},
        input_digest(),
        request,
        audit_policy,
        limits,
        {operation});
    const AuditLateralResult2 result =
        audit_deplete_segment(replay, motion, operation);
    require(
        result.request_digest().bytes() == request.digest().bytes()
            && result.cursor() == 0
            && result.motion_digest().bytes() == motion.digest().bytes()
            && result.authenticated_operation_digest().bytes() == operation.bytes()
            && result.decision_digest().bytes() == decision.digest().bytes()
            && result.depletion_witness_digest().bytes()
                == depletion.digest().bytes(),
        "result does not retain independently reproduced witnesses");
    const AuditResultDigest2 expected =
        AuditReplayTestAuthority2::lateral_result_digest(
            request,
            0,
            operation,
            motion.digest(),
            decision,
            depletion,
            result.pre_lineage(),
            result.post_lineage());
    require(
        expected.bytes() == result.digest().bytes(),
        "lateral result diverges from typed constituent authority");
    const AuditLineage2 foreign_lineage =
        AuditReplayTestAuthority2::seed_lineage(
            input_digest("foreign-lineage"), request);
    const std::vector<AuditResultDigest2> mutations{
        AuditReplayTestAuthority2::lateral_result_digest(
            changed_request, 0, operation, motion.digest(), decision, depletion,
            result.pre_lineage(), result.post_lineage()),
        AuditReplayTestAuthority2::lateral_result_digest(
            request, 1, operation, motion.digest(), decision, depletion,
            result.pre_lineage(), result.post_lineage()),
        AuditReplayTestAuthority2::lateral_result_digest(
            request, 0, operation_digest("result-changed"), motion.digest(),
            decision, depletion, result.pre_lineage(), result.post_lineage()),
        AuditReplayTestAuthority2::lateral_result_digest(
            request, 0, operation, changed_motion.digest(), decision, depletion,
            result.pre_lineage(), result.post_lineage()),
        AuditReplayTestAuthority2::lateral_result_digest(
            request, 0, operation, motion.digest(), changed_decision, depletion,
            result.pre_lineage(), result.post_lineage()),
        AuditReplayTestAuthority2::lateral_result_digest(
            request, 0, operation, motion.digest(), decision, changed_depletion,
            result.pre_lineage(), result.post_lineage()),
        AuditReplayTestAuthority2::lateral_result_digest(
            request, 0, operation, motion.digest(), decision, depletion,
            foreign_lineage.digest(), result.post_lineage()),
        AuditReplayTestAuthority2::lateral_result_digest(
            request, 0, operation, motion.digest(), decision, depletion,
            result.pre_lineage(), foreign_lineage.digest()),
    };
    for (const AuditResultDigest2& mutation : mutations) {
        require(
            mutation.bytes() != result.digest().bytes(),
            "lateral result digest ignores one typed component");
    }

    const AuditVerticalPlunge2 plunge_motion = plunge();
    const AuditVerticalPlunge2 changed_plunge = std::get<AuditVerticalPlunge2>(
        classify_audit_line(
            {2.0, 3.0, 5.0},
            {2.0, 3.0, 0.0},
            0.0,
            5.0,
            "plunge"));
    const AuditNativeRequestIdentity2 plunge_request = one_motion_request(
        boundary, audit_policy, limits, plunge_motion);
    const AuditNativeRequestIdentity2 changed_plunge_request =
        one_motion_request(boundary, audit_policy, limits, changed_plunge);
    Stock2 plunge_authority(boundary, {});
    Stock2 plunge_trial = plunge_authority.clone();
    const AuditDepletionWitness2 plunge_depletion =
        apply_audit_plunge_depletion_to_trial(
            plunge_authority,
            plunge_trial,
            plunge_motion,
            audit_policy);
    Stock2 changed_plunge_authority(boundary, {});
    Stock2 changed_plunge_trial = changed_plunge_authority.clone();
    const AuditDepletionWitness2 changed_plunge_depletion =
        apply_audit_plunge_depletion_to_trial(
            changed_plunge_authority,
            changed_plunge_trial,
            changed_plunge,
            audit_policy);
    const AuthenticatedOperationDigest2 plunge_operation =
        operation_digest("plunge-result");
    AuditReplay2 plunge_replay = one_motion_replay(
        boundary,
        audit_policy,
        limits,
        plunge_motion,
        plunge_operation);
    const AuditPlungeResult2 plunge_result = deplete_audit_plunge(
        plunge_replay, plunge_motion, plunge_operation);
    require(
        plunge_result.request_digest().bytes()
                == plunge_request.digest().bytes()
            && plunge_result.cursor() == 0
            && plunge_result.motion_digest().bytes()
                == plunge_motion.digest().bytes()
            && plunge_result.authenticated_operation_digest().bytes()
                == plunge_operation.bytes()
            && plunge_result.depletion_witness_digest().bytes()
                == plunge_depletion.digest().bytes(),
        "plunge result does not retain typed operation and depletion evidence");
    const AuditResultDigest2 expected_plunge =
        AuditReplayTestAuthority2::plunge_result_digest(
            plunge_request,
            0,
            plunge_operation,
            plunge_motion.digest(),
            plunge_depletion,
            plunge_result.pre_lineage(),
            plunge_result.post_lineage());
    require(
        expected_plunge.bytes() == plunge_result.digest().bytes(),
        "plunge result diverges from typed constituent authority");
    const std::vector<AuditResultDigest2> plunge_result_mutations{
        AuditReplayTestAuthority2::plunge_result_digest(
            changed_plunge_request, 0, plunge_operation, plunge_motion.digest(),
            plunge_depletion, plunge_result.pre_lineage(), plunge_result.post_lineage()),
        AuditReplayTestAuthority2::plunge_result_digest(
            plunge_request, 1, plunge_operation, plunge_motion.digest(),
            plunge_depletion, plunge_result.pre_lineage(), plunge_result.post_lineage()),
        AuditReplayTestAuthority2::plunge_result_digest(
            plunge_request, 0, operation_digest("changed-plunge-result"),
            plunge_motion.digest(), plunge_depletion, plunge_result.pre_lineage(),
            plunge_result.post_lineage()),
        AuditReplayTestAuthority2::plunge_result_digest(
            plunge_request, 0, plunge_operation, changed_plunge.digest(),
            plunge_depletion, plunge_result.pre_lineage(), plunge_result.post_lineage()),
        AuditReplayTestAuthority2::plunge_result_digest(
            plunge_request, 0, plunge_operation, plunge_motion.digest(),
            changed_plunge_depletion, plunge_result.pre_lineage(), plunge_result.post_lineage()),
        AuditReplayTestAuthority2::plunge_result_digest(
            plunge_request, 0, plunge_operation, plunge_motion.digest(),
            plunge_depletion, foreign_lineage.digest(), plunge_result.post_lineage()),
        AuditReplayTestAuthority2::plunge_result_digest(
            plunge_request, 0, plunge_operation, plunge_motion.digest(),
            plunge_depletion, plunge_result.pre_lineage(), foreign_lineage.digest()),
    };
    for (const AuditResultDigest2& mutation : plunge_result_mutations) {
        require(
            mutation.bytes() != plunge_result.digest().bytes(),
            "plunge result digest ignores one typed component");
    }

    const AuditVerticalRetract2 retract_motion = retract();
    const AuditClearanceTransport2 clearance_motion = clearance();
    const AuthenticatedOperationDigest2 retract_operation =
        operation_digest("nonengaging-result");
    const AuditNativeRequestIdentity2 retract_request = one_motion_request(
        boundary, audit_policy, limits, retract_motion);
    const AuditNativeRequestIdentity2 clearance_request = one_motion_request(
        boundary, audit_policy, limits, clearance_motion);
    AuditReplay2 retract_replay = one_motion_replay(
        boundary,
        audit_policy,
        limits,
        retract_motion,
        retract_operation);
    const AuditNonEngagingResult2 retract_result = record_audit_retract(
        retract_replay, retract_motion, retract_operation);
    require(
        retract_result.request_digest().bytes()
                == retract_request.digest().bytes()
            && retract_result.cursor() == 0
            && retract_result.motion_digest().bytes()
                == retract_motion.digest().bytes()
            && retract_result.authenticated_operation_digest().bytes()
            == retract_operation.bytes(),
        "non-engaging result does not retain authenticated operation digest");
    const AuditResultDigest2 expected_retract =
        AuditReplayTestAuthority2::non_engaging_result_digest(
            retract_request,
            0,
            retract_operation,
            retract_motion.digest(),
            AuditNonEngagingReason2::VERTICAL_RETRACT,
            retract_result.pre_lineage());
    require(
        expected_retract.bytes() == retract_result.digest().bytes(),
        "non-engaging result diverges from typed constituent authority");
    const std::vector<AuditResultDigest2> non_engaging_mutations{
        AuditReplayTestAuthority2::non_engaging_result_digest(
            clearance_request, 0, retract_operation, retract_motion.digest(),
            AuditNonEngagingReason2::VERTICAL_RETRACT, retract_result.pre_lineage()),
        AuditReplayTestAuthority2::non_engaging_result_digest(
            retract_request, 1, retract_operation, retract_motion.digest(),
            AuditNonEngagingReason2::VERTICAL_RETRACT, retract_result.pre_lineage()),
        AuditReplayTestAuthority2::non_engaging_result_digest(
            retract_request, 0, operation_digest("changed-nonengaging-result"),
            retract_motion.digest(), AuditNonEngagingReason2::VERTICAL_RETRACT,
            retract_result.pre_lineage()),
        AuditReplayTestAuthority2::non_engaging_result_digest(
            retract_request, 0, retract_operation, clearance_motion.digest(),
            AuditNonEngagingReason2::VERTICAL_RETRACT, retract_result.pre_lineage()),
        AuditReplayTestAuthority2::non_engaging_result_digest(
            retract_request, 0, retract_operation, retract_motion.digest(),
            AuditNonEngagingReason2::CLEARANCE_TRANSPORT,
            retract_result.pre_lineage()),
        AuditReplayTestAuthority2::non_engaging_result_digest(
            retract_request, 0, retract_operation, retract_motion.digest(),
            AuditNonEngagingReason2::VERTICAL_RETRACT, foreign_lineage.digest()),
    };
    for (const AuditResultDigest2& mutation : non_engaging_mutations) {
        require(
            mutation.bytes() != retract_result.digest().bytes(),
            "non-engaging result digest ignores one typed component");
    }

    const AuditReplayCompletion2 completion = finish_audit_replay(replay);
    const AuditInputDigest2 replay_input = input_digest();
    const AuditLineage2 seed =
        AuditReplayTestAuthority2::seed_lineage(replay_input, request);
    const AuditInputDigest2 changed_input = input_digest("changed-completion");
    const AuditLineage2 changed_seed =
        AuditReplayTestAuthority2::seed_lineage(changed_input, request);
    const AuditResultDigest2 expected_completion =
        AuditReplayTestAuthority2::completion_digest(
            replay_input, seed.digest(), request, 1, result.post_lineage());
    require(
        completion.input_digest().bytes() == replay_input.bytes()
            && completion.seed_lineage().bytes() == seed.digest().bytes()
            && completion.digest().bytes() == expected_completion.bytes(),
        "completion diverges from typed constituent authority");
    const std::string completion_v2_bytes = canonical_encode_tagged_union(
        "audit-replay-completion-v2",
        canonical_encode_component_map({
            {"audit-input-digest", replay_input.bytes()},
            {"native-request-digest", request.digest().bytes()},
            {"operation-count", canonical_audit_rational_bytes(Epeck::FT(1))},
            {"seed-lineage-digest", seed.digest().bytes()},
            {"terminal-lineage-digest", result.post_lineage().bytes()},
        }));
    require(
        sha256_bytes(completion_v2_bytes) == completion.digest().bytes(),
        "completion digest is not pinned to audit-replay-completion-v2");
    require(
        AuditReplayTestAuthority2::completion_digest(
            changed_input, seed.digest(), request, 1, result.post_lineage())
                .bytes()
                != completion.digest().bytes()
            && AuditReplayTestAuthority2::completion_digest(
                   replay_input, changed_seed.digest(), request, 1,
                   result.post_lineage())
                   .bytes()
                != completion.digest().bytes()
            && AuditReplayTestAuthority2::completion_digest(
                   replay_input, seed.digest(), changed_request, 1,
                   result.post_lineage())
                .bytes()
                != completion.digest().bytes()
            && AuditReplayTestAuthority2::completion_digest(
                   replay_input, seed.digest(), request, 2,
                   result.post_lineage())
                   .bytes()
                != completion.digest().bytes()
            && AuditReplayTestAuthority2::completion_digest(
                   replay_input, seed.digest(), request, 1,
                   foreign_lineage.digest())
                   .bytes()
                != completion.digest().bytes(),
        "completion v2 digest omits input, seed, request, count, or terminal lineage");
}

} // namespace

void audit_replay_identity_lineage_gate()
{
    limits_are_request_and_decision_identity_gate();
    request_preflight_gate();
    lineage_component_mutation_gate();
    no_material_and_nonengaging_lineage_gate();
    result_and_completion_component_gate();
}
