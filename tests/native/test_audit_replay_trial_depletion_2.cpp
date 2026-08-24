#include "audit_replay_fixtures_2.h"

#include "audit_depletion_witness_2.h"
#include "audit_stock_state_identity_2.h"
#include "exact_disk_region_2.h"
#include "stock_exact_depletion_2.h"

using namespace audit_replay_fixtures;

namespace {

std::string stock_state_digest_bytes(const Stock2& stock)
{
    return AuditStockStateIdentity2::build(stock).digest().bytes();
}

void require_stock_parity(
    const Stock2& actual,
    const Stock2& expected,
    const char* point_set_message,
    const char* digest_message,
    const char* representation_message)
{
    require(actual.exactly_equals(expected), point_set_message);
    require(
        stock_state_digest_bytes(actual) == stock_state_digest_bytes(expected),
        digest_message);
    require(
        actual.representation_is_valid() && expected.representation_is_valid(),
        representation_message);
}

void require_witness_transition(
    const AuditDepletionWitness2& witness,
    AuditDepletionKind2 kind,
    const Stock2& authority,
    const Stock2& trial,
    const char* message)
{
    require(
        witness.kind() == kind
            && witness.pre_stock_digest().bytes()
                == stock_state_digest_bytes(authority)
            && witness.post_stock_digest().bytes()
                == stock_state_digest_bytes(trial)
            && witness.digest().bytes().size() == 32,
        message);
}

void require_zero_trial_internal_work(const char* message)
{
    const AuditTrialDepletionInstrumentation2 metrics =
        audit_trial_depletion_instrumentation_for_test();
    require(
        metrics.internal_clone_count == 0 && metrics.internal_swap_count == 0,
        message);
}

void exact_disk_guard_gate()
{
    // Kills empty and exactly nonpositive disk-union construction.
    require_throws<ExactDiskRegionEmptyCentersError>(
        [&] {
            static_cast<void>(build_exact_disk_union_region_2({}, Epeck::FT(1)));
        },
        "empty exact disk region was accepted");
    require_throws<ExactDiskRegionRadiusError>(
        [&] {
            static_cast<void>(build_exact_disk_union_region_2(
                {EPoint(0, 0)}, Epeck::FT(0)));
        },
        "zero exact disk radius was accepted");
    require_throws<ExactDiskRegionRadiusError>(
        [&] {
            static_cast<void>(build_exact_disk_union_region_2(
                {EPoint(0, 0)}, Epeck::FT(-1)));
        },
        "negative exact disk radius was accepted");
}

void segment_trial_gate()
{
    // Kills authority mutation, incomplete witness binding, hidden work during
    // validation, and divergence from the independent legacy exact path.
    const compas::RowMatrixXd boundary = rectangle(-5.0, -5.0, 5.0, 5.0);
    const AuditPolicy2 audit_policy = policy();
    const AuditSegmentMotion2 motion = segment();
    const Stock2 authority(boundary, {});
    Stock2 pristine = authority.clone();
    Stock2 trial = authority.clone();
    Stock2 legacy = authority.clone();
    reset_audit_trial_depletion_instrumentation_for_test();
    const AuditDepletionWitness2 witness =
        apply_audit_segment_depletion_to_trial(
            authority, trial, motion, audit_policy);
    static_cast<void>(legacy.subtract_exact_segment(
        motion.xy(),
        audit_policy.tool_radius_mm(),
        audit_policy.depletion_chord_bound_mm(),
        audit_policy.center_count_limit()));
    require_stock_parity(
        authority,
        pristine,
        "segment applicator mutated authority point set",
        "segment applicator mutated authority state digest",
        "segment authority representation is invalid");
    require(!trial.exactly_equals(authority),
            "intersecting segment applicator was a no-op");
    require_stock_parity(
        trial,
        legacy,
        "segment applicator diverges from legacy point set",
        "segment applicator diverges from legacy stock-state digest",
        "segment trial or legacy representation is invalid");
    require_witness_transition(
        witness,
        AuditDepletionKind2::SEGMENT,
        authority,
        trial,
        "segment witness omits kind, pre/post stock, or digest");
    require(validate_audit_segment_depletion_witness(
                authority, trial, motion, audit_policy, witness),
            "fresh segment witness failed replay");

    Stock2 foreign_trial = authority.clone();
    foreign_trial.subtract_disk(4.0, 4.0, 0.25);
    Stock2 foreign_trial_pristine = foreign_trial.clone();
    require_throws<AuditTrialStockMismatchError>(
        [&] {
            static_cast<void>(apply_audit_segment_depletion_to_trial(
                authority, foreign_trial, motion, audit_policy));
        },
        "foreign segment trial was accepted");
    require_stock_parity(
        foreign_trial,
        foreign_trial_pristine,
        "rejected foreign segment trial was mutated",
        "rejected foreign segment trial digest changed",
        "rejected foreign segment trial representation is invalid");

    Stock2 foreign_authority = authority.clone();
    foreign_authority.subtract_disk(4.0, 4.0, 0.25);
    Stock2 foreign_post = trial.clone();
    foreign_post.subtract_disk(-4.0, -4.0, 0.25);
    require(
        !validate_audit_segment_depletion_witness(
            authority, trial, segment(1.0), audit_policy, witness),
        "segment witness replay accepted foreign motion");
    require(
        !validate_audit_segment_depletion_witness(
            authority, trial, motion, policy(2048), witness),
        "segment witness replay accepted foreign policy");
    require(
        !validate_audit_segment_depletion_witness(
            foreign_authority, trial, motion, audit_policy, witness),
        "segment witness replay accepted foreign authority");
    require(
        !validate_audit_segment_depletion_witness(
            authority, foreign_post, motion, audit_policy, witness),
        "segment witness replay accepted foreign post stock");

    Stock2 wrong_kind_trial = authority.clone();
    const AuditDepletionWitness2 wrong_kind =
        apply_audit_plunge_depletion_to_trial(
            authority, wrong_kind_trial, plunge(), audit_policy);
    require(
        !validate_audit_segment_depletion_witness(
            authority, trial, motion, audit_policy, wrong_kind),
        "segment validator accepted a plunge witness kind");
    require_zero_trial_internal_work(
        "segment apply/validation cloned or swapped internally");
}

void circle_trial_gate()
{
    // Kills full-circle divergence and omission of any replay-bound input.
    const compas::RowMatrixXd boundary = rectangle(-5.0, -5.0, 5.0, 5.0);
    const AuditPolicy2 audit_policy = policy();
    const AuditCircleMotion2 motion = circle();
    const Stock2 authority(boundary, {});
    Stock2 pristine = authority.clone();
    Stock2 trial = authority.clone();
    Stock2 legacy = authority.clone();
    reset_audit_trial_depletion_instrumentation_for_test();
    const AuditDepletionWitness2 witness =
        apply_audit_circle_depletion_to_trial(
            authority, trial, motion, audit_policy);
    static_cast<void>(legacy.subtract_exact_full_circle(
        motion.xy(),
        audit_policy.tool_radius_mm(),
        audit_policy.depletion_chord_bound_mm(),
        audit_policy.center_count_limit()));
    require_stock_parity(
        authority,
        pristine,
        "circle applicator mutated authority point set",
        "circle applicator mutated authority state digest",
        "circle authority representation is invalid");
    require(!trial.exactly_equals(authority),
            "intersecting circle applicator was a no-op");
    require_stock_parity(
        trial,
        legacy,
        "circle applicator diverges from legacy point set",
        "circle applicator diverges from legacy stock-state digest",
        "circle trial or legacy representation is invalid");
    require_witness_transition(
        witness,
        AuditDepletionKind2::FULL_CIRCLE,
        authority,
        trial,
        "circle witness omits kind, pre/post stock, or digest");
    require(validate_audit_circle_depletion_witness(
                authority, trial, motion, audit_policy, witness),
            "fresh circle witness failed replay");

    Stock2 foreign_trial = authority.clone();
    foreign_trial.subtract_disk(4.0, 4.0, 0.25);
    Stock2 foreign_trial_pristine = foreign_trial.clone();
    require_throws<AuditTrialStockMismatchError>(
        [&] {
            static_cast<void>(apply_audit_circle_depletion_to_trial(
                authority, foreign_trial, motion, audit_policy));
        },
        "foreign circle trial was accepted");
    require_stock_parity(
        foreign_trial,
        foreign_trial_pristine,
        "rejected foreign circle trial was mutated",
        "rejected foreign circle trial digest changed",
        "rejected foreign circle trial representation is invalid");

    Stock2 foreign_authority = authority.clone();
    foreign_authority.subtract_disk(4.0, 4.0, 0.25);
    Stock2 foreign_post = trial.clone();
    foreign_post.subtract_disk(-4.0, -4.0, 0.25);
    require(
        !validate_audit_circle_depletion_witness(
            authority, trial, circle(true), audit_policy, witness),
        "circle witness replay accepted foreign motion");
    require(
        !validate_audit_circle_depletion_witness(
            authority, trial, motion, policy(2048), witness),
        "circle witness replay accepted foreign policy");
    require(
        !validate_audit_circle_depletion_witness(
            foreign_authority, trial, motion, audit_policy, witness),
        "circle witness replay accepted foreign authority");
    require(
        !validate_audit_circle_depletion_witness(
            authority, foreign_post, motion, audit_policy, witness),
        "circle witness replay accepted foreign post stock");

    Stock2 wrong_kind_trial = authority.clone();
    const AuditDepletionWitness2 wrong_kind =
        apply_audit_segment_depletion_to_trial(
            authority, wrong_kind_trial, segment(), audit_policy);
    require(
        !validate_audit_circle_depletion_witness(
            authority, trial, motion, audit_policy, wrong_kind),
        "circle validator accepted a segment witness kind");
    require_zero_trial_internal_work(
        "circle apply/validation cloned or swapped internally");
}

void arc_trial_gate()
{
    // Kills Task3 trace bypass, incomplete witness replay, or legacy sweep use.
    const compas::RowMatrixXd boundary = rectangle(-5.0, -5.0, 5.0, 5.0);
    const AuditPolicy2 audit_policy = policy();
    const AuditArcMotion2 motion = arc();
    const Stock2 authority(boundary, {});
    Stock2 pristine = authority.clone();
    Stock2 trial = authority.clone();
    Stock2 legacy = authority.clone();
    reset_audit_trial_depletion_instrumentation_for_test();
    reset_legacy_arc_sweep_reachability_for_test();
    const AuditDepletionWitness2 witness = apply_audit_arc_depletion_to_trial(
        authority, trial, motion, audit_policy);
    static_cast<void>(legacy.subtract_exact_arc(
        motion,
        audit_policy.tool_radius_mm(),
        audit_policy.depletion_chord_bound_mm(),
        audit_policy.center_count_limit()));
    require_stock_parity(
        authority,
        pristine,
        "arc applicator mutated authority point set",
        "arc applicator mutated authority state digest",
        "arc authority representation is invalid");
    require(!trial.exactly_equals(authority),
            "intersecting arc applicator was a no-op");
    require_stock_parity(
        trial,
        legacy,
        "arc applicator diverges from Task3 exact point set",
        "arc applicator diverges from Task3 stock-state digest",
        "arc trial or legacy representation is invalid");
    require(legacy_arc_sweep_reachability_for_test() == 0,
            "authoritative exact arc path reached subtract_arc_sweep");
    require_witness_transition(
        witness,
        AuditDepletionKind2::ARC,
        authority,
        trial,
        "arc witness omits kind, pre/post stock, or digest");
    require(validate_audit_arc_depletion_witness(
                authority, trial, motion, audit_policy, witness),
            "fresh arc witness failed replay");

    Stock2 foreign_trial = authority.clone();
    foreign_trial.subtract_disk(4.0, 4.0, 0.25);
    Stock2 foreign_trial_pristine = foreign_trial.clone();
    require_throws<AuditTrialStockMismatchError>(
        [&] {
            static_cast<void>(apply_audit_arc_depletion_to_trial(
                authority, foreign_trial, motion, audit_policy));
        },
        "foreign arc trial was accepted");
    require_stock_parity(
        foreign_trial,
        foreign_trial_pristine,
        "rejected foreign arc trial was mutated",
        "rejected foreign arc trial digest changed",
        "rejected foreign arc trial representation is invalid");

    const AuditArcMotion2 foreign_motion = AuditArcMotion2::build(
        motion.center(),
        motion.zero_phase(),
        motion.guide_radius(),
        0.0,
        std::numbers::pi,
        motion.clockwise(),
        motion.cut_z());
    Stock2 foreign_authority = authority.clone();
    foreign_authority.subtract_disk(4.0, 4.0, 0.25);
    Stock2 foreign_post = trial.clone();
    foreign_post.subtract_disk(-4.0, -4.0, 0.25);
    require(
        !validate_audit_arc_depletion_witness(
            authority, trial, foreign_motion, audit_policy, witness),
        "arc witness replay accepted foreign motion");
    require(
        !validate_audit_arc_depletion_witness(
            authority, trial, motion, policy(2048), witness),
        "arc witness replay accepted foreign policy");
    require(
        !validate_audit_arc_depletion_witness(
            foreign_authority, trial, motion, audit_policy, witness),
        "arc witness replay accepted foreign authority");
    require(
        !validate_audit_arc_depletion_witness(
            authority, foreign_post, motion, audit_policy, witness),
        "arc witness replay accepted foreign post stock");

    Stock2 wrong_kind_trial = authority.clone();
    const AuditDepletionWitness2 wrong_kind =
        apply_audit_circle_depletion_to_trial(
            authority, wrong_kind_trial, circle(), audit_policy);
    require(
        !validate_audit_arc_depletion_witness(
            authority, trial, motion, audit_policy, wrong_kind),
        "arc validator accepted a circle witness kind");
    require_zero_trial_internal_work(
        "arc apply/validation cloned or swapped internally");
}

void plunge_trial_gate()
{
    // Kills approximate plunge geometry and omission of any witness-bound input.
    const compas::RowMatrixXd boundary = rectangle(-5.0, -5.0, 5.0, 5.0);
    const AuditPolicy2 audit_policy = policy();
    const AuditVerticalPlunge2 motion = plunge();
    const Stock2 authority(boundary, {});
    Stock2 pristine = authority.clone();
    Stock2 trial = authority.clone();
    Stock2 legacy = authority.clone();
    reset_audit_trial_depletion_instrumentation_for_test();
    const AuditDepletionWitness2 witness =
        apply_audit_plunge_depletion_to_trial(
            authority, trial, motion, audit_policy);
    // This binary64 reference is independent of the new exact-disk builder. All
    // values are dyadic, so legacy injection describes the identical exact disk.
    legacy.subtract_disk(3.0, 3.0, 0.5);
    require_stock_parity(
        authority,
        pristine,
        "plunge applicator mutated authority point set",
        "plunge applicator mutated authority state digest",
        "plunge authority representation is invalid");
    require(!trial.exactly_equals(authority),
            "intersecting plunge applicator was a no-op");
    require_stock_parity(
        trial,
        legacy,
        "plunge applicator diverges from legacy disk point set",
        "plunge applicator diverges from legacy disk state digest",
        "plunge trial or legacy representation is invalid");
    require_witness_transition(
        witness,
        AuditDepletionKind2::PLUNGE,
        authority,
        trial,
        "plunge witness omits kind, pre/post stock, or digest");
    require(validate_audit_plunge_depletion_witness(
                authority, trial, motion, audit_policy, witness),
            "fresh plunge witness failed replay");

    Stock2 foreign_trial = authority.clone();
    foreign_trial.subtract_disk(4.0, 4.0, 0.25);
    Stock2 foreign_trial_pristine = foreign_trial.clone();
    require_throws<AuditTrialStockMismatchError>(
        [&] {
            static_cast<void>(apply_audit_plunge_depletion_to_trial(
                authority, foreign_trial, motion, audit_policy));
        },
        "foreign plunge trial was accepted");
    require_stock_parity(
        foreign_trial,
        foreign_trial_pristine,
        "rejected foreign plunge trial was mutated",
        "rejected foreign plunge trial digest changed",
        "rejected foreign plunge trial representation is invalid");

    const AuditVerticalPlunge2 foreign_motion =
        std::get<AuditVerticalPlunge2>(classify_audit_line(
            {4.0, 3.0, 5.0},
            {4.0, 3.0, 0.0},
            0.0,
            5.0,
            "plunge"));
    Stock2 foreign_authority = authority.clone();
    foreign_authority.subtract_disk(4.0, 4.0, 0.25);
    Stock2 foreign_post = trial.clone();
    foreign_post.subtract_disk(-4.0, -4.0, 0.25);
    require(
        !validate_audit_plunge_depletion_witness(
            authority, trial, foreign_motion, audit_policy, witness),
        "plunge witness replay accepted foreign motion");
    require(
        !validate_audit_plunge_depletion_witness(
            authority, trial, motion, policy(2048), witness),
        "plunge witness replay accepted foreign policy");
    require(
        !validate_audit_plunge_depletion_witness(
            foreign_authority, trial, motion, audit_policy, witness),
        "plunge witness replay accepted foreign authority");
    require(
        !validate_audit_plunge_depletion_witness(
            authority, foreign_post, motion, audit_policy, witness),
        "plunge witness replay accepted foreign post stock");

    Stock2 wrong_kind_trial = authority.clone();
    const AuditDepletionWitness2 wrong_kind =
        apply_audit_segment_depletion_to_trial(
            authority, wrong_kind_trial, segment(), audit_policy);
    require(
        !validate_audit_plunge_depletion_witness(
            authority, trial, motion, audit_policy, wrong_kind),
        "plunge validator accepted a segment witness kind");
    require_zero_trial_internal_work(
        "plunge apply/validation cloned or swapped internally");
}

template <class Motion, class Apply, class Validate>
void require_no_material_removal(
    const Motion& motion,
    AuditDepletionKind2 kind,
    Apply&& apply,
    Validate&& validate,
    const char* parity_message,
    const char* witness_message,
    const char* validation_message,
    const char* instrumentation_message)
{
    const compas::RowMatrixXd boundary = rectangle(8.0, 8.0, 10.0, 10.0);
    const AuditPolicy2 audit_policy = policy();
    const Stock2 authority(boundary, {});
    Stock2 trial = authority.clone();
    reset_audit_trial_depletion_instrumentation_for_test();

    const AuditDepletionWitness2 witness =
        std::forward<Apply>(apply)(authority, trial, motion, audit_policy);
    require_stock_parity(
        trial,
        authority,
        parity_message,
        "no-material removal changed stock-state digest",
        "no-material removal invalidated stock representation");
    require_witness_transition(
        witness, kind, authority, trial, witness_message);
    require(
        std::forward<Validate>(validate)(
            authority, trial, motion, audit_policy, witness),
        validation_message);
    require_zero_trial_internal_work(instrumentation_message);
}

void no_material_removal_gate()
{
    // Kills an applicator that requires a nonempty boolean difference to emit a
    // complete, replayable witness for an authenticated mutating operation.
    require_no_material_removal(
        segment(),
        AuditDepletionKind2::SEGMENT,
        apply_audit_segment_depletion_to_trial,
        validate_audit_segment_depletion_witness,
        "clear segment changed stock point set",
        "clear segment witness lost unchanged transition",
        "clear segment witness failed replay",
        "clear segment apply/validation performed internal clone or swap");
    require_no_material_removal(
        circle(),
        AuditDepletionKind2::FULL_CIRCLE,
        apply_audit_circle_depletion_to_trial,
        validate_audit_circle_depletion_witness,
        "clear circle changed stock point set",
        "clear circle witness lost unchanged transition",
        "clear circle witness failed replay",
        "clear circle apply/validation performed internal clone or swap");
    require_no_material_removal(
        arc(),
        AuditDepletionKind2::ARC,
        apply_audit_arc_depletion_to_trial,
        validate_audit_arc_depletion_witness,
        "clear arc changed stock point set",
        "clear arc witness lost unchanged transition",
        "clear arc witness failed replay",
        "clear arc apply/validation performed internal clone or swap");
    require_no_material_removal(
        plunge(),
        AuditDepletionKind2::PLUNGE,
        apply_audit_plunge_depletion_to_trial,
        validate_audit_plunge_depletion_witness,
        "clear plunge changed stock point set",
        "clear plunge witness lost unchanged transition",
        "clear plunge witness failed replay",
        "clear plunge apply/validation performed internal clone or swap");
}

template <class Motion, class Apply, class Legacy>
void require_center_limit_preserved(
    const Motion& motion,
    Apply&& apply,
    Legacy&& legacy_apply,
    const char* trial_error_message,
    const char* legacy_error_message,
    const char* state_message,
    const char* instrumentation_message)
{
    const compas::RowMatrixXd boundary = rectangle(-5.0, -5.0, 5.0, 5.0);
    const AuditPolicy2 bounded_policy = policy(2);
    const Stock2 authority(boundary, {});
    Stock2 trial = authority.clone();
    Stock2 legacy = authority.clone();
    reset_audit_trial_depletion_instrumentation_for_test();

    require_throws<ExactDepletionCenterLimitError>(
        [&] {
            static_cast<void>(std::forward<Apply>(apply)(
                authority, trial, motion, bounded_policy));
        },
        trial_error_message);
    require_throws<ExactDepletionCenterLimitError>(
        [&] {
            static_cast<void>(std::forward<Legacy>(legacy_apply)(
                legacy, motion, bounded_policy));
        },
        legacy_error_message);
    require_stock_parity(
        trial,
        authority,
        state_message,
        "center-limit rejection changed trial stock-state digest",
        "center-limit rejection invalidated trial representation");
    require_stock_parity(
        legacy,
        authority,
        "legacy center-limit rejection changed stock point set",
        "legacy center-limit rejection changed stock-state digest",
        "legacy center-limit rejection invalidated representation");
    require_zero_trial_internal_work(instrumentation_message);
}

void construction_limit_gate()
{
    // Kills bypass of the existing exact construction bound and verifies that a
    // one-center plunge remains admissible under the same sealed policy field.
    require_center_limit_preserved(
        segment(),
        apply_audit_segment_depletion_to_trial,
        [](Stock2& stock,
           const AuditSegmentMotion2& motion,
           const AuditPolicy2& audit_policy) {
            return stock.subtract_exact_segment(
                motion.xy(),
                audit_policy.tool_radius_mm(),
                audit_policy.depletion_chord_bound_mm(),
                audit_policy.center_count_limit());
        },
        "segment applicator ignored exact center-count limit",
        "legacy segment unexpectedly fits exact center-count limit",
        "failed bounded segment applicator changed trial point set",
        "bounded segment applicator cloned or swapped internally");
    require_center_limit_preserved(
        circle(),
        apply_audit_circle_depletion_to_trial,
        [](Stock2& stock,
           const AuditCircleMotion2& motion,
           const AuditPolicy2& audit_policy) {
            return stock.subtract_exact_full_circle(
                motion.xy(),
                audit_policy.tool_radius_mm(),
                audit_policy.depletion_chord_bound_mm(),
                audit_policy.center_count_limit());
        },
        "circle applicator ignored exact center-count limit",
        "legacy circle unexpectedly fits exact center-count limit",
        "failed bounded circle applicator changed trial point set",
        "bounded circle applicator cloned or swapped internally");
    require_center_limit_preserved(
        arc(),
        apply_audit_arc_depletion_to_trial,
        [](Stock2& stock,
           const AuditArcMotion2& motion,
           const AuditPolicy2& audit_policy) {
            return stock.subtract_exact_arc(
                motion,
                audit_policy.tool_radius_mm(),
                audit_policy.depletion_chord_bound_mm(),
                audit_policy.center_count_limit());
        },
        "arc applicator ignored exact center-count limit",
        "Task3 arc unexpectedly fits exact center-count limit",
        "failed bounded arc applicator changed trial point set",
        "bounded arc applicator cloned or swapped internally");

    const compas::RowMatrixXd boundary = rectangle(-5.0, -5.0, 5.0, 5.0);
    const AuditPolicy2 one_center_policy = policy(1);
    const AuditVerticalPlunge2 motion = plunge();
    const Stock2 authority(boundary, {});
    Stock2 trial = authority.clone();
    Stock2 legacy = authority.clone();
    reset_audit_trial_depletion_instrumentation_for_test();
    const AuditDepletionWitness2 witness =
        apply_audit_plunge_depletion_to_trial(
            authority, trial, motion, one_center_policy);
    legacy.subtract_disk(3.0, 3.0, 0.5);
    require_stock_parity(
        trial,
        legacy,
        "one-center plunge diverges from legacy point set",
        "one-center plunge diverges from legacy stock-state digest",
        "one-center plunge invalidated representation");
    require(
        validate_audit_plunge_depletion_witness(
            authority, trial, motion, one_center_policy, witness),
        "one-center plunge witness failed replay");
    require_zero_trial_internal_work(
        "one-center plunge apply/validation cloned or swapped internally");
}

} // namespace

void audit_replay_trial_depletion_gate()
{
    exact_disk_guard_gate();
    segment_trial_gate();
    circle_trial_gate();
    arc_trial_gate();
    plunge_trial_gate();
    no_material_removal_gate();
    construction_limit_gate();
}
