#include "stock_exact_depletion_2.h"

#include "canonical_encoding.h"
#include "audit_replay_internal_2.h"
#include "exact_depletion_2.h"
#include "exact_disk_region_2.h"

#include <utility>
#include <vector>

namespace {

thread_local AuditTrialDepletionInstrumentation2 trial_metrics{0, 0};
thread_local std::size_t legacy_arc_sweep_count = 0;
thread_local std::size_t trial_primitive_scope_depth = 0;

class AuditTrialPrimitiveScope2 {
public:
    AuditTrialPrimitiveScope2() noexcept
    {
        ++trial_primitive_scope_depth;
    }

    ~AuditTrialPrimitiveScope2()
    {
        --trial_primitive_scope_depth;
    }

    AuditTrialPrimitiveScope2(const AuditTrialPrimitiveScope2&) = delete;
    AuditTrialPrimitiveScope2& operator=(
        const AuditTrialPrimitiveScope2&) = delete;
};

std::string exact_size_bytes(std::size_t value)
{
    return canonical_audit_rational_bytes(Epeck::FT(value));
}

std::string kind_bytes(AuditDepletionKind2 kind)
{
    switch (kind) {
    case AuditDepletionKind2::SEGMENT:
        return canonical_encode_bytes("segment");
    case AuditDepletionKind2::FULL_CIRCLE:
        return canonical_encode_bytes("full-circle");
    case AuditDepletionKind2::ARC:
        return canonical_encode_bytes("arc");
    case AuditDepletionKind2::PLUNGE:
        return canonical_encode_bytes("plunge");
    }
    throw AuditDepletionWitnessError("unknown audit depletion kind");
}

void validate_legacy_trace(const DepletionTrace& trace)
{
    if (trace.center_count != trace.center_parameters.size()
        || !trace.exact_incidence
        || !trace.exact_parameters_in_range
        || !trace.exact_anchors_present
        || !trace.exact_removal_radius_valid
        || !trace.exact_chord_bound_holds
        || (trace.cyclic && !trace.exact_seam_chord_bound_holds)
        || trace.strategy_version.empty()) {
        throw AuditDepletionWitnessError(
            "exact depletion construction trace is structurally invalid");
    }
}

std::string parameter_bytes(const ExactCenterParameter2& parameter)
{
    return canonical_encode_component_map({
        {"chart", canonical_audit_rational_bytes(Epeck::FT(parameter.chart))},
        {"denominator", exact_size_bytes(parameter.denominator)},
        {"numerator", exact_size_bytes(parameter.numerator)},
    });
}

std::string legacy_trace_bytes(
    AuditDepletionKind2 kind,
    const NativeMotionDigest2& motion_digest,
    const AuditPolicy2& policy,
    const DepletionTrace& trace)
{
    std::vector<std::string> parameters;
    parameters.reserve(trace.center_parameters.size());
    for (const ExactCenterParameter2& parameter : trace.center_parameters) {
        parameters.push_back(parameter_bytes(parameter));
    }
    return canonical_encode_tagged_union(
        "exact-depletion-trace-v1",
        canonical_encode_component_map({
            {"center-count", exact_size_bytes(trace.center_count)},
            {"center-count-limit", exact_size_bytes(policy.center_count_limit())},
            {"cyclic", canonical_encode_boolean(trace.cyclic)},
            {"kind", kind_bytes(kind)},
            {"max-chord", canonical_audit_rational_bytes(trace.max_chord)},
            {"motion-digest", motion_digest.bytes()},
            {"parameters", canonical_encode_sequence(parameters)},
            {"removal-radius", canonical_audit_rational_bytes(trace.removal_radius)},
            {"strategy", canonical_encode_bytes(trace.strategy_version)},
        }));
}

std::string plunge_trace_bytes(
    const AuditVerticalPlunge2& motion,
    const AuditPolicy2& policy)
{
    return canonical_encode_tagged_union(
        "exact-plunge-disk-v1",
        canonical_encode_component_map({
            {"center", canonical_encode_component_map({
                           {"x", canonical_audit_rational_bytes(
                                     motion.cut_endpoint().x())},
                           {"y", canonical_audit_rational_bytes(
                                     motion.cut_endpoint().y())},
                       })},
            {"motion-digest", motion.digest().bytes()},
            {"radius", canonical_audit_rational_bytes(policy.tool_radius_mm())},
            {"strategy", canonical_encode_bytes("exact-plunge-disk-v1")},
        }));
}

} // namespace

class AuditDepletionWitnessFactory2 {
public:
    static AuditDepletionWitness2 build(
        AuditDepletionKind2 kind,
        const Stock2& authority,
        const Stock2& post,
        const NativeMotionDigest2& motion_digest,
        const AuditPolicy2& policy,
        const std::string& construction_bytes,
        const std::string& strategy_version)
    {
        const AuditStockStateIdentity2 pre_identity =
            AuditStockStateIdentity2::build(authority);
        const AuditStockStateIdentity2 post_identity =
            AuditStockStateIdentity2::build(post);
        const ExactDepletionTraceDigest2 construction_digest =
            ExactDepletionTraceDigestAuthority2::hash_canonical(
                construction_bytes);
        std::string canonical = canonical_encode_tagged_union(
            "audit-depletion-witness-v1",
            canonical_encode_component_map({
                {"construction-digest", construction_digest.bytes()},
                {"kind", kind_bytes(kind)},
                {"motion-digest", motion_digest.bytes()},
                {"policy-digest", policy.digest().bytes()},
                {"post-stock-digest", post_identity.digest().bytes()},
                {"pre-stock-digest", pre_identity.digest().bytes()},
                {"strategy", canonical_encode_bytes(strategy_version)},
            }));
        return AuditDepletionWitness2(
            kind,
            pre_identity.digest(),
            post_identity.digest(),
            motion_digest,
            policy.digest(),
            construction_digest,
            strategy_version,
            canonical,
            DepletionWitnessDigestAuthority2::hash_canonical(canonical));
    }

    static bool matches(
        AuditDepletionKind2 kind,
        const Stock2& authority,
        const Stock2& post,
        const NativeMotionDigest2& motion_digest,
        const AuditPolicy2& policy,
        const std::string& construction_bytes,
        const std::string& strategy_version,
        const AuditDepletionWitness2& witness)
    {
        const AuditDepletionWitness2 expected = build(
            kind,
            authority,
            post,
            motion_digest,
            policy,
            construction_bytes,
            strategy_version);
        return witness.kind_ == expected.kind_
            && witness.pre_stock_digest_.bytes()
                == expected.pre_stock_digest_.bytes()
            && witness.post_stock_digest_.bytes()
                == expected.post_stock_digest_.bytes()
            && witness.motion_digest_.bytes()
                == expected.motion_digest_.bytes()
            && witness.policy_digest_.bytes()
                == expected.policy_digest_.bytes()
            && witness.construction_digest_.bytes()
                == expected.construction_digest_.bytes()
            && witness.strategy_version_ == expected.strategy_version_
            && witness.canonical_bytes_ == expected.canonical_bytes_
            && witness.digest_.bytes() == expected.digest_.bytes();
    }
};

namespace {

void require_equal_trial(const Stock2& authority, const Stock2& trial)
{
    if (!authority.representation_is_valid()
        || !trial.representation_is_valid()
        || !trial.exactly_equals(authority)
        || AuditStockStateIdentity2::build(trial).digest().bytes()
            != AuditStockStateIdentity2::build(authority).digest().bytes()) {
        throw AuditTrialStockMismatchError(
            "audit depletion trial and authority must be canonical exact equals");
    }
}

bool gps_sets_equal(const Gps& first, const Gps& second)
{
    Gps first_difference(first);
    first_difference.difference(second);
    if (!first_difference.is_empty()) {
        return false;
    }
    Gps second_difference(second);
    second_difference.difference(first);
    return second_difference.is_empty();
}

bool expected_post_matches(
    const Stock2& authority,
    const Stock2& post,
    Gps removal)
{
    if (!authority.representation_is_valid()
        || !post.representation_is_valid()) {
        return false;
    }
    Gps expected(authority.set());
    expected.difference(removal);
    return gps_sets_equal(expected, post.set());
}

ExactDepletionConstruction2 segment_construction(
    const AuditSegmentMotion2& motion,
    const AuditPolicy2& policy)
{
    ExactDepletionConstruction2 construction = construct_exact_segment_depletion(
        motion.xy(),
        policy.tool_radius_mm(),
        policy.depletion_chord_bound_mm(),
        policy.center_count_limit());
    validate_legacy_trace(construction.trace);
    return construction;
}

ExactDepletionConstruction2 circle_construction(
    const AuditCircleMotion2& motion,
    const AuditPolicy2& policy)
{
    ExactDepletionConstruction2 construction =
        construct_exact_full_circle_depletion(
            motion.xy(),
            policy.tool_radius_mm(),
            policy.depletion_chord_bound_mm(),
            policy.center_count_limit());
    validate_legacy_trace(construction.trace);
    return construction;
}

ExactArcDepletionConstruction2 arc_construction(
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy)
{
    ExactArcDepletionConstruction2 construction = construct_exact_arc_depletion(
        motion,
        policy.tool_radius_mm(),
        policy.depletion_chord_bound_mm(),
        policy.center_count_limit());
    if (!exact_arc_structural_density_holds(
            motion,
            policy.depletion_chord_bound_mm(),
            construction.trace.parameters())
        || !construction.trace.matches_exact_inputs(
            policy.tool_radius_mm(),
            policy.depletion_chord_bound_mm(),
            policy.center_count_limit())
        || !construction.trace.matches_motion(motion)) {
        throw AuditDepletionWitnessError(
            "exact arc depletion trace failed structural validation");
    }
    return construction;
}

} // namespace

AuditDepletionWitness2 apply_audit_segment_depletion_to_trial(
    const Stock2& authority,
    Stock2& trial,
    const AuditSegmentMotion2& motion,
    const AuditPolicy2& policy)
{
    const AuditTrialPrimitiveScope2 scope;
    require_equal_trial(authority, trial);
    ExactDepletionConstruction2 construction = segment_construction(motion, policy);
    const std::string trace = legacy_trace_bytes(
        AuditDepletionKind2::SEGMENT,
        motion.digest(),
        policy,
        construction.trace);
    Gps removal = build_exact_disk_union_region_2(
        construction.centers, policy.tool_radius_mm());
    trial.set().difference(removal);
    return AuditDepletionWitnessFactory2::build(
        AuditDepletionKind2::SEGMENT,
        authority,
        trial,
        motion.digest(),
        policy,
        trace,
        construction.trace.strategy_version);
}

AuditDepletionWitness2 apply_audit_circle_depletion_to_trial(
    const Stock2& authority,
    Stock2& trial,
    const AuditCircleMotion2& motion,
    const AuditPolicy2& policy)
{
    const AuditTrialPrimitiveScope2 scope;
    require_equal_trial(authority, trial);
    ExactDepletionConstruction2 construction = circle_construction(motion, policy);
    const std::string trace = legacy_trace_bytes(
        AuditDepletionKind2::FULL_CIRCLE,
        motion.digest(),
        policy,
        construction.trace);
    Gps removal = build_exact_disk_union_region_2(
        construction.centers, policy.tool_radius_mm());
    trial.set().difference(removal);
    return AuditDepletionWitnessFactory2::build(
        AuditDepletionKind2::FULL_CIRCLE,
        authority,
        trial,
        motion.digest(),
        policy,
        trace,
        construction.trace.strategy_version);
}

AuditDepletionWitness2 apply_audit_arc_depletion_to_trial(
    const Stock2& authority,
    Stock2& trial,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy)
{
    const AuditTrialPrimitiveScope2 scope;
    require_equal_trial(authority, trial);
    ExactArcDepletionConstruction2 construction = arc_construction(motion, policy);
    Gps removal = build_exact_disk_union_region_2(
        construction.centers, policy.tool_radius_mm());
    trial.set().difference(removal);
    return AuditDepletionWitnessFactory2::build(
        AuditDepletionKind2::ARC,
        authority,
        trial,
        motion.digest(),
        policy,
        construction.trace.canonical_bytes(),
        construction.trace.strategy_version());
}

AuditDepletionWitness2 apply_audit_plunge_depletion_to_trial(
    const Stock2& authority,
    Stock2& trial,
    const AuditVerticalPlunge2& motion,
    const AuditPolicy2& policy)
{
    const AuditTrialPrimitiveScope2 scope;
    require_equal_trial(authority, trial);
    const std::string trace = plunge_trace_bytes(motion, policy);
    Gps removal = build_exact_disk_union_region_2(
        {motion.cut_endpoint()}, policy.tool_radius_mm());
    trial.set().difference(removal);
    return AuditDepletionWitnessFactory2::build(
        AuditDepletionKind2::PLUNGE,
        authority,
        trial,
        motion.digest(),
        policy,
        trace,
        "exact-plunge-disk-v1");
}

bool validate_audit_segment_depletion_witness(
    const Stock2& authority,
    const Stock2& post,
    const AuditSegmentMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDepletionWitness2& witness)
{
    const AuditTrialPrimitiveScope2 scope;
    try {
        const ExactDepletionConstruction2 construction =
            segment_construction(motion, policy);
        Gps removal = build_exact_disk_union_region_2(
            construction.centers, policy.tool_radius_mm());
        return expected_post_matches(authority, post, std::move(removal))
            && AuditDepletionWitnessFactory2::matches(
            AuditDepletionKind2::SEGMENT,
            authority,
            post,
            motion.digest(),
            policy,
            legacy_trace_bytes(
                AuditDepletionKind2::SEGMENT,
                motion.digest(),
                policy,
                construction.trace),
            construction.trace.strategy_version,
            witness);
    } catch (const ExactDepletionConstructionError&) {
        return false;
    }
}

bool validate_audit_circle_depletion_witness(
    const Stock2& authority,
    const Stock2& post,
    const AuditCircleMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDepletionWitness2& witness)
{
    const AuditTrialPrimitiveScope2 scope;
    try {
        const ExactDepletionConstruction2 construction =
            circle_construction(motion, policy);
        Gps removal = build_exact_disk_union_region_2(
            construction.centers, policy.tool_radius_mm());
        return expected_post_matches(authority, post, std::move(removal))
            && AuditDepletionWitnessFactory2::matches(
            AuditDepletionKind2::FULL_CIRCLE,
            authority,
            post,
            motion.digest(),
            policy,
            legacy_trace_bytes(
                AuditDepletionKind2::FULL_CIRCLE,
                motion.digest(),
                policy,
                construction.trace),
            construction.trace.strategy_version,
            witness);
    } catch (const ExactDepletionConstructionError&) {
        return false;
    }
}

bool validate_audit_arc_depletion_witness(
    const Stock2& authority,
    const Stock2& post,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDepletionWitness2& witness)
{
    const AuditTrialPrimitiveScope2 scope;
    try {
        const ExactArcDepletionConstruction2 construction =
            arc_construction(motion, policy);
        Gps removal = build_exact_disk_union_region_2(
            construction.centers, policy.tool_radius_mm());
        return expected_post_matches(authority, post, std::move(removal))
            && AuditDepletionWitnessFactory2::matches(
            AuditDepletionKind2::ARC,
            authority,
            post,
            motion.digest(),
            policy,
            construction.trace.canonical_bytes(),
            construction.trace.strategy_version(),
            witness);
    } catch (const ExactDepletionConstructionError&) {
        return false;
    }
}

bool validate_audit_plunge_depletion_witness(
    const Stock2& authority,
    const Stock2& post,
    const AuditVerticalPlunge2& motion,
    const AuditPolicy2& policy,
    const AuditDepletionWitness2& witness)
{
    const AuditTrialPrimitiveScope2 scope;
    Gps removal = build_exact_disk_union_region_2(
        {motion.cut_endpoint()}, policy.tool_radius_mm());
    return expected_post_matches(authority, post, std::move(removal))
        && AuditDepletionWitnessFactory2::matches(
        AuditDepletionKind2::PLUNGE,
        authority,
        post,
        motion.digest(),
        policy,
        plunge_trace_bytes(motion, policy),
        "exact-plunge-disk-v1",
        witness);
}

void reset_audit_trial_depletion_instrumentation_for_test() noexcept
{
    trial_metrics = {0, 0};
}

void note_audit_trial_stock_clone_for_test() noexcept
{
    if (trial_primitive_scope_depth != 0) {
        ++trial_metrics.internal_clone_count;
    }
}

void note_audit_trial_stock_swap_for_test() noexcept
{
    if (trial_primitive_scope_depth != 0) {
        ++trial_metrics.internal_swap_count;
    }
}

AuditTrialDepletionInstrumentation2
audit_trial_depletion_instrumentation_for_test() noexcept
{
    return trial_metrics;
}

void reset_legacy_arc_sweep_reachability_for_test() noexcept
{
    legacy_arc_sweep_count = 0;
}

std::size_t legacy_arc_sweep_reachability_for_test() noexcept
{
    return legacy_arc_sweep_count;
}

void note_legacy_arc_sweep_reached_for_test() noexcept
{
    ++legacy_arc_sweep_count;
}
