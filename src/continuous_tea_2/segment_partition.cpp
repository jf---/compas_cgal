#include "segment_partition.h"

#include "../exact_algebraic_1.h"
#include "event_certificate.h"
#include "event_partition.h"
#include "segment_projection.h"
#include "segment_fibre.h"
#include "segment_strata.h"
#include "sha256.h"

#include <algorithm>
#include <optional>
#include <string_view>
#include <utility>

namespace {

std::string boolean_text(bool value)
{
    return value ? "1" : "0";
}

void normalize_tangent_root_multiplicities(
    EventPartitionCertificate2& certificate,
    const std::vector<ProjectionInput2>& projections)
{
    std::vector<ProjectionInput2> tangent_projections;
    for (const ProjectionInput2& projection :
         projections) {
        if (std::any_of(
                projection.events.begin(),
                projection.events.end(),
                [](const PartitionEvent2& event) {
                    return event.kind == "tangent";
                })) {
            tangent_projections.push_back(projection);
        }
    }
    if (tangent_projections.empty()) {
        return;
    }
    const EventPartitionCertificate2 tangent_partition =
        partition_projections(tangent_projections);
    for (AlgebraicRootRecord2& root :
         certificate.roots) {
        const auto tangent_root = std::find_if(
            tangent_partition.roots.begin(),
            tangent_partition.roots.end(),
            [&root](
                const AlgebraicRootRecord2& candidate) {
                return candidate.root_id
                    == root.root_id;
            });
        if (tangent_root
            != tangent_partition.roots.end()) {
            root.multiplicity =
                tangent_root->multiplicity;
        }
    }
}

std::string canonical_event(
    const PartitionEvent2& event)
{
    return encode_string_sequence(
        {
            event.kind,
            event.feature_id,
            event.support_id,
            event.trim_id,
            event.vertex_id,
            event.branch_id,
            event.disposition,
            std::to_string(event.left_active_count),
            std::to_string(event.right_active_count),
            boolean_text(
                event.incidence_permutation_rechecked),
            boolean_text(
                event.original_equations_rechecked),
            boolean_text(event.orientation_rechecked),
            event.trim_disposition,
            event.pair_sheet_id,
            event.first_feature_id,
            event.second_feature_id,
            event.first_chart_id,
            event.second_chart_id,
            event.first_branch_id,
            event.second_branch_id,
        });
}

std::string canonical_branch(
    const SegmentBoundaryBranch2& branch)
{
    return encode_string_sequence(
        {
            branch.branch_id,
            branch.feature_id,
            branch.support_id,
            branch.support_kind,
            branch.trim_id,
            branch.vertex_id,
            branch.material_side,
            branch.trim_disposition,
            branch.rim_chart_id,
            std::to_string(branch.rim_sheet_ordinal),
            branch.rim_root_id,
            encode_string_sequence(
                branch.rim_factor_coefficients),
            std::to_string(branch.rim_root_ordinal),
        });
}

std::string canonical_pair(
    const BranchPairDisposition2& pair)
{
    return encode_string_sequence(
        {
            pair.pair_sheet_id,
            pair.first_branch_id,
            pair.second_branch_id,
            pair.orientation_disposition,
            pair.cap_disposition,
        });
}

std::string canonical_pairs(
    const std::vector<BranchPairDisposition2>& values)
{
    std::vector<std::string> pairs;
    pairs.reserve(values.size());
    for (const BranchPairDisposition2& pair : values) {
        pairs.push_back(canonical_pair(pair));
    }
    return encode_string_sequence(pairs);
}

std::string canonical_stratum(
    const SegmentEventStratum2& stratum)
{
    std::vector<std::string> events;
    events.reserve(stratum.events.size());
    for (const PartitionEvent2& event :
         stratum.events) {
        events.push_back(canonical_event(event));
    }
    return encode_string_sequence(
        {
            stratum.kind,
            stratum.root_id,
            stratum.local_root_id,
            stratum.global_fibre_id,
            stratum.chart_id,
            stratum.witness_numerator,
            stratum.witness_denominator,
            encode_string_sequence(
                stratum.root_factor_coefficients),
            std::to_string(stratum.root_ordinal),
            encode_string_sequence(
                stratum.active_branch_ids),
            encode_string_sequence(events),
            canonical_pairs(
                stratum.left_pair_dispositions),
            canonical_pairs(
                stratum.pair_dispositions),
            canonical_pairs(
                stratum.right_pair_dispositions),
            boolean_text(
                stratum.algebraic_root_evaluated),
            boolean_text(
                stratum.original_equations_rechecked),
            boolean_text(
                stratum.orientation_rechecked),
            boolean_text(
                stratum.trim_predicates_rechecked),
        });
}

std::string canonical_segment_partition(
    const SegmentEventPartition2& partition)
{
    std::vector<std::string> projections;
    projections.reserve(partition.projections.size());
    for (const ProjectionRecord2& projection :
         partition.projections) {
        projections.push_back(
            encode_string_sequence(
                {
                    projection.projection_id,
                    projection.normalized_coefficient_bytes,
                    projection.degree_bound_id,
                }));
    }
    std::vector<std::string> branches;
    branches.reserve(partition.branches.size());
    for (const SegmentBoundaryBranch2& branch :
         partition.branches) {
        branches.push_back(canonical_branch(branch));
    }
    std::vector<std::string> strata;
    strata.reserve(partition.strata.size());
    for (const SegmentEventStratum2& stratum :
         partition.strata) {
        strata.push_back(canonical_stratum(stratum));
    }
    return encode_string_sequence(
        {
            "segment-event-partition-v1",
            partition.source.canonical_bytes(),
            encode_string_sequence(
                partition.boundary_feature_ids),
            encode_string_sequence(projections),
            encode_string_sequence(branches),
            encode_string_sequence(strata),
            partition.certificate.canonical_bytes,
        });
}

std::vector<std::string> branch_ids(
    const std::vector<SegmentBranchState2>& states)
{
    std::vector<std::string> result;
    result.reserve(states.size());
    for (const SegmentBranchState2& state : states) {
        if (std::find(
                result.begin(),
                result.end(),
                state.branch.branch_id)
            == result.end()) {
            result.push_back(state.branch.branch_id);
        }
    }
    return result;
}

std::vector<std::string> ordered_union_ids(
    const std::vector<std::string>& left,
    const std::vector<std::string>& right)
{
    std::vector<std::string> result = left;
    for (const std::string& branch_id : right) {
        if (std::find(
                result.begin(),
                result.end(),
                branch_id)
            == result.end()) {
            result.push_back(branch_id);
        }
    }
    return result;
}

std::vector<BranchPairDisposition2> union_pairs(
    const std::vector<BranchPairDisposition2>& left,
    const std::vector<BranchPairDisposition2>& right)
{
    std::vector<BranchPairDisposition2> result = left;
    for (const BranchPairDisposition2& pair : right) {
        const auto found = std::find_if(
            result.begin(),
            result.end(),
            [&pair](const BranchPairDisposition2& candidate) {
                    return candidate.pair_sheet_id
                        == pair.pair_sheet_id;
            });
        if (found == result.end()) {
            result.push_back(pair);
        }
    }
    return result;
}

std::vector<ChartSeam2> owned_seams(
    const std::vector<ParameterChart2>& charts)
{
    std::vector<ChartSeam2> seams;
    seams.reserve(charts.size());
    for (const ParameterChart2& chart : charts) {
        seams.push_back(
            {
                chart.start_seam_id,
                chart.chart_id,
            });
    }
    return seams;
}

std::optional<std::string> projection_event_class(
    const ProjectionRecord2& projection)
{
    std::vector<std::string> fields;
    try {
        fields = decode_string_sequence(
            projection.projection_id);
    } catch (const EventSubstrateError&) {
        return std::nullopt;
    }
    if (fields.empty()) {
        return std::nullopt;
    }
    if (fields.front() == "segment-tangency-v1") {
        return "support-tangency";
    }
    if (fields.front()
        == "segment-vertex-passage-v1") {
        return "vertex-passage";
    }
    if (fields.front()
        == "segment-support-overlap-v1") {
        return "support-overlap";
    }
    if (fields.front()
        == "segment-cap-crossing-v1") {
        return "cap-crossing";
    }
    if (fields.front()
        == "segment-endpoint-order-v1") {
        return "endpoint-order-resultant";
    }
    if (fields.front()
        == "segment-pair-orientation-v1") {
        return "pair-orientation-boundary";
    }
    if (fields.front() == "segment-pair-cap-v1") {
        return "cross-feature-cap-resultant";
    }
    if (fields.front()
        == "segment-pair-endpoint-v1") {
        return "cross-feature-endpoint-resultant";
    }
    return std::nullopt;
}

void finalize_segment_partition(
    SegmentEventPartition2& partition)
{
    finalize_event_partition(
        partition.certificate);
    partition.canonical_bytes =
        canonical_segment_partition(partition);
    partition.canonical_digest =
        sha256_bytes(partition.canonical_bytes);
}

std::optional<std::size_t> cell_with_upper_root(
    const EventPartitionCertificate2& certificate,
    const std::string& root_id)
{
    for (std::size_t index = 0;
         index < certificate.cells.size();
         ++index) {
        if (certificate.cells[index].upper_root_id
            == root_id) {
            return index;
        }
    }
    return std::nullopt;
}

std::optional<std::size_t> cell_with_lower_root(
    const EventPartitionCertificate2& certificate,
    const std::string& root_id)
{
    for (std::size_t index = 0;
         index < certificate.cells.size();
         ++index) {
        if (certificate.cells[index].lower_root_id
            == root_id) {
            return index;
        }
    }
    return std::nullopt;
}

} // namespace

SegmentEventPartition2 construct_segment_event_partition(
    const Stock2& stock,
    const SegmentEventSource2& source)
{
    const std::vector<BoundaryFeatureRecord2> records =
        extract_boundary_records(stock);
    if (records.empty()) {
        throw IncompleteSegmentPartitionError(
            "segment partition requires a nonempty stock boundary");
    }

    std::vector<std::string> feature_ids;
    feature_ids.reserve(records.size());
    for (const BoundaryFeatureRecord2& record : records) {
        feature_ids.push_back(record.feature_id);
    }
    if (!std::is_sorted(
            feature_ids.begin(),
            feature_ids.end())) {
        throw IncompleteSegmentPartitionError(
            "Task 4 boundary records are not canonically ordered");
    }

    SegmentProjectionBundle2 projections =
        derive_segment_projections(records, source);
    if (projections.event_projections.empty()) {
        throw IncompleteSegmentPartitionError(
            "segment partition has no exact event projections");
    }
    EventPartitionCertificate2 certificate =
        partition_projections(
            projections.event_projections);
    normalize_tangent_root_multiplicities(
        certificate,
        projections.event_projections);
    certificate.projections.insert(
        certificate.projections.end(),
        projections.constant_event_projections.begin(),
        projections.constant_event_projections.end());
    certificate.overlaps.insert(
        certificate.overlaps.end(),
        projections.overlaps.begin(),
        projections.overlaps.end());
    certificate.projections.insert(
        certificate.projections.begin(),
        projections.pullbacks.begin(),
        projections.pullbacks.end());
    certificate.seams =
        owned_seams(certificate.charts);
    certificate.source_kind =
        "segment-event-substrate-v1";
    certificate.source_payload =
        encode_string_sequence(
            {
                source.canonical_bytes(),
                encode_string_sequence(feature_ids),
            });

    std::vector<std::vector<SegmentBranchState2>>
        cell_states;
    std::vector<
        std::vector<BranchPairDisposition2>>
        cell_pairs;
    cell_states.reserve(certificate.cells.size());
    cell_pairs.reserve(certificate.cells.size());
    for (const ParameterCell2& cell :
         certificate.cells) {
        cell_states.push_back(
            segment_branches_at(
                records,
                source,
                cell.witness_numerator,
                cell.witness_denominator));
        cell_pairs.push_back(
            segment_branch_pair_dispositions(
                cell_states.back(),
                source));
    }
    for (ParameterCell2& cell : certificate.cells) {
        cell.disposition = "sign-invariant";
    }

    std::vector<SegmentBoundaryBranch2> branches;
    for (const auto& states : cell_states) {
        for (const SegmentBranchState2& state : states) {
            const auto found = std::find_if(
                branches.begin(),
                branches.end(),
                [&state](
                    const SegmentBoundaryBranch2& candidate) {
                    return candidate.branch_id
                        == state.branch.branch_id;
                });
            if (found == branches.end()) {
                branches.push_back(state.branch);
            }
        }
    }
    std::sort(
        branches.begin(),
        branches.end(),
        [](const SegmentBoundaryBranch2& first,
           const SegmentBoundaryBranch2& second) {
            return first.branch_id < second.branch_id;
        });

    std::vector<SegmentFibreEvaluation2>
        fibre_evaluations;
    fibre_evaluations.reserve(
        certificate.fibres.size());
    for (std::size_t fibre_index = 0;
         fibre_index < certificate.fibres.size();
         ++fibre_index) {
        const EventFibre2& fibre =
            certificate.fibres[fibre_index];
        const std::optional<std::size_t> left_index =
            cell_with_upper_root(
                certificate,
                fibre.root_id);
        const std::optional<std::size_t> right_index =
            cell_with_lower_root(
                certificate,
                fibre.root_id);
        const std::vector<SegmentBranchState2>
            no_states;
        SegmentFibreEvaluation2 evaluation =
            evaluate_segment_fibre(
                records,
                source,
                projections.pullbacks,
                projections.event_projections,
                certificate.roots[fibre_index],
                fibre.events,
                left_index.has_value()
                    ? cell_states[*left_index]
                    : no_states,
                right_index.has_value()
                    ? cell_states[*right_index]
                    : no_states);
        const bool incidence_rechecked =
            std::all_of(
                evaluation.events.begin(),
                evaluation.events.end(),
                [](const PartitionEvent2& event) {
                    return event.kind
                            != "endpoint-order"
                        || event
                               .incidence_permutation_rechecked;
                });
        const bool has_support_overlap =
            std::any_of(
                evaluation.events.begin(),
                evaluation.events.end(),
                [](const PartitionEvent2& event) {
                    return event.kind
                        == "support-overlap";
                });
        const std::size_t active_count =
            evaluation.active_branch_ids.size();
        const bool pairs_complete =
            active_count < 2
            || evaluation.pair_dispositions.size()
                == active_count * (active_count - 1)
            || has_support_overlap;
        if (!evaluation
                    .original_equations_rechecked
            || !evaluation.orientation_rechecked
            || !evaluation.trim_predicates_rechecked
            || !incidence_rechecked
            || !pairs_complete) {
            throw IncompleteSegmentPartitionError(
                "segment fibre replay did not prove complete exact evidence: events="
                + std::to_string(
                    evaluation.events.size())
                + "/"
                + std::to_string(fibre.events.size())
                + ", original="
                + boolean_text(
                    evaluation
                        .original_equations_rechecked)
                + ", orientation="
                + boolean_text(
                    evaluation.orientation_rechecked)
                + ", trim="
                + boolean_text(
                    evaluation
                        .trim_predicates_rechecked)
                + ", incidence="
                + boolean_text(incidence_rechecked)
                + ", active="
                + std::to_string(active_count)
                + ", pairs="
                + std::to_string(
                    evaluation
                        .pair_dispositions.size())
                + ", kind="
                + (
                    evaluation.events.empty()
                        ? std::string("none")
                        : evaluation.events.front().kind));
        }
        certificate.fibres[fibre_index].events =
            evaluation.events;
        for (const SegmentBranchState2& state :
             evaluation.branches) {
            const auto found = std::find_if(
                branches.begin(),
                branches.end(),
                [&state](
                    const SegmentBoundaryBranch2& candidate) {
                    return candidate.branch_id
                        == state.branch.branch_id;
                });
            if (found == branches.end()) {
                branches.push_back(state.branch);
            }
        }
        fibre_evaluations.push_back(
            std::move(evaluation));
    }
    std::sort(
        branches.begin(),
        branches.end(),
        [](const SegmentBoundaryBranch2& first,
           const SegmentBoundaryBranch2& second) {
            return first.branch_id < second.branch_id;
        });
    finalize_event_partition(certificate);

    const auto cell_stratum =
        [&certificate, &cell_states, &cell_pairs](
            std::size_t cell_index) {
        const ParameterCell2& cell =
            certificate.cells[cell_index];
        return SegmentEventStratum2{
            "cell",
            {},
            {},
            {},
            "segment-linear-v1",
            cell.witness_numerator,
            cell.witness_denominator,
            {},
            0,
            branch_ids(cell_states[cell_index]),
            {},
            {},
            cell_pairs[cell_index],
            {},
            false,
            false,
            false,
            false,
        };
    };

    std::vector<SegmentEventStratum2> strata;
    std::vector<bool> emitted_cells(
        certificate.cells.size(),
        false);
    for (std::size_t fibre_index = 0;
         fibre_index < certificate.fibres.size();
         ++fibre_index) {
        const EventFibre2& fibre =
            certificate.fibres[fibre_index];
        const std::optional<std::size_t> left_index =
            cell_with_upper_root(
                certificate,
                fibre.root_id);
        const std::optional<std::size_t> right_index =
            cell_with_lower_root(
                certificate,
                fibre.root_id);
        if (left_index.has_value()
            && !emitted_cells[*left_index]) {
            strata.push_back(
                cell_stratum(*left_index));
            emitted_cells[*left_index] = true;
        }
        const std::vector<std::string> before =
            left_index.has_value()
            ? branch_ids(cell_states[*left_index])
            : std::vector<std::string>{};
        const std::vector<std::string> after =
            right_index.has_value()
            ? branch_ids(cell_states[*right_index])
            : std::vector<std::string>{};
        const std::vector<BranchPairDisposition2>
            left_pairs =
                left_index.has_value()
                ? cell_pairs[*left_index]
                : std::vector<
                      BranchPairDisposition2>{};
        const std::vector<BranchPairDisposition2>
            right_pairs =
                right_index.has_value()
                ? cell_pairs[*right_index]
                : std::vector<
                      BranchPairDisposition2>{};
        const SegmentFibreEvaluation2& evaluation =
            fibre_evaluations[fibre_index];
        strata.push_back(
            {
                "fibre",
                fibre.root_id,
                fibre.root_id,
                encode_string_sequence(
                    {
                        "segment-global-fibre-v1",
                        source.canonical_bytes(),
                        fibre.root_id,
                    }),
                "segment-linear-v1",
                {},
                {},
                certificate.roots[fibre_index]
                    .factor_coefficients,
                certificate.roots[fibre_index]
                    .root_ordinal,
                evaluation.active_branch_ids,
                fibre.events,
                left_pairs,
                evaluation.pair_dispositions,
                right_pairs,
                evaluation.algebraic_root_evaluated,
                evaluation
                    .original_equations_rechecked,
                evaluation.orientation_rechecked,
                evaluation.trim_predicates_rechecked,
            });
        if (right_index.has_value()
            && !emitted_cells[*right_index]) {
            strata.push_back(
                cell_stratum(*right_index));
            emitted_cells[*right_index] = true;
        }
    }
    for (std::size_t cell_index = 0;
         cell_index < certificate.cells.size();
         ++cell_index) {
        if (!emitted_cells[cell_index]) {
            strata.push_back(
                cell_stratum(cell_index));
        }
    }

    SegmentEventPartition2 result{
        source,
        std::move(feature_ids),
        projections.pullbacks,
        std::move(branches),
        std::move(strata),
        std::move(certificate),
        {},
        {},
    };
    result.canonical_bytes =
        canonical_segment_partition(result);
    result.canonical_digest =
        sha256_bytes(result.canonical_bytes);
    return result;
}

SegmentEventPartition2 construct_segment_event_partition(
    const Stock2& stock,
    double x0,
    double y0,
    double x1,
    double y1,
    double tool_radius,
    double cap_chord_ratio)
{
    return construct_segment_event_partition(
        stock,
        SegmentEventSource2::from_binary64(
            x0,
            y0,
            x1,
            y1,
            tool_radius,
            cap_chord_ratio));
}

VerifiedSegmentEventPartition2
verify_segment_event_partition(
    const Stock2& stock,
    const SegmentEventSource2& source,
    const SegmentEventPartition2& candidate)
{
    try {
        const SegmentEventPartition2 expected =
            construct_segment_event_partition(
                stock,
                source);
        const std::string candidate_bytes =
            canonical_segment_partition(candidate);
        EventPartitionCertificate2
            candidate_certificate =
                candidate.certificate;
        const std::string claimed_certificate_bytes =
            candidate_certificate.canonical_bytes;
        const std::string claimed_certificate_digest =
            candidate_certificate.canonical_digest;
        finalize_event_partition(
            candidate_certificate);
        const bool nested_certificate_valid =
            claimed_certificate_bytes
                == candidate_certificate
                       .canonical_bytes
            && claimed_certificate_digest
                == candidate_certificate
                       .canonical_digest;
        const bool valid =
            nested_certificate_valid
            && candidate.source.canonical_bytes()
                == source.canonical_bytes()
            && candidate.canonical_bytes
                == candidate_bytes
            && candidate.canonical_digest
                == sha256_bytes(candidate_bytes)
            && candidate.canonical_bytes
                == expected.canonical_bytes
            && candidate.canonical_digest
                == expected.canonical_digest;
        return {
            valid
                ? ContinuousTeaVerdict::CERTIFIED
                : ContinuousTeaVerdict::
                      UNRESOLVED_DEGENERACY,
            candidate,
        };
    } catch (const EventSubstrateError&) {
        return {
            ContinuousTeaVerdict::
                UNRESOLVED_DEGENERACY,
            candidate,
        };
    }
}

SegmentEventPartition2 mutate_segment_event_partition(
    const SegmentEventPartition2& partition,
    const std::string& mutation)
{
    SegmentEventPartition2 result = partition;
    if (mutation == "delete-pair-disposition") {
        const auto found = std::find_if(
            result.strata.begin(),
            result.strata.end(),
            [](const SegmentEventStratum2& stratum) {
                return !stratum.pair_dispositions.empty();
            });
        if (found != result.strata.end()) {
            found->pair_dispositions.pop_back();
        }
        finalize_segment_partition(result);
        return result;
    }
    constexpr std::string_view projection_prefix =
        "delete-projection-";
    if (mutation.starts_with(projection_prefix)) {
        const std::string event_class =
            mutation.substr(
                projection_prefix.size());
        const auto new_end = std::remove_if(
            result.certificate.projections.begin(),
            result.certificate.projections.end(),
            [&event_class](
                const ProjectionRecord2& projection) {
                return projection_event_class(
                           projection)
                    == event_class;
            });
        if (new_end
            == result.certificate.projections.end()) {
            throw EventPartitionVerificationError(
                "segment projection class is absent");
        }
        result.certificate.projections.erase(
            new_end,
            result.certificate.projections.end());
        finalize_segment_partition(result);
        return result;
    }
    throw EventPartitionVerificationError(
        "unknown segment partition mutation");
}
