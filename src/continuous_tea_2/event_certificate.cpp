#include "event_certificate.h"

#include "../exact_algebraic_1.h"
#include "cap_partition.h"
#include "circle_source.h"
#include "sha256.h"

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string_view>
#include <tuple>

namespace {

std::string u64_record(std::size_t value)
{
    std::string result;
    for (int shift = 56; shift >= 0; shift -= 8) {
        result.push_back(
            static_cast<char>(
                (static_cast<std::uint64_t>(value) >> shift)
                & 0xffU));
    }
    return result;
}

std::size_t read_u64(
    const std::string& encoded,
    std::size_t& offset)
{
    if (encoded.size() - offset < 8) {
        throw EventPartitionVerificationError(
            "truncated canonical sequence length");
    }
    std::uint64_t value = 0;
    for (std::size_t index = 0; index < 8; ++index) {
        value = (value << 8)
            | static_cast<unsigned char>(
                encoded[offset + index]);
    }
    offset += 8;
    return static_cast<std::size_t>(value);
}

std::string record(
    std::string_view tag,
    const std::vector<std::string>& fields)
{
    return std::string(tag) + '\0'
        + encode_string_sequence(fields);
}

std::string boolean_text(bool value)
{
    return value ? "1" : "0";
}

std::string canonical_backend(
    const AlgebraicBackendEvidence2& evidence)
{
    return record(
        "algebraic-backend-evidence-v1",
        {
            evidence.cgal_version,
            evidence.integer_backend,
            evidence.algebraic_kernel_1,
            evidence.algebraic_kernel_2,
            evidence.arrangement_traits,
            encode_string_sequence(
                evidence.compile_definitions),
        });
}

std::string canonical_chart(const ParameterChart2& chart)
{
    return record(
        "parameter-chart-v1",
        {
            chart.chart_id,
            chart.family,
            chart.domain_low,
            chart.domain_high,
            chart.orientation,
            chart.start_seam_id,
            chart.end_seam_id,
            boolean_text(chart.owns_start_seam),
            boolean_text(chart.owns_end_seam),
        });
}

std::string canonical_projection(
    const ProjectionRecord2& projection)
{
    std::vector<std::string> rows;
    for (const auto& row : projection.coefficient_rows) {
        rows.push_back(encode_string_sequence(row));
    }
    std::vector<std::string> factors;
    for (const auto& factor :
         projection.factor_coefficients) {
        factors.push_back(
            encode_string_sequence(factor));
    }
    std::vector<std::string> fields{
        projection.projection_id,
        encode_string_sequence(rows),
        encode_string_sequence(factors),
        std::to_string(
            projection.actual_motion_degree),
        std::to_string(
            projection.actual_rim_degree),
        std::to_string(
            projection.bound_motion_degree),
        std::to_string(
            projection.bound_rim_degree),
        projection.degree_bound_id,
        projection.normalized_coefficient_bytes,
    };
    if (!projection
             .signed_predicate_coefficients.empty()) {
        fields.push_back(
            encode_string_sequence(
                projection
                    .signed_predicate_coefficients));
    }
    return record(
        projection.signed_predicate_coefficients
                .empty()
            ? "projection-record-v1"
            : "signed-projection-record-v2",
        fields);
}

std::string canonical_root(
    const AlgebraicRootRecord2& root)
{
    return record(
        "algebraic-root-record-v1",
        {
            root.root_id,
            encode_string_sequence(
                root.factor_coefficients),
            std::to_string(root.root_ordinal),
            std::to_string(root.multiplicity),
        });
}

std::string canonical_cell(const ParameterCell2& cell)
{
    return record(
        "parameter-cell-v1",
        {
            cell.lower_root_id,
            cell.upper_root_id,
            cell.witness_numerator,
            cell.witness_denominator,
            cell.disposition,
        });
}

std::string canonical_event(const PartitionEvent2& event)
{
    return record(
        "partition-event-v2",
        {
            event.kind,
            event.feature_id,
            event.support_id,
            event.trim_id,
            event.vertex_id,
            event.branch_id,
            event.disposition,
            event.endpoint_role,
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

std::string canonical_physical_incidence(
    const PhysicalIncidence2& incidence)
{
    return record(
        "physical-incidence-v1",
        {
            incidence.kind,
            incidence.feature_id,
            incidence.support_id,
            incidence.trim_id,
            incidence.vertex_id,
            incidence.endpoint_role,
            std::to_string(incidence.sheet_ordinal),
        });
}

std::string canonical_fibre(const EventFibre2& fibre)
{
    std::vector<std::string> events;
    for (const PartitionEvent2& event : fibre.events) {
        events.push_back(canonical_event(event));
    }
    if (fibre.left_active_branches.empty()
        && fibre.right_active_branches.empty()
        && fibre.local_root_ids.empty()
        && fibre.seam_id.empty()
        && fibre.local_event_witnesses.empty()
        && fibre.ccw_direction.empty()
        && fibre.cw_direction.empty()) {
        return record(
            "event-fibre-v1",
            {
                fibre.root_id,
                encode_string_sequence(events),
            });
    }
    const auto canonical_branches =
        [](const std::vector<ActiveBoundaryBranch2>&
               branches) {
            std::vector<std::string> values;
            values.reserve(branches.size());
            for (const ActiveBoundaryBranch2& branch :
                 branches) {
                values.push_back(
                    record(
                        "active-boundary-sheet-v3",
                        {
                            branch.branch_id,
                            branch.feature_id,
                            branch.support_id,
                            branch.trim_id,
                            branch.chart_id,
                            std::to_string(
                                branch.sheet_ordinal),
                            branch.root_id,
                            branch.physical_incidence
                                .incidence_id,
                            canonical_physical_incidence(
                                branch
                                    .physical_incidence),
                        }));
            }
            return encode_string_sequence(values);
        };
    std::vector<std::string> witnesses;
    witnesses.reserve(
        fibre.local_event_witnesses.size());
    for (const LocalEventWitness2& witness :
         fibre.local_event_witnesses) {
        witnesses.push_back(
            record(
                "local-event-witness-v2",
                {
                    witness.kind,
                    witness.feature_id,
                    witness.support_id,
                    witness.trim_id,
                    witness.vertex_id,
                    witness.endpoint_role,
                    witness.disposition,
                    witness.local_branch_id,
                    witness.local_root_id,
                    witness.source_projection_id,
                    witness.source_factor_id,
                    std::to_string(
                        witness.multiplicity),
                }));
    }
    std::vector<std::string> fields{
        fibre.root_id,
        encode_string_sequence(events),
        encode_string_sequence(witnesses),
        canonical_branches(
            fibre.left_active_branches),
        canonical_branches(
            fibre.right_active_branches),
        fibre.ccw_direction,
        fibre.cw_direction,
    };
    if (!fibre.local_root_ids.empty()) {
        fields.insert(
            fields.begin() + 1,
            encode_string_sequence(
                fibre.local_root_ids));
        fields.insert(
            fields.begin() + 2,
            fibre.seam_id);
    }
    return record(
        fibre.local_root_ids.empty()
            ? "evidenced-event-fibre-v3"
            : "cyclic-seam-event-fibre-v4",
        fields);
}

void validate_local_event_witnesses(
    const EventFibre2& fibre)
{
    for (const LocalEventWitness2& witness :
         fibre.local_event_witnesses) {
        if (witness.source_projection_id.empty()
            || witness.source_factor_id.empty()
            || witness.local_root_id.empty()
            || witness.multiplicity == 0) {
            throw EventPartitionVerificationError(
                "local event witness lacks exact source provenance");
        }
    }
    for (std::size_t first = 0;
         first < fibre.local_event_witnesses.size();
         ++first) {
        for (std::size_t second = first + 1;
             second
             < fibre.local_event_witnesses.size();
             ++second) {
            const LocalEventWitness2& left =
                fibre.local_event_witnesses[first];
            const LocalEventWitness2& right =
                fibre.local_event_witnesses[second];
            const auto physical_key =
                [](const LocalEventWitness2& witness) {
                    return std::tie(
                        witness.kind,
                        witness.feature_id,
                        witness.support_id,
                        witness.trim_id,
                        witness.vertex_id,
                        witness.endpoint_role);
                };
            if (physical_key(left)
                    == physical_key(right)
                && (left.multiplicity
                        != right.multiplicity
                    || left.disposition
                        != right.disposition)) {
                throw EventPartitionVerificationError(
                    "physical event has contradictory local branch evidence");
            }
        }
    }
}

void validate_phase_seam_incidence(
    const EventFibre2& fibre)
{
    if (fibre.local_root_ids.empty()) {
        return;
    }
    const auto validate_incidence =
        [](const PhysicalIncidence2& incidence) {
            if (incidence.incidence_id.empty()
                || incidence.incidence_id
                    != canonical_physical_incidence(
                        incidence)
                || incidence.vertex_id.empty()
                || incidence.endpoint_role.empty()) {
                throw EventPartitionVerificationError(
                    "phase-seam physical incidence identity is invalid");
            }
        };
    for (const ActiveBoundaryBranch2& branch :
         fibre.left_active_branches) {
        validate_incidence(
            branch.physical_incidence);
    }
    for (const ActiveBoundaryBranch2& branch :
         fibre.right_active_branches) {
        validate_incidence(
            branch.physical_incidence);
    }
}

void validate_event_evidence(
    const EventPartitionCertificate2& certificate)
{
    for (const EventFibre2& fibre :
         certificate.fibres) {
        validate_local_event_witnesses(fibre);
        validate_phase_seam_incidence(fibre);
    }
}

std::string canonical_overlap(
    const OverlapInterval2& overlap)
{
    if (overlap.feature_id.empty()
        && overlap.support_id.empty()
        && overlap.trim_id.empty()
        && overlap.branch_id.empty()) {
        return record(
            "overlap-interval-v1",
            {
                overlap.kind,
                overlap.domain_low,
                overlap.domain_high,
                overlap.witness_numerator,
                overlap.witness_denominator,
                overlap.orientation_disposition,
            });
    }
    return record(
        "owned-overlap-interval-v2",
        {
            overlap.kind,
            overlap.domain_low,
            overlap.domain_high,
            overlap.witness_numerator,
            overlap.witness_denominator,
            overlap.orientation_disposition,
            overlap.feature_id,
            overlap.support_id,
            overlap.trim_id,
            overlap.branch_id,
        });
}

std::string canonical_seam(const ChartSeam2& seam)
{
    return record(
        "chart-seam-v1",
        {
            seam.seam_id,
            seam.owner_chart_id,
        });
}

template <typename Value, typename Encoder>
std::string canonical_vector(
    const std::vector<Value>& values,
    Encoder encoder)
{
    std::vector<std::string> encoded;
    encoded.reserve(values.size());
    for (const Value& value : values) {
        encoded.push_back(encoder(value));
    }
    return encode_string_sequence(encoded);
}

std::string canonical_certificate(
    const EventPartitionCertificate2& certificate)
{
    return record(
        "event-partition-certificate-v1",
        {
            canonical_backend(
                certificate.build_evidence),
            canonical_vector(
                certificate.charts,
                canonical_chart),
            canonical_vector(
                certificate.projections,
                canonical_projection),
            canonical_vector(
                certificate.roots,
                canonical_root),
            canonical_vector(
                certificate.cells,
                canonical_cell),
            canonical_vector(
                certificate.fibres,
                canonical_fibre),
            canonical_vector(
                certificate.overlaps,
                canonical_overlap),
            canonical_vector(
                certificate.seams,
                canonical_seam),
            certificate.source_kind,
            certificate.source_payload,
        });
}

PartitionEvent2 decode_event(
    const std::vector<std::string>& fields)
{
    if (fields.size() != 13) {
        throw EventPartitionVerificationError(
            "cap source payload has wrong field count");
    }
    PartitionEvent2 event(
        fields[6],
        fields[7],
        fields[8],
        fields[9],
        fields[10],
        fields[11],
        fields[12]);
    return event;
}

EventPartitionCertificate2 reconstruct(
    const EventPartitionCertificate2& certificate)
{
    if (certificate.source_kind
        == "full-circle-uniform-v1") {
        const std::vector<std::string> fields =
            decode_string_sequence(
                certificate.source_payload);
        if (fields.size() != 2) {
            throw EventPartitionVerificationError(
                "full-circle uniform source payload is malformed");
        }
        return construct_full_circle_uniform_partition(
            fields[0],
            fields[1]);
    }
    if (certificate.source_kind
        == "full-circle-boundary-pullbacks-v2") {
        const std::vector<std::string> fields =
            decode_string_sequence(
                certificate.source_payload);
        if (fields.size() != 6) {
            throw EventPartitionVerificationError(
                "full-circle line source payload is malformed");
        }
        return construct_full_circle_boundary_pullback_partition(
            fields[0],
            decode_string_sequence(fields[1]),
            fields[2],
            fields[3],
            decode_string_sequence(fields[4]),
            decode_string_sequence(fields[5]));
    }
    if (certificate.source_kind != "cap-crossings-v1") {
        throw EventPartitionVerificationError(
            "unsupported certificate source kind");
    }
    const std::vector<std::string> fields =
        decode_string_sequence(
            certificate.source_payload);
    const PartitionEvent2 event = decode_event(fields);
    return partition_cap_crossings(
        decode_string_sequence(fields[0]),
        decode_string_sequence(fields[1]),
        decode_string_sequence(fields[2]),
        decode_string_sequence(fields[3]),
        fields[4],
        fields[5],
        event);
}

} // namespace

std::string encode_string_sequence(
    const std::vector<std::string>& values)
{
    std::string result = u64_record(values.size());
    for (const std::string& value : values) {
        result += u64_record(value.size());
        result += value;
    }
    return result;
}

std::vector<std::string> decode_string_sequence(
    const std::string& encoded)
{
    std::size_t offset = 0;
    const std::size_t count = read_u64(encoded, offset);
    std::vector<std::string> result;
    result.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        const std::size_t size = read_u64(encoded, offset);
        if (size > encoded.size() - offset) {
            throw EventPartitionVerificationError(
                "truncated canonical sequence payload");
        }
        result.push_back(encoded.substr(offset, size));
        offset += size;
    }
    if (offset != encoded.size()) {
        throw EventPartitionVerificationError(
            "canonical sequence has trailing bytes");
    }
    return result;
}

std::string encode_canonical_record(
    const std::string& tag,
    const std::vector<std::string>& fields)
{
    return record(tag, fields);
}

void finalize_event_partition(
    EventPartitionCertificate2& certificate)
{
    certificate.canonical_bytes =
        canonical_certificate(certificate);
    certificate.canonical_digest =
        sha256_bytes(certificate.canonical_bytes);
}

VerifiedEventPartition2 verify_event_partition(
    const EventPartitionCertificate2& certificate)
{
    EventPartitionCertificate2 candidate = certificate;
    const std::string claimed_bytes =
        candidate.canonical_bytes;
    const std::string claimed_digest =
        candidate.canonical_digest;
    try {
        finalize_event_partition(candidate);
        validate_algebraic_root_intervals(
            candidate.roots);
        validate_event_evidence(candidate);
        const EventPartitionCertificate2 expected =
            reconstruct(candidate);
        const bool valid =
            claimed_bytes == candidate.canonical_bytes
            && claimed_digest == candidate.canonical_digest
            && candidate.canonical_bytes
                == expected.canonical_bytes
            && candidate.canonical_digest
                == expected.canonical_digest;
        return {
            valid
                ? ContinuousTeaVerdict::CERTIFIED
                : ContinuousTeaVerdict::
                      UNRESOLVED_DEGENERACY,
            std::move(candidate),
        };
    } catch (const EventSubstrateError&) {
        return {
            ContinuousTeaVerdict::UNRESOLVED_DEGENERACY,
            std::move(candidate),
        };
    }
}

EventPartitionCertificate2 mutate_certificate_record(
    const EventPartitionCertificate2& certificate,
    const std::string& mutation)
{
    EventPartitionCertificate2 result = certificate;
    if (mutation == "delete-seam" && !result.seams.empty()) {
        result.seams.pop_back();
    } else if (
        mutation == "delete-root" && !result.roots.empty()) {
        result.roots.pop_back();
    } else if (
        mutation == "delete-cell" && !result.cells.empty()) {
        result.cells.pop_back();
    } else if (
        mutation == "delete-fibre" && !result.fibres.empty()) {
        result.fibres.pop_back();
    } else if (
        mutation == "delete-factor"
        && !result.projections.empty()
        && !result.projections.front()
                .factor_coefficients.empty()) {
        result.projections.front()
            .factor_coefficients.pop_back();
    } else if (
        mutation == "coalesce-roots"
        && result.roots.size() >= 2) {
        result.roots[1].root_id =
            result.roots[0].root_id;
    } else if (
        mutation == "witness-outside-cell"
        && !result.cells.empty()) {
        result.cells.front().witness_numerator = "2";
        result.cells.front().witness_denominator = "1";
    } else if (
        mutation == "alter-multiplicity"
        && !result.roots.empty()) {
        ++result.roots.front().multiplicity;
    } else if (
        mutation == "alter-event-identity"
        && !result.fibres.empty()
        && !result.fibres.front().events.empty()) {
        result.fibres.front()
            .events.front()
            .feature_id.push_back('\0');
    } else if (
        mutation == "alter-disposition"
        && !result.cells.empty()) {
        result.cells.front().disposition =
            "mutated-disposition";
    } else if (
        mutation == "alter-fibre-direction") {
        const auto fibre = std::find_if(
            result.fibres.begin(),
            result.fibres.end(),
            [](const EventFibre2& candidate) {
                return !candidate
                            .ccw_direction.empty();
            });
        if (fibre == result.fibres.end()) {
            throw EventPartitionVerificationError(
                "certificate has no directed fibre");
        }
        fibre->ccw_direction =
            fibre->ccw_direction == "merge"
            ? "split"
            : "merge";
    } else if (
        mutation == "delete-signed-predicate") {
        const auto projection = std::find_if(
            result.projections.begin(),
            result.projections.end(),
            [](const ProjectionRecord2& candidate) {
                return !candidate
                            .signed_predicate_coefficients
                            .empty();
            });
        if (projection == result.projections.end()) {
            throw EventPartitionVerificationError(
                "certificate has no signed predicate");
        }
        projection->signed_predicate_coefficients
            .clear();
    } else if (
        mutation == "delete-active-sheet") {
        const auto fibre = std::find_if(
            result.fibres.begin(),
            result.fibres.end(),
            [](const EventFibre2& candidate) {
                return !candidate
                            .left_active_branches.empty()
                    || !candidate
                            .right_active_branches.empty();
            });
        if (fibre == result.fibres.end()) {
            throw EventPartitionVerificationError(
                "certificate has no active sheet");
        }
        std::vector<ActiveBoundaryBranch2>& active =
            !fibre->left_active_branches.empty()
            ? fibre->left_active_branches
            : fibre->right_active_branches;
        active.pop_back();
    } else if (
        mutation == "alter-degree"
        && !result.projections.empty()) {
        ++result.projections.front()
              .actual_motion_degree;
    } else if (
        mutation == "alter-bound"
        && !result.projections.empty()) {
        ++result.projections.front()
              .bound_motion_degree;
    } else if (
        mutation == "double-own-seam"
        && !result.charts.empty()) {
        result.charts.front().owns_end_seam = true;
    } else if (
        mutation == "alter-seam-incidence-vertex") {
        const auto fibre = std::find_if(
            result.fibres.begin(),
            result.fibres.end(),
            [](const EventFibre2& candidate) {
                return !candidate
                            .local_root_ids.empty()
                    && (!candidate
                             .left_active_branches
                             .empty()
                        || !candidate
                                .right_active_branches
                                .empty());
            });
        if (fibre == result.fibres.end()) {
            throw EventPartitionVerificationError(
                "certificate has no active seam incidence");
        }
        std::vector<ActiveBoundaryBranch2>& active =
            !fibre->left_active_branches.empty()
            ? fibre->left_active_branches
            : fibre->right_active_branches;
        active.front()
            .physical_incidence.vertex_id
            .push_back('\0');
    } else if (
        mutation == "alter-seam-incidence-trim") {
        const auto fibre = std::find_if(
            result.fibres.begin(),
            result.fibres.end(),
            [](const EventFibre2& candidate) {
                return !candidate
                            .local_root_ids.empty()
                    && (!candidate
                             .left_active_branches
                             .empty()
                        || !candidate
                                .right_active_branches
                                .empty());
            });
        if (fibre == result.fibres.end()) {
            throw EventPartitionVerificationError(
                "certificate has no active seam incidence");
        }
        std::vector<ActiveBoundaryBranch2>& active =
            !fibre->left_active_branches.empty()
            ? fibre->left_active_branches
            : fibre->right_active_branches;
        active.front()
            .physical_incidence.trim_id
            .push_back('\0');
    } else if (
        mutation == "alter-seam-local-multiplicity") {
        const auto fibre = std::find_if(
            result.fibres.begin(),
            result.fibres.end(),
            [](const EventFibre2& candidate) {
                for (std::size_t first = 0;
                     first
                     < candidate
                           .local_event_witnesses
                           .size();
                     ++first) {
                    for (std::size_t second =
                             first + 1;
                         second
                         < candidate
                               .local_event_witnesses
                               .size();
                         ++second) {
                        const auto physical_key =
                            [](const LocalEventWitness2&
                                   witness) {
                                return std::tie(
                                    witness.kind,
                                    witness.feature_id,
                                    witness.support_id,
                                    witness.trim_id,
                                    witness.vertex_id,
                                    witness
                                        .endpoint_role);
                            };
                        if (physical_key(
                                candidate
                                    .local_event_witnesses
                                        [first])
                            == physical_key(
                                candidate
                                    .local_event_witnesses
                                        [second])) {
                            return true;
                        }
                    }
                }
                return false;
            });
        if (fibre == result.fibres.end()) {
            throw EventPartitionVerificationError(
                "certificate has no coincident local event evidence");
        }
        bool mutated = false;
        for (std::size_t first = 0;
             !mutated
             && first
             < fibre->local_event_witnesses.size();
             ++first) {
            for (std::size_t second = first + 1;
                 second
                 < fibre->local_event_witnesses.size();
                 ++second) {
                const auto physical_key =
                    [](const LocalEventWitness2&
                           witness) {
                        return std::tie(
                            witness.kind,
                            witness.feature_id,
                            witness.support_id,
                            witness.trim_id,
                            witness.vertex_id,
                            witness.endpoint_role);
                    };
                if (physical_key(
                        fibre->local_event_witnesses
                            [first])
                    == physical_key(
                        fibre->local_event_witnesses
                            [second])) {
                    ++fibre->local_event_witnesses
                          [second]
                              .multiplicity;
                    mutated = true;
                    break;
                }
            }
        }
    } else {
        throw EventPartitionVerificationError(
            "unknown or inapplicable certificate mutation");
    }
    finalize_event_partition(result);
    return result;
}
