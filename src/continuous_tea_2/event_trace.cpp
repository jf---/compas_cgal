#include "event_trace.h"

#include "event_certificate.h"
#include "sha256.h"

#include <algorithm>
#include <string_view>
#include <utility>

namespace {

void canonicalize_ids(
    std::vector<std::string>& identifiers,
    std::string_view role)
{
    if (identifiers.empty()) {
        throw EventTraceVerificationError(
            std::string(role)
            + " must not be empty");
    }
    if (std::any_of(
            identifiers.begin(),
            identifiers.end(),
            [](const std::string& value) {
                return value.empty();
            })) {
        throw EventTraceVerificationError(
            std::string(role)
            + " contains an empty identity");
    }
    std::sort(
        identifiers.begin(),
        identifiers.end());
    identifiers.erase(
        std::unique(
            identifiers.begin(),
            identifiers.end()),
        identifiers.end());
}

std::string verdict_text(
    ContinuousTeaVerdict verdict)
{
    switch (verdict) {
    case ContinuousTeaVerdict::CERTIFIED:
        return "certified";
    case ContinuousTeaVerdict::CAP_EXCEEDED:
        return "cap_exceeded";
    case ContinuousTeaVerdict::UNRESOLVED_DEGENERACY:
        return "unresolved";
    }
    throw EventTraceVerificationError(
        "unknown continuous TEA verdict");
}

bool has_root(
    const EventPartitionCertificate2& partition,
    const std::string& root_id)
{
    return std::any_of(
        partition.roots.begin(),
        partition.roots.end(),
        [&root_id](
            const AlgebraicRootRecord2& root) {
            return root.root_id == root_id;
        });
}

bool trace_event_less(
    const EventTraceEvent2& left,
    const EventTraceEvent2& right)
{
    if (left.motion_order != right.motion_order) {
        return left.motion_order < right.motion_order;
    }
    if (left.global_fibre_id
        != right.global_fibre_id) {
        return left.global_fibre_id
            < right.global_fibre_id;
    }
    return left.canonical_id < right.canonical_id;
}

void validate_decision_authority(
    const std::string& authority_bytes,
    const std::string& authority_digest)
{
    if (authority_bytes.empty()) {
        throw EventTraceVerificationError(
            "event trace requires canonical decision-authority bytes");
    }
    if (authority_digest.empty()
        || sha256_bytes(authority_bytes)
            != authority_digest) {
        throw EventTraceVerificationError(
            "event trace decision-authority digest does not match its bytes");
    }
}

} // namespace

const std::string& event_oracle_component_version()
{
    static const std::string version =
        "event-exact-motion-oracle-v2";
    return version;
}

EventTraceEvent2 make_event_trace_event(
    const std::string& root_id,
    const std::string& global_fibre_id,
    const std::string& kind,
    std::vector<std::string> feature_ids,
    std::vector<std::string> branch_ids,
    unsigned int multiplicity,
    const std::string& disposition,
    std::size_t motion_order)
{
    if (root_id.empty()) {
        throw EventTraceVerificationError(
            "trace event requires an algebraic root identity");
    }
    if (global_fibre_id.empty()) {
        throw EventTraceVerificationError(
            "trace event requires a global fibre identity");
    }
    if (kind.empty()) {
        throw EventTraceVerificationError(
            "trace event kind must not be empty");
    }
    if (disposition.empty()) {
        throw EventTraceVerificationError(
            "trace event disposition must not be empty");
    }
    if (multiplicity == 0) {
        throw EventTraceVerificationError(
            "trace event multiplicity must be positive");
    }
    canonicalize_ids(feature_ids, "trace feature IDs");
    canonicalize_ids(branch_ids, "trace branch IDs");

    EventTraceEvent2 event{
        root_id,
        global_fibre_id,
        kind,
        std::move(feature_ids),
        std::move(branch_ids),
        multiplicity,
        disposition,
        motion_order,
        {},
        {},
    };
    event.canonical_bytes =
        encode_canonical_record(
            "event-trace-event-v1",
            {
                event.root_id,
                event.global_fibre_id,
                event.kind,
                encode_string_sequence(
                    event.feature_ids),
                encode_string_sequence(
                    event.branch_ids),
                std::to_string(event.multiplicity),
                event.disposition,
            });
    event.canonical_id =
        sha256_bytes(event.canonical_bytes);
    return event;
}

EventTrace2 build_event_trace(
    const EventPartitionCertificate2& partition,
    const std::string& motion_chart_id,
    const std::string& motion_identity,
    const std::string& effective_cap_bytes,
    ContinuousTeaVerdict verdict,
    const std::string& whole_rim_disposition,
    const std::string& oracle_strategy_version,
    std::vector<EventTraceEvent2> events)
{
    const VerifiedEventPartition2 verified =
        verify_event_partition(partition);
    if (verified.verdict
        != ContinuousTeaVerdict::CERTIFIED) {
        throw EventTraceVerificationError(
            "event trace requires a reconstruct-verified partition");
    }
    if (motion_chart_id.empty()
        || motion_identity.empty()) {
        throw EventTraceVerificationError(
            "event trace requires exact motion identity");
    }
    if (effective_cap_bytes.empty()) {
        throw EventTraceVerificationError(
            "event trace requires canonical effective-cap bytes");
    }
    if (oracle_strategy_version.empty()) {
        throw EventTraceVerificationError(
            "event trace requires an oracle strategy version");
    }
    if (whole_rim_disposition != "clear"
        && whole_rim_disposition != "material"
        && whole_rim_disposition != "partial"
        && whole_rim_disposition != "unresolved") {
        throw EventTraceVerificationError(
            "event trace has an unknown whole-rim disposition");
    }
    for (const EventTraceEvent2& event : events) {
        if (!has_root(
                verified.partition,
                event.root_id)) {
            throw EventTraceVerificationError(
                "trace event root is absent from the verified partition");
        }
    }
    std::sort(
        events.begin(),
        events.end(),
        trace_event_less);

    std::vector<std::string> event_bytes;
    event_bytes.reserve(events.size());
    for (const EventTraceEvent2& event : events) {
        event_bytes.push_back(event.canonical_bytes);
    }

    validate_decision_authority(
        verified.partition.canonical_bytes,
        verified.partition.canonical_digest);
    EventTrace2 trace{
        verdict,
        verified.partition,
        verified.partition.canonical_bytes,
        verified.partition.canonical_digest,
        std::move(events),
        motion_chart_id,
        motion_identity,
        effective_cap_bytes,
        whole_rim_disposition,
        oracle_strategy_version,
        verified.partition.cells.size(),
        {},
        {},
    };
    trace.canonical_bytes =
        encode_canonical_record(
            "event-trace-v2",
            {
                trace.motion_chart_id,
                trace.motion_identity,
                trace.effective_cap_bytes,
                verdict_text(trace.verdict),
                trace.whole_rim_disposition,
                trace.oracle_strategy_version,
                trace.partition.canonical_digest,
                trace.decision_authority_digest,
                encode_string_sequence(event_bytes),
            });
    trace.canonical_digest =
        sha256_bytes(trace.canonical_bytes);
    return trace;
}
