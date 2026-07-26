#include "segment_oracle.h"

#include "../stock_2.h"
#include "event_certificate.h"
#include "segment_partition.h"
#include "segment_source.h"
#include "segment_strata.h"
#include "sha256.h"

#include <algorithm>
#include <map>
#include <type_traits>
#include <string_view>
#include <utility>

#include <CGAL/CORE/BigRat.h>

namespace {

using Integer = CORE::BigInt;
using Rational = CORE::BigRat;

Rational parse_rational(
    const std::string& text,
    std::string_view role)
{
    const std::size_t separator = text.find('/');
    try {
        if (separator == std::string::npos) {
            return Rational(Integer(text));
        }
        const Integer numerator(
            text.substr(0, separator));
        const Integer denominator(
            text.substr(separator + 1));
        if (denominator == 0) {
            throw IncompleteSegmentOracleError(
                std::string(role)
                + " has zero denominator");
        }
        return Rational(numerator, denominator);
    } catch (const EventSubstrateError&) {
        throw;
    } catch (const std::exception&) {
        throw IncompleteSegmentOracleError(
            std::string(role)
            + " is not an exact rational");
    }
}

Epeck::FT exact_ft(const Rational& value)
{
    static_assert(
        std::is_same_v<Rational, CGAL::Epeck_ft>);
    return Epeck::FT(value);
}

enum class CellDecision {
    CLEAR,
    MATERIAL,
    SAFE_PARTIAL,
    CAP_EXCEEDED,
    UNRESOLVED,
};

bool reference_is_material(
    const Stock2& stock,
    const SegmentEventSource2& source,
    const std::string& numerator,
    const std::string& denominator,
    bool& resolved)
{
    const Rational parameter{
        Integer(numerator),
        Integer(denominator)};
    const Rational x0 =
        parse_rational(source.x0().text(), "start x");
    const Rational y0 =
        parse_rational(source.y0().text(), "start y");
    const Rational x1 =
        parse_rational(source.x1().text(), "end x");
    const Rational y1 =
        parse_rational(source.y1().text(), "end y");
    const Rational radius =
        parse_rational(
            source.tool_radius().text(),
            "tool radius");
    const GpsPoint reference(
        exact_ft(x0 + parameter * (x1 - x0)),
        exact_ft(y0 + parameter * (y1 - y0) - radius));
    const CGAL::Oriented_side side =
        stock.set().oriented_side(reference);
    resolved = side != CGAL::ON_ORIENTED_BOUNDARY;
    return side == CGAL::ON_POSITIVE_SIDE;
}

const BranchPairDisposition2* pair_for(
    const SegmentEventStratum2& stratum,
    const std::string& first,
    const std::string& second)
{
    const auto found = std::find_if(
        stratum.pair_dispositions.begin(),
        stratum.pair_dispositions.end(),
        [&first, &second](
            const BranchPairDisposition2& pair) {
            return pair.first_branch_id == first
                && pair.second_branch_id == second;
        });
    return found == stratum.pair_dispositions.end()
        ? nullptr
        : &*found;
}

CellDecision classify_cell(
    const Stock2& stock,
    const SegmentEventSource2& source,
    const SegmentEventStratum2& stratum)
{
    if (stratum.kind != "cell"
        || stratum.witness_numerator.empty()
        || stratum.witness_denominator.empty()) {
        return CellDecision::UNRESOLVED;
    }
    bool reference_resolved = false;
    const bool reference_material =
        reference_is_material(
            stock,
            source,
            stratum.witness_numerator,
            stratum.witness_denominator,
            reference_resolved);
    if (!reference_resolved) {
        return CellDecision::UNRESOLVED;
    }
    const std::size_t count =
        stratum.active_branch_ids.size();
    if (count == 0) {
        return reference_material
            ? CellDecision::MATERIAL
            : CellDecision::CLEAR;
    }
    if (count % 2 != 0) {
        return CellDecision::UNRESOLVED;
    }

    std::vector<std::pair<std::size_t, std::size_t>>
        material_runs;
    if (reference_material) {
        material_runs.emplace_back(count - 1, 0);
        for (std::size_t index = 1;
             index + 1 < count;
             index += 2) {
            material_runs.emplace_back(
                index,
                index + 1);
        }
    } else {
        for (std::size_t index = 0;
             index + 1 < count;
             index += 2) {
            material_runs.emplace_back(
                index,
                index + 1);
        }
    }
    for (const auto [first_index, second_index] :
         material_runs) {
        const BranchPairDisposition2* pair =
            pair_for(
                stratum,
                stratum.active_branch_ids[first_index],
                stratum.active_branch_ids[second_index]);
        if (pair == nullptr) {
            return CellDecision::UNRESOLVED;
        }
        if (pair->cap_disposition == "above-cap") {
            return CellDecision::CAP_EXCEEDED;
        }
        if (pair->cap_disposition != "below-cap"
            && pair->cap_disposition != "equal-cap") {
            return CellDecision::UNRESOLVED;
        }
    }
    return CellDecision::SAFE_PARTIAL;
}

const AlgebraicRootRecord2& root_for(
    const EventPartitionCertificate2& certificate,
    const std::string& root_id)
{
    const auto found = std::find_if(
        certificate.roots.begin(),
        certificate.roots.end(),
        [&root_id](
            const AlgebraicRootRecord2& root) {
            return root.root_id == root_id;
        });
    if (found == certificate.roots.end()) {
        throw IncompleteSegmentOracleError(
            "segment trace fibre lacks its root");
    }
    return *found;
}

std::string trace_kind(
    const PartitionEvent2& event)
{
    if (event.kind == "cap-crossing") {
        return "cap";
    }
    if (event.kind == "support-overlap") {
        return "overlap";
    }
    if (event.kind == "endpoint-order"
        && event.left_active_count
            != event.right_active_count
        && !event.vertex_id.empty()) {
        return "vertex";
    }
    return event.kind;
}

std::string trace_disposition(
    const PartitionEvent2& event,
    const std::string& kind)
{
    if (kind == "cap") {
        return "equal";
    }
    if (kind == "vertex") {
        return "incident";
    }
    return event.disposition;
}

void append_unique(
    std::vector<std::string>& target,
    const std::string& value)
{
    if (!value.empty()
        && std::find(
               target.begin(),
               target.end(),
               value)
            == target.end()) {
        target.push_back(value);
    }
}

struct TraceAccumulator {
    std::string kind;
    std::string disposition;
    std::vector<std::string> feature_ids;
    std::vector<std::string> branch_ids;
    unsigned int multiplicity = 1;
};

std::string verdict_text(ContinuousTeaVerdict verdict)
{
    if (verdict == ContinuousTeaVerdict::CERTIFIED) {
        return "certified";
    }
    if (verdict == ContinuousTeaVerdict::CAP_EXCEEDED) {
        return "cap_exceeded";
    }
    return "unresolved";
}

EventTrace2 build_segment_trace(
    const VerifiedSegmentEventPartition2& verified,
    ContinuousTeaVerdict verdict,
    const std::string& whole_rim_disposition,
    std::vector<EventTraceEvent2> events)
{
    if (verified.verdict
        != ContinuousTeaVerdict::CERTIFIED) {
        throw IncompleteSegmentOracleError(
            "segment trace requires a reconstructed Task 5 partition");
    }
    std::sort(
        events.begin(),
        events.end(),
        [](const EventTraceEvent2& left,
           const EventTraceEvent2& right) {
            if (left.motion_order != right.motion_order) {
                return left.motion_order
                    < right.motion_order;
            }
            if (left.global_fibre_id
                != right.global_fibre_id) {
                return left.global_fibre_id
                    < right.global_fibre_id;
            }
            return left.canonical_id
                < right.canonical_id;
        });
    std::vector<std::string> event_bytes;
    event_bytes.reserve(events.size());
    for (const EventTraceEvent2& event : events) {
        event_bytes.push_back(event.canonical_bytes);
    }
    EventTrace2 trace{
        verdict,
        verified.partition.certificate,
        std::move(events),
        "segment-linear-0-1-v1",
        verified.partition.source.canonical_bytes(),
        verified.partition.source.cap_chord_ratio()
            .canonical_bytes(),
        whole_rim_disposition,
        "segment-event-exact-v1",
        verified.partition.certificate.cells.size(),
        {},
        {},
    };
    trace.canonical_bytes =
        encode_canonical_record(
            "event-trace-v1",
            {
                trace.motion_chart_id,
                trace.motion_identity,
                trace.effective_cap_bytes,
                verdict_text(trace.verdict),
                trace.whole_rim_disposition,
                trace.oracle_strategy_version,
                trace.partition.canonical_digest,
                encode_string_sequence(event_bytes),
            });
    trace.canonical_digest =
        sha256_bytes(trace.canonical_bytes);
    return trace;
}

std::vector<EventTraceEvent2> trace_events(
    const SegmentEventPartition2& partition)
{
    std::vector<EventTraceEvent2> result;
    for (std::size_t motion_order = 0;
         motion_order < partition.strata.size();
         ++motion_order) {
        const SegmentEventStratum2& stratum =
            partition.strata[motion_order];
        if (stratum.kind != "fibre") {
            continue;
        }
        const AlgebraicRootRecord2& root =
            root_for(
                partition.certificate,
                stratum.root_id);
        std::map<std::string, TraceAccumulator>
            accumulators;
        for (const PartitionEvent2& event :
             stratum.events) {
            const std::string kind =
                trace_kind(event);
            const std::string key =
                kind
                + (kind == "endpoint-order"
                       ? event.vertex_id
                       : std::string{});
            TraceAccumulator& accumulator =
                accumulators[key];
            accumulator.kind = kind;
            accumulator.disposition =
                trace_disposition(event, kind);
            accumulator.multiplicity =
                kind == "tangent"
                ? 2U
                : std::max(
                      accumulator.multiplicity,
                      root.multiplicity);
            append_unique(
                accumulator.feature_ids,
                event.feature_id);
            append_unique(
                accumulator.feature_ids,
                event.first_feature_id);
            append_unique(
                accumulator.feature_ids,
                event.second_feature_id);
            append_unique(
                accumulator.branch_ids,
                event.branch_id);
            append_unique(
                accumulator.branch_ids,
                event.first_branch_id);
            append_unique(
                accumulator.branch_ids,
                event.second_branch_id);
        }
        for (auto& [key, accumulator] :
             accumulators) {
            static_cast<void>(key);
            if (accumulator.feature_ids.empty()) {
                continue;
            }
            if (accumulator.branch_ids.empty()) {
                for (const std::string& branch_id :
                     stratum.active_branch_ids) {
                    const auto branch = std::find_if(
                        partition.branches.begin(),
                        partition.branches.end(),
                        [&branch_id](
                            const SegmentBoundaryBranch2& candidate) {
                            return candidate.branch_id
                                == branch_id;
                        });
                    if (branch != partition.branches.end()
                        && std::find(
                               accumulator.feature_ids.begin(),
                               accumulator.feature_ids.end(),
                               branch->feature_id)
                            != accumulator.feature_ids.end()) {
                        append_unique(
                            accumulator.branch_ids,
                            branch_id);
                    }
                }
            }
            if (accumulator.branch_ids.empty()) {
                continue;
            }
            result.push_back(
                make_event_trace_event(
                    stratum.root_id,
                    stratum.global_fibre_id,
                    accumulator.kind,
                    std::move(
                        accumulator.feature_ids),
                    std::move(
                        accumulator.branch_ids),
                    accumulator.multiplicity,
                    accumulator.disposition,
                    motion_order));
        }
    }
    return result;
}

} // namespace

SegmentTeaAudit2 audit_segment_tea_event_exact(
    const Stock2& stock,
    const SegmentEventSource2& source)
{
    const SegmentEventPartition2 partition =
        construct_segment_event_partition(
            stock,
            source);
    bool any_partial = false;
    bool any_clear = false;
    bool any_material = false;
    bool unresolved =
        !partition.certificate.overlaps.empty();
    bool cap_exceeded = false;
    for (const SegmentEventStratum2& stratum :
         partition.strata) {
        if (stratum.kind == "fibre") {
            unresolved = unresolved
                || !stratum.algebraic_root_evaluated
                || !stratum.original_equations_rechecked
                || !stratum.orientation_rechecked
                || !stratum.trim_predicates_rechecked;
            continue;
        }
        const CellDecision decision =
            classify_cell(stock, source, stratum);
        any_clear = any_clear
            || decision == CellDecision::CLEAR;
        any_material = any_material
            || decision == CellDecision::MATERIAL;
        any_partial = any_partial
            || decision == CellDecision::SAFE_PARTIAL
            || decision == CellDecision::CAP_EXCEEDED;
        cap_exceeded = cap_exceeded
            || decision == CellDecision::MATERIAL
            || decision == CellDecision::CAP_EXCEEDED;
        unresolved = unresolved
            || decision == CellDecision::UNRESOLVED;
    }
    const ContinuousTeaVerdict verdict =
        cap_exceeded
        ? ContinuousTeaVerdict::CAP_EXCEEDED
        : (
              unresolved
                  ? ContinuousTeaVerdict::
                        UNRESOLVED_DEGENERACY
                  : ContinuousTeaVerdict::CERTIFIED);
    const std::string whole_rim_disposition =
        any_partial
        ? "partial"
        : (
              any_material && !any_clear
                  ? "material"
                  : (
                        any_clear && !any_material
                            ? "clear"
                            : "unresolved"));
    const VerifiedSegmentEventPartition2 verified =
        verify_segment_event_partition(
            stock,
            source,
            partition);
    return {
        verdict,
        build_segment_trace(
            verified,
            verdict,
            whole_rim_disposition,
            trace_events(partition)),
    };
}

bool segment_station_cap_exceeded_exact(
    const Stock2& stock,
    const SegmentEventSource2& source,
    const std::string& witness_numerator,
    const std::string& witness_denominator)
{
    const SegmentCellStratum2 cell =
        construct_segment_cell_stratum(
            stock,
            source,
            witness_numerator,
            witness_denominator);
    const CellDecision decision =
        classify_cell(
            stock,
            source,
            cell.stratum);
    if (decision == CellDecision::UNRESOLVED) {
        throw IncompleteSegmentOracleError(
            "exact station disposition is unresolved");
    }
    return decision == CellDecision::MATERIAL
        || decision == CellDecision::CAP_EXCEEDED;
}
