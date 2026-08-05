#include "segment_oracle.h"

#include "../exact_stock_region_2.h"
#include "../exact_sweep_2.h"
#include "../stock_2.h"
#include "boundary_events.h"
#include "event_certificate.h"
#include "segment_partition.h"
#include "segment_source.h"
#include "segment_strata.h"
#include "sha256.h"
#include "station_classifier.h"

#include <algorithm>
#include <cstddef>
#include <map>
#include <string_view>
#include <utility>

#include <CGAL/CORE/BigRat.h>

namespace {

using Integer = CORE::BigInt;
using Rational = CORE::BigRat;

const std::string SWEPT_PREFIX_STRATEGY_VERSION =
    "swept-prefix-segment-tea-exact-v1";
const std::string SWEPT_PREFIX_THEOREM_VERSION =
    "translation-swept-prefix-forward-semicircle-v1";
constexpr std::size_t SWEPT_PREFIX_MOTION_STRATA = 2;

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

bool overlap_requires_resolution(
    const OverlapInterval2& overlap)
{
    return overlap.kind
        != "identically-equal-cap-interval";
}

StationEventSource2 station_source(
    const SegmentEventSource2& source,
    const std::string& numerator,
    const std::string& denominator)
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
    const Rational cap =
        parse_rational(
            source.cap_chord_ratio().text(),
            "cap chord ratio");
    const auto text =
        [](const Rational& value) {
            const Integer numerator =
                CORE::numerator(value);
            const Integer denominator =
                CORE::denominator(value);
            return denominator == 1
                ? numerator.convert_to<std::string>()
                : numerator.convert_to<std::string>()
                    + "/"
                    + denominator.convert_to<std::string>();
        };
    return StationEventSource2::build(
        text(x0 + parameter * (x1 - x0)),
        text(y0 + parameter * (y1 - y0)),
        text(radius),
        text(cap));
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

std::string segment_stock_boundary_digest(
    const Stock2& stock)
{
    return sha256_bytes(
        canonical_stock_boundary_identity(stock));
}

bool start_disk_has_no_material_interior(
    const Stock2& stock,
    const SegmentEventSource2& source)
{
    const ReachKernelPoint start(
        ReachFT(parse_rational(
            source.x0().text(),
            "swept-prefix start x")),
        ReachFT(parse_rational(
            source.y0().text(),
            "swept-prefix start y")));
    const ReachFT radius(parse_rational(
        source.tool_radius().text(),
        "swept-prefix tool radius"));
    ReachSet overlap;
    overlap.insert(reach_disk_polygon(start, radius));
    overlap.intersection(lift_exact_stock_region(stock));
    return overlap.is_empty();
}

bool cap_is_exact_pi(
    const SegmentEventSource2& source)
{
    // chord_ratio = 4 sin^2(theta/2), hence theta = pi iff ratio = 4.
    const Rational pi_cap_chord_ratio(4);
    return parse_rational(
               source.cap_chord_ratio().text(),
               "swept-prefix cap chord ratio")
        == pi_cap_chord_ratio;
}

std::string canonical_swept_prefix_audit(
    const std::string& source_canonical_bytes,
    const std::string& stock_boundary_digest,
    bool start_clear,
    bool pi_cap,
    ContinuousTeaVerdict verdict)
{
    return encode_canonical_record(
        "swept-prefix-segment-tea-proof-v1",
        {
            SWEPT_PREFIX_STRATEGY_VERSION,
            SWEPT_PREFIX_THEOREM_VERSION,
            stock_boundary_digest,
            source_canonical_bytes,
            start_clear
                ? "start-disk-clear"
                : "start-disk-overlaps-material",
            pi_cap
                ? "effective-cap-exact-pi"
                : "effective-cap-not-exact-pi",
            verdict_text(verdict),
            std::to_string(
                SWEPT_PREFIX_MOTION_STRATA),
        });
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
    if (verified.partition.canonical_bytes.empty()
        || verified.partition.canonical_digest.empty()
        || sha256_bytes(
               verified.partition.canonical_bytes)
            != verified.partition.canonical_digest) {
        throw IncompleteSegmentOracleError(
            "segment trace decision authority is not digest-valid");
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
        verified.partition.canonical_bytes,
        verified.partition.canonical_digest,
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

const std::string& swept_prefix_strategy_version()
{
    return SWEPT_PREFIX_STRATEGY_VERSION;
}

const std::string& swept_prefix_theorem_version()
{
    return SWEPT_PREFIX_THEOREM_VERSION;
}

bool swept_prefix_segment_tea_audit_is_self_consistent(
    const SweptPrefixSegmentTeaAudit2& audit)
{
    const ContinuousTeaVerdict expected_verdict =
        audit.start_disk_clear && audit.exact_pi_cap
        ? ContinuousTeaVerdict::CERTIFIED
        : ContinuousTeaVerdict::UNRESOLVED_DEGENERACY;
    const std::string expected_bytes =
        canonical_swept_prefix_audit(
            audit.source_canonical_bytes,
            audit.stock_boundary_digest,
            audit.start_disk_clear,
            audit.exact_pi_cap,
            expected_verdict);
    return audit.verdict == expected_verdict
        && audit.strategy_version
            == SWEPT_PREFIX_STRATEGY_VERSION
        && audit.theorem_version
            == SWEPT_PREFIX_THEOREM_VERSION
        && audit.motion_stratum_count
            == SWEPT_PREFIX_MOTION_STRATA
        && audit.canonical_bytes == expected_bytes
        && audit.canonical_digest
            == sha256_bytes(expected_bytes);
}

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
        std::any_of(
            partition.certificate.overlaps.begin(),
            partition.certificate.overlaps.end(),
            overlap_requires_resolution);
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
        const StationCellDecision decision =
            classify_station_cell(
                stock,
                station_source(
                    source,
                    stratum.witness_numerator,
                    stratum.witness_denominator),
                stratum)
                .decision;
        any_clear = any_clear
            || decision == StationCellDecision::CLEAR;
        any_material = any_material
            || decision == StationCellDecision::MATERIAL;
        any_partial = any_partial
            || decision
                == StationCellDecision::SAFE_PARTIAL
            || decision
                == StationCellDecision::CAP_EXCEEDED;
        cap_exceeded = cap_exceeded
            || decision == StationCellDecision::MATERIAL
            || decision
                == StationCellDecision::CAP_EXCEEDED;
        unresolved = unresolved
            || decision
                == StationCellDecision::UNRESOLVED;
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

SweptPrefixSegmentTeaAudit2
audit_swept_prefix_segment_tea_exact(
    const Stock2& stock,
    const SegmentEventSource2& source)
{
    bool start_clear;
    try {
        start_clear =
            start_disk_has_no_material_interior(stock, source);
    } catch (const ExactStockRegionLiftError& error) {
        throw IncompleteSegmentOracleError(
            std::string(
                "swept-prefix start-clearance lift is incomplete: ")
            + error.what());
    } catch (const ReachableDomainConstructionError& error) {
        throw IncompleteSegmentOracleError(
            std::string(
                "swept-prefix start-clearance construction is incomplete: ")
            + error.what());
    }
    const bool pi_cap = cap_is_exact_pi(source);
    const ContinuousTeaVerdict verdict =
        start_clear && pi_cap
        ? ContinuousTeaVerdict::CERTIFIED
        : ContinuousTeaVerdict::UNRESOLVED_DEGENERACY;
    const std::string stock_boundary_digest =
        segment_stock_boundary_digest(stock);
    const std::string source_canonical_bytes =
        source.canonical_bytes();
    const std::string canonical_bytes =
        canonical_swept_prefix_audit(
            source_canonical_bytes,
            stock_boundary_digest,
            start_clear,
            pi_cap,
            verdict);
    SweptPrefixSegmentTeaAudit2 audit{
        verdict,
        canonical_bytes,
        sha256_bytes(canonical_bytes),
        SWEPT_PREFIX_STRATEGY_VERSION,
        SWEPT_PREFIX_THEOREM_VERSION,
        source_canonical_bytes,
        stock_boundary_digest,
        SWEPT_PREFIX_MOTION_STRATA,
        start_clear,
        pi_cap,
    };
    if (!swept_prefix_segment_tea_audit_is_self_consistent(
            audit)) {
        throw IncompleteSegmentOracleError(
            "swept-prefix audit failed native self-verification");
    }
    return audit;
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
    const StationCellDecision decision =
        classify_station_cell(
            stock,
            station_source(
                source,
                witness_numerator,
                witness_denominator),
            cell.stratum)
            .decision;
    if (decision
        == StationCellDecision::UNRESOLVED) {
        throw IncompleteSegmentOracleError(
            "exact station disposition is unresolved");
    }
    return decision == StationCellDecision::MATERIAL
        || decision
            == StationCellDecision::CAP_EXCEEDED;
}
