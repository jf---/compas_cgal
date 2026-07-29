#include "circle_strata.h"

#include "../stock_2.h"
#include "boundary_events.h"
#include "circle_pair_projection.h"
#include "circle_projection.h"
#include "circle_source.h"
#include "event_certificate.h"
#include "segment_partition.h"
#include "segment_strata.h"
#include "sha256.h"
#include "station_classifier.h"
#include "station_source.h"

#include <algorithm>
#include <array>
#include <string_view>
#include <utility>
#include <vector>

#include <CGAL/CORE/BigRat.h>

namespace {

using Integer = CORE::BigInt;
using Rational = CORE::BigRat;

constexpr std::array<const char*, 4>
    CENTER_CHART_IDS{
        "center-quarter-0-v1",
        "center-quarter-1-v1",
        "center-quarter-2-v1",
        "center-quarter-3-v1",
    };

struct ChartWitness {
    std::size_t chart;
    Rational local_parameter;
};

Rational parse_rational(
    const std::string& text,
    std::string_view role)
{
    const std::size_t separator = text.find('/');
    try {
        if (separator == std::string::npos) {
            return Rational(Integer(text));
        }
        if (text.find('/', separator + 1)
            != std::string::npos) {
            throw IncompleteFullCircleCellAuthorityError(
                std::string(role)
                + " is not one exact rational");
        }
        const Integer numerator(
            text.substr(0, separator));
        const Integer denominator(
            text.substr(separator + 1));
        if (denominator == 0) {
            throw IncompleteFullCircleCellAuthorityError(
                std::string(role)
                + " has zero denominator");
        }
        return Rational(numerator, denominator);
    } catch (const EventSubstrateError&) {
        throw;
    } catch (const std::exception&) {
        throw IncompleteFullCircleCellAuthorityError(
            std::string(role)
            + " is not one exact rational");
    }
}

std::string rational_text(const Rational& value)
{
    const Integer numerator =
        CORE::numerator(value);
    const Integer denominator =
        CORE::denominator(value);
    return denominator == 1
        ? numerator.convert_to<std::string>()
        : numerator.convert_to<std::string>()
            + "/"
            + denominator.convert_to<std::string>();
}

ChartWitness chart_witness(const Rational& global)
{
    if (global < 0 || global >= 1) {
        throw IncompleteFullCircleCellAuthorityError(
            "full-circle cell witness lies outside the global chart domain");
    }
    const Rational scaled = Rational(4) * global;
    for (std::size_t chart = 0;
         chart < CENTER_CHART_IDS.size();
         ++chart) {
        if (scaled < Rational(chart + 1)) {
            return {
                chart,
                scaled - Rational(chart),
            };
        }
    }
    throw IncompleteFullCircleCellAuthorityError(
        "full-circle cell witness has no owned quarter chart");
}

std::pair<Rational, Rational> unit_direction(
    std::size_t chart,
    const Rational& parameter)
{
    const Rational denominator =
        Rational(1) + parameter * parameter;
    Rational x =
        (Rational(1) - parameter * parameter)
        / denominator;
    Rational y =
        Rational(2) * parameter / denominator;
    if (chart == 1) {
        return {-y, x};
    }
    if (chart == 2) {
        return {-x, -y};
    }
    if (chart == 3) {
        return {y, -x};
    }
    if (chart != 0) {
        throw IncompleteFullCircleCellAuthorityError(
            "full-circle cell uses an unknown quarter chart");
    }
    return {x, y};
}

StationEventSource2 station_source(
    double center_x,
    double center_y,
    double phase_dx,
    double phase_dy,
    double tool_radius,
    double cap_chord_ratio,
    const ChartWitness& witness)
{
    const Rational exact_center_x =
        parse_rational(
            exact_rational_text(
                Epeck::FT(center_x)),
            "circle center x");
    const Rational exact_center_y =
        parse_rational(
            exact_rational_text(
                Epeck::FT(center_y)),
            "circle center y");
    const Rational exact_phase_x =
        parse_rational(
            exact_rational_text(
                Epeck::FT(phase_dx)),
            "circle phase x");
    const Rational exact_phase_y =
        parse_rational(
            exact_rational_text(
                Epeck::FT(phase_dy)),
            "circle phase y");
    const auto [unit_x, unit_y] =
        unit_direction(
            witness.chart,
            witness.local_parameter);
    const Rational station_x =
        exact_center_x
        + exact_phase_x * unit_x
        - exact_phase_y * unit_y;
    const Rational station_y =
        exact_center_y
        + exact_phase_y * unit_x
        + exact_phase_x * unit_y;
    return StationEventSource2::build(
        rational_text(station_x),
        rational_text(station_y),
        exact_rational_text(
            Epeck::FT(tool_radius)),
        exact_rational_text(
            Epeck::FT(cap_chord_ratio)));
}

const SegmentBoundaryBranch2* branch_for(
    const SegmentCellStratum2& cell,
    const std::string& branch_id)
{
    const auto found = std::find_if(
        cell.branches.begin(),
        cell.branches.end(),
        [&branch_id](
            const SegmentBoundaryBranch2& branch) {
            return branch.branch_id == branch_id;
        });
    return found == cell.branches.end()
        ? nullptr
        : &*found;
}

std::vector<std::string> projection_fields(
    const std::string& encoded)
{
    try {
        return decode_string_sequence(encoded);
    } catch (const EventSubstrateError&) {
        return {};
    }
}

bool pair_projection_matches(
    const std::vector<std::string>& fields,
    const std::string& projection_tag,
    const SegmentBoundaryBranch2& first,
    const SegmentBoundaryBranch2& second,
    const std::string& chart_id)
{
    if (fields.size() != 8
        || fields[0] != projection_tag
        || fields[5] != chart_id) {
        return false;
    }
    const bool forward =
        fields[1] == first.feature_id
        && fields[2] == second.feature_id
        && fields[3] == first.support_id
        && fields[4] == second.support_id
        && fields[6] == first.rim_chart_id
        && fields[7] == second.rim_chart_id;
    const bool reverse =
        fields[1] == second.feature_id
        && fields[2] == first.feature_id
        && fields[3] == second.support_id
        && fields[4] == first.support_id
        && fields[6] == second.rim_chart_id
        && fields[7] == first.rim_chart_id;
    return forward || reverse;
}

std::vector<std::string> cap_projections_for(
    const EventPartitionCertificate2& partition,
    const SegmentCellStratum2& cell,
    const BranchPairDisposition2& run,
    const std::string& chart_id,
    bool cap_is_pi)
{
    const SegmentBoundaryBranch2* first =
        branch_for(cell, run.first_branch_id);
    const SegmentBoundaryBranch2* second =
        branch_for(cell, run.second_branch_id);
    if (first == nullptr || second == nullptr) {
        return {};
    }
    const bool same_support =
        first->support_id == second->support_id
        && first->support_kind
            == second->support_kind;
    if (!same_support) {
        if (run.cap_disposition == "equal-cap") {
            return {};
        }
        std::vector<std::string> result;
        for (const std::string& projection_tag :
             {
                 std::string(
                     "full-circle-pair-orientation-v1"),
                 std::string(
                     "full-circle-pair-cap-v2"),
             }) {
            if (cap_is_pi
                && projection_tag
                    == "full-circle-pair-cap-v2") {
                continue;
            }
            std::string matched;
            for (const ProjectionRecord2& projection :
                 partition.projections) {
                if (pair_projection_matches(
                        projection_fields(
                            projection.projection_id),
                        projection_tag,
                        *first,
                        *second,
                        chart_id)) {
                    matched =
                        projection.projection_id;
                    break;
                }
            }
            if (matched.empty()) {
                return {};
            }
            result.push_back(std::move(matched));
        }
        std::sort(result.begin(), result.end());
        return result;
    }
    const std::string prefix =
        first->support_kind == "line"
        ? "full-circle-line-cap-v1"
        : (
              first->support_kind == "circle"
                  ? "full-circle-circle-cap-v1"
                  : std::string{});
    if (prefix.empty()) {
        return {};
    }
    const auto owns_feature =
        [first, second](
            const std::string& feature_id) {
            return feature_id == first->feature_id
                || feature_id
                    == second->feature_id;
        };
    if (run.cap_disposition != "equal-cap") {
        for (const ProjectionRecord2& projection :
             partition.projections) {
            const std::vector<std::string> fields =
                projection_fields(
                    projection.projection_id);
            if (fields.size() == 3
                && fields[0] == prefix
                && owns_feature(fields[1])
                && fields[2] == chart_id) {
                return {
                    projection.projection_id,
                };
            }
        }
        return {};
    }
    for (const OverlapInterval2& overlap :
         partition.overlaps) {
        const std::vector<std::string> fields =
            projection_fields(overlap.branch_id);
        if (overlap.kind
                == "identically-equal-cap-interval"
            && fields.size() == 3
            && fields[0] == prefix
            && owns_feature(fields[1])
            && fields[2] == chart_id) {
            return {
                overlap.branch_id,
            };
        }
    }
    return {};
}

std::string decision_text(
    StationCellDecision decision)
{
    switch (decision) {
    case StationCellDecision::CLEAR:
        return "clear";
    case StationCellDecision::MATERIAL:
        return "material";
    case StationCellDecision::SAFE_PARTIAL:
        return "safe-partial";
    case StationCellDecision::CAP_EXCEEDED:
        return "cap-exceeded";
    case StationCellDecision::UNRESOLVED:
        return "unresolved";
    }
    throw IncompleteFullCircleCellAuthorityError(
        "full-circle cell has an unknown decision");
}

bool overlap_requires_resolution(
    const OverlapInterval2& overlap)
{
    return overlap.kind
        != "identically-equal-cap-interval";
}

} // namespace

std::vector<std::string>
derive_full_circle_pair_requests(
    const Stock2& stock,
    const EventPartitionCertificate2& topology_partition,
    double center_x,
    double center_y,
    double phase_dx,
    double phase_dy,
    double tool_radius,
    double cap_chord_ratio)
{
    const std::vector<BoundaryFeatureRecord2>
        boundary_records =
            extract_boundary_records(stock);
    if (boundary_records.empty()
        || topology_partition.cells.empty()) {
        throw IncompleteFullCircleCellAuthorityError(
            "full-circle pair demand requires nonempty topology cells and stock features");
    }
    std::vector<std::string> requests;
    for (const ParameterCell2& parameter_cell :
         topology_partition.cells) {
        const Rational global_witness(
            Integer(
                parameter_cell.witness_numerator),
            Integer(
                parameter_cell.witness_denominator));
        const ChartWitness witness =
            chart_witness(global_witness);
        const StationEventSource2 source =
            station_source(
                center_x,
                center_y,
                phase_dx,
                phase_dy,
                tool_radius,
                cap_chord_ratio,
                witness);
        const SegmentCellStratum2 cell =
            construct_station_cell_stratum(
                boundary_records,
                source);
        for (std::size_t first_index = 0;
             first_index < cell.branches.size();
             ++first_index) {
            for (std::size_t second_index =
                     first_index + 1;
                 second_index < cell.branches.size();
                 ++second_index) {
                const SegmentBoundaryBranch2& first =
                    cell.branches[first_index];
                const SegmentBoundaryBranch2& second =
                    cell.branches[second_index];
                if (first.support_kind
                        == second.support_kind
                    && first.support_id
                        == second.support_id) {
                    continue;
                }
                requests.push_back(
                    FullCirclePairRequest2::build(
                        CENTER_CHART_IDS[
                            witness.chart],
                        first.feature_id,
                        first.support_id,
                        first.support_kind,
                        first.rim_chart_id,
                        second.feature_id,
                        second.support_id,
                        second.support_kind,
                        second.rim_chart_id)
                        .canonical_source());
            }
        }
    }
    std::sort(
        requests.begin(),
        requests.end());
    requests.erase(
        std::unique(
            requests.begin(),
            requests.end()),
        requests.end());
    return requests;
}

EventPartitionCertificate2
construct_full_circle_pair_closed_partition(
    const Stock2& stock,
    const std::string& stock_identity,
    const std::vector<std::string>& motion_data,
    const std::string& cutter_radius,
    const std::string& cap_chord_ratio,
    const std::vector<std::string>& line_sources,
    const std::vector<std::string>& circle_sources,
    double center_x,
    double center_y,
    double phase_dx,
    double phase_dy,
    double tool_radius,
    double cap_chord_ratio_value)
{
    std::vector<std::string> pair_requests;
    while (true) {
        EventPartitionCertificate2 partition =
            construct_full_circle_boundary_pullback_partition(
                stock_identity,
                motion_data,
                cutter_radius,
                cap_chord_ratio,
                line_sources,
                circle_sources,
                pair_requests);
        const std::vector<std::string>
            observed_requests =
                derive_full_circle_pair_requests(
                    stock,
                    partition,
                    center_x,
                    center_y,
                    phase_dx,
                    phase_dy,
                    tool_radius,
                    cap_chord_ratio_value);
        std::vector<std::string> closed_requests;
        std::set_union(
            pair_requests.begin(),
            pair_requests.end(),
            observed_requests.begin(),
            observed_requests.end(),
            std::back_inserter(closed_requests));
        if (closed_requests == pair_requests) {
            return partition;
        }
        pair_requests = std::move(closed_requests);
    }
}

FullCircleCellAuthority2
construct_full_circle_cell_authority(
    const Stock2& stock,
    const VerifiedEventPartition2& verified_partition,
    double center_x,
    double center_y,
    double phase_dx,
    double phase_dy,
    double tool_radius,
    double cap_chord_ratio)
{
    if (verified_partition.verdict
        != ContinuousTeaVerdict::CERTIFIED) {
        throw IncompleteFullCircleCellAuthorityError(
            "full-circle cell authority requires a reconstructed event partition");
    }
    const EventPartitionCertificate2& partition =
        verified_partition.partition;
    const std::vector<BoundaryFeatureRecord2>
        boundary_records =
            extract_boundary_records(stock);
    std::vector<std::string> feature_ids;
    feature_ids.reserve(boundary_records.size());
    for (const BoundaryFeatureRecord2& record :
         boundary_records) {
        feature_ids.push_back(record.feature_id);
    }
    if (feature_ids.empty()
        || !std::is_sorted(
            feature_ids.begin(),
            feature_ids.end())) {
        throw IncompleteFullCircleCellAuthorityError(
            "full-circle cell authority requires canonical nonempty stock features");
    }

    bool any_clear = false;
    bool any_material = false;
    bool any_partial = false;
    bool cap_exceeded = false;
    bool unresolved =
        partition.cells.empty()
        || std::any_of(
               partition.overlaps.begin(),
               partition.overlaps.end(),
               overlap_requires_resolution);
    const bool cap_is_pi =
        parse_rational(
            exact_rational_text(
                Epeck::FT(cap_chord_ratio)),
            "circle cap chord ratio")
        == Rational(4);
    std::vector<std::string> cell_records;
    cell_records.reserve(partition.cells.size());
    for (const ParameterCell2& parameter_cell :
         partition.cells) {
        const Rational global_witness(
            Integer(
                parameter_cell.witness_numerator),
            Integer(
                parameter_cell.witness_denominator));
        const ChartWitness witness =
            chart_witness(global_witness);
        const StationEventSource2 source =
            station_source(
                center_x,
                center_y,
                phase_dx,
                phase_dy,
                tool_radius,
                cap_chord_ratio,
                witness);
        const SegmentCellStratum2 cell =
            construct_station_cell_stratum(
                boundary_records,
                source);
        StationCellClassification2 classification =
            classify_station_cell(
                stock,
                source,
                cell.stratum);
        std::vector<std::string>
            cap_projection_ids;
        if (classification.decision
                == StationCellDecision::SAFE_PARTIAL
            || classification.decision
                == StationCellDecision::CAP_EXCEEDED) {
            for (const BranchPairDisposition2& run :
                 classification.material_runs) {
                std::vector<std::string>
                    projection_ids =
                        cap_projections_for(
                            partition,
                            cell,
                            run,
                            CENTER_CHART_IDS[
                                witness.chart],
                            cap_is_pi);
                if (projection_ids.empty()) {
                    classification.decision =
                        StationCellDecision::UNRESOLVED;
                    break;
                }
                cap_projection_ids.insert(
                    cap_projection_ids.end(),
                    std::make_move_iterator(
                        projection_ids.begin()),
                    std::make_move_iterator(
                        projection_ids.end()));
            }
        }
        std::sort(
            cap_projection_ids.begin(),
            cap_projection_ids.end());
        cap_projection_ids.erase(
            std::unique(
                cap_projection_ids.begin(),
                cap_projection_ids.end()),
            cap_projection_ids.end());

        any_clear = any_clear
            || classification.decision
                == StationCellDecision::CLEAR;
        any_material = any_material
            || classification.decision
                == StationCellDecision::MATERIAL;
        any_partial = any_partial
            || classification.decision
                == StationCellDecision::SAFE_PARTIAL
            || classification.decision
                == StationCellDecision::CAP_EXCEEDED;
        cap_exceeded = cap_exceeded
            || classification.decision
                == StationCellDecision::MATERIAL
            || classification.decision
                == StationCellDecision::CAP_EXCEEDED;
        unresolved = unresolved
            || classification.decision
                == StationCellDecision::UNRESOLVED;

        cell_records.push_back(
            encode_string_sequence(
                {
                    "full-circle-cell-decision-v1",
                    parameter_cell.lower_root_id,
                    parameter_cell.upper_root_id,
                    parameter_cell
                        .witness_numerator,
                    parameter_cell
                        .witness_denominator,
                    CENTER_CHART_IDS[
                        witness.chart],
                    rational_text(
                        witness.local_parameter),
                    source.canonical_bytes(),
                    canonical_segment_cell_stratum(
                        cell),
                    decision_text(
                        classification.decision),
                    classification.reference_resolved
                        ? "reference-resolved"
                        : "reference-unresolved",
                    classification.reference_material
                        ? "reference-material"
                        : "reference-clear",
                    encode_string_sequence(
                        cap_projection_ids),
                }));
    }

    const ContinuousTeaVerdict verdict =
        unresolved
        ? ContinuousTeaVerdict::
              UNRESOLVED_DEGENERACY
        : (
              cap_exceeded
                  ? ContinuousTeaVerdict::
                        CAP_EXCEEDED
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
    const auto verdict_text =
        [verdict]() {
            if (verdict
                == ContinuousTeaVerdict::CERTIFIED) {
                return "certified";
            }
            if (verdict
                == ContinuousTeaVerdict::
                    CAP_EXCEEDED) {
                return "cap-exceeded";
            }
            return "unresolved";
        };
    const std::string canonical_bytes =
        encode_string_sequence(
            {
                "full-circle-cell-authority-v1",
                partition.canonical_bytes,
                encode_string_sequence(feature_ids),
                exact_rational_text(
                    Epeck::FT(center_x)),
                exact_rational_text(
                    Epeck::FT(center_y)),
                exact_rational_text(
                    Epeck::FT(phase_dx)),
                exact_rational_text(
                    Epeck::FT(phase_dy)),
                exact_rational_text(
                    Epeck::FT(tool_radius)),
                exact_rational_text(
                    Epeck::FT(cap_chord_ratio)),
                encode_string_sequence(cell_records),
                verdict_text(),
                whole_rim_disposition,
            });
    return {
        verdict,
        whole_rim_disposition,
        canonical_bytes,
        sha256_bytes(canonical_bytes),
    };
}

bool verify_full_circle_cell_authority(
    const Stock2& stock,
    const VerifiedEventPartition2& verified_partition,
    double center_x,
    double center_y,
    double phase_dx,
    double phase_dy,
    double tool_radius,
    double cap_chord_ratio,
    const FullCircleCellAuthority2& candidate)
{
    try {
        const FullCircleCellAuthority2 expected =
            construct_full_circle_cell_authority(
                stock,
                verified_partition,
                center_x,
                center_y,
                phase_dx,
                phase_dy,
                tool_radius,
                cap_chord_ratio);
        return candidate.verdict == expected.verdict
            && candidate.whole_rim_disposition
                == expected.whole_rim_disposition
            && candidate.canonical_bytes
                == expected.canonical_bytes
            && candidate.canonical_digest
                == expected.canonical_digest
            && sha256_bytes(
                   candidate.canonical_bytes)
                == candidate.canonical_digest;
    } catch (const EventSubstrateError&) {
        return false;
    }
}
