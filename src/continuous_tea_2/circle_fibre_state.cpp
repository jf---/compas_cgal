#include "circle_fibre_state.h"

#include "event_certificate.h"

#include <algorithm>
#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include <CGAL/CORE/BigRat.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/number_utils.h>

namespace {

using Integer = CORE::BigInt;
using Polynomial = CGAL::Polynomial<Integer>;

CGAL::Sign signed_projection_sign(
    const ProjectionRecord2& projection,
    const ParameterCell2& cell)
{
    if (projection
            .signed_predicate_coefficients.empty()) {
        throw EventPartitionVerificationError(
            "full-circle event state lacks its signed radial predicate");
    }
    std::vector<Integer> coefficients;
    coefficients.reserve(
        projection
            .signed_predicate_coefficients.size());
    for (const std::string& coefficient :
         projection.signed_predicate_coefficients) {
        coefficients.emplace_back(coefficient);
    }
    using Traits =
        CGAL::Polynomial_traits_d<Polynomial>;
    const Polynomial predicate =
        typename Traits::Construct_polynomial()(
            coefficients.begin(),
            coefficients.end());
    const Integer value =
        typename Traits::Evaluate_homogeneous()(
            predicate,
            Integer(cell.witness_numerator),
            Integer(cell.witness_denominator));
    const CGAL::Sign sign = CGAL::sign(value);
    if (sign == CGAL::ZERO) {
        throw EventPartitionVerificationError(
            "full-circle cell witness lies on an event projection");
    }
    return sign;
}

const ParameterCell2& adjacent_cell(
    const EventPartitionCertificate2& certificate,
    const std::string& root_id,
    bool left)
{
    const auto found = std::find_if(
        certificate.cells.begin(),
        certificate.cells.end(),
        [&root_id, left](
            const ParameterCell2& cell) {
            return left
                ? cell.upper_root_id == root_id
                : cell.lower_root_id == root_id;
        });
    if (found != certificate.cells.end()) {
        return *found;
    }
    if (certificate.cells.empty()) {
        throw EventPartitionVerificationError(
            "full-circle event fibre lacks adjacent cells");
    }
    return left
        ? certificate.cells.back()
        : certificate.cells.front();
}

const ProjectionRecord2& governing_projection(
    const EventPartitionCertificate2& certificate,
    const PartitionEvent2& event)
{
    const std::vector<std::string> fields =
        decode_string_sequence(event.branch_id);
    if (fields.size() != 4
        || fields[0]
            != "full-circle-vertex-passage-v1") {
        throw EventPartitionVerificationError(
            "full-circle endpoint event has malformed branch identity");
    }
    const std::string projection_id =
        encode_string_sequence(
            {
                fields[0],
                fields[1],
                fields[2],
            });
    const auto found = std::find_if(
        certificate.projections.begin(),
        certificate.projections.end(),
        [&projection_id](
            const ProjectionRecord2& projection) {
            return projection.projection_id
                == projection_id;
        });
    if (found == certificate.projections.end()) {
        throw EventPartitionVerificationError(
            "full-circle endpoint event lacks its governing projection");
    }
    return *found;
}

ActiveBoundaryBranch2 active_branch(
    const PartitionEvent2& event,
    const std::string& root_id,
    std::size_t sheet_ordinal)
{
    const std::vector<std::string> fields =
        decode_string_sequence(event.branch_id);
    if (fields.size() != 4) {
        throw EventPartitionVerificationError(
            "full-circle active sheet has malformed event identity");
    }
    return {
        encode_string_sequence(
            {
                "full-circle-active-boundary-sheet-v2",
                event.feature_id,
                event.support_id,
                event.trim_id,
                fields[2],
                std::to_string(sheet_ordinal),
                root_id,
            }),
        event.feature_id,
        event.support_id,
        event.trim_id,
        fields[2],
        sheet_ordinal,
        root_id,
    };
}

void canonicalize(
    std::vector<ActiveBoundaryBranch2>& branches)
{
    std::sort(
        branches.begin(),
        branches.end(),
        [](const ActiveBoundaryBranch2& first,
           const ActiveBoundaryBranch2& second) {
            return first.branch_id
                < second.branch_id;
        });
}

bool same_branches(
    const std::vector<ActiveBoundaryBranch2>& first,
    const std::vector<ActiveBoundaryBranch2>& second)
{
    return first.size() == second.size()
        && std::equal(
            first.begin(),
            first.end(),
            second.begin(),
            [](const ActiveBoundaryBranch2& left,
               const ActiveBoundaryBranch2& right) {
                return left.branch_id == right.branch_id
                    && left.feature_id
                        == right.feature_id
                    && left.support_id == right.support_id;
            });
}

void coalesce_phase_seam(
    EventPartitionCertificate2& certificate)
{
    std::optional<std::size_t> start_index;
    std::optional<std::size_t> end_index;
    for (std::size_t index = 0;
         index < certificate.fibres.size();
         ++index) {
        const std::string& root_id =
            certificate.fibres[index].root_id;
        const bool has_left = std::any_of(
            certificate.cells.begin(),
            certificate.cells.end(),
            [&root_id](const ParameterCell2& cell) {
                return cell.upper_root_id == root_id;
            });
        const bool has_right = std::any_of(
            certificate.cells.begin(),
            certificate.cells.end(),
            [&root_id](const ParameterCell2& cell) {
                return cell.lower_root_id == root_id;
            });
        if (!has_left) {
            start_index = index;
        }
        if (!has_right) {
            end_index = index;
        }
    }
    if (!start_index.has_value()
        || !end_index.has_value()
        || start_index == end_index) {
        return;
    }
    EventFibre2& start =
        certificate.fibres[*start_index];
    EventFibre2& end =
        certificate.fibres[*end_index];
    start.local_root_ids = {
        start.root_id,
        end.root_id,
    };
    const auto anchor = std::find_if(
        certificate.charts.begin(),
        certificate.charts.end(),
        [](const ParameterChart2& chart) {
            return chart.chart_id
                == "center-quarter-0-v1";
        });
    if (anchor == certificate.charts.end()
        || !anchor->owns_start_seam) {
        throw EventPartitionVerificationError(
            "full-circle phase seam lacks its anchor chart");
    }
    start.seam_id = anchor->start_seam_id;
    start.events.insert(
        start.events.end(),
        end.events.begin(),
        end.events.end());
    const auto event_key =
        [](const PartitionEvent2& event) {
            return std::tie(
                event.kind,
                event.feature_id,
                event.support_id,
                event.trim_id,
                event.vertex_id,
                event.branch_id,
                event.disposition);
        };
    std::sort(
        start.events.begin(),
        start.events.end(),
        [&event_key](
            const PartitionEvent2& first,
            const PartitionEvent2& second) {
            return event_key(first)
                < event_key(second);
        });
    start.events.erase(
        std::unique(
            start.events.begin(),
            start.events.end(),
            [&event_key](
                const PartitionEvent2& first,
                const PartitionEvent2& second) {
                return event_key(first)
                    == event_key(second);
            }),
        start.events.end());
    certificate.fibres.erase(
        certificate.fibres.begin()
            + static_cast<std::ptrdiff_t>(
                *end_index));
}

} // namespace

void classify_full_circle_endpoint_fibres(
    EventPartitionCertificate2& certificate)
{
    coalesce_phase_seam(certificate);
    for (EventFibre2& fibre : certificate.fibres) {
        std::vector<const PartitionEvent2*> endpoints;
        for (const PartitionEvent2& event :
             fibre.events) {
            if (event.kind == "endpoint-order") {
                endpoints.push_back(&event);
            }
        }
        std::sort(
            endpoints.begin(),
            endpoints.end(),
            [](const PartitionEvent2* first,
               const PartitionEvent2* second) {
                return std::tie(
                           first->feature_id,
                           first->support_id,
                           first->trim_id,
                           first->branch_id)
                    < std::tie(
                           second->feature_id,
                           second->support_id,
                           second->trim_id,
                           second->branch_id);
            });
        if (endpoints.size() < 2) {
            continue;
        }

        const ParameterCell2& left =
            adjacent_cell(
                certificate,
                fibre.root_id,
                true);
        const ParameterCell2& right =
            adjacent_cell(
                certificate,
                fibre.root_id,
                false);
        std::vector<ActiveBoundaryBranch2>
            incident_branches;
        incident_branches.reserve(endpoints.size());
        std::map<std::string, std::size_t>
            sheet_ordinals;
        for (const PartitionEvent2* event :
             endpoints) {
            const std::vector<std::string> fields =
                decode_string_sequence(
                    event->branch_id);
            if (fields.size() != 4) {
                throw EventPartitionVerificationError(
                    "full-circle active sheet has malformed event identity");
            }
            const std::string sheet_key =
                encode_string_sequence(
                    {
                        event->feature_id,
                        event->support_id,
                        event->trim_id,
                        fields[2],
                    });
            const std::size_t sheet_ordinal =
                sheet_ordinals[sheet_key]++;
            const ActiveBoundaryBranch2 branch =
                active_branch(
                    *event,
                    fibre.root_id,
                    sheet_ordinal);
            incident_branches.push_back(branch);
            const ProjectionRecord2& projection =
                governing_projection(
                    certificate,
                    *event);
            const bool cyclic_seam =
                fibre.local_root_ids.size() == 2;
            if ((!cyclic_seam
                    || fields[2]
                        == "center-quarter-3-v1")
                && signed_projection_sign(
                    projection,
                    left)
                == CGAL::NEGATIVE) {
                fibre.left_active_branches.push_back(
                    branch);
            }
            if ((!cyclic_seam
                    || fields[2]
                        == "center-quarter-0-v1")
                && signed_projection_sign(
                    projection,
                    right)
                == CGAL::NEGATIVE) {
                fibre.right_active_branches.push_back(
                    branch);
            }
        }
        canonicalize(incident_branches);
        canonicalize(fibre.left_active_branches);
        canonicalize(fibre.right_active_branches);

        const auto contained_in_incidents =
            [&incident_branches](
                const std::vector<ActiveBoundaryBranch2>&
                    active) {
                return std::all_of(
                    active.begin(),
                    active.end(),
                    [&incident_branches](
                        const ActiveBoundaryBranch2& branch) {
                        return std::any_of(
                            incident_branches.begin(),
                            incident_branches.end(),
                            [&branch](
                                const ActiveBoundaryBranch2&
                                    incident) {
                                return branch.branch_id
                                    == incident.branch_id;
                            });
                    });
            };
        const bool cyclic_seam =
            fibre.local_root_ids.size() == 2;
        if (
            (
                fibre.left_active_branches.empty()
                && same_branches(
                    fibre.right_active_branches,
                    incident_branches))
            || (
                cyclic_seam
                && fibre.left_active_branches.empty()
                && !fibre.right_active_branches.empty()
                && contained_in_incidents(
                    fibre.right_active_branches))) {
            fibre.ccw_direction = "split";
            fibre.cw_direction = "merge";
        } else if (
            (
                fibre.right_active_branches.empty()
                && same_branches(
                    fibre.left_active_branches,
                    incident_branches))
            || (
                cyclic_seam
                && fibre.right_active_branches.empty()
                && !fibre.left_active_branches.empty()
                && contained_in_incidents(
                    fibre.left_active_branches))) {
            fibre.ccw_direction = "merge";
            fibre.cw_direction = "split";
        } else {
            fibre.ccw_direction = "unchanged";
            fibre.cw_direction = "unchanged";
        }
    }
}
