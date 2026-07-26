#include "circle_fibre_state.h"

#include "event_certificate.h"

#include <algorithm>
#include <cstddef>
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
    branches.erase(
        std::unique(
            branches.begin(),
            branches.end(),
            [](const ActiveBoundaryBranch2& first,
               const ActiveBoundaryBranch2& second) {
                return first.branch_id
                    == second.branch_id;
            }),
        branches.end());
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
    const auto append_witnesses =
        [&certificate, &start](
            const EventFibre2& local) {
        const auto root = std::find_if(
            certificate.roots.begin(),
            certificate.roots.end(),
            [&local](
                const AlgebraicRootRecord2& candidate) {
                return candidate.root_id
                    == local.root_id;
            });
        if (root == certificate.roots.end()) {
            throw EventPartitionVerificationError(
                "phase-seam event lacks its local root");
        }
        for (const PartitionEvent2& event :
             local.events) {
            start.local_event_witnesses.push_back(
                {
                    event.kind,
                    event.feature_id,
                    event.support_id,
                    event.trim_id,
                    event.vertex_id,
                    event.branch_id,
                    local.root_id,
                    root->multiplicity,
                });
        }
    };
    const EventFibre2 start_local = start;
    append_witnesses(start_local);
    append_witnesses(end);
    std::sort(
        start.local_event_witnesses.begin(),
        start.local_event_witnesses.end(),
        [](const LocalEventWitness2& first,
           const LocalEventWitness2& second) {
            return std::tie(
                       first.kind,
                       first.feature_id,
                       first.support_id,
                       first.trim_id,
                       first.vertex_id,
                       first.local_root_id,
                       first.local_branch_id)
                < std::tie(
                       second.kind,
                       second.feature_id,
                       second.support_id,
                       second.trim_id,
                       second.vertex_id,
                       second.local_root_id,
                       second.local_branch_id);
        });
    for (std::size_t index = 1;
         index < start.local_event_witnesses.size();
         ++index) {
        const LocalEventWitness2& previous =
            start.local_event_witnesses[index - 1];
        const LocalEventWitness2& current =
            start.local_event_witnesses[index];
        if (std::tie(
                previous.kind,
                previous.feature_id,
                previous.support_id,
                previous.trim_id,
                previous.vertex_id)
                == std::tie(
                    current.kind,
                    current.feature_id,
                    current.support_id,
                    current.trim_id,
                    current.vertex_id)
            && previous.multiplicity
                != current.multiplicity) {
            throw EventPartitionVerificationError(
                "phase-seam local event multiplicities disagree");
        }
    }
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
                event.vertex_id);
        };
    std::sort(
        start.events.begin(),
        start.events.end(),
        [&event_key](
            const PartitionEvent2& first,
            const PartitionEvent2& second) {
            if (event_key(first)
                != event_key(second)) {
                return event_key(first)
                    < event_key(second);
            }
            return first.branch_id
                < second.branch_id;
        });
    for (std::size_t index = 1;
         index < start.events.size();
         ++index) {
        if (event_key(start.events[index - 1])
                == event_key(start.events[index])
            && start.events[index - 1].disposition
                != start.events[index].disposition) {
            throw EventPartitionVerificationError(
                "phase-seam local event dispositions disagree");
        }
    }
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
        std::vector<PartitionEvent2>
            local_endpoint_events;
        if (!fibre.seam_id.empty()) {
            for (const LocalEventWitness2& witness :
                 fibre.local_event_witnesses) {
                if (witness.kind != "endpoint-order") {
                    continue;
                }
                const auto physical = std::find_if(
                    fibre.events.begin(),
                    fibre.events.end(),
                    [&witness](
                        const PartitionEvent2& event) {
                        return event.kind == witness.kind
                            && event.feature_id
                                == witness.feature_id
                            && event.support_id
                                == witness.support_id
                            && event.trim_id
                                == witness.trim_id
                            && event.vertex_id
                                == witness.vertex_id;
                    });
                if (physical == fibre.events.end()) {
                    throw EventPartitionVerificationError(
                        "phase-seam witness lacks its physical incidence");
                }
                local_endpoint_events.push_back(
                    *physical);
                local_endpoint_events.back().branch_id =
                    witness.local_branch_id;
            }
        }
        std::vector<const PartitionEvent2*> endpoints;
        const std::vector<PartitionEvent2>&
            endpoint_source =
                fibre.seam_id.empty()
            ? fibre.events
            : local_endpoint_events;
        for (const PartitionEvent2& event :
             endpoint_source) {
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
        for (const PartitionEvent2* event :
             endpoints) {
            const std::vector<std::string> fields =
                decode_string_sequence(
                    event->branch_id);
            if (fields.size() != 4) {
                throw EventPartitionVerificationError(
                    "full-circle active sheet has malformed event identity");
            }
            std::string sheet_root_id =
                fibre.root_id;
            if (!fibre.seam_id.empty()) {
                const auto witness = std::find_if(
                    fibre.local_event_witnesses.begin(),
                    fibre.local_event_witnesses.end(),
                    [event](
                        const LocalEventWitness2& candidate) {
                        return candidate.local_branch_id
                            == event->branch_id;
                    });
                if (witness
                    == fibre.local_event_witnesses.end()) {
                    throw EventPartitionVerificationError(
                        "phase-seam active sheet lacks local root provenance");
                }
                sheet_root_id =
                    witness->local_root_id;
            }
            const ActiveBoundaryBranch2 branch =
                active_branch(
                    *event,
                    sheet_root_id,
                    0);
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

        const auto physical_keys =
            [](const std::vector<ActiveBoundaryBranch2>&
                   branches) {
                std::vector<std::string> result;
                result.reserve(branches.size());
                for (const ActiveBoundaryBranch2& branch :
                     branches) {
                    result.push_back(
                        encode_string_sequence(
                            {
                                "full-circle-physical-incidence-v1",
                                branch.feature_id,
                                branch.support_id,
                                branch.trim_id,
                                std::to_string(
                                    branch.sheet_ordinal),
                            }));
                }
                std::sort(
                    result.begin(),
                    result.end());
                result.erase(
                    std::unique(
                        result.begin(),
                        result.end()),
                    result.end());
                return result;
            };
        const std::vector<std::string> left_keys =
            physical_keys(
                fibre.left_active_branches);
        const std::vector<std::string> right_keys =
            physical_keys(
                fibre.right_active_branches);
        const std::vector<std::string> incident_keys =
            physical_keys(incident_branches);
        const auto contains =
            [](const std::vector<std::string>& whole,
               const std::vector<std::string>& part) {
                return std::includes(
                    whole.begin(),
                    whole.end(),
                    part.begin(),
                    part.end());
            };
        std::vector<std::string> active_union =
            left_keys;
        active_union.insert(
            active_union.end(),
            right_keys.begin(),
            right_keys.end());
        std::sort(
            active_union.begin(),
            active_union.end());
        active_union.erase(
            std::unique(
                active_union.begin(),
                active_union.end()),
            active_union.end());
        if (left_keys != right_keys
            && contains(right_keys, left_keys)
            && contains(incident_keys, right_keys)) {
            fibre.ccw_direction = "split";
            fibre.cw_direction = "merge";
        } else if (
            left_keys != right_keys
            && contains(left_keys, right_keys)
            && contains(incident_keys, left_keys)) {
            fibre.ccw_direction = "merge";
            fibre.cw_direction = "split";
        } else if (
            left_keys == right_keys
            && contains(incident_keys, left_keys)) {
            fibre.ccw_direction = "unchanged";
            fibre.cw_direction = "unchanged";
        } else if (
            left_keys.size() == right_keys.size()
            && active_union == incident_keys) {
            fibre.ccw_direction = "unchanged";
            fibre.cw_direction = "unchanged";
        } else {
            throw EventPartitionVerificationError(
                "full-circle active sheets do not prove a cyclic merge or split");
        }
    }
}
