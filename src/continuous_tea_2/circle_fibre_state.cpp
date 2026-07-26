#include "circle_fibre_state.h"

#include "event_certificate.h"

#include <algorithm>
#include <string>
#include <vector>

#include <CGAL/CORE/BigRat.h>

namespace {

using Integer = CORE::BigInt;
using Rational = CORE::BigRat;

Rational projection_value(
    const ProjectionRecord2& projection,
    const ParameterCell2& cell)
{
    if (projection.coefficient_rows.size() != 1) {
        throw EventPartitionVerificationError(
            "full-circle event state requires a univariate projection");
    }
    const Rational witness(
        Integer(cell.witness_numerator),
        Integer(cell.witness_denominator));
    Rational value = 0;
    Rational power = 1;
    for (const std::string& coefficient :
         projection.coefficient_rows.front()) {
        value += Rational(Integer(coefficient))
            * power;
        power *= witness;
    }
    if (value == 0) {
        throw EventPartitionVerificationError(
            "full-circle cell witness lies on an event projection");
    }
    return value;
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
    const PartitionEvent2& event)
{
    return {
        encode_string_sequence(
            {
                "full-circle-active-boundary-branch-v1",
                event.feature_id,
                event.support_id,
            }),
        event.support_id,
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
                    && left.support_id == right.support_id;
            });
}

} // namespace

void classify_full_circle_endpoint_fibres(
    EventPartitionCertificate2& certificate)
{
    for (EventFibre2& fibre : certificate.fibres) {
        std::vector<const PartitionEvent2*> endpoints;
        for (const PartitionEvent2& event :
             fibre.events) {
            if (event.kind == "endpoint-order") {
                endpoints.push_back(&event);
            }
        }
        std::vector<std::string> supports;
        for (const PartitionEvent2* event :
             endpoints) {
            supports.push_back(event->support_id);
        }
        std::sort(supports.begin(), supports.end());
        supports.erase(
            std::unique(
                supports.begin(),
                supports.end()),
            supports.end());
        if (supports.size() < 2) {
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
            const ActiveBoundaryBranch2 branch =
                active_branch(*event);
            incident_branches.push_back(branch);
            const ProjectionRecord2& projection =
                governing_projection(
                    certificate,
                    *event);
            if (projection_value(projection, left) < 0) {
                fibre.left_active_branches.push_back(
                    branch);
            }
            if (projection_value(projection, right) < 0) {
                fibre.right_active_branches.push_back(
                    branch);
            }
        }
        canonicalize(incident_branches);
        canonicalize(fibre.left_active_branches);
        canonicalize(fibre.right_active_branches);

        if (fibre.left_active_branches.empty()
            && same_branches(
                fibre.right_active_branches,
                incident_branches)) {
            fibre.ccw_direction = "split";
            fibre.cw_direction = "merge";
        } else if (
            fibre.right_active_branches.empty()
            && same_branches(
                fibre.left_active_branches,
                incident_branches)) {
            fibre.ccw_direction = "merge";
            fibre.cw_direction = "split";
        } else {
            fibre.ccw_direction = "unchanged";
            fibre.cw_direction = "unchanged";
        }
    }
}
