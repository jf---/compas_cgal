#include "station_classifier.h"

#include "../stock_2.h"

#include <algorithm>
#include <string_view>
#include <type_traits>
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
            throw InvalidStationSourceError(
                std::string(role)
                + " has zero denominator");
        }
        return Rational(numerator, denominator);
    } catch (const EventSubstrateError&) {
        throw;
    } catch (const std::exception&) {
        throw InvalidStationSourceError(
            std::string(role)
            + " is not one exact rational");
    }
}

Epeck::FT exact_ft(const Rational& value)
{
    static_assert(
        std::is_same_v<Rational, CGAL::Epeck_ft>);
    return Epeck::FT(value);
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

} // namespace

StationCellClassification2 classify_station_cell(
    const Stock2& stock,
    const StationEventSource2& source,
    const SegmentEventStratum2& stratum)
{
    if (stratum.kind != "cell") {
        return {
            StationCellDecision::UNRESOLVED,
            false,
            false,
            {},
        };
    }
    const Epeck::FT center_x = exact_ft(
        parse_rational(
            source.center_x().text(),
            "station center x"));
    const Epeck::FT center_y = exact_ft(
        parse_rational(
            source.center_y().text(),
            "station center y"));
    const Epeck::FT radius = exact_ft(
        parse_rational(
            source.tool_radius().text(),
            "station tool radius"));
    const GpsPoint reference(
        center_x,
        center_y - radius);
    const CGAL::Oriented_side side =
        stock.set().oriented_side(reference);
    const bool reference_resolved =
        side != CGAL::ON_ORIENTED_BOUNDARY;
    const bool reference_material =
        side == CGAL::ON_POSITIVE_SIDE;
    if (!reference_resolved) {
        return {
            StationCellDecision::UNRESOLVED,
            reference_material,
            false,
            {},
        };
    }

    const std::size_t count =
        stratum.active_branch_ids.size();
    if (count == 0) {
        return {
            reference_material
                ? StationCellDecision::MATERIAL
                : StationCellDecision::CLEAR,
            reference_material,
            true,
            {},
        };
    }
    if (count % 2 != 0) {
        return {
            StationCellDecision::UNRESOLVED,
            reference_material,
            true,
            {},
        };
    }

    std::vector<std::pair<std::size_t, std::size_t>>
        run_indices;
    if (reference_material) {
        run_indices.emplace_back(count - 1, 0);
        for (std::size_t index = 1;
             index + 1 < count;
             index += 2) {
            run_indices.emplace_back(
                index,
                index + 1);
        }
    } else {
        for (std::size_t index = 0;
             index + 1 < count;
             index += 2) {
            run_indices.emplace_back(
                index,
                index + 1);
        }
    }

    std::vector<BranchPairDisposition2>
        material_runs;
    material_runs.reserve(run_indices.size());
    StationCellDecision decision =
        StationCellDecision::SAFE_PARTIAL;
    for (const auto [first_index, second_index] :
         run_indices) {
        const BranchPairDisposition2* pair =
            pair_for(
                stratum,
                stratum.active_branch_ids[first_index],
                stratum.active_branch_ids[second_index]);
        if (pair == nullptr
            || (
                pair->cap_disposition != "below-cap"
                && pair->cap_disposition != "equal-cap"
                && pair->cap_disposition != "above-cap")) {
            return {
                StationCellDecision::UNRESOLVED,
                reference_material,
                true,
                std::move(material_runs),
            };
        }
        material_runs.push_back(*pair);
        if (pair->cap_disposition == "above-cap") {
            decision =
                StationCellDecision::CAP_EXCEEDED;
        }
    }
    return {
        decision,
        reference_material,
        true,
        std::move(material_runs),
    };
}
