#include "exact_region_2.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <utility>

namespace {

void append_u64(std::string& target, std::uint64_t value)
{
    for (int shift = 56; shift >= 0; shift -= 8) {
        target.push_back(static_cast<char>((value >> shift) & 0xffU));
    }
}

} // namespace

std::string reach_u64_record(std::size_t value)
{
    std::string result;
    append_u64(result, static_cast<std::uint64_t>(value));
    return result;
}

std::string reach_length_prefixed(const std::string& value)
{
    std::string result = reach_u64_record(value.size());
    result += value;
    return result;
}

std::string reach_tagged_record(
    std::string_view tag,
    const std::vector<std::string>& fields)
{
    std::string result(tag);
    result.push_back('\0');
    result += reach_u64_record(fields.size());
    for (const std::string& field : fields) {
        result += reach_length_prefixed(field);
    }
    return result;
}

std::string reach_binary64_record(double value)
{
    if (!std::isfinite(value)) {
        throw ReachableDomainConstructionError(
            "structural binary64 records require finite values.");
    }
    const double canonical = value == 0.0 ? 0.0 : value;
    const std::uint64_t bits = std::bit_cast<std::uint64_t>(canonical);
    std::string result;
    append_u64(result, bits);
    return result;
}

bool reach_exact_subset(const ReachSet& subset, const ReachSet& superset)
{
    ReachSet difference(subset);
    difference.difference(superset);
    return difference.is_empty();
}

bool reach_exact_equal(const ReachSet& left, const ReachSet& right)
{
    return reach_exact_subset(left, right)
        && reach_exact_subset(right, left);
}

std::size_t reach_component_count(const ReachSet& set)
{
    std::vector<ReachPolygonWithHoles> components;
    set.polygons_with_holes(std::back_inserter(components));
    return components.size();
}

std::vector<std::string> reach_component_records(
    const ReachSet& set,
    const std::string& record_tag)
{
    std::vector<ReachPolygonWithHoles> components;
    set.polygons_with_holes(std::back_inserter(components));
    struct OrderedComponent {
        ReachPoint minimum;
        std::size_t outer_curve_count;
        std::size_t hole_count;
    };
    std::vector<OrderedComponent> ordered;
    ordered.reserve(components.size());
    for (const ReachPolygonWithHoles& component : components) {
        const ReachPolygon& outer = component.outer_boundary();
        if (outer.curves_begin() == outer.curves_end()) {
            throw ReachableDomainConstructionError(
                "exact component has an empty outer boundary.");
        }
        ReachPoint minimum = outer.curves_begin()->source();
        for (auto curve = outer.curves_begin();
             curve != outer.curves_end();
             ++curve) {
            if (ReachTraits().compare_xy_2_object()(
                    curve->source(),
                    minimum)
                == CGAL::SMALLER) {
                minimum = curve->source();
            }
            if (ReachTraits().compare_xy_2_object()(
                    curve->target(),
                    minimum)
                == CGAL::SMALLER) {
                minimum = curve->target();
            }
        }
        ordered.push_back({
            minimum,
            outer.size(),
            static_cast<std::size_t>(
                std::distance(
                    component.holes_begin(),
                    component.holes_end())),
        });
    }
    std::sort(
        ordered.begin(),
        ordered.end(),
        [](const OrderedComponent& left, const OrderedComponent& right) {
            return ReachTraits().compare_xy_2_object()(
                       left.minimum,
                       right.minimum)
                == CGAL::SMALLER;
        });
    std::vector<std::string> records;
    records.reserve(ordered.size());
    for (std::size_t index = 0; index < ordered.size(); ++index) {
        records.push_back(reach_tagged_record(
            record_tag,
            {
                reach_u64_record(index),
                reach_u64_record(ordered[index].outer_curve_count),
                reach_u64_record(ordered[index].hole_count),
            }));
    }
    return records;
}

ExactRegion2 ExactRegion2::build(
    ReachSet set,
    ExactRegionRole2 role,
    std::string recipe_record)
{
    return ExactRegion2(
        std::make_shared<const ReachSet>(std::move(set)),
        role,
        std::move(recipe_record));
}

ExactRegion2::ExactRegion2(
    std::shared_ptr<const ReachSet> set,
    ExactRegionRole2 role,
    std::string recipe_record)
    : set_(std::move(set))
    , role_(role)
    , recipe_record_(std::move(recipe_record))
{
    if (!set_) {
        throw ReachableDomainConstructionError(
            "exact region requires owned native storage.");
    }
}

ExactRegion2 ExactRegion2::clone() const
{
    return *this;
}

bool ExactRegion2::contains(double x, double y) const
{
    if (!std::isfinite(x) || !std::isfinite(y)) {
        throw ReachableDomainConstructionError(
            "exact region membership requires finite binary64 coordinates.");
    }
    return set_->oriented_side(ReachPoint(ReachFT(x), ReachFT(y)))
        != CGAL::ON_NEGATIVE_SIDE;
}

bool ExactRegion2::is_empty() const
{
    return set_->is_empty();
}

std::size_t ExactRegion2::component_count() const
{
    return reach_component_count(*set_);
}

bool ExactRegion2::is_subset_of(const ExactRegion2& other) const
{
    return reach_exact_subset(*set_, *other.set_);
}

bool ExactRegion2::exactly_equals(const ExactRegion2& other) const
{
    return reach_exact_equal(*set_, *other.set_);
}

bool ExactRegion2::shares_storage_with_for_audit(
    const ExactRegion2& other) const
{
    return set_.get() == other.set_.get();
}

const ReachSet& ExactRegion2::set() const
{
    return *set_;
}

ExactRegionRole2 ExactRegion2::role() const
{
    return role_;
}

const std::string& ExactRegion2::recipe_record() const
{
    return recipe_record_;
}
