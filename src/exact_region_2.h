#pragma once

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Gps_circle_segment_traits_2.h>

using ReachKernel =
    CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;
using ReachTraits = CGAL::Gps_circle_segment_traits_2<ReachKernel>;
using ReachSet = CGAL::General_polygon_set_2<ReachTraits>;
using ReachPolygon = ReachTraits::Polygon_2;
using ReachPolygonWithHoles = ReachTraits::Polygon_with_holes_2;
using ReachCurve = ReachTraits::Curve_2;
using ReachXCurve = ReachTraits::X_monotone_curve_2;
using ReachPoint = ReachTraits::Point_2;
using ReachKernelPoint = ReachKernel::Point_2;
using ReachKernelVector = ReachKernel::Vector_2;
using ReachFT = ReachKernel::FT;

class ReachableDomainConstructionError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

enum class ExactRegionRole2 {
    Design,
    CenterDomain,
    ReachableMaterial,
    UnreachableResidual,
    AccumulatedSweeps,
    CoverageResidual,
};

class ExactRegion2 {
public:
    static ExactRegion2 build(
        ReachSet set,
        ExactRegionRole2 role,
        std::string recipe_record);
    ExactRegion2 clone() const;
    bool contains(double x, double y) const;
    bool is_empty() const;
    std::size_t component_count() const;
    bool is_subset_of(const ExactRegion2& other) const;
    bool exactly_equals(const ExactRegion2& other) const;
    bool shares_storage_with_for_audit(const ExactRegion2& other) const;
    const ReachSet& set() const;
    ExactRegionRole2 role() const;
    const std::string& recipe_record() const;

private:
    ExactRegion2(
        std::shared_ptr<const ReachSet> set,
        ExactRegionRole2 role,
        std::string recipe_record);

    std::shared_ptr<const ReachSet> set_;
    ExactRegionRole2 role_;
    std::string recipe_record_;
};

bool reach_exact_subset(const ReachSet& subset, const ReachSet& superset);
bool reach_exact_equal(const ReachSet& left, const ReachSet& right);
std::size_t reach_component_count(const ReachSet& set);
std::vector<std::string> reach_component_records(
    const ReachSet& set,
    const std::string& record_tag);

std::string reach_tagged_record(
    std::string_view tag,
    const std::vector<std::string>& fields);
std::string reach_binary64_record(double value);
std::string reach_length_prefixed(const std::string& value);
std::string reach_u64_record(std::size_t value);
