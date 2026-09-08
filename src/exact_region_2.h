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

// A set was offered to ExactRegion2::build whose arrangement reads geometry
// traits that no live ReachSet owns. There is no repair at this point: the
// caller has to build the set on traits something keeps alive.
class ExactRegionTraitsUnownedError
    : public ReachableDomainConstructionError {
public:
    using ReachableDomainConstructionError::ReachableDomainConstructionError;
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
    // Adopt `set` as this region's exact storage. The set is taken, never
    // copied: the region reads the very arrangement the caller built.
    //
    // `traits_owner_candidates` names the sets whose geometry traits `set`'s
    // arrangement may borrow -- pass each source region's traits_owner() when
    // `set` was seeded by copy from another region's set, and an empty list when
    // `set` was built from scratch and owns its own traits. `build` resolves
    // which one the arrangement actually reads and keeps it alive; a set that
    // matches none is rejected with ExactRegionTraitsUnownedError rather than
    // adopted with a pointer into freed memory.
    static ExactRegion2 build(
        std::shared_ptr<const ReachSet> set,
        ExactRegionRole2 role,
        std::string recipe_record,
        std::vector<std::shared_ptr<const ReachSet>> traits_owner_candidates);
    ExactRegion2 clone() const;
    bool contains(double x, double y) const;
    bool is_empty() const;
    std::size_t component_count() const;
    bool is_subset_of(const ExactRegion2& other) const;
    bool exactly_equals(const ExactRegion2& other) const;
    bool shares_storage_with_for_audit(const ExactRegion2& other) const;
    // True when the geometry traits this region's arrangement reads belong to a
    // ReachSet this region keeps alive.
    //
    // CGAL's Gps copy constructor gives the copy a fresh Traits_2 of its own but
    // builds the copy's arrangement as Aos_2(*(ps.m_arr)), and
    // Arrangement_on_surface_2::assign propagates a BORROWED traits pointer
    // verbatim (m_geom_traits = arr.m_own_traits ? new Traits_adaptor_2
    // : arr.m_geom_traits). Every Gps builds its arrangement in borrow mode, so
    // a copied set reads the traits of the ROOT set it was copied from -- and
    // its own freshly allocated traits is never used. This accessor states the
    // lifetime the object graph assumes: whoever owns that root must outlive
    // this region.
    bool arrangement_traits_are_owned_for_audit() const;
    const ReachSet& set() const;
    // The set that owns the geometry traits this region's arrangement reads --
    // this region's own set when that set is a root. Hand it to build() as a
    // candidate whenever a new set is seeded by copy from set().
    const std::shared_ptr<const ReachSet>& traits_owner() const;
    ExactRegionRole2 role() const;
    const std::string& recipe_record() const;

private:
    ExactRegion2(
        std::shared_ptr<const ReachSet> traits_owner,
        std::shared_ptr<const ReachSet> set,
        ExactRegionRole2 role,
        std::string recipe_record);

    // Declared BEFORE set_ so it is destroyed AFTER the arrangement that reads
    // it: members are destroyed in reverse declaration order. Aliases set_ when
    // set_ is a root.
    std::shared_ptr<const ReachSet> traits_owner_;
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
