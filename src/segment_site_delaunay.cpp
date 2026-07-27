#include "segment_site_delaunay.h"

#include <algorithm>
#include <cstddef>
#include <set>
#include <string>
#include <utility>

#include <CGAL/enum.h>
#include <CGAL/Kernel/global_functions_2.h>

namespace {

CGAL::Comparison_result compare_points(
    const MatTraits::Point_2& left,
    const MatTraits::Point_2& right)
{
    return CGAL::compare_xy(left, right);
}

bool point_less(
    const MatTraits::Point_2& left,
    const MatTraits::Point_2& right)
{
    return compare_points(left, right)
        == CGAL::SMALLER;
}

std::pair<
    MatTraits::Point_2,
    MatTraits::Point_2>
ordered_segment_endpoints(
    MatTraits::Point_2 first,
    MatTraits::Point_2 second)
{
    if (point_less(second, first)) {
        std::swap(first, second);
    }
    return {
        std::move(first),
        std::move(second),
    };
}

} // namespace

bool CanonicalMatSiteGeometryIndex2::
    segment_key_less(
        const SegmentRecord& left,
        const SegmentRecord& right)
{
    const CGAL::Comparison_result first =
        compare_points(left.first, right.first);
    return first == CGAL::SMALLER
        || (first == CGAL::EQUAL
            && point_less(
                left.second,
                right.second));
}

bool CanonicalMatSiteGeometryIndex2::
    segment_key_equal(
        const SegmentRecord& left,
        const SegmentRecord& right)
{
    return compare_points(
               left.first,
               right.first)
            == CGAL::EQUAL
        && compare_points(
               left.second,
               right.second)
            == CGAL::EQUAL;
}

CanonicalMatSiteGeometryIndex2::
    CanonicalMatSiteGeometryIndex2(
        const CanonicalMatSiteCatalog2& catalog)
{
    points_.reserve(
        catalog.sites().size() / 2);
    segments_.reserve(
        catalog.sites().size() / 2);
    for (const CanonicalMatSite2& site :
         catalog.sites()) {
        if (site.site().is_point()) {
            points_.push_back(
                {
                    site.site().point(),
                    site.stable_id(),
                });
            continue;
        }
        const auto [first, second] =
            ordered_segment_endpoints(
                site.site().source(),
                site.site().target());
        segments_.push_back(
            {
                first,
                second,
                site.stable_id(),
            });
    }

    std::sort(
        points_.begin(),
        points_.end(),
        [](const PointRecord& left,
           const PointRecord& right) {
            return point_less(
                left.point,
                right.point);
        });
    if (std::adjacent_find(
            points_.begin(),
            points_.end(),
            [](const PointRecord& left,
               const PointRecord& right) {
                return compare_points(
                           left.point,
                           right.point)
                    == CGAL::EQUAL;
            })
        != points_.end()) {
        throw DuplicateCanonicalMatSiteGeometryError(
            "canonical MAT point sites have duplicate exact geometry");
    }

    std::sort(
        segments_.begin(),
        segments_.end(),
        segment_key_less);
    if (std::adjacent_find(
            segments_.begin(),
            segments_.end(),
            segment_key_equal)
        != segments_.end()) {
        throw DuplicateCanonicalMatSiteGeometryError(
            "canonical MAT segment sites have duplicate exact geometry");
    }
}

std::size_t
CanonicalMatSiteGeometryIndex2::size() const noexcept
{
    return points_.size() + segments_.size();
}

const std::string&
CanonicalMatSiteGeometryIndex2::stable_id(
    const MatTraits::Site_2& site) const
{
    if (site.is_point()) {
        const MatTraits::Point_2 point =
            site.point();
        const auto found = std::lower_bound(
            points_.begin(),
            points_.end(),
            point,
            [](const PointRecord& record,
               const MatTraits::Point_2& key) {
                return point_less(
                    record.point,
                    key);
            });
        if (found != points_.end()
            && compare_points(
                   found->point,
                   point)
                == CGAL::EQUAL) {
            return found->stable_id;
        }
        throw UnknownCanonicalMatSiteGeometryError(
            "point generator is absent from the canonical MAT site index");
    }

    const auto [first, second] =
        ordered_segment_endpoints(
            site.source(),
            site.target());
    const SegmentRecord key{
        first,
        second,
        {},
    };
    const auto found = std::lower_bound(
        segments_.begin(),
        segments_.end(),
        key,
        segment_key_less);
    if (found != segments_.end()
        && segment_key_equal(*found, key)) {
        return found->stable_id;
    }
    throw UnknownCanonicalMatSiteGeometryError(
        "segment generator is absent from the canonical MAT site index");
}

CanonicalMatDelaunaySource2
CanonicalMatDelaunaySource2::build(
    const CanonicalMatSiteCatalog2& catalog)
{
    return CanonicalMatDelaunaySource2(
        catalog);
}

CanonicalMatDelaunaySource2::
    CanonicalMatDelaunaySource2(
        const CanonicalMatSiteCatalog2& catalog)
    : site_index_(catalog)
{
    for (const CanonicalMatSite2& site :
         catalog.sites()) {
        if (site.site().is_segment()) {
            delaunay_.insert(
                site.site().source(),
                site.site().target());
        }
    }
    if (!delaunay_.is_valid()) {
        throw CanonicalMatDelaunayBijectionError(
            "catalog-fed segment-Delaunay graph is invalid");
    }

    std::set<std::string> matched_ids;
    for (auto vertex =
             delaunay_.finite_vertices_begin();
         vertex
         != delaunay_.finite_vertices_end();
         ++vertex) {
        if (!matched_ids.insert(
                site_index_.stable_id(
                    vertex->site()))
                 .second) {
            throw CanonicalMatDelaunayBijectionError(
                "multiple live generators resolve to one canonical MAT site");
        }
    }
    if (matched_ids.size()
            != site_index_.size()
        || matched_ids.size()
            != catalog.sites().size()) {
        throw CanonicalMatDelaunayBijectionError(
            "canonical MAT site catalog and live generators are not bijective");
    }
}

const SegmentSiteDelaunay2&
CanonicalMatDelaunaySource2::delaunay() const noexcept
{
    return delaunay_;
}

const CanonicalMatSiteGeometryIndex2&
CanonicalMatDelaunaySource2::site_index() const noexcept
{
    return site_index_;
}
