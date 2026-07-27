#include "segment_site_catalog.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

static_assert(std::is_same_v<ReachKernel, MatKernel>);
static_assert(
    std::is_same_v<
        ReachKernelPoint,
        MatTraits::Point_2>);

std::int64_t checked_index(
    const std::size_t value,
    const char* name)
{
    if (value
        > static_cast<std::uintmax_t>(
            std::numeric_limits<std::int64_t>::max())) {
        throw MatSiteCatalogOverflowError(
            std::string(name)
            + " exceeds int64 range");
    }
    return static_cast<std::int64_t>(value);
}

std::int64_t global_ring_index(
    const CanonicalReachRing2& ring)
{
    if (ring.outer) {
        return 0;
    }
    if (ring.canonical_ordinal
        >= static_cast<std::uintmax_t>(
            std::numeric_limits<std::int64_t>::max())) {
        throw MatSiteCatalogOverflowError(
            "canonical hole ordinal exceeds MAT ring range");
    }
    return static_cast<std::int64_t>(
        ring.canonical_ordinal + 1);
}

std::string stable_site_identity(
    const MatSiteProvenance2& provenance)
{
    return reach_tagged_record(
        "mat-site-v1",
        {
            reach_u64_record(
                static_cast<std::size_t>(
                    provenance.kind)),
            reach_u64_record(
                static_cast<std::size_t>(
                    provenance.ring)),
            reach_u64_record(
                static_cast<std::size_t>(
                    provenance.feature)),
        });
}

void require_canonical_catalog(
    const std::vector<CanonicalMatSite2>& sites)
{
    if (sites.empty()) {
        throw InvalidCanonicalMatSiteCatalogError(
            "canonical MAT site catalog is empty");
    }
    for (std::size_t index = 0;
         index < sites.size();
         ++index) {
        const CanonicalMatSite2& site =
            sites[index];
        if (site.stable_id().empty()
            || (index > 0
                && site.stable_id()
                    <= sites[index - 1].stable_id())) {
            throw InvalidCanonicalMatSiteCatalogError(
                "canonical MAT site identities are not strictly ordered");
        }
        if ((site.provenance().kind
                    == MatSiteKind2::Point)
                != site.site().is_point()
            || (site.provenance().kind
                    == MatSiteKind2::Point)
                == site.endpoint_site_ids().has_value()) {
            throw InvalidCanonicalMatSiteCatalogError(
                "canonical MAT site kind contradicts exact geometry");
        }
    }
}

} // namespace

CanonicalMatSite2
CanonicalMatSite2::build_point(
    const std::int64_t ring,
    const std::int64_t feature,
    const MatTraits::Point_2& point)
{
    if (ring < 0 || feature < 0) {
        throw InvalidCanonicalMatSiteCatalogError(
            "MAT point provenance is negative");
    }
    const MatSiteProvenance2 provenance{
        MatSiteKind2::Point,
        ring,
        feature,
    };
    return CanonicalMatSite2(
        stable_site_identity(provenance),
        provenance,
        MatTraits::Site_2::construct_site_2(
            point),
        std::nullopt);
}

CanonicalMatSite2
CanonicalMatSite2::build_open_segment(
    const std::int64_t ring,
    const std::int64_t feature,
    const std::int64_t target_feature,
    const MatTraits::Point_2& source,
    const MatTraits::Point_2& target)
{
    if (ring < 0
        || feature < 0
        || target_feature < 0) {
        throw InvalidCanonicalMatSiteCatalogError(
            "MAT segment provenance is negative");
    }
    if (source == target) {
        throw InvalidCanonicalMatSiteCatalogError(
            "MAT segment site endpoints coincide");
    }
    const MatSiteProvenance2 provenance{
        MatSiteKind2::OpenSegment,
        ring,
        feature,
    };
    const MatSiteProvenance2 source_provenance{
        MatSiteKind2::Point,
        ring,
        feature,
    };
    const MatSiteProvenance2 target_provenance{
        MatSiteKind2::Point,
        ring,
        target_feature,
    };
    return CanonicalMatSite2(
        stable_site_identity(provenance),
        provenance,
        MatTraits::Site_2::construct_site_2(
            source,
            target),
        std::array<std::string, 2>{
            stable_site_identity(
                source_provenance),
            stable_site_identity(
                target_provenance),
        });
}

CanonicalMatSite2::CanonicalMatSite2(
    std::string stable_id,
    MatSiteProvenance2 provenance,
    MatTraits::Site_2 site,
    std::optional<std::array<std::string, 2>>
        endpoint_site_ids)
    : stable_id_(std::move(stable_id)),
      provenance_(provenance),
      site_(std::move(site)),
      endpoint_site_ids_(
          std::move(endpoint_site_ids))
{
}

const std::string&
CanonicalMatSite2::stable_id() const noexcept
{
    return stable_id_;
}

const MatSiteProvenance2&
CanonicalMatSite2::provenance() const noexcept
{
    return provenance_;
}

const MatTraits::Site_2&
CanonicalMatSite2::site() const noexcept
{
    return site_;
}

const std::optional<std::array<std::string, 2>>&
CanonicalMatSite2::endpoint_site_ids() const noexcept
{
    return endpoint_site_ids_;
}

CanonicalMatSiteCatalog2::CanonicalMatSiteCatalog2(
    std::vector<CanonicalMatSite2> sites)
    : sites_(std::move(sites))
{
}

const std::vector<CanonicalMatSite2>&
CanonicalMatSiteCatalog2::sites() const noexcept
{
    return sites_;
}

std::vector<std::array<std::int64_t, 3>>
CanonicalMatSiteCatalog2::site_provenance() const
{
    std::vector<std::array<std::int64_t, 3>>
        provenance;
    provenance.reserve(sites_.size());
    for (const CanonicalMatSite2& site : sites_) {
        provenance.push_back(
            {
                static_cast<std::int64_t>(
                    site.provenance().kind),
                site.provenance().ring,
                site.provenance().feature,
            });
    }
    return provenance;
}

std::int64_t CanonicalMatSiteCatalog2::index_of(
    const std::string& stable_id) const
{
    const auto found = std::lower_bound(
        sites_.begin(),
        sites_.end(),
        stable_id,
        [](const CanonicalMatSite2& site,
           const std::string& identity) {
            return site.stable_id() < identity;
        });
    if (found == sites_.end()
        || found->stable_id() != stable_id) {
        throw UnknownMatSiteIdentityError(
            "MAT graph references a site outside its canonical catalog");
    }
    return static_cast<std::int64_t>(
        std::distance(
            sites_.begin(),
            found));
}

CanonicalMatSiteCatalog2
canonical_mat_site_catalog(
    const CanonicalReachInput2& input)
{
    validate_canonical_reach_input(input);
    std::size_t point_count =
        input.outer.points.size();
    for (const CanonicalReachRing2& hole :
         input.holes) {
        if (point_count
            > std::numeric_limits<std::size_t>::max()
                - hole.points.size()) {
            throw MatSiteCatalogOverflowError(
                "MAT site count overflows size_t");
        }
        point_count += hole.points.size();
    }
    const std::uintmax_t maximum_sites =
        static_cast<std::uintmax_t>(
            std::numeric_limits<std::int64_t>::max());
    if (static_cast<std::uintmax_t>(point_count)
        > maximum_sites / 2) {
        throw MatSiteCatalogOverflowError(
            "MAT site catalog exceeds int64 range");
    }

    std::vector<CanonicalMatSite2> sites;
    sites.reserve(point_count * 2);
    const auto append_points =
        [&sites](
            const CanonicalReachRing2& ring) {
            const std::int64_t ring_index =
                global_ring_index(ring);
            for (std::size_t feature = 0;
                 feature < ring.points.size();
                 ++feature) {
                sites.push_back(
                    CanonicalMatSite2::build_point(
                        ring_index,
                        checked_index(
                            feature,
                            "MAT point feature ordinal"),
                        ring.points[feature]));
            }
        };
    const auto append_segments =
        [&sites](
            const CanonicalReachRing2& ring) {
            const std::int64_t ring_index =
                global_ring_index(ring);
            for (std::size_t feature = 0;
                 feature < ring.points.size();
                 ++feature) {
                const std::size_t target =
                    (feature + 1)
                    % ring.points.size();
                sites.push_back(
                    CanonicalMatSite2::
                        build_open_segment(
                            ring_index,
                            checked_index(
                                feature,
                                "MAT segment feature ordinal"),
                            checked_index(
                                target,
                                "MAT segment target feature ordinal"),
                            ring.points[feature],
                            ring.points[target]));
            }
        };
    append_points(input.outer);
    for (const CanonicalReachRing2& hole :
         input.holes) {
        append_points(hole);
    }
    append_segments(input.outer);
    for (const CanonicalReachRing2& hole :
         input.holes) {
        append_segments(hole);
    }
    require_canonical_catalog(sites);
    return CanonicalMatSiteCatalog2(
        std::move(sites));
}
