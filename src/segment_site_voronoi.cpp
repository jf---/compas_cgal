#include "segment_site_voronoi.h"

#include <cstddef>
#include <iterator>
#include <set>
#include <string>
#include <utility>

CanonicalMatVoronoiSource2
CanonicalMatVoronoiSource2::build(
    CanonicalMatDelaunaySource2&& source)
{
    if (source.delaunay_
            .finite_vertices_begin()
        == source.delaunay_
               .finite_vertices_end()) {
        throw ConsumedCanonicalMatDelaunaySourceError(
            "canonical MAT segment-Delaunay source was already consumed");
    }
    return CanonicalMatVoronoiSource2(
        std::move(source));
}

CanonicalMatVoronoiSource2::
    CanonicalMatVoronoiSource2(
        CanonicalMatDelaunaySource2&& source)
    : voronoi_(source.delaunay_, true),
      site_index_(
          std::move(source.site_index_))
{
    if (source.delaunay_
            .finite_vertices_begin()
        != source.delaunay_
               .finite_vertices_end()) {
        throw InvalidCanonicalMatVoronoiTransferError(
            "segment-Delaunay source retained vertices after swap transfer");
    }
    if (!voronoi_.is_valid()) {
        throw InvalidCanonicalMatVoronoiTransferError(
            "catalog-fed Voronoi adaptor is invalid after swap transfer");
    }

    std::set<std::string> matched_ids;
    const SegmentSiteDelaunay2& dual =
        voronoi_.dual();
    for (auto vertex =
             dual.finite_vertices_begin();
         vertex != dual.finite_vertices_end();
         ++vertex) {
        if (!matched_ids.insert(
                site_index_.stable_id(
                    vertex->site()))
                 .second) {
            throw InvalidCanonicalMatVoronoiTransferError(
                "Voronoi adaptor duplicates a canonical MAT generator");
        }
    }
    if (matched_ids.size()
            != site_index_.size()) {
        throw InvalidCanonicalMatVoronoiTransferError(
            "Voronoi adaptor lost canonical MAT generators during transfer");
    }
}

const SegmentSiteVoronoi2&
CanonicalMatVoronoiSource2::voronoi() const noexcept
{
    return voronoi_;
}

const CanonicalMatSiteGeometryIndex2&
CanonicalMatVoronoiSource2::site_index() const noexcept
{
    return site_index_;
}
