#include "segment_site_voronoi.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

static_assert(
    !std::is_aggregate_v<
        CanonicalMatVoronoiSource2>);
static_assert(
    !std::is_default_constructible_v<
        CanonicalMatVoronoiSource2>);
static_assert(
    !std::is_copy_constructible_v<
        CanonicalMatVoronoiSource2>);
static_assert(
    !std::is_move_constructible_v<
        CanonicalMatVoronoiSource2>);

using RawPair =
    std::tuple<int, std::string, std::string>;

compas::RowMatrixXd rectangle(
    const bool transformed)
{
    const std::vector<std::array<double, 2>> points =
        transformed
        ? std::vector<std::array<double, 2>>{
              {4.0, 2.0},
              {4.0, -2.0},
              {-4.0, -2.0},
              {-4.0, 2.0},
          }
        : std::vector<std::array<double, 2>>{
              {-4.0, -2.0},
              {4.0, -2.0},
              {4.0, 2.0},
              {-4.0, 2.0},
          };
    compas::RowMatrixXd result(
        static_cast<Eigen::Index>(points.size()),
        2);
    for (std::size_t index = 0;
         index < points.size();
         ++index) {
        result(
            static_cast<Eigen::Index>(index),
            0) = points[index][0];
        result(
            static_cast<Eigen::Index>(index),
            1) = points[index][1];
    }
    return result;
}

CanonicalMatSiteCatalog2 catalog(
    const bool transformed)
{
    return canonical_mat_site_catalog(
        canonical_reach_input(
            rectangle(transformed),
            {},
            1.0));
}

std::vector<std::string> catalog_ids(
    const CanonicalMatSiteCatalog2& sites)
{
    std::vector<std::string> ids;
    ids.reserve(sites.sites().size());
    for (const CanonicalMatSite2& site :
         sites.sites()) {
        ids.push_back(site.stable_id());
    }
    return ids;
}

std::vector<std::string> live_ids(
    const CanonicalMatVoronoiSource2& source)
{
    std::vector<std::string> ids;
    const SegmentSiteDelaunay2& dual =
        source.voronoi().dual();
    for (auto vertex =
             dual.finite_vertices_begin();
         vertex != dual.finite_vertices_end();
         ++vertex) {
        ids.push_back(
            source.site_index().stable_id(
                vertex->site()));
    }
    std::sort(ids.begin(), ids.end());
    return ids;
}

std::vector<RawPair> raw_pair_signature(
    const CanonicalMatVoronoiSource2& source)
{
    std::vector<RawPair> signature;
    for (auto halfedge =
             source.voronoi().halfedges_begin();
         halfedge
         != source.voronoi().halfedges_end();
         ++halfedge) {
        const MatTraits::Site_2 up =
            halfedge->up()->site();
        const MatTraits::Site_2 down =
            halfedge->down()->site();
        std::string up_id =
            source.site_index().stable_id(up);
        std::string down_id =
            source.site_index().stable_id(
                down);
        if (down_id < up_id) {
            continue;
        }
        const int kind =
            up.is_segment()
                && down.is_segment()
            ? 2
            : up.is_point()
                    != down.is_point()
                ? 1
                : 0;
        signature.emplace_back(
            kind,
            std::move(up_id),
            std::move(down_id));
    }
    std::sort(
        signature.begin(),
        signature.end());
    return signature;
}

bool swap_transfer_is_complete()
{
    const CanonicalMatSiteCatalog2 sites =
        catalog(false);
    CanonicalMatDelaunaySource2 delaunay =
        CanonicalMatDelaunaySource2::build(
            sites);
    const CanonicalMatVoronoiSource2 voronoi =
        CanonicalMatVoronoiSource2::build(
            std::move(delaunay));
    const std::vector<RawPair> signature =
        raw_pair_signature(voronoi);
    const std::size_t point_segment_count =
        static_cast<std::size_t>(
            std::count_if(
                signature.begin(),
                signature.end(),
                [](const RawPair& pair) {
                    return std::get<0>(pair)
                        == 1;
                }));
    const std::size_t segment_segment_count =
        static_cast<std::size_t>(
            std::count_if(
                signature.begin(),
                signature.end(),
                [](const RawPair& pair) {
                    return std::get<0>(pair)
                        == 2;
                }));
    return std::distance(
               delaunay.delaunay()
                   .finite_vertices_begin(),
               delaunay.delaunay()
                   .finite_vertices_end())
            == 0
        && voronoi.voronoi().dual().is_valid()
        && live_ids(voronoi)
            == catalog_ids(sites)
        && signature.size() == 13
        && point_segment_count == 8
        && segment_segment_count == 5;
}

bool adaptor_is_input_symmetry_invariant()
{
    const CanonicalMatSiteCatalog2 first_sites =
        catalog(false);
    const CanonicalMatSiteCatalog2 second_sites =
        catalog(true);
    CanonicalMatDelaunaySource2 first_source =
        CanonicalMatDelaunaySource2::build(
            first_sites);
    CanonicalMatDelaunaySource2 second_source =
        CanonicalMatDelaunaySource2::build(
            second_sites);
    const CanonicalMatVoronoiSource2 first =
        CanonicalMatVoronoiSource2::build(
            std::move(first_source));
    const CanonicalMatVoronoiSource2 second =
        CanonicalMatVoronoiSource2::build(
            std::move(second_source));
    return raw_pair_signature(first)
        == raw_pair_signature(second);
}

bool second_consumption_is_rejected()
{
    const CanonicalMatSiteCatalog2 sites =
        catalog(false);
    CanonicalMatDelaunaySource2 delaunay =
        CanonicalMatDelaunaySource2::build(
            sites);
    const CanonicalMatVoronoiSource2 first =
        CanonicalMatVoronoiSource2::build(
            std::move(delaunay));
    static_cast<void>(first);
    try {
        static_cast<void>(
            CanonicalMatVoronoiSource2::build(
                std::move(delaunay)));
    } catch (
        const ConsumedCanonicalMatDelaunaySourceError&) {
        return true;
    }
    return false;
}

} // namespace

bool catalog_voronoi_gate()
{
    return swap_transfer_is_complete()
        && adaptor_is_input_symmetry_invariant()
        && second_consumption_is_rejected();
}
