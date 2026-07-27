#include "segment_site_delaunay.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <string>
#include <type_traits>
#include <vector>

namespace {

static_assert(
    !std::is_aggregate_v<
        CanonicalMatDelaunaySource2>);
static_assert(
    !std::is_default_constructible_v<
        CanonicalMatDelaunaySource2>);
static_assert(
    !std::is_copy_constructible_v<
        CanonicalMatDelaunaySource2>);
static_assert(
    !std::is_move_constructible_v<
        CanonicalMatDelaunaySource2>);

compas::RowMatrixXd matrix(
    const std::vector<std::array<double, 2>>& points)
{
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

CanonicalMatSiteCatalog2 rectangle_catalog(
    const bool transformed)
{
    const compas::RowMatrixXd boundary =
        transformed
        ? matrix(
              {
                  {4.0, 2.0},
                  {4.0, -2.0},
                  {-4.0, -2.0},
                  {-4.0, 2.0},
              })
        : matrix(
              {
                  {-4.0, -2.0},
                  {4.0, -2.0},
                  {4.0, 2.0},
                  {-4.0, 2.0},
              });
    return canonical_mat_site_catalog(
        canonical_reach_input(
            boundary,
            {},
            1.0));
}

std::vector<std::string> catalog_ids(
    const CanonicalMatSiteCatalog2& catalog)
{
    std::vector<std::string> ids;
    ids.reserve(catalog.sites().size());
    for (const CanonicalMatSite2& site :
         catalog.sites()) {
        ids.push_back(site.stable_id());
    }
    return ids;
}

std::vector<std::string> live_ids(
    const CanonicalMatDelaunaySource2& source)
{
    std::vector<std::string> ids;
    for (auto vertex =
             source.delaunay()
                 .finite_vertices_begin();
         vertex
         != source.delaunay()
                .finite_vertices_end();
         ++vertex) {
        ids.push_back(
            source.site_index().stable_id(
                vertex->site()));
    }
    std::sort(ids.begin(), ids.end());
    return ids;
}

bool rectangle_source_is_bijective()
{
    const CanonicalMatSiteCatalog2 catalog =
        rectangle_catalog(false);
    const CanonicalMatDelaunaySource2 source =
        CanonicalMatDelaunaySource2::build(
            catalog);
    const CanonicalMatSite2& segment =
        catalog.sites()[4];
    const MatTraits::Site_2 reversed_segment =
        MatTraits::Site_2::construct_site_2(
            segment.site().target(),
            segment.site().source());
    return source.delaunay().is_valid()
        && source.site_index().size() == 8
        && static_cast<std::size_t>(
               std::distance(
                   source.delaunay()
                       .finite_vertices_begin(),
                   source.delaunay()
                       .finite_vertices_end()))
            == catalog.sites().size()
        && live_ids(source)
            == catalog_ids(catalog)
        && source.site_index().stable_id(
               reversed_segment)
            == segment.stable_id();
}

bool source_is_input_symmetry_invariant()
{
    const CanonicalMatSiteCatalog2 canonical =
        rectangle_catalog(false);
    const CanonicalMatSiteCatalog2 transformed =
        rectangle_catalog(true);
    const CanonicalMatDelaunaySource2
        canonical_source =
            CanonicalMatDelaunaySource2::build(
                canonical);
    const CanonicalMatDelaunaySource2
        transformed_source =
            CanonicalMatDelaunaySource2::build(
                transformed);
    return live_ids(canonical_source)
            == live_ids(transformed_source)
        && live_ids(canonical_source)
            == catalog_ids(canonical);
}

bool holes_are_inserted_once()
{
    const compas::RowMatrixXd outer =
        matrix(
            {
                {0.0, 0.0},
                {20.0, 0.0},
                {20.0, 20.0},
                {0.0, 20.0},
            });
    const compas::RowMatrixXd first =
        matrix(
            {
                {2.0, 2.0},
                {2.0, 4.0},
                {4.0, 4.0},
                {4.0, 2.0},
            });
    const compas::RowMatrixXd second =
        matrix(
            {
                {10.0, 10.0},
                {12.0, 10.0},
                {12.0, 12.0},
                {10.0, 12.0},
            });
    const CanonicalMatSiteCatalog2 catalog =
        canonical_mat_site_catalog(
            canonical_reach_input(
                outer,
                {second, first},
                1.0));
    const CanonicalMatDelaunaySource2 source =
        CanonicalMatDelaunaySource2::build(
            catalog);
    return source.site_index().size() == 24
        && live_ids(source)
            == catalog_ids(catalog);
}

bool unknown_geometry_is_rejected()
{
    const CanonicalMatSiteCatalog2 catalog =
        rectangle_catalog(false);
    const CanonicalMatDelaunaySource2 source =
        CanonicalMatDelaunaySource2::build(
            catalog);
    bool point_rejected = false;
    try {
        static_cast<void>(
            source.site_index().stable_id(
                MatTraits::Site_2::
                    construct_site_2(
                        MatTraits::Point_2(
                            99,
                            99))));
    } catch (
        const UnknownCanonicalMatSiteGeometryError&) {
        point_rejected = true;
    }
    bool segment_rejected = false;
    try {
        static_cast<void>(
            source.site_index().stable_id(
                MatTraits::Site_2::
                    construct_site_2(
                        MatTraits::Point_2(
                            99,
                            99),
                        MatTraits::Point_2(
                            100,
                            100))));
    } catch (
        const UnknownCanonicalMatSiteGeometryError&) {
        segment_rejected = true;
    }
    return point_rejected
        && segment_rejected;
}

bool duplicate_geometry_is_rejected()
{
    const compas::RowMatrixXd outer =
        matrix(
            {
                {0.0, 0.0},
                {10.0, 0.0},
                {10.0, 10.0},
                {0.0, 10.0},
            });
    const compas::RowMatrixXd touching_hole =
        matrix(
            {
                {0.0, 0.0},
                {1.0, 2.0},
                {2.0, 0.0},
            });
    const CanonicalMatSiteCatalog2 catalog =
        canonical_mat_site_catalog(
            canonical_reach_input(
                outer,
                {touching_hole},
                1.0));
    try {
        static_cast<void>(
            CanonicalMatDelaunaySource2::build(
                catalog));
    } catch (
        const DuplicateCanonicalMatSiteGeometryError&) {
        return true;
    }
    return false;
}

} // namespace

bool catalog_delaunay_gate()
{
    return rectangle_source_is_bijective()
        && source_is_input_symmetry_invariant()
        && holes_are_inserted_once()
        && unknown_geometry_is_rejected()
        && duplicate_geometry_is_rejected();
}
