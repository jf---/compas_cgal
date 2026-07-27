#include "reachable_errors_2.h"
#include "segment_site_catalog.h"
#include "segment_site_graph_csr.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

static_assert(
    !std::is_aggregate_v<CanonicalMatSite2>);
static_assert(
    !std::is_default_constructible_v<
        CanonicalMatSiteCatalog2>);

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

compas::RowMatrixXd rectangle()
{
    return matrix(
        {
            {-4.0, -2.0},
            {4.0, -2.0},
            {4.0, 2.0},
            {-4.0, 2.0},
        });
}

compas::RowMatrixXd rotated_reversed_rectangle()
{
    return matrix(
        {
            {4.0, 2.0},
            {4.0, -2.0},
            {-4.0, -2.0},
            {-4.0, 2.0},
        });
}

bool exact_site_equal(
    const MatTraits::Site_2& left,
    const MatTraits::Site_2& right)
{
    if (left.is_point() != right.is_point()) {
        return false;
    }
    if (left.is_point()) {
        return left.point() == right.point();
    }
    return (left.source() == right.source()
            && left.target() == right.target())
        || (left.source() == right.target()
            && left.target() == right.source());
}

bool catalog_equal(
    const CanonicalMatSiteCatalog2& left,
    const CanonicalMatSiteCatalog2& right)
{
    if (left.sites().size() != right.sites().size()
        || left.site_provenance()
            != right.site_provenance()) {
        return false;
    }
    for (std::size_t index = 0;
         index < left.sites().size();
         ++index) {
        const CanonicalMatSite2& left_site =
            left.sites()[index];
        const CanonicalMatSite2& right_site =
            right.sites()[index];
        if (left_site.stable_id()
                != right_site.stable_id()
            || left_site.provenance()
                != right_site.provenance()
            || left_site.endpoint_site_ids()
                != right_site.endpoint_site_ids()
            || !exact_site_equal(
                left_site.site(),
                right_site.site())) {
            return false;
        }
    }
    return true;
}

bool rectangle_catalog_is_canonical()
{
    const CanonicalMatSiteCatalog2 catalog =
        canonical_mat_site_catalog(
            canonical_reach_input(
                rectangle(),
                {},
                1.0));
    const std::vector<std::array<std::int64_t, 3>>
        expected_provenance{
            {0, 0, 0},
            {0, 0, 1},
            {0, 0, 2},
            {0, 0, 3},
            {1, 0, 0},
            {1, 0, 1},
            {1, 0, 2},
            {1, 0, 3},
        };
    if (catalog.site_provenance()
            != expected_provenance
        || catalog.sites().size() != 8) {
        return false;
    }
    const std::array<MatTraits::Point_2, 4>
        expected_points{
            MatTraits::Point_2(-4, -2),
            MatTraits::Point_2(4, -2),
            MatTraits::Point_2(4, 2),
            MatTraits::Point_2(-4, 2),
        };
    const std::string identity_tag(
        "mat-site-v1\0",
        12);
    for (std::size_t index = 0;
         index < 4;
         ++index) {
        const CanonicalMatSite2& point =
            catalog.sites()[index];
        const CanonicalMatSite2& segment =
            catalog.sites()[index + 4];
        if (!point.site().is_point()
            || point.site().point()
                != expected_points[index]
            || point.endpoint_site_ids().has_value()
            || !point.stable_id().starts_with(
                identity_tag)
            || !segment.site().is_segment()
            || !segment.endpoint_site_ids().has_value()
            || !segment.stable_id().starts_with(
                identity_tag)
            || (*segment.endpoint_site_ids())[0]
                != catalog.sites()[index].stable_id()
            || (*segment.endpoint_site_ids())[1]
                != catalog.sites()[
                    (index + 1) % 4].stable_id()
            || !((segment.site().source()
                          == expected_points[index]
                      && segment.site().target()
                          == expected_points[
                              (index + 1) % 4])
                || (segment.site().target()
                          == expected_points[index]
                    && segment.site().source()
                        == expected_points[
                            (index + 1) % 4]))) {
            return false;
        }
    }
    return catalog.sites().front().stable_id()
            < catalog.sites().back().stable_id()
        && catalog.index_of(
               catalog.sites()[6].stable_id())
            == 6;
}

bool input_symmetries_preserve_catalog()
{
    const CanonicalMatSiteCatalog2 canonical =
        canonical_mat_site_catalog(
            canonical_reach_input(
                rectangle(),
                {},
                1.0));
    const CanonicalMatSiteCatalog2 transformed =
        canonical_mat_site_catalog(
            canonical_reach_input(
                rotated_reversed_rectangle(),
                {},
                1.0));
    return catalog_equal(canonical, transformed);
}

bool radius_is_not_site_identity()
{
    const CanonicalMatSiteCatalog2 small_radius =
        canonical_mat_site_catalog(
            canonical_reach_input(
                rectangle(),
                {},
                0.5));
    const CanonicalMatSiteCatalog2 large_radius =
        canonical_mat_site_catalog(
            canonical_reach_input(
                rectangle(),
                {},
                2.0));
    return catalog_equal(
        small_radius,
        large_radius);
}

bool malformed_canonical_input_is_rejected()
{
    CanonicalReachInput2 malformed =
        canonical_reach_input(
            rectangle(),
            {},
            1.0);
    malformed.outer.canonical_ordinal = 1;
    try {
        static_cast<void>(
            canonical_mat_site_catalog(
                malformed));
    } catch (
        const InvalidReachableDomainInputError&) {
        return true;
    }
    return false;
}

bool hole_order_is_canonical()
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
    const CanonicalMatSiteCatalog2 forward =
        canonical_mat_site_catalog(
            canonical_reach_input(
                outer,
                {first, second},
                1.0));
    const CanonicalMatSiteCatalog2 reordered =
        canonical_mat_site_catalog(
            canonical_reach_input(
                outer,
                {second, first},
                1.0));
    if (!catalog_equal(forward, reordered)
        || forward.sites().size() != 24) {
        return false;
    }
    const std::vector<std::array<std::int64_t, 3>>
        provenance = forward.site_provenance();
    return provenance.front()
            == std::array<std::int64_t, 3>{0, 0, 0}
        && provenance[4]
            == std::array<std::int64_t, 3>{0, 1, 0}
        && provenance[8]
            == std::array<std::int64_t, 3>{0, 2, 0}
        && provenance[12]
            == std::array<std::int64_t, 3>{1, 0, 0}
        && provenance.back()
            == std::array<std::int64_t, 3>{1, 2, 3};
}

bool numeric_csr_maps_canonical_identities()
{
    const CanonicalMatSiteCatalog2 catalog =
        canonical_mat_site_catalog(
            canonical_reach_input(
                rectangle(),
                {},
                1.0));
    const MatExactGraph2 graph{
        {
            {
                "node-a",
                {},
                {
                    catalog.sites()[0].stable_id(),
                    catalog.sites()[4].stable_id(),
                },
                {},
            },
            {
                "node-b",
                {},
                {
                    catalog.sites()[2].stable_id(),
                    catalog.sites()[6].stable_id(),
                    catalog.sites()[7].stable_id(),
                },
                {},
            },
        },
        {},
        0,
        8,
    };
    const MatNumericNodeSiteCsr2 csr =
        numeric_node_site_csr(graph, catalog);
    const MatNumericNodeSiteCsr2 empty =
        numeric_node_site_csr(
            {{}, {}, 0, 0},
            catalog);
    bool unknown_rejected = false;
    try {
        MatExactGraph2 unknown = graph;
        unknown.nodes[0].generator_site_ids[0] =
            "unknown-site";
        std::sort(
            unknown.nodes[0].generator_site_ids.begin(),
            unknown.nodes[0].generator_site_ids.end());
        static_cast<void>(
            numeric_node_site_csr(
                unknown,
                catalog));
    } catch (const UnknownMatSiteIdentityError&) {
        unknown_rejected = true;
    }
    return csr.node_site_offsets
            == std::vector<std::int64_t>{0, 2, 5}
        && csr.node_site_ids
            == std::vector<std::int64_t>{
                0,
                4,
                2,
                6,
                7,
            }
        && empty.node_site_offsets
            == std::vector<std::int64_t>{0}
        && empty.node_site_ids.empty()
        && unknown_rejected;
}

} // namespace

bool site_catalog_gate()
{
    return rectangle_catalog_is_canonical()
        && input_symmetries_preserve_catalog()
        && radius_is_not_site_identity()
        && malformed_canonical_input_is_rejected()
        && hole_order_is_canonical()
        && numeric_csr_maps_canonical_identities();
}
