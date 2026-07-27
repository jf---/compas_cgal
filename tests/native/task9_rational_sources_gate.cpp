#include "segment_site_rational_sources.h"

#include <array>
#include <cstddef>
#include <string>
#include <type_traits>
#include <vector>

namespace {

static_assert(
    !std::is_aggregate_v<
        CanonicalMatRationalPointSource2>);
static_assert(
    !std::is_aggregate_v<
        CanonicalMatRationalOpenSegmentSource2>);
static_assert(
    !std::is_default_constructible_v<
        CanonicalMatRationalSources2>);

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

bool rational_sources_equal(
    const CanonicalMatRationalSources2& left,
    const CanonicalMatRationalSources2& right)
{
    if (left.points().size()
            != right.points().size()
        || left.segments().size()
            != right.segments().size()) {
        return false;
    }
    for (std::size_t index = 0;
         index < left.points().size();
         ++index) {
        const auto& first = left.points()[index];
        const auto& second = right.points()[index];
        if (first.stable_site_id
                != second.stable_site_id
            || first.x != second.x
            || first.y != second.y) {
            return false;
        }
    }
    for (std::size_t index = 0;
         index < left.segments().size();
         ++index) {
        const auto& first =
            left.segments()[index];
        const auto& second =
            right.segments()[index];
        if (first.stable_site_id
                != second.stable_site_id
            || first.source_point_id
                != second.source_point_id
            || first.target_point_id
                != second.target_point_id
            || first.support.line_a
                != second.support.line_a
            || first.support.line_b
                != second.support.line_b
            || first.support.line_c
                != second.support.line_c) {
            return false;
        }
    }
    return true;
}

bool rectangle_sources_are_exact()
{
    const CanonicalMatSiteCatalog2 catalog =
        rectangle_catalog(false);
    const CanonicalMatRationalSources2 sources =
        CanonicalMatRationalSources2::build(
            catalog);
    const std::array<std::array<CORE::BigRat, 2>, 4>
        expected_points{
            std::array<CORE::BigRat, 2>{-4, -2},
            std::array<CORE::BigRat, 2>{4, -2},
            std::array<CORE::BigRat, 2>{4, 2},
            std::array<CORE::BigRat, 2>{-4, 2},
        };
    const std::array<std::array<CORE::BigRat, 3>, 4>
        expected_lines{
            std::array<CORE::BigRat, 3>{0, 1, 2},
            std::array<CORE::BigRat, 3>{1, 0, -4},
            std::array<CORE::BigRat, 3>{0, 1, -2},
            std::array<CORE::BigRat, 3>{1, 0, 4},
        };
    if (sources.points().size() != 4
        || sources.segments().size() != 4) {
        return false;
    }
    for (std::size_t index = 0;
         index < 4;
         ++index) {
        const auto& point =
            sources.points()[index];
        const auto& segment =
            sources.segments()[index];
        if (point.stable_site_id
                != catalog.sites()[index]
                       .stable_id()
            || point.x
                != expected_points[index][0]
            || point.y
                != expected_points[index][1]
            || segment.stable_site_id
                != catalog.sites()[index + 4]
                       .stable_id()
            || segment.source_point_id
                != catalog.sites()[index]
                       .stable_id()
            || segment.target_point_id
                != catalog.sites()[
                    (index + 1) % 4]
                       .stable_id()
            || segment.support.stable_site_id
                != segment.stable_site_id
            || segment.support.line_a
                != expected_lines[index][0]
            || segment.support.line_b
                != expected_lines[index][1]
            || segment.support.line_c
                != expected_lines[index][2]) {
            return false;
        }
    }
    return true;
}

bool binary64_coordinates_are_not_decimalized()
{
    const CanonicalMatSiteCatalog2 catalog =
        canonical_mat_site_catalog(
            canonical_reach_input(
                matrix(
                    {
                        {0.1, 0.2},
                        {1.0, 0.2},
                        {1.0, 1.0},
                        {0.1, 1.0},
                    }),
                {},
                0.25));
    const CanonicalMatRationalSources2 sources =
        CanonicalMatRationalSources2::build(
            catalog);
    const CORE::BigInt numerator(
        "3602879701896397");
    const CORE::BigRat exact_point_one(
        numerator,
        CORE::BigInt(1) << 55);
    const CORE::BigRat exact_point_two(
        numerator,
        CORE::BigInt(1) << 54);
    return sources.points().front().x
            == exact_point_one
        && sources.points().front().y
            == exact_point_two;
}

bool radical_coordinate_is_rejected()
{
    try {
        static_cast<void>(
            exact_mat_input_rational(
                CORE::sqrt(
                    CORE::Expr(2))));
    } catch (
        const NonRationalCanonicalMatCoordinateError&) {
        return true;
    }
    return false;
}

bool source_projection_is_input_symmetry_invariant()
{
    const CanonicalMatRationalSources2 first =
        CanonicalMatRationalSources2::build(
            rectangle_catalog(false));
    const CanonicalMatRationalSources2 second =
        CanonicalMatRationalSources2::build(
            rectangle_catalog(true));
    return rational_sources_equal(
        first,
        second);
}

} // namespace

bool rational_sources_gate()
{
    return rectangle_sources_are_exact()
        && binary64_coordinates_are_not_decimalized()
        && radical_coordinate_is_rejected()
        && source_projection_is_input_symmetry_invariant();
}
