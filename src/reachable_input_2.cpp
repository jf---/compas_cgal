#include "reachable_input_2.h"

#include "canonical_rotation_2.h"
#include "reachable_errors_2.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include <CGAL/Polygon_2.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>

namespace {

// Booth's algorithm performs at most three comparisons per ring element.
constexpr std::size_t MINIMAL_ROTATION_COMPARISON_FACTOR = 3;

bool point_equal(
    const ReachKernelPoint& left,
    const ReachKernelPoint& right)
{
    return CGAL::compare_xy(left, right) == CGAL::EQUAL;
}

std::string ring_record(
    bool outer,
    const std::vector<std::array<double, 2>>& points)
{
    std::vector<std::string> fields;
    fields.reserve(points.size() + 1);
    fields.emplace_back(outer ? "outer" : "hole");
    for (const std::array<double, 2>& point : points) {
        fields.push_back(reach_tagged_record(
            "binary64-point-v1",
            {
                reach_binary64_record(point[0]),
                reach_binary64_record(point[1]),
            }));
    }
    return reach_tagged_record("input-ring-v1", fields);
}

std::string reach_input_recipe_record(
    const std::string& outer_record,
    const std::vector<std::string>& hole_records,
    double binary64_radius)
{
    std::vector<std::string> fields{
        outer_record,
        reach_binary64_record(binary64_radius),
    };
    fields.insert(
        fields.end(),
        hole_records.begin(),
        hole_records.end());
    return reach_tagged_record("reachable-input-recipe-v1", fields);
}

std::string reach_input_integrity_record(
    const std::string& recipe_record,
    std::size_t input_vertex_count,
    std::size_t ring_rotation_comparisons)
{
    return reach_tagged_record(
        "canonical-reach-integrity-v1",
        {
            recipe_record,
            reach_u64_record(input_vertex_count),
            reach_u64_record(ring_rotation_comparisons),
        });
}

void rotate_ring(
    CanonicalReachRing2& ring,
    std::size_t& comparison_count)
{
    const std::size_t first = minimal_rotation_index(
        ring.points,
        [](const ReachKernelPoint& left, const ReachKernelPoint& right) {
            return CGAL::compare_xy(left, right);
        },
        comparison_count);
    std::rotate(
        ring.points.begin(),
        ring.points.begin() + static_cast<std::ptrdiff_t>(first),
        ring.points.end());
    std::rotate(
        ring.binary64_points.begin(),
        ring.binary64_points.begin() + static_cast<std::ptrdiff_t>(first),
        ring.binary64_points.end());
}

CanonicalReachRing2 canonical_ring(
    Eigen::Ref<const compas::RowMatrixXd> vertices,
    bool outer,
    std::size_t& comparison_count)
{
    if (vertices.cols() < 2 || vertices.rows() < 3) {
        throw InvalidReachableDomainInputError(
            "reachable-domain rings require at least three XY vertices");
    }
    CanonicalReachRing2 ring{outer, 0, {}, {}, {}};
    ring.points.reserve(static_cast<std::size_t>(vertices.rows()));
    ring.binary64_points.reserve(
        static_cast<std::size_t>(vertices.rows()));
    for (Eigen::Index row = 0; row < vertices.rows(); ++row) {
        const double x = vertices(row, 0);
        const double y = vertices(row, 1);
        if (!std::isfinite(x) || !std::isfinite(y)) {
            throw InvalidReachableDomainInputError(
                "reachable-domain vertices must be finite binary64 values");
        }
        const double canonical_x = x == 0.0 ? 0.0 : x;
        const double canonical_y = y == 0.0 ? 0.0 : y;
        ring.points.emplace_back(canonical_x, canonical_y);
        ring.binary64_points.push_back({canonical_x, canonical_y});
    }
    if (ring.points.size() > 3
        && point_equal(ring.points.front(), ring.points.back())) {
        ring.points.pop_back();
        ring.binary64_points.pop_back();
    }
    if (ring.points.size() < 3) {
        throw InvalidReachableDomainInputError(
            "reachable-domain rings require three distinct vertices");
    }
    for (std::size_t index = 0; index < ring.points.size(); ++index) {
        if (point_equal(
                ring.points[index],
                ring.points[(index + 1) % ring.points.size()])) {
            throw InvalidReachableDomainInputError(
                "reachable-domain edges must be nondegenerate");
        }
    }

    CGAL::Polygon_2<ReachKernel> polygon(
        ring.points.begin(),
        ring.points.end());
    if (!polygon.is_simple()) {
        throw InvalidReachableDomainInputError(
            "reachable-domain rings must be exact simple polygons");
    }
    const CGAL::Orientation orientation = polygon.orientation();
    if (orientation == CGAL::COLLINEAR) {
        throw InvalidReachableDomainInputError(
            "reachable-domain rings must have nonzero exact area");
    }
    if ((orientation == CGAL::COUNTERCLOCKWISE) != outer) {
        std::reverse(ring.points.begin(), ring.points.end());
        std::reverse(
            ring.binary64_points.begin(),
            ring.binary64_points.end());
    }
    rotate_ring(ring, comparison_count);
    ring.record = ring_record(outer, ring.binary64_points);
    return ring;
}

std::string validate_canonical_ring(
    const CanonicalReachRing2& ring,
    bool expected_outer,
    std::size_t expected_ordinal)
{
    if (ring.outer != expected_outer
        || ring.canonical_ordinal != expected_ordinal) {
        throw InvalidReachableDomainInputError(
            "canonical reach ring role or ordinal is invalid");
    }
    if (ring.points.size() < 3
        || ring.points.size() != ring.binary64_points.size()) {
        throw InvalidReachableDomainInputError(
            "canonical reach ring point mirrors are incomplete");
    }
    for (std::size_t index = 0; index < ring.points.size(); ++index) {
        const double x = ring.binary64_points[index][0];
        const double y = ring.binary64_points[index][1];
        if (!std::isfinite(x) || !std::isfinite(y)
            || CGAL::compare_xy(
                   ring.points[index],
                   ReachKernelPoint(x, y))
                != CGAL::EQUAL) {
            throw InvalidReachableDomainInputError(
                "canonical reach exact point diverges from its binary64 mirror");
        }
        if (point_equal(
                ring.points[index],
                ring.points[(index + 1) % ring.points.size()])) {
            throw InvalidReachableDomainInputError(
                "canonical reach ring contains a degenerate edge");
        }
    }
    CGAL::Polygon_2<ReachKernel> polygon(
        ring.points.begin(),
        ring.points.end());
    if (!polygon.is_simple()) {
        throw InvalidReachableDomainInputError(
            "canonical reach ring is not exact simple geometry");
    }
    const CGAL::Orientation expected_orientation =
        expected_outer
        ? CGAL::COUNTERCLOCKWISE
        : CGAL::CLOCKWISE;
    if (polygon.orientation() != expected_orientation) {
        throw InvalidReachableDomainInputError(
            "canonical reach ring orientation contradicts its role");
    }
    std::size_t validation_comparisons = 0;
    if (minimal_rotation_index(
            ring.points,
            [](const ReachKernelPoint& left,
               const ReachKernelPoint& right) {
                return CGAL::compare_xy(left, right);
            },
            validation_comparisons)
        != 0) {
        throw InvalidReachableDomainInputError(
            "canonical reach ring is not minimally rotated");
    }
    const std::string recomputed_record =
        ring_record(expected_outer, ring.binary64_points);
    if (ring.record != recomputed_record) {
        throw InvalidReachableDomainInputError(
            "canonical reach ring record does not match structural input");
    }
    return recomputed_record;
}

} // namespace

CanonicalReachInput2 canonical_reach_input(
    Eigen::Ref<const compas::RowMatrixXd> boundary,
    const std::vector<compas::RowMatrixXd>& holes,
    double tool_radius)
{
    if (!std::isfinite(tool_radius) || !(tool_radius > 0.0)) {
        throw InvalidReachableDomainInputError(
            "reachable-domain tool radius must be finite and positive");
    }
    std::size_t comparison_count = 0;
    CanonicalReachInput2 input;
    input.outer = canonical_ring(boundary, true, comparison_count);
    input.radius = ReachFT(tool_radius);
    input.binary64_radius = tool_radius;
    input.holes.reserve(holes.size());
    for (const compas::RowMatrixXd& hole : holes) {
        input.holes.push_back(
            canonical_ring(hole, false, comparison_count));
    }
    std::sort(
        input.holes.begin(),
        input.holes.end(),
        [](const CanonicalReachRing2& left,
           const CanonicalReachRing2& right) {
            return left.record < right.record;
        });
    for (std::size_t index = 0; index < input.holes.size(); ++index) {
        if (index != 0
            && input.holes[index - 1].record
                == input.holes[index].record) {
            throw InvalidReachableDomainInputError(
                "reachable-domain holes must have distinct exact identities");
        }
        input.holes[index].canonical_ordinal = index;
    }
    input.input_vertex_count_ = input.outer.points.size();
    for (const CanonicalReachRing2& hole : input.holes) {
        input.input_vertex_count_ += hole.points.size();
    }
    input.ring_rotation_comparisons_ = comparison_count;
    std::vector<std::string> hole_records;
    hole_records.reserve(input.holes.size());
    for (const CanonicalReachRing2& hole : input.holes) {
        hole_records.push_back(hole.record);
    }
    input.recipe_record =
        reach_input_recipe_record(
            input.outer.record,
            hole_records,
            input.binary64_radius);
    input.integrity_record_ = reach_input_integrity_record(
        input.recipe_record,
        input.input_vertex_count_,
        input.ring_rotation_comparisons_);
    return input;
}

void validate_canonical_reach_input(
    const CanonicalReachInput2& input)
{
    if (!std::isfinite(input.binary64_radius)
        || !(input.binary64_radius > 0.0)
        || CGAL::compare(input.radius, ReachFT(0)) != CGAL::LARGER
        || CGAL::compare(
               input.radius,
               ReachFT(input.binary64_radius))
            != CGAL::EQUAL) {
        throw InvalidReachableDomainInputError(
            "canonical reach radius must be a positive exact binary64 injection");
    }
    std::size_t vertex_count = input.outer.points.size();
    const std::string outer_record =
        validate_canonical_ring(input.outer, true, 0);
    std::vector<std::string> hole_records;
    hole_records.reserve(input.holes.size());
    for (std::size_t index = 0; index < input.holes.size(); ++index) {
        const CanonicalReachRing2& hole = input.holes[index];
        if (vertex_count
            > std::numeric_limits<std::size_t>::max()
                - hole.points.size()) {
            throw InvalidReachableDomainInputError(
                "canonical reach input vertex count overflows");
        }
        vertex_count += hole.points.size();
        std::string record =
            validate_canonical_ring(hole, false, index);
        if (!hole_records.empty()
            && !(hole_records.back() < record)) {
            throw InvalidReachableDomainInputError(
                "canonical reach holes are not strictly identity-sorted");
        }
        hole_records.push_back(std::move(record));
    }
    const std::string recomputed_recipe =
        reach_input_recipe_record(
            outer_record,
            hole_records,
            input.binary64_radius);
    if (input.recipe_record != recomputed_recipe
        || input.input_vertex_count_ != vertex_count
        || (vertex_count
                <= std::numeric_limits<std::size_t>::max()
                    / MINIMAL_ROTATION_COMPARISON_FACTOR
            && input.ring_rotation_comparisons_
                > MINIMAL_ROTATION_COMPARISON_FACTOR * vertex_count)
        || input.integrity_record_
            != reach_input_integrity_record(
                recomputed_recipe,
                input.input_vertex_count_,
                input.ring_rotation_comparisons_)) {
        throw InvalidReachableDomainInputError(
            "canonical reach input integrity does not match factory construction");
    }
}
