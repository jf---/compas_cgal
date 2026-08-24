#include "audit_stock_identity_2.h"

#include "canonical_encoding.h"

#include <CGAL/Polygon_2.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

struct AuthoredVertex2 {
    double x;
    double y;
    EPoint exact;
    std::string canonical_bytes;
};

std::string canonical_authored_point(double x, double y)
{
    return canonical_encode_tagged_union(
        "point2-world-xy-v1",
        canonical_encode_sequence({
            canonical_encode_binary64(x),
            canonical_encode_binary64(y),
        }));
}

std::string canonical_ring(
    Eigen::Ref<const compas::RowMatrixXd> matrix,
    bool outer)
{
    if (matrix.cols() != 2) {
        throw AuditNativeStockShapeError(
            "native stock rings must have exactly two world-XY columns");
    }
    std::vector<AuthoredVertex2> vertices;
    vertices.reserve(static_cast<std::size_t>(matrix.rows()));
    for (Eigen::Index index = 0; index < matrix.rows(); ++index) {
        double x = matrix(index, 0);
        double y = matrix(index, 1);
        if (!std::isfinite(x) || !std::isfinite(y)) {
            throw AuditNativeStockNonFiniteInputError(
                "native stock ring coordinates must be finite binary64 values");
        }
        if (x == 0.0) {
            x = 0.0;
        }
        if (y == 0.0) {
            y = 0.0;
        }
        vertices.push_back({
            x,
            y,
            EPoint(Epeck::FT(x), Epeck::FT(y)),
            canonical_authored_point(x, y),
        });
    }
    if (vertices.size() >= 2
        && vertices.front().canonical_bytes == vertices.back().canonical_bytes) {
        vertices.pop_back();
    }
    if (vertices.size() < 3) {
        throw AuditNativeStockShapeError(
            "native stock rings require at least three open vertices");
    }
    std::unordered_set<std::string> unique;
    for (const AuthoredVertex2& vertex : vertices) {
        if (!unique.insert(vertex.canonical_bytes).second) {
            throw AuditNativeStockRingError(
                "native stock ring contains a repeated interior vertex");
        }
    }
    CGAL::Polygon_2<Epeck> polygon;
    for (const AuthoredVertex2& vertex : vertices) {
        polygon.push_back(vertex.exact);
    }
    if (!polygon.is_simple()) {
        throw AuditNativeStockRingError(
            "native stock ring must be exactly simple");
    }
    const CGAL::Orientation orientation = polygon.orientation();
    if (orientation == CGAL::COLLINEAR) {
        throw AuditNativeStockRingError(
            "native stock ring must have nonzero exact area");
    }
    const bool should_reverse = outer
        ? orientation != CGAL::COUNTERCLOCKWISE
        : orientation != CGAL::CLOCKWISE;
    if (should_reverse) {
        std::reverse(vertices.begin(), vertices.end());
    }
    const auto canonical_start = std::min_element(
        vertices.begin(),
        vertices.end(),
        [](const AuthoredVertex2& lhs, const AuthoredVertex2& rhs) {
            return lhs.canonical_bytes < rhs.canonical_bytes;
        });
    std::rotate(vertices.begin(), canonical_start, vertices.end());
    std::vector<std::string> point_bytes;
    point_bytes.reserve(vertices.size());
    for (const AuthoredVertex2& vertex : vertices) {
        point_bytes.push_back(vertex.canonical_bytes);
    }
    return canonical_encode_tagged_union(
        outer ? "outer-ring-ccw-v1" : "hole-ring-cw-v1",
        canonical_encode_sequence(point_bytes));
}

} // namespace

AuditNativeStockIdentity2 AuditNativeStockIdentity2::build(
    Eigen::Ref<const compas::RowMatrixXd> boundary,
    const std::vector<compas::RowMatrixXd>& holes)
{
    const std::string boundary_bytes = canonical_ring(boundary, true);
    std::vector<std::string> hole_bytes;
    hole_bytes.reserve(holes.size());
    for (const compas::RowMatrixXd& hole : holes) {
        hole_bytes.push_back(canonical_ring(hole, false));
    }
    std::sort(hole_bytes.begin(), hole_bytes.end());
    if (std::adjacent_find(hole_bytes.begin(), hole_bytes.end())
        != hole_bytes.end()) {
        throw AuditNativeStockDuplicateHoleError(
            "native stock identity requires unique canonical holes");
    }
    std::string canonical = canonical_encode_tagged_union(
        "audit-native-stock-v1",
        canonical_encode_component_map({
            {"boundary", boundary_bytes},
            {"holes", canonical_encode_sequence(hole_bytes)},
        }));
    AuditNativeStockDigest2 digest =
        AuditNativeStockDigestAuthority2::hash_canonical(canonical);
    return AuditNativeStockIdentity2(
        std::move(canonical), std::move(digest));
}

AuditNativeStockIdentity2::AuditNativeStockIdentity2(
    std::string canonical_bytes,
    AuditNativeStockDigest2 digest)
    : canonical_bytes_(std::move(canonical_bytes)),
      digest_(std::move(digest))
{
}

const std::string& AuditNativeStockIdentity2::canonical_bytes() const noexcept
{
    return canonical_bytes_;
}

const AuditNativeStockDigest2& AuditNativeStockIdentity2::digest() const noexcept
{
    return digest_;
}
