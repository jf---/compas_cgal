#pragma once

#include "exact_build_audit_2.h"
#include "exact_region_2.h"

#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arrangement_2.h>
#include <Eigen/Core>

namespace compas {

using RowMatrixXd =
    Eigen::Matrix<
        double,
        Eigen::Dynamic,
        Eigen::Dynamic,
        Eigen::RowMajor>;

} // namespace compas

struct ReachCurveLabels2 {
    std::vector<std::string> source_piece_ids;
    std::vector<std::string> primitive_ids;
    bool operator==(const ReachCurveLabels2&) const = default;
};

struct MergeReachCurveLabels2 {
    ReachCurveLabels2 operator()(
        const ReachCurveLabels2& left,
        const ReachCurveLabels2& right) const;
};

using ReachDataTraits2 = CGAL::Arr_curve_data_traits_2<
    ReachTraits,
    ReachCurveLabels2,
    MergeReachCurveLabels2>;

struct ReachFaceState2 {
    bool classified = false;
    bool selected = false;
    bool outer_active = false;
    std::size_t active_holes = 0;
    std::size_t active_forbidden = 0;
};

using ReachFaceDcel2 =
    CGAL::Arr_face_extended_dcel<ReachDataTraits2, ReachFaceState2>;
using ReachArrangement2 =
    CGAL::Arrangement_2<ReachDataTraits2, ReachFaceDcel2>;

enum class ReachPrimitiveKind2 {
    Outer,
    Hole,
    Forbidden,
};

using ReachPrimitiveKinds2 =
    std::map<std::string, ReachPrimitiveKind2>;

struct CanonicalReachRing2 {
    bool outer;
    std::size_t canonical_ordinal;
    std::vector<ReachKernelPoint> points;
    std::vector<std::array<double, 2>> binary64_points;
    std::string record;
};

struct ReachableArrangement2;

struct CanonicalReachInput2 {
    CanonicalReachRing2 outer;
    std::vector<CanonicalReachRing2> holes;
    ReachFT radius;
    double binary64_radius;
    std::string recipe_record;

private:
    std::size_t input_vertex_count_ = 0;
    std::size_t ring_rotation_comparisons_ = 0;

    friend CanonicalReachInput2 canonical_reach_input(
        Eigen::Ref<const compas::RowMatrixXd> boundary,
        const std::vector<compas::RowMatrixXd>& holes,
        double tool_radius);
    friend ReachableArrangement2 build_reachable_arrangement(
        CanonicalReachInput2 input);
};

struct ReachableArrangement2 {
    ReachArrangement2 arrangement;
    CanonicalReachInput2 input;
    std::vector<std::string> source_records;
    ReachPolygonWithHoles design_polygon;
    ReachPolygonWithHoles center_polygon;
    ReachableDomainBuildAudit2 audit;
};

CanonicalReachInput2 canonical_reach_input(
    Eigen::Ref<const compas::RowMatrixXd> boundary,
    const std::vector<compas::RowMatrixXd>& holes,
    double tool_radius);

ReachableArrangement2 build_reachable_arrangement(
    CanonicalReachInput2 input);

void classify_faces_by_primitive_parity(
    ReachArrangement2& arrangement,
    const ReachPrimitiveKinds2& primitive_kinds);
