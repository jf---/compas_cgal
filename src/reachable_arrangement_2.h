#pragma once

#include "exact_build_audit_2.h"
#include "exact_region_2.h"
#include "reachable_input_2.h"

#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arrangement_2.h>

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

struct ReachableArrangement2 {
    ReachArrangement2 arrangement;
    CanonicalReachInput2 input;
    std::vector<std::string> source_records;
    ReachPolygonWithHoles design_polygon;
    ReachPolygonWithHoles center_polygon;
    ReachableDomainBuildAudit2 audit;
};

ReachableArrangement2 build_reachable_arrangement(
    CanonicalReachInput2 input);

void classify_faces_by_primitive_parity(
    ReachArrangement2& arrangement,
    const ReachPrimitiveKinds2& primitive_kinds);
