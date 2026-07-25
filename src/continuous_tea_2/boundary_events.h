#pragma once

#include "../stock_2.h"

#include <cstddef>
#include <string>
#include <vector>

struct BoundaryFeatureRecord2 {
    GpsXCurve curve;
    std::string support_kind;
    std::vector<std::string> support_coefficients;
    std::vector<std::string> primitive_coefficients;
    std::string support_id;
    std::string source_exact;
    std::string target_exact;
    std::string source_vertex_id;
    std::string target_vertex_id;
    std::string material_side;
    std::string trim_predicate;
    std::string feature_id;
    std::size_t overlap_multiplicity;
};

struct BoundaryEvent2 {
    std::string kind;
    std::string first_feature_id;
    std::string second_feature_id;
    std::string vertex_id;
    std::string exact_overlap_record;
    unsigned int multiplicity;
};

std::vector<BoundaryFeatureRecord2> extract_boundary_records(
    const Stock2& stock);

std::vector<BoundaryEvent2> classify_boundary_pair(
    const Stock2& stock,
    std::size_t first_index,
    std::size_t second_index);

std::vector<BoundaryEvent2> extract_boundary_events(
    const Stock2& stock);
