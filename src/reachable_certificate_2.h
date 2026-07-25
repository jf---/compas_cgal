#pragma once

#include "reachable_arrangement_2.h"

#include <string>
#include <vector>

struct ReachableDomainCertificate2 {
    std::string strategy_version;
    std::vector<std::string> source_curve_records;
    std::vector<std::string> selected_cell_records;
    std::vector<std::string> component_records;
    bool exact_cell_selection = false;
    bool complete_source_provenance = false;
    bool reachable_subset_of_design = false;
    std::string input_recipe_record;

    bool matches_exact_inputs(
        Eigen::Ref<const compas::RowMatrixXd> boundary,
        const std::vector<compas::RowMatrixXd>& holes,
        double tool_radius) const;
};

ReachableDomainCertificate2 build_reachable_certificate(
    ReachableArrangement2& reachable,
    bool reachable_subset_of_design);
