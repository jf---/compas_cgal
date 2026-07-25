#pragma once

#include "compas.h"
#include "exact_region_2.h"

#include <string>
#include <vector>

class PocketNotMachinableError : public ReachableDomainConstructionError {
public:
    using ReachableDomainConstructionError::ReachableDomainConstructionError;
};

struct ReachableDomainCertificate2 {
    std::string strategy_version;
    std::vector<std::string> source_curve_records;
    std::vector<std::string> selected_cell_records;
    std::vector<std::string> component_records;
    bool exact_disk_containment;
    bool exact_reconstruction;
    bool reachable_subset_of_design;
    bool residual_matches_difference;
    std::string input_recipe_record;

    bool matches_exact_inputs(
        Eigen::Ref<const compas::RowMatrixXd> boundary,
        const std::vector<compas::RowMatrixXd>& holes,
        double tool_radius) const;
};

class ReachableDomain2 {
public:
    ReachableDomain2(
        Eigen::Ref<const compas::RowMatrixXd> design_boundary,
        const std::vector<compas::RowMatrixXd>& holes,
        double tool_radius);

    ExactRegion2 design_region() const;
    ExactRegion2 center_domain() const;
    ExactRegion2 reachable_material() const;
    ExactRegion2 unreachable_residual() const;
    ReachableDomainCertificate2 certificate() const;

private:
    ExactRegion2 design_;
    ExactRegion2 center_domain_;
    ExactRegion2 reachable_material_;
    ExactRegion2 unreachable_residual_;
    ReachableDomainCertificate2 certificate_;
};
