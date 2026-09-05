#pragma once

#include "exact_build_audit_2.h"
#include "exact_region_2.h"
#include "reachable_input_2.h"

#include <memory>
#include <vector>

bool reach_curve_within_radius(
    const ReachKernelPoint& query,
    const ReachXCurve& curve,
    const ReachFT& radius);

class ReachableMaterialPredicate2 {
public:
    static ReachableMaterialPredicate2 build(
        Eigen::Ref<const compas::RowMatrixXd> design_boundary,
        const std::vector<compas::RowMatrixXd>& holes,
        double tool_radius);

    bool contains(double x, double y) const;
    const ReachableDomainBuildAudit2&
        build_audit_for_native_gate() const;

private:
    explicit ReachableMaterialPredicate2(
        std::shared_ptr<const struct ReachableMaterialPredicateStorage2>
            storage);

    std::shared_ptr<const struct ReachableMaterialPredicateStorage2>
        storage_;
};
