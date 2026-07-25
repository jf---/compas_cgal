#pragma once

#include "exact_region_2.h"
#include "reachable_certificate_2.h"

#include <vector>

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
    const ReachableDomainBuildAudit2&
        build_audit_for_native_gate() const;

private:
    struct State {
        ExactRegion2 design;
        ExactRegion2 center_domain;
        ExactRegion2 reachable_material;
        ExactRegion2 unreachable_residual;
        ReachableDomainCertificate2 certificate;
        ReachableDomainBuildAudit2 audit;
    };

    explicit ReachableDomain2(State state);
    static State build_state(
        Eigen::Ref<const compas::RowMatrixXd> design_boundary,
        const std::vector<compas::RowMatrixXd>& holes,
        double tool_radius);

    State state_;
};
