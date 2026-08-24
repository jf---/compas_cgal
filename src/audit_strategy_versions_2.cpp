#include "audit_strategy_versions_2.h"

#include "canonical_encoding.h"
#include "exact_depletion_2.h"

#include <string>

const std::string& audit_native_decision_contract_version()
{
    static const std::string version = "audit-native-decision-contract-v1";
    return version;
}

const std::string& audit_native_depletion_contract_version()
{
    static const std::string version = canonical_encode_tagged_union(
        "audit-native-depletion-contract-v1",
        canonical_encode_component_map({
            {"arc", exact_arc_depletion_strategy_version()},
            {"segment-circle", exact_depletion_strategy_version()},
        }));
    return version;
}
