#pragma once

#include "partition_certificate.h"

#include <string>
#include <vector>

class InvalidFullCirclePairRequestError
    : public EventSubstrateError {
public:
    using EventSubstrateError::EventSubstrateError;
};

class FullCirclePairRequest2 {
public:
    static FullCirclePairRequest2 build(
        std::string center_chart_id,
        std::string first_feature_id,
        std::string first_support_id,
        std::string first_support_kind,
        std::string first_rim_chart_id,
        std::string second_feature_id,
        std::string second_support_id,
        std::string second_support_kind,
        std::string second_rim_chart_id);

    static FullCirclePairRequest2 decode(
        const std::string& canonical_source);

    const std::string& center_chart_id() const;
    const std::string& first_feature_id() const;
    const std::string& first_support_id() const;
    const std::string& first_support_kind() const;
    const std::string& first_rim_chart_id() const;
    const std::string& second_feature_id() const;
    const std::string& second_support_id() const;
    const std::string& second_support_kind() const;
    const std::string& second_rim_chart_id() const;
    std::string canonical_source() const;

private:
    FullCirclePairRequest2(
        std::string center_chart_id,
        std::string first_feature_id,
        std::string first_support_id,
        std::string first_support_kind,
        std::string first_rim_chart_id,
        std::string second_feature_id,
        std::string second_support_id,
        std::string second_support_kind,
        std::string second_rim_chart_id);

    std::string center_chart_id_;
    std::string first_feature_id_;
    std::string first_support_id_;
    std::string first_support_kind_;
    std::string first_rim_chart_id_;
    std::string second_feature_id_;
    std::string second_support_id_;
    std::string second_support_kind_;
    std::string second_rim_chart_id_;
};

struct FullCirclePairProjectionBundle2 {
    std::vector<ProjectionInput2> projections;
    std::vector<ProjectionRecord2> constants;
};

std::vector<std::string>
full_circle_exhaustive_pair_requests(
    const std::vector<std::string>& line_sources,
    const std::vector<std::string>& circle_sources);

FullCirclePairProjectionBundle2
derive_full_circle_pair_cap_projections(
    const std::vector<std::string>& line_sources,
    const std::vector<std::string>& circle_sources,
    const std::string& cap_chord_ratio,
    const std::vector<ProjectionRecord2>& pullbacks,
    const std::vector<std::string>& pair_requests);
