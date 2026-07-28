#include "segment_site_mat_proposal_table.h"

#include "segment_site_catalog.h"
#include "segment_site_rational_sources.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace {

double evaluate_reported(const RationalPolynomial &coefficients,
                         const double parameter) {
  double result = 0.0;
  for (auto coefficient = coefficients.rbegin();
       coefficient != coefficients.rend(); ++coefficient) {
    result = result * parameter + CGAL::to_double(*coefficient);
  }
  return result;
}

} // namespace

MatToolRadiusMm2::MatToolRadiusMm2(const double value,
                                   CORE::BigRat squared_exact)
    : value_(value), squared_exact_(std::move(squared_exact)) {}

MatToolRadiusMm2 MatToolRadiusMm2::build(const double value) {
  if (!std::isfinite(value) || value <= 0.0) {
    throw InvalidMatProposalToolRadiusError(
        "MAT proposal tool radius must be finite and positive");
  }
  const CORE::BigRat exact = exact_mat_input_rational(MatTraits::FT(value));
  return MatToolRadiusMm2(value, exact * exact);
}

double MatToolRadiusMm2::value() const noexcept { return value_; }

const CORE::BigRat &MatToolRadiusMm2::squared_exact() const noexcept {
  return squared_exact_;
}

MatNumericSampleTable2
numeric_sample_table(const MatProposalSamplingGraph2 &sampled,
                     const MatToolRadiusMm2 &tool_radius) {
  if (!sampled.clearance_radius_squared().has_value() ||
      sampled.sample_exact_flags().size() != sampled.samples().size()) {
    throw UnverifiedMatProposalSamplingGraphError(
        "MAT proposal samples have no exact production verdict");
  }
  if (*sampled.clearance_radius_squared() != tool_radius.squared_exact()) {
    throw MismatchedMatProposalToolRadiusError(
        "MAT proposal tool radius does not match its exact clearance bound");
  }

  MatNumericSampleTable2 table;
  table.edge_sample_offsets = sampled.sample_offsets();
  table.sample_flags = sampled.sample_exact_flags();
  table.sample_centers.reserve(sampled.samples().size());
  table.sample_clearance.reserve(sampled.samples().size());
  table.sample_guide_radius.reserve(sampled.samples().size());
  table.sample_parameter.reserve(sampled.samples().size());

  const std::array<std::int64_t, 2> verified{1, 1};
  for (std::size_t edge_index = 0;
       edge_index < sampled.profile_graph().profiles().size(); ++edge_index) {
    const auto &profile = sampled.profile_graph().profiles()[edge_index];
    const std::size_t begin =
        static_cast<std::size_t>(sampled.sample_offsets()[edge_index]);
    const std::size_t end =
        static_cast<std::size_t>(sampled.sample_offsets()[edge_index + 1]);
    for (std::size_t sample_index = begin; sample_index < end; ++sample_index) {
      const auto &sample = sampled.samples()[sample_index];
      if (sampled.sample_exact_flags()[sample_index] != verified) {
        throw UnverifiedMatProposalSamplingGraphError(
            "MAT proposal sample exact verdict is not complete");
      }
      const double squared_clearance = evaluate_reported(
          profile.squared_clearance().coefficients(), sample.parameter());
      if (!std::isfinite(squared_clearance) || squared_clearance < 0.0) {
        throw InvalidMatProposalReportingDataError(
            "MAT proposal sample has invalid reported squared clearance");
      }
      const double clearance = std::sqrt(squared_clearance);
      const double guide_radius = (clearance - tool_radius.value()) / 2.0;
      if (!std::isfinite(clearance) || clearance < tool_radius.value() ||
          !std::isfinite(guide_radius) || guide_radius < 0.0) {
        throw InvalidMatProposalReportingDataError(
            "MAT proposal sample violates its reported clearance bound");
      }
      table.sample_centers.push_back({sample.x_mm(), sample.y_mm(), 0.0});
      table.sample_clearance.push_back(clearance);
      table.sample_guide_radius.push_back(guide_radius);
      table.sample_parameter.push_back(sample.parameter());
    }
  }
  return table;
}

MatNumericProposalTable2 canonical_l_shape_mat_numeric_proposal_table(
    const CanonicalReachInput2 &input,
    const MatStationSpacingMm2 &station_spacing,
    const MatSagittaBoundMm2 &max_sagitta,
    const std::size_t max_refinement_depth) {
  const MatToolRadiusMm2 tool_radius =
      MatToolRadiusMm2::build(input.binary64_radius);
  const MatProposalSamplingGraph2 sampled =
      canonical_l_shape_mat_proposal_graph(input, tool_radius.squared_exact(),
                                           station_spacing, max_sagitta,
                                           max_refinement_depth);
  return {
      numeric_graph_table(sampled.profile_graph().graph(),
                          canonical_mat_site_catalog(input)),
      numeric_sample_table(sampled, tool_radius),
  };
}
