#include "segment_site_mat_proposal_table.h"

#include "segment_site_catalog.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace {

compas::RowMatrixXd matrix(const std::vector<std::array<double, 2>> &points) {
  compas::RowMatrixXd result(static_cast<Eigen::Index>(points.size()), 2);
  for (std::size_t index = 0; index < points.size(); ++index) {
    result(static_cast<Eigen::Index>(index), 0) = points[index][0];
    result(static_cast<Eigen::Index>(index), 1) = points[index][1];
  }
  return result;
}

CanonicalReachInput2 l_shape_input(const bool reversed) {
  return canonical_reach_input(reversed ? matrix({
                                              {0.0, 6.0},
                                              {2.0, 6.0},
                                              {2.0, 2.0},
                                              {6.0, 2.0},
                                              {6.0, 0.0},
                                              {0.0, 0.0},
                                          })
                                        : matrix({
                                              {0.0, 0.0},
                                              {6.0, 0.0},
                                              {6.0, 2.0},
                                              {2.0, 2.0},
                                              {2.0, 6.0},
                                              {0.0, 6.0},
                                          }),
                               {}, 1.0);
}

MatProposalSamplingGraph2 sampled_l_shape() {
  return canonical_l_shape_mat_proposal_graph(
      l_shape_input(false), CORE::BigRat(1), MatStationSpacingMm2::build(0.75),
      MatSagittaBoundMm2::build(0.02), 32);
}

std::vector<MatProposalSamplingRun2>
runs(const MatProposalSamplingGraph2 &sampled) {
  std::vector<MatProposalSamplingRun2> result;
  const auto &graph = sampled.profile_graph().graph();
  for (std::size_t edge_index = 0; edge_index < graph.edges.size();
       ++edge_index) {
    const std::size_t begin =
        static_cast<std::size_t>(sampled.sample_offsets()[edge_index]);
    const std::size_t end =
        static_cast<std::size_t>(sampled.sample_offsets()[edge_index + 1]);
    result.push_back(MatProposalSamplingRun2::build(
        graph.edges[edge_index].edge_id,
        std::vector<MatWorldXYProposalSample2>(
            sampled.samples().begin() + static_cast<std::ptrdiff_t>(begin),
            sampled.samples().begin() + static_cast<std::ptrdiff_t>(end))));
  }
  return result;
}

bool positive_radius_sample_fields_are_complete() {
  const MatProposalSamplingGraph2 sampled = sampled_l_shape();
  const MatToolRadiusMm2 radius = MatToolRadiusMm2::build(1.0);
  const MatNumericSampleTable2 table = numeric_sample_table(sampled, radius);
  const std::size_t count = sampled.samples().size();
  if (table.sample_centers.size() != count ||
      table.sample_clearance.size() != count ||
      table.sample_guide_radius.size() != count ||
      table.sample_flags.size() != count ||
      table.edge_sample_offsets != sampled.sample_offsets() ||
      table.sample_parameter.size() != count) {
    return false;
  }

  const std::array<std::int64_t, 2> verified{1, 1};
  for (std::size_t index = 0; index < count; ++index) {
    const auto &sample = sampled.samples()[index];
    if (table.sample_centers[index] !=
            std::array<double, 3>{sample.x_mm(), sample.y_mm(), 0.0} ||
        table.sample_parameter[index] != sample.parameter() ||
        table.sample_flags[index] != verified ||
        !std::isfinite(table.sample_clearance[index]) ||
        table.sample_clearance[index] < radius.value() ||
        table.sample_guide_radius[index] < 0.0 ||
        table.sample_guide_radius[index] !=
            (table.sample_clearance[index] - radius.value()) / 2.0) {
      return false;
    }
  }
  return true;
}

bool fixed_native_proposal_table_composes_graph_and_samples() {
  const CanonicalReachInput2 input = l_shape_input(false);
  const MatNumericProposalTable2 table =
      canonical_l_shape_mat_numeric_proposal_table(
          input, MatStationSpacingMm2::build(0.75),
          MatSagittaBoundMm2::build(0.02), 32);
  const MatProposalSamplingGraph2 sampled =
      canonical_l_shape_mat_proposal_graph(input, CORE::BigRat(1),
                                           MatStationSpacingMm2::build(0.75),
                                           MatSagittaBoundMm2::build(0.02), 32);
  return table.graph ==
             numeric_graph_table(sampled.profile_graph().graph(),
                                 canonical_mat_site_catalog(input)) &&
         table.samples ==
             numeric_sample_table(sampled, MatToolRadiusMm2::build(1.0));
}

bool fixed_native_proposal_table_is_reversal_invariant() {
  const MatNumericProposalTable2 first =
      canonical_l_shape_mat_numeric_proposal_table(
          l_shape_input(false), MatStationSpacingMm2::build(0.75),
          MatSagittaBoundMm2::build(0.02), 32);
  const MatNumericProposalTable2 reversed =
      canonical_l_shape_mat_numeric_proposal_table(
          l_shape_input(true), MatStationSpacingMm2::build(0.75),
          MatSagittaBoundMm2::build(0.02), 32);
  return first == reversed;
}

bool unverified_or_mismatched_sample_tables_fail_loudly() {
  const MatProposalSamplingGraph2 verified = sampled_l_shape();

  bool unverified_rejected = false;
  try {
    const MatProposalSamplingGraph2 unverified =
        MatProposalSamplingGraph2::build(verified.profile_graph(),
                                         runs(verified));
    static_cast<void>(
        numeric_sample_table(unverified, MatToolRadiusMm2::build(1.0)));
  } catch (const UnverifiedMatProposalSamplingGraphError &) {
    unverified_rejected = true;
  }

  bool radius_rejected = false;
  try {
    static_cast<void>(
        numeric_sample_table(verified, MatToolRadiusMm2::build(2.0)));
  } catch (const MismatchedMatProposalToolRadiusError &) {
    radius_rejected = true;
  }

  bool invalid_radius_rejected = false;
  try {
    static_cast<void>(
        MatToolRadiusMm2::build(std::numeric_limits<double>::infinity()));
  } catch (const InvalidMatProposalToolRadiusError &) {
    invalid_radius_rejected = true;
  }
  return unverified_rejected && radius_rejected && invalid_radius_rejected;
}

} // namespace

bool mat_proposal_table_gate() {
  return positive_radius_sample_fields_are_complete() &&
         fixed_native_proposal_table_composes_graph_and_samples() &&
         fixed_native_proposal_table_is_reversal_invariant() &&
         unverified_or_mismatched_sample_tables_fail_loudly();
}
