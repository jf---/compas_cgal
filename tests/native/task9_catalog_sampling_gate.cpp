#include "segment_site_catalog_sampling.h"

#include "segment_site_catalog.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
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

MatProposalSamplingGraph2 sampled_l_shape(const bool reversed,
                                          const double spacing,
                                          const double sagitta) {
  return canonical_l_shape_mat_proposal_graph(
      l_shape_input(reversed), CORE::BigRat(0),
      MatStationSpacingMm2::build(spacing), MatSagittaBoundMm2::build(sagitta),
      32);
}

bool production_samples_bind_every_canonical_edge() {
  const MatProposalSamplingGraph2 sampled = sampled_l_shape(false, 0.75, 0.02);
  const MatExactGraph2 &graph = sampled.profile_graph().graph();
  const auto &offsets = sampled.sample_offsets();
  const auto &samples = sampled.samples();
  if (graph.edges.size() != 9 ||
      sampled.profile_graph().profiles().size() != graph.edges.size() ||
      offsets.size() != graph.edges.size() + 1 || offsets.front() != 0 ||
      offsets.back() != static_cast<std::int64_t>(samples.size())) {
    return false;
  }

  std::size_t lines = 0;
  std::size_t parabolas = 0;
  for (std::size_t edge_index = 0; edge_index < graph.edges.size();
       ++edge_index) {
    const MatExactGraphEdge2 &edge = graph.edges[edge_index];
    lines += edge.primitive_kind == "LINE";
    parabolas += edge.primitive_kind == "PARABOLA";
    const std::int64_t begin = offsets[edge_index];
    const std::int64_t end = offsets[edge_index + 1];
    if (end - begin < 2) {
      return false;
    }
    for (std::int64_t sample_index = begin; sample_index < end;
         ++sample_index) {
      const auto &sample = samples[static_cast<std::size_t>(sample_index)];
      if (sample.edge_id() != edge.edge_id ||
          sample.exact_parameter_id().find(edge.edge_id) == std::string::npos ||
          !std::isfinite(sample.parameter()) || !std::isfinite(sample.x_mm()) ||
          !std::isfinite(sample.y_mm()) ||
          (sample_index > begin &&
           samples[static_cast<std::size_t>(sample_index - 1)].parameter() >=
               sample.parameter())) {
        return false;
      }
    }
  }
  return lines == 7 && parabolas == 2;
}

bool production_sampling_is_deterministic_and_reversal_invariant() {
  const MatProposalSamplingGraph2 first = sampled_l_shape(false, 0.75, 0.02);
  const MatProposalSamplingGraph2 repeated = sampled_l_shape(false, 0.75, 0.02);
  const MatProposalSamplingGraph2 reversed = sampled_l_shape(true, 0.75, 0.02);
  return first.sample_offsets() == repeated.sample_offsets() &&
         first.samples() == repeated.samples() &&
         first.sample_offsets() == reversed.sample_offsets() &&
         first.samples() == reversed.samples();
}

bool policy_refinement_changes_only_proposals() {
  const MatProposalSamplingGraph2 coarse = sampled_l_shape(false, 10.0, 10.0);
  const MatProposalSamplingGraph2 refined = sampled_l_shape(false, 0.5, 0.01);
  const MatExactGraph2 &coarse_graph = coarse.profile_graph().graph();
  const MatExactGraph2 &refined_graph = refined.profile_graph().graph();
  if (refined.samples().size() <= coarse.samples().size() ||
      coarse_graph.edges.size() != refined_graph.edges.size()) {
    return false;
  }
  for (std::size_t index = 0; index < coarse_graph.edges.size(); ++index) {
    if (coarse_graph.edges[index].edge_id !=
            refined_graph.edges[index].edge_id ||
        coarse.profile_graph().profiles()[index].edge_id() !=
            refined.profile_graph().profiles()[index].edge_id()) {
      return false;
    }
  }
  return true;
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

bool malformed_sampling_graph_fails_loudly() {
  const MatProposalSamplingGraph2 valid = sampled_l_shape(false, 0.75, 0.02);

  bool missing_rejected = false;
  try {
    auto missing = runs(valid);
    missing.pop_back();
    static_cast<void>(MatProposalSamplingGraph2::build(valid.profile_graph(),
                                                       std::move(missing)));
  } catch (const IncompleteMatProposalSamplingGraphError &) {
    missing_rejected = true;
  }

  bool order_rejected = false;
  try {
    auto reordered = runs(valid);
    std::swap(reordered[0], reordered[1]);
    static_cast<void>(MatProposalSamplingGraph2::build(valid.profile_graph(),
                                                       std::move(reordered)));
  } catch (const IncompleteMatProposalSamplingGraphError &) {
    order_rejected = true;
  }

  bool owner_rejected = false;
  try {
    auto valid_runs = runs(valid);
    static_cast<void>(MatProposalSamplingRun2::build(valid_runs[1].edge_id(),
                                                     valid_runs[0].samples()));
  } catch (const InvalidMatProposalSamplingRunError &) {
    owner_rejected = true;
  }
  return missing_rejected && order_rejected && owner_rejected;
}

} // namespace

bool catalog_sampling_gate() {
  return production_samples_bind_every_canonical_edge() &&
         production_sampling_is_deterministic_and_reversal_invariant() &&
         policy_refinement_changes_only_proposals() &&
         malformed_sampling_graph_fails_loudly();
}
