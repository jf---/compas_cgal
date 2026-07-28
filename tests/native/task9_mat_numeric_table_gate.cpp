#include "segment_site_mat_numeric_table.h"

#include "segment_site_catalog.h"
#include "segment_site_catalog_graph.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
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

CanonicalReachInput2 l_shape_input(const bool reversed,
                                   const double tool_radius = 0.5) {
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
                               {}, tool_radius);
}

bool native_table_fills_fields_one_through_nineteen() {
  const CanonicalReachInput2 input = l_shape_input(false);
  const MatNumericMatTable2 table = canonical_l_shape_mat_numeric_table(
      input, MatStationSpacingMm2::build(0.75), MatSagittaBoundMm2::build(0.02),
      32);
  if (table.nodes.size() != table.proposal.graph.node_ids.size() ||
      table.nodes.size() != 10 || table.proposal.graph.edges.size() != 9 ||
      table.neck_evidence.size() != 2 ||
      table.neck_cut_offsets != std::vector<std::int64_t>{0, 8, 16} ||
      table.neck_cut_edge_ids.size() != 16 ||
      table.center_domain_digest.size() != 32) {
    return false;
  }
  for (const auto &node : table.nodes) {
    if (!std::isfinite(node[0]) || !std::isfinite(node[1]) || node[2] != 0.0) {
      return false;
    }
  }
  for (std::size_t neck = 0; neck < table.neck_evidence.size(); ++neck) {
    if (table.neck_evidence[neck].find("mat-neck-plateau-v1") ==
        std::string::npos) {
      return false;
    }
    const auto begin =
        table.neck_cut_edge_ids.begin() + table.neck_cut_offsets[neck];
    const auto end =
        table.neck_cut_edge_ids.begin() + table.neck_cut_offsets[neck + 1];
    if (!std::is_sorted(begin, end) || std::adjacent_find(begin, end) != end) {
      return false;
    }
  }
  return true;
}

bool threshold_tool_radius_has_no_separating_neck() {
  const MatNumericMatTable2 table = canonical_l_shape_mat_numeric_table(
      l_shape_input(false, 1.0), MatStationSpacingMm2::build(0.75),
      MatSagittaBoundMm2::build(0.02), 32);
  return table.nodes.size() == 6 && table.proposal.graph.edges.size() == 5 &&
         table.neck_evidence.empty() &&
         table.neck_cut_offsets == std::vector<std::int64_t>{0} &&
         table.neck_cut_edge_ids.empty();
}

bool native_table_composes_existing_proposal_projection() {
  const CanonicalReachInput2 input = l_shape_input(false);
  const MatNumericMatTable2 table = canonical_l_shape_mat_numeric_table(
      input, MatStationSpacingMm2::build(0.75), MatSagittaBoundMm2::build(0.02),
      32);
  return table.proposal == canonical_l_shape_mat_numeric_proposal_table(
                               input, MatStationSpacingMm2::build(0.75),
                               MatSagittaBoundMm2::build(0.02), 32);
}

bool native_table_is_reversal_invariant() {
  const MatNumericMatTable2 first = canonical_l_shape_mat_numeric_table(
      l_shape_input(false), MatStationSpacingMm2::build(0.75),
      MatSagittaBoundMm2::build(0.02), 32);
  const MatNumericMatTable2 reversed = canonical_l_shape_mat_numeric_table(
      l_shape_input(true), MatStationSpacingMm2::build(0.75),
      MatSagittaBoundMm2::build(0.02), 32);
  return first == reversed;
}

bool mismatched_neck_graph_fails_loudly() {
  const MatClearanceProfileGraph2 l_bundle =
      canonical_l_shape_mat_clearance_graph(l_shape_input(false),
                                            CORE::BigRat(1, 4));
  const auto evidence = exact_neck_evidence_v1(l_bundle);
  try {
    static_cast<void>(
        numeric_neck_cut_table(MatExactGraph2{{}, {}, 0, 0}, evidence));
  } catch (const UnknownMatNumericNeckCutEdgeError &) {
    return true;
  }
  return false;
}

} // namespace

bool mat_numeric_table_gate() {
  return native_table_fills_fields_one_through_nineteen() &&
         threshold_tool_radius_has_no_separating_neck() &&
         native_table_composes_existing_proposal_projection() &&
         native_table_is_reversal_invariant() &&
         mismatched_neck_graph_fails_loudly();
}
