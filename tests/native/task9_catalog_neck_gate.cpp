#include "segment_site_catalog_neck.h"

#include "segment_site_catalog.h"
#include "segment_site_catalog_graph.h"
#include "segment_site_graph_csr.h"
#include "segment_site_neck_evidence.h"
#include "segment_site_neck_evidence_bytes.h"

#include <algorithm>
#include <array>
#include <cstddef>
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

CanonicalReachInput2 l_shape_input(const bool transformed) {
  return canonical_reach_input(transformed ? matrix({
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

bool endpoints_equal(const MatParameterEndpoint2 &lhs,
                     const MatParameterEndpoint2 &rhs) {
  return lhs.parameter.has_value() == rhs.parameter.has_value() &&
         (!lhs.parameter.has_value() ||
          algebraic_root_identity_v1(*lhs.parameter) ==
              algebraic_root_identity_v1(*rhs.parameter)) &&
         lhs.provenance_ids == rhs.provenance_ids &&
         lhs.exact_evidence == rhs.exact_evidence;
}

bool profiles_equal(const MatClearanceEdgeProfile2 &lhs,
                    const MatClearanceEdgeProfile2 &rhs) {
  return lhs.edge_id() == rhs.edge_id() &&
         lhs.defining_site_ids() == rhs.defining_site_ids() &&
         endpoints_equal(lhs.lower(), rhs.lower()) &&
         endpoints_equal(lhs.upper(), rhs.upper()) &&
         lhs.squared_clearance().coefficients() ==
             rhs.squared_clearance().coefficients();
}

bool l_shape_profiles_are_production_bound() {
  const CanonicalReachInput2 input = l_shape_input(false);
  const CanonicalMatSiteCatalog2 catalog = canonical_mat_site_catalog(input);
  const MatClearanceProfileGraph2 bundle =
      canonical_l_shape_mat_clearance_graph(input, CORE::BigRat(0));
  const MatExactGraph2 graph =
      canonical_l_shape_mat_graph(input, CORE::BigRat(0));
  if (numeric_graph_table(bundle.graph(), catalog) !=
      numeric_graph_table(graph, catalog)) {
    return false;
  }
  std::size_t constants = 0;
  std::size_t quadratics = 0;
  std::size_t quartics = 0;
  std::size_t strict_minima = 0;
  for (std::size_t index = 0; index < bundle.profiles().size(); ++index) {
    const auto &profile = bundle.profiles()[index];
    if (profile.edge_id() != bundle.graph().edges[index].edge_id) {
      return false;
    }
    const std::size_t degree =
        profile.squared_clearance().coefficients().size() - 1;
    if ((degree == 0 &&
         profile.squared_clearance().coefficients() != RationalPolynomial{1}) ||
        (degree == 2 &&
         profile.squared_clearance().coefficients() != RationalPolynomial{
                                                           0,
                                                           0,
                                                           1,
                                                       })) {
      return false;
    }
    constants += degree == 0;
    quadratics += degree == 2;
    quartics += degree == 4;
    strict_minima += strict_edge_clearance_minima(profile).size();
  }
  return bundle.profiles().size() == 9 && constants == 2 && quadratics == 5 &&
         quartics == 2 && strict_minima == 0;
}

bool positive_radius_profiles_match_clips() {
  const CanonicalReachInput2 input = l_shape_input(false);
  const CanonicalMatSiteCatalog2 catalog = canonical_mat_site_catalog(input);
  const MatClearanceProfileGraph2 bundle =
      canonical_l_shape_mat_clearance_graph(input, CORE::BigRat(1));
  if (numeric_graph_table(bundle.graph(), catalog) !=
      numeric_graph_table(canonical_l_shape_mat_graph(input, CORE::BigRat(1)),
                          catalog)) {
    return false;
  }
  std::size_t clearance_endpoints = 0;
  ExactAlgebraicKernel1 kernel;
  for (const MatClearanceEdgeProfile2 &profile : bundle.profiles()) {
    for (const MatParameterEndpoint2 *endpoint :
         {&profile.lower(), &profile.upper()}) {
      if (!endpoint->exact_evidence.clearance_root) {
        continue;
      }
      ++clearance_endpoints;
      const MatExactSquaredWidth2 width = MatExactSquaredWidth2::from_clearance(
          profile.squared_clearance(), *endpoint->parameter);
      if (kernel.compare_1_object()(width.value(), CORE::BigRat(4)) !=
              CGAL::EQUAL ||
          width.root_id() != algebraic_root_id_v1({-4, 1}, 0)) {
        return false;
      }
    }
  }
  return clearance_endpoints > 0;
}

bool l_shape_profiles_are_reversal_invariant() {
  const MatClearanceProfileGraph2 first = canonical_l_shape_mat_clearance_graph(
      l_shape_input(false), CORE::BigRat(0));
  const MatClearanceProfileGraph2 reversed =
      canonical_l_shape_mat_clearance_graph(l_shape_input(true),
                                            CORE::BigRat(0));
  if (first.profiles().size() != reversed.profiles().size()) {
    return false;
  }
  for (std::size_t index = 0; index < first.profiles().size(); ++index) {
    if (!profiles_equal(first.profiles()[index], reversed.profiles()[index])) {
      return false;
    }
  }
  return true;
}

bool l_shape_necks_are_two_exact_plateaus() {
  const MatClearanceProfileGraph2 first_bundle =
      canonical_l_shape_mat_clearance_graph(l_shape_input(false),
                                            CORE::BigRat(0));
  const MatClearanceProfileGraph2 reversed_bundle =
      canonical_l_shape_mat_clearance_graph(l_shape_input(true),
                                            CORE::BigRat(0));
  const auto first = exact_neck_evidence_v1(first_bundle);
  const auto reversed = exact_neck_evidence_v1(reversed_bundle);
  if (first.size() != 2 || reversed.size() != first.size()) {
    return false;
  }
  ExactAlgebraicKernel1 kernel;
  std::vector<std::string> records;
  for (std::size_t index = 0; index < first.size(); ++index) {
    const MatExactNeckEvidence2 &first_evidence = first[index].evidence();
    const MatExactNeckEvidence2 &reversed_evidence = reversed[index].evidence();
    const auto *location =
        std::get_if<MatPlateauNeckLocation2>(&first_evidence.location());
    const auto *reversed_location =
        std::get_if<MatPlateauNeckLocation2>(&reversed_evidence.location());
    if (location == nullptr || reversed_location == nullptr ||
        location->edge_ids().size() != 1 || location->node_ids().size() != 2 ||
        first_evidence.owner_id() != reversed_evidence.owner_id() ||
        first_evidence.defining_site_ids() !=
            reversed_evidence.defining_site_ids() ||
        first_evidence.squared_width().root_id() !=
            reversed_evidence.squared_width().root_id() ||
        kernel.compare_1_object()(first_evidence.squared_width().value(),
                                  CORE::BigRat(4)) != CGAL::EQUAL ||
        first_evidence.separating_cut() != reversed_evidence.separating_cut() ||
        location->edge_ids() != reversed_location->edge_ids() ||
        location->node_ids() != reversed_location->node_ids() ||
        first_evidence.separating_cut().edge_partitions().size() < 2 ||
        first[index].canonical_bytes() != reversed[index].canonical_bytes() ||
        first[index].canonical_digest() != reversed[index].canonical_digest()) {
      return false;
    }
    records.push_back(first[index].canonical_bytes());
  }
  return verify_neck_evidence_v1(first_bundle, records).size() ==
             first.size() &&
         verify_neck_evidence_v1(reversed_bundle, records).size() ==
             reversed.size();
}

bool malformed_profile_graph_is_rejected() {
  const MatClearanceProfileGraph2 valid = canonical_l_shape_mat_clearance_graph(
      l_shape_input(false), CORE::BigRat(0));
  bool missing_rejected = false;
  try {
    std::vector<MatClearanceEdgeProfile2> missing = valid.profiles();
    missing.pop_back();
    static_cast<void>(
        MatClearanceProfileGraph2::build(valid.graph(), std::move(missing)));
  } catch (const IncompleteMatClearanceProfileGraphError &) {
    missing_rejected = true;
  }

  bool order_rejected = false;
  try {
    std::vector<MatClearanceEdgeProfile2> reversed = valid.profiles();
    std::swap(reversed[0], reversed[1]);
    static_cast<void>(
        MatClearanceProfileGraph2::build(valid.graph(), std::move(reversed)));
  } catch (const IncompleteMatClearanceProfileGraphError &) {
    order_rejected = true;
  }
  bool graph_order_rejected = false;
  try {
    MatExactGraph2 graph = valid.graph();
    std::vector<MatClearanceEdgeProfile2> profiles = valid.profiles();
    std::swap(graph.edges[0], graph.edges[1]);
    std::swap(profiles[0], profiles[1]);
    static_cast<void>(MatClearanceProfileGraph2::build(std::move(graph),
                                                       std::move(profiles)));
  } catch (const IncompleteMatClearanceProfileGraphError &) {
    graph_order_rejected = true;
  }
  return missing_rejected && order_rejected && graph_order_rejected;
}

} // namespace

bool catalog_neck_gate() {
  return l_shape_profiles_are_production_bound() &&
         positive_radius_profiles_match_clips() &&
         l_shape_profiles_are_reversal_invariant() &&
         l_shape_necks_are_two_exact_plateaus() &&
         malformed_profile_graph_is_rejected();
}
