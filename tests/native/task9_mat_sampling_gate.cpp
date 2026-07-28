#include "segment_site_mat_sampling.h"

#include <cmath>
#include <iterator>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include <CGAL/Polynomial_traits_d.h>

namespace {

using Algebraic = ExactAlgebraicKernel1::Algebraic_real_1;

MatParameterEndpoint2 endpoint(const CORE::BigRat &parameter) {
  ExactAlgebraicKernel1 kernel;
  const auto algebraic = kernel.construct_algebraic_real_1_object()(parameter);
  return {
      algebraic,
      {
          algebraic_root_identity_v1(algebraic),
      },
      {
          true,
          false,
          false,
          {},
      },
  };
}

MatParameterEndpoint2 endpoint(const Algebraic &parameter) {
  return {
      parameter,
      {
          algebraic_root_identity_v1(parameter),
      },
      {
          true,
          false,
          false,
          {},
      },
  };
}

MatExactGraphEdge2 edge(const std::string &edge_id,
                        const std::string &primitive_kind,
                        const CORE::BigRat &lower, const CORE::BigRat &upper) {
  return {
      edge_id,
      primitive_kind,
      "dual-" + edge_id,
      "source-" + edge_id,
      "target-" + edge_id,
      endpoint(lower),
      endpoint(upper),
      {
          edge_id + "-site-a",
          edge_id + "-site-b",
      },
      {
          edge_id + "-parent-a",
          edge_id + "-parent-b",
      },
      0,
      true,
  };
}

MatExactGraphEdge2 algebraic_edge(const std::string &edge_id,
                                  const std::string &primitive_kind,
                                  const Algebraic &lower,
                                  const Algebraic &upper) {
  return {
      edge_id,
      primitive_kind,
      "dual-" + edge_id,
      "source-" + edge_id,
      "target-" + edge_id,
      endpoint(lower),
      endpoint(upper),
      {
          edge_id + "-site-a",
          edge_id + "-site-b",
      },
      {
          edge_id + "-parent-a",
          edge_id + "-parent-b",
      },
      0,
      true,
  };
}

MatClearanceEdgeProfile2 profile(const MatExactGraphEdge2 &edge,
                                 RationalPolynomial coefficients) {
  return MatClearanceEdgeProfile2::build(
      edge, MatSquaredClearanceFunction2::build(std::move(coefficients)));
}

MatRationalSamplingCurve2 line_curve(const std::string &edge_id = "line") {
  const MatExactGraphEdge2 graph_edge = edge(edge_id, "LINE", 0, 4);
  return MatRationalSamplingCurve2::build(graph_edge, profile(graph_edge, {1}),
                                          {0, 1}, {0});
}

MatRationalSamplingCurve2 parabola_curve() {
  const MatExactGraphEdge2 graph_edge = edge("parabola", "PARABOLA", -1, 1);
  return MatRationalSamplingCurve2::build(
      graph_edge, profile(graph_edge, {1, 0, 1}), {0, 1}, {0, 0, 1});
}

MatRationalSamplingCurve2 algebraic_endpoint_parabola_curve() {
  using Polynomial = ExactAlgebraicKernel1::Polynomial_1;
  const std::vector<ExactAlgebraicInteger1> coefficients{-2, 0, 1};
  const Polynomial polynomial =
      typename CGAL::Polynomial_traits_d<Polynomial>::Construct_polynomial()(
          coefficients.begin(), coefficients.end());
  ExactAlgebraicKernel1 kernel;
  std::vector<Algebraic> roots;
  kernel.solve_1_object()(polynomial, true, std::back_inserter(roots));
  if (roots.size() != 2) {
    throw InvalidMatSamplingCurveError(
        "algebraic sampling fixture requires both roots of t^2 - 2");
  }
  const MatExactGraphEdge2 graph_edge = algebraic_edge(
      "algebraic-parabola", "PARABOLA", roots.front(), roots.back());
  return MatRationalSamplingCurve2::build(graph_edge, profile(graph_edge, {2}),
                                          {0, 1}, {0, 0, 1});
}

bool line_samples_are_exact_barycentric() {
  const MatRationalSamplingCurve2 curve = line_curve();
  const auto samples = proposal_samples(curve, MatStationSpacingMm2::build(1.5),
                                        MatSagittaBoundMm2::build(0.25), 4);
  if (samples.size() != 4 ||
      samples.front().exact_parameter_id().find(
          "mat-line-barycentric-parameter-v1") == std::string::npos ||
      samples.back().exact_parameter_id().find(
          "mat-line-barycentric-parameter-v1") == std::string::npos) {
    return false;
  }
  const std::vector<double> expected_x{
      0.0,
      4.0 / 3.0,
      8.0 / 3.0,
      4.0,
  };
  for (std::size_t index = 0; index < samples.size(); ++index) {
    if (samples[index].x_mm() != expected_x[index] ||
        samples[index].y_mm() != 0.0 ||
        samples[index].parameter() != expected_x[index]) {
      return false;
    }
  }
  return proposal_samples(curve, MatStationSpacingMm2::build(1.5),
                          MatSagittaBoundMm2::build(0.25), 4) == samples;
}

bool line_parameter_identity_binds_edge() {
  const auto first =
      proposal_samples(line_curve("line-a"), MatStationSpacingMm2::build(10.0),
                       MatSagittaBoundMm2::build(1.0), 0);
  const auto second =
      proposal_samples(line_curve("line-b"), MatStationSpacingMm2::build(10.0),
                       MatSagittaBoundMm2::build(1.0), 0);
  return first.size() == 2 && second.size() == 2 &&
         first.front().exact_parameter_id().find(
             "mat-line-barycentric-parameter-v1") != std::string::npos &&
         first.front().exact_parameter_id() !=
             second.front().exact_parameter_id() &&
         first.back().exact_parameter_id() !=
             second.back().exact_parameter_id();
}

bool parabola_refinement_is_proposal_only() {
  const MatRationalSamplingCurve2 curve = parabola_curve();
  const auto coarse = proposal_samples(curve, MatStationSpacingMm2::build(10.0),
                                       MatSagittaBoundMm2::build(2.0), 4);
  const auto refined =
      proposal_samples(curve, MatStationSpacingMm2::build(10.0),
                       MatSagittaBoundMm2::build(0.3), 4);
  if (coarse.size() != 2 || refined.size() != 3 ||
      refined.front().parameter() != -1.0 ||
      refined.back().parameter() != 1.0 || curve.edge_id() != "parabola" ||
      refined.front().exact_parameter_id().find(
          "mat-parabola-algebraic-parameter-v1") == std::string::npos ||
      refined[1].exact_parameter_id().find(
          "mat-parabola-rational-parameter-v1") == std::string::npos ||
      refined.back().exact_parameter_id().find(
          "mat-parabola-algebraic-parameter-v1") == std::string::npos) {
    return false;
  }
  for (std::size_t index = 1; index < refined.size(); ++index) {
    if (refined[index - 1].parameter() >= refined[index].parameter()) {
      return false;
    }
  }
  return refined[1].parameter() == 0.0 && refined[1].x_mm() == 0.0 &&
         refined[1].y_mm() == 0.0 &&
         proposal_samples(curve, MatStationSpacingMm2::build(10.0),
                          MatSagittaBoundMm2::build(0.3), 4) == refined;
}

bool parabola_station_spacing_refines_independently() {
  const MatRationalSamplingCurve2 curve = parabola_curve();
  const MatStationSpacingMm2 spacing = MatStationSpacingMm2::build(1.0);
  const auto samples =
      proposal_samples(curve, spacing, MatSagittaBoundMm2::build(10.0), 4);
  if (samples.size() <= 2) {
    return false;
  }
  for (std::size_t index = 1; index < samples.size(); ++index) {
    if (std::hypot(samples[index].x_mm() - samples[index - 1].x_mm(),
                   samples[index].y_mm() - samples[index - 1].y_mm()) >
        spacing.value()) {
      return false;
    }
  }
  return true;
}

bool algebraic_endpoint_refinement_respects_reported_sagitta() {
  const MatRationalSamplingCurve2 curve = algebraic_endpoint_parabola_curve();
  const MatSagittaBoundMm2 sagitta_bound = MatSagittaBoundMm2::build(0.4);
  const auto samples = proposal_samples(
      curve, MatStationSpacingMm2::build(10.0), sagitta_bound, 6);
  for (std::size_t index = 1; index < samples.size(); ++index) {
    const auto &lower = samples[index - 1];
    const auto &upper = samples[index];
    const double midpoint_parameter =
        (lower.parameter() + upper.parameter()) / 2.0;
    const double midpoint_x = midpoint_parameter;
    const double midpoint_y = midpoint_parameter * midpoint_parameter;
    const double chord_midpoint_x = (lower.x_mm() + upper.x_mm()) / 2.0;
    const double chord_midpoint_y = (lower.y_mm() + upper.y_mm()) / 2.0;
    if (std::hypot(midpoint_x - chord_midpoint_x,
                   midpoint_y - chord_midpoint_y) > sagitta_bound.value()) {
      return false;
    }
  }
  return true;
}

bool conic_depth_cap_fails_loudly() {
  try {
    static_cast<void>(proposal_samples(parabola_curve(),
                                       MatStationSpacingMm2::build(10.0),
                                       MatSagittaBoundMm2::build(0.3), 0));
  } catch (const ConicSamplingLimitError &) {
    return true;
  }
  return false;
}

bool line_cardinality_fails_loudly() {
  try {
    static_cast<void>(proposal_samples(
        line_curve(),
        MatStationSpacingMm2::build(std::numeric_limits<double>::denorm_min()),
        MatSagittaBoundMm2::build(1.0), 0));
  } catch (const MatSamplingCardinalityError &) {
    return true;
  }
  return false;
}

bool malformed_proposal_sample_fails_loudly() {
  bool empty_identity_rejected = false;
  try {
    static_cast<void>(MatWorldXYProposalSample2::build("", 0.0, 0.0, 0.0));
  } catch (const InvalidMatProposalSampleError &) {
    empty_identity_rejected = true;
  }

  bool nonfinite_coordinate_rejected = false;
  try {
    static_cast<void>(MatWorldXYProposalSample2::build(
        "exact-parameter", 0.0, std::numeric_limits<double>::infinity(), 0.0));
  } catch (const InvalidMatProposalSampleError &) {
    nonfinite_coordinate_rejected = true;
  }
  return empty_identity_rejected && nonfinite_coordinate_rejected;
}

bool malformed_sampling_contracts_fail_loudly() {
  const auto invalid_policy = [](const double spacing, const double sagitta) {
    try {
      static_cast<void>(MatStationSpacingMm2::build(spacing));
      static_cast<void>(MatSagittaBoundMm2::build(sagitta));
    } catch (const InvalidMatSamplingPolicyError &) {
      return true;
    }
    return false;
  };
  if (!invalid_policy(0.0, 1.0) || !invalid_policy(1.0, 0.0) ||
      !invalid_policy(std::numeric_limits<double>::infinity(), 1.0) ||
      !invalid_policy(1.0, std::numeric_limits<double>::quiet_NaN())) {
    return false;
  }

  bool kind_rejected = false;
  try {
    const MatExactGraphEdge2 graph_edge = edge("kind", "CIRCLE", 0, 1);
    static_cast<void>(MatRationalSamplingCurve2::build(
        graph_edge, profile(graph_edge, {1}), {0, 1}, {0}));
  } catch (const InvalidMatSamplingCurveError &) {
    kind_rejected = true;
  }

  bool degenerate_rejected = false;
  try {
    const MatExactGraphEdge2 graph_edge = edge("degenerate", "LINE", 0, 1);
    static_cast<void>(MatRationalSamplingCurve2::build(
        graph_edge, profile(graph_edge, {1}), {0}, {0}));
  } catch (const InvalidMatSamplingCurveError &) {
    degenerate_rejected = true;
  }

  bool mismatch_rejected = false;
  try {
    const MatExactGraphEdge2 first = edge("first", "LINE", 0, 1);
    const MatExactGraphEdge2 second = edge("second", "LINE", 0, 1);
    static_cast<void>(MatRationalSamplingCurve2::build(
        first, profile(second, {1}), {0, 1}, {0}));
  } catch (const InvalidMatSamplingCurveError &) {
    mismatch_rejected = true;
  }
  return kind_rejected && degenerate_rejected && mismatch_rejected;
}

} // namespace

bool mat_sampling_gate() {
  return line_samples_are_exact_barycentric() &&
         line_parameter_identity_binds_edge() &&
         parabola_refinement_is_proposal_only() &&
         parabola_station_spacing_refines_independently() &&
         algebraic_endpoint_refinement_respects_reported_sagitta() &&
         conic_depth_cap_fails_loudly() && line_cardinality_fails_loudly() &&
         malformed_proposal_sample_fails_loudly() &&
         malformed_sampling_contracts_fail_loudly();
}
