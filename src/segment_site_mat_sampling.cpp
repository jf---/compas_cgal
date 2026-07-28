#include "segment_site_mat_sampling.h"

#include "canonical_encoding.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace {

using Algebraic = ExactAlgebraicKernel1::Algebraic_real_1;

bool finite_positive(const double value) {
  return std::isfinite(value) && value > 0.0;
}

void require_finite_sample(const double parameter, const double x_mm,
                           const double y_mm) {
  if (!std::isfinite(parameter) || !std::isfinite(x_mm) ||
      !std::isfinite(y_mm)) {
    throw InvalidMatProposalSampleError(
        "MAT sampling curve has no finite reporting coordinates");
  }
}

double evaluate_reported(const RationalPolynomial &coefficients,
                         const double parameter) {
  double result = 0.0;
  for (auto coefficient = coefficients.rbegin();
       coefficient != coefficients.rend(); ++coefficient) {
    result = result * parameter + CGAL::to_double(*coefficient);
  }
  return result;
}

MatWorldXYProposalSample2
algebraic_sample(const MatRationalSamplingCurve2 &curve,
                 const Algebraic &parameter) {
  const double reported_parameter = parameter.to_double();
  const double x_mm =
      evaluate_reported(curve.x_coefficients(), reported_parameter);
  const double y_mm =
      evaluate_reported(curve.y_coefficients(), reported_parameter);
  require_finite_sample(reported_parameter, x_mm, y_mm);
  return MatWorldXYProposalSample2::build(
      curve.edge_id(),
      canonical_encode_tagged_union(
          "mat-parabola-algebraic-parameter-v1",
          canonical_encode_component_map({
              {
                  "edge-id",
                  curve.edge_id(),
              },
              {
                  "parameter-root-id",
                  algebraic_root_identity_v1(parameter),
              },
          })),
      reported_parameter, x_mm, y_mm);
}

MatWorldXYProposalSample2
rational_parabola_sample(const MatRationalSamplingCurve2 &curve,
                         const CORE::BigRat &parameter) {
  const double reported_parameter = CGAL::to_double(parameter);
  const double x_mm = CGAL::to_double(
      evaluate_rational_polynomial(curve.x_coefficients(), parameter));
  const double y_mm = CGAL::to_double(
      evaluate_rational_polynomial(curve.y_coefficients(), parameter));
  require_finite_sample(reported_parameter, x_mm, y_mm);
  return MatWorldXYProposalSample2::build(
      curve.edge_id(),
      canonical_encode_tagged_union(
          "mat-parabola-rational-parameter-v1",
          canonical_encode_component_map({
              {
                  "edge-id",
                  curve.edge_id(),
              },
              {
                  "parameter",
                  canonical_encode_rational(parameter),
              },
          })),
      reported_parameter, x_mm, y_mm);
}

MatWorldXYProposalSample2 line_sample(const MatRationalSamplingCurve2 &curve,
                                      const CORE::BigRat &barycentric_parameter,
                                      const MatWorldXYProposalSample2 &lower,
                                      const MatWorldXYProposalSample2 &upper) {
  const double fraction = CGAL::to_double(barycentric_parameter);
  const double parameter =
      lower.parameter() + (upper.parameter() - lower.parameter()) * fraction;
  const double x_mm = lower.x_mm() + (upper.x_mm() - lower.x_mm()) * fraction;
  const double y_mm = lower.y_mm() + (upper.y_mm() - lower.y_mm()) * fraction;
  require_finite_sample(parameter, x_mm, y_mm);
  return MatWorldXYProposalSample2::build(
      curve.edge_id(),
      canonical_encode_tagged_union(
          "mat-line-barycentric-parameter-v1",
          canonical_encode_component_map({
              {
                  "edge-id",
                  curve.edge_id(),
              },
              {
                  "fraction",
                  canonical_encode_rational(barycentric_parameter),
              },
          })),
      parameter, x_mm, y_mm);
}

double distance(const MatWorldXYProposalSample2 &first,
                const MatWorldXYProposalSample2 &second) {
  return std::hypot(second.x_mm() - first.x_mm(), second.y_mm() - first.y_mm());
}

// For a quadratic coordinate map, the parameter-matched chord-deviation
// vector is proportional to lambda * (1 - lambda), so its norm peaks at the
// parameter midpoint. The exact rational subdivision separator is not
// generally that midpoint.
double
reported_quadratic_chord_deviation(const MatRationalSamplingCurve2 &curve,
                                   const MatWorldXYProposalSample2 &lower,
                                   const MatWorldXYProposalSample2 &upper) {
  const double midpoint_parameter =
      std::midpoint(lower.parameter(), upper.parameter());
  const double midpoint_x =
      evaluate_reported(curve.x_coefficients(), midpoint_parameter);
  const double midpoint_y =
      evaluate_reported(curve.y_coefficients(), midpoint_parameter);
  const double chord_midpoint_x = std::midpoint(lower.x_mm(), upper.x_mm());
  const double chord_midpoint_y = std::midpoint(lower.y_mm(), upper.y_mm());
  require_finite_sample(midpoint_parameter, midpoint_x, midpoint_y);
  return std::hypot(midpoint_x - chord_midpoint_x,
                    midpoint_y - chord_midpoint_y);
}

CORE::BigRat exact_split_parameter(const Algebraic &lower,
                                   const Algebraic &upper) {
  if (lower.is_rational() && upper.is_rational()) {
    return (lower.rational() + upper.rational()) / 2;
  }
  ExactAlgebraicKernel1 kernel;
  return kernel.bound_between_1_object()(lower, upper);
}

void refine_parabola(const MatRationalSamplingCurve2 &curve,
                     const Algebraic &lower_parameter,
                     const MatWorldXYProposalSample2 &lower,
                     const Algebraic &upper_parameter,
                     const MatWorldXYProposalSample2 &upper,
                     const MatStationSpacingMm2 &station_spacing,
                     const MatSagittaBoundMm2 &max_sagitta,
                     const std::size_t max_refinement_depth,
                     const std::size_t depth,
                     std::vector<MatWorldXYProposalSample2> &samples) {
  if (distance(lower, upper) <= station_spacing.value() &&
      reported_quadratic_chord_deviation(curve, lower, upper) <=
          max_sagitta.value()) {
    samples.push_back(upper);
    return;
  }
  if (depth >= max_refinement_depth) {
    throw ConicSamplingLimitError(
        "parabola sampling exceeded refinement depth on edge " +
        curve.edge_id());
  }
  const CORE::BigRat split_parameter =
      exact_split_parameter(lower_parameter, upper_parameter);
  const MatWorldXYProposalSample2 middle =
      rational_parabola_sample(curve, split_parameter);
  ExactAlgebraicKernel1 kernel;
  const Algebraic algebraic_split =
      kernel.construct_algebraic_real_1_object()(split_parameter);
  refine_parabola(curve, lower_parameter, lower, algebraic_split, middle,
                  station_spacing, max_sagitta, max_refinement_depth, depth + 1,
                  samples);
  refine_parabola(curve, algebraic_split, middle, upper_parameter, upper,
                  station_spacing, max_sagitta, max_refinement_depth, depth + 1,
                  samples);
}

std::vector<MatWorldXYProposalSample2>
line_samples(const MatRationalSamplingCurve2 &curve,
             const MatStationSpacingMm2 &station_spacing) {
  const MatWorldXYProposalSample2 lower =
      algebraic_sample(curve, curve.lower_parameter());
  const MatWorldXYProposalSample2 upper =
      algebraic_sample(curve, curve.upper_parameter());
  const long double raw_segment_count =
      std::ceil(static_cast<long double>(distance(lower, upper)) /
                static_cast<long double>(station_spacing.value()));
  if (!std::isfinite(raw_segment_count) ||
      raw_segment_count > static_cast<long double>(
                              std::numeric_limits<std::size_t>::max() - 1)) {
    throw MatSamplingCardinalityError(
        "line proposal sample count exceeds native address space");
  }
  const std::size_t segment_count =
      std::max<std::size_t>(1, static_cast<std::size_t>(raw_segment_count));
  std::vector<MatWorldXYProposalSample2> samples;
  samples.reserve(segment_count + 1);
  for (std::size_t index = 0; index <= segment_count; ++index) {
    samples.push_back(
        line_sample(curve,
                    CORE::BigRat(ExactAlgebraicInteger1(index),
                                 ExactAlgebraicInteger1(segment_count)),
                    lower, upper));
  }
  return samples;
}

std::vector<MatWorldXYProposalSample2>
parabola_samples(const MatRationalSamplingCurve2 &curve,
                 const MatStationSpacingMm2 &station_spacing,
                 const MatSagittaBoundMm2 &max_sagitta,
                 const std::size_t max_refinement_depth) {
  const MatWorldXYProposalSample2 lower =
      algebraic_sample(curve, curve.lower_parameter());
  const MatWorldXYProposalSample2 upper =
      algebraic_sample(curve, curve.upper_parameter());
  std::vector<MatWorldXYProposalSample2> samples{
      lower,
  };
  refine_parabola(curve, curve.lower_parameter(), lower,
                  curve.upper_parameter(), upper, station_spacing, max_sagitta,
                  max_refinement_depth, 0, samples);
  return samples;
}

} // namespace

MatStationSpacingMm2::MatStationSpacingMm2(const double value)
    : value_(value) {}

MatStationSpacingMm2 MatStationSpacingMm2::build(const double value) {
  if (!finite_positive(value)) {
    throw InvalidMatSamplingPolicyError(
        "MAT station spacing must be finite and positive");
  }
  return MatStationSpacingMm2(value);
}

double MatStationSpacingMm2::value() const noexcept { return value_; }

MatSagittaBoundMm2::MatSagittaBoundMm2(const double value) : value_(value) {}

MatSagittaBoundMm2 MatSagittaBoundMm2::build(const double value) {
  if (!finite_positive(value)) {
    throw InvalidMatSamplingPolicyError(
        "MAT sagitta bound must be finite and positive");
  }
  return MatSagittaBoundMm2(value);
}

double MatSagittaBoundMm2::value() const noexcept { return value_; }

MatRationalSamplingCurve2::MatRationalSamplingCurve2(
    std::string edge_id, std::string primitive_kind, Algebraic lower_parameter,
    Algebraic upper_parameter, RationalPolynomial x_coefficients,
    RationalPolynomial y_coefficients)
    : edge_id_(std::move(edge_id)), primitive_kind_(std::move(primitive_kind)),
      lower_parameter_(std::move(lower_parameter)),
      upper_parameter_(std::move(upper_parameter)),
      x_coefficients_(std::move(x_coefficients)),
      y_coefficients_(std::move(y_coefficients)) {}

MatRationalSamplingCurve2 MatRationalSamplingCurve2::build(
    const MatExactGraphEdge2 &edge,
    const MatClearanceEdgeProfile2 &clearance_profile,
    RationalPolynomial x_coefficients, RationalPolynomial y_coefficients) {
  if (!edge.source_endpoint.parameter.has_value() ||
      !edge.target_endpoint.parameter.has_value() ||
      clearance_profile.edge_id() != edge.edge_id ||
      clearance_profile.defining_site_ids() != edge.generator_site_ids ||
      algebraic_root_identity_v1(*clearance_profile.lower().parameter) !=
          algebraic_root_identity_v1(*edge.source_endpoint.parameter) ||
      algebraic_root_identity_v1(*clearance_profile.upper().parameter) !=
          algebraic_root_identity_v1(*edge.target_endpoint.parameter)) {
    throw InvalidMatSamplingCurveError(
        "MAT sampling curve does not match its exact edge/profile");
  }
  if (x_coefficients.empty() || y_coefficients.empty()) {
    throw InvalidMatSamplingCurveError(
        "MAT sampling curve requires two coordinate polynomials");
  }
  trim(x_coefficients);
  trim(y_coefficients);
  const bool is_line = edge.primitive_kind == "LINE";
  const bool is_parabola = edge.primitive_kind == "PARABOLA";
  const bool line_shape =
      x_coefficients.size() <= 2 && y_coefficients.size() <= 2 &&
      (x_coefficients.size() == 2 || y_coefficients.size() == 2);
  const bool parabola_shape =
      x_coefficients.size() <= 3 && y_coefficients.size() <= 3 &&
      (x_coefficients.size() == 3 || y_coefficients.size() == 3);
  if ((!is_line && !is_parabola) || (is_line && !line_shape) ||
      (is_parabola && !parabola_shape)) {
    throw InvalidMatSamplingCurveError(
        "MAT sampling curve kind disagrees with coordinate degree");
  }
  ExactAlgebraicKernel1 kernel;
  if (kernel.compare_1_object()(*edge.source_endpoint.parameter,
                                *edge.target_endpoint.parameter) !=
      CGAL::SMALLER) {
    throw InvalidMatSamplingCurveError(
        "MAT sampling curve parameter interval is not increasing");
  }
  return MatRationalSamplingCurve2(
      edge.edge_id, edge.primitive_kind, *edge.source_endpoint.parameter,
      *edge.target_endpoint.parameter, std::move(x_coefficients),
      std::move(y_coefficients));
}

const std::string &MatRationalSamplingCurve2::edge_id() const noexcept {
  return edge_id_;
}

const std::string &MatRationalSamplingCurve2::primitive_kind() const noexcept {
  return primitive_kind_;
}

const Algebraic &MatRationalSamplingCurve2::lower_parameter() const noexcept {
  return lower_parameter_;
}

const Algebraic &MatRationalSamplingCurve2::upper_parameter() const noexcept {
  return upper_parameter_;
}

const RationalPolynomial &
MatRationalSamplingCurve2::x_coefficients() const noexcept {
  return x_coefficients_;
}

const RationalPolynomial &
MatRationalSamplingCurve2::y_coefficients() const noexcept {
  return y_coefficients_;
}

MatWorldXYProposalSample2::MatWorldXYProposalSample2(
    std::string edge_id, std::string exact_parameter_id, const double parameter,
    const double x_mm, const double y_mm)
    : edge_id_(std::move(edge_id)),
      exact_parameter_id_(std::move(exact_parameter_id)), parameter_(parameter),
      x_mm_(x_mm), y_mm_(y_mm) {}

MatWorldXYProposalSample2
MatWorldXYProposalSample2::build(std::string exact_parameter_id,
                                 const double parameter, const double x_mm,
                                 const double y_mm) {
  if (exact_parameter_id.empty()) {
    throw InvalidMatProposalSampleError(
        "MAT proposal sample has no exact parameter identity");
  }
  require_finite_sample(parameter, x_mm, y_mm);
  return MatWorldXYProposalSample2("", std::move(exact_parameter_id), parameter,
                                   x_mm, y_mm);
}

MatWorldXYProposalSample2 MatWorldXYProposalSample2::build(
    std::string edge_id, std::string exact_parameter_id, const double parameter,
    const double x_mm, const double y_mm) {
  if (edge_id.empty()) {
    throw InvalidMatProposalSampleError(
        "MAT proposal sample has no exact edge owner");
  }
  if (exact_parameter_id.empty()) {
    throw InvalidMatProposalSampleError(
        "MAT proposal sample has no exact parameter identity");
  }
  require_finite_sample(parameter, x_mm, y_mm);
  return MatWorldXYProposalSample2(
      std::move(edge_id), std::move(exact_parameter_id), parameter, x_mm, y_mm);
}

const std::string &MatWorldXYProposalSample2::edge_id() const noexcept {
  return edge_id_;
}

const std::string &
MatWorldXYProposalSample2::exact_parameter_id() const noexcept {
  return exact_parameter_id_;
}

double MatWorldXYProposalSample2::parameter() const noexcept {
  return parameter_;
}

double MatWorldXYProposalSample2::x_mm() const noexcept { return x_mm_; }

double MatWorldXYProposalSample2::y_mm() const noexcept { return y_mm_; }

MatProposalSamplingRun2::MatProposalSamplingRun2(
    std::string edge_id, std::vector<MatWorldXYProposalSample2> samples)
    : edge_id_(std::move(edge_id)), samples_(std::move(samples)) {}

MatProposalSamplingRun2
MatProposalSamplingRun2::build(std::string edge_id,
                               std::vector<MatWorldXYProposalSample2> samples) {
  if (edge_id.empty() || samples.size() < 2) {
    throw InvalidMatProposalSamplingRunError(
        "MAT proposal sampling run requires an edge and both endpoints");
  }
  for (std::size_t index = 0; index < samples.size(); ++index) {
    if (samples[index].edge_id() != edge_id ||
        (index > 0 &&
         samples[index - 1].parameter() >= samples[index].parameter())) {
      throw InvalidMatProposalSamplingRunError(
          "MAT proposal sampling run is not bound and parameter ordered");
    }
  }
  return MatProposalSamplingRun2(std::move(edge_id), std::move(samples));
}

const std::string &MatProposalSamplingRun2::edge_id() const noexcept {
  return edge_id_;
}

const std::vector<MatWorldXYProposalSample2> &
MatProposalSamplingRun2::samples() const noexcept {
  return samples_;
}

std::vector<MatWorldXYProposalSample2>
proposal_samples(const MatRationalSamplingCurve2 &curve,
                 const MatStationSpacingMm2 &station_spacing,
                 const MatSagittaBoundMm2 &max_sagitta,
                 const std::size_t max_refinement_depth) {
  if (curve.primitive_kind() == "LINE") {
    return line_samples(curve, station_spacing);
  }
  return parabola_samples(curve, station_spacing, max_sagitta,
                          max_refinement_depth);
}

MatProposalSamplingRun2
proposal_sampling_run(const MatRationalSamplingCurve2 &curve,
                      const MatStationSpacingMm2 &station_spacing,
                      const MatSagittaBoundMm2 &max_sagitta,
                      const std::size_t max_refinement_depth) {
  return MatProposalSamplingRun2::build(
      curve.edge_id(), proposal_samples(curve, station_spacing, max_sagitta,
                                        max_refinement_depth));
}
