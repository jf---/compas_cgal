#pragma once

#include "segment_site_neck_clearance.h"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

class InvalidMatSamplingPolicyError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class InvalidMatSamplingCurveError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class InvalidMatProposalSampleError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class MatSamplingCardinalityError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class ConicSamplingLimitError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class MatStationSpacingMm2 {
public:
  static MatStationSpacingMm2 build(double value);
  double value() const noexcept;

private:
  explicit MatStationSpacingMm2(double value);
  double value_;
};

class MatSagittaBoundMm2 {
public:
  static MatSagittaBoundMm2 build(double value);
  double value() const noexcept;

private:
  explicit MatSagittaBoundMm2(double value);
  double value_;
};

class MatRationalSamplingCurve2 {
public:
  static MatRationalSamplingCurve2
  build(const MatExactGraphEdge2 &edge,
        const MatClearanceEdgeProfile2 &clearance_profile,
        RationalPolynomial x_coefficients, RationalPolynomial y_coefficients);

  const std::string &edge_id() const noexcept;
  const std::string &primitive_kind() const noexcept;
  const ExactAlgebraicKernel1::Algebraic_real_1 &
  lower_parameter() const noexcept;
  const ExactAlgebraicKernel1::Algebraic_real_1 &
  upper_parameter() const noexcept;
  const RationalPolynomial &x_coefficients() const noexcept;
  const RationalPolynomial &y_coefficients() const noexcept;

private:
  MatRationalSamplingCurve2(
      std::string edge_id, std::string primitive_kind,
      ExactAlgebraicKernel1::Algebraic_real_1 lower_parameter,
      ExactAlgebraicKernel1::Algebraic_real_1 upper_parameter,
      RationalPolynomial x_coefficients, RationalPolynomial y_coefficients);

  std::string edge_id_;
  std::string primitive_kind_;
  ExactAlgebraicKernel1::Algebraic_real_1 lower_parameter_;
  ExactAlgebraicKernel1::Algebraic_real_1 upper_parameter_;
  RationalPolynomial x_coefficients_;
  RationalPolynomial y_coefficients_;
};

class MatWorldXYProposalSample2 {
public:
  static MatWorldXYProposalSample2 build(std::string exact_parameter_id,
                                         double parameter, double x_mm,
                                         double y_mm);

  const std::string &exact_parameter_id() const noexcept;
  double parameter() const noexcept;
  double x_mm() const noexcept;
  double y_mm() const noexcept;

  bool operator==(const MatWorldXYProposalSample2 &) const = default;

private:
  MatWorldXYProposalSample2(std::string exact_parameter_id, double parameter,
                            double x_mm, double y_mm);

  std::string exact_parameter_id_;
  double parameter_;
  double x_mm_;
  double y_mm_;
};

std::vector<MatWorldXYProposalSample2>
proposal_samples(const MatRationalSamplingCurve2 &curve,
                 const MatStationSpacingMm2 &station_spacing,
                 const MatSagittaBoundMm2 &max_sagitta,
                 std::size_t max_refinement_depth);
