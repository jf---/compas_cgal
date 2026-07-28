#include "segment_site_neck_classification.h"

#include "canonical_encoding.h"
#include "continuous_tea_2/sha256.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace {

std::vector<CORE::BigRat>
decode_boundaries(const std::vector<std::string> &boundary_bytes) {
  if (boundary_bytes.size() >
      static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max())) {
    throw InvalidMatNeckWidthBoundariesError(
        "neck width boundary count exceeds int64 class identity");
  }
  std::vector<CORE::BigRat> boundaries;
  boundaries.reserve(boundary_bytes.size());
  for (std::size_t index = 0; index < boundary_bytes.size(); ++index) {
    try {
      boundaries.push_back(canonical_decode_rational(boundary_bytes[index]));
    } catch (const InvalidCanonicalEncodingError &) {
      throw InvalidMatNeckWidthBoundariesError(
          "neck width boundary at index " + std::to_string(index) +
          " is not canonical ExactRationalV1");
    }
    if (boundaries.back() < 0) {
      throw InvalidMatNeckWidthBoundariesError(
          "neck squared-width boundary at index " + std::to_string(index) +
          " is negative");
    }
    if (index > 0 && boundaries[index - 1] >= boundaries[index]) {
      throw InvalidMatNeckWidthBoundariesError(
          "neck squared-width boundaries are not strictly increasing at "
          "index " +
          std::to_string(index));
    }
  }
  return boundaries;
}

std::int64_t comparison_id(const CGAL::Comparison_result comparison) {
  if (comparison == CGAL::SMALLER) {
    return -1;
  }
  if (comparison == CGAL::LARGER) {
    return 1;
  }
  return 0;
}

std::string boundary_comparison_bytes(const std::string &boundary,
                                      const std::int64_t comparison) {
  return canonical_encode_tagged_union(
      "mat-neck-width-boundary-comparison-v1",
      canonical_encode_component_map({
          {
              "boundary",
              boundary,
          },
          {
              "comparison",
              canonical_encode_integer(comparison),
          },
      }));
}

std::pair<std::int64_t, std::string>
classify_one(const MatNeckEvidenceV1 &evidence,
             const std::string &evidence_bytes,
             const std::vector<std::string> &boundary_bytes,
             const std::vector<CORE::BigRat> &boundaries) {
  ExactAlgebraicKernel1 kernel;
  std::vector<std::string> comparisons;
  comparisons.reserve(boundaries.size());
  std::int64_t class_id = static_cast<std::int64_t>(boundaries.size());
  for (std::size_t index = 0; index < boundaries.size(); ++index) {
    const CGAL::Comparison_result comparison = kernel.compare_1_object()(
        evidence.evidence().squared_width().value(), boundaries[index]);
    comparisons.push_back(boundary_comparison_bytes(boundary_bytes[index],
                                                    comparison_id(comparison)));
    if (class_id == static_cast<std::int64_t>(boundaries.size()) &&
        comparison != CGAL::LARGER) {
      class_id = static_cast<std::int64_t>(index);
    }
  }
  const std::string certificate = canonical_encode_tagged_union(
      "mat-neck-width-comparison-v1",
      canonical_encode_component_map({
          {
              "neck-evidence-digest",
              sha256_bytes(evidence_bytes),
          },
          {
              "squared-width-root-id",
              evidence.evidence().squared_width().root_id(),
          },
          {
              "boundary-comparisons",
              canonical_encode_sequence(comparisons),
          },
          {
              "class-id",
              canonical_encode_integer(class_id),
          },
      }));
  return {
      class_id,
      certificate,
  };
}

} // namespace

std::pair<std::vector<std::int64_t>, std::vector<std::string>>
validate_and_classify_necks_v1(
    const MatClearanceProfileGraph2 &bundle,
    const std::vector<std::string> &neck_evidence,
    const std::vector<std::string> &squared_width_boundaries) {
  const std::vector<CORE::BigRat> boundaries =
      decode_boundaries(squared_width_boundaries);
  const std::vector<MatNeckEvidenceV1> verified =
      verify_neck_evidence_v1(bundle, neck_evidence);
  std::vector<std::int64_t> classes;
  std::vector<std::string> certificates;
  classes.reserve(verified.size());
  certificates.reserve(verified.size());
  for (std::size_t index = 0; index < verified.size(); ++index) {
    auto [class_id, certificate] =
        classify_one(verified[index], neck_evidence[index],
                     squared_width_boundaries, boundaries);
    classes.push_back(class_id);
    certificates.push_back(std::move(certificate));
  }
  return {
      std::move(classes),
      std::move(certificates),
  };
}
