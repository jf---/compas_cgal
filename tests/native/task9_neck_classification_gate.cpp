#include "canonical_encoding.h"
#include "continuous_tea_2/sha256.h"
#include "segment_site_neck_classification.h"
#include "segment_site_neck_evidence_bytes.h"

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace {

MatExactGraphNode2 node(const std::string &node_id) {
  return {
      node_id,
      {"provenance-" + node_id},
      {"generator-" + node_id},
      {"parent-" + node_id},
  };
}

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

MatExactGraphEdge2 edge(const std::string &edge_id,
                        const std::string &source_node_id,
                        const std::string &target_node_id) {
  return {
      edge_id,
      "LINE",
      "dual-" + edge_id,
      source_node_id,
      target_node_id,
      endpoint(0),
      endpoint(1),
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

MatClearanceProfileGraph2 two_width_bundle() {
  std::vector<MatExactGraphEdge2> edges{
      edge("e0", "n0", "n1"),
      edge("e1", "n2", "n3"),
  };
  std::vector<MatClearanceEdgeProfile2> profiles{
      MatClearanceEdgeProfile2::build(edges[0],
                                      MatSquaredClearanceFunction2::build({
                                          CORE::BigRat(5, 4),
                                          -1,
                                          1,
                                      })),
      MatClearanceEdgeProfile2::build(edges[1],
                                      MatSquaredClearanceFunction2::build({
                                          CORE::BigRat(9, 4),
                                          -1,
                                          1,
                                      })),
  };
  return MatClearanceProfileGraph2::build(
      {
          {
              node("n0"),
              node("n1"),
              node("n2"),
              node("n3"),
          },
          std::move(edges),
          0,
          0,
      },
      std::move(profiles));
}

std::vector<std::string>
record_bytes(const std::vector<MatNeckEvidenceV1> &evidence) {
  std::vector<std::string> result;
  result.reserve(evidence.size());
  for (const MatNeckEvidenceV1 &record : evidence) {
    result.push_back(record.canonical_bytes());
  }
  return result;
}

bool rational_bytes_match_python_golden() {
  const unsigned char golden_bytes[]{
      0x43, 0x43, 0x41, 0x4e, 0x00, 0x01, 0x52, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x22, 0x43, 0x43, 0x41, 0x4e, 0x00,
      0x01, 0x49, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02,
      0x01, 0x03, 0x43, 0x43, 0x41, 0x4e, 0x00, 0x01, 0x49, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0x00, 0x05,
  };
  const std::string golden(reinterpret_cast<const char *>(golden_bytes),
                           sizeof(golden_bytes));
  return canonical_encode_rational(CORE::BigRat(-3, 5)) == golden &&
         canonical_decode_rational(golden) == CORE::BigRat(-3, 5);
}

bool malformed_rational_bytes_fail_loudly() {
  const std::string valid = canonical_encode_rational(CORE::BigRat(-3, 5));
  if (valid.size() != 49) {
    return false;
  }
  const auto rejected = [](const std::string &record) {
    try {
      static_cast<void>(canonical_decode_rational(record));
    } catch (const InvalidCanonicalEncodingError &) {
      return true;
    }
    return false;
  };

  std::string truncated = valid;
  truncated.pop_back();
  std::string trailing = valid;
  trailing.push_back('\0');
  std::string wrong_child_kind = valid;
  wrong_child_kind[21] = 'B';
  std::string invalid_sign = valid;
  invalid_sign[30] = static_cast<char>(2);
  std::string nonminimal_magnitude = valid;
  nonminimal_magnitude[30] = 0;
  nonminimal_magnitude[31] = 0;
  std::string negative_denominator = valid;
  negative_denominator[47] = 1;
  std::string unreduced = valid;
  unreduced[31] = 5;
  return rejected(truncated) && rejected(trailing) &&
         rejected(wrong_child_kind) && rejected(invalid_sign) &&
         rejected(nonminimal_magnitude) && rejected(negative_denominator) &&
         rejected(unreduced);
}

std::string boundary_comparison(const std::string &boundary,
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

bool exact_classes_and_certificates_are_stable() {
  const MatClearanceProfileGraph2 bundle = two_width_bundle();
  const std::vector<MatNeckEvidenceV1> evidence =
      exact_neck_evidence_v1(bundle);
  const std::vector<std::string> records = record_bytes(evidence);
  const std::vector<std::string> boundaries{
      canonical_encode_rational(4),
      canonical_encode_rational(6),
      canonical_encode_rational(8),
  };
  const auto [classes, certificates] =
      validate_and_classify_necks_v1(bundle, records, boundaries);
  const auto [repeated_classes, repeated_certificates] =
      validate_and_classify_necks_v1(bundle, records, boundaries);
  if (classes != std::vector<std::int64_t>{0, 2} ||
      repeated_classes != classes || certificates != repeated_certificates ||
      certificates.size() != records.size()) {
    return false;
  }

  const std::string expected_first = canonical_encode_tagged_union(
      "mat-neck-width-comparison-v1",
      canonical_encode_component_map({
          {
              "neck-evidence-digest",
              sha256_bytes(records[0]),
          },
          {
              "squared-width-root-id",
              evidence[0].evidence().squared_width().root_id(),
          },
          {
              "boundary-comparisons",
              canonical_encode_sequence({
                  boundary_comparison(boundaries[0], 0),
                  boundary_comparison(boundaries[1], -1),
                  boundary_comparison(boundaries[2], -1),
              }),
          },
          {
              "class-id",
              canonical_encode_integer(0),
          },
      }));
  return certificates[0] == expected_first &&
         certificates[1].find("mat-neck-width-comparison-v1") !=
             std::string::npos &&
         certificates[1].find(sha256_bytes(records[1])) != std::string::npos;
}

bool empty_boundaries_form_one_class() {
  const MatClearanceProfileGraph2 bundle = two_width_bundle();
  const std::vector<std::string> records =
      record_bytes(exact_neck_evidence_v1(bundle));
  const auto [classes, certificates] =
      validate_and_classify_necks_v1(bundle, records, {});
  return classes == std::vector<std::int64_t>{0, 0} &&
         certificates.size() == records.size();
}

bool malformed_boundaries_fail_loudly() {
  const MatClearanceProfileGraph2 bundle = two_width_bundle();
  const std::vector<std::string> records =
      record_bytes(exact_neck_evidence_v1(bundle));
  const auto rejected = [&bundle,
                         &records](const std::vector<std::string> &boundaries) {
    try {
      static_cast<void>(
          validate_and_classify_necks_v1(bundle, records, boundaries));
    } catch (const InvalidMatNeckWidthBoundariesError &) {
      return true;
    }
    return false;
  };
  return rejected({canonical_encode_rational(-1)}) &&
         rejected({
             canonical_encode_rational(4),
             canonical_encode_rational(4),
         }) &&
         rejected({
             canonical_encode_rational(6),
             canonical_encode_rational(4),
         }) &&
         rejected({canonical_encode_integer(4)});
}

} // namespace

bool neck_classification_gate() {
  return rational_bytes_match_python_golden() &&
         malformed_rational_bytes_fail_loudly() &&
         exact_classes_and_certificates_are_stable() &&
         empty_boundaries_form_one_class() &&
         malformed_boundaries_fail_loudly();
}
