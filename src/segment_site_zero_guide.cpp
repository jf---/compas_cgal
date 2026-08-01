#include "segment_site_zero_guide.h"

#include "canonical_encoding.h"

#include <algorithm>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr std::size_t SHA256_DIGEST_SIZE = 32;

void require_digest(const std::string &digest, const char *name) {
  if (digest.size() != SHA256_DIGEST_SIZE) {
    throw InvalidMatZeroGuideRecordError(std::string(name) +
                                         " must contain 32 bytes");
  }
}

} // namespace

MatZeroGuideRecordV1::MatZeroGuideRecordV1(std::string edge_id,
                                           std::string canonical_bytes)
    : edge_id_(std::move(edge_id)),
      canonical_bytes_(std::move(canonical_bytes)) {}

MatZeroGuideRecordV1
MatZeroGuideRecordV1::build(const std::string &mat_certificate_digest,
                            const std::string &center_domain_digest,
                            const MatClearanceEdgeProfile2 &profile,
                            const CORE::BigRat &tool_radius_squared) {
  require_digest(mat_certificate_digest, "MAT certificate digest");
  require_digest(center_domain_digest, "MAT center-domain digest");
  if (tool_radius_squared <= CORE::BigRat(0)) {
    throw InvalidMatZeroGuideRecordError(
        "MAT zero-guide tool-radius square must be positive");
  }
  const RationalPolynomial &coefficients =
      profile.squared_clearance().coefficients();
  if (coefficients.size() != 1 || coefficients.front() != tool_radius_squared) {
    throw InvalidMatZeroGuideRecordError(
        "MAT zero-guide profile must be the constant tool-radius square");
  }

  const CORE::BigRat exact_tool_radius_squared = tool_radius_squared;
  std::string bytes = canonical_encode_tagged_union(
      "segment-site-mat-zero-guide-record-v1",
      canonical_encode_component_map({
          {
              "center-domain-digest",
              canonical_encode_bytes(center_domain_digest),
          },
          {
              "identity-strategy",
              canonical_encode_bytes(
                  "constant-clearance-square-equals-tool-radius-square-v1"),
          },
          {
              "identity-verdict",
              canonical_encode_boolean(true),
          },
          {
              "mat-certificate-digest",
              canonical_encode_bytes(mat_certificate_digest),
          },
          {
              "profile",
              canonical_mat_clearance_profile_v1(profile),
          },
          {
              "tool-radius-squared",
              canonical_encode_rational(exact_tool_radius_squared),
          },
      }));
  return MatZeroGuideRecordV1(profile.edge_id(), std::move(bytes));
}

const std::string &MatZeroGuideRecordV1::edge_id() const noexcept {
  return edge_id_;
}

const std::string &MatZeroGuideRecordV1::canonical_bytes() const noexcept {
  return canonical_bytes_;
}

std::vector<MatZeroGuideRecordV1>
certified_mat_zero_guide_records_v1(const MatClearanceProfileGraph2 &exact,
                                    const MatCertificateV1 &mat_certificate,
                                    const std::string &center_domain_digest,
                                    const MatToolRadiusMm2 &tool_radius) {
  std::vector<MatZeroGuideRecordV1> records;
  for (const MatClearanceEdgeProfile2 &profile : exact.profiles()) {
    const RationalPolynomial &coefficients =
        profile.squared_clearance().coefficients();
    if (coefficients.size() == 1 &&
        coefficients.front() == tool_radius.squared_exact()) {
      records.push_back(MatZeroGuideRecordV1::build(
          mat_certificate.canonical_digest(), center_domain_digest, profile,
          tool_radius.squared_exact()));
    }
  }
  std::sort(
      records.begin(), records.end(),
      [](const MatZeroGuideRecordV1 &left, const MatZeroGuideRecordV1 &right) {
        return left.edge_id() < right.edge_id();
      });
  if (std::adjacent_find(records.begin(), records.end(),
                         [](const MatZeroGuideRecordV1 &left,
                            const MatZeroGuideRecordV1 &right) {
                           return left.edge_id() == right.edge_id();
                         }) != records.end()) {
    throw DuplicateMatZeroGuideEdgeError(
        "MAT zero-guide inventory contains a duplicate edge identity");
  }
  return records;
}

std::vector<std::string> validate_mat_zero_guide_records_v1(
    const MatClearanceProfileGraph2 &exact,
    const CanonicalMatSiteCatalog2 &catalog,
    const std::string &center_domain_digest,
    const MatToolRadiusMm2 &tool_radius, const std::string &mat_certificate,
    const std::vector<std::pair<std::string, std::string>> &records) {
  const MatCertificateV1 replayed = replay_mat_certificate_v1(
      exact, catalog, center_domain_digest, mat_certificate);
  const std::vector<MatZeroGuideRecordV1> expected =
      certified_mat_zero_guide_records_v1(exact, replayed, center_domain_digest,
                                          tool_radius);

  std::map<std::string, std::string> provided_by_edge;
  for (const auto &[edge_id, record] : records) {
    if (!provided_by_edge.emplace(edge_id, record).second) {
      throw DuplicateMatZeroGuideEdgeError(
          "MAT zero-guide records contain a duplicate edge identity");
    }
  }

  std::vector<std::string> verified_edge_ids;
  verified_edge_ids.reserve(expected.size());
  for (const MatZeroGuideRecordV1 &record : expected) {
    const auto provided = provided_by_edge.find(record.edge_id());
    if (provided == provided_by_edge.end()) {
      throw MissingMatZeroGuideEdgeError(
          "MAT zero-guide records are missing a proved edge identity");
    }
    if (provided->second != record.canonical_bytes()) {
      throw MismatchedMatZeroGuideRecordError(
          "MAT zero-guide record bytes disagree with exact replay");
    }
    verified_edge_ids.push_back(record.edge_id());
  }
  if (provided_by_edge.size() != expected.size()) {
    throw MismatchedMatZeroGuideRecordError(
        "MAT zero-guide records contain a foreign edge identity");
  }
  return verified_edge_ids;
}
