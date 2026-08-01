#pragma once

#include "segment_site_mat_certificate.h"
#include "segment_site_mat_proposal_table.h"

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

class InvalidMatZeroGuideRecordError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class MissingMatZeroGuideEdgeError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class DuplicateMatZeroGuideEdgeError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class MismatchedMatZeroGuideRecordError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class MatZeroGuideRecordV1 {
public:
  static MatZeroGuideRecordV1 build(const std::string &mat_certificate_digest,
                                    const std::string &center_domain_digest,
                                    const MatClearanceEdgeProfile2 &profile,
                                    const CORE::BigRat &tool_radius_squared);

  const std::string &edge_id() const noexcept;
  const std::string &canonical_bytes() const noexcept;

private:
  MatZeroGuideRecordV1(std::string edge_id, std::string canonical_bytes);

  std::string edge_id_;
  std::string canonical_bytes_;
};

std::vector<MatZeroGuideRecordV1>
certified_mat_zero_guide_records_v1(const MatClearanceProfileGraph2 &exact,
                                    const MatCertificateV1 &mat_certificate,
                                    const std::string &center_domain_digest,
                                    const MatToolRadiusMm2 &tool_radius);

std::vector<std::string> validate_mat_zero_guide_records_v1(
    const MatClearanceProfileGraph2 &exact,
    const CanonicalMatSiteCatalog2 &catalog,
    const std::string &center_domain_digest,
    const MatToolRadiusMm2 &tool_radius, const std::string &mat_certificate,
    const std::vector<std::pair<std::string, std::string>> &records);
