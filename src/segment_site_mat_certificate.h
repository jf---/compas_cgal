#pragma once

#include "segment_site_catalog.h"
#include "segment_site_catalog_neck.h"
#include "segment_site_mat_numeric_table.h"

#include <stdexcept>
#include <string>
#include <vector>

class InvalidMatCenterDomainDigestError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class InvalidMatCertificateReplayError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class InvalidMatCertificateEncodingError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class MatCertificateV1;
struct MatCertifiedExactProjection2;

std::string canonical_mat_clearance_profile_v1(
    const MatClearanceEdgeProfile2 &profile);

class MatCertificateV1 {
public:
  const std::string &canonical_bytes() const noexcept;
  const std::string &canonical_digest() const noexcept;
  std::string release_canonical_bytes() && noexcept;

private:
  MatCertificateV1(std::string canonical_bytes, std::string canonical_digest);

  std::string canonical_bytes_;
  std::string canonical_digest_;

  friend MatCertifiedExactProjection2
  certified_mat_exact_projection_v1(const MatClearanceProfileGraph2 &exact,
                                    const CanonicalMatSiteCatalog2 &catalog,
                                    const std::string &center_domain_digest);
};

struct MatCertifiedExactProjection2 {
  MatNumericGraphTable2 graph;
  MatNumericNeckCutTable2 necks;
  std::vector<MatNeckEvidenceV1> neck_evidence;
  MatCertificateV1 certificate;
};

MatCertifiedExactProjection2
certified_mat_exact_projection_v1(const MatClearanceProfileGraph2 &exact,
                                  const CanonicalMatSiteCatalog2 &catalog,
                                  const std::string &center_domain_digest);

MatCertificateV1
replay_mat_certificate_v1(const MatClearanceProfileGraph2 &exact,
                          const CanonicalMatSiteCatalog2 &catalog,
                          const std::string &center_domain_digest,
                          const std::string &record);
