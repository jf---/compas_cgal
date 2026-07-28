#pragma once

#include "segment_site_neck_evidence.h"

#include <stdexcept>
#include <string>
#include <vector>

class InvalidMatNeckEvidenceBytesError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class MatNeckEvidenceV1 {
public:
  static MatNeckEvidenceV1 build(MatExactNeckEvidence2 evidence);

  const MatExactNeckEvidence2 &evidence() const noexcept;
  const std::string &canonical_bytes() const noexcept;
  const std::string &canonical_digest() const noexcept;

private:
  MatNeckEvidenceV1(MatExactNeckEvidence2 evidence, std::string canonical_bytes,
                    std::string canonical_digest);

  MatExactNeckEvidence2 evidence_;
  std::string canonical_bytes_;
  std::string canonical_digest_;
};

std::vector<MatNeckEvidenceV1>
exact_neck_evidence_v1(const MatClearanceProfileGraph2 &bundle);

std::vector<MatNeckEvidenceV1>
verify_neck_evidence_v1(const MatClearanceProfileGraph2 &bundle,
                        const std::vector<std::string> &records);
